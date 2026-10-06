import { spawn } from 'node:child_process';
import type { ChildProcessWithoutNullStreams } from 'node:child_process';
import { existsSync } from 'node:fs';
import path from 'node:path';
import readline from 'node:readline';

/**
 * The Electron side of the Python service boundary.
 *
 * Runs `python -m molecular_diagnosis.service` as a CHILD PROCESS and talks to
 * it in newline-delimited JSON. Nothing about the scientific code is embedded
 * in Electron: this module knows method names and plain data, nothing else.
 *
 * Why a child process rather than an embedded interpreter:
 *   - a long DMC search cannot block the Electron main process or the renderer;
 *   - a crash in the analysis cannot take the window down;
 *   - the same service can be driven from a terminal for debugging.
 *
 * stdout carries the protocol and nothing else. The child redirects its own
 * `print` to stderr for exactly this reason; on this side, stderr is collected
 * as diagnostics and attached to failures instead of being parsed.
 *
 * Three kinds of line arrive on stdout, distinguished by shape and never by
 * position: a response carries `ok`, a progress notification carries
 * `type: "progress"`, and the one-off startup handshake has id `ready`. A
 * progress line leaves its request PENDING — that is the whole point of it —
 * so the framing is unchanged and a reader that ignores progress still works.
 *
 * stderr MUST keep being read. The child writes structured diagnostics there
 * for every request; a parent that pipes stderr and stops draining it fills
 * the OS pipe buffer, and the child then blocks mid-write and never answers.
 * That failure looks exactly like a hung analysis, which is why the readline
 * interface below is not optional.
 */

export interface BackendError {
  readonly code: string;
  readonly message: string;
  readonly detail?: string;
  /** Python traceback, kept for logs. Never required for the UI to behave. */
  readonly traceback?: string;
}

export type BackendResult<T> =
  | { readonly ok: true; readonly result: T }
  | { readonly ok: false; readonly error: BackendError };

interface PendingRequest {
  readonly resolve: (value: BackendResult<unknown>) => void;
  readonly method: string;
}

/** A progress notification, with the id of the request it belongs to. */
export interface BridgeProgress {
  readonly requestId: string;
  readonly method: string;
  readonly progress: Record<string, unknown>;
}

/** Lines kept from stderr so a failure can carry recent diagnostics. */
const STDERR_RING_SIZE = 60;

/** Startup handshake budget. The analysis itself is deliberately untimed. */
const READY_TIMEOUT_MS = 20_000;

export class PythonBridge {
  private child: ChildProcessWithoutNullStreams | null = null;
  private ready: Promise<void> | null = null;
  private readonly pending = new Map<string, PendingRequest>();
  private readonly stderrLines: string[] = [];
  private nextId = 1;

  /** Called for every progress notification. Set by the owner. */
  private onProgress: ((progress: BridgeProgress) => void) | null = null;

  constructor(
    private readonly repoRoot: string,
    private readonly log: (message: string) => void = () => undefined,
  ) {}

  /**
   * Listen for progress from long-running requests.
   *
   * One listener, owned by the process that created the bridge. Progress is
   * pushed as it arrives; the request it belongs to is still pending and will
   * resolve separately.
   */
  setProgressListener(listener: ((progress: BridgeProgress) => void) | null): void {
    this.onProgress = listener;
  }

  /**
   * Resolve the interpreter to run.
   *
   * Order: explicit override, then the project virtualenv, then whatever
   * `python` is on PATH. The venv is preferred because it is the environment
   * the project's dependencies (openpyxl) are installed into.
   */
  private resolvePython(): string {
    const override = process.env.MOLECULAR_TOOL_PYTHON;
    if (override && existsSync(override)) return override;

    const candidates =
      process.platform === 'win32'
        ? [path.join(this.repoRoot, '.venv', 'Scripts', 'python.exe')]
        : [path.join(this.repoRoot, '.venv', 'bin', 'python3'), path.join(this.repoRoot, '.venv', 'bin', 'python')];

    for (const candidate of candidates) {
      if (existsSync(candidate)) return candidate;
    }

    return process.platform === 'win32' ? 'python' : 'python3';
  }

  private start(): Promise<void> {
    if (this.ready) return this.ready;

    this.ready = new Promise<void>((resolve, reject) => {
      const python = this.resolvePython();
      this.log(`starting backend: ${python} -m molecular_diagnosis.service (cwd ${this.repoRoot})`);

      let child: ChildProcessWithoutNullStreams;
      try {
        child = spawn(python, ['-u', '-m', 'molecular_diagnosis.service'], {
          cwd: this.repoRoot,
          stdio: ['pipe', 'pipe', 'pipe'],
          env: { ...process.env, PYTHONIOENCODING: 'utf-8' },
        });
      } catch (error) {
        reject(error instanceof Error ? error : new Error(String(error)));
        return;
      }

      this.child = child;

      const stdout = readline.createInterface({ input: child.stdout, crlfDelay: Infinity });
      stdout.on('line', (line) => this.handleLine(line, resolve));

      const stderr = readline.createInterface({ input: child.stderr, crlfDelay: Infinity });
      stderr.on('line', (line) => {
        // Diagnostics channel. Never parsed, only recorded and logged.
        this.stderrLines.push(line);
        if (this.stderrLines.length > STDERR_RING_SIZE) this.stderrLines.shift();
        this.log(`[python] ${line}`);
      });

      child.on('error', (error) => {
        this.teardown(`failed to start: ${error.message}`);
        reject(error);
      });

      child.on('exit', (code, signal) => {
        this.teardown(`backend exited (code ${code ?? 'null'}, signal ${signal ?? 'none'})`);
      });

      const timer = setTimeout(() => {
        reject(new Error('The analysis backend did not start in time.'));
      }, READY_TIMEOUT_MS);

      // The child announces itself with an id of "ready" before anything else.
      this.readyResolvers.push(() => {
        clearTimeout(timer);
        resolve();
      });
    });

    // A failed start must not be cached, or every later call fails too.
    this.ready = this.ready.catch((error) => {
      this.ready = null;
      throw error;
    });

    return this.ready;
  }

  private readonly readyResolvers: Array<() => void> = [];

  private handleLine(line: string, _resolveReady: () => void): void {
    const trimmed = line.trim();
    if (!trimmed) return;

    let message: unknown;
    try {
      message = JSON.parse(trimmed);
    } catch {
      // Not protocol JSON. The child guards against this, so treat it as a
      // diagnostic rather than failing the in-flight request.
      this.log(`[python:non-json] ${trimmed}`);
      return;
    }

    if (typeof message !== 'object' || message === null) return;
    const envelope = message as {
      id?: unknown;
      ok?: unknown;
      result?: unknown;
      error?: unknown;
      type?: unknown;
      progress?: unknown;
    };

    if (envelope.id === 'ready') {
      while (this.readyResolvers.length > 0) this.readyResolvers.pop()?.();
      return;
    }

    if (typeof envelope.id !== 'string') return;
    const pending = this.pending.get(envelope.id);

    /*
     * A progress notification. It does NOT resolve the request: the run is
     * still going, and its response will arrive on a later line.
     */
    if (envelope.type === 'progress') {
      if (!pending) return; // stale: the request already finished or was lost
      const progress = envelope.progress;
      if (typeof progress !== 'object' || progress === null) return;
      this.onProgress?.({
        requestId: envelope.id,
        method: pending.method,
        progress: progress as Record<string, unknown>,
      });
      return;
    }

    if (!pending) return;
    this.pending.delete(envelope.id);

    if (envelope.ok === true) {
      pending.resolve({ ok: true, result: envelope.result });
    } else {
      const error = (envelope.error ?? {}) as Partial<BackendError>;
      pending.resolve({
        ok: false,
        error: {
          code: error.code ?? 'UNKNOWN',
          message: error.message ?? 'The analysis backend reported an error.',
          detail: error.detail,
          traceback: error.traceback,
        },
      });
    }
  }

  /** Fail everything in flight; a dead child will never answer them. */
  private teardown(reason: string): void {
    this.log(reason);
    const diagnostics = this.stderrLines.slice(-12).join('\n');

    for (const [id, pending] of this.pending) {
      this.pending.delete(id);
      pending.resolve({
        ok: false,
        error: {
          code: 'BACKEND_UNAVAILABLE',
          message: 'The analysis backend stopped unexpectedly.',
          detail: `${reason} (during ${pending.method})`,
          traceback: diagnostics || undefined,
        },
      });
    }

    this.child = null;
    this.ready = null;
    this.readyResolvers.length = 0;
  }

  /**
   * Send one request and await its response.
   *
   * Never rejects for a backend-level failure: those come back as
   * `{ ok: false, error }` so callers handle them as data rather than as
   * exceptions. It only rejects if the child cannot be started at all.
   */
  async call<T = unknown>(method: string, params: Record<string, unknown>): Promise<BackendResult<T>> {
    try {
      await this.start();
    } catch (error) {
      return {
        ok: false,
        error: {
          code: 'BACKEND_UNAVAILABLE',
          message:
            'The analysis backend could not be started. Check that Python and the project ' +
            'dependencies are installed.',
          detail: error instanceof Error ? error.message : String(error),
          traceback: this.stderrLines.slice(-12).join('\n') || undefined,
        },
      };
    }

    const child = this.child;
    if (!child || child.stdin.destroyed) {
      return {
        ok: false,
        error: { code: 'BACKEND_UNAVAILABLE', message: 'The analysis backend is not running.' },
      };
    }

    const id = String(this.nextId++);

    return new Promise<BackendResult<T>>((resolve) => {
      this.pending.set(id, {
        method,
        resolve: resolve as (value: BackendResult<unknown>) => void,
      });
      child.stdin.write(`${JSON.stringify({ id, method, params })}\n`);
    });
  }

  dispose(): void {
    const child = this.child;
    this.teardown('backend shutting down');
    if (child && !child.killed) {
      child.stdin.end();
      child.kill();
    }
  }
}

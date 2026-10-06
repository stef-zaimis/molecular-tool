import { useEffect, useState } from 'react';
import type {
  DiagnosisProgress,
  DiagnosisProgressStage,
  DiagnosisResumeState,
} from '../../backendContract';
import type { DiagnosisRunState } from '../../app/state/projectState';
import type { RunGate } from '../../app/state/runGate';
import './DiagnosisRunPanel.css';

/**
 * What each backend stage is called on screen.
 *
 * EXHAUSTIVE over `DiagnosisProgressStage`, deliberately: a stage added in
 * `molecular_diagnosis/progress.py` and to the union will not compile until it
 * has wording here, so the status line can never go blank because the two
 * sides drifted.
 */
const STAGE_LABELS: Record<DiagnosisProgressStage, string> = {
  starting: 'starting',
  verifying_sources: 'verifying FASTA',
  loading_alignment: 'reading alignment',
  building_scope: 'combining alignments',
  validating_focal: 'checking focal entries',
  consensus: 'focal consensus',
  dmc_search: 'DMC search',
  five_site: '5-site optimisation',
  writing_report: 'writing report',
  writing_consensus: 'writing consensus report',
  writing_workbook: 'writing comparison workbook',
  finishing: 'preparing output files',
};

/** `01:23`, or `1:02:03` once a run passes an hour. */
function formatElapsed(ms: number): string {
  const total = Math.max(0, Math.floor(ms / 1000));
  const seconds = String(total % 60).padStart(2, '0');
  const minutes = Math.floor(total / 60) % 60;
  const hours = Math.floor(total / 3600);
  if (hours > 0) return `${hours}:${String(minutes).padStart(2, '0')}:${seconds}`;
  return `${String(minutes).padStart(2, '0')}:${seconds}`;
}

/**
 * The middle of the status line: what the backend is doing, with its counts.
 *
 * Formatting lives here rather than in Python because it is wording, and the
 * numbers arrive as numbers. A stage with nothing to count is just its label.
 */
function describeProgress(progress: DiagnosisProgress | null): string {
  if (!progress) return 'starting';

  const label = STAGE_LABELS[progress.stage];
  const { current, total, detail, stage } = progress;

  if (stage === 'dmc_search' && detail) {
    const size = `${label}, size ${detail}`;
    return total !== null && current !== null
      ? `${size} — ${current.toLocaleString()} / ${total.toLocaleString()}`
      : size;
  }

  if (current !== null && total !== null) {
    return `${label} ${current.toLocaleString()}/${total.toLocaleString()}`;
  }
  if (current !== null) return `${label} ${current.toLocaleString()}`;
  return label;
}

/**
 * A once-a-second clock, live only while a run is going.
 *
 * The elapsed time has to keep moving even when no progress arrives — a long
 * stage that reports nothing is exactly when the user most needs to see that
 * the app is still alive.
 */
function useElapsed(since: number | null): number {
  const [now, setNow] = useState(() => Date.now());

  useEffect(() => {
    if (since === null) return;
    setNow(Date.now());
    const timer = setInterval(() => setNow(Date.now()), 1000);
    return () => clearInterval(timer);
  }, [since]);

  return since === null ? 0 : Math.max(0, now - since);
}

interface DiagnosisRunPanelProps {
  readonly run: DiagnosisRunState;
  /**
   * Whether a run may start, and why not.
   *
   * Computed in `runGate.ts` from the same facts the backend checks, so the
   * button explains itself before the click instead of the user discovering
   * the refusal afterwards.
   */
  readonly gate: RunGate;
  readonly onRun: () => void;
  readonly onContinue: (resume: DiagnosisResumeState) => void;
  readonly onDismissContinuation: () => void;
  readonly onReveal: (filePath: string) => void;
}

function fileName(filePath: string): string {
  const parts = filePath.split(/[\\/]/);
  return parts[parts.length - 1] || filePath;
}

/**
 * Run controls and results for Molecular Diagnosis.
 *
 * Not in the reference design, which shows no way to start an analysis. It is
 * required to reach the pipeline at all, so it is styled in the existing visual
 * language rather than introducing a new one.
 *
 * The continuation prompt lives here because the *decision* belongs to the user,
 * while the continuation state itself is produced and consumed by the backend —
 * this component only hands the opaque `resume` blob back.
 */
export function DiagnosisRunPanel({
  run,
  gate,
  onRun,
  onContinue,
  onDismissContinuation,
  onReveal,
}: DiagnosisRunPanelProps): JSX.Element {
  const running = run.status === 'running';
  const canRun = gate.canRun;
  const blockedReason = gate.reason ?? undefined;
  const elapsed = useElapsed(run.status === 'running' ? run.startedAt : null);

  return (
    <section className="run-panel" aria-label="Molecular Diagnosis run">
      <div className="run-panel__actions">
        <button
          type="button"
          className="run-panel__run"
          onClick={onRun}
          disabled={!canRun}
          title={blockedReason}
        >
          {running ? 'Running...' : 'Run Molecular Diagnosis'}
        </button>

        {run.status === 'running' && (
          /*
           * The one live line. It says the same three things throughout —
           * running, which stage with its counts, how long so far — so a run
           * that looks stuck can be told apart from one that is simply long.
           */
          <span className="run-panel__status run-panel__status--running" role="status">
            <span className="run-panel__progress-state">
              {run.continuing ? 'Continuing' : 'Running'}
            </span>
            <span className="run-panel__progress-stage">{describeProgress(run.progress)}</span>
            <span className="run-panel__progress-clock">{formatElapsed(elapsed)}</span>
          </span>
        )}

        {blockedReason && !running && <span className="run-panel__status">{blockedReason}</span>}

        {run.status === 'succeeded' && (
          <span className="run-panel__status">
            {run.result.sequenceCount} sequences · alignment length{' '}
            {run.result.alignmentLength ?? '—'}
          </span>
        )}
      </div>

      {run.status === 'failed' && (
        <div className="run-panel__result run-panel__result--error" role="alert">
          <p className="run-panel__headline">{run.error.message}</p>
          {run.error.detail && <p className="run-panel__detail">{run.error.detail}</p>}
        </div>
      )}

      {run.status === 'succeeded' && (
        <div className="run-panel__result" role="status">
          <p className="run-panel__headline">
            {summarise(run.result.dmc.combinationsByLength)} · {run.result.dmc.uniqueSites.length}{' '}
            unique diagnostic sites · stopped at size {run.result.dmc.stoppedAtLength}
          </p>

          <ul className="run-panel__outputs">
            {[
              run.result.outputs.reportTxt,
              run.result.outputs.workbookXlsx,
              run.result.outputs.consensusTxt,
            ]
              .filter((value): value is string => typeof value === 'string' && value.length > 0)
              .map((filePath) => (
                <li key={filePath}>
                  <button
                    type="button"
                    className="run-panel__output"
                    onClick={() => onReveal(filePath)}
                    title={filePath}
                  >
                    {fileName(filePath)}
                  </button>
                </li>
              ))}
          </ul>

          {run.pendingContinuation && (
            <div className="run-panel__continuation">
              <p className="run-panel__detail">
                The search reached the maximum candidate size ({run.result.dmc.maxCombinationLength})
                without a stopping result. It can continue from size{' '}
                {run.pendingContinuation.startCombinationLength}, keeping everything already tested.
              </p>
              <div className="run-panel__continuation-actions">
                <button
                  type="button"
                  className="run-panel__run run-panel__run--small"
                  onClick={() => onContinue(run.pendingContinuation as DiagnosisResumeState)}
                >
                  Continue from size {run.pendingContinuation.startCombinationLength}
                </button>
                <button
                  type="button"
                  className="run-panel__dismiss"
                  onClick={onDismissContinuation}
                >
                  Stop here
                </button>
              </div>
              <p className="run-panel__hint">
                Continuing uses the maximum candidate size set above, so raise it first if you want
                the search to go further than one more size.
              </p>
            </div>
          )}
        </div>
      )}
    </section>
  );
}

/** "2 one-site, 5 two-site" style summary of what the search found. */
function summarise(byLength: Readonly<Record<string, readonly (readonly number[])[]>>): string {
  const parts = Object.entries(byLength)
    .map(([length, combos]) => [Number(length), combos.length] as const)
    .sort((a, b) => a[0] - b[0])
    .map(([length, count]) => `${count} × ${length}-site`);

  return parts.length > 0 ? `Found ${parts.join(', ')}` : 'No diagnostic combinations found';
}

# Molecular Diagnosis Tool

Identifies molecular diagnostic characters (DMCs) from aligned FASTA files.

The repository holds three layers:

- **the scientific core** — `molecular_diagnosis/`, the DMC search, punishment
  scoring, consensus generation and the report/Excel writers;
- **the project layer** — `molecular_diagnosis/project/` and
  `molecular_diagnosis/service/`: a SQLite-backed project format (linked FASTA
  sources, persistent focal sets) exposed over a newline-delimited JSON service;
- **the desktop app** — `desktop/`, an Electron + React UI that talks to that
  service and never imports the scientific code itself.

`REPO_MAP.md` maps the whole tree; `desktop/UI_NOTES.md` records the UI's
decisions and open questions.

## Features

- Reads and validates aligned FASTA files
- Identifies focal-set fixed sites and filters globally conserved ones
- Reports n-site diagnostic characters, with a resumable search
- Searches 5-site combinations for similarity optimisation
- Focal-only punishment/anomaly scoring and focal consensus generation
- Exports a text report and a formatted Excel workbook
- Persistent projects: linked FASTA sources with live status, and saved focal
  sets

## Installation

Create and activate a virtual environment, then install the Python
dependencies:

```bash
pip install -r requirements.txt
```

For the desktop app, install its own dependencies:

```bash
cd desktop && npm install
```

## Running

The desktop app (current UI):

```bash
cd desktop && npm start
```

It spawns `python -m molecular_diagnosis.service` as a child process. The
interpreter is `MOLECULAR_TOOL_PYTHON` if set, otherwise the project `.venv`,
otherwise `python`/`python3` on PATH.

The legacy Tkinter UI is still runnable and has no project database:

```bash
python main.py
```

## Tests

```bash
python -m pytest
```

```bash
cd desktop && npm test && npm run typecheck && npm run lint
```

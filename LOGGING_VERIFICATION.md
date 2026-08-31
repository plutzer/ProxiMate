# Logging and provenance: verification checklist

ProxiMate emits four layers of diagnostics. This document records what each
layer guarantees, what the automated suite already covers, and the checks that
need a running container.

| Layer | Destination | Who writes it |
| --- | --- | --- |
| Operational stream | stderr (`docker logs`) + `$PROXIMATE_LOG_DIR/proximate-server.log`, rotated | `log_config` |
| Provenance log | `<dataset>/proximate.log` | `log_config.add_file_handler` / `dataset_log` |
| Run manifest | `<dataset>/run.json` | `provenance.stage` |
| User notifications | Shiny toasts | `notify()` in `GUI/app.py` |

A **run ID** minted once per action is inherited by every child process through
`PROXIMATE_RUN_ID` and stamped on every record, so one run can be followed across
the GUI and the three stages it launches.

## Automated coverage

Run from the repository root:

```
python -m pytest
```

| File | Covers |
| --- | --- |
| `Scripts/tests/test_log_config.py` | logger scoping, handler de-duplication and ref-counting, `dataset_log`, run IDs and `run_context`, rotation, level handling |
| `Scripts/tests/test_provenance.py` | manifest schema, stage success/failure/`SystemExit`, reruns, checksums, atomic and corrupt-file behavior |
| `Scripts/tests/test_pipeline_logging.py` | static guards: no `print` in pipeline modules, every one uses `log_config`, `ui.notification_show` only inside `notify`, no bare `except`, no `traceback.print_exc` |

`GUI/app.py` cannot be imported by the suite — it needs `shiny`, which the
analysis environment does not carry, and `app.py` hardcodes `sys.path.append('/Scripts')`.
Everything GUI-side is therefore checked by the static guards above plus the
manual list below.

## Manual checks

Build with a version stamp, since the image carries no git repository:

```
docker build --build-arg PROXIMATE_VERSION=$(git rev-parse --short HEAD) -t proximate:log .
docker run -p 3838:3838 --mount type=bind,source=$PWD/out,target=/Outputs proximate:log
```

### Startup and the operational log

| # | Step | Expected |
| --- | --- | --- |
| 1 | Watch `docker logs -f` at startup | Lines shaped `<time> [INFO   ] <run-id> pid=<n> proximate.app: ...` |
| 2 | `ls out/` on the host | `proximate-server.log` exists and grows |
| 3 | `docker restart` the container, then `ls -l out/proximate-server.log` | Appended to, not truncated — earlier lines survive |
| 4 | Run with `-e LOG_LEVEL=DEBUG` | ProxiMate DEBUG lines appear, with **no** matplotlib / urllib3 / py4cytoscape flood |
| 5 | Run with `-e LOG_LEVEL=NONSENSE` | One warning naming `NONSENSE`, then normal INFO logging; the app still starts |

### Parsing

| # | Step | Expected |
| --- | --- | --- |
| 6 | Parse a MaxQuant dataset `A` | Lines stream to the terminal during the parse; `out/A/proximate.log` and `out/A/run.json` both exist |
| 7 | Inspect `out/A/run.json` | One run, one `parse` stage, `status: "ok"`, input `sha256`s, `n_experiments` / `n_controls`, and a non-null `proximate.version` |
| 8 | Parse a second dataset `B`, then re-check `out/A/proximate.log` | **A's log gained no new lines.** This is the cross-dataset contamination check — before `dataset_log`, every later GUI line was written into every previously parsed dataset's log |
| 9 | Upload a malformed ED file (e.g. missing `Type`) | Actionable toast **and** a matching ERROR in `docker logs` |
| 10 | Upload a bait TSV with the wrong column count | Toast naming the problem; the app does not silently keep an empty table |
| 11 | Parse with a corrupt proteinGroups file | Toast names the run ID and points at `proximate.log`, *not* "check the console"; the full traceback is in the log |

### Scoring

| # | Step | Expected |
| --- | --- | --- |
| 12 | Score dataset `A`, watching `docker logs -f` | Output appears **during** the run, not in one burst at the end. With imputation 2 or 3, periodic `imputed N/M preys` lines |
| 13 | Compare `docker logs` against `out/A/proximate.log` | Each `score.py` line appears **once**, not twice |
| 14 | Inspect `out/A/run.json` after scoring | Three stages (`parse`, `score`, `annotate`) sharing one `run_id`; `score.extra.saint_binary` names the binary that ran; `cli_args.seed` present |
| 15 | Score with quant type `Spectral Counts`; then Intensity + imputation 2; then Intensity + imputation 0 | `saint_binary` is `/bin/SAINTexpress-spc`, then `/bin/SAINTexpress-int`, then `/bin/SAINTexpress-int_default` — all three branches of `_build_saint_cmd` |
| 16 | Delete `out/A/to_CompPASS.csv`, then score | Toast shows *this run's* log tail naming the missing file; `run.json` gains a `score` stage with `status: "error"` and `exit_code: 1` |
| 17 | Make scoring fail before its log handler attaches (e.g. `chmod 000` the dataset dir) | Toast falls back to "exit N; check the container logs" rather than showing a previous run's text |
| 18 | Score the same dataset three times in one container session | No duplicated lines in `proximate.log` — the handler de-duplication check |

### Downloads and housekeeping

| # | Step | Expected |
| --- | --- | --- |
| 18b | Compare `scored_rows_in` against `annotated_rows` in `run.json` | Equal. Annotation joins reference databases by gene name and accession; a duplicated key there multiplies interaction rows, which reads downstream as extra evidence rather than as a join fault |
| 19 | Download the batch ZIP for a scored dataset | Contains `run.json`, `proximate.log` and `build_info.txt` alongside the existing seven files |
| 20 | Click "Clear datasets" | Dataset directories go; `proximate-server.log` **survives** and the GUI keeps logging normally afterwards |
| 21 | Select a dataset that has not been scored | Selects reset as before; at `LOG_LEVEL=DEBUG` a line explains why |
| 22 | Open the Cytoscape tab with Cytoscape not running | Same status text as before, plus a logged line |

### Command-line pipeline

| # | Step | Expected |
| --- | --- | --- |
| 23 | `./run_pipeline.sh --format maxquant` with a **nonexistent** PG file | Stops at the parse stage with a non-zero `$?`; scoring and annotation never run. **Previously the run continued and produced plausible-looking garbage** |
| 24 | A valid run | `log.txt` opens with a run banner carrying the timestamp, run ID and full command line, and now contains the `=== Parsing ===` stage banners |
| 25 | Run twice into the same output directory | Two separated banners in `log.txt`; `run.json` holds two runs, neither clobbering the other |
| 26 | `--format saint` and `--format maxquant` | Both now produce `proximate.log` and `run.json`; neither did before |
| 27 | `--format maxquant` with a malformed ED | Rejected at parse. The CLI MaxQuant path now runs the same validation as the GUI, so inputs it previously accepted silently are refused |
| 28 | Run twice with the same `--seed` | Identical `WD_pval` columns in `compPASS.csv`, and the seed recorded in `run.json` both times |

## Behavior changes to expect

These are intended, and visible:

- **`run_pipeline.sh` now stops at the first failing stage.** Runs that used to
  limp to completion on missing inputs will now abort. This surfaces failures
  that were previously silent, so a dataset that "worked" before may now fail —
  it was producing junk.
- **The MaxQuant CLI path validates its inputs**, because it now calls
  `parse_ed_pg` instead of duplicating a subset of it. It also copies
  `proteinGroups.txt` into the output directory, as the GUI path always did.
- **A parse that fails validation now leaves an output directory** containing a
  `run.json` recording the failure, where previously it left nothing.
- **Error toasts are sourced from the dataset log tail**, so their wording differs.
- **The WDFDR warning can now appear as a read error.** That check previously
  swallowed a failing `read_csv` and rendered as "nothing to warn about".
- **Log lines gained a run ID and PID**, and logger names are prefixed
  `proximate.` — anything that greps the logs needs updating.

## Out of scope

- `Scripts/setup_datasets.py` keeps its own `log()`. It runs at image build time
  with no output directory and no run context, where unbuffered stderr is the
  right behavior.
- `Scripts/saint_normalization.py`'s `print` calls are inside an
  `if __name__ == "__main__"` demonstration block and are not on any pipeline path.
- `SAINTexpress-int_oldimp` is still built by the Dockerfile but referenced by no
  code. Removing it is a separate cleanup.
- `upload_session` extracts a user-supplied zip with `extractall`, which is a
  path-traversal sink. Hardening it is a separate change.

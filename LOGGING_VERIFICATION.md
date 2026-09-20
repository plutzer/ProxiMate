# Logging and provenance: verification checklist

ProxiMate emits four layers of diagnostics.

| Layer | Destination | Who writes it |
| --- | --- | --- |
| Operational stream | stderr (`docker logs`) + `$PROXIMATE_LOG_DIR/proximate-server.log`, rotated | `log_config` |
| Provenance log | `<dataset>/proximate.log` | `log_config.add_file_handler` / `dataset_log` |
| Run manifest | `<dataset>/run.json` | `provenance.stage` |
| User notifications | Shiny toasts | `notify()` in `GUI/app.py` |

A **run ID** minted once per action is inherited by every child process through
`PROXIMATE_RUN_ID` and stamped on every record, so one run can be followed across
the GUI and the three stages it launches.

## Manual checks

These need a running container. Build with a version stamp, since the image
carries no git repository:

```
docker build --build-arg PROXIMATE_VERSION=$(git rev-parse --short HEAD) -t proximate:log .
docker run -p 3838:3838 --mount type=bind,source=$PWD/out,target=/Outputs proximate:log
```

| # | Step | Expected |
| --- | --- | --- |
| 8 | Parse dataset `A`, then dataset `B`, then re-check `out/A/proximate.log` | A's log gained no new lines: `dataset_log` detaches each dataset's handler when its action ends |
| 13 | Score `A` and compare `docker logs` against `out/A/proximate.log` | Each `score.py` line appears once, not twice |
| 18 | Score the same dataset three times in one container session | No duplicated lines in `proximate.log`: one handler per file, however many actions attach it |
| 18b | Compare `scored_rows_in` against `annotated_rows` in `run.json` | Equal. Annotation joins reference databases by gene name and accession; a duplicated key there multiplies interaction rows, which reads downstream as extra evidence rather than as a join fault |
| 23 | `./run_pipeline.sh --format maxquant` with a nonexistent PG file | Stops at the parse stage with a non-zero `$?`; scoring and annotation never run |
| 25 | Run the pipeline twice into the same output directory | Two separated banners in `log.txt`; `run.json` holds two runs, neither clobbering the other |

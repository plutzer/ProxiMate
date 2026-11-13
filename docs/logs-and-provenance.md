# Logs and provenance

Every run records what it did. **Attach `run.json` to a bug report.**

## Files

| File | Scope | Contents |
| --- | --- | --- |
| `<output_dir>/run.json` | one dataset | parameters, input checksums, row counts, which SAINTexpress binary ran, package versions, pass/fail per stage |
| `<output_dir>/proximate.log` | one dataset | timestamped log from parsing, scoring and annotation |
| `<output_dir>/log.txt` | one dataset, command line only | raw console transcript, including output the logger never sees |
| `$PROXIMATE_LOG_DIR/proximate-server.log` | all datasets | the server's operational log, rotated at 10 MB |

## Run ID

Every log line and manifest entry carries a **run ID**. The GUI and the parse,
score and annotate processes it launches share it, so one run can be traced
across all four files.

## Environment variables

| Variable | Default | Effect |
| --- | --- | --- |
| `LOG_LEVEL` | `INFO` | `DEBUG` adds per-prey imputation detail. Affects ProxiMate's own loggers only. |
| `PROXIMATE_LOG_DIR` | `/Outputs` (container), the output directory (command line) | Where the operational log is written |
| `PROXIMATE_OUTPUT_DIR` | `/Outputs` | Where the GUI stores datasets |
| `PROXIMATE_RUN_ID` | minted per run | Set it to correlate an external job with a ProxiMate run |

Cytoscape variables are in [Cytoscape setup](cytoscape.md); port variables are in
[Agent access (MCP)](mcp.md#ports).

# Command line

`run_pipeline.sh` parses, scores and annotates one dataset without the GUI. Use it
on a cluster or in a script.

## Run the pipeline

```
docker run --mount type=bind,source=<data_dir>,target=/data plutzer/proximate \
  /bin/bash /run_pipeline.sh [OPTIONS] --format <format> <arguments>
```

Paths in the arguments are paths inside the container, so they start with `/data`.

| `--format` | Arguments, in order |
| --- | --- |
| `maxquant` | `ED_file PG_file quant_type output_dir n_iterations imputation` |
| `diann` | `ED_file matrix_file output_dir n_iterations imputation` |
| `pioneer` | `ED_file protein_groups_wide.tsv output_dir n_iterations imputation` |
| `fragpipe` | `ED_file FP_file quant_type output_dir n_iterations imputation` |
| `msstats` | `ED_file ProteinLevelData.csv output_dir n_iterations imputation` |
| `saint` | `bait_file prey_file interaction_file quant_type output_dir n_iterations imputation` |

Example, MaxQuant LFQ with 1000 WDFDR iterations and no imputation:

```
docker run --mount type=bind,source=/home/me/screen,target=/data plutzer/proximate \
  /bin/bash /run_pipeline.sh --format maxquant \
  /data/ED.csv /data/proteinGroups.txt LFQ /data/out 1000 0
```

The pipeline stops at the first stage that fails. Run `/run_pipeline.sh` with no
arguments to print the full usage.

## Arguments

| Argument | Values |
| --- | --- |
| `quant_type` | `Intensity`, `LFQ` or `"Spectral Counts"`. DIA-NN and Pioneer always use Intensity. |
| `n_iterations` | CompPASS permutations for the WD false discovery rate. `0` skips WDFDR. |
| `imputation` | `0` none, `1` prey-specific AFT, `2` two-component AFT, `3` one-component AFT |

MSstats `ProteinLevelData.csv` is already log2-transformed, normalized and
imputed, so score it with `imputation` `0`.

## Options

Options go before `--format`.

| Option | Effect |
| --- | --- |
| `--organism` | `human` (default), `mouse` or `yeast` |
| `--exclude-hcm` | Annotate against BioGRID without Human Cell Map evidence; human only. See [Databases](databases.md#excluding-human-cell-map-evidence). |
| `--seed` | CompPASS permutation seed, so WD p-values reproduce between runs |
| `--pi-method` | `weighted_average` (default) or `single_bait`; `imputation` `2` only |
| `--pi-bait` | Control bait name; required with `--pi-method single_bait` |
| `--aft-min-obs` | Preys observed in fewer than N runs take the dataset median per-prey SD as their sigma lower bound (`imputation` `2` or `3`). Default `4`; `0` disables. |

## Output

`output_dir` receives the scores (`merged.csv`, `annotated_scores.csv`), the SAINT
input files, and the run's logs and manifest. See
[Logs and provenance](logs-and-provenance.md).

## Sample data

[`examples/`](../examples/README.md) holds small public datasets to try the
pipeline on.

## Open a shell in the container

To run the Python scripts, GOGO or the SAINTexpress binaries by hand:

```
docker run -it --mount type=bind,source=<data_dir>,target=/data plutzer/proximate /bin/bash
```

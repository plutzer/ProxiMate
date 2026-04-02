# ProxiMate
All-in-one GUI and scripts for analyzing proximity labelling data.

## How to run the tool
### Running the GUI locally (easiest)
1. Download and install Docker: https://www.docker.com/products/docker-desktop/
2. Run the pre-built container from dockerhub, exposing the 3838 port
    `docker run -p 3838:3838 plutzer/proximate`
3. Access the GUI through a web browswer at localhost:3838

### Running the docker container interactively (experienced users)
The backend of the application can be accessed interactively by overriding the command to start the shiny app:
1. Start the docker container interactively. In order to analyze files, you'll need to mount a directory to the container:
    `docker run -it --mount type=bind,source=<native_path_to_data_directory>,target=<working_directory_within_container> plutzer/proximate /bin/bash`
2. Within the docker container, python, R, perl, or SAINT scripts can be run manually.

### Running scripts on a high-performance computing cluster (experienced users)
The docker container contains a shell script to run the entire pipeline. Use the `--format` flag to specify your input type:

**MaxQuant:**
```
docker run --mount type=bind,source=<data_dir>,target=<container_dir> plutzer/proximate \
  /bin/bash /run_pipeline.sh --format maxquant \
  <ED_file> <PG_file> <quant_type> <output_dir> <n_iterations> <imputation>
```

**DIA-NN:**
```
docker run --mount type=bind,source=<data_dir>,target=<container_dir> plutzer/proximate \
  /bin/bash /run_pipeline.sh --format diann \
  <ED_file> <matrix_file> <output_dir> <n_iterations> <imputation>
```

**FragPipe:**
```
docker run --mount type=bind,source=<data_dir>,target=<container_dir> plutzer/proximate \
  /bin/bash /run_pipeline.sh --format fragpipe \
  <ED_file> <FP_file> <quant_type> <output_dir> <n_iterations> <imputation>
```

**SAINT:**
```
docker run --mount type=bind,source=<data_dir>,target=<container_dir> plutzer/proximate \
  /bin/bash /run_pipeline.sh --format saint \
  <bait_file> <prey_file> <interaction_file> <quant_type> <output_dir> <n_iterations> <imputation>
```

**Options:**
- `--organism`: `human` (default), `mouse`, or `yeast` — add before `--format` if needed
- `quant_type`: `Intensity`, `LFQ`, or `Spectral Counts` (DIA-NN always uses Intensity; FragPipe supports all three)
- `imputation`: `0` (none), `1` (prey-specific AFT), or `2` (refactored AFT)

## Annotations and Databases:
I will periodically push newer versions of the tool with updated databases. 

Current versions of the databases:

BioGRID: March 21, 2026

UniProt: March 21, 2026

Human Protein Atlas: March 21, 2026

CORUM: March 21, 2026



### Updating databases manually (experienced users):
Databases can be downloaded automatically by running `python3 Scripts/setup_datasets.py --output-dir Datasets`. Use `--skip` to exclude specific databases (e.g., `--skip corum`). See `python3 Scripts/setup_datasets.py --help` for all options.

Alternatively, you can manually assemble the following files in a `/Datasets` subdirectory inside the ProxiMate parent directory and re-build the docker container. File names will need to match or be changed in the Dockerfile before building.


`BIOGRID-ALL.tab3.txt` - downloaded from [BioGRID](https://downloads.thebiogrid.org/BioGRID)

`BIOGRID-MV-Physical.tab3.txt` - downloaded from [BioGRID](https://downloads.thebiogrid.org/BioGRID)

`uniprot_anns.tsv` - tsv-formatted annotations for human proteins from [UniProt](https://www.uniprot.org/uniprotkb?query=%28proteome%3AUP000005640%29&facets=reviewed%3Atrue)

`subcellular_location.tsv` - downloaded from [Human Protein Atlas](https://www.proteinatlas.org/humanproteome/subcellular/data#locations)

`corum_humanComplexes.txt` - downloaded from [CORUM](https://mips.helmholtz-muenchen.de/corum/download): Human Complexes.




## Common Errors

### SAINT
Permission Denied:
    Can occur when python tries to run a subprocess on a directory rather than a file. Is the SAINT directory correct?

Invalid delimiter:
    SAINT throws this error when the file path is incorrect

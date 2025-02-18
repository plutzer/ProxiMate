# ProxiMate
All-in-one GUI and scripts for analyzing proximity labelling data.

## How to run the tool
### Running the GUI locally (easiest)
1. Download and install Docker: https://www.docker.com/products/docker-desktop/
2. Run the pre-built container from dockerhub, exposing the 3838 port
    `docker run -p 3838:3838 plutzer/score_apms`
3. Access the GUI through a web browswer at localhost:3838

### Running the docker container interactively (experienced users)
The backend of the application can be accessed interactively by overriding the command to start the shiny app:
1. Start the docker container interactively. In order to analyze files, you'll need to mount a directory to the container:
    `docker run -it --mount type=bind,source=<native_path_to_data_directory>,target=<working_directory_within_container> plutzer/score_apms /bin/bash`
2. Within the docker container, python, R, perl, or SAINT scripts can be run manually.

### Running scripts on a high-performance computing cluster (experienced users)
The docker container contains a shell script to run the entire pipeline. This can be run by overriding the command to start the shiny app:
1. Run the command:
       `docker run --mount type=bind,source=<native_path_to_data_directory>,target=<working_directory_within_container> plutzer/score_apms /bin/bash /run_pipeline.sh <path_to_ED_file> <path_to_PG_file> <quant_type (Intensity, Spectral Counts, or LFQ)> <path_to_output_directory> <n_iterations_for_WD_scoring> <imputation (1 for yes or 0 for no)>`

## Common Errors

### SAINT
Permission Denied:
    Can occur when python tries to run a subprocess on a directory rather than a file. Is the SAINT directory correct?

Invalid delimiter:
    SAINT throws this error when the file path is incorrect

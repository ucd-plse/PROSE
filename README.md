# Toward Automated Precision Tuning of Weather and Climate Models: A Case Study

This is the executable artifact accompanying the paper: __Toward Automated Precision Tuning of Weather and Climate Models: A Case Study__. The artifact is intended to increase reproducibility and replicability by providing:

1. All of the data resulting from all of the experiments conducted in the evaluation section of the paper.
2. A means of generating interactive html versions of all of the plots in the evaluation section of the paper from the aforementioned data.
3. A means of demonstrating the execution of our methodology by running the same floating-point precision tuning we applied in the evaluation section of the paper to the `funarc` motivating example described in the preliminaries section of the paper within a docker container.
4. A means of reproducing our results by rerunning the exact same precision tuning experiments described in the evaluation section of the paper provided that one has access to the machinery on which we executed the experiments.

## Table of Contents
- __[0]__ A Note on Reproducibility
- __[1]__ The Docker Artifact
    - __[1.1]__ Requirements
    - __[1.2]__ Run the Docker Artifact (~20 minutes)
    - __[1.3]__ Inspecting the Interactive Plots and the Accompanying Data
    - __[1.4]__ Inspecting the `funarc` Tuning Results
- __[2]__ The Supercomputer Artifact
    - __[2.1]__ Requirements
    - __[2.2]__ Run the Supercomputer Artifact (~12 hours)
    - __[2.3]__ Inspecting the Interactive Plots and the Accompanying Data

## [0] A Note on Reproducibilty

In order to fully reproduce the results described in the evaluation section of this paper, one needs access to the same machine on which we conducted these experiments. There are two main reasons for this. __First__, managing a single environment that can seamlessly switch between the long list of dependencies required by each of the three weather and climate models targeted in the paper is a daunting task. It would not have been possible without the sysadmins of the supercomputing facility we had the privilege of working with. __Second__, the experiments in this paper are the result of >50K core hours of computation. It would not have been possible without the amount of hardware made available to us. _If the artifact reviewing committee wants to reproduce the results of this paper, we can provide a contact for the machine's admin that can set up both a user account and a core-hour allocation. Once this is done, Section [2] describes the necessary steps._

Otherwise, we provide a docker artifact described in Section [1] which includes all of the data from the supercomputer experiments, generates interactive html versions of all of the plots in the paper from this data, and also demonstrates the execution of our methodology by running the same floating-point precision tuning we applied to the weather and climate models to the `funarc` motivating example described in the preliminary section of the paper.


## [1] The Docker Artifact

### [1.1] Requirements
- Git (tested with v2.25.1)
- Docker (tested with v24.0.5)
- GNU bash (tested with v5.0.17)
- An application with which to view html files
- ~27 GB of disk space (15.7 GB docker image, 11GB paper data)
- <1 hour (~20 minutes of execution time, the rest for manual inspection)

### [1.2] Run the Docker Artifact (~20 minutes)
First clone the repo via the following:
```
git clone -b paper_artifact git@github.com:ucd-plse/PROSE
```

Then, from the root of the artifact repo, execute:
```
execute_artifact_docker.sh
```

This will:
1. Acquire the docker image
2. Copy the paper data from the image to the root of the artifact repo
3. Run a `generate_figure.py` script for each experiment's data. (If desired, these scripts can be found by running `find submodules/*/experiments -name "generate_figures.py"` from the root of the artifact repo)
4. Run the same floating-point precision tuning we applied in the evaluation section of the paper to the `funarc` motivating example described in the preliminaries section of the paper.

### [1.3] Inspecting the Interactive Plots and the Accompanying Data
The plots will be saved in the root of the artifact repo. Their names will be of the form `figure_#_model` where `#` is the figure number from the paper and `model` is one of `MOM6`, `MPAS-A`, or `ADCIRC`. These interactive plots are generated using the same data used for the figures in the paper. __Accordingly, the markers plotted will be identical to those plotted in the corresponding figure from the paper.__ These html plots are interactive, allowing zooming, panning, and hover-info. Notably, the hover-info allows for the inclusion of extra data that is not possible in the static figures such as the exact percentage of 32-bit values in each variant and the variant number (labeled as "Configuration Number"). In particular, the variant number allows for further exploration of any particular variants :

The relevant paper data is organized thusly:
```
paper_results/
├── {MODEL}_{SECTION_ID}
.   └── prose_logs                  (directory of all generated variants)
.       ├── search_log.txt              (summary of variants explored)
.       ├── {VARIANT_NUMBER}            (directory of variant src)
.       .   ├── config_{LABEL}              (list of var precisions)
.       .   ├── outlog.txt                  (stdout if variant was executable)
.       .   ├── gptl_timing.tar.gz          (raw hotspot timing info if available)
.       .   ├── gptl_subset_info.pckl       (summarized hotspot timing info if available)
.       .   ├── {SRC_FILE_NAME}.F[90]       (transformed src)
.       .   .
.       .   .
.       .   .
```

### [1.4] Inspecting the `funarc` Tuning Results
An interactive plot analogous to those generated above will be saved to the root of the artifact repo as `funarc_search_results.html`. While the motivating example in the paper was an exhaustive search exploring 256 variants, this search uses the same delta-debugging inspired search used for the weather and climate models in the experimental evaluation. In our tests, this search generates something on the order of only 12 variants before discovering an optimal variant. Note in particular:
1. Variant 0 (uniform 64-bit) plotted at the origin, i.e., 0 relative error and 1x speedup, the baseline for the search
2. Variant 1 (uniform 32-bit) plotted near the upper right, i.e., very high relative error and very high speedup
3. A number of optimal variants on the right side of the plot but under the horizontal dotted line representing the tolerable error threshold

As above, this freshly generated data is available for inspection. It is structured similarly to the above:
```
submodules/
├── funarc
.   └── prose_logs                  (directory of all generated variants)
.       ├── search_log.txt              (summary of variants explored)
.       ├── {VARIANT_NUMBER}            (directory of variant src)
.       .   ├── config_{LABEL}              (list of var precisions)
.       .   ├── outlog.txt                  (stdout if variant was executable)
.       .   ├── {SRC_FILE_NAME}.f90         (transformed src)
.       .   .
.       .   .
.       .   .
```

## [2] The Supercomputer Artifact
### [2.1] Requirements
- A user account and a core-hour allocation must be established on the supercomputer we used for our case study. If desired, we can provide a contact for the machine's admin that can set these up for reproducibilty reviewers.
- ~12 hours for the execution + the time for establishing a user account
- An application with which to view html files


### [2.2] Run the Supercomputer Artifact (~12 hours)
Once logged into the machine, clone the repo via the following:
```
git clone -b paper_artifact git@github.com:ucd-plse/PROSE
```

Then, from the root of the artifact repo, execute:
```
execute_artifact_supercomputer.sh
```

This will run the exact same precision tuning experiments described in the evaluation section of the paper and will generate interactive plots corresponding to all figures in the evaluation section of the paper.

### [2.3] Inspecting the Interactive Plots and the Accompanying Data
The plots will be saved in the root of the artifact repo. Their names will be of the form `figure_#_model` where `#` is the figure number from the paper and `model` is one of `MOM6`, `MPAS-A`, or `ADCIRC`. These html plots are interactive, allowing zooming, panning, and hover-info. Notably, the hover-info allows for the inclusion of extra data that is not possible in the static figures such as the exact percentage of 32-bit values in each variant and the variant number (labeled as "Configuration Number"). In particular, the variant number allows for further exploration of any particular variants :

The relevant paper data is organized thusly:
```
submodules/
├── {MODEL}_{SECTION_ID}
.   ├── experiments
.   .   ├── {EXPERIMENT_NAME}
.   .   .   └── prose_logs                  (directory of all generated variants)
.   .   .       ├── search_log.txt              (summary of variants explored)
.   .   .       ├── {VARIANT_NUMBER}            (directory of variant src)
.   .   .       .   ├── config_{LABEL}              (list of var precisions)
.   .   .       .   ├── outlog.txt                  (stdout if variant was executable)
.   .   .       .   ├── gptl_timing.tar.gz          (raw hotspot timing info if available)
.   .   .       .   ├── gptl_subset_info.pckl       (summarized hotspot timing info if available)
.   .   .       .   ├── {SRC_FILE_NAME}.F[90]       (transformed src)
.   .   .       .   .
.   .   .       .   .
.   .   .       .   .
```

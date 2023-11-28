# Pipeline for sorting mixed samples into detailed groups
Nov. 28, 2023

---

## Data preparation
The input data should be a mixed sample or a testing artificial mixed sample. The target sequencing should use Angiosperm353 as the baits. The testing artificial sample from the Kew should use the method PAFTOL/SRA(Need check)/GAP/Brassiworld.

If the Kew data was used, the `.fastp.gz` file should be merged into one before trimming.

Git clone this directory to the working directory. Then, go to `full_run` folder to view the code and pipeline.

The selected panel can be generated from the following command lines:
```bash
# prepare a gene.list.txt to the name of the original per gene alignment
DIRECTORY=
NAMELIST=
cd full_run
mkdir ref
while read line; do python pick_match_list.py ${DIRECTORY}/${line} ref/${line} ${NAMELIST}; done < gene.list.txt
```

## Run the pipeline
In the `full_run` folder, there are two submission scripts. There are some constants need to be edited before submission. Also, check the python environment if any packages need to be installed. If everything is prepared and ready to go, then just simply submit as:
```bash
sbatch 01-target.sh
sbatch 02-ref_panel.sh
```
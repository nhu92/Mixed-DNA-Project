# Pipeline for sorting mixed samples into detailed groups
Nov. 28, 2023

---

## Data preparation
The input data should be a mixed sample or a testing artificial mixed sample. The target sequencing should use Angiosperm353 as the baits. The testing artificial sample from the Kew should use the method PAFTOL/SRA(Need check)/GAP/Brassiworld.

If the Kew data was used, the `.fastp.gz` file should be merged into one before trimming.

Git clone this directory to the working directory. Then, go to `full_run` folder to view the code and pipeline.
```bash
git clone -n  https://github.com/gudusanjiao/Mixed-DNA-Project.git --depth=1
cd Mixed-DNA-Project/
git checkout HEAD pipeline/full_run
cd ..
cp -r Mixed-DNA-Project/pipeline/full_run/ ./
rm -rf Mixed-DNA-Project/
```

The selected panel can be generated from the following command lines:
```bash
# prepare a gene.list.txt to the name of the original per gene alignment (included as default)
# prepare a list file that contains the species selected from the large reference (70 or 700+)
DIRECTORY=
NAMELIST=
mkdir ref
while read line; do python pick_match_list.py ${DIRECTORY}/${line} ref/${line} ${NAMELIST}; done < gene.list.txt
```

## Run the pipeline
In the `full_run` folder, there are two submission scripts. There are some constants need to be edited before submission. Also, check the python environment if any packages need to be installed. If everything is prepared and ready to go, then just simply submit as:
```bash
mkdir log
sbatch 01-target.sh
sbatch 02-ref_panel.sh
```

One thing need to be make sure is that the `02-ref_panel.sh` must finish before the `01-target.sh` (normally will) because the later one has the commandline to merge the output and generate the contour map. It only happens when job submission system is jammed thus the target pipeline runs way ahead than the reference panel pipeline.


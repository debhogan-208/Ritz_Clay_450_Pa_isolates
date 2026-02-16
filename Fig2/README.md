### Kenneth B. Hoehn
### 2/10/26

Rscript formatDNA.R

Set model parameters using BEAUti with sample_full.fa as input, create sample_full_relaxed.xml

sbatch runAnalyses.sh

sbatch logCombiner.sh

Rscript treeAnalyses.R


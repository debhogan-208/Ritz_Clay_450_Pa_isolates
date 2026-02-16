#!/usr/bin/bash
#SBATCH --ntasks=1 --nodes=1
#SBATCH --nodelist=t01
#SBATCH --account=hoehnlab
#SBATCH --qos=lab_priority
#SBATCH --partition=preempt_t01
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=5gb 
#SBATCH --time=125:00:00
#SBATCH --mail-type=END,FAIL

../beast/bin/beast -overwrite -threads 24 sample_full_relaxed.xml
#../beast/bin/beast -resume -threads 24 sample_full_relaxed.xml


#!/bin/bash -l
#PBS -N @@jobName@@
#PBS -l walltime=@@walltime@@
#PBS -l mem=@@mem@@
#PBS -l ncpus=@@cpus@@
#PBS -P @@rpid@@@@extraDirectives@@

cd $PBS_O_WORKDIR


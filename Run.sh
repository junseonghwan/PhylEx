#!/bin/bash

#SBATCH -A snic2020-5-280
#SBATCH -n 1
#SBATCH -t 48:00:00
#SBATCH --mem 8G
#SBATCH -J scRNA

./../run -c $1

#!/bin/sh
#PBS -lwalltime=01:00:00
#PBS -lselect=1:ncpus=1:mem=4gb

module load anaconda3/personal
source activate py39
conda install netCDF4
pip install pyrealm


echo "running the isoprene model for month ${INDEX}/143"

python isoprene_model_HPC.py $INDEX

echo "done!"

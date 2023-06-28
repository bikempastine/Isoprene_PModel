#!/bin/sh

module load cdo/1.9.0

echo "start submission"

for index in $(seq 110 116); do
	echo "${index}/143"
	qsub -N iso_job_${index} -o $HOME/iso_model_messages/${index}.out -e $HOME/iso_model_messages/${index}.err -v "INDEX=${index}" $HOME/run_iso_model_at_x.sh
done

echo "done!"

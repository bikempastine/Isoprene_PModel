#!/bin/sh

# List of folders to process
folders=("2005" "2006" "2007" "2008" "2009" "2010" "2011" "2012" "2013" "2014" "2015" "2016")

echo "Submitting jobs for years: ${folders[@]} "


# SWdown
echo "Start submssion for SWdown"​

for folder in "${folders[@]}"; do
    qsub -N job_SW_${folder} -o $HOME/averager_messages/messages_SW/${folder}.out -e $HOME/averager_messages/messages_SW/${folder}.err -v "YEAR=${folder}" $HOME/averager_scripts/SWdown_averager.sh 
done


# Tair
echo "Start submssion for Tair"​

for folder in "${folders[@]}"; do
    qsub -N job_Tair_${folder} -o $HOME/averager_messages/messages_Tair/${folder}.out -e $HOME/averager_messages/messages_Tair/${folder}.err -v "YEAR=${folder}" $HOME/averager_scripts/Tair_averager.sh
done

# FPAR
echo "Start submssion for FPAR"​

for folder in "${folders[@]}"; do
    qsub -N job_fpar_${folder} -o $HOME/averager_messages/messages_fpar/${folder}.out -e $HOME/averager_messages/messages_fpar/${folder}.err -v "YEAR=${folder}" $HOME/averager_scripts/fpar_averager.sh 
done

#VAP
echo "Averaging VAP"
cdo -b F32 -monmean -selvar,vap $HOME/average_vap_full/vap/cru_ts4.04.2001.2010.vap.dat.nc $HOME/average_vap_full/monthly/file1_vc.nc
cdo -b F32 -monmean -selvar,vap $HOME/average_vap_full/vap/cru_ts4.04.2011.2019.vap.dat.nc $HOME/average_vap_full/monthly/file2_vc.nc

echo "Script complete"​



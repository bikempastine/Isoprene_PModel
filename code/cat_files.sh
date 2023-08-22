#!/bin/sh

module load cdo/1.9.0
module load anaconda3/personal

rm $HOME/full_averaged/*

echo "starting concatenation"

echo "Tair"
cdo cat $HOME/average_Tair_full/monthly/* $HOME/full_averaged/Tair_monthly_K.nc
cdo addc,-273.15 $HOME/full_averaged/Tair_monthly_K.nc $HOME/full_averaged/Tair_monthly.nc
rm $HOME/full_averaged/Tair_monthly_K.nc

echo "SWdown"
cdo cat $HOME/average_swdown_full/monthly/* $HOME/full_averaged/swdown_monthly.nc

echo "FPAR"
cdo cat $HOME/average_fpar_full/monthly/* $HOME/full_averaged/fpar_monthly.nc

echo "VAP"
rm $HOME/average_vap_full/sandbox/*
cdo cat $HOME/average_vap_full/monthly/* $HOME/average_vap_full/sandbox/vap_monthly_all.nc
cdo selyear,2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016 $HOME/average_vap_full/sandbox/vap_monthly_all.nc $HOME/full_averaged/vap_monthly.nc

echo "done concatenation"

echo "merge the files"
cdo merge $HOME/full_averaged/*_monthly.nc $HOME/full_averaged/Pmodel_4_of_6_data.nc
echo "merged"

echo "remove intermediate files"
rm $HOME/full_averaged/*_monthly.nc

echo "Adding the two other two variables"
python < $HOME/full_add_CO2_elev.py

echo "Done! Saved to: $HOME/full_averaged/Pmodel_data.nc"
rm $HOME/full_averaged/Pmodel_4_of_6_data.nc
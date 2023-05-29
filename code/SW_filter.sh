#!/bin/bash

# Define the input directories for SWdown and Tair files
SW_dir="../data/SWdown_2005"

function daytime_month_averager{
}
# Loop over all the SWdown files in the directory
for x in "$SW_dir"/*.nc
do
    # get the SWdown file name from above
    SWdown=${x}
    echo $SWdown

    #replace the word SWdown with Tair to get the corresponding Tair file
    Tair=$(echo "$SWdown" | sed 's/SWdown/Tair/g')
    echo $Tair

    # Define the output file name
    out_file="../temporary/filtered_mean_$(basename "$Tair")"
    echo $out_file

    # if $SW -gtc(greater than constant)0 (assumed to be daytime)
    # -timmean: get the mean across time of $Tair_file
    # save to the results
    cdo -timmean -ifthen -gtc,0 $SW $Tair $out_file

done

# concatenate all of the monthly mean files using CDO and save the 
# result to a single file in the results folder
cdo cat ../temporary/mean*.nc ../results/$1.nc

# remove all the temporary, monthly mean files
rm ../temporary/filtered_mean*.nc

# display info about the results
cdo -info ../results/$1.nc


# look at data
cdo ntime $out_file
cdo griddes $SW
cdo -info $out_file
cdo showname $out_file



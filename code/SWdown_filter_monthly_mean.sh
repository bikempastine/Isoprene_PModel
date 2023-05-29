
function month_averager {


    ### Loop over all netCDF files in the folder specified in the data folder
    for x in ../data/$1/*.nc; 
    do
        # get the file name from the path given 
        SWdown=${x}

        # replace the word SWdown with Tair to get the corresponding Tair file
        Tair=$(echo "$SWdown" | sed 's/SWdown/Tair/g')
        
        # create a file name for the output file in the temporary folder
        out_file="../temporary/filtered_mean_$(basename "$Tair")"

        # get the mean over time when SWdown is > 0 (meaning it is daytime)
        cdo -timmean -ifthen -gtc,0 $SWdown $Tair $out_file #-timemean = get the mean over time of Tair; -gtc = greater than constant


    done

    # concatenate all of the monthly mean files using CDO and save the 
    # result to a single file in the results folder
    cdo cat ../temporary/filtered_mean*.nc ../results/Tair_filtered_2005.nc

    # remove all the temporary, monthly mean files
    rm ../temporary/filtered_mean*.nc

    # display info about the results
    cdo -info ../results/Tair_filtered_2005.nc

}


# Run the function on the folder "Tair_2005"
month_averager "../data/SWdown_2005"

#convert from K to C
# cdo -addc,-273.15 ../results/Tair_2005.nc ../results/Tair_2005_C.nc

# # display info about the Celcius result
# cdo -info ../results/Tair_2005_C.nc

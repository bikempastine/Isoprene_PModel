#!/bin/bash

# This function takes a single argument, which is the name of the folder 
# containing the netCDF files to be averaged, and calculates the monthly 
# mean for each file in the folder. It then concatenates the resulting 
# monthly mean files and saves the output to a single file in the 
# specified directory.


function month_averager {


    ### Loop over all netCDF files in the folder specified in the data folder
    for x in ../data/$1/*.nc; 
    do
        # get the file name from the patch given 
        data=${x}

        # create a file name for the output file in the temporary folder
        out_file="../temporary/mean_"$(basename $data)

        # Calculate the monthly mean of the file using CDO
        cdo -b F64 -ymonmean $data $out_file

        ### Perform a check to make sure each monthly file returns only one mean 

        # get the number of time steps in the output file
        ntime=$(cdo ntime $out_file)
        # remove the characters from the string return after 
        # the first non-numeric character
        ntime=${ntime%%[^0-9]*} 
        # convert the string result to an integer
        ntime=$((ntime))  

        # check if the output file contains 1 value and give feedback to the user
        if [ "$ntime" -eq 1 ] #if the number of time steps is 1
        then
            echo "sucsessful averaging, monthly mean found!"
        else
            echo "Error: More than 1 time point returned. Loop ended"
            break  #break the for loop
        fi

    done

    # concatenate all of the monthly mean files using CDO and save the 
    # result to a single file in the results folder
    cdo cat ../temporary/mean*.nc ../results/$1.nc

    # remove all the temporary, monthly mean files
    rm ../temporary/mean*.nc

    # display info about the results
    cdo -info ../results/$1.nc

}


# # Run the function on the folder "Tair_2005"
# month_averager "Tair_2005"

# #convert from K to C
# cdo -addc,-273.15 ../results/Tair_2005.nc ../results/Tair_2005_C.nc

# # display info about the Celcius result
# cdo -info ../results/Tair_2005_C.nc

# cdo -info ../results/Tair_2005.nc



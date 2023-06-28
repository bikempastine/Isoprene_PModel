#!/bin/sh
#PBS -lwalltime=03:00:00
#PBS -lselect=1:ncpus=1:mem=1gb

module load cdo/1.9.0

# Print the folder name
echo "Processing SW folder: $YEAR"

### Loop over all netCDF files in the folder specified in the data folder
for x in $HOME/average_swdown_full/$YEAR/*.nc; 
    do
        # get the file name from the patch given 
        data=${x}
        echo "starting ${x}"

        # create a file name for the output file in the temporary folder
        out_file="$HOME/average_swdown_full/sandbox/${YEAR}_mean_"$(basename $data)

        # Calculate the monthly mean of the file using CDO
        #cdo -b F32 -ymonmean $data $out_file
	cdo timmean -setvrange,0,100000 $data $out_file

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

echo "Done averaging SWDown for $YEAR"

# concatenate all of the monthly mean files using CDO and save the 
# result to a single file in the results folder
cdo cat $HOME/average_swdown_full/sandbox/${YEAR}_mean_*.nc $HOME/average_swdown_full/monthly/SWdown_${YEAR}_output.nc

# remove all the temporary, monthly mean files
rm  $HOME/average_swdown_full/sandbox/${YEAR}_mean*.nc

echo "Done!"

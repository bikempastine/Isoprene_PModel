
### Contents
| Name  |      Objective      |  Input |  Output |
|----------|:-------------:|------:|------:|
| clean_OMI_MODEL_data.py |  Clean up the OMI and Model isoprene data. For example, flip the omi data, filter out extreme values or fill values etc. | '../data/model.nc' , '../data/model_monthly_average.nc', '../data/OMI.nc', '../data/OMI_monthly_average.nc' | `../data/model_clean.nc`,`../data/model_monthly_average_clean.nc`, `../data/OMI_clean.nc`, `../data/OMI_monthly_average_clean.nc` |
| factor_isoprene_scatter_plot.py | Calculate F using the full dataset and plot it in scatter plots against isoprene levels | '../data/OMI_clean.nc' , '../data/model_clean.nc' | two scatter plots |
| monthly_avg_omi_model_scatter.py | plot the monthly averaged omi and model data against each other wih a line with the slope of the averaged F | `../data/model_monthly_average_clean.nc`, `../data/OMI_monthly_average_clean.nc` | one scatter plot |

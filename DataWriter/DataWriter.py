import netCDF4
import numpy as np
from DataExtraction.ExtractCMAQ import CMAQGridInfo
import shutil

def combineFusionResults(data_fusion_results):
    """

    :param data_fusion_results: results dictionary from data fusion process
    :return: a tensor contains all years data fusion results (time dimension, 1, x_dimension, y_dimension)
    """
    year_list = []
    for key in data_fusion_results.keys():
        year_list.append(key)
    year_list.sort()
    all_year_results = []
    for current_year in year_list:
        current_fusion_res = data_fusion_results[current_year]
        # add back layer dimension
        current_fusion_res = current_fusion_res[:, np.newaxis, :, :]
        all_year_results.append(current_fusion_res)
    all_year_results = np.concatenate(all_year_results, axis=0)
    return all_year_results


def writeToNetCDF(data_fusion_results, pollutant_name, input_file_name, output_file_name):
    combined_results = combineFusionResults(data_fusion_results)
    shutil.copyfile(input_file_name, output_file_name)
    dset = netCDF4.Dataset(output_file_name, 'r+')
    dset[pollutant_name][:] = combined_results
    dset.close()


def writeToNetCDFwithGeo(data_fusion_results, pollutant_name, input_file_name, output_file_name):
    combined_results = combineFusionResults(data_fusion_results)
    # Generate spatial information
    grid_info = CMAQGridInfo(input_file_name)
    with netCDF4.Dataset(input_file_name) as src, netCDF4.Dataset(output_file_name, "w") as dst:
        # copy global attributes all at once via dictionary
        dst.setncatts(src.__dict__)
        # copy dimensions
        for name, dimension in src.dimensions.items():
            dst.createDimension(name, len(dimension))
        # copy all file data except for the excluded
        for name, variable in src.variables.items():
            if name == pollutant_name:
                x = dst.createVariable(name, variable.datatype, variable.dimensions)
                dst[name][:] = combined_results
                # copy variable attributes all at once via dictionary
                dst[name].setncatts(src[name].__dict__)

                # Lat
                x = dst.createVariable("Lat", variable.datatype, ('ROW', 'COL'))
                dst["Lat"][:] = grid_info["Lat"]
                description = {
                    "long_name": "Latitude",
                    "units": "degree",
                    "var_desc": "Latitude coordinate for CMAQ grids"
                }
                dst["Lat"].setncatts(description)

                # Lon
                x = dst.createVariable("Lon", variable.datatype, ('ROW', 'COL'))
                dst["Lon"][:] = grid_info["Lon"]
                description = {
                    "long_name": "Longitude",
                    "units": "degree",
                    "var_desc": "Longitude coordinate for CMAQ grids"
                }
                dst["Lon"].setncatts(description)

            if name == "TFLAG":
                x = dst.createVariable(name, variable.datatype, variable.dimensions)
                dst[name][:] = src[name][:]
                # copy variable attributes all at once via dictionary
                dst[name].setncatts(src[name].__dict__)
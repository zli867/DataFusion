import netCDF4
import numpy as np


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
            if name == "TFLAG":
                x = dst.createVariable(name, variable.datatype, variable.dimensions)
                dst[name][:] = src[name][:]
                # copy variable attributes all at once via dictionary
                dst[name].setncatts(src[name].__dict__)
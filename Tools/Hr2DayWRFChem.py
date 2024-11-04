import sys
import os
from pathlib import Path
project_path = Path(__file__).parent.parent
sys.path.append(os.path.abspath(project_path))
import numpy as np
import netCDF4 as nc
from datetime import datetime, timedelta
from util.DailyMetrics import HR_MAX, AVG, MDA8
from DataExtraction.ExtractCMAQ import WRFGridInfo
from util.GeoProcess import UTC_offset


start_process_time = datetime(2017, 1, 1)
end_process_time = datetime(2017, 1, 1)
file_dir = "/Volumes/DataStorage/WRFChem"
output_file = os.path.join("/Volumes/DataStorage/WRFChem/", "WRFChem_Combine_2017.nc")

# extract static info
first_chem_file = os.path.join(file_dir, "wrfout_d01_%s_00:00:00" % start_process_time.strftime("%Y-%m-%d"))
static_info = WRFGridInfo(first_chem_file)
static_info["UTC_Offset"] = UTC_offset(static_info["Lon"], static_info["Lat"])
utc_offset = static_info["UTC_Offset"]
lat = static_info["Lat"]
lon = static_info["Lon"]
m, n = utc_offset.shape
time_array = []
variables_info = {
    "NO2_AVG": {"vars": "no2", "func": AVG, "values": []},
    "O3_MDA8": {"vars": "o3", "func": MDA8, "values": []},
    "PM25_AVG": {"vars": "PM2_5_DRY", "func": AVG, "values": []}
    }
current_time = start_process_time

while current_time <= end_process_time:
    print("processing %s" % current_time.strftime("%Y%m%d"))
    cur_chem = os.path.join(file_dir, "wrfout_d01_%s_00:00:00" % current_time.strftime("%Y-%m-%d"))
    nxt_chem = os.path.join(file_dir, "wrfout_d01_%s_00:00:00" % (current_time + timedelta(days=1)).strftime("%Y-%m-%d"))
    if os.path.exists(cur_chem) and os.path.exists(nxt_chem):
        cur_chem_ds = nc.Dataset(cur_chem)
        nxt_chem_ds = nc.Dataset(nxt_chem)
        time_array.append(cur_chem_ds["Times"][0, :])
        for variable_name in variables_info.keys():
            vars = []
            input_var_name = variables_info[variable_name]["vars"]
            cur_chem_var_value = cur_chem_ds[input_var_name][:][:, 0, :, :]
            nxt_chem_var_value = nxt_chem_ds[input_var_name][:][:, 0, :, :]
            chem_var_value = np.concatenate((cur_chem_var_value, nxt_chem_var_value), axis=0)
            input_var_value = np.zeros((24, m, n))
            # process to extract local date variable values
            for i in range(0, m):
                for j in range(0, n):
                    # only for US
                    start_idx = int(utc_offset[i, j] * -1)
                    input_var_value[:, i, j] = chem_var_value[start_idx: start_idx + 24, i, j]
            vars.append(input_var_value)
            # function processed
            daily_value = variables_info[variable_name]["func"](*vars)
            # save data to dict
            variables_info[variable_name]["values"].append(daily_value[np.newaxis, np.newaxis, :, :])
    else:
        print("Missing Chem files: %s" % (current_time.strftime("%Y%m%d")))

    current_time += timedelta(days=1)

# processed data and wrap to csv
chem_vars = {}
for variable_name in variables_info.keys():
    variables_info[variable_name]["values"] = np.concatenate(variables_info[variable_name]["values"], axis=0)
    chem_vars[variables_info[variable_name]["vars"]] = variable_name

# print("Hello")
with nc.Dataset(first_chem_file) as src, nc.Dataset(output_file, "w") as dst:
    dst.setncatts(src.__dict__)
    # copy dimensions
    for name, dimension in src.dimensions.items():
        if name == "Time":
            dst.createDimension(name, len(time_array))
        if name == "bottom_top":
             dst.createDimension(name, 1)
        elif name in ["DateStrLen", "west_east", "south_north"]:
            dst.createDimension(name, len(dimension))
        # else:
        #     print(name)
        #     dst.createDimension(name, len(dimension))
    # copy all file data except for the excluded
    for name, variable in src.variables.items():
        if name == "Times":
            var_tmp = dst.createVariable(name, variable.datatype, variable.dimensions)
            dst[name][:] = time_array
            dst[name].setncatts(src[name].__dict__)
        elif name in ["XLAT", "XLONG"]:
            var_tmp = dst.createVariable(name, variable.datatype, variable.dimensions)
            src_val = np.squeeze(src[name][0, :, :])
            dst[name][:] = np.repeat(src_val[np.newaxis, :, :], len(time_array), axis=0)
            dst[name].setncatts(src[name].__dict__)
        elif name in chem_vars.keys():
            var_tmp = dst.createVariable(chem_vars[name], variable.datatype, variable.dimensions)
            attrs = src[name].__dict__
            dst[chem_vars[name]].setncatts(src[name].__dict__)
            dst[chem_vars[name]][:] = variables_info[chem_vars[name]]["values"]
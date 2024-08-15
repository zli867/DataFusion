import sys
import os
from pathlib import Path
project_path = Path(__file__).parent.parent
sys.path.append(os.path.abspath(project_path))
import numpy as np
import netCDF4 as nc
from datetime import datetime, timedelta
from util.DailyMetrics import HR_MAX, AVG, MDA8
from DataExtraction.ExtractCMAQ import GCGridInfo
from util.GeoProcess import UTC_offset


start_process_time = datetime(2015, 1, 1)
end_process_time = datetime(2015, 12, 30)
file_dir = "/Volumes/Expansion/DataFusionData/GC/2015"
output_file = os.path.join("/Volumes/Expansion/DataFusionData/GC/combined", "GC_Combine_2015.nc4")

# extract static info
first_gc_file = os.path.join(file_dir, "GEOSChem.SpeciesConc.%s_0000z.nc4" % start_process_time.strftime("%Y%m%d"))
static_info = GCGridInfo(first_gc_file)
static_info["UTC_Offset"] = UTC_offset(static_info["Lon"], static_info["Lat"])
utc_offset = static_info["UTC_Offset"]
lat = static_info["Lat"]
lon = static_info["Lon"]
m, n = utc_offset.shape
time_array = []
variables_info = {"NO2_AVG": {"vars": "SpeciesConc_NO2", "func": AVG, "values": []},
                  "O3_MDA8": {"vars": "SpeciesConc_O3", "func": MDA8, "values": []}}
current_time = start_process_time

while current_time <= end_process_time:
    print("processing %s" % current_time.strftime("%Y%m%d"))
    cur_gc = os.path.join(file_dir, "GEOSChem.SpeciesConc.%s_0000z.nc4" % current_time.strftime("%Y%m%d"))
    nxt_gc = os.path.join(file_dir, "GEOSChem.SpeciesConc.%s_0000z.nc4" % (current_time + timedelta(days=1)).strftime("%Y%m%d"))
    if os.path.exists(cur_gc) and os.path.exists(nxt_gc):
        cur_gc_ds = nc.Dataset(cur_gc)
        nxt_gc_ds = nc.Dataset(nxt_gc)
        time_array.append((current_time - start_process_time).total_seconds() // 60)
        for variable_name in variables_info.keys():
            vars = []
            input_var_name = variables_info[variable_name]["vars"]
            cur_gc_var_value = cur_gc_ds[input_var_name][:][:, 0, :, :]
            nxt_gc_var_value = nxt_gc_ds[input_var_name][:][:, 0, :, :]
            gc_var_value = np.concatenate((cur_gc_var_value, nxt_gc_var_value), axis=0)
            input_var_value = np.zeros((24, m, n))
            # process to extract local date variable values
            for i in range(0, m):
                for j in range(0, n):
                    # only for US
                    start_idx = int(utc_offset[i, j] * -1)
                    input_var_value[:, i, j] = gc_var_value[start_idx: start_idx + 24, i, j]
            vars.append(input_var_value)
            # function processed
            daily_value = variables_info[variable_name]["func"](*vars)
            # save data to dict
            variables_info[variable_name]["values"].append(daily_value[np.newaxis, np.newaxis, :, :])
    else:
        print("Missing GC files: %s" % (current_time.strftime("%Y%m%d")))

    current_time += timedelta(days=1)

# processed data and wrap to csv
gc_vars = {}
for variable_name in variables_info.keys():
    variables_info[variable_name]["values"] = np.concatenate(variables_info[variable_name]["values"], axis=0)
    gc_vars[variables_info[variable_name]["vars"]] = variable_name

with nc.Dataset(first_gc_file) as src, nc.Dataset(output_file, "w") as dst:
    dst.setncatts(src.__dict__)
    # copy dimensions
    for name, dimension in src.dimensions.items():
        if name == "time":
            dst.createDimension(name, len(time_array))
        elif name == "lev":
            dst.createDimension(name, 1)
        elif name == "ilev":
            dst.createDimension(name, 2)
        else:
            dst.createDimension(name, len(dimension))

    # copy all file data except for the excluded
    for name, variable in src.variables.items():
        if name == "time":
            var_tmp = dst.createVariable(name, variable.datatype, variable.dimensions)
            dst[name][:] = time_array
            dst[name].setncatts(src[name].__dict__)
        elif name in ["lev", "hyam", "hybm"]:
            var_tmp = dst.createVariable(name, variable.datatype, variable.dimensions)
            dst[name][:] = src[name][0:1]
            dst[name].setncatts(src[name].__dict__)
        elif name in ["ilev", "hyai", "hybi"]:
            var_tmp = dst.createVariable(name, variable.datatype, variable.dimensions)
            dst[name][:] = src[name][0:2]
            dst[name].setncatts(src[name].__dict__)
        elif name in ["lat", "lon", "Area", "P0"]:
            var_tmp = dst.createVariable(name, variable.datatype, variable.dimensions)
            dst[name][:] = src[name][:]
            dst[name].setncatts(src[name].__dict__)
        elif name in gc_vars.keys():
            var_tmp = dst.createVariable(gc_vars[name], variable.datatype, variable.dimensions)
            attrs = src[name].__dict__
            dst[gc_vars[name]].setncatts(src[name].__dict__)
            dst[gc_vars[name]][:] = variables_info[gc_vars[name]]["values"]
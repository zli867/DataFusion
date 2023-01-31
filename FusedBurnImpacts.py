import netCDF4
from DataExtraction.ExtractCMAQ import CMAQGridInfo
import numpy as np
from util.FormatConverter import timeArrayToTFLAG

CMAQ_fire_file = "/Users/zongrunli/Desktop/work_dir/CCTM.ACONC.combined.FIRE.PM25.hires4.2015_2020.nc"
CMAQ_no_fire_file = "/Users/zongrunli/Desktop/work_dir/CCTM.ACONC.combined.NOFIRE.PM25.hires4.2015_2020.nc"
CMAQ_fused_file = "/Users/zongrunli/Desktop/work_dir/CCTM.ACONC.combined.FIRE.PM25.hires4.2015_2020_fused.nc"
output_filename = "/Users/zongrunli/Desktop/work_dir/CCTM.ACONC.combined.FIRE.PM25.hires4.2015_2020_fused_impact.nc"
pollutant_name = "PM25_TOT_AVG"

# read time for each file
fire_info = CMAQGridInfo(CMAQ_fire_file)
no_fire_info = CMAQGridInfo(CMAQ_no_fire_file)
fused_info = CMAQGridInfo(CMAQ_fused_file)

# common time
common_time = list(set(fire_info["time"]) & set(no_fire_info["time"]) & set(fused_info["time"]))
common_time.sort()

# valid pollutant index for each field
fire_valid_idx = []
no_fire_valid_idx = []
fused_valid_idx = []
for current_time in common_time:
    fire_valid_idx.append(fire_info["time"].index(current_time))
    no_fire_valid_idx.append(no_fire_info["time"].index(current_time))
    fused_valid_idx.append(fused_info["time"].index(current_time))

# calculate burn impact
fire_pollutant = netCDF4.Dataset(CMAQ_fire_file)[pollutant_name][fire_valid_idx, :, :, :]
no_fire_pollutant = netCDF4.Dataset(CMAQ_no_fire_file)[pollutant_name][no_fire_valid_idx, :, :, :]
fused_pollutant = netCDF4.Dataset(CMAQ_fused_file)[pollutant_name][fused_valid_idx, :, :, :]
fire_fused_impact = ((fire_pollutant - no_fire_pollutant)/fire_pollutant) * fused_pollutant


# Write NetCDF for Burn Impact
with netCDF4.Dataset(CMAQ_fused_file) as src, netCDF4.Dataset(output_filename, "w") as dst:
    # copy global attributes all at once via dictionary
    dst.setncatts(src.__dict__)
    # copy dimensions
    for name, dimension in src.dimensions.items():
        if name == "TSTEP":
            dst.createDimension(name, len(common_time))
        elif name == "VAR":
            dst.createDimension(name, 1)
        else:
            dst.createDimension(name, len(dimension))
    # copy all file data except for the excluded
    for name, variable in src.variables.items():
        if name == pollutant_name:
            x = dst.createVariable(name, variable.datatype, variable.dimensions)
            # calculate burn impact
            dst[name][:] = fire_fused_impact
            # copy variable attributes all at once via dictionary
            dst[name].setncatts(src[name].__dict__)
        if name == "TFLAG":
            x = dst.createVariable(name, variable.datatype, variable.dimensions)
            dst[name][:] = timeArrayToTFLAG(common_time, 1)
            # copy variable attributes all at once via dictionary
            dst[name].setncatts(src[name].__dict__)

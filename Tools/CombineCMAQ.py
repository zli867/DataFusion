import netCDF4
import os

pollutant_names = ["O3_MDA8", "NO2_AVG", "PM25_AVG"]
dir = "/Volumes/Expansion/DataFusionData/CMAQ/test/"
output_file = "/Volumes/Expansion/DataFusionData/CCTM.ACONC.combined.2015.nc"

filenames = os.listdir(dir)
#  sort the filename to make sure the CMAQ data is from previous year to current year
filenames.sort()
filepaths = []
for filename in filenames:
    if ".nc" not in filename or filename[0] == ".":
        continue
    filepaths.append(dir + filename)

combined_dataset = netCDF4.MFDataset(filepaths)

with netCDF4.Dataset(filepaths[0]) as src, netCDF4.Dataset(output_file, "w",
                                                           format='NETCDF3_64BIT_OFFSET',
                                                           diskless=True,
                                                           persist=True) as dst:
    globle_dict = src.__dict__
    dst.setncatts(globle_dict)
    # copy dimensions
    for name, dimension in combined_dataset.dimensions.items():
        dst.createDimension(name, len(dimension))
    # copy all file data except for the excluded
    for name, variable in src.variables.items():
        if name == 'TFLAG' or name in pollutant_names:
            var_tmp = dst.createVariable(name, variable.datatype, variable.dimensions)
            dst[name][:] = combined_dataset[name][:]
            dst[name].setncatts(src[name].__dict__)
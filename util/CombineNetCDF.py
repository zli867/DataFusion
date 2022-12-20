import netCDF4
import os

dir = "/Volumes/SERDP Modeling/CDC_CMAQ/2016_2020/FIRE/4km/"
filenames = os.listdir(dir)
#  sort the filename to make sure the CMAQ data is from previous year to current year
filenames.sort()
output_file = "../data/CCTM.ACONC.combined.FIRE.hires4.2016_2020.nc"
filepaths = []
for filename in filenames:
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
        if name == 'TFLAG' or name == 'PM25_TOT_AVG':
            var_tmp = dst.createVariable(name, variable.datatype, variable.dimensions)
            dst[name][:] = combined_dataset[name][:]
            dst[name].setncatts(src[name].__dict__)
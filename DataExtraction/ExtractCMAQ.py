import netCDF4 as nc
import numpy as np
import pyproj
from datetime import datetime, timedelta


def CMAQGridInfo(CMAQ_file):
    cmaq_ds = nc.Dataset(CMAQ_file)
    time_data = cmaq_ds['TFLAG'][:]
    lat_1 = cmaq_ds.getncattr('P_ALP')
    lat_2 = cmaq_ds.getncattr('P_BET')
    lat_0 = cmaq_ds.getncattr('YCENT')
    lon_0 = cmaq_ds.getncattr('XCENT')
    crs = pyproj.Proj("+proj=lcc +a=6370000.0 +b=6370000.0 +lat_1=" + str(lat_1)
                      + " +lat_2=" + str(lat_2) + " +lat_0=" + str(lat_0) +
                      " +lon_0=" + str(lon_0))
    xcell = cmaq_ds.getncattr('XCELL')
    ycell = cmaq_ds.getncattr('YCELL')
    xorig = cmaq_ds.getncattr('XORIG')
    yorig = cmaq_ds.getncattr('YORIG')

    ncols = cmaq_ds.getncattr('NCOLS')
    nrows = cmaq_ds.getncattr('NROWS')

    # > for X, Y cell centers
    x_center_range = np.linspace(xorig + xcell / 2, (xorig + xcell / 2) + xcell * (ncols - 1), ncols)
    y_center_range = np.linspace(yorig + ycell / 2, (yorig + ycell / 2) + ycell * (nrows - 1), nrows)

    Xcenters, Ycenters = np.meshgrid(x_center_range, y_center_range)

    # > for X, Y cell boundaries (i.e., cell corners)
    x_bound_range = np.linspace(xorig, xorig + xcell * ncols, ncols + 1)
    y_bound_range = np.linspace(yorig, yorig + ycell * nrows, nrows + 1)

    Xbounds, Ybounds = np.meshgrid(x_bound_range, y_bound_range)

    x_max = np.max(Xbounds)
    x_min = np.min(Xbounds)
    y_max = np.max(Ybounds)
    y_min = np.min(Ybounds)

    cmaq_time_array = []
    for i in range(0, time_data.shape[0]):
        time_data_tmp = time_data[i, 0, :]
        time_str = str(time_data_tmp[0]) + str(time_data_tmp[1]).rjust(6, '0')
        parsed = datetime.strptime(time_str, '%Y%j%H%M%S')
        cmaq_time_array.append(parsed)

    lon_ctr, lat_ctr = crs(Xcenters, Ycenters, inverse=True)
    lat_min = np.min(lat_ctr)
    lat_max = np.max(lat_ctr)
    lon_min = np.min(lon_ctr)
    lon_max = np.max(lon_ctr)
    res_dict = {"crs": crs, "X": Xcenters, "Y": Ycenters, "X_bdry": [x_min, x_max], "Y_bdry": [y_min, y_max],
                "time": cmaq_time_array, "Lat": lat_ctr, "Lon": lon_ctr, "Lat_bdry": [lat_min, lat_max],
                "Lon_bdry": [lon_min, lon_max], "XCELL": xcell, "YCELL": ycell}
    return res_dict


def searchTimeIndex(time, year):
    """

    :param time: a list of datetime
    :param year: int number, the year need to be searched
    :return: (start_idx, end_idx)
             start_idx: the index of the first day of the year in the list
             end_idx: the index of the end day of the year in the list
    """
    start_idx = None
    end_idx = None
    start_time = datetime(year, 1, 1)
    end_time = datetime(year, 12, 31)
    current_time = start_time
    # find start_idx
    while current_time <= end_time:
        if current_time in time:
            start_idx = time.index(current_time)
            break
        else:
            current_time += timedelta(days=1)
    # find end_idx
    current_time = end_time
    while current_time >= start_time:
        if current_time in time:
            end_idx = time.index(current_time)
            break
        else:
            current_time -= timedelta(days=1)
    return start_idx, end_idx


def extractCMAQ(pollutant, filename):
    """
    :param pollutant: the variable name need to do the data fusion
    :param filename: CMAQ file name
    :return: a list of a dictionary. The dictionary contains daily, yearly CMAQ data and unit of the variable
    """
    file_data = nc.Dataset(filename)
    polltant_conc = file_data.variables[pollutant]
    conc_surface_layer = polltant_conc[:][:, 0, :, :]
    grid_info = CMAQGridInfo(filename)
    year_set = []
    for i in range(0, len(grid_info["time"])):
        year_set.append(grid_info["time"][i].year)
    year_set = list(set(year_set))
    # CMAQ combined dict
    cmaq_combined_dict = {
        "Time": grid_info["time"],
        "Daily": conc_surface_layer,
        "unit": polltant_conc.units
    }
    # CMAQ yearly dict
    cmaq_dict = {}
    for year in year_set:
        start_time_idx, end_time_idx = searchTimeIndex(grid_info["time"], year)
        current_year_conc = conc_surface_layer[start_time_idx: end_time_idx + 1, :, :]
        current_time_list = grid_info["time"][start_time_idx: end_time_idx + 1]
        cmaq_yearly_dict = {
            "Time": current_time_list,
            "Daily": current_year_conc,
            "Yearly": np.mean(current_year_conc, axis=0),
            "unit": polltant_conc.units,
            "year": year
        }
        cmaq_dict[year] = cmaq_yearly_dict

    return [cmaq_combined_dict, cmaq_dict]

# #
# filename = "/Users/zongrunli/Desktop/DataFusion_v1.0/util/out.nc"
# print(extractCMAQ("PM25_TOT_AVG", filename)[0]["Time"])
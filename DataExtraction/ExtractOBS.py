import pandas as pd
import numpy as np
from datetime import datetime
import sys


def columnCheck(pollutant, dataframe):
    """

    :param dataframe: dataframe read from observation files
    :return: True: the data is valid, False: the data is invalid
    """
    valid = True
    required_fields = ["Time", "siteCode", "Latitude", "Longitude", pollutant]
    for required_field in required_fields:
        valid = valid and (required_field in dataframe.columns)
    # check datetime format
    return valid


def extractOBS(pollutant, filename, grid_info):
    """
    We firstly remove the observations out of the CMAQ boundary, remove the date we do not have any CMAQ simulations (it
    means date range of observation <= date range of CMAQ), the data points which have negative pollutant concentration.
    Then, the concentration is rearraged into the format: (timeSeries * site). Notice that, if we calculate
    np.nanmean(conc, axis=0) or np.nanmean(conc, axis=1), there will be no nan values.
    :param pollutant: pollutant name in the csv file
    :param filename: csv filename which includes the observation information
    :param grid_info:
    :return: a list of a dictionary. The dictionary contains daily, yearly observation data
    """
    dateparse = lambda x: datetime.strptime(x, '%Y-%m-%d')
    df = pd.read_csv(filename, dtype={"Time": str, "siteCode": str, "Latitude": float,
                                      "Longitude": float, pollutant: float}, parse_dates=['Time'], date_parser=dateparse)

    column_check = columnCheck(pollutant, df)
    if not column_check:
        sys.exit("The observation data format is not correct, check the columns.")

    # Remove the negative pollutant data (really simple quality check, quality check should be
    # user's responsibility)
    df[df[pollutant] < 0] = np.nan
    df = df.dropna()
    df = df.reset_index(drop=True)

    crs = grid_info["crs"]
    cmaq_time = grid_info["time"]
    # add X and Y to the dataframe
    loc_x, loc_y = crs(df["Longitude"].to_numpy(), df["Latitude"].to_numpy(), inverse=False)
    df["X"] = loc_x
    df["Y"] = loc_y
    # filtered by the CMAQ boundary
    x_bdry = grid_info["X_bdry"]
    y_bdry = grid_info["Y_bdry"]
    df = df[(df["X"] >= x_bdry[0]) & (df["X"] <= x_bdry[1]) & (df["Y"] >= y_bdry[0]) & (df["Y"] <= y_bdry[1])]
    df = df.reset_index(drop=True)

    # filtered by the CMAQ time, it means date range of observation <= date range of CMAQ
    df = df[df['Time'].isin(cmaq_time)]
    date_series = list(set(list(df["Time"])))
    site_code_set = list(set(df["siteCode"]))
    date_series.sort()
    site_code_set.sort()
    year_set = []
    for date_data in date_series:
        year_set.append(date_data.year)
    year_set = list(set(year_set))

    lat = np.empty((len(date_series), len(site_code_set)))
    lat[:] = np.NaN
    lon = np.empty((len(date_series), len(site_code_set)))
    lon[:] = np.NaN
    X = np.empty((len(date_series), len(site_code_set)))
    X[:] = np.NaN
    Y = np.empty((len(date_series), len(site_code_set)))
    Y[:] = np.NaN
    conc = np.empty((len(date_series), len(site_code_set)))
    conc[:] = np.NaN

    for idx, row in df.iterrows():
        spatial_index = site_code_set.index(row["siteCode"])
        time_index = date_series.index(row["Time"])
        lat[time_index, spatial_index] = row["Latitude"]
        lon[time_index, spatial_index] = row["Longitude"]
        X[time_index, spatial_index] = row["X"]
        Y[time_index, spatial_index] = row["Y"]
        conc[time_index, spatial_index] = row[pollutant]

    # Combined obs dict
    obs_combined_dict = {
        "siteCode": site_code_set,
        "dateSeries": date_series,
        "Lat": np.nanmean(lat, axis=0),
        "Lon": np.nanmean(lon, axis=0),
        "X": np.nanmean(X, axis=0),
        "Y": np.nanmean(Y, axis=0),
        "Conc": conc
    }

    # Yearly obs dict
    obs_dict = {}
    for year in year_set:
        current_year_df = df[(df["Time"] >= datetime(year, 1, 1)) & (df["Time"] <= datetime(year, 12, 31))]
        current_date_series = list(set(list(current_year_df["Time"])))
        current_site_code = list(set(current_year_df["siteCode"]))
        current_date_series.sort()
        current_site_code.sort()
        # timeSeries * spatial
        # Conc[0] first day of all sites
        lat = np.empty((len(current_date_series), len(current_site_code)))
        lat[:] = np.NaN
        lon = np.empty((len(current_date_series), len(current_site_code)))
        lon[:] = np.NaN
        X = np.empty((len(current_date_series), len(current_site_code)))
        X[:] = np.NaN
        Y = np.empty((len(current_date_series), len(current_site_code)))
        Y[:] = np.NaN
        conc = np.empty((len(current_date_series), len(current_site_code)))
        conc[:] = np.NaN

        for idx, row in current_year_df.iterrows():
            spatial_index = current_site_code.index(row["siteCode"])
            time_index = current_date_series.index(row["Time"])
            lat[time_index, spatial_index] = row["Latitude"]
            lon[time_index, spatial_index] = row["Longitude"]
            X[time_index, spatial_index] = row["X"]
            Y[time_index, spatial_index] = row["Y"]
            conc[time_index, spatial_index] = row[pollutant]

        dict = {
            "year": year,
            "siteCode": current_site_code,
            "dateSeries": current_date_series,
            "Lat": np.nanmean(lat, axis=0),
            "Lon": np.nanmean(lon, axis=0),
            "X": np.nanmean(X, axis=0),
            "Y": np.nanmean(Y, axis=0),
            "Conc": conc
        }
        obs_dict[year] = dict
    return [obs_combined_dict, obs_dict]


# from ExtractCMAQ import CMAQGridInfo
# grid_info = CMAQGridInfo("/Users/zongrunli/Desktop/DataFusion_v1.0/util/out.nc")
# x, y = extractOBS("PM25", "/Users/zongrunli/Desktop/DataFusion_v1.0/util/test.csv", grid_info)
# print(x)
# print(y)

import pandas as pd
import os

data_dir = "/Volumes/Expansion/DataFusionData/obs/PM25/"
filenames = os.listdir(data_dir)
for filename in filenames:
    if ".csv" not in filename or filename[0] == ".":
        filenames.remove(filename)

dfs = []
for filename in filenames:
    filepath = data_dir + filename
    df = pd.read_csv(filepath)
    print(filename)
    df = df[["State Code", "County Code", "Site Num", "Latitude", "Longitude", "Date Local", "Arithmetic Mean",
                     "Units of Measure"]]
    df["Time"] = df["Date Local"]
    df["PM25"] = df["Arithmetic Mean"]
    # add site code
    siteCode = []
    for idx, row in df.iterrows():
        stateCode_temp = str(row[0])
        countyCode_temp = str(row[1]).rjust(3, '0')
        siteNum_temp = str(row[2]).rjust(4, '0')
        siteCode_temp = "US" + stateCode_temp + countyCode_temp + siteNum_temp
        siteCode.append(siteCode_temp)
    df["siteCode"] = siteCode
    df = df[["Time", "PM25", "siteCode", "Latitude", "Longitude"]]
    dfs.append(df)
dfs = pd.concat(dfs, ignore_index=True)
dfs = dfs.sort_values(by='Time')
dfs.to_csv("/Volumes/Expansion/DataFusionData/obs/combined/PM25_obs_2010_2019.csv", index=False)
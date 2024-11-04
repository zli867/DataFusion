from DataExtraction.ExtractCMAQ import extractCMAQ, CMAQGridInfo, GCGridInfo, WRFGridInfo
from DataExtraction.ExtractOBS import extractOBS
from DataFusion.DataFusionDriver import dataFusion
from DataWriter.DataWriter import combineFusionResults
import netCDF4
import pandas as pd
from sklearn.model_selection import KFold
from Evaluation.StatisticalMetrics import stats_metrics
import os
from datetime import datetime

dateparse = lambda x: datetime.strptime(x, '%Y-%m-%d')

# I want to just save the concentration I evaluated, so I overwrite the method.
def writeToNetCDF(data_fusion_results, pollutant_name, input_file_name, output_file_name):
    combined_results = combineFusionResults(data_fusion_results)
    # Generate spatial information
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
                dst[name].setncatts(src[name].__dict__)
                dst[name][:] = combined_results
            if name in ["Times", "XLAT", "XLONG"]:
                x = dst.createVariable(name, variable.datatype, variable.dimensions)
                dst[name].setncatts(src[name].__dict__)
                dst[name][:] = src[name][:]

# random withhoding
work_dir =  "/Volumes/DataStorage/DataFusionData/validation/WRFChem/O3"
CMAQ_file = "/Volumes/DataStorage/DataFusionData/WRFChem/WRFChem_Combine_2017.nc"
obs_file =  "/Volumes/DataStorage/DataFusionData/WRFChem/O3_obs_2017_CONUS_ppm.csv"
CMAQ_pollutant = "O3_MDA8"
obs_pollutant = "O3"
k_fold = 5

geo = WRFGridInfo(CMAQ_file)
CMAQList = extractCMAQ(CMAQ_pollutant, CMAQ_file, geo)

obs_df = pd.read_csv(obs_file, parse_dates=["Time"], date_parser=dateparse)
obs_df = obs_df[(obs_df["Time"] >= datetime(2017, 1, 1)) & (obs_df["Time"] <= datetime(2017, 12, 31))]
kf = KFold(n_splits=k_fold, shuffle=True, random_state=100)
obs_train_filenames = []
obs_test_filenames = []
for i, (train_index, test_index) in enumerate(kf.split(obs_df)):
    obs_train_filename = os.path.join(work_dir, "train_obs_rfold_" + str(i) + ".csv")
    obs_test_filename = os.path.join(work_dir, "test_obs_rfold_" + str(i) + ".csv")
    train_df = obs_df.iloc[train_index]
    test_df = obs_df.iloc[test_index]
    train_df = train_df.sort_values(by='Time')
    test_df = test_df.sort_values(by='Time')
    train_df.to_csv(obs_train_filename, index=False)
    test_df.to_csv(obs_test_filename, index=False)
    obs_train_filenames.append(obs_train_filename)
    obs_test_filenames.append(obs_test_filename)

for i in range(0, k_fold):
    data_fused_filename = os.path.join(work_dir, "WRFChem_Combine_2017_o3_fused_fold_" + str(i) + ".nc4")
    train_obs_list = extractOBS(obs_pollutant, obs_train_filenames[i], geo)
    datafusion_results = dataFusion(CMAQList, train_obs_list, geo, 'default')
    writeToNetCDF(datafusion_results, CMAQ_pollutant, CMAQ_file, data_fused_filename)


# site wide withholding
work_dir =  "/Volumes/DataStorage/DataFusionData/validation/WRFChem/O3"
CMAQ_file = "/Volumes/DataStorage/DataFusionData/WRFChem/WRFChem_Combine_2017.nc"
obs_file =  "/Volumes/DataStorage/DataFusionData/WRFChem/O3_obs_2017_CONUS_ppm.csv"
CMAQ_pollutant = "O3_MDA8"
obs_pollutant = "O3"
k_fold = 5

geo = WRFGridInfo(CMAQ_file)
CMAQList = extractCMAQ(CMAQ_pollutant, CMAQ_file, geo)

obs_df = pd.read_csv(obs_file, parse_dates=["Time"], date_parser=dateparse)
obs_df = obs_df[(obs_df["Time"] >= datetime(2017, 1, 1)) & (obs_df["Time"] <= datetime(2017, 12, 31))]
kf = KFold(n_splits=k_fold, shuffle=True, random_state=100)
obs_train_filenames = []
obs_test_filenames = []
obs_site_id = obs_df["siteCode"].unique()

for i, (train_index, test_index) in enumerate(kf.split(obs_site_id)):
    obs_train_filename = os.path.join(work_dir, "train_obs_sfold_" + str(i) + ".csv")
    obs_test_filename = os.path.join(work_dir, "test_obs_sfold_" + str(i) + ".csv")
    train_site, test_site = obs_site_id[train_index], obs_site_id[test_index]
    train_df = obs_df[obs_df['siteCode'].isin(train_site)]
    test_df = obs_df[obs_df['siteCode'].isin(test_site)]
    train_df = train_df.sort_values(by='Time')
    test_df = test_df.sort_values(by='Time')
    train_df.to_csv(obs_train_filename, index=False)
    test_df.to_csv(obs_test_filename, index=False)
    obs_train_filenames.append(obs_train_filename)
    obs_test_filenames.append(obs_test_filename)

for i in range(0, k_fold):
    data_fused_filename = os.path.join(work_dir, "WRFChem_Combine_2017_o3_fused_fold_site_" + str(i) + ".nc")
    train_obs_list = extractOBS(obs_pollutant, obs_train_filenames[i], geo)
    datafusion_results = dataFusion(CMAQList, train_obs_list, geo, 'default')
    writeToNetCDF(datafusion_results, CMAQ_pollutant, CMAQ_file, data_fused_filename)
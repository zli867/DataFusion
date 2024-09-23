from DataExtraction.ExtractCMAQ import extractCMAQ, CMAQGridInfo
from DataExtraction.ExtractOBS import extractOBS
from DataFusion.DataFusionDriver import dataFusion
from DataWriter.DataWriter import combineFusionResults
import netCDF4
import pandas as pd
from sklearn.model_selection import KFold
from Evaluation.StatisticalMetrics import stats_metrics
import os

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
            if name in ["TFLAG"]:
                x = dst.createVariable(name, variable.datatype, variable.dimensions)
                dst[name].setncatts(src[name].__dict__)
                dst[name][:] = src[name][:]

# random withhoding
work_dir =  "/home/zli867/DFManuscript/DataFusion/data/validation/PM25"
CMAQ_file = "/home/zli867/DataFusion/data/CMAQ/CCTM_ACONC_combined_EQUATES_2010_2019.nc"
obs_file =  "/home/zli867/DataFusion/data/CMAQ/PM25_obs_2010_2019_CONUS_ug.csv"
CMAQ_pollutant = "PM25_AVG"
obs_pollutant = "PM25"
k_fold = 5

geo = CMAQGridInfo(CMAQ_file)
CMAQList = extractCMAQ(CMAQ_pollutant, CMAQ_file, geo)

obs_df = pd.read_csv(obs_file)
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
    print("Fuse %d Fold" % i)
    data_fused_filename = os.path.join(work_dir, "CTM_ACONC_fused_fold_" + obs_pollutant + "_" + str(i) + ".nc4")
    train_obs_list = extractOBS(obs_pollutant, obs_train_filenames[i], geo)
    datafusion_results = dataFusion(CMAQList, train_obs_list, geo, 'default')
    writeToNetCDF(datafusion_results, CMAQ_pollutant, CMAQ_file, data_fused_filename)


# # site wide withholding
# work_dir =  "/home/zli867/DFManuscript/DataFusion/data/validation/PM25"
# CMAQ_file = "/home/zli867/DataFusion/data/CMAQ/CCTM_ACONC_combined_EQUATES_2010_2019.nc"
# obs_file =  "/home/zli867/DataFusion/data/CMAQ/PM25_obs_2010_2019_CONUS_ug.csv"
# CMAQ_pollutant = "PM25_AVG"
# obs_pollutant = "PM25"
# k_fold = 5

# geo = CMAQGridInfo(CMAQ_file)
# CMAQList = extractCMAQ(CMAQ_pollutant, CMAQ_file, geo)

# obs_df = pd.read_csv(obs_file)
# kf = KFold(n_splits=k_fold, shuffle=True, random_state=100)
# obs_train_filenames = []
# obs_test_filenames = []
# obs_site_id = obs_df["siteCode"].unique()

# for i, (train_index, test_index) in enumerate(kf.split(obs_site_id)):
#     obs_train_filename = os.path.join(work_dir, "train_obs_sfold_" + str(i) + ".csv")
#     obs_test_filename = os.path.join(work_dir, "test_obs_sfold_" + str(i) + ".csv")
#     train_site, test_site = obs_site_id[train_index], obs_site_id[test_index]
#     train_df = obs_df[obs_df['siteCode'].isin(train_site)]
#     test_df = obs_df[obs_df['siteCode'].isin(test_site)]
#     train_df = train_df.sort_values(by='Time')
#     test_df = test_df.sort_values(by='Time')
#     train_df.to_csv(obs_train_filename, index=False)
#     test_df.to_csv(obs_test_filename, index=False)
#     obs_train_filenames.append(obs_train_filename)
#     obs_test_filenames.append(obs_test_filename)

# for i in range(0, k_fold):
#     print("Fuse %d Fold" % i)
#     data_fused_filename = os.path.join(work_dir, "CTM_ACONC_fused_fold_site_" + obs_pollutant + "_" + str(i) + ".nc")
#     train_obs_list = extractOBS(obs_pollutant, obs_train_filenames[i], geo)
#     datafusion_results = dataFusion(CMAQList, train_obs_list, geo, 'default')
#     writeToNetCDF(datafusion_results, CMAQ_pollutant, CMAQ_file, data_fused_filename)
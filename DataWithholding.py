from DataExtraction.ExtractCMAQ import extractCMAQ, CMAQGridInfo
from DataExtraction.ExtractOBS import extractOBS
from DataFusion.DataFusionDriver import dataFusion
from DataWriter.DataWriter import writeToNetCDF
import pandas as pd
from sklearn.model_selection import KFold
from Evaluation.StatisticalMetrics import stats_metrics

work_dir = "./work_dir/"
CMAQ_file = "./data/CCTM.ACONC.combined.FIRE.hires4.2016_2020.nc"
obs_file = "./data/obs_2016_2020.csv"
CMAQ_pollutant = "PM25_TOT_AVG"
obs_pollutant = "PM25"
k_fold = 10

CMAQList = extractCMAQ(CMAQ_pollutant, CMAQ_file)
geo = CMAQGridInfo(CMAQ_file)
OBSList = extractOBS(obs_pollutant, obs_file, geo)
cmaq_performance = stats_metrics(CMAQList[0], OBSList[0], geo)
print("CMAQ performance: ")
print(cmaq_performance)


obs_df = pd.read_csv(obs_file)
kf = KFold(n_splits=k_fold, shuffle=True, random_state=100)
data_fused_filenames = []
obs_train_filenames = []
obs_test_filenames = []
for i, (train_index, test_index) in enumerate(kf.split(obs_df)):
    obs_train_filename = work_dir + "train_obs_fold_" + str(i) + ".csv"
    obs_test_filename = work_dir + "test_obs_fold_" + str(i) + ".csv"
    train_df = obs_df.iloc[train_index]
    test_df = obs_df.iloc[test_index]
    train_df = train_df.sort_values(by='Time')
    test_df = test_df.sort_values(by='Time')
    train_df.to_csv(obs_train_filename, index=False)
    test_df.to_csv(obs_test_filename, index=False)
    obs_train_filenames.append(obs_train_filename)
    obs_test_filenames.append(obs_test_filename)

for i in range(0, k_fold):
    data_fused_filename = work_dir + "CCTM.ACONC.combined.fused_fold_" + str(i) + ".nc"
    train_obs_list = extractOBS(obs_pollutant, obs_train_filenames[i], geo)
    datafusion_results = dataFusion(CMAQList, train_obs_list, geo)
    writeToNetCDF(datafusion_results, CMAQ_pollutant, CMAQ_file, data_fused_filename)
    # evaluation
    test_obs_list = extractOBS(obs_pollutant, obs_test_filenames[i], geo)
    fused_CMAQ_List = extractCMAQ(CMAQ_pollutant, data_fused_filename)
    fused_cmaq_performance = stats_metrics(fused_CMAQ_List[0], test_obs_list[0], geo)
    print("Data fusion performance " + str(i) + "-fold: ")
    print(fused_cmaq_performance)
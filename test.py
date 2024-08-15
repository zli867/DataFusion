# from sklearn.model_selection import KFold
# import numpy as np
# from DataFusion.adjustCMAQ import adjustParameters
# from DataFusion.adjustCMAQ import adjustParamentersLinear
# from DataFusion.adjustCMAQ import adjustCMAQ
# from Evaluation.StatisticalMetrics import RMSE

# obs_yearly_mean_all, CMAQ_at_obs_yearly_mean_all = np.random.rand((1234)), np.random.rand((1234))
# # cross validation
# non_linear_score, linear_score = [], []
# kf = KFold(n_splits=10, random_state=100, shuffle=True)
# for i, (train_index, test_index) in enumerate(kf.split(obs_yearly_mean_all)):
#     cur_obs_train, cur_cmaq_train = obs_yearly_mean_all[train_index], CMAQ_at_obs_yearly_mean_all[train_index]
#     cur_obs_test, cur_cmaq_test = obs_yearly_mean_all[test_index], CMAQ_at_obs_yearly_mean_all[test_index]
#     # non linear
#     alpha_nonlinear, beta_nonlinear, _ = adjustParameters(cur_cmaq_train,cur_obs_train)
#     nonlinear_predict = adjustCMAQ(cur_cmaq_test, alpha_nonlinear, beta_nonlinear)
#     cur_nonlinear_score = RMSE(nonlinear_predict, cur_obs_test)
#     non_linear_score.append(cur_nonlinear_score)
#     # linear
#     alpha_linear, beta_linear, _ = adjustParamentersLinear(cur_cmaq_train, cur_obs_train)
#     linear_predict = adjustCMAQ(cur_cmaq_test, alpha_linear, beta_linear)
#     cur_linear_score = RMSE(linear_predict, cur_obs_test)
#     linear_score.append(cur_linear_score)
# print(np.mean(non_linear_score))
# print(np.mean(linear_score))

import netCDF4 as nc
import numpy as np
filename = "/Volumes/DataStorage/DataFusionData/GC/GC_Combine_2015.nc4"
# filename = "/Volumes/Expansion/DataFusionData/GC/2015/GEOSChem.SpeciesConc.20150621_0000z.nc4"
ds = nc.Dataset(filename)
o3 = ds["O3_MDA8"][:] * 1e9
# # o3 = ds["SpeciesConc_O3"][0:10, 0, :, :] * 1e9
import matplotlib.pyplot as plt

for i in range(0, len(o3)):
    cur_o3 = np.squeeze(o3[i, :, :, :])
    plt.pcolor(cur_o3)
    plt.colorbar()
    plt.show()

# # # import pandas as pd
# filename = "/Volumes/Expansion/DataFusionData/GC/2015/GEOSChem.SpeciesConc.20150621_0000z.nc4"
# ds = nc.Dataset(filename)
# o3 =  ds["SpeciesConc_O3"][:] * 1e9
# print(np.max(o3))
# file_format = ds.file_format
# print(file_format)
# filename = "/Volumes/DataStorage/DataFusionData/obs/O3_obs_2010_2019_CONUS_ppb.csv"
# df = pd.read_csv(filename)
# df["O3"] = df["O3"] * 1e-9
# df.to_csv("/Volumes/DataStorage/DataFusionData/obs/O3_obs_2010_2019_CONUS_GC.csv", index=False)
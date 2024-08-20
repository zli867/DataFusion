# # from sklearn.model_selection import KFold
# # import numpy as np
# # from DataFusion.adjustCMAQ import adjustParameters
# # from DataFusion.adjustCMAQ import adjustParamentersLinear
# # from DataFusion.adjustCMAQ import adjustCMAQ
# # from Evaluation.StatisticalMetrics import RMSE

# # obs_yearly_mean_all, CMAQ_at_obs_yearly_mean_all = np.random.rand((1234)), np.random.rand((1234))
# # # cross validation
# # non_linear_score, linear_score = [], []
# # kf = KFold(n_splits=10, random_state=100, shuffle=True)
# # for i, (train_index, test_index) in enumerate(kf.split(obs_yearly_mean_all)):
# #     cur_obs_train, cur_cmaq_train = obs_yearly_mean_all[train_index], CMAQ_at_obs_yearly_mean_all[train_index]
# #     cur_obs_test, cur_cmaq_test = obs_yearly_mean_all[test_index], CMAQ_at_obs_yearly_mean_all[test_index]
# #     # non linear
# #     alpha_nonlinear, beta_nonlinear, _ = adjustParameters(cur_cmaq_train,cur_obs_train)
# #     nonlinear_predict = adjustCMAQ(cur_cmaq_test, alpha_nonlinear, beta_nonlinear)
# #     cur_nonlinear_score = RMSE(nonlinear_predict, cur_obs_test)
# #     non_linear_score.append(cur_nonlinear_score)
# #     # linear
# #     alpha_linear, beta_linear, _ = adjustParamentersLinear(cur_cmaq_train, cur_obs_train)
# #     linear_predict = adjustCMAQ(cur_cmaq_test, alpha_linear, beta_linear)
# #     cur_linear_score = RMSE(linear_predict, cur_obs_test)
# #     linear_score.append(cur_linear_score)
# # print(np.mean(non_linear_score))
# # print(np.mean(linear_score))

# # import netCDF4 as nc
# # import numpy as np
# # filename = "/Volumes/DataStorage/DataFusionData/GC/GC_Combine_2015.nc4"
# # # filename = "/Volumes/Expansion/DataFusionData/GC/2015/GEOSChem.SpeciesConc.20150621_0000z.nc4"
# # ds = nc.Dataset(filename)
# # o3 = ds["NO2_AVG"][:] * 1e9
# # # # o3 = ds["SpeciesConc_O3"][0:10, 0, :, :] * 1e9
# # import matplotlib.pyplot as plt

# # for i in range(0, len(o3)):
# #     cur_o3 = np.squeeze(o3[i, :, :, :])
# #     plt.pcolor(cur_o3)
# #     plt.colorbar()
# #     plt.show()

# # # # import pandas as pd
# # filename = "/Volumes/Expansion/DataFusionData/GC/2015/GEOSChem.SpeciesConc.20150621_0000z.nc4"
# # ds = nc.Dataset(filename)
# # o3 =  ds["SpeciesConc_O3"][:] * 1e9
# # print(np.max(o3))
# # file_format = ds.file_format
# # print(file_format)
# # filename = "/Volumes/DataStorage/DataFusionData/obs/O3_obs_2010_2019_CONUS_ppb.csv"
# # df = pd.read_csv(filename)
# # df["O3"] = df["O3"] * 1e-9
# # df.to_csv("/Volumes/DataStorage/DataFusionData/obs/O3_obs_2010_2019_CONUS_GC.csv", index=False)

# # import pandas as pd
# # filename = "/Volumes/DataStorage/DataFusionData/CMAQ/O3_obs_2010_2019_CONUS_ppb.csv"
# # df = pd.read_csv(filename)
# # df = df.sample(frac=0.1).reset_index(drop=True)
# # df.to_csv("/Volumes/DataStorage/DataFusionData/CMAQ/O3_obs_2010_2019_sample_ppb.csv", index=False)

# # import pandas as pd
# # filename = "/Volumes/DataStorage/DataFusionData/obs/NO2_obs_2010_2019_CONUS_ppb.csv"
# # df = pd.read_csv(filename)
# # df["NO2"] = df["NO2"] * 1e-9
# # df.to_csv("/Volumes/DataStorage/DataFusionData/obs/NO2_obs_2010_2019_CONUS_GC.csv", index=False)

# # import concurrent.futures
# # import numpy as np
# # t = np.zeros((10000, 500, 500))

# # def process_date(i):
# #     t[i, :, :] = np.ones((500, 500))
    

# # # List of processed years
# # processed_years = np.arange(365)

# # # Parallelize the loop
# # with concurrent.futures.ThreadPoolExecutor() as executor:
# #     executor.map(process_date(processed_years), processed_years)

# # print(np.min(t))

# # SuperFastPython.com
# # example of calling map and processing results
# # from time import sleep
# # from random import random
# # from concurrent.futures import ThreadPoolExecutor
 
# # # custom task that will sleep for a variable amount of time
# # def task(name):
# #     # sleep for less than a second
# #     sleep(random())
# #     return f'Task: {name} done.'
 
# # # start the thread pool
# # with ThreadPoolExecutor(10) as executor:
# #     # execute tasks concurrently and process results in order
# #     for result in executor.map(task, range(10)):
# #         # report the result
# #         print(result)
# # print("hello")

# # import numpy as np
# # from DataFusion.KrigingOBS import krigingOBS
# # from scipy.optimize import curve_fit
# # import concurrent.futures
# # from multiprocessing import 

# # def dataFusionOne(CMAQDict, OBSDict, geoDict, alpha_yearly, beta):
# #     # Use CMAQAdjust function to calculate adjusted CMAQ
# #     FCData = FC(CMAQDict["Yearly"], alpha_yearly, beta)
# #     result = np.empty(CMAQDict["Daily"].shape)
# #     result[:] = np.NaN
# #     CMAQTime = CMAQDict["Time"]
# #     dateRange = result.shape[0]
# #     # Use kriging to calculate ratio
# #     OBSConc = OBSDict["Conc"]
# #     OBSConcMean = np.nanmean(OBSConc, axis=0)
# #     OBSTime = OBSDict["dateSeries"]
# #     obsX = OBSDict["X"]
# #     obsY = OBSDict["Y"]
# #     predictX = geoDict["X"]
# #     predictY = geoDict["Y"]
# #     # TODO: Use parallel or not?
# #     def fusion_daily(i):
# #         print(i)
# #         current_time = CMAQTime[i]
# #         if current_time in OBSTime:
# #             valid_krig = True
# #             obs_time_idx = OBSTime.index(current_time)
# #             OBSNormalized = OBSConc[obs_time_idx, :] / OBSConcMean
# #             # Krig method to interpolate the observation
# #             # delete the Nan number in array
# #             nanIndex = np.isnan(OBSNormalized)
# #             OBSNormalized = OBSNormalized[~nanIndex]
# #             siteX = obsX[~nanIndex]
# #             siteY = obsY[~nanIndex]
# #             # If there is less than one observation, do not do kriging
# #             # if the data variance is zero (constant), do not do kriging
# #             if OBSNormalized.shape[0] == 0 or OBSNormalized.shape[0] == 1 or np.var(OBSNormalized) == 0:
# #                 valid_krig = False

# #             if valid_krig:
# #                 krigRatio = krigingOBS(siteX, siteY, OBSNormalized, predictX, predictY)
# #                 resultDaily = krigRatio * FCData
# #             else:
# #                 resultDaily = FCData
# #         else:
# #             resultDaily = FCData
# #         result[i, :, :] = resultDaily
# #     all_dates_indices = np.arange(0, dateRange)
# #     # Parallelize the loop
# #     with concurrent.futures.ThreadPoolExecutor(max_workers=2) as executor:
# #         executor.map(fusion_daily, all_dates_indices)
        
# #     # for i in range(0, dateRange):
# #     #     current_time = CMAQTime[i]
# #     #     if current_time in OBSTime:
# #     #         valid_krig = True
# #     #         obs_time_idx = OBSTime.index(current_time)
# #     #         OBSNormalized = OBSConc[obs_time_idx, :] / OBSConcMean
# #     #         # Krig method to interpolate the observation
# #     #         # delete the Nan number in array
# #     #         nanIndex = np.isnan(OBSNormalized)
# #     #         OBSNormalized = OBSNormalized[~nanIndex]
# #     #         siteX = obsX[~nanIndex]
# #     #         siteY = obsY[~nanIndex]
# #     #         # If there is less than one observation, do not do kriging
# #     #         # if the data variance is zero (constant), do not do kriging
# #     #         if OBSNormalized.shape[0] == 0 or OBSNormalized.shape[0] == 1 or np.var(OBSNormalized) == 0:
# #     #             valid_krig = False

# #     #         if valid_krig:
# #     #             krigRatio = krigingOBS(siteX, siteY, OBSNormalized, predictX, predictY)
# #     #             resultDaily = krigRatio * FCData
# #     #         else:
# #     #             resultDaily = FCData
# #     #     else:
# #     #         resultDaily = FCData
# #     #     result[i, :, :] = resultDaily
# #     return result



# # import concurrent.futures
# # import numpy as np
# # import multiprocessing as mp 
# # t = np.zeros((10000, 500, 500))

# # def process_date(i):
# #     t[i, :, :] = np.ones((500, 500))
    

# # # List of processed years
# # processed_years = np.arange(365)
# # pool = mp.Pool(mp.cpu_count()) 
# # result = pool.map(process_date, processed_years)


# # import concurrent.futures
# # import numpy as np
# # x = np.array([[1, 2], [3, 4]])
# # arr_reshaped = x[np.newaxis, :, :]

# # # Repeat the reshaped array 10 times along the first dimension
# # arr_3d = np.repeat(arr_reshaped, repeats=3, axis=0)
# # print(arr_3d)

# import multiprocess as mp

# # Define the function to process data for each year
# def process_year(year):
#     def yest(k):
#         return k + 1
#     print(f"Data fusion for {year}")
#     return (year, yest(year))

# # Set up multiprocessing pool

# # List of years to process
# processed_years = [2021, 2022]  # Your list of years

# # Create a pool of workers
# with mp.Pool(mp.cpu_count()) as pool:
#     print(mp.cpu_count())
#     # Map the process_year function to the years
#     results = pool.map(process_year, processed_years)
# print(results)
# # Combine results into a dictionary
# fused_conc = dict(results)

# print(fused_conc)

import netCDF4 as nc
import numpy as np

species = "NO2_AVG"
filename1 = "/Volumes/DataStorage/DataFusionData/fused/GC_Combine_2015_fused.nc4"
ds1 = nc.Dataset(filename1)
filename2 = "/Volumes/DataStorage/DataFusionData/fused/GC_Combine_2015_fused_test.nc4"
ds2 = nc.Dataset(filename2)
value_1 = ds1[species][:]
value_2 = ds2[species][:]
diff = np.abs(value_1 - value_2)
print(np.max(diff) * 1e9)
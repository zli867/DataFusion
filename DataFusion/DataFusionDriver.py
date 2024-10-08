from DataFusion.CMAQatOBS import CMAQatOBSDaily
import numpy as np
from DataFusion.adjustCMAQ import adjustParamentersYearly
from DataFusion.adjustCMAQ import adjustParameters
from DataFusion.adjustCMAQ import adjustParamentersLinear
from DataFusion.CorrelationParameter import ROBSData
from DataFusion.CorrelationParameter import ROBS
from DataFusion.CorrelationParameter import RCMAQ
from DataFusion.dataFusion import dataFusionOne
from DataFusion.dataFusion import dataFusionTwo
from DataFusion.dataFusion import temporalCorrection
from DataFusion.dataFusion import FC
from DataFusion.CorrelationParameter import R1
from DataFusion.dataFusion import weightFactor
import sys, os
from sklearn.model_selection import KFold
from DataFusion.adjustCMAQ import adjustCMAQ
from Evaluation.StatisticalMetrics import RMSE
from concurrent.futures import ProcessPoolExecutor, as_completed


def dataFusion(CMAQList, OBSList, geo, method='default', worker=1):
    # Generate Alpha and Beta
    combined_CMAQ_dict = CMAQList[0]
    combined_obs_dict = OBSList[0]
    yearly_CMAQ_dict = CMAQList[1]
    yearly_obs_dict = OBSList[1]

    processed_years = list(yearly_CMAQ_dict.keys())
    processed_years.sort()

    for obs_year in yearly_obs_dict.keys():
        if obs_year not in processed_years:
            sys.exit("Missing observation data for " + str(obs_year) + ", please add.")

    obs_yearly_mean_all = []
    CMAQ_at_obs_yearly_mean_all = []
    CMAQ_at_obs_yearly_mean = {}
    CMAQ_at_obs_spatial_mean = {}
    obs_mean_yearly = {}
    obs_mean_spatial = {}
    alpha_yearly = {}
    for select_year in processed_years:
        current_obs_dict = yearly_obs_dict[select_year]
        current_CMAQ_dict = yearly_CMAQ_dict[select_year]
        current_CMAQ_at_obs_conc = CMAQatOBSDaily(current_obs_dict, geo, current_CMAQ_dict)
        # Temporal mean values
        current_obs_yearly_mean = np.nanmean(current_obs_dict["Conc"], axis=0)
        current_CMAQ_at_obs_yearly_mean = np.nanmean(current_CMAQ_at_obs_conc, axis=0)
        CMAQ_at_obs_yearly_mean[select_year] = current_CMAQ_at_obs_yearly_mean
        obs_mean_yearly[select_year] = current_obs_yearly_mean
        # Spatial mean values
        current_CMAQ_at_obs_spatial_mean = np.nanmean(current_CMAQ_at_obs_conc, axis=1)
        current_obs_spatial_mean = np.nanmean(current_obs_dict["Conc"], axis=1)
        CMAQ_at_obs_spatial_mean[select_year] = current_CMAQ_at_obs_spatial_mean
        obs_mean_spatial[select_year] = current_obs_spatial_mean

        obs_yearly_mean_all.append(current_obs_yearly_mean)
        CMAQ_at_obs_yearly_mean_all.append(current_CMAQ_at_obs_yearly_mean)

    obs_yearly_mean_all = np.concatenate(obs_yearly_mean_all)
    CMAQ_at_obs_yearly_mean_all = np.concatenate(CMAQ_at_obs_yearly_mean_all)
    alpha_nonlinear, beta_nonlinear, score_nonlinear = adjustParameters(CMAQ_at_obs_yearly_mean_all,
                                                                        obs_yearly_mean_all)
    alpha_linear, beta_linear, score_linear = adjustParamentersLinear(CMAQ_at_obs_yearly_mean_all, obs_yearly_mean_all)
    # calculate alpha and beta based on method
    if method == 'linear':
        alpha = alpha_linear
        beta = beta_linear
    elif method == 'nonlinear':
        alpha = alpha_nonlinear
        beta = beta_nonlinear
    else:
        # cross validation for model selection if we have enough data
        if len(CMAQ_at_obs_yearly_mean_all) >= 10:
            non_linear_score, linear_score = [], []
            kf = KFold(n_splits=10, random_state=100, shuffle=True)
            for i, (train_index, test_index) in enumerate(kf.split(obs_yearly_mean_all)):
                cur_obs_train, cur_cmaq_train = obs_yearly_mean_all[train_index], CMAQ_at_obs_yearly_mean_all[train_index]
                cur_obs_test, cur_cmaq_test = obs_yearly_mean_all[test_index], CMAQ_at_obs_yearly_mean_all[test_index]
                # non linear
                alpha_nonlinear, beta_nonlinear, _ = adjustParameters(cur_cmaq_train,cur_obs_train)
                nonlinear_predict = adjustCMAQ(cur_cmaq_test, alpha_nonlinear, beta_nonlinear)
                cur_nonlinear_score = RMSE(nonlinear_predict, cur_obs_test)
                non_linear_score.append(cur_nonlinear_score)
                # linear
                alpha_linear, beta_linear, _ = adjustParamentersLinear(cur_cmaq_train, cur_obs_train)
                linear_predict = adjustCMAQ(cur_cmaq_test, alpha_linear, beta_linear)
                cur_linear_score = RMSE(linear_predict, cur_obs_test)
                linear_score.append(cur_linear_score)
            score_linear = np.mean(linear_score)
            score_nonlinear = np.mean(non_linear_score)
        if score_linear >= score_nonlinear:
            print("Select linear adjustment for CMAQ, RMSE score is: " + str(score_linear))
            alpha = alpha_linear
            beta = beta_linear
        else:
            print("Select nonlinear adjustment for CMAQ, RMSE score is: " + str(score_nonlinear))
            alpha = alpha_nonlinear
            beta = beta_nonlinear

    for select_year in processed_years:
        current_alpha = adjustParamentersYearly(beta, CMAQ_at_obs_yearly_mean[select_year],
                                                obs_mean_yearly[select_year])
        alpha_yearly[select_year] = current_alpha
        print("alpha " + str(select_year) + " :" + str(current_alpha))

    print("alpha: " + str(alpha))
    print("beta: " + str(beta))

    # Calculate Correlation Parameters
    # Calculate parameters
    [corrCoef, distance] = ROBSData(combined_obs_dict["Conc"], combined_obs_dict["X"], combined_obs_dict["Y"])
    [r, Rcoll] = ROBS(corrCoef, distance)

    print("r: " + str(r))
    print("Rcoll: " + str(Rcoll))

    CMAQ_at_obs_combined = CMAQatOBSDaily(combined_obs_dict, geo, combined_CMAQ_dict)
    R2 = RCMAQ(combined_obs_dict["Conc"], CMAQ_at_obs_combined)
    print("R2: " + str(R2))

    all_obs_time_series = []
    all_obs_adjusted_CMAQ_ratio = []
    # Calculate A and t_max
    for year in processed_years:
        current_alpha = alpha_yearly[year]
        yearly_obs = yearly_obs_dict[year]
        # Calculate A and t_max
        current_obs_time = yearly_obs["dateSeries"]
        current_CMAQ_at_obs_spatial_mean = CMAQ_at_obs_spatial_mean[year]
        current_obs_spatial_mean = obs_mean_spatial[year]
        current_adjusted_CMAQ_at_obs_spatial_mean = FC(current_CMAQ_at_obs_spatial_mean, current_alpha, beta)
        current_obs_adjusted_CMAQ_ratio = current_obs_spatial_mean/current_adjusted_CMAQ_at_obs_spatial_mean
        all_obs_time_series.extend(current_obs_time)
        all_obs_adjusted_CMAQ_ratio.append(current_obs_adjusted_CMAQ_ratio)
    all_obs_adjusted_CMAQ_ratio = np.concatenate(all_obs_adjusted_CMAQ_ratio, axis=0)
    A, t_max = temporalCorrection(all_obs_adjusted_CMAQ_ratio, all_obs_time_series)
    print("A:" + str(A) + " , t_max: " + str(t_max))

    global_params = {
        "A": A,
        "t_max": t_max,
        "beta": beta,
        "Rcoll": Rcoll,
        "r": r,
        "R2": R2
    }
    
    fused_conc = {}
    if worker == 1:
        for year in processed_years:
            cur_conc = yearly_fusion(year, alpha_yearly, yearly_obs_dict, yearly_CMAQ_dict, geo, global_params)
            fused_conc[year] = cur_conc
    else:
        input_data = [(year, alpha_yearly, yearly_obs_dict, yearly_CMAQ_dict, geo, global_params) for year in processed_years]
        # Use ProcessPoolExecutor for parallel processing
        default_workers = min(32, os.cpu_count() + 4)
        worker_count = min(worker, default_workers)
        with ProcessPoolExecutor(max_workers=worker_count) as executor:
            futures = {executor.submit(yearly_fusion, *args): args[0] for args in input_data}
            for future in as_completed(futures):
                year = futures[future]
                fused_conc[year] = future.result()
    return fused_conc


def yearly_fusion(year, alpha_yearly, yearly_obs_dict, yearly_CMAQ_dict, geo, global_params):
    print("Data fusion for " + str(year))
    current_alpha = alpha_yearly[year]
    yearly_obs = yearly_obs_dict[year]
    yearly_CMAQ = yearly_CMAQ_dict[year]
    FC_1 = dataFusionOne(yearly_CMAQ, yearly_obs, geo, current_alpha, global_params["beta"])
    FC_2 = dataFusionTwo(yearly_CMAQ, current_alpha, global_params["beta"], global_params["A"], global_params["t_max"])

    R_1 = R1(yearly_obs, geo, global_params["Rcoll"], global_params["r"], year)
    W = weightFactor(R_1, global_params["R2"])
    FC_opt = W * FC_1 + (1 - W) * FC_2
    yearly_fused_conc = FC_opt
    return yearly_fused_conc
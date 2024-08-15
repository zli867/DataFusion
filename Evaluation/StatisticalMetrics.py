import scipy.stats
from DataFusion.CMAQatOBS import CMAQatOBSDaily
import numpy as np
from sklearn.metrics import mean_squared_error

# reference: https://www.tandfonline.com/doi/full/10.1080/10962247.2016.1265027
def remove_nan_values(prediction, observation):
    valid_idx = (~np.isnan(prediction)) & (~np.isnan(observation))
    valid_p = prediction[valid_idx]
    valid_o = observation[valid_idx]
    return valid_p, valid_o


def MB(prediction, observation):
    return np.sum(prediction - observation) / len(prediction)


def ME(prediction, observation):
    return np.sum(np.abs(prediction - observation)) / len(prediction)


def RMSE(prediction, observation):
    return np.sqrt(np.sum((prediction - observation) ** 2) / len(prediction))


def CRMSE(prediction, observation):
    p_mean, o_mean = np.mean(prediction), np.mean(observation)
    return np.sqrt((1 / len(prediction) * np.sum(((prediction - p_mean) - (observation - o_mean)) ** 2)))


def NMB(prediction, observation):
    return (np.sum(prediction - observation) / np.sum(observation)) * 100


def NME(prediction, observation):
    return ((np.sum(np.abs(prediction - observation))) / np.sum(observation)) * 100


def MNB(prediction, observation):
    cur_predict, cur_obs = prediction[observation > 0], observation[observation > 0]
    return (1 / len(cur_predict)) * np.sum((cur_predict - cur_obs) / cur_obs) * 100


def MNE(prediction, observation):
    cur_predict, cur_obs = prediction[observation > 0], observation[observation > 0]
    return (1 / len(cur_predict)) * np.sum(np.abs(cur_predict - cur_obs) / cur_obs) * 100


def FB(prediction, observation):
    return (2 / len(prediction)) * np.sum((prediction - observation) / (prediction + observation)) * 100


def FE(prediction, observation):
    return (2 / len(prediction)) * np.sum(np.abs(prediction - observation) / (prediction + observation)) * 100


def IOA(prediction, observation):
    numerator = np.sum((prediction - observation) ** 2)
    prediction_shift = np.abs((prediction - np.mean(observation)))
    observation_shift = np.abs((observation - np.mean(observation)))
    denominator = np.sum((prediction_shift + observation_shift) ** 2)
    return 1 - numerator / denominator


def pearson_r(prediction, observation):
    r = scipy.stats.pearsonr(prediction, observation)
    return r[0]


def stats_metrics(combined_CMAQ_dict, combined_obs_dict, geo, metrics_dict):
    obs_conc = combined_obs_dict["Conc"]
    cmaq_at_obs_conc = CMAQatOBSDaily(combined_obs_dict, geo, combined_CMAQ_dict)
    siteSize = obs_conc.shape[1]
    obs_vals, model_vals = [], []
    for i in range(0, siteSize):
        CMAQConcTmp = cmaq_at_obs_conc[:, i]
        OBSConcTmp = obs_conc[:, i]
        valid_idx = (~np.isnan(CMAQConcTmp)) & (~np.isnan(OBSConcTmp))
        # Delete the nan number
        CMAQConcTmp = CMAQConcTmp[valid_idx]
        OBSConcTmp = OBSConcTmp[valid_idx]
        obs_vals.append(OBSConcTmp)
        model_vals.append(CMAQConcTmp)
    obs_vals = np.concatenate(obs_vals)
    model_vals = np.concatenate(model_vals)
    performance = {"metrics": [], "values": []}
    for metric_name in metrics_dict.keys():
        performance["metrics"].append(metric_name)
        performance["values"].append(metrics_dict[metric_name]["func"](model_vals, obs_vals))
    return performance, model_vals, obs_vals


def monitor_stats_metrics(combined_CMAQ_dict, combined_obs_dict, geo, metrics_dict):
    obs_conc = combined_obs_dict["Conc"]
    cmaq_at_obs_conc = CMAQatOBSDaily(combined_obs_dict, geo, combined_CMAQ_dict)
    siteSize = obs_conc.shape[1]
    performance = {"Lat": [], "Lon": []}
    for metric_name in metrics_dict.keys():
        performance[metric_name] = []
    
    for i in range(0, siteSize):
        CMAQConcTmp = cmaq_at_obs_conc[:, i]
        OBSConcTmp = obs_conc[:, i]
        valid_idx = (~np.isnan(CMAQConcTmp)) & (~np.isnan(OBSConcTmp))
        # Delete the nan number
        CMAQConcTmp = CMAQConcTmp[valid_idx]
        OBSConcTmp = OBSConcTmp[valid_idx]
        if len(CMAQConcTmp) > 1:
            valid_flag, cur_metrics = True, {}
            for metric_name in metrics_dict.keys():
                cur_val = metrics_dict[metric_name]["func"](CMAQConcTmp, OBSConcTmp)
                if np.isnan(cur_val):
                    valid_flag = False
                    break
                else:
                    cur_metrics[metric_name] = cur_val
            if valid_flag:
                for metric_name in metrics_dict.keys():
                    performance[metric_name].append(cur_metrics[metric_name])
                performance["Lat"].append(combined_obs_dict["Lat"][i])
                performance["Lon"].append(combined_obs_dict["Lon"][i])
    return performance
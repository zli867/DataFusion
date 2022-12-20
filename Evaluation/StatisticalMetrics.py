import scipy.stats
from DataFusion.CMAQatOBS import CMAQatOBSDaily
import numpy as np
from sklearn.metrics import mean_squared_error


def stats_metrics(combined_CMAQ_dict, combined_obs_dict, geo):
    obs_conc = combined_obs_dict["Conc"]
    cmaq_at_obs_conc = CMAQatOBSDaily(combined_obs_dict, geo, combined_CMAQ_dict)
    siteSize = obs_conc.shape[1]

    corr_coef = []
    rmse = []
    for i in range(0, siteSize):
        CMAQConcTmp = cmaq_at_obs_conc[:, i]
        OBSConcTmp = obs_conc[:, i]
        valid_idx = (~np.isnan(CMAQConcTmp)) & (~np.isnan(OBSConcTmp))
        # Delete the nan number
        CMAQConcTmp = CMAQConcTmp[valid_idx]
        OBSConcTmp = OBSConcTmp[valid_idx]
        if len(CMAQConcTmp) <= 1 or len(OBSConcTmp) <= 1:
            continue
        pearson_r = scipy.stats.pearsonr(CMAQConcTmp, OBSConcTmp)
        rmse_value = mean_squared_error(OBSConcTmp, CMAQConcTmp, squared=False)
        rmse.append(rmse_value)
        corr_coef.append(pearson_r[0])

    performance = {"Person r": np.mean(corr_coef),
                   "RMSE": np.mean(rmse)}
    return performance
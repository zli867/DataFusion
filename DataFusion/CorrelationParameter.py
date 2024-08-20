import numpy as np
import scipy.stats
from scipy.optimize import curve_fit


def exponential(x, A, b):
    return A * np.exp(-b * x)


def RCMAQ(OBSConc, CMAQConc):
    """
    Calculate pearson correlation between CMAQ and observations
    :param OBSConc:
    :param CMAQConc:
    :return:
    """
    siteSize = OBSConc.shape[1]
    corrCoef = []
    for i in range(0, siteSize):
        CMAQConcTmp = CMAQConc[:, i]
        OBSConcTmp = OBSConc[:, i]
        valid_idx = (~np.isnan(CMAQConcTmp)) & (~np.isnan(OBSConcTmp))
        # Delete the nan number
        CMAQConcTmp = CMAQConcTmp[valid_idx]
        OBSConcTmp = OBSConcTmp[valid_idx]
        if len(CMAQConcTmp) <= 1 or len(OBSConcTmp) <= 1:
            continue
        pearson_r = scipy.stats.pearsonr(CMAQConcTmp, OBSConcTmp)
        # There are some reason make pearson r as valid. For example, all data are zeros
        if not np.isnan(pearson_r[0]):
            corrCoef.append(pearson_r[0])

    R2 = np.mean(corrCoef)
    return R2


def ROBSData(combinedOBS, X, Y):
    siteNum = X.shape[0]
    corrCoef = []
    distance = []
    for i in range(0, siteNum - 1):
        for j in range(i+1, siteNum):
            OBS1 = combinedOBS[:, i]
            OBS2 = combinedOBS[:, j]
            # delete the NaN data and maintain the non NaN data for the same day
            naNIndex = np.isnan(OBS1 + OBS2)
            OBS1 = OBS1[~naNIndex]
            OBS2 = OBS2[~naNIndex]

            # we cannot calculate correlation if data less than 2 or all values are zeros <=> mean is zero
            if OBS1.shape[0] <= 1 or OBS2.shape[0] <= 1 or np.mean(OBS1) == 0 or np.mean(OBS2) == 0:
                continue
            # we cannot calculate correlation if data is a constant
            variance_obs1 = (OBS1 - np.mean(OBS1)) / np.mean(OBS1)
            variance_obs2 = (OBS2 - np.mean(OBS2)) / np.mean(OBS2)
            denominator_pearson_r = np.sum(variance_obs1 ** 2) * np.sum(variance_obs2 ** 2)
            if denominator_pearson_r <= 1e-20:
                continue

            pearson_r = scipy.stats.pearsonr(OBS1, OBS2)
            # There are some reason make pearson r as valid. For example, all data are zeros, we do not include for R2 training
            if not np.isnan(pearson_r[0]):
                corrCoef.append(pearson_r[0])
                # calculate the distance change meter to kilometer
                X1 = X[i]/1000
                Y1 = Y[i]/1000
                X2 = X[j]/1000
                Y2 = Y[j]/1000
                d = np.sqrt((X1 - X2) * (X1 - X2) + (Y1 - Y2) * (Y1 - Y2))
                distance.append(d)
    corrCoef = np.array(corrCoef)
    corrCoef = corrCoef.reshape(-1, 1)
    distance = np.array(distance)
    distance = distance.reshape(-1, 1)
    return [corrCoef, distance]


def ROBS(corrCoef, distance):
    ydata = corrCoef.reshape(1, corrCoef.shape[0])
    xdata = distance.reshape(1, distance.shape[0])
    popt, pcov = curve_fit(exponential, xdata[0], ydata[0], p0=(1, 0.01))
    Rcoll = popt[0]
    r = 1/popt[1]
    return [r, Rcoll]

# TODO Need to improve
def R1(obsDict, geo, Rcoll, r, year):
    # If there is no observation, R1 = 0, FC = FC2
    # return size: [CMAQ_yearly_range, m, n]
    # A bug found, 2024/08/15
    CMAQTime = []
    for t in geo["time"]:
        if t.year == year:
            CMAQTime.append(t)
    OBSTime = obsDict["dateSeries"]
    obsConc = obsDict["Conc"]
    obsX = obsDict["X"]
    obsY = obsDict["Y"]
    spatialShape = geo["X"].shape
    dateRange = len(CMAQTime)
    result = np.zeros((dateRange, spatialShape[0], spatialShape[1]))
    for i in range(0, dateRange):
        cur_time = CMAQTime[i]
        if cur_time in OBSTime:
            obs_time_idx = OBSTime.index(cur_time)
            # Find each day's valid site
            obsXDaily, obsYDaily = obsX.copy(), obsY.copy()
            valid_idx = ~np.isnan(obsConc[obs_time_idx, :])
            obsXDaily, obsYDaily = obsXDaily[valid_idx], obsYDaily[valid_idx]
            distance = np.sqrt(np.power((obsXDaily - geo["X"][:, :, np.newaxis]), 2) + np.power((obsYDaily - geo["Y"][:, :, np.newaxis]), 2))
            minDistance = np.min(distance, axis=2) / 1000
            result[i, :, :] = Rcoll * np.exp(-minDistance/r)
        else:
            result[i, :, :] = 0
    return result

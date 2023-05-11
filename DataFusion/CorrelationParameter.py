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
            
            # we cannot calculate correlation if data less than 2
            if OBS1.shape[0] <= 1 or OBS2.shape[0] <= 1:
                continue
            # we cannot calculate correlation if data is a constant
            variance_obs1 = (OBS1 - np.mean(OBS1)) / np.mean(OBS1)
            variance_obs2 = (OBS2 - np.mean(OBS2)) / np.mean(OBS2)
            denominator_pearson_r = np.sum(variance_obs1 ** 2) * np.sum(variance_obs2 ** 2)
            if denominator_pearson_r <= 1e-20:
                continue
            
            pearson_r = scipy.stats.pearsonr(OBS1, OBS2)
            corrCoef.append(pearson_r[0])
            # calculate the distance change meter to kilometer
            X1 = X[i]/1000
            Y1 = Y[i]/1000
            X2 = X[j]/1000
            Y2 = Y[j]/1000
            d = np.sqrt((X1 - X2) * (X1 - X2) + (Y1 - Y2) * (Y1 - Y2))
            d = np.asscalar(d)
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


def R1(obsDict, geo, Rcoll, r):
    obsConc = obsDict["Conc"]
    obsX = obsDict["X"]
    obsY = obsDict["Y"]
    spatialShape = geo["X"].shape
    dateRange = obsConc.shape[0]
    siteSize = len(obsX)
    result = np.empty((dateRange, spatialShape[0], spatialShape[1]))
    result[:] = np.NaN
    for i in range(0, dateRange):
        # Find each day's valid site
        invalidSite = np.isnan(obsConc[i, :])
        obsXDaily = obsX.copy()
        obsYDaily = obsY.copy()
        obsXDaily[invalidSite] = np.NaN
        obsYDaily[invalidSite] = np.NaN
        # Find the nearest site for each grids
        for j in range(0, spatialShape[0]):
            for k in range(0, spatialShape[1]):
                xTmp = geo["X"][j, k]
                yTmp = geo["Y"][j, k]
                distance = np.sqrt(np.power((obsXDaily - xTmp), 2) + np.power((obsYDaily - yTmp), 2))
                # calculate the distance change meter to kilometer
                if np.isnan(distance).all():
                    result[i, j, k] = 0
                else:
                    minDistance = np.nanmin(distance)/1000
                    result[i, j, k] = Rcoll * np.exp(-minDistance/r)
    return result

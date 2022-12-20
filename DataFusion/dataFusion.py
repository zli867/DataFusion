import numpy as np
from DataFusion.KrigingOBS import krigingOBS
from scipy.optimize import curve_fit


def dataFusionOne(CMAQDict, OBSDict, geoDict, alpha_yearly, beta):
    # Use CMAQAdjust function to calculate adjusted CMAQ
    FCData = FC(CMAQDict["Yearly"], alpha_yearly, beta)

    result = np.empty(CMAQDict["Daily"].shape)
    result[:] = np.NaN
    dateRange = result.shape[0]
    # Use kriging to calculate ratio
    OBSConc = OBSDict["Conc"]
    OBSConcMean = np.nanmean(OBSConc, axis=0)
    obsX = OBSDict["X"]
    obsY = OBSDict["Y"]
    predictX = geoDict["X"]
    predictY = geoDict["Y"]
    for i in range(0, dateRange):
        OBSNormalized = OBSConc[i, :] / OBSConcMean
        # Krig method to interpolate the observation
        # delete the Nan number in array
        nanIndex = np.isnan(OBSNormalized)
        OBSNormalized = OBSNormalized[~nanIndex]
        siteX = obsX[~nanIndex]
        siteY = obsY[~nanIndex]
        # If there is less than one observation, do not do kriging
        # if the data variance is zero (constant), do not do kriging
        valid_krig = True
        if OBSNormalized.shape[0] == 0 or OBSNormalized.shape[0] == 1:
            valid_krig = False
        else:
            variance_data = np.var(OBSNormalized)
            if variance_data == 0:
                valid_krig = False
        if valid_krig:
            krigRatio = krigingOBS(siteX, siteY, OBSNormalized, predictX, predictY)
            resultDaily = krigRatio * FCData
        else:
            resultDaily = FCData
        result[i, :, :] = resultDaily
    return result


def temporal_func(t, A, t_max):
    return np.exp(A * np.cos(((2 * np.pi)/365.25) * (t - t_max)))


def is_leap_year(year):
    result = False
    if (year % 4 == 0 and year % 100 != 0) or (year % 400 == 0):
        result = True
    return result


def time_series_generator(time_list):
    results = []
    for current_time in time_list:
        tt = current_time.timetuple()
        results.append(tt.tm_yday)
    return np.array(results)


def temporalCorrection(obs_adjusted_CMAQ_ratio, obs_time_series):
    time_series = time_series_generator(obs_time_series)
    popt, pcov = curve_fit(temporal_func, time_series, obs_adjusted_CMAQ_ratio, p0=(0, 180))
    # plt.scatter(time_series, obs_adjusted_CMAQ_ratio)
    # plt.show()
    A = popt[0]
    t_max = popt[1]
    return A, t_max


def dataFusionTwo(CMAQDict, alpha_yearly, beta, A, tmax):
    # If there is no saesonal correction, A = 0
    FCData = FC(CMAQDict["Yearly"], alpha_yearly, beta) / CMAQDict["Yearly"]
    CMAQDaily = CMAQDict["Daily"]
    julian_dates = time_series_generator(CMAQDict["Time"])
    correction_ratio = temporal_func(julian_dates, A, tmax)
    result = CMAQDaily * FCData * correction_ratio[:, np.newaxis, np.newaxis]
    return result


def FC(CMAQConcYearly, alpha_yearly, beta):
    return alpha_yearly * np.power(CMAQConcYearly, beta)


def weightFactor(R1, R2):
    return (R1 * (1 - R2)) / (R1 * (1 - R2) + R2 * (1 - R1))

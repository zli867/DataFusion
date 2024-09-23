from sklearn import linear_model
import numpy as np
from scipy.optimize import curve_fit
from Evaluation.StatisticalMetrics import RMSE


def powerFunc(x, alpha, beta):
    return alpha * (x ** beta)


def adjustParameters(CMAQatObsConc, obsMeanConc):
    """
    obs = alpha*CMAQ^beta => log(obs) = log(alpha) + beta*log(CMAQ)
    :param CMAQatObsConc: CMAQ values at obs location
    :param obsMeanConc: observed values
    :return: alpha, beta and performance r2 score
    """
    if 0 not in CMAQatObsConc and 0 not in obsMeanConc:
        log_CMAQatOBSConc = np.log(CMAQatObsConc)
        log_obs = np.log(obsMeanConc)
        log_CMAQatOBSConc = log_CMAQatOBSConc.reshape(-1, 1)
        log_obs = log_obs.reshape(-1, 1)
        reg = linear_model.LinearRegression(fit_intercept=True).fit(log_CMAQatOBSConc, log_obs)
        beta = reg.coef_[0][0]
        log_alpha = reg.intercept_
        alpha = np.exp(log_alpha)
    else:
        # previous version, use scipy do the optimization
        popt, pcov = curve_fit(powerFunc, CMAQatObsConc, obsMeanConc, p0=(1, 1))
        alpha = popt[0]
        beta = popt[1]

    adjusted_CMAQ = adjustCMAQ(CMAQatObsConc, alpha, beta)
    score = RMSE(adjusted_CMAQ, obsMeanConc)
    return alpha, beta, score

# OBS = alpha * CMAQ^beta => regress for OBS and CMAQ^beta
def adjustParamentersYearly(beta, CMAQatObsConc, obsMeanConc):
    CMAQ_beta = np.power(CMAQatObsConc, beta)
    CMAQ_beta = CMAQ_beta.reshape(-1, 1)
    obsMeanConc = obsMeanConc.reshape(-1, 1)
    reg = linear_model.LinearRegression(fit_intercept=False).fit(CMAQ_beta, obsMeanConc)
    alpha = reg.coef_[0][0]
    return alpha


def adjustParamentersLinear(CMAQatObsConc, obsMeanConc):
    CMAQ_beta = CMAQatObsConc.reshape(-1, 1)
    obsMeanConc = obsMeanConc.reshape(-1, 1)
    reg = linear_model.LinearRegression(fit_intercept=False).fit(CMAQ_beta, obsMeanConc)
    alpha = reg.coef_[0][0]
    # evaluate performance by R2
    adjusted_CMAQ = reg.predict(CMAQ_beta)
    score = RMSE(adjusted_CMAQ, obsMeanConc)
    return alpha, 1.0, score


def adjustCMAQ(CMAQYearlyConc, alpha, beta):
    result = alpha * np.power(CMAQYearlyConc, beta)
    return result
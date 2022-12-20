import numpy as np


def CMAQatOBSDaily(dictOBS, dictGeo, dictCMAQ):
    """

    :param dictOBS: dictionary of observation
    :param dictGeo: static geographic definition of CMAQ
    :param dictCMAQ: dictionary of CMAQ
    :return: a concentration matrix (m, n) where m is the size of date range for observation, n is the number of sites.
             The concentration denotes the CMAQ results at the observation sites. If the value is nan at (i, j), it
             means there is no observation at i date and j site.
    """
    xCMAQ = dictGeo["X"]
    yCMAQ = dictGeo["Y"]
    siteCode = dictOBS["siteCode"]
    siteSize = len(siteCode)
    CMAQ_time = dictCMAQ["Time"]
    obs_time = dictOBS["dateSeries"]
    dateSize = len(obs_time)
    xOBS = dictOBS["X"]
    yOBS = dictOBS["Y"]

    geoShape = xCMAQ.shape
    row = geoShape[0]
    col = geoShape[1]
    ConcDaily = np.empty((dateSize, siteSize))
    ConcDaily[:] = np.NaN

    # calculate the corresponding position in CMAQ
    positionCMAQ = []
    for i in range(0, siteSize):
        siteX = xOBS[i]
        siteY = yOBS[i]
        deltaX = xCMAQ - siteX
        deltaY = yCMAQ - siteY
        distanceSquare = np.power(deltaX, 2) + np.power(deltaY, 2)
        pos = np.argmin(distanceSquare)
        x = int(pos / col)
        y = pos % col
        coordinator = (x, y)
        positionCMAQ.append(coordinator)

    for i in range(0, dateSize):
        current_time = obs_time[i]
        # Notice: we just extract the observation when we have CMAQ data, observation time series is a subset of CMAQ
        # time series.
        CMAQ_time_idx = CMAQ_time.index(current_time)
        for j in range(0, siteSize):
            x, y = positionCMAQ[j]
            if not np.isnan(dictOBS["Conc"][i, j]):
                ConcDaily[i, j] = dictCMAQ["Daily"][CMAQ_time_idx, x, y]
    return ConcDaily
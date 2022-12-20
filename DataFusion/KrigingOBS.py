import pykrige
import numpy as np


def krigingOBS(obsX, obsY, obsConc, predictX, predictY):
    model = pykrige.ok.OrdinaryKriging(obsX, obsY, obsConc, variogram_model='exponential')
    geoSize = predictX.shape
    X = np.reshape(predictX, predictX.size)
    Y = np.reshape(predictY, predictY.size)
    [h, t] = model.execute('points', X, Y)
    h = np.reshape(h, (geoSize[0], geoSize[1]))
    return h
import numpy as np


def timeArrayToTFLAG(time_array, dimension):
    res = []
    for current_time in time_array:
        current_time_str = current_time.strftime("%Y%j")
        current_time_date = int(current_time_str)
        current_time_hour = current_time.hour
        res.append([current_time_date, current_time_hour])
    res = np.array(res)
    res = res[:, np.newaxis, :]
    res = np.repeat(res, dimension, axis=1)
    return res

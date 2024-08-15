import numpy as np


def HR_MAX(vals):
    return np.max(vals, axis=0)


def AVG(vals):
    return np.mean(vals, axis=0)


def MDA8(vals):
    candidates = []
    for i in range(0, len(vals) - 8):
        cur_mean = np.mean(vals[i: i + 8, :, :], axis=0)
        candidates.append(cur_mean[np.newaxis, :, :])
    combined_val = np.concatenate(candidates, axis=0)
    max_val = np.max(combined_val, axis=0)
    return max_val
import pandas as pd
import numpy as np


def standardization(matrix):
    matrix_sum = np.sum(matrix, axis=0)
    matrix = matrix / (matrix_sum[np.newaxis, :]+1e-16)
    return matrix


def unknown_initialize(sources, sinks):
    sinks = sinks.reshape((sinks.shape[0],))
    unknown_pre = sinks - np.sum(sources, axis=1)
    unknown_zero = np.zeros((unknown_pre.shape[0],))
    unknown_init = np.vstack((unknown_pre, unknown_zero))
    unknown_init = np.max(unknown_init, axis=0)+1e-16
    unknown_init = unknown_init.reshape((-1, 1))
    unknown_init = standardization(unknown_init)
    return unknown_init


def parameters_initialize(sources, sinks, unknown_init, A):
    y = sources
    n = sources.shape[0]
    k = sources.shape[1]
    x = sinks  # only 1 need to change
    w = y.copy()
    h = np.zeros((k + 1, 1)) + 1 / (k + 1)
    a = np.ones((n, k)) * A  # previous K columns in the weighted matrix are one
    zeros = np.zeros((n, 1), dtype=y.dtype)
    a = np.hstack((a, zeros))  # while (K+1)-th column has zero values
    y = np.hstack((y, zeros))
    w = np.hstack((w, unknown_init))
    w_plus = w.copy()
    h_plus = h.copy()
    alpha_w = np.zeros((n, k + 1))
    alpha_h = np.zeros((k + 1, 1))
    i = np.identity(k + 1)  # Let I be the identity matrix
    return x, y, w, a, i, w_plus, h_plus, alpha_w, alpha_h






import sys
import numpy as np
import concurrent.futures
import multiprocessing as mp
from numpy.linalg import inv
from data_plot import iteration_performance
from scipy.optimize import root_scalar
from tqdm import tqdm
# mp.set_start_method("fork")
# thread = int(sys.argv[1])


def w_update(w_plus, alpha_w, h, x, i, rho):  # w_renew = w_update()
    return (x.dot(h.T) + rho * w_plus - alpha_w).dot(inv(h.dot(h.T) + rho * i))


def h_update(w_renew, x, i, h_plus, alpha_h, rho):  # h_renew = h_update()
    return inv(w_renew.T.dot(w_renew) + rho * i).dot(w_renew.T.dot(x) + rho * h_plus - alpha_h)


def w_add(w_renew, d):  # d = np.diag(w_renew)
    return w_renew.dot(inv(d))  # w_latest = w_add()


def h_add(h_renew, d):  # d = np.diag(w_renew)
    return d.dot(h_renew)  # h_latest = h_add()


def water_filling_w(beta_w, w, alpha_w, y, a, rho):
    w1 = np.power(a, 2) * y + alpha_w + rho * w - beta_w
    w2 = np.power(a, 2) + rho
    return np.sum(np.maximum(0, w1 / w2)) - 1


def water_filling_h(beta_h, h, alpha_h, rho):
    return np.sum(np.maximum(0, (alpha_h + rho * h - beta_h) / rho)) - 1


def root_scalar_w(w, alpha_w, y, a, rho):
    sol = root_scalar(water_filling_w, args=(w, alpha_w, y, a, rho), method='bisect', bracket=[-10000, 10000])  #
    return sol.root


def root_scalar_h(h, alpha_h, rho):
    sol = root_scalar(water_filling_h, args=(h, alpha_h, rho), method='bisect', bracket=[-10000, 10000])  #
    return sol.root


def call_beta(args):
    i = args[0]
    root = root_scalar_w(args[1][:, i], args[2][:, i], args[3][:, i], args[4][:, i], args[5])
    return root


def w_plus_update(a, y, w_latest, alpha_w, rho):
    thread = 16
    query_list = [(i, w_latest, alpha_w, y, a, rho) for i in range(w_latest.shape[1])]
    # query_list = [i for i in range(w_latest.shape[1])]
    with concurrent.futures.ProcessPoolExecutor(thread) as executor:
        beta = list(tqdm(executor.map(call_beta, query_list), total=len(query_list)))
    w1 = np.power(a, 2) * y + alpha_w + rho * w_latest - beta
    w2 = np.power(a, 2) + rho
    w_plus = np.maximum(0, w1 / w2)
    return w_plus
    # return np.maximum((np.multiply(np.multiply(a, y), a) + alpha_w + rho * w_latest) / (np.multiply(a, a) + rho), 0)


def h_plus_update(h_latest, alpha_h, rho):
    root = root_scalar_h(h_latest, alpha_h, rho)
    h_plus = np.maximum(0, (alpha_h + rho * h_latest - root) / rho)
    return h_plus
    # return np.maximum(h_latest + alpha_h / rho, 0)


def alpha_w_update(alpha_w, w_latest, w_plus, rho):
    return alpha_w + rho * (w_latest - w_plus)


def alpha_h_update(alpha_h, h_latest, h_plus, rho):
    return alpha_h + rho * (h_latest - h_plus)


def standardization(matrix):
    matrix_sum = np.sum(matrix, axis=0)
    matrix = matrix / matrix_sum[np.newaxis, :]
    return matrix


def simluated_data(noise_data, sink_data, unknown_init, rho):
    y = noise_data
    n = noise_data.shape[0]
    k = noise_data.shape[1]
    x = sink_data  # only 1 need to change
    w = y.copy()
    h = np.zeros((k + 1, 1)) + 1 / (k + 1)
    a = np.ones((n, k))  # previous K columns in the weighted matrix are one, while (K+1)-th column has zero values
    zeros = np.zeros((n, 1), dtype=y.dtype)
    # random = np.zeros((n, 1), dtype=y.dtype) + 1e-05
    unknown_init = standardization(unknown_init)
    a = np.hstack((a, zeros))
    y = np.hstack((y, zeros))
    w = np.hstack((w, unknown_init))
    w_plus = w.copy()
    h_plus = h.copy()
    alpha_w = np.zeros((n, k + 1))
    alpha_h = np.zeros((k + 1, 1))
    i = np.identity(k + 1)  # Let I be the identity matrix
    return x, y, w, h, a, i, w_plus, h_plus, alpha_w, alpha_h, rho


def admm_step(w, x, y, a, i, w_plus, h_plus, alpha_w, alpha_h, rho):
    h_renew = h_update(w, x, i, h_plus, alpha_h, rho)
    w_renew = w_update(w_plus, alpha_w, h_renew, x, i, rho)
    d = np.diag(np.sum(w_renew, axis=0))
    h_renew = h_add(h_renew, d)
    w_renew = w_add(w_renew, d)
    h_plus = h_plus_update(h_renew, alpha_h, rho)
    w_plus = w_plus_update(a, y, w_renew, alpha_w, rho)
    alpha_w = alpha_w_update(alpha_w, w_renew, w_plus, rho)
    alpha_h = alpha_h_update(alpha_h, h_renew, h_plus, rho)
    return w_renew, h_renew, w_plus, h_plus, alpha_w, alpha_h


def lagrangian(x, y, w, h, a, w_plus, h_plus, alpha_w, alpha_h, rho):
    l1 = 0.5 * np.linalg.norm(x - w.dot(h))
    l2 = 0.5 * np.linalg.norm(np.multiply(a, (w_plus - y)))
    l3 = np.sum(np.diag(alpha_w.T.dot(w - w_plus)))
    l4 = (rho / 2) * np.linalg.norm(w - w_plus)
    l5 = np.sum(np.diag(alpha_h.T.dot(h - h_plus)))
    l6 = (rho / 2) * np.linalg.norm(h - h_plus)
    return l1 + l2 + l3 + l4 + l5 + l6


def run_admm(iteration, sink_abundance, noise_data, sink_data, original_data, unknown_init, rho, value):
    x, y, w, h, a, i, w_plus, h_plus, alpha_w, alpha_h, rho = simluated_data(noise_data, sink_data, unknown_init, rho)
    loss = []
    diff_xwh = []
    jsd = []
    diff_wo = []
    diff_un = []
    w_renew, h_renew, w_plus, h_plus, alpha_w, alpha_h = admm_step(w, x, y, a, i, w_plus, h_plus, alpha_w, alpha_h, rho)
    l_init = lagrangian(x, y, w_renew, h_renew, a, w_plus, h_plus, alpha_w, alpha_h, rho)
    loss.append(l_init)
    for m in range(1, iteration):
        w_renew, h_renew, w_plus, h_plus, alpha_w, alpha_h = admm_step(w_renew, x, y, a, i, w_plus, h_plus, alpha_w, alpha_h, rho)
        l_update = lagrangian(x, y, w_renew, h_renew, a, w_plus, h_plus, alpha_w, alpha_h, rho)
        diff_xwh_iters, jsd_iters, diff_worigin_iters, diff_un_iters = iteration_performance(x, w_plus, h_plus, sink_abundance, original_data)
        loss.append(l_update)
        diff_xwh.append(diff_xwh_iters)
        jsd.append(jsd_iters)
        diff_wo.append(diff_worigin_iters)
        diff_un.append(diff_un_iters)
        if m == 1 and abs(l_update - l_init) / l_init < value:
            break
        if m > 1 and abs(l_update - loss[m - 1]) / loss[m - 1] < value:
            break
        else:
            continue
    h_plus = np.reshape(h_plus, (h_plus.shape[0]))
    loss_value = abs(loss[-1] - loss[-2]) / loss[-1]
    return h_plus, w_plus, w_renew, h_renew, m, loss, loss_value, diff_xwh, jsd, diff_wo, diff_un
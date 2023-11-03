import numpy as np
import concurrent.futures
import multiprocessing as mp
from tqdm import tqdm
from numpy.linalg import inv
from data_initialization import standardization
from scipy.optimize import root_scalar
mp.set_start_method("fork")


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


def w_plus_update(a, y, w_latest, alpha_w, rho, thread):
    query_list = [(i, w_latest, alpha_w, y, a, rho) for i in range(w_latest.shape[1])]
    with concurrent.futures.ProcessPoolExecutor(thread) as executor:
        beta = list(tqdm(executor.map(call_beta, query_list), total=len(query_list)))
    w1 = np.power(a, 2) * y + alpha_w + rho * w_latest - beta
    w2 = np.power(a, 2) + rho
    w_plus = np.maximum(0, w1 / w2)
    return w_plus


def h_plus_update(h_latest, alpha_h, rho):
    root = root_scalar_h(h_latest, alpha_h, rho)
    h_plus = np.maximum(0, (alpha_h + rho * h_latest - root) / rho)
    return h_plus


def alpha_w_update(alpha_w, w_latest, w_plus, rho):
    return alpha_w + rho * (w_latest - w_plus)


def alpha_h_update(alpha_h, h_latest, h_plus, rho):
    return alpha_h + rho * (h_latest - h_plus)


def admm_step(w, x, y, a, i, w_plus, h_plus, alpha_w, alpha_h, rho, thread):
    h_renew = h_update(w, x, i, h_plus, alpha_h, rho)
    w_renew = w_update(w_plus, alpha_w, h_renew, x, i, rho)
    d = np.diag(np.sum(w_renew, axis=0))
    h_renew = h_add(h_renew, d)
    w_renew = w_add(w_renew, d)
    h_plus = h_plus_update(h_renew, alpha_h, rho)
    w_plus = w_plus_update(a, y, w_renew, alpha_w, rho, thread)
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




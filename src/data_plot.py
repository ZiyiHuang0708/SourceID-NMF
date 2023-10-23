import numpy as np
from scipy.spatial.distance import jensenshannon
import matplotlib.pyplot as plt
import pylab as pl
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


def iteration_performance(x, w_plus, h_plus, sink_abundance, original_data):
    h_plus = np.reshape(h_plus, (h_plus.shape[0],))
    x = np.reshape(x, (x.shape[0],))
    norm_xwh_iters = np.linalg.norm(x - np.dot(w_plus, h_plus))
    diff_xwh_iters = np.sum(abs(x - np.dot(w_plus, h_plus)))
    jsd_iters = jensenshannon(sink_abundance, h_plus)
    diff_worigin_iters = np.sum(abs(w_plus - original_data))
    diff_un_iters = np.sum(abs(w_plus[:, w_plus.shape[1]-1] - original_data[:, original_data.shape[1]-1]))
    return diff_xwh_iters, jsd_iters, diff_worigin_iters, diff_un_iters


def plot_loss_function(loss1, loss2, loss3, loss4, path):
    iters1 = list(range(1, len(loss1) + 1))
    iters2 = list(range(1, len(loss2) + 1))
    iters3 = list(range(1, len(loss3) + 1))
    iters4 = list(range(1, len(loss4) + 1))
    plt.figure(1)
    ax1 = plt.figure(9).add_subplot(1, 1, 1)
    plt.plot(iters1, loss1, 'r-', label='sink1 loss function', color='red', linewidth=0.75)
    plt.legend()
    plt.plot(iters2, loss2, 'r-', label='sink2 loss function', color='deepskyblue', linewidth=0.75)
    plt.legend()
    plt.plot(iters3, loss3, 'r-', label='sink3 loss function', color='green', linewidth=0.75)
    plt.legend()
    plt.plot(iters4, loss4, 'r-', label='sink4 loss function', color='yellow', linewidth=0.75)
    plt.legend()
    plt.title('Compare loss for models in training')
    plt.xlabel('iters')
    plt.ylabel('loss')
    tx0 = 0
    tx1 = 50
    ty0 = 0
    ty1 = 0.1
    sx = [tx0, tx1, tx1, tx0, tx0]
    sy = [ty0, ty0, ty1, ty1, ty0]
    plt.plot(sx, sy, "purple", linewidth=0.75)
    axins1 = inset_axes(ax1, width=1.5, height=1.5, loc='right')
    axins1.plot(iters1, loss1, 'r-', color='red', linewidth=0.75)
    axins1.plot(iters2, loss2, 'r-', color='deepskyblue', linewidth=0.75)
    axins1.plot(iters3, loss3, 'r-', color='green', linewidth=0.75)
    axins1.plot(iters4, loss4, 'r-', color='yellow', linewidth=0.75)
    axins1.axis([tx0, tx1, ty0, ty1])
    plt.savefig(path)


def plot_jsd_iteration(jsd1, jsd2, jsd3, jsd4, path):
    iters1 = list(range(1, len(jsd1) + 1))
    iters2 = list(range(1, len(jsd2) + 1))
    iters3 = list(range(1, len(jsd3) + 1))
    iters4 = list(range(1, len(jsd4) + 1))
    plt.figure(2)
    plt.plot(iters1, jsd1, 'r-', label='sink1 jsd_iteration', color='red', linewidth=0.75)
    plt.plot(iters2, jsd2, 'r-', label='sink2 jsd_iteration', color='deepskyblue', linewidth=0.75)
    plt.plot(iters3, jsd3, 'r-', label='sink3 jsd_iteration', color='green', linewidth=0.75)
    plt.plot(iters4, jsd4, 'r-', label='sink4 jsd_iteration', color='yellow', linewidth=0.75)
    plt.legend()
    plt.xlabel('iters')
    plt.ylabel('jsd_iteration')
    plt.title('Compare jsd_iteration for models in training')
    plt.savefig(path)


def plot_diff_xwh_iteration(xwh1, xwh2, xwh3, xwh4, path):
    iters1 = list(range(1, len(xwh1) + 1))
    iters2 = list(range(1, len(xwh2) + 1))
    iters3 = list(range(1, len(xwh3) + 1))
    iters4 = list(range(1, len(xwh4) + 1))
    plt.figure(3)
    plt.plot(iters1, xwh1, 'r-', label='sink1 diff_xwh_iteration', color='red', linewidth=0.75)
    plt.plot(iters2, xwh2, 'r-', label='sink2 diff_xwh_iteration', color='deepskyblue', linewidth=0.75)
    plt.plot(iters3, xwh3, 'r-', label='sink3 diff_xwh_iteration', color='green', linewidth=0.75)
    plt.plot(iters4, xwh4, 'r-', label='sink4 diff_xwh_iteration', color='yellow', linewidth=0.75)
    plt.legend()
    plt.xlabel('iters')
    plt.ylabel('diff_xwh_iteration')
    plt.title('Compare diff_xwh_iteration for models in training')
    plt.savefig(path)


def plot_w_origin_iteration(wo1, wo2, wo3, wo4, path):
    iters1 = list(range(1, len(wo1) + 1))
    iters2 = list(range(1, len(wo2) + 1))
    iters3 = list(range(1, len(wo3) + 1))
    iters4 = list(range(1, len(wo4) + 1))
    plt.figure(4)
    plt.plot(iters1, wo1, 'r-', label='sink1 w_origin_iteration', color='red', linewidth=0.75)
    plt.plot(iters2, wo2, 'r-', label='sink2 w_origin_iteration', color='deepskyblue', linewidth=0.75)
    plt.plot(iters3, wo3, 'r-', label='sink3 w_origin_iteration', color='green', linewidth=0.75)
    plt.plot(iters4, wo4, 'r-', label='sink4 w_origin_iteration', color='yellow', linewidth=0.75)
    plt.legend()
    plt.xlabel('iters')
    plt.ylabel('w_origin_iteration')
    plt.title('Compare w_origin_iteration for models in training')
    plt.savefig(path)


def plot_unknown_iteration(unknown1, unknown2, unknown3, unknown4, path):
    iters1 = list(range(1, len(unknown1) + 1))
    iters2 = list(range(1, len(unknown2) + 1))
    iters3 = list(range(1, len(unknown3) + 1))
    iters4 = list(range(1, len(unknown4) + 1))
    plt.figure(5)
    plt.plot(iters1, unknown1, 'r-', label='sink1 unknown_iteration', color='red', linewidth=0.75)
    plt.plot(iters2, unknown2, 'r-', label='sink2 unknown_iteration', color='deepskyblue', linewidth=0.75)
    plt.plot(iters3, unknown3, 'r-', label='sink3 unknown_iteration', color='green', linewidth=0.75)
    plt.plot(iters4, unknown4, 'r-', label='sink4 unknown_iteration', color='yellow', linewidth=0.75)
    plt.legend()
    plt.xlabel('iters')
    plt.ylabel('unknown_iteration')
    plt.title('Compare unknown_iteration for models in training')
    plt.savefig(path)









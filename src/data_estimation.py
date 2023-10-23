import numpy as np
from scipy.stats import pearsonr
from scipy.spatial.distance import jensenshannon
import math
# returns an S1 by S2 matrix P, where S1 is the number sinks and S2 is the number of sources
# Each row in matrix P sums to 1. Pij is the contribution of source j to sink i
# estimated_proportions = ...
# true_proportions = ...


def source_jsd_estimation(source_data):
    jsd_ave = []
    for i in range(source_data.shape[1] - 1):
        for j in range(i + 1, source_data.shape[1]):
            jsd = jensenshannon(source_data[:, i], source_data[:, j])
            jsd_ave.append(jsd)
    jensen_shannon_mean = np.mean(jsd_ave)
    return jensen_shannon_mean


def noise_source_jsd_estimation(noise_data, source_data):
    jensen_shannon = []
    for i in range(noise_data.shape[1]):
        noise_pre = noise_data[:, i]
        source_pre = source_data[:, i]
        j = jensenshannon(noise_pre, source_pre)
        jensen_shannon.append(j)
    jensen_shannon_mean = np.mean(jensen_shannon)
    return jensen_shannon_mean


# Calculate the squared Pearson correlation (r2) between the estimated and the true mixing proportions per source
# and average across sources
def pearson_correlation(estimated_proportions, true_proportions):
    pearson = []
    for i in range(estimated_proportions.shape[1]):
        estimated_data = estimated_proportions[:, i]
        true_data = true_proportions[:, i]
        print(estimated_data.shape)
        print(true_data.shape)
        r2, p_value = pearsonr(estimated_data, true_data)
        pearson.append(abs(r2))
    pearson_mean = np.mean(pearson)
    return pearson_mean, pearson


# Calculate the average Jensen-Shannon divergence of m (based on the pairwise Jensen-Shannon divergence)
def jensen_shannon_divergence(estimated_proportions, true_proportions):
    jensen_shannon = []
    for i in range(estimated_proportions.shape[0]):
        mixing_data = estimated_proportions[i, :]
        true_data = true_proportions[i, :]
        j = jensenshannon(mixing_data, true_data)
        jensen_shannon.append(j)
    jensen_shannon_mean = np.mean(jensen_shannon)
    return jensen_shannon_mean, jensen_shannon


def perf_admm(estimated_proportions, sink_abundance):
    # pearson, _ = pearson_correlation(sink_abundance, estimated_proportions)
    jensen_shannon_mean, jensen_shannon = jensen_shannon_divergence(sink_abundance, estimated_proportions)
    difference = np.sum(abs(sink_abundance - estimated_proportions), axis=1)
    return jensen_shannon_mean, jensen_shannon, difference


def diff_admm(sink_data, original_data, w_plus, w, h):
    x = sink_data
    num = original_data.shape[1] - 1
    unknown = original_data[:, num]
    diff_xwh = np.sum(abs(x - np.dot(w, h)))
    norm_xwh = np.linalg.norm(x - np.dot(w, h))
    diff_unknown = np.sum(abs(w_plus[:, num] - unknown))
    jsd_unknown = jensenshannon(w_plus[:, num], unknown)
    return diff_xwh, norm_xwh, diff_unknown, jsd_unknown


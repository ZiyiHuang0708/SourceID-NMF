import argparse
import numpy as np
import pandas as pd
from scipy.spatial.distance import jensenshannon
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


def jsd_estimation(noise_data, source_data):
    jensen_shannon = []
    for i in range(noise_data.shape[1]):
        noise_pre = noise_data[:, i]
        source_pre = source_data[:, i]
        j = jensenshannon(noise_pre, source_pre)
        jensen_shannon.append(j)
    jensen_shannon_mean = np.mean(jensen_shannon)
    return jensen_shannon_mean


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


def perf_admm(estimated_proportions, true_proportions):
    jensen_shannon_mean, jensen_shannon = jensen_shannon_divergence(estimated_proportions, true_proportions)
    difference = np.sum(abs(estimated_proportions - true_proportions), axis=1)
    return jensen_shannon_mean, jensen_shannon, np.mean(difference), difference


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Evaluation of SourceID-NMF Estimated Proportions.""")
    parser.add_argument('--estimated', '-e', help='input the estimated proportions', default='./estimated_proportions.txt')
    parser.add_argument('--true', '-t', help='input the true proportions', default='./true_proportions.txt')
    inputs = parser.parse_args()

    estimated_proportions_path = inputs.estimated
    true_proportions_path      = inputs.true
    estimated_proportions = pd.read_csv(estimated_proportions_path, sep=" ", header=0, index_col=0)
    true_proportions = pd.read_csv(true_proportions_path, sep="	", header=0, index_col=0)
    row_labels = list(estimated_proportions.index.values)
    estimated_proportions = np.array(estimated_proportions)
    true_proportions = np.array(true_proportions)

    perf = []
    jsd_ave, jsd, diff_ave, diff = perf_admm(estimated_proportions, true_proportions)
    print("The average JSD between estimated and true proportion is", jsd_ave)
    print("The average Difference between estimated and true proportion is", diff_ave)
    perf.append(jsd)
    perf.append(diff)
    perf = np.array(perf)
    nmf_perf_label = ['jsd', 'difference']
    perf = pd.DataFrame(perf, index=nmf_perf_label, columns=row_labels)
    perf.to_csv('./proportion_perf.txt', sep=' ')




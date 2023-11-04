import argparse
import numpy as np
import pandas as pd


def sink_proportion_output(sources, unknown_proportion):
    sink_abundance = []
    unknown = unknown_proportion
    m = np.random.pareto(2, sources[:, :-1].shape[1])  # Add the mix ratio of unknown data
    m = m / np.sum(m) * (1 - unknown)
    m = np.append(m, unknown)
    for j in range(sources.shape[0]):
        abundance = np.dot(sources[j, :], m)
        abundance = np.round(abundance)
        sink_abundance.append(abundance)
    return sink_abundance, m


def sink_data_output(sources, unknown_proportion, value):
    sink_abundances = []
    true_proportions = []
    noise_data = []
    for i in range(len(unknown_proportion)):
        sink_abundance, m = sink_proportion_output(sources, unknown_proportion[i])
        sink_abundances.append(sink_abundance)
        true_proportions.append(m)
    sink_abundances = np.array(sink_abundances)
    sink_abundances = sink_abundances.T
    true_proportions = np.array(true_proportions)
    for i in range(sources.shape[1] - 1):
        count_list = list(sources[:, i])
        probs = count_list / np.sum(count_list)
        multinomial_data = np.random.multinomial(value, probs)
        noise_data.append(multinomial_data)
    noise_data = np.array(noise_data).T
    nmf_data = np.hstack((noise_data, sink_abundances))
    return nmf_data, true_proportions


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Generate simulated data.""")
    parser.add_argument('--input', '-i', help='input the original data', default='./0.70jsd.txt')
    parser.add_argument('--sample', '-s', help='Number of samples for polynomial distribution', type=int, default=100000)
    inputs = parser.parse_args()

    sources_path = inputs.input
    sample_path = inputs.sample
    sources = np.array(pd.read_csv(sources_path, sep="	", header=0))
    sources_num = sources.shape[1]
    unknown_proportions = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    sinks_num = len(unknown_proportions)
    nmf_data, true_proportions = sink_data_output(sources, unknown_proportions, sample_path)

    taxon_list = []
    for i in range(nmf_data.shape[0]):
        taxon_list.append("taxon_" + str(i))
    sample_list = []
    for j in range(1, nmf_data.shape[1]+1):
        sample_list.append("D" + str(j))
    source_list = []
    for k in range(1, sources_num):
        source_list.append("D" + str(k))
    source_list.append("unknown")
    sink_list = []
    for n in range(sources_num, sources_num+sinks_num):
        sink_list.append("D" + str(n))

    nmf_data = pd.DataFrame(nmf_data, index=taxon_list, columns=sample_list)
    true_proportions = pd.DataFrame(true_proportions, index=sink_list, columns=source_list)
    nmf_data.to_csv('./nmf_data.txt', sep=' ')
    true_proportions.to_csv('./true_proportions.txt', sep=' ')


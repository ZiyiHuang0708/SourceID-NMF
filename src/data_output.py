import numpy as np
import pandas as pd


def sink_proportion_output(source_data, i, proportion_path):
    sink_abundance = []
    proportion = np.array(pd.read_excel(proportion_path, header=0, sheet_name='true_proportions'))
    for j in range(source_data.shape[0]):
        abundance = np.dot(source_data[j, :], proportion[i, :])
        sink_abundance.append(abundance)
    return sink_abundance


def sink_data_output(project_data, sink_number, proportion_path, value):
    sink_abundances = []
    probs_sum = []
    noise_data = []

    for i in range(sink_number):
        sink_abundance = sink_proportion_output(project_data, i, proportion_path)
        sink_abundances.append(sink_abundance)
    sink_abundances = np.array(sink_abundances)
    sink_abundances = sink_abundances.T

    for i in range(project_data.shape[1]):
        count_list = list(project_data[:, i])
        probs = count_list / np.sum(count_list)
        multinomial_data = np.random.multinomial(value, probs)
        probs_sum.append(probs)
        noise_data.append(multinomial_data)

    project_data = np.hstack((project_data, sink_abundances))
    noise_data = np.array(noise_data)
    noise_data = noise_data.T
    noise_sink_data = np.hstack((noise_data, sink_abundances))

    return project_data, noise_sink_data




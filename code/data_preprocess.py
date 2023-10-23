import pandas as pd
import numpy as np
from data_output import sink_data_output
from data_cluster import jsd_correlation_matrix
from data_cluster import data_cluster
from data_estimation import source_jsd_estimation
from data_estimation import noise_source_jsd_estimation


def data_to_excel(workbook, sheet_name, iters, matrix):
    worksheet = workbook.create_sheet(sheet_name, index=iters)
    for j in range(1, matrix.shape[0] + 1):
        for k in range(1, matrix.shape[1] + 1):
            worksheet.cell(row=j + 1, column=k, value=matrix[j - 1][k - 1])
    return worksheet


def standardization(matrix):
    matrix_sum = np.sum(matrix, axis=0)
    matrix = matrix / matrix_sum[np.newaxis, :]
    return matrix


def ground_truth_exchange(estimated_proportion, path):
    original_data = np.array(pd.read_csv(path, header=0, index_col=0))
    source_data = original_data[:, 0:19]
    sink_data = original_data[:, 20:23]

    source_data_sum = []
    for i in range(source_data.shape[1]):
        original_source = np.sum(source_data[:, i], axis=0)
        source_data_sum.append(original_source)

    sink_data_sum = []
    for i in range(sink_data.shape[1]):
        original_sink = np.sum(sink_data[:, i], axis=0)
        sink_data_sum.append(original_sink)

    result_sum = []
    for i in range(estimated_proportion.shape[0]):
        result = estimated_proportion[i, 0:11] * sink_data_sum[i] / source_data_sum
        unknown = 1 - np.sum(result)
        result = np.append(result, unknown)
        result_sum.append(result)
    result_sum = np.array(result_sum)
    return result_sum


def abundance_exchange(proportion_path, path, sources_num):
    project_data = np.array(pd.read_excel(path, header=0, sheet_name='Sheet1'))
    source_data = project_data[:, 0:sources_num]
    true_proportion = np.array(pd.read_excel(proportion_path, header=0, sheet_name='true_proportions'))

    source_data_sum = []
    for i in range(source_data.shape[1]):
        original_source = np.sum(source_data[:, i], axis=0)
        source_data_sum.append(original_source)

    result_sum = []
    for i in range(true_proportion.shape[0]):
        result = true_proportion[i, :] * source_data_sum
        result = result / np.sum(result)
        result_sum.append(result)
    result_sum = np.array(result_sum)
    return result_sum


def data_output(path1, sink_number, proportion_path, value):
    project_data = np.array(pd.read_excel(path1, header=0, sheet_name='Sheet1'))
    sources_num = project_data.shape[1]
    project_data, noise_sink_data = sink_data_output(project_data, sink_number, proportion_path, value)
    feast_data = np.hstack((noise_sink_data[:, 0:sources_num-1], noise_sink_data[:, sources_num:sources_num+sink_number]))
    source_data = project_data[:, 0:sources_num-1]
    source_jsd = source_jsd_estimation(source_data)
    noise_source_jsd = noise_source_jsd_estimation(noise_sink_data[:, 0:sources_num], project_data[:, 0:sources_num])
    return project_data, feast_data, noise_sink_data, source_jsd, sources_num, noise_source_jsd


def data_clustering(sink_number, noise_sink_data, sink_abundance, sources_num, jsd_value):
    noise_pre = noise_sink_data[:, 0:sources_num-1]
    noise_sink_pre = noise_sink_data[:, 0:sources_num-1]
    sink_abundance_pre = sink_abundance[:, 0:sources_num-1]
    noise_pre = standardization(noise_pre)
    corr_matrix = jsd_correlation_matrix(noise_pre)
    source_index, noise_sink_pre, sink_abundance_pre = data_cluster(noise_pre, corr_matrix, noise_sink_pre, sink_abundance_pre, jsd_value)

    unknown_pre = np.reshape(noise_sink_data[:, sources_num - 1], (noise_sink_data[:, sources_num - 1].shape[0], 1))
    sink_pre = noise_sink_data[:, sources_num:sources_num+sink_number]
    noise_sink_data = np.concatenate((noise_sink_pre, unknown_pre, sink_pre), axis=1)
    unknownpro_pre = np.reshape(sink_abundance[:, sources_num - 1], (sink_abundance[:, sources_num - 1].shape[0], 1))
    sink_abundance = np.concatenate((sink_abundance_pre, unknownpro_pre), axis=1)
    sources_num = noise_sink_data[:, 0:noise_sink_pre.shape[1]+1].shape[1]
    return noise_sink_data, sink_abundance, source_index, sources_num


def data_preprocess(iters, noise_sink_data, sources_num):
    sink_data_pre = np.reshape(noise_sink_data[:, sources_num+iters], (noise_sink_data[:, sources_num+iters].shape[0], 1))
    preprocess_data = np.hstack((noise_sink_data[:, 0:sources_num], sink_data_pre))
    preprocess_pre = np.hstack((noise_sink_data[:, 0:sources_num-1], sink_data_pre))
    remove_pre = preprocess_data[[not np.all(preprocess_pre[:, 0:preprocess_pre.shape[1]][i] == 0) for i in range(preprocess_pre.shape[0])], :]
    remove_data = standardization(remove_pre)

    noise_data = remove_data[:, 0:sources_num-1]
    sink_data = remove_data[:, sources_num]
    original_data = remove_data[:, 0:sources_num]

    unknown_pre = sink_data - np.sum(noise_data, axis=1)
    unknown_zero = np.zeros((unknown_pre.shape[0],))
    unknown_init = np.vstack((unknown_pre, unknown_zero))
    unknown_init = np.max(unknown_init, axis=0)
    sink_data = np.reshape(sink_data, (sink_data.shape[0], 1))
    unknown_init = np.reshape(unknown_init, (unknown_init.shape[0], 1))
    return remove_data, noise_data, sink_data, unknown_init, original_data






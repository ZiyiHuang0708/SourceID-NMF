import numpy as np
import pandas as pd
from scipy.spatial.distance import jensenshannon


def jsd_correlation_matrix(source_data):
    num = len(source_data[0])
    record = []
    jsd = np.zeros([num, num])
    for i in range(num - 1):
        for j in range(i + 1, num):
            value = jensenshannon(source_data[:, i], source_data[:, j])
            jsd[i, j] = value
            jsd[j, i] = value
            record.append(value)
    record = np.array(jsd)
    return record


def min_value(corr_matrix):
    min_val = corr_matrix[0][0]
    row_index, col_index = 0, 0
    for i in range(0, corr_matrix.shape[0]):
        for j in range(0, corr_matrix.shape[1]):
            if corr_matrix[j][i] <= min_val:
                min_val = corr_matrix[i][j]
                row_index = i
                col_index = j
    return min_val, row_index, col_index


def index_update(source_index, row_index, col_index):
    if isinstance(source_index[row_index], int) is True:
        source_index_update = [source_index[row_index], source_index[col_index]]
    else:
        source_index[row_index].append(source_index[col_index])
        source_index_update = source_index[row_index]
    source_index.pop(row_index)
    source_index.pop(col_index)
    source_index.append(source_index_update)
    return source_index


def update_matrix(source_data, row_index, col_index):
    cluster_data = np.sum(source_data.take([row_index, col_index], axis=1), axis=1)
    cluster_data = np.reshape(cluster_data, (cluster_data.shape[0], 1))
    delete_data = np.delete(source_data, [row_index, col_index], axis=1)
    source_data_update = np.concatenate((delete_data, cluster_data), axis=1)
    corr_matrix_update = jsd_correlation_matrix(source_data_update)
    np.fill_diagonal(corr_matrix_update, 1)
    return source_data_update, corr_matrix_update


def update_proportion(sink_abundance, row_index, col_index):
    cluster_proportion = np.sum(sink_abundance.take([row_index, col_index], axis=1), axis=1)
    cluster_proportion = np.reshape(cluster_proportion, (cluster_proportion.shape[0], 1))
    delete_proportion = np.delete(sink_abundance, [row_index, col_index], axis=1)
    sink_abundance_update = np.concatenate((delete_proportion, cluster_proportion), axis=1)
    return sink_abundance_update


def data_cluster(noise_data, corr_matrix, noise_sink_data, sink_abundance, jsd_value):
    source_index = list(range(noise_data.shape[1]))
    np.fill_diagonal(corr_matrix, 1)
    for i in range(1, noise_data.shape[1]):
        min_val, row_index, col_index = min_value(corr_matrix)
        if min_val > jsd_value:
            break
        source_index = index_update(source_index, row_index, col_index)
        noise_data, corr_matrix = update_matrix(noise_data, row_index, col_index)
        noise_sink_data, _ = update_matrix(noise_sink_data, row_index, col_index)
        sink_abundance = update_proportion(sink_abundance, row_index, col_index)
    return source_index, noise_sink_data, sink_abundance
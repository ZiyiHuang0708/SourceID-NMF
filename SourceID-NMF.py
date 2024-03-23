import argparse
import warnings
import numpy as np
import pandas as pd
from admm import lagrangian
from admm import admm_step
from data_clustering import jsd_correlation_matrix
from data_clustering import data_cluster
from data_estimation import jsd_estimation
from data_estimation import source_jsd_estimation
from data_initialization import standardization
from data_initialization import unknown_initialize
from data_initialization import parameters_initialize
warnings.filterwarnings('ignore')

parser = argparse.ArgumentParser(description="""Main script of SourceID-NMF.""")
parser.add_argument('--input', '-i', help='Input the original source data and sink data', default='./nmf_data.txt')
parser.add_argument('--name', '-n', help='Input the source data and sink data labels', default='./name.txt')
parser.add_argument('--output', '-o', help='Output the result of the mixing proportions', default='./estimated_proportions.txt')
parser.add_argument('--mode', '-m', help='Whether to optimize using clustering algorithms', default='normal')
parser.add_argument('--cutoff', '-f', help='Thresholding for clustering algorithms', type=float, default=0.25)
parser.add_argument('--thread', '-t', help='Number of Threads for multiprocessing', type=int, default=20)
parser.add_argument('--iter', '-e', help='Number of iterations per round for the model', type=int, default=2000)
parser.add_argument('--rho', '-r', help='rho', type=int, default=1)
parser.add_argument('--A', '-a', help='weighting matrix factor', type=int, default=1)
parser.add_argument('--threshold', '-c', help='threshold of convergence of Lagrangian functions', type=float, default=1e-06)
parser.add_argument('--perf', '-p', help='Whether you need to export the performance of the model', default=None)
inputs = parser.parse_args()

data_path   = inputs.input
name_path   = inputs.name
output_path = inputs.output
mode        = inputs.mode
cutoff      = inputs.cutoff
thread      = inputs.thread
iteration   = inputs.iter
rho         = inputs.rho
A           = inputs.A
threshold   = inputs.threshold
perf_mode   = inputs.perf

# output sources / sinks
data = pd.read_csv(data_path, sep="	", header=0, index_col=0)
name = pd.read_csv(name_path, sep="	", header=0)
name_id = list(name.loc[:, 'SampleID'])
source_sink_id = list(name.loc[:, 'SourceSink'])
sources = []
sinks = []
sources_label = []
sinks_label = []
for i in range(len(source_sink_id)):
    if source_sink_id[i] == 'Source':
        source = data[name_id[i]].values
        sources.append(source)
        sources_label.append(name_id[i])
    elif source_sink_id[i] == 'Sink':
        sink = data[name_id[i]].values
        sinks.append(sink)
        sinks_label.append(name_id[i])
    else:
        print('Error occurred: abnormal category in data preprocessing')
sources = np.array(sources).T  # organize and summarize all sources to prepare inputs for the model
sinks = np.array(sinks).T  # organize and summarize all sinks to prepare inputs for the model
sources_label.append('unknown')  # adding unknown sources to the source's label
sources_jsd = source_jsd_estimation(sources)
print("The Jensen Shannon divergence between the original source data is", sources_jsd)

if mode == 'cluster':
    corr_matrix = jsd_correlation_matrix(sources)
    source_index, sources = data_cluster(sources, corr_matrix, cutoff)
    print("The newly generated combination of sources after clustering is", source_index)
    print("The number of sources after clustering is", sources.shape[1])

# data preprocessing and training
estimated_proportions = []  # output the final estimated proportion
jsd_wy = []  # output the jsd between the w_plus and the y
diff_xwh = []  # output the difference between the x and np.dot(w, h)
nmf_perf = []  # output the final performance including jsd_wy and diff_xwh
for i in range(sinks.shape[1]):
    print("The model is running for the", i+1, "sink and all sources")
    # remove non-contributing sources
    sources_sums = np.sum(sources, axis=0)
    zero_sum_columns = np.where(sources_sums == 0)[0]
    sources = np.delete(sources, zero_sum_columns, axis=1)

    # remove the taxon of 0 between sources with sinks and standardize
    merging_data = np.hstack((sources, sinks[:, i].reshape((-1, 1))))
    remove_data = merging_data[[not np.all(merging_data[:, 0:merging_data.shape[1]][i] == 0) for i in range(merging_data.shape[0])], :]
    remove_data = standardization(remove_data)
    source = remove_data[:, :-1]
    sink = remove_data[:, -1].reshape((-1, 1))

    # initialize the parameters of the input model
    unknown_init = unknown_initialize(source, sink)
    x, y, w, a, i, w_plus, h_plus, alpha_w, alpha_h = parameters_initialize(source, sink, unknown_init, A)

    # NMF model running
    loss = []
    w_renew, h_renew, w_plus, h_plus, alpha_w, alpha_h = admm_step(w, x, y, a, i, w_plus, h_plus, alpha_w, alpha_h, rho, thread)
    l_init = lagrangian(x, y, w_renew, h_renew, a, w_plus, h_plus, alpha_w, alpha_h, rho)
    loss.append(l_init)
    for m in range(1, iteration):
        w_renew, h_renew, w_plus, h_plus, alpha_w, alpha_h = admm_step(w_renew, x, y, a, i, w_plus, h_plus, alpha_w, alpha_h, rho, thread)
        l_update = lagrangian(x, y, w_renew, h_renew, a, w_plus, h_plus, alpha_w, alpha_h, rho)
        loss.append(l_update)
        if m == 1 and abs(l_update - l_init) / l_init < threshold:
            break
        if m > 1 and abs(l_update - loss[m - 1]) / loss[m - 1] < threshold:
            break
        else:
            continue
    jsd_wy_iter = jsd_estimation(w_plus[:, :-1], y)
    jsd_wy.append(jsd_wy_iter)
    diff_xwh_iter = np.sum(abs(x - np.dot(w_renew, h_renew)))
    diff_xwh.append(diff_xwh_iter)
    estimated_proportions.append(h_plus.reshape((h_plus.shape[0])))

# output the final estimated proportions
estimated_proportions = np.array(estimated_proportions)
estimated_proportions = pd.DataFrame(estimated_proportions, index=sinks_label, columns=sources_label)
estimated_proportions.to_csv(output_path, sep=' ')

# output the final performance
if perf_mode == 'perf_needed':
    nmf_perf.append(jsd_wy)
    nmf_perf.append(diff_xwh)
    nmf_perf_label = ['jsd_wy', 'diff_xwh']
    nmf_perf = np.array(nmf_perf).T
    nmf_perf = pd.DataFrame(nmf_perf, index=sinks_label, columns=nmf_perf_label)
    nmf_perf.to_csv('./nmf_perf.txt', sep=' ')

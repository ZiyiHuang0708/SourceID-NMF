# SourceID-NMF
SourceID-NMF: Towards more accurate microbial source tracking via non-negative matrix factorization.

A major challenge of analyzing the compositional structure of microbiome data is identifying its potential origins. Here, we introduce a novel tool called SourceID-NMF for precise microbial source tracking. SourceID-NMF utilizes a non-negative matrix factorization (NMF) algorithm to trace the microbial sources contributing to a target sample, without relying on specific probability distributions.

## Support
For support using SourceID-NMF, please email: zhuang82-c@my.cityu.edu.hk

## Required Dependencies
Detailed package information can be found in SourceID-NMF.yaml. The main environment configuration we need includes:
* Conda
* Python >=3.8.13
* numpy >=1.24.3
* pandas >=2.0.3
* tqdm >=4.66.1

We suggest to install SourceID-NMF's environment SourceID-NMF.yaml by using anaconda after cloning the repository. This will install all the required packages in cpu mode.

The command is: 
```
conda env create -f SourceID-NMF.yaml -n nmf
```

Alternatively, you can use conda to install all packages, as shown on the command line below:
```
# create a new environment nmf with python 3.8 using conda
conda create -n nmf python=3.8.13

# install some basic conda packages to support the modeling operations
conda install -c anaconda numpy
conda install -c anaconda pandas
conda install -c anaconda scipy

# multiprocessing tqdm progress bar display
conda install -c conda-forge tqdm
```

After installing, you only need to activate the "nmf" environment before using SourceID-NMF for the next time.
```
conda activate nmf
```


## Usage

#### Command
```
python SourceID-NMF.py -i ./nmf_data.txt -n ./name.txt -o ./estimated_proportions.txt -t 20 -e 20 -r 1 -a 1 -c 1e-06
```

#### Parameters

As input, SourceID-NMF provides a number of parameters for different purposes. Some of the parameters exist as inputs and outputs of data, including:

```
Options
-i | --input:        A path to an input.txt file: Input of count table.
-n | --name:         A path to an name.txt file: Data labels for input data. 
-o | --output:       A path to an output.txt file: Output of estimated proportions.
-m | --perf:         A path to an output.txt file: Output of model performance including jsd_wy and diff_xwh.
```
The input to SourceID-NMF is composed of two txt files:

* Input

The input count table containing sources and sinks (M by N). where M is the number samples and N is the number of taxa. Row names are the sample ids ('SampleID'). Column names are the taxa ids. Every consecutive column contains read counts for each sample.

The specific input table case is shown below:

| | D1 | D2 | D3 | ... | D27 | D28 |
| ------------- | ------------- |------------- |------------- |------------- |------------- |------------- |
| taxon_1  |  0 | 5 | 0 | ... | 20 | 5 |
| taxon_2  |  20 | 5 | 0 | ... | 0 | 11 |
| taxon_3  |  0 | 13 | 210 | ... | 0 | 20 |
| taxon_4  |  80 | 6 | 0 | ... | 0 | 0 |
| taxon_5  |  4 | 38 | 0 | ... | 14 | 0 |
| ... | ... | ... | ... | ... | ... | ... |
| taxon_n  |  24 | 25 | 0 | ... | 0 | 14 |

* Name

The name table containing four columns, 'SampleID', 'Env' and 'SourceSink'. The 'SampleID' column describes the labels for each source data or sink data. The 'Env' column describes the environment to which each source or sink belongs, e.g. the first row Env = 'Electronics' means that the source was collected from Electronics. This 'SourceSink' column describes the source or sink to which the data belongs. 

The specific name table case is shown below:

| SampleID | Env |SourceSink |
| ------------- | ------------- |------------- |
| D1 | Electronic | Source |
| D2 | Hand | Source |
| D3 | Incubator | Source|
| D4 | Surface | Source|
| ... | ... | ... |
| D27 | fecal8 | Sink |
| D28 | fecal9 | Sink |

The output to SourceID-NMF is composed of one txt files:

* Output

The count table contains all the estimated proportions (K by S). where K is the number sinks and S is the number of sources (including an unknown source). The specific value in this table represents the contribution of each source to each sink. The sum of the proportions in each row is 1.

The specific output table case is shown below:

| | D1 | D2 | D3 | ... | D19 | Unknown |
| ------------- | ------------- |------------- |------------- |------------- |------------- |------------- |
| D20 | 0.004057049 | 0.000486668	| 0.014082136 | ... | 0.068413642	| 0.048694265 |
| D21 | 0.055386974 | 0.013956499	| 0.011004204 | ... | 0.001979577	| 0.155531113 |
| D22 | 0.127182548 | 0.000429836	| 0.017506031 | ... | 0.016176403	| 0.280087424 |
| D23 | 0.114428201 | 0.068275637	| 0.042867499 | ... | 0.069508314	| 0.378383698 |
| ... | ... | ... | ... | ... | ... | ... |
| D29 | 0.038159300 | 0.001559895	| 0.003974004 | ... | 0.005082442	| 0.875595719 |

* Perf

This is an optional parameter for the user. When this parameter = 'perf_needed' it means that the user needs to output the final result of the model iteration. This output includes the similarity between W_plus and Y as measured by Jensen-Shannon divergence, and the difference between X and WH. With these results we can evaluate the performance of the model iteration. The specific commands that can be added are as follows:

```
python SourceID-NMF.py -p pref_needed
```

Meanwhile, SourceID-NMF also provides some parameters that may affect the final result, including:
```
Options
-t | --thread:       Max workers for multiprocessing operation.
-e | --iter:         Maximum number of iterations for the NMF model.
-r | --rho:          The penalty parameter.
-a | --A:            The weighting matrix coefficients. 
-c | --threshold:    The convergence threshold.
```

## Demo
Here, we provide some datasets for SourceID-NMF to test. The /data folder contains sections for simulated data and real data. The section on simulated data contains three sets of data available for testing. Each data set contains two txt input files 'nmf_data.txt' and 'name.txt'. We can run it on simulated data or real data by running the following command:

```
Command:
# In simulated data
cd ./data/simulated_data
cd ./0.60jsd / cd ./0.70jsd / cd ./0.80jsd
python SourceID-NMF.py -i ./nmf_data.txt -n ./name.txt -o ./estimated_proportions.txt -t 20 -e 20 -r 1 -a 1 -c 1e-06

# In true data
cd ./data/true_data
python SourceID-NMF.py -i ./nmf_data.txt -n ./name.txt -o ./estimated_proportions.txt -t 20 -e 20 -r 1 -a 1 -c 1e-06
```

After running the code, you can find the file 'estimated_proportions.txt' in the folder corresponding to that dataset, which contains the results of the model run, i.e., the contributions of the sources to the sinks. 

                     




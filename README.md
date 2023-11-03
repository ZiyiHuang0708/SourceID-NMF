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

The command is: `conda env create -f SourceID-NMF.yaml -n nmf`

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

## Usage

#### Command
```
python main.py -i ./nmf_data.txt -n ./name.txt -o ./estimated_proportions.txt -t 20 -e 20 -r 1 -a 1 -c 1e-06
```

#### Parameters

As input, SourceID-NMF provides a number of parameters for different purposes. Some of the parameters exist as inputs and outputs of data, including:

```
Options
-i | --input:        A path to an input.txt file.
-n | --name:         Data labels and tags for input data. 
-o | --output:       A path to an output.txt file.
```
Input: A count table containing sources and sinks (M by N). where M is the number samples and N is the number of taxa. Row names are the sample ids ('SampleID'). Column names are the taxa ids. Every consecutive column contains read counts for each sample.











                     




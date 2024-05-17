#!/usr/bin/env python


import numpy as np
import random
from plinkio import plinkfile
import matplotlib.pyplot as plt
import pandas as pd
import sys
import os
import glob
from statistics import mean, median
from operator import index
from math import sqrt
import statistics
import time
from optparse import OptionParser

VERSION='v2.00'

Usage = """

Options when running READv2:

\t -i, --input_file <val>\tInput file prefix (required). The current READ version only supports Plink bed/bim/fam files.
\t -n, --norm_method <val>\tNormalization method (either 'mean', 'median', 'max' or 'value').
\t\t\tmedian (default) - assuming that most pairs of compared individuals are unrelated, READ uses the median across all pairs for normalization.
\t\t\tmean - READ uses the mean across all pairs for normalization, this would be more sensitive to outliers in the data (e.g. recent migrants or identical twins)
\t\t\tmax - READ uses the maximum across all pairs for normalization. This should be used to test trios where both parents are supposed to be unrelated but the offspring is a first degree relative to both others.
\t\t\tvalue - READ uses a user-defined value for normalization. This can be used if the value from another population should be used for normalization. That would be useful if only two individuals are available for the test population. Normalization value needs to be provided through --norm_value
\t --window_size <val>\tChange window size for block jackknife or for window-based P0 estimates (as in READv1), default: 5000000
\t --window_est\tWindow based estimate of P0 (as opposed to the genome-wide estimate, default in READv2)
\t -h, --help\tPrint this message
\t -v, --version\tPrint version



"""

# turn of interactive plotting, does not show the plot unless plt.show() is used.
plt.ioff()


def Pair_Calc(window, out_matrix):
    missing = np.count_nonzero(window == -1)
    full_call = window.size
    full_call_var_corrected = full_call - missing
    if (full_call_var_corrected != 0):
        valor = np.count_nonzero(window == 1)
        IBS0 = np.count_nonzero(window == 0)
        P0 = float(IBS0)/full_call_var_corrected
        P2 = float(valor)/full_call_var_corrected
        Line_info = [pair, valor, IBS0, P2, P0, full_call_var_corrected]
        out_matrix.append(Line_info)


def Perc_Calc(x):  # a function to calculate the proportion of windows that does not fall between deg_thresholds[2] and deg_thresholds[0] for each pair
    length_windows = len(x)
    cnt = 0
    for el in x:
        if (el > deg_thresholds[1] or el < deg_thresholds[3]):
             cnt += 1
    return (cnt/length_windows)


def standard_error(x):
    return (statistics.stdev(x) / sqrt(len(x)))

def blockJN_SE(x):
    x=list(x)
    blockJN_means = []
    num_blocks=len(x)

    for i in range(num_blocks-1):
        block_mean = np.mean(x[:i]+x[i+1:])
        blockJN_means.append(block_mean)

    grand_mean=np.mean(blockJN_means)

    variance = np.sum((blockJN_means - grand_mean)**2)* ((num_blocks-1) / num_blocks)

    return (sqrt(variance))


Arguments = sys.argv[1:]
norm_method = "median"
window_size = 5000000 #default value chosen to correspond to commonly used blockJN block sizes, Readv1 default was 1000000
window_size_1stdeg = 20000000
deg_thresholds = [0.953125,0.90625,0.8125,0.625]
effectiveSNP_cutoff_first = 10000
effectiveSNP_cutoff_third = 5000
window_based = False
Check_1stdeg_type = True
list_all_individuals = []
previous_window = 0
Missing_info = 0
IBS2 = 0
IBS0 = 0
Total_snps_window = 0
full_call_var = 0
full_call_var_corrected = 0
ind_dict = {}  # a dict. to hold ind. snp info
pair_dict = {}  # a dict. to hold pair snp info
pair_matrix = []
pair_matrix_blockJN = []
pair_matrix_1stdeg = []
List_individuals = []


previous_window = 0
snp_count = 0

parser = OptionParser()
parser.add_option("-i","--input",action="store",type='string',dest='infile',help="Prefix of input Plink bed/bim/fam files")
parser.add_option("-n","--norm_method",action="store",type='string',dest='norm_method',help="Normalization method (either 'mean', 'median', 'max' or 'value')", default='median')
parser.add_option("--norm_value", action="store",type="float", dest="norm_value",help="User-specified normalization value",default=0.0)
parser.add_option("--window_size", action="store",type="int", dest="window_size",help="Window size (default 5000000)",default=5000000)
parser.add_option("--window_est", action="store_true", dest="window_based",help="Window based estimate of P0 (as opposed to the genome-wide estimate, default in READv2)",default=False)
parser.add_option("--2pow", action="store_true", dest="twopow_thresh",help="Use alternative classification thresholds",default=False)
parser.add_option("-v","--version", action="store_true", dest="print_vers",help="Print version",default=False)

(options, args) = parser.parse_args()

if options.print_vers:
    print(VERSION)
    sys.exit(0)


print("    ===Thank you for using Relationship Estimation from Ancient DNA version 2 (READv2)===    ", end='\n\n')

norm_method=options.norm_method
norm_value=options.norm_value
infile=options.infile
window_based=options.window_based
window_size=options.window_size

if norm_method not in ["median", "mean", "max", "value"]:
    print("No valid standardization method specified! Using the default method (i.e. median).")
    norm_method='median'

if norm_value:
    norm_method='value'
    print('User-specified normalization value provided.')
    
if norm_method=='value' and norm_value<=0:
    print('User-specified normalization value chosen but no valid value provided!!!')
    sys.exit(1)

if options.twopow_thresh: #use updated degree thresholds
    deg_thresholds = [1-1/(2**4.5),1-1/(2**3.5),1-1/(2**2.5),1-1/(2**1.5)]

plink_file = plinkfile.open(infile)
locus_list = plink_file.get_loci()
sample_list = plink_file.get_samples()

#Check for pseuodhaploidy, if any heterozygous site is detected, the program quits with an error message.
for row in plink_file:
    if(1 in row):
        print('Heterozygous site detected. Please make sure your input files are pseudo-haploid.')
        sys.exit(1)
    
for sample in sample_list:
    list_all_individuals.append(sample.iid)

if len(sample_list)<2:
    print("At least two individuals are required to run READ")
    sys.exit(1)

if len(sample_list)<4 and norm_method=="median":
	print('Warning: Normalization using the median can be problematic for small sample sizes.')

# locus.bp_position -> base pair position, locus.allele* gives the alleles in that position.

print('Window size: %s' %window_size)
if window_based:
     print('Calculating P0 in windows across the genome and using the mean P0 for classification, the reported standard error is the SE of that mean.')
else:
     print('Calculating one genome-wide P0 value and using a block-jackknife to estimate the standard error.')

start_time = time.time()
tmp_list = []
tmp_count = 0

plink_file.transpose(infile + "_test")
plink_file.close()
plink_file = plinkfile.open(infile + "_test")
i = 0
# prep. np arrays for each individual
locus_list = plink_file.get_loci()

locus_pos = np.array(list(map(lambda d: d.bp_position, locus_list)))
locus_pos = np.floor_divide(locus_pos, window_size)
windows = np.where(locus_pos[:-1] != locus_pos[1:])[0]

locus_pos_1stdeg = np.array(list(map(lambda d: d.bp_position, locus_list)))
locus_pos_1stdeg = np.floor_divide(locus_pos_1stdeg, window_size_1stdeg)
windows_1stdeg = np.where(locus_pos_1stdeg[:-1] != locus_pos_1stdeg[1:])[0]

for row in plink_file:
    ind_dict[list_all_individuals[i]] = np.array(row)
    i += 1

keys = ind_dict.keys()
num_pairs=len(keys)*(len(keys)-1)/2
pair_count=0
for key in keys:
    for key2 in keys:
        if (not (key == key2) and not ((key2+','+key) in pair_dict.keys())):
            nan_indices = np.union1d(np.where(
                ind_dict[key] == 3), np.where(
                ind_dict[key2] == 3))
            pairwise = 1 * (ind_dict[key] == ind_dict[key2])
            pairwise[nan_indices] = -1

            pair=key+","+key2

            pair_dict[pair]=0

            if len(nan_indices)==len(pairwise):
                print("!!! Warning: no overlapping SNPs for the pair %s. They won't be part of the output files." %pair)

            arr = pairwise
            if not window_based:
                Pair_Calc(arr, pair_matrix)
                pre_ind = 0
                for i in windows:
                    window = arr[pre_ind:i+1]
                    Pair_Calc(window, pair_matrix_blockJN)
                    pre_ind = i
                window = arr[pre_ind:]  # for the last window
                Pair_Calc(window, pair_matrix_blockJN)
                pre_ind = 0

            else:
                pre_ind = 0
                for i in windows:
                    window = arr[pre_ind:i+1]
                    Pair_Calc(window, pair_matrix)
                    pre_ind = i
                window = arr[pre_ind:]  # for the last window
                Pair_Calc(window, pair_matrix)
                pre_ind = 0

            if Check_1stdeg_type:
                pre_ind = 0
                for i in windows_1stdeg:
                    window = arr[pre_ind:i+1]
                    Pair_Calc(window, pair_matrix_1stdeg)
                    pre_ind = i
                window = arr[pre_ind:]  # for the last window
                Pair_Calc(window, pair_matrix_1stdeg)
                pre_ind = 0

            pair_count+=1
            if (pair_count % 1000 == 0):
                print("Pair %s out of %s processed." %(pair_count,int(num_pairs)))
            

#pair_matrix = np.array(pair_matrix)


#start_time = time.time()
df_pair = pd.DataFrame(pair_matrix, columns=[
                       'PairIndividuals', 'IBS2', 'IBS0', 'P1', 'P0', 'NSNPs'])

df_pair = df_pair.astype(dtype={
                         'PairIndividuals': str, 'IBS2': float, 'IBS0': float, 'P1': float, 'P0': float, 'NSNPs': int})


df_pair_1stdeg = pd.DataFrame(pair_matrix_1stdeg, columns=[
                       'PairIndividuals', 'IBS2', 'IBS0', 'P1', 'P0', 'NSNPs'])

df_pair_1stdeg = df_pair_1stdeg.astype(dtype={
                         'PairIndividuals': str, 'IBS2': float, 'IBS0': float, 'P1': float, 'P0': float, 'NSNPs': int})


del locus_list, plink_file, pair_matrix, ind_dict, pair_dict

df_pair['TwoAlleleDiffMeanPerc'] = df_pair.groupby(
    ['PairIndividuals'])['P0'].transform('mean')

zero_pairs=df_pair[df_pair.TwoAlleleDiffMeanPerc==0]['PairIndividuals']
for zp in zero_pairs:
    print('Warning: The pair %s seems completely identical. Please check if this is a duplicate!' %zp)

if (norm_method == "median"):
    norm_value = median(df_pair['TwoAlleleDiffMeanPerc'])
elif (norm_method == "mean"):
    norm_value = mean(df_pair['TwoAlleleDiffMeanPerc'])
elif (norm_method == "max"):
    norm_value = max(df_pair['TwoAlleleDiffMeanPerc'])

if norm_value==0.0:
    print("!!! Normalization value is zero. Are many individuals duplicates? Terminating normalization, please fix and rerun READ !!!")
    sys.exit(1)

start_time = time.time()
df_pair['Norm2AlleleDiff'] = df_pair['P0'] / norm_value

df_pair_1stdeg['Norm2AlleleDiff'] = df_pair_1stdeg['P0'] / norm_value

Perc_cal_pair = df_pair_1stdeg.groupby(['PairIndividuals'])[
    'Norm2AlleleDiff'].apply(Perc_Calc).to_frame() #now always on 20mbp windows unless changed

nsnps_pair = df_pair.groupby(
    ['PairIndividuals'])['NSNPs'].sum().to_frame()

means2Allele_Diff_Normalized = df_pair.groupby(
    ['PairIndividuals'])['Norm2AlleleDiff'].mean().to_frame()
means2Allele_Diff_Normalized['Perc_P0'] = Perc_cal_pair['Norm2AlleleDiff']
means2Allele_Diff_Normalized['OverlapNSNPs'] = nsnps_pair['NSNPs']
means2Allele_Diff_Normalized['NSNPsXNorm'] = means2Allele_Diff_Normalized['OverlapNSNPs']*norm_value
df_P0 = df_pair.groupby(['PairIndividuals'])[
    'P0'].mean().to_frame()  # unnormalized P0

if not window_based: #calculate SEs using a block jackknife
    df_pair_blocks = pd.DataFrame(pair_matrix_blockJN, columns=[
                       'PairIndividuals', 'IBS2', 'IBS0', 'P1', 'P0', 'NSNPs'])
    df_pair_blocks = df_pair_blocks.astype(dtype={
                         'PairIndividuals': str, 'IBS2': float, 'IBS0': float, 'P1': float, 'P0': float, 'NSNPs': int})
    df_pair_blocks['Norm2AlleleDiff'] = df_pair_blocks['P0'] / norm_value
    StError_2Allele_Norm = df_pair_blocks.groupby(['PairIndividuals'])['Norm2AlleleDiff'].apply(blockJN_SE).to_frame()

else: #SE across windows
    StError_2Allele_Norm = df_pair.groupby(['PairIndividuals'])['Norm2AlleleDiff'].apply(standard_error).to_frame()  # normalized P0 standard error

means2Allele_Diff_Normalized['StError_2Allele_Norm'] = StError_2Allele_Norm['Norm2AlleleDiff']
Nonnormalized_P0_serr = StError_2Allele_Norm * norm_value #df_pair.groupby(['PairIndividuals'])['P0'].apply(standard_error).to_frame()  # nonNormalized P0 standard error
means2Allele_Diff_Normalized['Nonnormalized_P0'] = df_P0['P0']
means2Allele_Diff_Normalized['Nonnormalized_P0_serr'] = Nonnormalized_P0_serr
df_P0['Error'] = np.real(Nonnormalized_P0_serr)
df_P0.sort_values(by=['P0'],inplace = True)
df_P0 = df_P0.reset_index()


means2Allele_Diff_Normalized['theta'] = 1 - means2Allele_Diff_Normalized['Norm2AlleleDiff']


#means2Allele_Diff_Normalized = means2Allele_Diff_Normalized.astype(
#    dtype={'Norm2AlleleDiff': np.float16, 'StError_2Allele_Norm': np.float16, 'Nonnormalized_P0': np.float16, 'Nonnormalized_P0_serr': np.float16})
means2Allele_Diff_Normalized[['Norm2AlleleDiff', 'StError_2Allele_Norm', 'Nonnormalized_P0',
                              'Nonnormalized_P0_serr','OverlapNSNPs']].to_csv("meansP0_AncientDNA_normalized_READv2", index=True, sep='\t')

if num_pairs<=1000:
    start_time = time.time()
    # boxplot=df_P0.reset_index().boxplot(by = 'PairIndividuals',rot=90, figsize=((8+0.1*(len(df_P0.index)),6)))
    boxplot = df_P0.plot(x='PairIndividuals', y='P0', kind="scatter", yerr='Error',
                     s=6, rot=90, figsize=((8+0.1*(len(df_P0.index)), 8)), fontsize=8)

    # hor. lines for 2nd degree
    boxplot.axhline(norm_value*deg_thresholds[1], color="purple", alpha=0.2, lw=10)
    boxplot.axhline(norm_value*deg_thresholds[1], color="black", alpha=0.5, linestyle="--")
    # why use normalized values while you have nonnormalized values in the graph.
    boxplot.text(0.55, (norm_value*deg_thresholds[1])*0.95,
             "2nd degree", size=9, color="purple")
    # hor. lines for 1st degree
    boxplot.axhline(norm_value*deg_thresholds[2], color="purple", alpha=0.2, lw=10)
    boxplot.axhline(norm_value*deg_thresholds[2], color="black", alpha=0.5, linestyle="--")
    boxplot.text(0.55, (norm_value*deg_thresholds[2])*0.95,
             "1st degree", size=9, color="purple")
    # hor.lines for Identical/twins
    boxplot.axhline(norm_value*deg_thresholds[3], color="purple", alpha=0.2, lw=10)
    boxplot.axhline(norm_value*deg_thresholds[3], color="black", alpha=0.5, linestyle="--")
    boxplot.text(0.55, (norm_value*deg_thresholds[3])*0.95,
             "Identical/Twin", size=9, color="purple")
    fig = boxplot.get_figure()
    fig.savefig('READ.pdf', format="pdf")
    # print("Plot created!")
else:
    print("No plotting for more than 1000 pairs.")


start_time = time.time()
filters = [
    (means2Allele_Diff_Normalized.Norm2AlleleDiff <= deg_thresholds[0]) & (means2Allele_Diff_Normalized.Norm2AlleleDiff >= deg_thresholds[1]) & (means2Allele_Diff_Normalized.NSNPsXNorm >= effectiveSNP_cutoff_third),
    (means2Allele_Diff_Normalized.Norm2AlleleDiff <= deg_thresholds[1]) & (means2Allele_Diff_Normalized.Norm2AlleleDiff >= deg_thresholds[2]),
    (means2Allele_Diff_Normalized.Norm2AlleleDiff <= deg_thresholds[2]) & (means2Allele_Diff_Normalized.Norm2AlleleDiff >= deg_thresholds[3]),
    means2Allele_Diff_Normalized.Norm2AlleleDiff <= deg_thresholds[3]
]
values_Rel = ["Third Degree", "Second Degree", "First Degree",'IdenticalTwins/SameIndividual']
values_Zup = [
    np.absolute((deg_thresholds[0]-means2Allele_Diff_Normalized.Norm2AlleleDiff) /
                StError_2Allele_Norm.Norm2AlleleDiff),
    np.absolute((deg_thresholds[1]-means2Allele_Diff_Normalized.Norm2AlleleDiff) /
                StError_2Allele_Norm.Norm2AlleleDiff),
    np.absolute((deg_thresholds[2]-means2Allele_Diff_Normalized.Norm2AlleleDiff) /
                StError_2Allele_Norm.Norm2AlleleDiff),
    np.absolute((deg_thresholds[3]-means2Allele_Diff_Normalized.Norm2AlleleDiff) /
                StError_2Allele_Norm.Norm2AlleleDiff)
]
values_Zdown = [
    np.absolute((deg_thresholds[1]-means2Allele_Diff_Normalized.Norm2AlleleDiff) /
                StError_2Allele_Norm.Norm2AlleleDiff),
    np.absolute((deg_thresholds[2]-means2Allele_Diff_Normalized.Norm2AlleleDiff) /
                StError_2Allele_Norm.Norm2AlleleDiff),
    np.absolute((deg_thresholds[3]-means2Allele_Diff_Normalized.Norm2AlleleDiff) /
                StError_2Allele_Norm.Norm2AlleleDiff),
    'NA'
]

#var_df = df_pair.groupby(['PairIndividuals'])[
#    'Norm2AlleleDiff'].apply(statistics.variance).to_frame()


df_to_print = pd.DataFrame(
    data={'PairIndividuals': means2Allele_Diff_Normalized.index.get_level_values(0).values, 'P0_mean': means2Allele_Diff_Normalized['Norm2AlleleDiff'],
           'Perc_Win_1stdeg_P0': means2Allele_Diff_Normalized['Perc_P0'], 'OverlapNSNPs': means2Allele_Diff_Normalized['OverlapNSNPs'], 'Nonnormalized_P0':means2Allele_Diff_Normalized['Nonnormalized_P0'], 'Nonnormalized_P0_serr':means2Allele_Diff_Normalized['Nonnormalized_P0_serr'], 'NSNPsXNorm': means2Allele_Diff_Normalized['NSNPsXNorm'], 'KinshipCoefficient': means2Allele_Diff_Normalized['theta'] })


df_to_print['Rel'] = np.select(
    filters, values_Rel, default="Unrelated")
#df_to_print['Zup'] = np.select(filters, values_Zup, default=np.abs(
#    (deg_thresholds[3]-means2Allele_Diff_Normalized.Norm2AlleleDiff)/StError_2Allele_Norm.Norm2AlleleDiff))
#df_to_print['Zdown'] = np.select(filters, values_Zdown, default='NA')
df_to_print['Zup'] = np.select(filters, values_Zup, default='NA')
df_to_print['Zdown'] = np.select(filters, values_Zdown, default=np.abs((deg_thresholds[1]-means2Allele_Diff_Normalized.Norm2AlleleDiff)/StError_2Allele_Norm.Norm2AlleleDiff))


filters_first_Deg = [
    (means2Allele_Diff_Normalized.Perc_P0 >= 0.6) & (
        df_to_print.Rel == "First Degree") & (means2Allele_Diff_Normalized.NSNPsXNorm>=effectiveSNP_cutoff_first),
    (means2Allele_Diff_Normalized.Perc_P0 >= 0.35) & (
        df_to_print.Rel == "First Degree") & (means2Allele_Diff_Normalized.NSNPsXNorm>=effectiveSNP_cutoff_first),
    (means2Allele_Diff_Normalized.Perc_P0 <= 0.3) & (
        df_to_print.Rel == "First Degree") & (means2Allele_Diff_Normalized.NSNPsXNorm>=effectiveSNP_cutoff_first)
]
values_first_deg = [
    "NotApplicable", "Siblings", "Parent-offspring"
]

df_to_print['1st_Type'] = np.select(
    filters_first_Deg, values_first_deg, default="N/A")


df_to_print[['PairIndividuals', 'Rel', 'Zup', 'Zdown','P0_mean','Nonnormalized_P0', 'Nonnormalized_P0_serr', '1st_Type', 'Perc_Win_1stdeg_P0','OverlapNSNPs','NSNPsXNorm','KinshipCoefficient']].to_csv(
    'Read_Results.tsv', index=False, sep='\t')

fname = "./"+infile + "_test*"
for filename in glob.glob(fname):
    os.remove(filename)
print("READ analysis finished. Please check Read_Results.tsv for results!")

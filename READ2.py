#!/usr/bin/python
# PlinkIO


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


Usage = """

Current READ version only supports .tped files.

READ2.py <InputFile> <normalization> <normalization_value>

A normalization is only required when a user-defined values is used instead of the median, mean or maximum of the test population. All normalization settings are described below:

    median (default) - assuming that most pairs of compared individuals are unrelated, READ uses the median across all pairs for normalization.
    mean - READ uses the mean across all pairs for normalization, this would be more sensitive to outliers in the data (e.g. recent migrants or identical twins)
    max - READ uses the maximum across all pairs for normalization. This should be used to test trios where both parents are supposed to be unrelated but the offspring is a first degree relative to both others.
    value <val> - READ uses a user-defined value for normalization. This can be used if the value from another population should be used for normalization. That would be useful if only two individuals are available for the test population.

Optionally, one can add --window_size <value> at the end in order to modify the window size used (default: 5000000).

"""

print("    ===Thank you for using Relationship Estimation from Ancient DNA (READ)===    ", end='\n\n')


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


def Perc_Calc(x):  # a function to calculate the proportion of windows that does not fall between 0.625 and 0.90625 for each pair
    length_windows = len(x)
    cnt = 0
    for el in x:
        if (el > 0.90625 or el < 0.625):
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
norm_value = ''
window_size = 5000000 #default value chosen to correspond to commonly used blockJN block sizes, Readv1 default was 1000000
window_size_1stdeg = 20000000
Single_Sites = False
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

# Input check!
if (len(Arguments) == 0):
    sys.exit(Usage)

if (len(Arguments) >= 2 and Arguments[1] in ["median", "mean", "max", "value"]):
    norm_method = Arguments[1]
elif (len(Arguments) >= 2 and Arguments[1] not in ["median", "mean", "max", "value"]):
    print("No valid standardization method specified! Using the default method.")
if (len(Arguments) >= 3 and Arguments[1] == "value"):
    norm_value = float(Arguments[2])
if ("--window_size" in Arguments):
    ws_index = Arguments.index("--window_size")
    if (len(Arguments) > (ws_index+1)):
        window_size = int(Arguments[ws_index+1])
    else:
        sys.exit("No window size specified!")
if ("--blocks" in Arguments):
    Single_Sites = True

plink_file = plinkfile.open(Arguments[0])

sample_list = plink_file.get_samples()
for sample in sample_list:
    list_all_individuals.append(sample.iid)

# locus.bp_position -> base pair position, locus.allele* gives the alleles in that position.

start_time = time.time()
tmp_list = []
tmp_count = 0

plink_file.transpose(Arguments[0] + "_test")
plink_file.close()
plink_file = plinkfile.open(Arguments[0] + "_test")
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
        if (not (key == key2) and not (key2+key in pair_dict.keys())):
            nan_indices = np.union1d(np.where(
                ind_dict[key] == 3), np.where(
                ind_dict[key2] == 3))
            pairwise = 1 * (ind_dict[key] == ind_dict[key2])
            pairwise[nan_indices] = -1

            pair=key+key2

            pair_dict[pair]=0

            arr = pairwise
            if Single_Sites:
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

df_pair['2AlleleDiffMeanPerc'] = df_pair.groupby(
    ['PairIndividuals'])['P0'].transform('mean')

if (norm_method == "median"):
    norm_value = median(df_pair['2AlleleDiffMeanPerc'])
elif (norm_method == "mean"):
    norm_value = mean(df_pair['2AlleleDiffMeanPerc'])
elif (norm_method == "max"):
    norm_value = max(df_pair['2AlleleDiffMeanPerc'])

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

if Single_Sites: #calculate SEs using a block jackknife
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
    boxplot.axhline(norm_value*0.90625, color="purple", alpha=0.2, lw=10)
    boxplot.axhline(norm_value*0.90625, color="black", alpha=0.5, linestyle="--")
    # why use normalized values while you have nonnormalized values in the graph.
    boxplot.text(0.55, (norm_value*0.90625)-0.005,
             "2nd degree", size=9, color="purple")
    # hor. lines for 1st degree
    boxplot.axhline(norm_value*0.8125, color="purple", alpha=0.2, lw=10)
    boxplot.axhline(norm_value*0.8125, color="black", alpha=0.5, linestyle="--")
    boxplot.text(0.55, (norm_value*0.8125)-0.005,
             "1st degree", size=9, color="purple")
    # hor.lines for Identical/twins
    boxplot.axhline(norm_value*0.625, color="purple", alpha=0.2, lw=10)
    boxplot.axhline(norm_value*0.625, color="black", alpha=0.5, linestyle="--")
    boxplot.text(0.55, (norm_value*0.625)+0.005,
             "Identical/Twin", size=9, color="purple")
    fig = boxplot.get_figure()
    fig.savefig('READ.pdf', format="pdf")
    # print("Plot created!")
else:
    print("No plotting for more than 1000 pairs.")


start_time = time.time()
filters = [
    means2Allele_Diff_Normalized.Norm2AlleleDiff >= 0.90625,
    means2Allele_Diff_Normalized.Norm2AlleleDiff >= 0.8125,
    means2Allele_Diff_Normalized.Norm2AlleleDiff >= 0.625
]
values_Rel = ["Unrelated", "Second Degree", "First Degree"]
values_Zup = [
    "NA",
    np.absolute((0.90625-means2Allele_Diff_Normalized.Norm2AlleleDiff) /
                StError_2Allele_Norm.Norm2AlleleDiff),
    np.absolute((0.8125-means2Allele_Diff_Normalized.Norm2AlleleDiff) /
                StError_2Allele_Norm.Norm2AlleleDiff),
]
values_Zdown = [
    np.absolute((0.90625-means2Allele_Diff_Normalized.Norm2AlleleDiff) /
                StError_2Allele_Norm.Norm2AlleleDiff),
    np.absolute((0.8125-means2Allele_Diff_Normalized.Norm2AlleleDiff) /
                StError_2Allele_Norm.Norm2AlleleDiff),
    np.absolute((0.625-means2Allele_Diff_Normalized.Norm2AlleleDiff) /
                StError_2Allele_Norm.Norm2AlleleDiff)
]

#var_df = df_pair.groupby(['PairIndividuals'])[
#    'Norm2AlleleDiff'].apply(statistics.variance).to_frame()


df_to_print = pd.DataFrame(
    data={'PairIndividuals': means2Allele_Diff_Normalized.index.get_level_values(0).values, 'P0_mean': means2Allele_Diff_Normalized['Norm2AlleleDiff'],
           'Perc_Win_1stdeg_P0': means2Allele_Diff_Normalized['Perc_P0'], 'OverlapNSNPs': means2Allele_Diff_Normalized['OverlapNSNPs'], 'Nonnormalized_P0':means2Allele_Diff_Normalized['Nonnormalized_P0'], 'Nonnormalized_P0_serr':means2Allele_Diff_Normalized['Nonnormalized_P0_serr'], 'NSNPsXNorm': means2Allele_Diff_Normalized['NSNPsXNorm'] })


df_to_print['Rel'] = np.select(
    filters, values_Rel, default="IdenticalTwins/SameIndividual")
df_to_print['Zup'] = np.select(filters, values_Zup, default=np.abs(
    (0.625-means2Allele_Diff_Normalized.Norm2AlleleDiff)/StError_2Allele_Norm.Norm2AlleleDiff))
df_to_print['Zdown'] = np.select(filters, values_Zdown, default='NA')

filters_first_Deg = [
    (means2Allele_Diff_Normalized.Perc_P0 >= 0.6) & (
        df_to_print.Rel == "First Degree"),
    (means2Allele_Diff_Normalized.Perc_P0 >= 0.35) & (
        df_to_print.Rel == "First Degree"),
    (means2Allele_Diff_Normalized.Perc_P0 <= 0.3) & (
        df_to_print.Rel == "First Degree")
]
values_first_deg = [
    "NotApplicable", "Siblings", "Parent-offspring"
]

df_to_print['1st_Type'] = np.select(
    filters_first_Deg, values_first_deg, default="N/A")


df_to_print[['PairIndividuals', 'Rel', 'Zup', 'Zdown','P0_mean','Nonnormalized_P0', 'Nonnormalized_P0_serr', '1st_Type', 'Perc_Win_1stdeg_P0','OverlapNSNPs','NSNPsXNorm']].to_csv(
    'Read_Results.tsv', index=False, sep='\t')

fname = "./"+Arguments[0] + "_test*"
for filename in glob.glob(fname):
    os.remove(filename)
print("READ analysis finished. Please check Read_Results.tsv for results!")

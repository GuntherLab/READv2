#!/usr/bin/python
# PlinkIO


import numpy as np
import asyncio
from plinkio import plinkfile
import matplotlib.pyplot as plt
import pandas as pd
import sys
import os
import glob
from statistics import mean, median
from operator import index
from cmath import sqrt
import statistics
import time


def standard_error(x):
    return (statistics.stdev(x) / sqrt(len(x)))


Usage = """

Current READ version only supports .tped files.

READ.py <InputFile> <normalization> <normalization_value>

A normalization is only required when a user-defined values is used instead of the median, mean or maximum of the test population. All normalization settings are described below:

    median (default) - assuming that most pairs of compared individuals are unrelated, READ uses the median across all pairs for normalization.
    mean - READ uses the mean across all pairs for normalization, this would be more sensitive to outliers in the data (e.g. recent migrants or identical twins)
    max - READ uses the maximum across all pairs for normalization. This should be used to test trios where both parents are supposed to be unrelated but the offspring is a first degree relative to both others.
    value <val> - READ uses a user-defined value for normalization. This can be used if the value from another population should be used for normalization. That would be useful if only two individuals are available for the test population.

Optionally, one can add --window_size <value> at the end in order to modify the window size used (default: 1000000).

"""

print("    ===Thank you for using Relationship Estimation from Ancient DNA (READ)===    ", end='\n\n')


# turn of interactive plotting, does not show the plot unless plt.show() is used.
plt.ioff()

'''
def Ind_list():  # adds the individual pairs that exists in the sample in List_individiuals list
    for idx1, j in enumerate(list_all_individuals):
        for idx2, i in enumerate(list_all_individuals):
            # if same individuals or this comparison is done before: ignore.
            if (i == j or idx1 >= idx2):
                continue
            else:
                dictionary_pair_individuals["%s%s" % (j, i)] = 0
                dictionary_pair_individuals_missinginfo["%s%s" % (j, i)] = 0
                List_individuals.append("%s%s" % (j, i))
'''


def Pair_Calc(window):
    # print(window)
    missing = np.count_nonzero(window == -1)
    full_call = window.size
    # print("Window size:", window.size)
    # print("Missing: ", missing)
    full_call_var_corrected = full_call - missing
    # print("Full call var: ", full_call_var_corrected)
    if (full_call_var_corrected != 0):
        valor = np.count_nonzero(window == 1)
        IBS0 = np.count_nonzero(window == 0)
        # print("Full call corr:", full_call_var_corrected, "valor: ", valor, "IBS0: ", IBS0)
        P0 = float(IBS0)/full_call_var_corrected
        P2 = float(valor)/full_call_var_corrected
        # print("P0: ", P0, "P2:", P2)
        Line_info = [pair, valor, IBS0, P2, P0]
        pair_matrix.append(Line_info)


Arguments = sys.argv[1:]
norm_method = "median"
norm_value = ''
window_size = 100000
# dictionary_pair_individuals = {}
# dictionary_pair_individuals_missinginfo = {}
list_all_individuals = []
# list_2_all_individuals=[]
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
List_individuals = []

# print("Starting pairwise comparisons in windows", end='\n')


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

plink_file = plinkfile.open(Arguments[0])

sample_list = plink_file.get_samples()
# locus_list = plink_file.get_loci()
# print(locus_list.bp_position)
for sample in sample_list:
    list_all_individuals.append(sample.iid)


# Ind_list()
# del sample_list
# locus.bp_position -> base pair position, locus.allele* gives the alleles in that position.

start_time = time.time()
tmp_list = []
tmp_count = 0

# rows = np.fromiter(plink_file, int)
plink_file.transpose(Arguments[0] + "_test")
plink_file.close()
plink_file = plinkfile.open(Arguments[0] + "_test")
i = 0
# prep. np arrays for each individual
locus_list = plink_file.get_loci()
locus_pos = np.array(list(map(lambda d: d.bp_position, locus_list)))
locus_pos = np.floor_divide(locus_pos, window_size)
# print(locus_pos[-50:])
windows = np.where(locus_pos[:-1] != locus_pos[1:])[0]
# print(windows[-50:])
for row in plink_file:
    # row_array = np.array(row)
    # res_row = np.where(row_array == 3, np.nan, row_array)
    ind_dict[list_all_individuals[i]] = np.array(row)
    i += 1

keys = ind_dict.keys()
for key in keys:
    for key2 in keys:
        if (not (key == key2) and not (key2+key in pair_dict.keys())):
            nan_indices = np.union1d(np.where(
                ind_dict[key] == 3), np.where(
                ind_dict[key2] == 3))
            # print(type(nan_indices), nan_indices.size)
            pair_dict[key+key2] = 1 * (ind_dict[key] == ind_dict[key2])
            pair_dict[key+key2][nan_indices] = -1
            # print(-1 in pair_dict[key+key2])  # float nans?

for key in pair_dict.keys():
    pair = key
    arr = pair_dict[pair]
    pre_ind = 0


print("PlinkIO Pairwise comparison with Numpy takes: --- %s seconds ---" %
      (time.time() - start_time))
pair_matrix = np.array(pair_matrix)  # is not size of (length,6) fix!
# print(pair_matrix[:15])


# print("Create matrix!")

start_time = time.time()
# try to write everything directly to dataframe in the beginning
df_pair = pd.DataFrame(pair_matrix, columns=[
                       'PairIndividuals', 'IBS2', 'IBS0', 'P1', 'P0'])

df_pair = df_pair.astype(dtype={
                         'PairIndividuals': str, 'IBS2': float, 'IBS0': float, 'P1': float, 'P0': float})

print("Dataframe creation takes: --- %s seconds ---" %
      (time.time() - start_time))

# print("Dataframe created.")
start_time = time.time()
# print("")
del locus_list, plink_file, pair_matrix, ind_dict, pair_dict
# sample_list=plink_file.get_samples()
# print("Normalization started", end='\n')

print("Deletion of unused variables takes: --- %s seconds ---" %
      (time.time() - start_time))

# print(df_pair['P0'][:15])
# mean_2Allele_Diff_perc = df_pair.groupby(['PairIndividuals'])['P0'].mean()
#df_pair['2AlleleDiffMeanPerc'] = df_pair.groupby(['PairIndividuals'])['P0'].transform('mean')

if (norm_method == "median"):
    norm_value = median(df_pair['P0'])
elif (norm_method == "mean"):
    norm_value = mean(df_pair['P0'])
elif (norm_method == "max"):
    norm_value = max(df_pair['P0'])

start_time = time.time()
df_pair['Norm2AlleleDiff'] = df_pair['P0'] / norm_value
# print("P0 Normalized")
means2Allele_Diff_Normalized = df_pair.groupby(
    ['PairIndividuals'])['Norm2AlleleDiff'].mean().to_frame()
df_P0 = df_pair.groupby(['PairIndividuals'])[
    'P0'].mean().to_frame()  # unnormalized P0
StError_2Allele_Norm = df_pair.groupby(['PairIndividuals'])[
    'Norm2AlleleDiff'].apply(standard_error).to_frame()  # normalized P0 standard error
means2Allele_Diff_Normalized['StError_2Allele_Norm'] = StError_2Allele_Norm['Norm2AlleleDiff']
Nonnormalized_P0_serr = df_pair.groupby(['PairIndividuals'])[
    'P0'].apply(standard_error).to_frame()  # nonNormalized P0 standard error
means2Allele_Diff_Normalized['Nonnormalized_P0'] = df_P0['P0']
means2Allele_Diff_Normalized['Nonnormalized_P0_serr'] = Nonnormalized_P0_serr["P0"]
df_P0['Error'] = Nonnormalized_P0_serr['P0'].values
df_P0 = df_P0.reset_index()

print("Normalization takes: --- %s seconds ---" %
      (time.time() - start_time))

means2Allele_Diff_Normalized = means2Allele_Diff_Normalized.astype(
    dtype={'Norm2AlleleDiff': np.float16, 'StError_2Allele_Norm': np.float16, 'Nonnormalized_P0': np.float16, 'Nonnormalized_P0_serr': np.float16})
means2Allele_Diff_Normalized[['Norm2AlleleDiff', 'StError_2Allele_Norm', 'Nonnormalized_P0',
                              'Nonnormalized_P0_serr']].to_csv("meansP0_AncientDNA_normalized_READv2", index=True, sep='\t')
# print("Plotting!")


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

print("Plotting takes: --- %s seconds ---" %
      (time.time() - start_time))

# print("Estimating degree of relationships")


# Paper mentions these cutoffs are expressed as multiples of the standard error of the mean (Z).

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

df_to_print = pd.DataFrame(
    data={'PairIndividuals': means2Allele_Diff_Normalized.index.get_level_values(0).values})

df_to_print['Rel'] = np.select(
    filters, values_Rel, default="IdenticalTwins/SameIndividual")
df_to_print['Zup'] = np.select(filters, values_Zup, default=np.abs(
    (0.625-means2Allele_Diff_Normalized.Norm2AlleleDiff)/StError_2Allele_Norm.Norm2AlleleDiff))
df_to_print['Zdown'] = np.select(filters, values_Zdown, default='NA')

df_to_print[['PairIndividuals', 'Rel', 'Zup', 'Zdown']].to_csv(
    'Read_Results.tsv', index=False, sep='\t')

print("Printing takes: --- %s seconds ---" %
      (time.time() - start_time))


# print('\n')
###
fname = "./"+Arguments[0] + "_test*"
for filename in glob.glob(fname):
    # print(filename)
    os.remove(filename)
print("READ analysis finished. Please check Read_Results.tsv for results!")

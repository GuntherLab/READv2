#!/usr/bin/python
# PlinkIO


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
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
    value <val> - READ uses a user-defined value for normalization. This can be used if the value from another population should be used for normalization. That would be useful if only two individuals are available for the test population. A value can be obtained from the NonNormalizedP0 column of the file meansP0_AncientDNA_normalized from a previous run of READ.

Optionally, one can add --window_size <value> at the end in order to modify the window size used (default: 1000000).

"""

print("   ===Thank you for using Relationship Estimation from Ancient DNA (READ)===", end='\n\n')


# turn of interactive plotting, does not show the plot unless plt.show() is used.
plt.ioff()


def Ind_list():  # adds the individual pairs that exists in the sample in List_individiuals list

    for l in InFile_tfam:  # open tfam file.
        split_line = l.split()
        # create the list of individuals
        list_all_individuals.append(split_line[1])
    # for i in range(1,Number_individuals+1):
    #	list_all_individuals.append('Ind%d' % i)
    # list_2_all_individuals.extend(list_all_individuals)
    for idx1, j in enumerate(list_all_individuals):
        for idx2, i in enumerate(list_all_individuals):
            # if same individuals or this comparison is done before: ignore.
            if (i == j or idx1 >= idx2):
                continue
            else:
                dictionary_pair_individuals["%s%s" % (j, i)] = 0
                dictionary_pair_individuals_missinginfo["%s%s" % (j, i)] = 0
                List_individuals.append("%s%s" % (j, i))


Arguments = sys.argv[1:]
norm_method = "median"
norm_value = ''
Filetype = 'TPED'
window_size = 100000
dictionary_pair_individuals = {}
dictionary_pair_individuals_missinginfo = {}
list_all_individuals = []
# list_2_all_individuals=[]
previous_window = 0
Missing_info = 0
IBS2 = 0
IBS0 = 0
Total_snps_window = 0
full_call_var = 0
full_call_var_corrected = 0
pair_matrix = []  # a dict. to hold pair snp info, window info etc.
List_individuals = []

print("Starting pairwise comparisons in windows", end='\n')


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


InFileName = Arguments[0] + '.tped'
try:
    InFile = open(InFileName, 'r')
except IOError:
    print("the file %s.tped does not exist" % Arguments[0])
    sys.exit(1)

InfileName_tfam = Arguments[0] + '.tfam'
try:
    InFile_tfam = open(InfileName_tfam, 'r')
except IOError:
    print("the file %s.tfam does not exist" % Arguments[0])
    sys.exit(1)


Ind_list()
#del sample_list
# locus.bp_position -> base pair position, locus.allele* gives the alleles in that position.
snp_count = 0

start_time = time.time()
for Line in InFile:  # for each line in .tped file
    snp_count += 1
    Line = Line.strip("\n")
    ElementList = Line.split()
    Chromosome = int(ElementList[0])
    Position = int(ElementList[3])
    Genotype = ElementList[4:]
    window_index = Position/window_size
    # if (snp_count%5000)==0:
    #	print("Current position: Chromosome %s, bp %s" %(Chromosome,Position))
    Genotype = "".join(Genotype)
    # Divide the genotype into dublets to acquire allele info.
    Alleles_individuals = [Genotype[i:i+2] for i in range(0, len(Genotype), 2)]
    if window_index == previous_window:  # checks if we are still in the same window
        full_call_var += 1  # number of SNPs in the window
        for idx1, i in enumerate(list_all_individuals):
            for idx2, j in enumerate(list_all_individuals):
                if (i == j or idx1 >= idx2):
                    continue
                # if any of the ind. have 00 allele in corresponding index
                elif (Alleles_individuals[idx2] == "00") or (Alleles_individuals[idx1] == "00"):
                    dictionary_pair_individuals_missinginfo["%s%s" % (
                        list_all_individuals[idx1], list_all_individuals[idx2])] += 1  # increase missing info by one
                # if they have the same allele
                elif Alleles_individuals[idx1] == Alleles_individuals[idx2]:
                    dictionary_pair_individuals["%s%s" % (
                        list_all_individuals[idx1], list_all_individuals[idx2])] += 1

    else:
        IBS2 = [dictionary_pair_individuals[pair]
                for pair in List_individuals]  # IBS2 --> two alleles are shared
        Missing_info = [dictionary_pair_individuals_missinginfo[pair]
                        for pair in List_individuals]
        # parallel iteration
        for valor, missing, pair in zip(IBS2, Missing_info, List_individuals):
            full_call_var_corrected = full_call_var - missing
            # corrected SNP number in the window - IBS2 of that pair?
            IBS0 = float(full_call_var_corrected - valor)
            if (full_call_var_corrected != 0):
                # P1 referred as P2 here?
                P2 = float(valor)/full_call_var_corrected
                P0 = float(IBS0/full_call_var_corrected)
                Line_info = [pair, Chromosome, valor, IBS0, P2, P0]
                pair_matrix.extend(Line_info)
                # if(pair not in pair_dict.keys()):
                # pair_dict[pair]=[Line_info]
                # else:
                # pair_dict[pair].append(Line_info)

        # reinitialize the occurances of "two alleles shared" and missingness to zero
        for key in dictionary_pair_individuals:
            dictionary_pair_individuals[key] = 0
            dictionary_pair_individuals_missinginfo[key] = 0

        full_call_var = 1
        previous_window = window_index

        for idx1, i in enumerate(list_all_individuals):
            for idx2, j in enumerate(list_all_individuals):
                if (i == j or idx1 >= idx2):
                    continue
                elif (Alleles_individuals[idx2] == "00") or (Alleles_individuals[idx1] == "00"):
                    dictionary_pair_individuals_missinginfo["%s%s" % (
                        list_all_individuals[idx1], list_all_individuals[idx2])] += 1
                elif Alleles_individuals[idx1] == Alleles_individuals[idx2]:
                    dictionary_pair_individuals["%s%s" % (
                        list_all_individuals[idx1], list_all_individuals[idx2])] += 1

print("No Plink IO Pairwise comparison takes: --- %s seconds ---" %
      (time.time() - start_time))

'''
df_pair = pd.DataFrame(pair_matrix, columns=[
                       'PairIndividuals', 'Chromosome', 'IBS2', 'IBS0', 'P1', 'P0'])
print("Dataframe created.")
del pair_matrix
#mean_2Allele_Diff_perc = df_pair.groupby(['PairIndividuals'])['P0'].mean()
df_pair['2AlleleDiffMeanPerc'] = df_pair.groupby(
    ['PairIndividuals'])['P0'].transform('mean')

if (norm_method == "median"):
    norm_value = median(df_pair['2AlleleDiffMeanPerc'])
elif (norm_method == "mean"):
    norm_value = mean(df_pair['2AlleleDiffMeanPerc'])
elif (norm_method == "max"):
    norm_value = max(df_pair['2AlleleDiffMeanPerc'])

df_pair['Norm2AlleleDiff'] = df_pair['P0'] / norm_value
print("P0 Normalized")
means2Allele_Diff_Normalized = df_pair.groupby(
    ['PairIndividuals'])['Norm2AlleleDiff'].mean().to_frame()
df_P0 = df_pair.groupby(['PairIndividuals'])['P0'].mean().to_frame()
StError_2Allele_Norm = df_pair.groupby(['PairIndividuals'])[
    'Norm2AlleleDiff'].apply(standard_error).to_frame()
Nonnormalized_P0_serr = df_pair.groupby(['PairIndividuals'])[
    'P0'].apply(standard_error).to_frame()
df_P0['Error'] = Nonnormalized_P0_serr['P0'].values
df_P0 = df_P0.reset_index()

print("Plotting!")

#df_fig=df_pair[['PairIndividuals', 'P0']]
print("a")
# df_fig=df_fig.set_index('PairIndividuals').T
print("b")
#boxplot=df_P0.reset_index().boxplot(by = 'PairIndividuals',rot=90, figsize=((8+0.1*(len(df_P0.index)),6)))
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
print("????")


print("Estimating degree of relationships")

filters = [
    means2Allele_Diff_Normalized.Norm2AlleleDiff >= 0.90625,
    means2Allele_Diff_Normalized.Norm2AlleleDiff >= 0.8125,
    means2Allele_Diff_Normalized.Norm2AlleleDiff >= 0.625
]
values_Rel = ["Unrelated", "Second Degree", "First Degree"]
values_Zup = [
    "NA",
    (0.90625-pd.to_numeric(means2Allele_Diff_Normalized.Norm2AlleleDiff)) /
    pd.to_numeric(StError_2Allele_Norm.Norm2AlleleDiff),
    (0.8125-pd.to_numeric(means2Allele_Diff_Normalized.Norm2AlleleDiff)) /
    pd.to_numeric(StError_2Allele_Norm.Norm2AlleleDiff),
]
values_Zdown = [
    (0.90625-pd.to_numeric(means2Allele_Diff_Normalized.Norm2AlleleDiff)) /
    pd.to_numeric(StError_2Allele_Norm.Norm2AlleleDiff),
    (0.8125-pd.to_numeric(means2Allele_Diff_Normalized.Norm2AlleleDiff)) /
    pd.to_numeric(StError_2Allele_Norm.Norm2AlleleDiff),
    (0.625-pd.to_numeric(means2Allele_Diff_Normalized.Norm2AlleleDiff)) /
    pd.to_numeric(StError_2Allele_Norm.Norm2AlleleDiff)
]

df_to_print = pd.DataFrame(
    data={'PairIndividuals': means2Allele_Diff_Normalized.index.get_level_values(0).values})

df_to_print['Rel'] = np.select(
    filters, values_Rel, default="IdenticalTwins/SameIndividual")
df_to_print['Zup'] = np.select(filters, values_Zup, default=(
    0.625-pd.to_numeric(means2Allele_Diff_Normalized.Norm2AlleleDiff))/pd.to_numeric(StError_2Allele_Norm.Norm2AlleleDiff))
df_to_print['Zdown'] = np.select(filters, values_Zdown, default='NA')

df_to_print[['Zup', 'Zdown']] = df_to_print[[
    'Zup', 'Zdown']].apply(pd.to_numeric, errors='coerce')
df_to_print[['PairIndividuals', 'Rel', 'Zup', 'Zdown']].to_csv(
    'Read_Results.tsv', index=False, sep='\t')


print('\n)')
###

print("READ analysis finished. Please check READ_results for results!")
'''

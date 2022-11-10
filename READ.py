#!/usr/bin/python
#PlinkIO


import statistics


def standard_error(x):
	return( statistics.stdev(x) / sqrt(len(x)))


Usage =""" 

Current READ version only supports .tped files.

READ.py <InputFile> <normalization> <normalization_value>

A normalization is only required when a user-defined values is used instead of the median, mean or maximum of the test population. All normalization settings are described below:

    median (default) - assuming that most pairs of compared individuals are unrelated, READ uses the median across all pairs for normalization.
    mean - READ uses the mean across all pairs for normalization, this would be more sensitive to outliers in the data (e.g. recent migrants or identical twins)
    max - READ uses the maximum across all pairs for normalization. This should be used to test trios where both parents are supposed to be unrelated but the offspring is a first degree relative to both others.
    value <val> - READ uses a user-defined value for normalization. This can be used if the value from another population should be used for normalization. That would be useful if only two individuals are available for the test population. A value can be obtained from the NonNormalizedP0 column of the file meansP0_AncientDNA_normalized from a previous run of READ.

Optionally, one can add --window_size <value> at the end in order to modify the window size used (default: 1000000).

"""

print("   ===Thank you for using Relationship Estimation from Ancient DNA (READ)===" , end='\n\n')


from cmath import sqrt
from operator import index
from statistics import mean, median, stdev
import sys
import pandas as pd
import numpy as np


path_Rscript='/usr/bin/Rscript'
Arguments = sys.argv[1:]
norm_method="median"
norm_value=''
Filetype='TPED'
window_size=100000


def Ind_list(): #adds the individual pairs that exists in the sample in List_individiuals list

	for l in InFile_tfam: #open tfam file.
		split_line=l.split()
		list_all_individuals.append(split_line[1]) #create the list of individuals
	#for i in range(1,Number_individuals+1):
	#	list_all_individuals.append('Ind%d' % i)
	#list_2_all_individuals.extend(list_all_individuals)
	for idx1, j in enumerate(list_all_individuals):
		for idx2, i in enumerate(list_all_individuals):
			if (i == j or idx1 >= idx2): #if same individuals or this comparison is done before: ignore.
				continue
			else:
				dictionary_pair_individuals["%s%s" % (j,i)]=0
				dictionary_pair_individuals_missinginfo["%s%s" % (j,i)]=0
				List_individuals.append("%s%s" % (j,i))


#Input check!
if(len(Arguments)==0):
	sys.exit(Usage)

if(len(Arguments)>=2 and Arguments[1] in ["median","mean","max","value"]):
	norm_method=Arguments[1]
elif(len(Arguments)>=2 and Arguments[1] not in ["median","mean","max","value"]):
	print("No valid standardization method specified! Using the default method.")
if(len(Arguments)>=3 and Arguments[1]=="value"):
	norm_value=float(Arguments[2])
if("--window_size" in Arguments):
	ws_index=Arguments.index("--window_size")
	if(len(Arguments)>(ws_index+1)):
		window_size=int(Arguments[ws_index+1])
	else:
		sys.exit("No window size specified!")
	
	
	
InFileName= Arguments[0] + '.tped'
try:
	InFile = open(InFileName,'r')
except IOError:
	print("the file %s.tped does not exist" % Arguments[0])
	sys.exit(1)

InfileName_tfam=Arguments[0] + '.tfam'
try:
	InFile_tfam = open(InfileName_tfam,'r')
except IOError:
	print("the file %s.tfam does not exist" % Arguments[0])
	sys.exit(1)




dictionary_pair_individuals={}
dictionary_pair_individuals_missinginfo={}
list_all_individuals=[]
#list_2_all_individuals=[]
previous_window=0
Missing_info=0
IBS2=0
IBS0=0
Total_snps_window=0
full_call_var=0
full_call_var_corrected=0
pair_dict={} #a dict. to hold pair snp info, window info etc.
List_individuals=[]

print("Starting pairwise comparisons in windows", end='\n')

Ind_list()
previous_window=0
snp_count=0

for Line in InFile: #for each line in .tped file
	snp_count+=1
	Line = Line.strip("\n")
	ElementList = Line.split()
	Chromosome = int(ElementList[0])
	Position = int(ElementList[3])
	Genotype= ElementList[4:]
	window_index=Position/window_size 

	'''
	if (snp_count%5000)==0:
		print("Current position: Chromosome %s, bp %s" %(Chromosome,Position))
	'''
	Genotype="".join(Genotype)
	#Divide the genotype into dublets to acquire allele info.
	Alleles_individuals = [Genotype[i:i+2] for i in range(0, len(Genotype), 2)]
	if window_index == previous_window:  #checks if we are still in the same window
		full_call_var+=1 #number of SNPs in the window
		for idx1, i in enumerate(list_all_individuals):
			for idx2, j in enumerate(list_all_individuals): 
					if (i == j or idx1 >= idx2): 
						continue
					elif (Alleles_individuals[idx2]=="00") or (Alleles_individuals[idx1]=="00"): #if any of the ind. have 00 allele in corresponding index
						dictionary_pair_individuals_missinginfo["%s%s" % (list_all_individuals[idx1],list_all_individuals[idx2])]+=1 #increase missing info by one
					elif Alleles_individuals[idx1] == Alleles_individuals[idx2]: #if they have the same allele
						dictionary_pair_individuals["%s%s" % (list_all_individuals[idx1],list_all_individuals[idx2])]+=1

	else: 
		IBS2= [dictionary_pair_individuals[pair] for pair in List_individuals] #IBS2 --> two alleles are shared
		Missing_info = [dictionary_pair_individuals_missinginfo[pair] for pair in List_individuals]
		for valor, missing, pair in zip(IBS2,Missing_info,List_individuals): #parallel iteration
			full_call_var_corrected = full_call_var - missing
			IBS0= float(full_call_var_corrected - valor) #corrected SNP number in the window - IBS2 of that pair?
			if (full_call_var_corrected!=0):
				P2 = float(valor)/full_call_var_corrected #P1 referred as P2 here?
				P0 = float(IBS0/full_call_var_corrected)
				Line_info= [Chromosome, window_index, full_call_var, valor, IBS0,P2,P0,missing,full_call_var_corrected]
				if(pair not in pair_dict.keys()):
					pair_dict[pair]=[Line_info]
				else:
					pair_dict[pair].append(Line_info)

		
		#reinitialize the occurances of "two alleles shared" and missingness to zero
		for key in dictionary_pair_individuals:
			dictionary_pair_individuals[key]=0	
		for key in dictionary_pair_individuals_missinginfo:
			dictionary_pair_individuals_missinginfo[key]=0
		full_call_var=1
		previous_window= window_index	
			
		for idx1, i in enumerate(list_all_individuals):
			for idx2, j in enumerate(list_all_individuals):
				if (i == j or idx1 >= idx2):
					continue
				elif (Alleles_individuals[idx2]=="00") or (Alleles_individuals[idx1]=="00"):
					dictionary_pair_individuals_missinginfo["%s%s" % (list_all_individuals[idx1],list_all_individuals[idx2])]+=1
				elif Alleles_individuals[idx1] == Alleles_individuals[idx2]:
					dictionary_pair_individuals["%s%s" % (list_all_individuals[idx1],list_all_individuals[idx2])]+=1



InFile.close()



print("")

print("Normalization started", end='\n')

pair_matrix=[]
for key in pair_dict.keys():
	for line in pair_dict[key]:
		pair_matrix.append([key]+line)
del pair_dict,dictionary_pair_individuals, dictionary_pair_individuals_missinginfo

df_pair = pd.DataFrame(pair_matrix, columns= ['PairIndividuals','Chromosome','WindowIndex','SNVperWindow','IBS2','IBS0','P1','P0','Missing', 'SNVperWindowCorrected'])
del pair_matrix
#mean_2Allele_Diff_perc = df_pair.groupby(['PairIndividuals'])['P0'].mean()
df_pair['2AlleleDiffMeanPerc']=df_pair.groupby(['PairIndividuals'])['P0'].transform('mean')

if(norm_method == "median"):
	norm_value = median(df_pair['2AlleleDiffMeanPerc'])
elif(norm_method == "mean"):
	norm_value = mean(df_pair['2AlleleDiffMeanPerc'])
elif(norm_method == "max"):
	norm_value = max(df_pair['2AlleleDiffMeanPerc'])

df_pair['Norm2AlleleDiff'] = df_pair['P0'] / norm_value

#does not work as intended like this, try to convert it back.
means2Allele_Diff_Normalized = df_pair.groupby(['PairIndividuals'])['Norm2AlleleDiff'].mean().to_frame()
#df_pair['Norm2AlleleDiffMeans']=df_pair.groupby(['PairIndividuals'])['Norm2AlleleDiff'].transform('mean')
#df_pair['StandardError']=df_pair.groupby(['PairIndividuals'])['Norm2AlleleDiff'].transform(standard_error)
StError_2Allele_Norm = df_pair.groupby(['PairIndividuals'])['Norm2AlleleDiff'].apply(standard_error).to_frame()
#df_pair["NonNormalizedP0StandardError"]=df_pair.groupby(['PairIndividuals'])['P0'].transform(standard_error)
Nonnormalized_P0_serr = df_pair.groupby(['PairIndividuals'])['P0'].apply(standard_error).to_frame()




'''
if(len(df_pair['2AlleleDiffMeanPerc']) % 2 == 1):
	seb2_nn = df_pair[[df_pair['2AlleleDiffMeanPerc'] == norm_value]]['Norm2AlleleDiffMeans']  #get the index in mean_2Allele_Diiff_perc where P0 val. is equal to normalization value
	
else:
	sorted_nn_df = df_pair.sort_values('NonNormalizedP0StandardError')
	sorted_nn = sorted_nn_df['NonNormalizedP0StandardError']
	se1_nn = sorted_nn[len(df_pair['2AlleleDiffMeanPerc'])/2]
	se2_nn = sorted_nn[1+len(df_pair['2AlleleDiffMeanPerc']/2)]
	seb2_nn = mean([se1_nn,se2_nn])
'''

#indices=which(abs(means_2AlleleDifference_percentage$P0-normalization_val)== min(abs(means_2AlleleDifference_percentage$P0-normalization_val),na.rm=TRUE))


print("Estimating degree of relationships")

###Predict relationships

#OutFileName="READ_results"
#OutFile=open(OutFileName, 'w')

#OutFile.write("PairIndividuals\tRelationship\tZ_upper\tZ_lower\n")



#Paper mentions these cutoffs are expressed as multiples of the standard error of the mean (Z).

filters=[
	means2Allele_Diff_Normalized.Norm2AlleleDiff >= 0.90625,
	means2Allele_Diff_Normalized.Norm2AlleleDiff >= 0.8125,
	means2Allele_Diff_Normalized.Norm2AlleleDiff >= 0.625
]
values_Rel=[ "Unrelated", "Second Degree", "First Degree" ]
values_Zup=[
	"NA", 
	(0.90625-means2Allele_Diff_Normalized.Norm2AlleleDiff)/StError_2Allele_Norm.Norm2AlleleDiff,
	(0.8125-means2Allele_Diff_Normalized.Norm2AlleleDiff)/StError_2Allele_Norm.Norm2AlleleDiff,
	]
values_Zdown=[
	(0.90625-means2Allele_Diff_Normalized.Norm2AlleleDiff)/StError_2Allele_Norm.Norm2AlleleDiff,
	(0.8125-means2Allele_Diff_Normalized.Norm2AlleleDiff)/StError_2Allele_Norm.Norm2AlleleDiff,
	(0.625-means2Allele_Diff_Normalized.Norm2AlleleDiff)/StError_2Allele_Norm.Norm2AlleleDiff
]

df_to_print = pd.DataFrame(data = {'PairIndividuals':means2Allele_Diff_Normalized.index.get_level_values(0).values})

df_to_print['Rel']=np.select(filters, values_Rel, default="IdenticalTwins/SameIndividual")
df_to_print['Zup']=np.select(filters, values_Zup, default=(0.625-means2Allele_Diff_Normalized.Norm2AlleleDiff)/StError_2Allele_Norm.Norm2AlleleDiff)
df_to_print['Zdown']=np.select(filters, values_Zdown, default = 'NA')

df_to_print[['PairIndividuals', 'Rel', 'Zup', 'Zdown']].to_csv('Read_Results.tsv', index=False, sep='\t')

'''
for i in df_pair.index.get_level_values(0).values: #Printing results
	P0 = df_pair.iloc[i]['Norm2AlleleDiff']
	pair= df_pair.iloc[i]['PairIndividuals']
	if (P0 >= 0.90625):
		Relationship="Unrelated"
		Zup='NA'
		Zdown=(0.90625-P0)
	elif (P0 >= 0.8125):
		Relationship="Second Degree"
		Zup=(0.90625-P0)
		Zdown=(0.8125-P0)
	elif (P0 >= 0.625):
		Relationship="First Degree"
		Zup=(0.8125-P0)
		Zdown=(0.625-P0)
	else:
		Relationship="IdenticalTwins/SameIndividual"
		Zup=(0.625-P0)
		Zdown='NA'
	
	OutString="%s\t%s\t%s\t%s" % (pair,Relationship,Zup,Zdown)
	OutFile.write(OutString+ '\n')
'''

#OutFile.close()

print('\n)')
###
print("READ analysis finished. Please check READ_results for results!")
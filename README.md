# Relationship Estimation from Ancient DNA (READv2) #

# PLEASE NOTE: This repository has been set to public to have a version of record for Erkin Alaçamlı's MSc thesis. READv2 development is still ongoing and we do not recommend to use it at this point. The README below might also be incomplete. #

## Requirements ##

* Python 3.7 or higher
* Python pandas library
* Python NumPy library
* PLINKIO library
* Unix-like operating system

To make sure all these dependencies are available, we suggest to set up a conda environment for READv2. The required packages and libraries can be downloaded and installed with the following conda command:

    conda create -n readv2 python=3.9 pandas=1.3.1 numpy=1.21.1 pip=22.3.1
and the environment can be activated by the

    conda activate readv2
and PLINKIO library can be installed with the following pip command afterwards:

    pip install plinkio
This way, the PLINKIO library will be added to the environment, and READv2 will be ready to use.

Quit the conda environment with

    conda deactivate
  

## How to use READv2? ##

### Running READv2 ###

READv2 assumes pseudohaploid data which means that all individuals in the input file are coded as completely homozygous. It is a quite common approach in aDNA research to obtain such data by randomly sampling a sequencing read per SNP site (using a pre-defined panel of SNP sites). Implementations of this approach can be found in [ANGSD](http://www.popgen.dk/angsd/index.php/Haploid_calling) and [SequenceTools](https://github.com/stschiff/sequenceTools). Common file formats for genotype data are VCF, Eigenstrat or various Plink formats. [Plink](http://pngu.mgh.harvard.edu/~purcell/plink/], [vcftools](http://vcftools.sourceforge.net/man_latest.html], [PGDSpider](http://www.cmpg.unibe.ch/software/PGDSpider/) or [convertf](https://github.com/argriffing/eigensoft/tree/master/CONVERTF) can be used to convert between these formats.

The input for READ is a pair of files in Plink's TPED/TFAM format containing information on the individuals and their genotypes. Please note that you should use a single file pair per test population (you can use Plink's --keep to restrict the individuals to your population of interest). Plink can be used to filter and manipulate files. Make sure to exclude non-polymorphic and low frequency variants (e.g. --maf 0.01).

Assume you have downloaded READv2 and a group of files example.bed, example.bim and example.fam in your current directory, then you can simply run READv2 like this:

    python READ2.py -i example 

This runs the READv2 script in default settings. The results of your READv2 analysis can be found in the two files Read_Results.tsv and meansP0_AncientDNA_normalized. The main result file is Read_Results.tsv. It contains a number of columns: 
 * the pair of individuals, the predicted relationship
 * two columns showing how many standard errors that normalized mean P0 score is from a higher class of genetic difference (Z_upper, third column = distance to lower degrees of relationship) and a lower class of genetic difference (Z_lower, fourth column = distance to higher degrees of relationship) These values can be used to assess the certainty of the classification. We observed in our simulations that false classifications were enriched when the normalized mean P0 score were less than one standard error from the nearest threshold (i.e. |Z|<1).
 * the normalized mean P0 score for the pair  
 * the non-normalized mean P0 score for the pair
 * the non-normalized standard error of the mean P0
 * for 1st degree relatives whether they are siblings or parent-offspring (if sufficient data is included)
 * The percentage of 20 Mbp windows not classified as first or second degree.
 * The number of overlapping SNPs per pair
 * The number of overlapping SNPs per pair (*number of overlapping SNPs* times *the normalization value assumed to represent an unrelated pair for the population*)
 * The kinship coefficient theta (1 - *normalized mean P0*)

Additionally, a graphical representation of the results is produced (READ_results_plot.pdf) showing the results as well as uncertainties of individual estimates (plots are only produced for less than 1000 pairs of individuals). meansP0_AncientDNA_normalized is mainly for READ's internal use but it can be used for normalization with a user-defined value (see below).

#### Command line options ####

Options when running READv2:

* `-i`, `--input_file` *val* -- Input file prefix (required). The current READ version only supports Plink bed/bim/fam files.
* `-n`, `--norm_method` *val* -- Normalization method (either 'mean', 'median', 'max' or 'value').
   * `median` (default) -- assuming that most pairs of compared individuals are unrelated, READ uses the median across all pairs for normalization.
   * `mean` -- READ uses the mean across all pairs for normalization, this would be more sensitive to outliers in the data (e.g. recent migrants or identical twins)
   * `max` -- READ uses the maximum across all pairs for normalization. This should be used to test trios where both parents are supposed to be unrelated but the offspring is a first degree relative to both others.
   * `value` -- READ uses a user-defined value for normalization. This can be used if the value from another population should be used for normalization. That would be useful if only two individuals are available for the test population. Normalization value needs to be provided through --norm_value
* `--norm_value` *val* -- Provide a user-defined normalization value
* `--window_size` *val* -- Change window size for block jackknife or for window-based P0 estimates (as in READv1), default: 5000000
* `--window_est` -- Window based estimate of P0 (as opposed to the genome-wide estimate, default in READv2)
* `-h`, `--help` -- Print help message
* `-v`, `--version` -- Print version

#### Normalization ####

To scale the classification cutoffs for the expected distances between two unrelated individuals from the same population and using the same SNP panel, READ performs a normalization step using an estimate for this pairwise distance of unrelated individual. This helps to overcome potential biases arising from e.g. general population diversity or SNP ascertainment. Obtaining a good estimate is cruicial for the accuracy of READ but it is also difficult to formulate these expectations as large population samples are usually not available in aDNA studies. In default setting, READ estimates the normalization value from the input data which worked quite well in our test runs for sample sizes of 4 or more. READ calculates all pairwise distances between individuals. Based on the assumption that most pairs of individuals in a sample are unrelated, READ then uses the median of all pairwise differences as an estimate of the expected distance between two unrelated individuals.

This default setting (median of the sample -- **median**) may not be the best choice in all cases. Some studies may expect that the majority of all pairs are in fact related (e.g. when analyzing parents-child trios or larger pedigrees) or the sample size may be too small (e.g. 2). For trios, using the maximum value (**max**) of all three pairs may be a sensible choice as two of the three pairs would be expected to be first degrees (parent-child) while only one pairwise comparison (mother-father) would represent unrelated individuals. In other cases, such as sample sizes of only two, a normalization value can be obtained from a second (but genetically similar) population. To use this last approach, READ needs to be run twice: first on a dataset containing the genotype information from the population used to estimate the normalization value and second on the actual test population. For the first run, READ can be run in default settings. READ will output a file called meansP0_AncientDNA_normalized. The user can obtain a value for normalization from the NonNormalizedP0 column (e.g. the median). For the test population, READv2 can then be run with `--norm_value` to input the user-defined normalization value.

#### Window-based versus genome-wide approach ####

[READv1](https://bitbucket.org/tguenther/read/) used a default window size of 1000000 bp which means that the genome was split into non-overlapping windows of that size. P0 was then calculated per window and the genome-wide average across all windows was used as test statistic for classification. Standard errors were calculated as the standard error of that mean. During the development of READv2, we compared performance for different window sizes as well as a genome-wide estimate (i.e. calculating a single value across all sites without splitting the genome into windows). We observed that the genome-wide approach performed slightly better than the different window sizes which is why we decided to make it default for READv2. To estimate standard errors, we use a block-jackknife approach and blocks of 5000000 bp as commonly used in human population genetics. This block size can be changes with `--window_size`. Users who with to use the window-based approach can do so by using the command line options `--window_est` and `--window_size`.

Note that for the differentiation between parent-offspring and siblings, we use a different window size of 20000000 bp.

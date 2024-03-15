
import pandas as pd
import matplotlib.pyplot as plt
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen # matrix generation
from SigProfilerAssignment import Analyzer as Analyze # signature analysis
import sigProfilerPlotting as sigPlt # signatures plots
import glob

# Install reference genome
#from SigProfilerMatrixGenerator import install as genInstall
#genInstall.install('GRCh37', bash=True)




## High enhancers ###
# Extract variants of interest for "high" CtIP and GRHL clusters - SNVs
abel = "high"
clust_high_ctip = ["CtIP_cluster_2_Enh", "CtIP_cluster_3_Enh", "CtIP_cluster_5_Enh"]
clust_high_grhl = ["GRHL_cluster_1_Enh", "GRHL_cluster_2_Enh", "GRHL_cluster_4_Enh"]

SSMs_ctip = pd.read_csv("/Users/ieo6983/Desktop/fragile_enhancer_clinical/results/ICGC/enhancers_SSMs_overlaps/data/Table_enh_SSMs_CtIP.all_overlaps.tsv", sep = "\t") 
SSMs_ctip = SSMs_ctip[SSMs_ctip['cluster'].isin(clust_high_ctip) & SSMs_ctip['width_sbj'] == 1]
SSMs_grhl = pd.read_csv("/Users/ieo6983/Desktop/fragile_enhancer_clinical/results/ICGC/enhancers_SSMs_overlaps/data/Table_enh_SSMs_GRHL.all_overlaps.tsv", sep = "\t") 
SSMs_grhl = SSMs_grhl[SSMs_grhl['cluster'].isin(clust_high_grhl) & SSMs_grhl['width_sbj'] == 1]

# Convert to appropriate format - Info on input file txt format: https://osf.io/xfphr
format_ssms_ctip = {
    'Project' : SSMs_ctip['project_code'],
    'Sample' : SSMs_ctip['icgc_sample_id'],
    'ID' : SSMs_ctip['icgc_donor_id'],
    'Genome' : "GRCh37",
    'mut_type' : "SNP", 
    'chrom' : SSMs_ctip['seqnames_sbj'], 
    'pos_start' : SSMs_ctip['start_sbj'], 
    'pos_end' : SSMs_ctip['end_sbj'], 
    'ref' : SSMs_ctip['reference_genome_allele'], 
    'alt' : SSMs_ctip['mutated_to_allele'], 
    'Type' : "SOMATIC"
}
format_ssms_ctip = pd.DataFrame(format_ssms_ctip)

format_ssms_grhl = {
    'Project' : SSMs_grhl['project_code'],
    'Sample' : SSMs_grhl['icgc_sample_id'],
    'ID' : SSMs_grhl['icgc_donor_id'],
    'Genome' : "GRCh37",
    'mut_type' : "SNP", 
    'chrom' : SSMs_grhl['seqnames_sbj'], 
    'pos_start' : SSMs_grhl['start_sbj'], 
    'pos_end' : SSMs_grhl['end_sbj'], 
    'ref' : SSMs_grhl['reference_genome_allele'], 
    'alt' : SSMs_grhl['mutated_to_allele'], 
    'Type' : "SOMATIC"
}
format_ssms_grhl = pd.DataFrame(format_ssms_grhl)

# Save formatted file - required for matrix generation. The folder with input file will correspond to the output folder
OUT_FOLDER = "/Users/ieo6983/Desktop/fragile_enhancer_clinical/results/ICGC/mutational_signature_analysis/"
#format_ssms_ctip.to_csv(OUT_FOLDER+'SNVs_CtIP_{}/input.SNVs_ctip_{}.txt'.format(label), index=None, sep = '\t')
#format_ssms_grhl.to_csv(OUT_FOLDER+'SNVs_GRHL_{}/input.SNVs_grhl_{}.txt'.format(label), index=None, sep = '\t')

# Save version with only 1 sample name
format_ssms_ctip['Sample'] = "ONE_SAMPLE"
format_ssms_grhl['Sample'] = "ONE_SAMPLE"
#format_ssms_ctip.to_csv(OUT_FOLDER+'SNVs_CtIP_{}/input.SNVs_ctip_{}.ONE_SAMPLE.txt'.format(label), index=None, sep = '\t')
#format_ssms_grhl.to_csv(OUT_FOLDER+'SNVs_GRHL_{}/input.SNVs_grhl_{}.ONE_SAMPLE.txt'.format(label), index=None, sep = '\t')

# Save version with only Samples = Patients 
format_ssms_ctip['Sample'] = SSMs_ctip['icgc_donor_id']
format_ssms_grhl['Sample'] = SSMs_grhl['icgc_donor_id']
#format_ssms_ctip.to_csv(OUT_FOLDER+'SNVs_CtIP_{}.Patients/input.SNVs_ctip_{}.Patients.txt'.format(label), index=None, sep = '\t')
#format_ssms_grhl.to_csv(OUT_FOLDER+'SNVs_GRHL_{}.Patients/input.SNVs_grhl_{}.Patients.txt'.format(label), index=None, sep = '\t')


## Generate mutational matrix - For each sample
path_to_file = "/Users/ieo6983/Desktop/fragile_enhancer_clinical/results/ICGC/mutational_signature_analysis/" # Full path of the saved input files in the desired output folder.

matrices = matGen.SigProfilerMatrixGeneratorFunc("SNVs_CtIP_{}".format(label), "GRCh37", path_to_file+'SNVs_CtIP_{}/'.format(label), 
                                                 exome=False, bed_file=None, chrom_based=False, 
                                                 plot=True, tsb_stat=False, seqInfo=True)
matrices2 = matGen.SigProfilerMatrixGeneratorFunc("SNVs_GRHL_{}".format(label), "GRCh37", path_to_file+'SNVs_GRHL_{}/'.format(label), 
                                                 exome=False, bed_file=None, chrom_based=False, 
                                                 plot=True, tsb_stat=False, seqInfo=True)

## Generate mutational matrix - For ONE_SAMPLE only (unified sample contaning all of them)
path_to_file = "/Users/ieo6983/Desktop/fragile_enhancer_clinical/results/ICGC/mutational_signature_analysis/" # Full path of the saved input files in the desired output folder.

matrices = matGen.SigProfilerMatrixGeneratorFunc("SNVs_CtIP_{}".format(label), "GRCh37", path_to_file+'SNVs_CtIP_{}.ONE_SAMPLE/'.format(label), 
                                                 exome=False, bed_file=None, chrom_based=False, 
                                                 plot=True, tsb_stat=False, seqInfo=True)
matrices2 = matGen.SigProfilerMatrixGeneratorFunc("SNVs_GRHL_{}".format(label), "GRCh37", path_to_file+'SNVs_GRHL_{}.ONE_SAMPLE/'.format(label), 
                                                 exome=False, bed_file=None, chrom_based=False, 
                                                 plot=True, tsb_stat=False, seqInfo=True)

## Generate mutational matrix - For each icgc_donor_id
path_to_file = "/Users/ieo6983/Desktop/fragile_enhancer_clinical/results/ICGC/mutational_signature_analysis/" # Full path of the saved input files in the desired output folder.

matrices = matGen.SigProfilerMatrixGeneratorFunc("SNVs_CtIP_{}".format(label), "GRCh37", path_to_file+'SNVs_CtIP_{}.Patients/'.format(label), 
                                                 exome=False, bed_file=None, chrom_based=False, 
                                                 plot=True, tsb_stat=False, seqInfo=True)
matrices2 = matGen.SigProfilerMatrixGeneratorFunc("SNVs_GRHL_{}".format(label), "GRCh37", path_to_file+'SNVs_GRHL_{}.Patients/'.format(label), 
                                                 exome=False, bed_file=None, chrom_based=False, 
                                                 plot=True, tsb_stat=False, seqInfo=True)

'''
Actually, for each donor, we have only one sample. Hence in this case the first and third analyses correspond.
'''



## Plot SBSs frequencies for all samples
'''
SBS96 files represnt the frequency of each mutation x sample (C>A, C>G,...,T>G), with the trinucleotide contexts of each mutation.
- SNVs_CtIP_high.SBS96.all : Single Base Substitution-96 (SBS-96) mutational matrix 
- SBS_96_plots_SNVs_GRHL_high.pdf : Single Base Substitution-96 (SBS-96) numerical plot graphs for each sample

For other plots, such as plot of mutational frequencies for all samples, refer to: https://osf.io/2aj6t/wiki/3.%20Plotting%20Substitutions/.
'''

matrix_path_SBS = path_to_file+'SNVs_CtIP_high/output/SBS/SNVs_CtIP_high.SBS96.all'
path_to_plots = path_to_file+'SNVs_CtIP_high/output/plots/SBS_plots/'

sigPlt.plotSBS(matrix_path_SBS, path_to_plots, "SNVs_CtIP_high", "96", )




# Use SigProfilerAssignment's main function for mutational signatures assignment
# COSMIC v3.3 mutational signatures are used by default as the input reference signatures
MARKERS = ["CtIP", "GRHL"]

file_pattern = "*SBS96.all"
for marker in MARKERS:
    print(marker)

    dir_matrices = path_to_file + "SNVs_{}_{}/output/SBS/".format(marker, label)        
    file_list = glob.glob(dir_matrices+file_pattern)[0]

    Analyze.cosmic_fit(samples=file_list, # path to the input somatic mutations file (if using segmentation file/mutational matrix)
                   output=path_to_file+"SNVs_{}_{}/SigProfilerAssignment_output".format(marker, label), 
                   input_type="matrix")


# Extract activity values for each signature
SBS_freqs_ctip = pd.read_csv(path_to_file+"SNVs_CtIP_{}/SigProfilerAssignment_output/Assignment_Solution/Activities/Assignment_Solution_Activities.txt".format(label), 
                             sep = "\t", index_col='Samples')
total_freqs_ctip = SBS_freqs_ctip.sum().sort_values(ascending=False)[0:5]
# most frequent, aka present in most samples, although not necessarily with more mutations x signature
(SBS_freqs_ctip != 0).sum().sort_values(ascending=False)[0:5] 

SBS_freqs_grhl = pd.read_csv(path_to_file+"SNVs_GRHL_{}/SigProfilerAssignment_output/Assignment_Solution/Activities/Assignment_Solution_Activities.txt".format(label), 
                             sep = "\t", index_col='Samples')
total_freqs_grhl = SBS_freqs_grhl.sum().sort_values(ascending=False)[0:5]
(SBS_freqs_grhl != 0).sum().sort_values(ascending=False)[0:5] 





### Enhancers low ###
label = "low"
clust_low_ctip = ["CtIP_cluster_6.2_Enh", "CtIP_cluster_6.3_Enh", "CtIP_cluster_6.1_Enh", "CtIP_cluster_6.0"]
clust_low_grhl = ["GRHL_cluster_5.2_Enh", "GRHL_cluster_5.3_Enh", "GRHL_cluster_5.1_Enh", "GRHL_cluster_5.0"]

## High enhancers
SSMs_ctip = pd.read_csv("/Users/ieo6983/Desktop/fragile_enhancer_clinical/results/ICGC/enhancers_SSMs_overlaps/data/Table_enh_SSMs_CtIP.all_overlaps.tsv", sep = "\t") 
SSMs_ctip = SSMs_ctip[SSMs_ctip['cluster'].isin(clust_low_ctip) & SSMs_ctip['width_sbj'] == 1]
SSMs_grhl = pd.read_csv("/Users/ieo6983/Desktop/fragile_enhancer_clinical/results/ICGC/enhancers_SSMs_overlaps/data/Table_enh_SSMs_GRHL.all_overlaps.tsv", sep = "\t") 
SSMs_grhl = SSMs_grhl[SSMs_grhl['cluster'].isin(clust_low_grhl) & SSMs_grhl['width_sbj'] == 1]

# Convert to appropriate format - Info on input file txt format: https://osf.io/xfphr
format_ssms_ctip = {
    'Project' : SSMs_ctip['project_code'],
    'Sample' : SSMs_ctip['icgc_sample_id'],
    'ID' : SSMs_ctip['icgc_donor_id'],
    'Genome' : "GRCh37",
    'mut_type' : "SNP", 
    'chrom' : SSMs_ctip['seqnames_sbj'], 
    'pos_start' : SSMs_ctip['start_sbj'], 
    'pos_end' : SSMs_ctip['end_sbj'], 
    'ref' : SSMs_ctip['reference_genome_allele'], 
    'alt' : SSMs_ctip['mutated_to_allele'], 
    'Type' : "SOMATIC"
}
format_ssms_ctip = pd.DataFrame(format_ssms_ctip)

format_ssms_grhl = {
    'Project' : SSMs_grhl['project_code'],
    'Sample' : SSMs_grhl['icgc_sample_id'],
    'ID' : SSMs_grhl['icgc_donor_id'],
    'Genome' : "GRCh37",
    'mut_type' : "SNP", 
    'chrom' : SSMs_grhl['seqnames_sbj'], 
    'pos_start' : SSMs_grhl['start_sbj'], 
    'pos_end' : SSMs_grhl['end_sbj'], 
    'ref' : SSMs_grhl['reference_genome_allele'], 
    'alt' : SSMs_grhl['mutated_to_allele'], 
    'Type' : "SOMATIC"
}
format_ssms_grhl = pd.DataFrame(format_ssms_grhl)

# Save formatted file - required for matrix generation. The folder with input file will correspond to the output folder
OUT_FOLDER = "/Users/ieo6983/Desktop/fragile_enhancer_clinical/results/ICGC/mutational_signature_analysis/"
#format_ssms_ctip.to_csv(OUT_FOLDER+'SNVs_CtIP_low/input.SNVs_CtIP_low.txt', index=None, sep = '\t')
#format_ssms_grhl.to_csv(OUT_FOLDER+'SNVs_GRHL_low/input.SNVs_GRHL_low.txt', index=None, sep = '\t')

# Save version with only 1 sample name
format_ssms_ctip['Sample'] = "ONE_SAMPLE"
format_ssms_grhl['Sample'] = "ONE_SAMPLE"
#format_ssms_ctip.to_csv(OUT_FOLDER+'SNVs_CtIP_low.ONE_SAMPLE/input.SNVs_ctip_low.ONE_SAMPLE.txt', index=None, sep = '\t')
#format_ssms_grhl.to_csv(OUT_FOLDER+'SNVs_GRHL_low.ONE_SAMPLE/input.SNVs_grhl_low.ONE_SAMPLE.txt', index=None, sep = '\t')


## Generate mutational matrix - For each sample
path_to_file = "/Users/ieo6983/Desktop/fragile_enhancer_clinical/results/ICGC/mutational_signature_analysis/" # Full path of the saved input files in the desired output folder.

matrices = matGen.SigProfilerMatrixGeneratorFunc("SNVs_CtIP_{}".format(label), "GRCh37", path_to_file+'SNVs_CtIP_{}/'.format(label), 
                                                 exome=False, bed_file=None, chrom_based=False, 
                                                 plot=True, tsb_stat=False, seqInfo=True)
matrices2 = matGen.SigProfilerMatrixGeneratorFunc("SNVs_GRHL_{}".format(label), "GRCh37", path_to_file+'SNVs_GRHL_{}/'.format(label), 
                                                 exome=False, bed_file=None, chrom_based=False, 
                                                 plot=True, tsb_stat=False, seqInfo=True)

## Generate mutational matrix - For ONE_SAMPLE only (unified sample contaning all of them)
path_to_file = "/Users/ieo6983/Desktop/fragile_enhancer_clinical/results/ICGC/mutational_signature_analysis/" # Full path of the saved input files in the desired output folder.

matrices = matGen.SigProfilerMatrixGeneratorFunc("SNVs_CtIP_{}".format(label), "GRCh37", path_to_file+'SNVs_CtIP_{}.ONE_SAMPLE/'.format(label), 
                                                 exome=False, bed_file=None, chrom_based=False, 
                                                 plot=True, tsb_stat=False, seqInfo=True)
matrices2 = matGen.SigProfilerMatrixGeneratorFunc("SNVs_GRHL_{}".format(label), "GRCh37", path_to_file+'SNVs_GRHL_{}.ONE_SAMPLE/'.format(label), 
                                                 exome=False, bed_file=None, chrom_based=False, 
                                                 plot=True, tsb_stat=False, seqInfo=True)


# Use SigProfilerAssignment's main function for mutational signatures assignment
# COSMIC v3.3 mutational signatures are used by default as the input reference signatures
MARKERS = ["CtIP", "GRHL"]

file_pattern = "*SBS96.all"
for marker in MARKERS:
    print(marker)

    dir_matrices = path_to_file + "SNVs_{}_{}/output/SBS/".format(marker, label)        
    file_list = glob.glob(dir_matrices+file_pattern)[0]

    Analyze.cosmic_fit(samples=file_list, # path to the input somatic mutations file (if using segmentation file/mutational matrix)
                   output=path_to_file+"SNVs_{}_{}/SigProfilerAssignment_output".format(marker, label), 
                   input_type="matrix")


# Extract activity values for each signature
SBS_freqs_ctip = pd.read_csv(path_to_file+"SNVs_CtIP_{}/SigProfilerAssignment_output/Assignment_Solution/Activities/Assignment_Solution_Activities.txt".format(label), 
                             sep = "\t", index_col='Samples')
total_freqs_ctip = SBS_freqs_ctip.sum().sort_values(ascending=False)[0:5]
print(total_freqs_ctip)
# most frequent, aka present in most samples, although not necessarily with more mutations x signature
(SBS_freqs_ctip != 0).sum().sort_values(ascending=False)[0:5] 

SBS_freqs_grhl = pd.read_csv(path_to_file+"SNVs_GRHL_{}/SigProfilerAssignment_output/Assignment_Solution/Activities/Assignment_Solution_Activities.txt".format(label), 
                             sep = "\t", index_col='Samples')
total_freqs_grhl = SBS_freqs_grhl.sum().sort_values(ascending=False)[0:5]
print(total_freqs_grhl)
(SBS_freqs_grhl != 0).sum().sort_values(ascending=False)[0:5] 



import pandas as pd
import matplotlib.pyplot as plt
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen # matrix generation
from SigProfilerAssignment import Analyzer as Analyze # signature analysis
import sigProfilerPlotting as sigPlt # signatures plots
import glob
import os

# Install reference genome
#from SigProfilerMatrixGenerator import install as genInstall
#genInstall.install('GRCh37', bash=True)

WIN = 500
UNIT = "bp"

IN_FOLDER = "/Users/ieo6983/Desktop/fragile_enhancer_clinical/results/ICGC/enhancers_SSMs_overlaps/data/"
#OUT_FOLDER = "/Users/ieo6983/Desktop/fragile_enhancer_clinical/results/ICGC/mutational_signature_analysis.{}{}_WIN/".format(WIN, UNIT)
OUT_FOLDER = "/Users/ieo6983/Desktop/fragile_enhancer_clinical/results/ICGC/mutational_signature_analysis.{}{}_WIN.BRCA_D/".format(WIN, UNIT)
if not os.path.exists(OUT_FOLDER):
    os.makedirs(OUT_FOLDER)


##


# Read variants 
#SSMs_grhl = pd.read_csv(IN_FOLDER+"Table_enh_SSMs_GRHL.all_overlaps.{}{}_WIN.tsv".format(WIN,UNIT), sep = "\t") 
SSMs_grhl = pd.read_csv(IN_FOLDER+"Table_enh_SSMs_GRHL.all_overlaps.{}{}_WIN.BRCA_D.tsv".format(WIN,UNIT), sep = "\t") 
#SSMs_grhl = pd.read_csv(IN_FOLDER+"Table_enh_SSMs_GRHL.all_overlaps.tsv", sep = "\t") 

# Define 'high' and 'low' clusters
clust_high_grhl = ["GRHL_cluster_1_Enh", "GRHL_cluster_2_Enh", "GRHL_cluster_4_Enh"]
#clust_high_ctip = ["CtIP_cluster_2_Enh", "CtIP_cluster_3_Enh", "CtIP_cluster_5_Enh"]
clust_low_grhl = ["GRHL_cluster_5.2_Enh", "GRHL_cluster_5.3_Enh", "GRHL_cluster_5.1_Enh", "GRHL_cluster_5.0"]
#clust_low_ctip = ["CtIP_cluster_6.2_Enh", "CtIP_cluster_6.3_Enh", "CtIP_cluster_6.1_Enh", "CtIP_cluster_6.0"]

clusters = {
    "high" : clust_high_grhl, 
    "low" : clust_low_grhl
    }


##


# Generate mutational matrices 
for label in clusters.keys():
    print(label)

    SSMs_grhl_sub = SSMs_grhl[SSMs_grhl['cluster'].isin(clusters[label]) & SSMs_grhl['width_sbj'] == 1]

    # Convert to appropriate format 
    format_ssms_grhl = {
        'Project' : SSMs_grhl_sub['project_code'],
        'Sample' : SSMs_grhl_sub['icgc_sample_id'],
        'ID' : SSMs_grhl_sub['icgc_donor_id'],
        'Genome' : "GRCh37",
        'mut_type' : "SNP", 
        'chrom' : SSMs_grhl_sub['seqnames_sbj'], 
        'pos_start' : SSMs_grhl_sub['start_sbj'], 
        'pos_end' : SSMs_grhl_sub['end_sbj'], 
        'ref' : SSMs_grhl_sub['reference_genome_allele'], 
        'alt' : SSMs_grhl_sub['mutated_to_allele'], 
        'Type' : "SOMATIC"
    }
    format_ssms_grhl = pd.DataFrame(format_ssms_grhl)

    # Save formatted file - required for matrix generation. The folder with input file will correspond to the output folder
    dir_results = OUT_FOLDER+'SNVs_GRHL_{}/'.format(label)
    if not os.path.exists(dir_results):
        os.makedirs(dir_results)
    format_ssms_grhl.to_csv(dir_results+'input.SNVs_grhl_{}.txt'.format(label), index=None, sep = '\t')

    #

    ## Generate mutational matrix - For each sample
    path_to_file = OUT_FOLDER # Full path of the saved input files in the desired output folder.
    
    print("Generating mutational matrix for GRHL2 SNVs - Cluster {}".format(label))
    matrices = matGen.SigProfilerMatrixGeneratorFunc("SNVs_GRHL_{}".format(label), "GRCh37", path_to_file+'SNVs_GRHL_{}/'.format(label), 
                                                    exome=False, bed_file=None, chrom_based=False, 
                                                    plot=True, tsb_stat=False, seqInfo=True)
    

##



# Use SigProfilerAssignment's main function for mutational signatures assignment
# COSMIC v3.3 mutational signatures are used by default as the input reference signatures

for label in clusters.keys():
    file_pattern = "*SBS96.all"

    dir_matrices = path_to_file + "SNVs_GRHL_{}/output/SBS/".format(label)        
    file_list = glob.glob(dir_matrices+file_pattern)[0]

    dir_results_signatures = path_to_file+"SNVs_GRHL_{}/SigProfilerAssignment_output".format(label) 
    if not os.path.exists(dir_results_signatures):
        os.makedirs(dir_results_signatures)

    Analyze.cosmic_fit(samples = file_list, # path to the input somatic mutations file (if using segmentation file/mutational matrix)
                   output = dir_results_signatures, 
                   input_type = "matrix")


##


# Extract activity values for each signature
top_n = 5
for label in clusters.keys():
    dir_results_signatures_activities = OUT_FOLDER+"SNVs_GRHL_{}/SigProfilerAssignment_output/Assignment_Solution/Activities/".format(label) 

    SBS_freqs_grhl = pd.read_csv(dir_results_signatures_activities+"Assignment_Solution_Activities.txt", 
                             sep = "\t", index_col='Samples')

    # most frequent, aka present in most samples, although not necessarily with more mutations x signature
    total_freqs_grhl = SBS_freqs_grhl.sum().sort_values(ascending=False)[0:top_n]

    print(f'\n{label} - GRHL2 enhancers\n')

    print("Most active Signatures - Signatures with most mutations")
    print(
        total_freqs_grhl
    )

    print("Most frequent Signatures - Signatures present in most samples")
    print(
        (SBS_freqs_grhl != 0).sum().sort_values(ascending=False)[0:top_n]
    )


##


## De-novo signature analysis
from SigProfilerExtractor import sigpro as sig # for de-novo extraction of signatures
from SigProfilerAssignment import Analyzer as Analyze # to assign SSMs to a set of de-novo signatures

# Extract signatures de-novo from mutational matrix 
for label in clusters.keys():
    file_pattern = "*SBS96.all"

    dir_matrices = OUT_FOLDER + "SNVs_GRHL_{}/output/SBS/".format(label)        
    input_file = glob.glob(dir_matrices+file_pattern)[0]

    out_folder_sig = OUT_FOLDER+"SNVs_GRHL_{}/de_novo/sig_extraction".format(label) 
    if not os.path.exists(out_folder_sig):
        os.makedirs(out_folder_sig)

    sig.sigProfilerExtractor(input_data=input_file, input_type="matrix", output=out_folder_sig, 
                             reference_genome="GRCh37", minimum_signatures=1, maximum_signatures=25)



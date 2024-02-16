import os

### This script creates all the folders of the pipeline and works to link to the modules used by the pipeline.

### Initially, when running the script your working directory is the Script folder
### the folders function navigates to the base_folder and then Pipeline
### and creates all necessary folders and save paths as variables

base_folder = os.getcwd()


if 'Pipeline' not in os.listdir():
    os.mkdir('Pipeline')
pipeline_folder = os.getcwd() + '/Pipeline'

if 'Cured_folder' not in os.listdir(pipeline_folder):
    os.mkdir(os.getcwd() + '/Pipeline/Cured_folder')
cured_folder = os.getcwd() + '/Pipeline/Cured_folder/'

if 'Peak_folder_0' not in os.listdir(pipeline_folder):
    os.mkdir(os.getcwd() + '/Pipeline/Peak_folder_0')
peak_folder_0 = os.getcwd() + '/Pipeline/Peak_folder_0/'

if 'Peak_folder_0_zscores' not in os.listdir(pipeline_folder):
    os.mkdir(os.getcwd() + '/Pipeline/Peak_folder_0_zscores')
peak_folder_0_zscores = os.getcwd() + '/Pipeline/Peak_folder_0_zscores/'

if 'Peak_folder_1' not in os.listdir(pipeline_folder):
    os.mkdir(os.getcwd() + '/Pipeline/Peak_folder_1')
peak_folder_1 = os.getcwd() + '/Pipeline/Peak_folder_1/'

if 'Peak_folder_2' not in os.listdir(pipeline_folder):
    os.mkdir(os.getcwd() + '/Pipeline/Peak_folder_2')
peak_folder_2 = os.getcwd() + '/Pipeline/Peak_folder_2/'

if 'Peak_folder_3' not in os.listdir(pipeline_folder):
    os.mkdir(os.getcwd() + '/Pipeline/Peak_folder_3')
peak_folder_3 = os.getcwd() + '/Pipeline/Peak_folder_3/'

if 'GFF3_folder_0' not in os.listdir(pipeline_folder):
    os.mkdir(os.getcwd() + '/Pipeline/GFF3_folder_0')
GFF3_folder_0 = os.getcwd() + '/Pipeline/GFF3_folder_0/'

if 'GFF3_folder_1' not in os.listdir(pipeline_folder):
    os.mkdir(os.getcwd() + '/Pipeline/GFF3_folder_1')
GFF3_folder_1 = os.getcwd() + '/Pipeline/GFF3_folder_1/'

if 'GFF3_folder_2' not in os.listdir(pipeline_folder):
    os.mkdir(os.getcwd() + '/Pipeline/GFF3_folder_2')
GFF3_folder_2 = os.getcwd() + '/Pipeline/GFF3_folder_2/'

if 'missing_peaks' not in os.listdir(pipeline_folder):
    os.mkdir(os.getcwd() + '/Pipeline/missing_peaks')
missing_peaks = os.getcwd() + '/Pipeline/missing_peaks/'

if 'GFF3_folder_3' not in os.listdir(pipeline_folder):
    os.mkdir(os.getcwd() + '/Pipeline/GFF3_folder_3')
GFF3_folder_3 = os.getcwd() + '/Pipeline/GFF3_folder_3/'

if 'Annotation_files' not in os.listdir(pipeline_folder):
    os.mkdir(os.getcwd() + '/Pipeline/Annotation_files')
GFF3_folder_4 = os.getcwd() + '/Pipeline/Annotation_files/'

if 'Operon_folder' not in os.listdir(pipeline_folder):
    os.mkdir(os.getcwd() + '/Pipeline/Operon_folder')
operon_folder = os.getcwd() + '/Pipeline/Operon_folder/'

if 'Analysis' not in os.listdir(pipeline_folder):
    os.mkdir(os.getcwd() + '/Pipeline/Analysis')
analysis_folder = os.getcwd() + '/Pipeline/Analysis/'

if 'Distance_to_old_annotation' not in os.listdir(analysis_folder):
    os.mkdir(analysis_folder + '/Distance_to_old_annotation')
distance_folder = analysis_folder + '/Distance_to_old_annotation'

if 'differentially_expressed_genes' not in os.listdir(analysis_folder):
    os.mkdir(analysis_folder + '/differentially_expressed_genes')
differentially_expressed_genes = analysis_folder + '/differentially_expressed_genes'

# if 'differentially_expressed_genes_All' not in os.listdir(differentially_expressed_genes):
#     os.mkdir(differentially_expressed_genes + '/differentially_expressed_genes_All')
# differentially_expressed_genes_All = differentially_expressed_genes + '/differentially_expressed_genes_All'

# if 'differentially_expressed_genes_HighConfidence' not in os.listdir(differentially_expressed_genes):
#     os.mkdir(differentially_expressed_genes + '/differentially_expressed_genes_HighConfidence')
# differentially_expressed_genes_HighConfidence = differentially_expressed_genes + '/differentially_expressed_genes_HighConfidence'

# if 'differentially_expressed_genes_HighDensity' not in os.listdir(differentially_expressed_genes):
#     os.mkdir(differentially_expressed_genes + '/differentially_expressed_genes_HighDensity')
# differentially_expressed_genes_HighDensity = differentially_expressed_genes + '/differentially_expressed_genes_HighDensity'

if 'differentially_expressed_genes_Overlap' not in os.listdir(differentially_expressed_genes):
    os.mkdir(differentially_expressed_genes + '/differentially_expressed_genes_Overlap')
differentially_expressed_genes_Overlap = differentially_expressed_genes + '/differentially_expressed_genes_Overlap'

# if 'differentially_expressed_genes_Overlap_lessStringent' not in os.listdir(differentially_expressed_genes):
#     os.mkdir(differentially_expressed_genes + '/differentially_expressed_genes_Overlap_lessStringent')
# differentially_expressed_genes_Overlap_lessStringent = differentially_expressed_genes + '/differentially_expressed_genes_Overlap_lessStringent'



### Part of the pipeline is done with multiprocessing for faster analysis. This is done with seperate scripts incorporated into this pipeline.
curing_script = os.getcwd() + '/Scripts/Curing.py'
peak_calling_script = os.getcwd() + '/Scripts/Peak_calling.py'
operon_annotation_script = os.getcwd() + '/Scripts/Make_RendSeq_operon_annotation.py'
assign_and_GFF3_script = os.getcwd() + '/Scripts/Draw_GFF3.py'
delete_redundant_script = os.getcwd() + '/Scripts/delete_redundant.py'
fdr_script = os.getcwd() + '/Scripts/fdr.py'
GFF_to_GTF_script = os.getcwd() + '/Scripts/GFF_to_GTF.py'
GTF_to_GFF_script = os.getcwd() + '/Scripts/GTF_to_GFF.py'
supplyCapTerm_script = os.getcwd() + '/Scripts/supplyCapTerm.py'

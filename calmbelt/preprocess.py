import os
import pickle
import nbformat
from nbconvert.preprocessors import ExecutePreprocessor

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("genome_dir", help="Directory containing input genomes (500-10000 genomes) in fasta format")
parser.add_argument("reference_fname", help="A reference file")
parser.add_argument("reference_gene_loci", help="A tsv file containing a list of genes ane their loci")
parser.add_argument("predefined_clade", help="A tsv file containing a list of SNPs for each GISAID/custom clade")
parser.add_argument("predefined_label", help="A tsv file containing WHO label for each PANGO lineage")
parser.add_argument("preprocess_dir", help="Output directory for preprocessed data")
parser.add_argument("world_metadata_fname", help="A metadata for samples around the world (only be used for "
                                                 "illustrating overall trends)")
parser.add_argument("country_name", help="A country name for country-specific stats, e.g. Singapore")
parser.add_argument("subsampling_n_samples", nargs='?', default=-1, help=">= 500; Default = -1, i.e. all samples")

args = parser.parse_args()
print(args)

# Save arguments
pickle.dump(args, open('arguments.pickle', 'wb'))

# Create preprocess_dir if not exist
if not os.path.exists(args.preprocess_dir):
    os.mkdir(args.preprocess_dir)


def run_ipynb(notebook_fname, script_path='/'.join(os.path.realpath(__file__).split('/')[:-1])):
    notebook_path = os.path.join(script_path, notebook_fname)
    notebook_output_path = notebook_path.replace('.ipynb', '.nbconvert.ipynb')

    with open(notebook_path) as f:
        nb = nbformat.read(f, as_version=4)
    ep = ExecutePreprocessor(timeout=600, kernel_name='python3')
    ep.preprocess(nb, {'metadata': {'path': '.'}})
    with open(notebook_output_path, 'w', encoding='utf-8') as f:
        nbformat.write(nb, f)


print('Process metadata of patients around the world ...')
# Rolling average for world_metadata_fname
run_ipynb('Step0_plot_metadata.ipynb')
# Summarize number of lineage per month
run_ipynb('Step12_meta_data.ipynb')

print('Blast all genomes in genome_dir to reference, and then find SNPs and indel ...')
run_ipynb('Step1_split_blast.ipynb')
run_ipynb('Step2_align_delete_insertion.ipynb')
run_ipynb('Step3_create_genome_array.ipynb')
run_ipynb('Step3.1_insertion_position.ipynb')
run_ipynb('Step3.2_create_insertion_array.ipynb')

print('Identify mutation at protein level ...')
run_ipynb('Step9_change_protein.ipynb')

print('Identify mutations for each month ...')
run_ipynb('Step11_alarm_more_patient.ipynb')

print('Select positions based on entropy, perform clustering, and calculate pairwise distance ...')
run_ipynb('Step4_entropy.ipynb')
run_ipynb('Step6_prepare_X_cluster.ipynb')
run_ipynb('Step7_pairwise_distance.ipynb')
run_ipynb('Step8_clustering.ipynb')

print('Create dendrogram ...')
run_ipynb('Step10_dendro.ipynb')
print(f'Done. All preprocessing files were saved in {args.preprocess_dir}')

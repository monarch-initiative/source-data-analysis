
from utils import os, sys, read_g2d_files, read_mapping_files, prep_g2d_gencc, prep_g2d_omim, prep_g2d_hpoa, prep_g2d_orpha, venn4_plot_g2d

#parent_directory = os.path.abspath(os.path.join(os.getcwd(), ".."))
#input_data_folder_path = os.path.join(parent_directory, 'data/input/g2d')
#output_data_folder_path = os.path.join(parent_directory, 'data/output/g2d')

input_data_folder_path = './data/input/g2d'
output_data_folder_path = './data/output/g2d'

# Read files
gencc_df, omim_df, hpoa_df, orpha_tree = read_g2d_files(input_data_folder_path)
omim2mondo_df, orpha2mondo_df = read_mapping_files(input_data_folder_path)

# Prep files
gencc = prep_g2d_gencc(gencc_df)
omim = prep_g2d_omim(omim_df, omim2mondo_df, orpha2mondo_df)
hpoa = prep_g2d_hpoa(hpoa_df, omim2mondo_df, orpha2mondo_df)
orpha = prep_g2d_orpha(orpha_tree, omim2mondo_df, orpha2mondo_df)

# Plot output
old_stdout = sys.stdout # backup current stdout
sys.stdout = open(os.devnull, "w")
venn4_plot_g2d(gencc, omim, hpoa, orpha, output_data_folder_path)
sys.stdout = old_stdout # reset old stdout

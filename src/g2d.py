import os
from utils import os, read_g2d_files, read_mapping_files, prep_g2d_gencc, prep_g2d_omim

parent_directory = os.path.abspath(os.path.join(os.getcwd(), ".."))
data_folder_path = os.path.join(parent_directory, 'data')
output_folder_path = os.path.join(parent_directory, 'output')

# Read files
gencc_df, omim_df, hpoa_df, orpha_tree = read_g2d_files(data_folder_path)
omim2mondo_df, orpha2mondo_df = read_mapping_files(data_folder_path)

# Prep files
gencc = prep_g2d_gencc(gencc_df)
omim = prep_g2d_omim(omim_df, omim2mondo_df, orpha2mondo_df)
#hpoa = prep_g2d_hpoa(hpoa_df)
#orpha = prep_g2d_orpha(orpha_tree)

print(omim.head())

# Plot output


# Write output


# Prep summary





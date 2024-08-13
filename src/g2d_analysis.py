import os
import sys
from itertools import combinations
import xml.etree.ElementTree as ET

import pandas as pd
from venn import venn
import matplotlib.pyplot as plt
from upsetplot import UpSet


def read_g2d_files(folder_path):
    files = os.listdir(folder_path)
    gencc_df, omim_df, hpoa_df, orpha_tree, medgen_df = None, None, None, None, None
    missing_files = []

    # Define the expected files and their corresponding variables
    expected_files = {
        'gencc-submissions.tsv': 'gencc_df',
        'morbidmap.txt': 'omim_df',
        'genes_to_disease.txt': 'hpoa_df',
        'en_product6.xml': 'orpha_tree',
        'mim2gene_medgen.txt': 'medgen_df'
    }

    # Check for missing files
    for expected_file, var_name in expected_files.items():
        if expected_file not in files:
            # Remove suffixes starting with an underscore
            readable_var_name = var_name.split('_')[0]
            missing_files.append(f'{readable_var_name}: {expected_file}')

    # Print missing files
    if missing_files:
        print("Missing files:")
        for file in missing_files:
            print(f"- {file}")

    # Read the files if they are present
    for file in files:
        file_path = os.path.join(folder_path, file)
        if file == 'gencc-submissions.tsv':
            gencc_df = pd.read_csv(file_path, sep='\t')
            gencc_df = gencc_df[gencc_df['classification_title'].isin(['Moderate', 'Strong', 'Definitive'])]
        elif file == 'morbidmap.txt':
            omim_df = pd.read_csv(file_path, sep='\t', skiprows=3)
        elif file == 'genes_to_disease.txt':
            hpoa_df = pd.read_csv(file_path, sep='\t')
        elif file == 'en_product6.xml':
            orpha_tree = ET.parse(file_path)
        elif file == 'mim2gene_medgen.txt':
            medgen_df = pd.read_csv(file_path, sep='\t')
            medgen_df = medgen_df[medgen_df['type'].isin(['phenotype'])]
    return gencc_df, omim_df, hpoa_df, orpha_tree, medgen_df, missing_files

def read_mapping_files(folder_path):
    files = os.listdir(folder_path)
    disease_mappings_df, gene_mappings_df = None, None

    # Define the expected files and their corresponding variables
    expected_files = {
        'mondo.sssom.tsv': 'disease_mappings_df',
        'gene_mappings.sssom.tsv': 'gene_mappings_df'
    }

    # Check for missing files
    for expected_file, var_name in expected_files.items():
        if expected_file not in files:
            # Remove suffixes starting with an underscore
            readable_var_name = var_name.split('_')[0]
            print(f"Missing mapping file (required): {readable_var_name}: {expected_file}")
            sys.exit(1)  # Exit the script with a non-zero status code

    # Read the files if they are present
    for file in files:
        file_path = os.path.join(folder_path, file)
        if file == 'mondo.sssom.tsv':
            disease_mappings_df = pd.read_csv(file_path, sep='\t', skiprows=50)
            disease_mappings_df['object_id'] = disease_mappings_df['object_id'].str.replace('Orphanet:', 'ORPHA:', regex=False)
        elif file == 'gene_mappings.sssom.tsv':
            gene_mappings_df = pd.read_csv(file_path, sep='\t', skiprows=25)
            gene_mappings_df = gene_mappings_df[gene_mappings_df['subject_id'].str.startswith('HGNC:')]

    return disease_mappings_df, gene_mappings_df

def prep_g2d_gencc(gencc_df):
    print("Starting GENCC data preparation...")

    # Ensure the DataFrame is not a view
    gencc_df = gencc_df.copy()

    # Remove rows with NaN in 'disease_curie' and 'gene_curie'
    gencc_df = gencc_df[gencc_df['disease_curie'].notna()]
    gencc_df = gencc_df[gencc_df['gene_curie'].notna()]

    before_mapping = {
        'source': "GENCC",
        'disease_count': gencc_df['disease_curie'].notna().sum(),
        'disease_count_unique': len(gencc_df['disease_curie'].unique()),
        'gene_count': gencc_df['gene_curie'].notna().sum(),
        'gene_count_unique': len(gencc_df['gene_curie'].unique()),
        'unique_g2d_association_count': len(gencc_df.groupby(['disease_curie', 'gene_curie']).size())
    }

    after_mapping = before_mapping
    print("GENCC data preparation completed")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    return gencc_df, before_mapping, after_mapping

def prep_g2d_hpoa(hpoa_df, disease_mappings_df, gene_mappings_df):
    print("Starting HPO data preparation...")

    # Ensure the DataFrame is not a view
    hpoa_df = hpoa_df.copy()

    # Remove rows with NaN in 'disease_id' and 'ncbi_gene_id'
    hpoa_df = hpoa_df[hpoa_df['disease_id'].notna()]
    hpoa_df = hpoa_df[hpoa_df['ncbi_gene_id'].notna()]

    before_mapping = {
        'source': "HPO",
        'disease_count': hpoa_df['disease_id'].notna().sum(),
        'disease_count_unique': len(hpoa_df['disease_id'].unique()),
        'gene_count': hpoa_df['ncbi_gene_id'].notna().sum(),
        'gene_count_unique': len(hpoa_df['ncbi_gene_id'].unique()),
        'unique_g2d_association_count': len(hpoa_df.groupby(['disease_id', 'ncbi_gene_id']).size())
    }

    # Disease mapping
    print("\tHPO disease mapping in progress")
    hpoa_df["disease_mondo"] = disease_mondo_mapping(hpoa_df["disease_id"], disease_mappings_df)
    missing_diseases = hpoa_df[hpoa_df["disease_mondo"].isna()]["disease_id"].unique().tolist()
    print("\t\tDisease # with no mappings:", len(missing_diseases))
    print("\tHPO disease mapping complete")

    # Gene mapping
    print("\tHPO gene mapping in progress")
    hpoa_df["gene_hgnc"] = gene_hgnc_mapping(hpoa_df["ncbi_gene_id"].tolist(), gene_mappings_df)
    missing_genes = hpoa_df[hpoa_df["gene_hgnc"].isna()]["ncbi_gene_id"].unique().tolist()
    print("\t\tGene # with no mappings:", len(missing_genes))
    print("\tHPO gene mapping complete")

    after_mapping = {
        'source': "HPO",
        'disease_count': hpoa_df['disease_mondo'].notna().sum(),
        'disease_count_unique': len(hpoa_df['disease_mondo'].unique()),
        'gene_count': hpoa_df['gene_hgnc'].notna().sum(),
        'gene_count_unique': len(hpoa_df['gene_hgnc'].unique()),
        'unique_g2d_association_count': len(hpoa_df.groupby(['disease_mondo', 'gene_hgnc']).size())
    }
    print("HPO data preparation completed")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    hpoa_df = hpoa_df.dropna(subset=['disease_mondo', 'gene_hgnc'])
    return hpoa_df, before_mapping, after_mapping, missing_diseases, missing_genes

def prep_g2d_omim(omim_df, disease_mappings_df, gene_mappings_df):
    print("Starting OMIM data preparation...")

    # Ensure the DataFrame is not a view
    omim_df = omim_df.copy()

    # Clean and prepare '# Phenotype' column
    omim_df["# Phenotype"] = omim_df["# Phenotype"].str.split(',').str[-1]
    omim_df["# Phenotype"] = omim_df["# Phenotype"].str[1:7]
    omim_df['# Phenotype'] = pd.to_numeric(omim_df['# Phenotype'], errors='coerce')

    # Remove rows with NaN in '# Phenotype'
    omim_df = omim_df[omim_df['# Phenotype'].notna()]

    # Convert '# Phenotype' to integer
    omim_df['# Phenotype'] = omim_df['# Phenotype'].astype(int)

    # Remove rows with NaN in 'MIM Number' and convert to integer
    omim_df = omim_df[omim_df['MIM Number'].notna()]
    omim_df['MIM Number'] = omim_df['MIM Number'].astype(int)

    before_mapping = {
        'source': "OMIM",
        'disease_count': omim_df['# Phenotype'].notna().sum(),
        'disease_count_unique': len(omim_df['# Phenotype'].unique()),
        'gene_count': omim_df['MIM Number'].notna().sum(),
        'gene_count_unique': len(omim_df['MIM Number'].unique()),
        'unique_g2d_association_count': len(omim_df.groupby(['# Phenotype', 'MIM Number']).size())
    }

    # Prefix OMIM
    omim_df["# Phenotype"] = 'OMIM:' + omim_df["# Phenotype"].astype(str)
    omim_df["MIM Number"] = 'OMIM:' + omim_df["MIM Number"].astype(str)

    # Disease mapping
    print("\tOMIM disease mapping in progress")
    omim_df["disease_mondo"] = disease_mondo_mapping(omim_df["# Phenotype"], disease_mappings_df)
    missing_diseases = omim_df[omim_df["disease_mondo"].isna()]["# Phenotype"].unique().tolist()
    print(f"\t\tDisease # with no mappings: {len(missing_diseases)}")
    print("\tOMIM disease mapping complete")

    # Gene mapping
    print("\tOMIM gene mapping in progress")
    omim_df["gene_hgnc"] = gene_hgnc_mapping(omim_df["MIM Number"].tolist(), gene_mappings_df)
    missing_genes = omim_df[omim_df["gene_hgnc"].isna()]["MIM Number"].unique().tolist()
    print(f"\t\tGene # with no mappings: {len(missing_genes)}")
    print("\tOMIM gene mapping complete")

    after_mapping = {
        'source': "OMIM",
        'disease_count': omim_df['disease_mondo'].notna().sum(),
        'disease_count_unique': len(omim_df['disease_mondo'].unique()),
        'gene_count': omim_df['gene_hgnc'].notna().sum(),
        'gene_count_unique': len(omim_df['gene_hgnc'].unique()),
        'unique_g2d_association_count': len(omim_df.groupby(['disease_mondo', 'gene_hgnc']).size())
    }

    print("OMIM data preparation completed")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    omim_df = omim_df.dropna(subset=['disease_mondo', 'gene_hgnc'])
    return omim_df, before_mapping, after_mapping, missing_diseases, missing_genes

def prep_g2d_medgen(medgen_df, disease_mappings_df, gene_mappings_df):
    print("Starting MEDGEN data preparation...")

    # Ensure the DataFrame is not a view
    medgen_df = medgen_df.copy()

    medgen_df = medgen_df[medgen_df['#MIM number'].notna()]
    medgen_df = medgen_df[medgen_df['GeneID'].notna()]
    medgen_df = medgen_df[~medgen_df['GeneID'].str.contains('-')]

    before_mapping = {
        'source': "MEDGEN",
        'disease_count': medgen_df['#MIM number'].notna().sum(),
        'disease_count_unique': len(medgen_df['#MIM number'].unique()),
        'gene_count': medgen_df['GeneID'].notna().sum(),
        'gene_count_unique': len(medgen_df['GeneID'].unique()),
        'unique_g2d_association_count': len(medgen_df.groupby(['#MIM number', 'GeneID']).size())
    }

    medgen_df["#MIM number"] = 'OMIM:' + medgen_df["#MIM number"].astype(str)
    medgen_df["GeneID"] = 'NCBIGene:' + medgen_df["GeneID"].astype(str)

    # Disease mapping
    print("\tMEDGEN disease mapping in progress")
    medgen_df["disease_mondo"] = disease_mondo_mapping(medgen_df["#MIM number"], disease_mappings_df)
    missing_diseases = medgen_df[medgen_df["disease_mondo"].isna()]["#MIM number"].unique().tolist()
    print(f"\t\tDisease # with no mappings: {len(missing_diseases)}")
    print("\tMEDGEN disease mapping complete")

    # Gene mapping
    print("\tMEDGEN gene mapping in progress")
    medgen_df["gene_hgnc"] = gene_hgnc_mapping(medgen_df["GeneID"].tolist(), gene_mappings_df)
    missing_genes = medgen_df[medgen_df["gene_hgnc"].isna()]["GeneID"].unique().tolist()
    print(f"\t\tGene # with no mappings: {len(missing_genes)}")
    print("\tMEDGEN gene mapping complete")

    after_mapping = {
        'source': "MEDGEN",
        'disease_count': medgen_df['disease_mondo'].notna().sum(),
        'disease_count_unique': len(medgen_df['disease_mondo'].unique()),
        'gene_count': medgen_df['gene_hgnc'].notna().sum(),
        'gene_count_unique': len(medgen_df['gene_hgnc'].unique()),
        'unique_g2d_association_count': len(medgen_df.groupby(['disease_mondo', 'gene_hgnc']).size())
    }

    print("MEDGEN data preparation completed")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    medgen_df = medgen_df.dropna(subset=['disease_mondo', 'gene_hgnc'])
    return medgen_df, before_mapping, after_mapping, missing_diseases, missing_genes


def prep_g2d_orpha(orpha_tree, disease_mappings_df):
    print("Starting ORPHA data preparation...")

    df_cols = ["Disorder_ID", "Disorder_Name", "Orpha_Code", "gene_hgnc"]
    rows = []

    for node in orpha_tree.iter('Disorder'):
        disorder_id = node.attrib.get('id')
        disorder_name = node.find("Name").text
        orphacode = node.find("OrphaCode").text

        # Iterate through the 'DisorderGeneAssociationList' elements under the current 'Disorder' element
        for nodea in node.iter('DisorderGeneAssociation'):
            gene_hgnc_element = nodea.find(".//Gene/ExternalReferenceList/ExternalReference[Source='HGNC']/Reference")

            # Check if gene_hgnc_element is not None before accessing its text attribute
            if gene_hgnc_element is not None:
                gene_hgnc = gene_hgnc_element.text
                rows.append({
                    "Disorder_ID": disorder_id,
                    "Disorder_Name": disorder_name,
                    "Orpha_Code": orphacode,
                    "gene_hgnc": gene_hgnc
                })

    orpha_df = pd.DataFrame(rows, columns=df_cols)
    orpha_df = orpha_df[orpha_df['Disorder_ID'].notna()]
    orpha_df["Orpha_Code"] = 'ORPHA:' + orpha_df["Orpha_Code"].astype(str)
    orpha_df = orpha_df[orpha_df['gene_hgnc'].notna()]
    orpha_df["gene_hgnc"] = 'HGNC:' + orpha_df["gene_hgnc"].astype(str)

    before_mapping = {
        'disease_count': orpha_df['Orpha_Code'].notna().sum(),
        'disease_count_unique': len(orpha_df['Orpha_Code'].unique()),
        'gene_count': orpha_df['gene_hgnc'].notna().sum(),
        'gene_count_unique': len(orpha_df['gene_hgnc'].unique()),
        'unique_g2d_association_count': len(orpha_df.groupby(['Orpha_Code', 'gene_hgnc']).size())
    }

    # Disease mapping
    print("\tORPHA disease mapping in progress")
    orpha_df["disease_mondo"] = disease_mondo_mapping(orpha_df["Orpha_Code"], disease_mappings_df)
    missing_diseases = orpha_df[orpha_df["disease_mondo"].isna()]["Orpha_Code"].unique().tolist()
    print("\t\tDisease # with no mappings:", len(missing_diseases))
    print("\tORPHA disease mapping complete")

    after_mapping = {
        'source': "HPO",
        'disease_count': orpha_df['disease_mondo'].notna().sum(),
        'disease_count_unique': len(orpha_df['disease_mondo'].unique()),
        'gene_count': orpha_df['gene_hgnc'].notna().sum(),
        'gene_count_unique': len(orpha_df['gene_hgnc'].unique()),
        'unique_g2d_association_count': len(orpha_df.groupby(['disease_mondo', 'gene_hgnc']).size())
    }

    print("ORPHA data preparation completed")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    missing_genes = []
    orpha_df = orpha_df.dropna(subset=['disease_mondo', 'gene_hgnc'])
    return orpha_df, before_mapping, after_mapping, missing_diseases, missing_genes


def map_ids(id_list, mappings_df, id_col='object_id', mapping_col='subject_id'):
    # Create DataFrame from id_list
    id_df = pd.DataFrame({id_col: id_list})

    # Perform the merge
    merged_df = id_df.merge(mappings_df, on=id_col, how='left')

    # Handle missing values by filling NaN with None
    mapped_list = merged_df[mapping_col].apply(lambda x: x if pd.notna(x) else None).tolist()

    return mapped_list

def disease_mondo_mapping(disease_list, disease_mappings_df):
    return map_ids(disease_list, disease_mappings_df)

def gene_hgnc_mapping(gene_list, gene_mappings_df):
    return map_ids(gene_list, gene_mappings_df)

def analysis_summary(before, after):
    summary = {
        'disease_count': before['disease_count'],
        'disease_count_unique': before['disease_count_unique'],
        'gene_count': before['gene_count'],
        'gene_count_unique': before['gene_count_unique'],
        #'total_all_g2d_association_count': before['all_g2d_association_count'],
        #'total_unique_g2d_association_count': before['unique_g2d_association_count'],
        'disease_count (MONDO mapped)': after['disease_count'],
        'disease_count_unique (MONDO mapped)': after['disease_count_unique'],
        'gene_count (HGNC mapped)': after['gene_count'],
        'gene_count_unique (HGNC mapped) ': after['gene_count_unique'],
        #'mapped_all_g2d_association_count': after['all_g2d_association_count'],
        #'mapped_unique_g2d_association_count': after['unique_g2d_association_count']
        'g2d_edges_count (final set)': after['unique_g2d_association_count']
    }
    return summary

def save_analysis_summaries(summaries, output_path):
    with open(output_path, 'w') as file:
        for source, summary in summaries.items():
            file.write(f"Summary for {source.upper()}:\n")
            file.write(f"=============================\n")
            for key, value in summary.items():
                file.write(f"{key}: {value}\n")
            file.write("\n")

def save_missing_mappings(missing_mappings, output_path):
    with pd.ExcelWriter(output_path) as writer:
        for source, mappings in missing_mappings.items():
            # Ensure at least one row is written to the sheet
            if mappings['missing_diseases'] or mappings['missing_genes']:
                max_len = max(len(mappings['missing_diseases']), len(mappings['missing_genes']))
                missing_diseases = mappings['missing_diseases'] + [None] * (max_len - len(mappings['missing_diseases']))
                missing_genes = mappings['missing_genes'] + [None] * (max_len - len(mappings['missing_genes']))
                df = pd.DataFrame({'missing_disease': missing_diseases, 'missing_gene': missing_genes})
                df.to_excel(writer, sheet_name=source.upper(), index=False)

def save_g2d_venn(sets, g2d_venn_file_path):
    plt.figure(figsize=(10, 10))
    venn(sets, fontsize=9, legend_loc="upper right")
    plt.title('G2D VENN')
    plt.savefig(g2d_venn_file_path)
    #plt.show()



def g2d_edges_differences(sets, file_path):
    with pd.ExcelWriter(file_path) as writer:
        # Calculate and save differences
        for (name1, set1), (name2, set2) in combinations(sets.items(), 2):
            diff1 = set1 - set2
            diff2 = set2 - set1

            # Save differences where set1 - set2
            diff1_name = f'{name1[:15]}_minus_{name2[:15]}'  # Truncate set names to 15 characters
            df1 = pd.DataFrame(list(diff1), columns=['Disease', 'Gene'])
            df1.to_excel(writer, sheet_name=diff1_name[:31], index=False)  # Ensure sheet name is within 31 characters

            # Save differences where set2 - set1
            diff2_name = f'{name2[:15]}_minus_{name1[:15]}'  # Truncate set names to 15 characters
            df2 = pd.DataFrame(list(diff2), columns=['Disease', 'Gene'])
            df2.to_excel(writer, sheet_name=diff2_name[:31], index=False)  # Ensure sheet name is within 31 characters

def plot_upset(set_dict, output_file_path):
    # Union of all sets
    all_elems = set().union(*set_dict.values())
    #print("Total unique elements:", len(all_elems))

    # Create DataFrame indicating membership in each set
    df = pd.DataFrame(
        [[e in set_dict['HPO'], e in set_dict['MEDGEN'], e in set_dict['ORPHA'],
          e in set_dict['OMIM'], e in set_dict['GENCC']]
         for e in all_elems],
        columns=['HPO', 'MEDGEN', 'ORPHA', 'OMIM', 'GENCC']
    )
    #print(df.head())

    # Count occurrences of each combination of set memberships
    df_up = df.groupby(['HPO', 'MEDGEN', 'ORPHA', 'OMIM', 'GENCC']).size()
    #print("Grouped Data:")
    #print(df_up.head())

    # Check if df_up is empty
    if df_up.empty:
        print("The grouped DataFrame is empty. Check the input sets for valid data.")
        return

    # Create and save the UpSet plot
    plt.figure(figsize=(10, 6))
    upset = UpSet(df_up, orientation='horizontal', show_counts=True, sort_by='input')
    upset.plot()
    plt.title('G2D UpSet Plot')
    plt.savefig(output_file_path)
    #plt.show()

def consolidate_g2d_edges(sets, file_path):
    consolidated_data = []

    for name, data in sets.items():
        for disease_id, gene_id in data:
            consolidated_data.append({
                'disease_id': disease_id,
                'gene_id': gene_id,
                'sources': name
            })

    # Create DataFrame
    consolidated_df = pd.DataFrame(consolidated_data)

    # Group by disease_id and gene_id, then concatenate sources
    consolidated_df = (consolidated_df.groupby(['disease_id', 'gene_id'])
                       .agg({'sources': lambda x: '|'.join(sorted(set(x)))})
                       .reset_index())

    # Save to Excel
    consolidated_df.to_excel(file_path, index=False)

def g2d_edges_by_sources(sets, file_path):
    with pd.ExcelWriter(file_path) as writer:
        for name, data in sets.items():
            df = pd.DataFrame(list(data), columns=['Disease', 'Gene'])
            df.to_excel(writer, sheet_name=name[:31], index=False)


def consolidate_unique_g2d_edges(sets, file_path):
    # Dictionary to store unique edges for each source
    unique_sets = {}

    # Collect all edges across all sources
    all_edges = set()
    for source, edges in sets.items():
        all_edges.update(edges)

    # Identify unique edges for each source
    for source, edges in sets.items():
        unique_edges = edges - {edge for other_source, other_edges in sets.items() if other_source != source for edge in
                                other_edges}
        unique_sets[source] = unique_edges

    # Save unique edges to Excel file, one sheet per source
    with pd.ExcelWriter(file_path) as writer:
        for source, edges in unique_sets.items():
            unique_df = pd.DataFrame(list(edges), columns=['disease_id', 'gene_id'])
            unique_df.to_excel(writer, sheet_name=source[:31], index=False)


def main():
    global gencc_before, gencc_after, HPO_tuples, MEDGEN_tuples, OMIM_tuples, ORPHA_tuples, GENCC_tuples
    analysis_summaries = {}
    missing_mappings = {}
    sets = {}

    try:
        #parent_directory = os.path.abspath(os.path.join(os.getcwd(), ".."))
        parent_directory = os.path.abspath(os.path.join(os.getcwd(), "..", "source-data-analysis"))

        data_folder_path = os.path.join(parent_directory, 'data/input/g2d/')
        output_folder_path = os.path.join(parent_directory, 'data/output/g2d/')

        # Create the output directory if it does not exist
        if not os.path.exists(output_folder_path):
            os.makedirs(output_folder_path, exist_ok=True)

        folder_path = data_folder_path  # Replace with your folder path
        print(folder_path)
        try:
            gencc_df, omim_df, hpoa_df, orpha_tree, medgen_df, missing_files = read_g2d_files(folder_path)
            disease_mappings_df, gene_mappings_df = read_mapping_files(folder_path)
        except Exception as e:
            print(f"Error reading G2D files or mapping files: {e}")
            return

        if hpoa_df is not None and not hpoa_df.empty:
            try:
                hpoa_df, before_mapping, after_mapping, missing_diseases, missing_genes = prep_g2d_hpoa(hpoa_df,
                                                                                                        disease_mappings_df,
                                                                                                        gene_mappings_df)
                analysis_summaries['hpo'] = analysis_summary(before_mapping, after_mapping)
                missing_mappings['hpo'] = {'missing_diseases': missing_diseases, 'missing_genes': missing_genes}
                HPO_tuples = set(hpoa_df[['disease_mondo', 'gene_hgnc']].apply(tuple, axis=1))
                sets['HPO'] = HPO_tuples
            except Exception as e:
                print(f"Error processing HPOA data: {e}")

        if medgen_df is not None and not medgen_df.empty:
            try:
                medgen_df, before_mapping, after_mapping, missing_diseases, missing_genes = prep_g2d_medgen(medgen_df,
                                                                                                            disease_mappings_df,
                                                                                                            gene_mappings_df)
                analysis_summaries['medgen'] = analysis_summary(before_mapping, after_mapping)
                missing_mappings['medgen'] = {'missing_diseases': missing_diseases, 'missing_genes': missing_genes}
                MEDGEN_tuples = set(medgen_df[['disease_mondo', 'gene_hgnc']].apply(tuple, axis=1))
                sets['MEDGEN'] = MEDGEN_tuples
            except Exception as e:
                print(f"Error processing MEDGEN data: {e}")

        if orpha_tree is not None:
            try:
                orpha_df, before_mapping, after_mapping, missing_diseases, missing_genes = prep_g2d_orpha(orpha_tree,
                                                                                                          disease_mappings_df)
                analysis_summaries['orpha'] = analysis_summary(before_mapping, after_mapping)
                missing_mappings['orpha'] = {'missing_diseases': missing_diseases, 'missing_genes': missing_genes}
                ORPHA_tuples = set(orpha_df[['disease_mondo', 'gene_hgnc']].apply(tuple, axis=1))
                sets['ORPHA'] = ORPHA_tuples
            except Exception as e:
                print(f"Error processing ORPHA data: {e}")

        if omim_df is not None and not omim_df.empty:
            try:
                omim_df, before_mapping, after_mapping, missing_diseases, missing_genes = prep_g2d_omim(omim_df,
                                                                                                        disease_mappings_df,
                                                                                                        gene_mappings_df)
                analysis_summaries['omim'] = analysis_summary(before_mapping, after_mapping)
                missing_mappings['omim'] = {'missing_diseases': missing_diseases, 'missing_genes': missing_genes}
                OMIM_tuples = set(omim_df[['disease_mondo', 'gene_hgnc']].apply(tuple, axis=1))
                sets['OMIM'] = OMIM_tuples
            except Exception as e:
                print(f"Error processing OMIM data: {e}")

        if gencc_df is not None and not gencc_df.empty:
            try:
                gencc_df, gencc_before, gencc_after = prep_g2d_gencc(gencc_df)
                analysis_summaries['gencc'] = analysis_summary(gencc_before, gencc_after)
                GENCC_tuples = set(gencc_df[['disease_curie', 'gene_curie']].apply(tuple, axis=1))
                sets['GENCC'] = GENCC_tuples
            except Exception as e:
                print(f"Error processing GENCC data: {e}")

        try:
            analysis_summaries_file_path = os.path.join(output_folder_path, 'g2d_analysis_summary.txt')
            save_analysis_summaries(analysis_summaries, analysis_summaries_file_path)
            print(f"G2D Analysis summaries saved to {analysis_summaries_file_path}")
        except Exception as e:
            print(f"Error saving analysis summaries: {e}")

        try:
            missing_mappings_file_path = os.path.join(output_folder_path, 'g2d_missing_mappings.xlsx')
            save_missing_mappings(missing_mappings, missing_mappings_file_path)
            print(f"G2D missing mappings saved to {missing_mappings_file_path}")
        except Exception as e:
            print(f"Error saving missing mappings: {e}")

        try:
            g2d_venn_file_path = os.path.join(output_folder_path, 'g2d_venn.png')
            save_g2d_venn(sets, g2d_venn_file_path)
            print(f"G2D VENN plot saved to {g2d_venn_file_path}")
        except Exception as e:
            print(f"Error saving G2D VENN plot: {e}")

        try:
            g2d_edges_by_sources_file_path = os.path.join(output_folder_path, 'g2d_edges_by_sources.xlsx')
            g2d_edges_by_sources(sets, g2d_edges_by_sources_file_path)
            print(f"G2D edges by sources saved to {g2d_edges_by_sources_file_path}")
        except Exception as e:
            print(f"Error saving G2D edges by sources: {e}")

        try:
            g2d_differences_file_path = os.path.join(output_folder_path, 'g2d_differences.xlsx')
            g2d_edges_differences(sets, g2d_differences_file_path)
            print(f"G2D Differences saved to {g2d_differences_file_path}")
        except Exception as e:
            print(f"Error saving G2D differences: {e}")

        try:
            consolidated_file_path = os.path.join(output_folder_path, 'g2d_consolidated_edges.xlsx')
            consolidate_g2d_edges(sets, consolidated_file_path)
            print(f"G2D Consolidated edges saved to {consolidated_file_path}")
        except Exception as e:
            print(f"Error saving consolidated edges: {e}")

        try:
            upset_plot_file_path = os.path.join(output_folder_path, 'g2d_upset_plot.png')
            plot_upset(sets, upset_plot_file_path)
            print(f"G2D UpSet plot saved to {upset_plot_file_path}")
        except Exception as e:
            print(f"Error saving UpSet plot: {e}")

        try:
            unique_edges_file_path = os.path.join(output_folder_path, 'g2d_unique_edges.xlsx')
            consolidate_unique_g2d_edges(sets, unique_edges_file_path)
            print(f"G2D Unique edges saved to {unique_edges_file_path}")
        except Exception as e:
            print(f"Error saving unique edges: {e}")

    except Exception as e:
        print(f"An unexpected error occurred: {e}")


if __name__ == "__main__":
    main()

import pandas as pd
import xml.etree.ElementTree as ET
import os
import shutil
import sys

from oaklib import get_adapter
from venny4py.venny4py import *
#from venn import venn2, venn4
import matplotlib.pyplot as plt

def read_g2d_files(folder_path):
    files = os.listdir(folder_path)
    for file in files:
        file_path = os.path.join(folder_path, file)
        if file == 'gencc-submissions.tsv':
            gencc_df = pd.read_csv(file_path, sep='\t')
        elif file == 'omim-morbidmap.txt':
            omim_df = pd.read_csv(file_path, sep='\t', skiprows=3)
        elif file == 'hpoa-genes_to_disease.txt':
            hpoa_df = pd.read_csv(file_path, sep='\t')
        elif file == 'orpha-en_product6.xml':
            orpha_tree = ET.parse(file_path)
    return gencc_df, omim_df, hpoa_df, orpha_tree

def read_mapping_files(folder_path):
    files = os.listdir(folder_path)
    for file in files:
        file_path = os.path.join(folder_path, file)
        if file == 'mondo_exactmatch_omim.sssom.tsv':
            omim2mondo_df = pd.read_csv(file_path, sep='\t', skiprows=17)
        elif file == 'mondo_exactmatch_orphanet.sssom.tsv':
            orpha2mondo_df = pd.read_csv(file_path, sep='\t', skiprows=14)
            orpha2mondo_df["object_id"] = orpha2mondo_df["object_id"].str.replace('Orphanet', 'ORPHA')
    return omim2mondo_df, orpha2mondo_df

def prep_g2d_gencc(gencc_df):
    # Remove rows with null values
    gencc_df = gencc_df[gencc_df['disease_curie'].notna()]
    gencc_df = gencc_df[gencc_df['gene_curie'].notna()]
    print("GENCC Summary")
    print("~~~~~~~~~~~~~~")
    print("Disease count:", len(gencc_df['disease_curie']))
    print("Disease count (unique):", len(gencc_df['disease_curie'].unique()))
    print("Gene count:", len(gencc_df['gene_curie']))
    print("Gene count (unique):", len(gencc_df['gene_curie'].unique()))
    print("G2D association count:", len(gencc_df.groupby(['disease_curie', 'gene_curie']).size()))
    print("###########################################")
    return gencc_df

def prep_g2d_omim(omim_df, omim2mondo_df, orpha2mondo_df):
    ## parse columns "# Phenotype"  and "Gene Symbols" to get disease id and gene symbols ###
    omim_df["# Phenotype"] = omim_df["# Phenotype"].str.split(',').str[-1]
    omim_df["# Phenotype"] = omim_df["# Phenotype"].str[1:7]
    omim_df["Gene Symbols"] = omim_df["Gene Symbols"].str.split(',').str[0]

    ### Remove rows with null values ###
    omim_df['# Phenotype'] = pd.to_numeric(omim_df['# Phenotype'], errors='coerce')
    omim_df = omim_df[omim_df['# Phenotype'].notna()]
    omim_df = omim_df.copy()
    omim_df['# Phenotype'] = omim_df['# Phenotype'].astype(int)
    omim_df = omim_df[omim_df['Gene Symbols'].notna()]

    print("OMIM Summary")
    print("~~~~~~~~~~~~~~")
    print("Before mapping:")
    print("Disease count:", len(omim_df['# Phenotype']))
    print("Disease count (unique):", len(omim_df['# Phenotype'].unique()))
    print("Gene count:", len(omim_df['Gene Symbols']))
    print("Gene count (unique):", len(omim_df['Gene Symbols'].unique()))
    print("G2D association count:", len(omim_df.groupby(['# Phenotype', 'Gene Symbols']).size()))

    # prefix OMIM
    omim_df["# Phenotype"] = 'OMIM:' + omim_df["# Phenotype"].astype(str)

    # disease mapping
    print("OMIM disease mapping in progress")
    omim_df["disease_mondo"] = disease_mondo_mapping(omim_df["# Phenotype"], omim2mondo_df, orpha2mondo_df)
    print("OMIM disease mapping complete")

    # gene mapping
    print("OMIM gene mapping in progress")

    omim_df["gene_hgnc"] = gene_hgnc_mapping(omim_df["Gene Symbols"].tolist())
    omim_df["gene_hgnc"] = omim_df["gene_hgnc"].str.replace('hgnc', 'HGNC')
    print("OMIM gene mapping complete")

    ### Remove rows with null values (after mapping)###
    omim_df = omim_df[omim_df['disease_mondo'].notna()]
    omim_df = omim_df[omim_df['gene_hgnc'].notna()]

    print("After mapping:")
    print("Disease count:", len(omim_df['disease_mondo']))
    print("Disease count (unique):", len(omim_df['disease_mondo'].unique()))
    print("Gene count:", len(omim_df['gene_hgnc']))
    print("Gene count (unique):", len(omim_df['gene_hgnc'].unique()))
    print("G2D association count:", len(omim_df.groupby(['disease_mondo', 'gene_hgnc']).size()))
    print("###########################################")
    return omim_df

def prep_g2d_hpoa(hpoa_df, omim2mondo_df, orpha2mondo_df):

    ### Remove rows with null values ###
    hpoa_df = hpoa_df[hpoa_df['disease_id'].notna()]
    hpoa_df = hpoa_df[hpoa_df['gene_symbol'].notna()]

    print("HPOA Summary")
    print("~~~~~~~~~~~~~~")
    print("Before mapping:")
    print("Disease count:", len(hpoa_df['disease_id']))
    print("Disease count (unique):", len(hpoa_df['disease_id'].unique()))
    print("Gene count:", len(hpoa_df['gene_symbol']))
    print("Gene count (unique):", len(hpoa_df['gene_symbol'].unique()))
    print("G2D association count:", len(hpoa_df.groupby(['disease_id', 'gene_symbol']).size()))

    # disease mapping
    print("HPOA disease mapping in progress")
    hpoa_df["disease_mondo"] = disease_mondo_mapping(hpoa_df["disease_id"], omim2mondo_df, orpha2mondo_df)
    print("HPOA disease mapping complete")

    # gene mapping
    print("HPOA gene mapping in progress")

    hpoa_df["gene_hgnc"] = gene_hgnc_mapping(hpoa_df["gene_symbol"].tolist())
    hpoa_df["gene_hgnc"] = hpoa_df["gene_hgnc"].str.replace('hgnc', 'HGNC')
    print("HPOA gene mapping complete")

    ### Remove rows with null values (after mapping)###
    hpoa_df = hpoa_df[hpoa_df['disease_mondo'].notna()]
    hpoa_df = hpoa_df[hpoa_df['gene_hgnc'].notna()]

    print("After mapping:")
    print("Disease count:", len(hpoa_df['disease_mondo']))
    print("Disease count (unique):", len(hpoa_df['disease_mondo'].unique()))
    print("Gene count:", len(hpoa_df['gene_hgnc']))
    print("Gene count (unique):", len(hpoa_df['gene_hgnc'].unique()))
    print("G2D association count:", len(hpoa_df.groupby(['disease_mondo', 'gene_hgnc']).size()))
    print("###########################################")

    return hpoa_df

def prep_g2d_orpha(orpha_tree, omim2mondo_df, orpha2mondo_df):

    ### Parsing XML file -> dataframe ###
    #root = orpha_tree.getroot()
    df_cols = ["Disorder_ID", "Disorder_Name", "Orpha_Code", "gene_hgnc"]
    rows = []
    for node in orpha_tree.iter('Disorder'):
        disorder_id = node.attrib.get('id')
        disorder_name = node.find("Name").text
        orphacode = node.find("OrphaCode").text
        if (node.find("DisorderGeneAssociationList/DisorderGeneAssociation/Gene/ExternalReferenceList/ExternalReference[Source='HGNC']")):
            genehgnc = node.find("DisorderGeneAssociationList/DisorderGeneAssociation/Gene/ExternalReferenceList/ExternalReference[Source='HGNC']/Reference").text
        rows.append({"Disorder_ID": disorder_id, "Disorder_Name": disorder_name, "Orpha_Code": orphacode, "gene_hgnc": genehgnc})
    orpha_df = pd.DataFrame(rows, columns=df_cols)

    ### Remove rows with null values ###
    orpha_df = orpha_df[orpha_df['Orpha_Code'].notna()]
    orpha_df = orpha_df[orpha_df['gene_hgnc'].notna()]

    print("ORPHA Summary")
    print("~~~~~~~~~~~~~~")
    print("Before mapping:")
    print("Disease count:", len(orpha_df['Orpha_Code']))
    print("Disease count (unique):", len(orpha_df['Orpha_Code'].unique()))
    print("Gene count:", len(orpha_df['gene_hgnc']))
    print("Gene count (unique):", len(orpha_df['gene_hgnc'].unique()))
    print("G2D association count:", len(orpha_df.groupby(['Orpha_Code', 'gene_hgnc']).size()))

    ### prefix HGNC and ORPHA ###
    orpha_df["gene_hgnc"] = 'HGNC:' + orpha_df["gene_hgnc"].astype(str)
    orpha_df["Orpha_Code"] = 'ORPHA:' + orpha_df["Orpha_Code"].astype(str)

    # disease mapping
    print("ORPHA disease mapping in progress")
    orpha_df["disease_mondo"] = disease_mondo_mapping(orpha_df["Orpha_Code"], omim2mondo_df, orpha2mondo_df)
    print("ORPHA disease mapping complete")

    ### Remove rows with null values (after mapping) ###
    orpha_df = orpha_df[orpha_df['disease_mondo'].notna()]
    orpha_df = orpha_df[orpha_df['gene_hgnc'].notna()]

    print("After mapping:")
    print("Disease count:", len(orpha_df['disease_mondo']))
    print("Disease count (unique):", len(orpha_df['disease_mondo'].unique()))
    print("G2D association count:", len(orpha_df.groupby(['disease_mondo', 'gene_hgnc']).size()))
    print("###########################################")

    return orpha_df


def disease_mondo_mapping(ids_to_map, omim2mondo_df, orpha2mondo_df):
    mondo_mappings = []
    for id_to_map in ids_to_map:
        if (omim2mondo_df.loc[omim2mondo_df['object_id'] == id_to_map, 'subject_id'].tolist()):
            mondo_mappings.append(
                *omim2mondo_df.loc[omim2mondo_df['object_id'] == id_to_map, 'subject_id'].tolist())
        elif (orpha2mondo_df.loc[orpha2mondo_df['object_id'] == id_to_map, 'subject_id'].tolist()):
            mondo_mappings.append(
                *orpha2mondo_df.loc[orpha2mondo_df['object_id'] == id_to_map, 'subject_id'].tolist())
        else:
            mondo_mappings.append(None)
    return mondo_mappings

def gene_hgnc_mapping(genes_to_map):
    adapter = get_adapter("sqlite:obo:hgnc")
    hgnc_mappings = []
    for gene in genes_to_map:
        if adapter.curies_by_label(gene):
            for curie in adapter.basic_search(gene):
                hgnc_mappings.append(curie)
        else:
            hgnc_mappings.append(None)
    return hgnc_mappings

def venn4_plot_g2d(gencc, omim, hpoa, orpha, output_data_folder_path):
    sets = {
        'GENCC-submissions.tsv': set(gencc[['disease_curie', 'gene_curie']].apply(tuple, axis=1)),
        'OMIM-morbidmap.txt': set(omim[['disease_mondo', 'gene_hgnc']].apply(tuple, axis=1)),
        'HPOA-genes_to_disease.txt': set(hpoa[['disease_mondo', 'gene_hgnc']].apply(tuple, axis=1)),
        'ORPHA-en_product6.xml': set(orpha[['disease_mondo', 'gene_hgnc']].apply(tuple, axis=1))}
    venny4py(sets=sets)
    shutil.move('Intersections_4.txt', os.path.join(output_data_folder_path, 'Intersections_4.txt'))
    shutil.move('Venn_4.png', os.path.join(output_data_folder_path, 'Venn_4.png'))
    return plt.show()
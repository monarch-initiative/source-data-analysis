# source-data-analysis
Repository of analyses of data sources



# Gene to Disease (G2D) Data Source Analysis 

This project is focused on analyzing gene to disease data sources. This README file provides instructions on how to set up and run the project locally, as well as details about the required data sources and the expected output. 

## Installation and Setup 

To install and activate the virtual environment locally, please follow these steps: 
1. Clone the project repository from GitHub. 
2. Navigate to the project directory in your terminal. 
3. Create a virtual environment using the command: `python -m venv venv`. 
4. Activate the virtual environment: 
- For Windows: `venv\Scripts\activate` 
- For Mac/Linux: `source venv/bin/activate` 

## Package Installation 

To install the required packages for the project, use the following command: 

pip install -r requirements.txt

## Data Source 

For the G2D analysis, four files have been identified and should be placed in the `input/g2d` folder within the project: 
- `gencc-submission.tsv` 
- `hpoa-genes_to_disease.txt` 
- `omim-morbidmap.txt` 
- `orphan-en_product6.xml` 

Please ensure that the files are named exactly as mentioned above. 

## Mapping 

Some of the data sources require disease and gene mappings as they are not normalized to MONDO and HGNC. For disease mappings, the `mondo_exactmatch_omim.sssom.tsv` and `mondo_exactmatch_orphanet.sssom.tsv` files from SSSOM are used. For gene mappings, the Oaklib adapter is used. 

## Execution and Summary 

To run the project and obtain a summary for each data source, execute the following command in the project terminal: 

python src/g2d.py

This will parse the data sources and provide information on the number of diseases, genes, and G2D associations. 

## Output 

The project generates a Venn diagram illustrating the uniqueness of G2D associations. Additionally, the intersection results can be found in the `output/g2d` folder.

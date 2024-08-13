# Directory
INPUT_DIR := data/input/g2d

# Download GENCC data
$(INPUT_DIR)/gencc-submissions.tsv: $(INPUT_DIR)
	curl --silent --output $(INPUT_DIR)/gencc-submissions.tsv https://search.thegencc.org/download/action/submissions-export-tsv

# Download MONDO SSSOM mappings
$(INPUT_DIR)/mondo.sssom.tsv: $(INPUT_DIR)
	curl --silent --output $(INPUT_DIR)/mondo.sssom.tsv https://data.monarchinitiative.org/mappings/latest/mondo.sssom.tsv

# Download Gene Mappings SSSOM
$(INPUT_DIR)/gene_mappings.sssom.tsv: $(INPUT_DIR)
	curl --silent --output $(INPUT_DIR)/gene_mappings.sssom.tsv https://data.monarchinitiative.org/mappings/latest/gene_mappings.sssom.tsv

# Download MIM2GENE data
$(INPUT_DIR)/mim2gene_medgen.txt: $(INPUT_DIR)
	curl --silent --output $(INPUT_DIR)/mim2gene_medgen.txt https://ftp.ncbi.nih.gov/gene/DATA/mim2gene_medgen

# Download ORPHADATA
$(INPUT_DIR)/en_product6.xml: $(INPUT_DIR)
	curl --silent --output $(INPUT_DIR)/en_product6.xml https://www.orphadata.com/data/xml/en_product6.xml

# Create input directory if it does not exist
$(INPUT_DIR):
	mkdir -p $(INPUT_DIR)

# Rule to download all files
download_g2d_sources: $(INPUT_DIR)/gencc-submissions.tsv $(INPUT_DIR)/mondo.sssom.tsv $(INPUT_DIR)/gene_mappings.sssom.tsv $(INPUT_DIR)/mim2gene_medgen.txt $(INPUT_DIR)/en_product6.xml

# Run g2d_analysis.py script
analyze_g2d_sources:
	python src/g2d_analysis.py

.PHONY: download_g2d_sources analyze_g2d_sources

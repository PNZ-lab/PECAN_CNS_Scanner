# PECAN_CNS_Scanner
This script is fed PeCan expression data and a file with associated clinical data. It then calculates and plots the level of expression as stratified by the Central Nervous System (CNS) invasion category:
- CNS1 – no blasts in the CSF
- CNS2 – WBC count is less than 5/mL with blasts in the CSF
- CNS3 – WBC count is 5/mL or greater with blasts in the CSF or there are signs that leukemia has spread to the CNS

## Plot a single gene
Run cell "3. Plot a single gene" to analayze a gene and generate a graph for it.

There are three ways to collapse the data:
1. All three groups separate
2. CNS1 and CNS2 collapsed
3. CNS2 and CNS3 collapsed

PECAN_CNS_scanner.py can do all depending on how the Grapher function at the bottom of cell 3 is executed:
1. Grapher(gene, collapse_type=None)
2. Grapher(gene, collapse_type='collapse_1_2')
3. Grapher(gene, collapse_type='collapse_2_3')

#### collapse_type='None': 

<img width="543" alt="image" src="https://github.com/user-attachments/assets/a973a131-972b-4c46-abec-e25edeb83154" />

#### collapse_type='collapse_1_2':

<img width="543" alt="image" src="https://github.com/user-attachments/assets/3bf7ca20-9cdd-46f6-90d0-5b72edab6e30" />

#### collapse_type='collapse_2_3':

<img width="543" alt="image" src="https://github.com/user-attachments/assets/490b7b95-c8a8-4544-bae5-3dda9c7b8ed0" />

## Analyze all genes and plot most significant
Run cell "4. Plot the top n most significant genes" to run the analyze for every gene in the dataset, plot the most significant, and create a csv file containing significance for all genes.
1. significant_genes = AnalyzeAllGenes("CNS1_vs_CNS3")
2. significant_genes = AnalyzeAllGenes("CNS1_2_vs_CNS3")
3. significant_genes = AnalyzeAllGenes("CNS1_vs_CNS2_3")

#### csv file:

<img width="543" alt="image" src="https://github.com/user-attachments/assets/39390bba-626f-41b8-b953-826ccb3d0c0e" />






    

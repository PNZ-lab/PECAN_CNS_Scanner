#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ttest_ind

# Load the data files
PECAN_in  = '/Volumes/cmgg_pnlab/Kasper/Data/Interesting_Lists/Expression_Pecan.txt'  # PECAN FPKM data
clin_data = '/Users/kasperthorhaugechristensen/Library/CloudStorage/OneDrive-UGent/PNG lab/Kasper/PECAN_scanner/PECAN_clinical_reference.tsv'
out_dir   = '/Users/kasperthorhaugechristensen/Desktop/Dumpbox/CNS3_sign_and_substantial/'  # Output directory
df_PECAN  = pd.read_csv(PECAN_in, sep='\t')
clin_df   = pd.read_csv(clin_data, sep='\t')

#%%
# # Optional parameters
subanalysis_col  = 'CNS_at_Dx'

def WriteFile(name):
    plt.savefig(os.path.join(out_dir, name))
    print('File created: %s' % (os.path.join(out_dir, name)))

#%%
from tqdm import tqdm
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
import pandas as pd
import numpy as np

def AnalyzeAllGenes(comparison_type="CNS1_2_vs_CNS3"):
    gene_list = df_PECAN['Gene'].tolist()  # Get all gene names from the PECAN data
    gene_stats = []  # List to store gene statistics

    # Use tqdm to track progress over the list of genes
    for gene in tqdm(gene_list, desc="Analyzing Genes", unit="gene"):
        values = df_PECAN.loc[df_PECAN['Gene'] == gene].iloc[0, 1:].tolist()
        pecan_samples = df_PECAN.columns[1:].tolist()

        data = {
            'Sample': [],
            'Expression': [],
            'Group': []
        }

        # Filter out samples with a '.' in CNS_at_Dx
        clin_df_filtered = clin_df[clin_df['CNS_at_Dx'] != '.']

        for sample in pecan_samples:
            # Check if the sample exists in the filtered clinical data
            match = clin_df_filtered[clin_df_filtered['RNAseq_id_D'] == sample]
            if not match.empty:
                stage = match[subanalysis_col].values[0]

                # Determine the group based on the comparison type
                group = None
                if comparison_type == "CNS1_2_vs_CNS3":
                    if stage in ['CNS 1', 'CNS 2']:
                        group = 'CNS 1/2'
                    elif stage == 'CNS 3':
                        group = 'CNS 3'
                elif comparison_type == "CNS1_vs_CNS2_3":
                    if stage == 'CNS 1':
                        group = 'CNS 1'
                    elif stage in ['CNS 2', 'CNS 3']:
                        group = 'CNS 2/3'
                elif comparison_type == "CNS1_vs_CNS3":
                    if stage == 'CNS 1':
                        group = 'CNS 1'
                    elif stage == 'CNS 3':
                        group = 'CNS 3'

                # Only add to data if the group is assigned, excluding unwanted stages
                if group is not None:
                    data['Sample'].append(sample)
                    data['Expression'].append(values[pecan_samples.index(sample)])
                    data['Group'].append(group)

        # Convert the collected data into a DataFrame
        plot_df = pd.DataFrame(data)

        # Retrieve values for the groups
        group_1_values = plot_df[plot_df['Group'] == 'CNS 1/2']['Expression'] if comparison_type == "CNS1_2_vs_CNS3" else plot_df[plot_df['Group'] == 'CNS 1']['Expression']
        group_2_values = plot_df[plot_df['Group'] == 'CNS 3']['Expression'] if comparison_type == "CNS1_2_vs_CNS3" or comparison_type == "CNS1_vs_CNS3" else plot_df[plot_df['Group'] == 'CNS 2/3']['Expression']

        # Perform t-test if both groups have values
        if not group_1_values.empty and not group_2_values.empty:
            t_stat, p_value = ttest_ind(group_1_values, group_2_values, equal_var=False)
            
            # Calculate means for both groups
            mean_group_1 = group_1_values.mean()
            mean_group_2 = group_2_values.mean()
            
            # Calculate the ratio between the means
            mean_ratio = mean_group_1 / mean_group_2 if mean_group_2 != 0 else float('inf')

            # Append gene statistics to the list
            gene_stats.append({
                'Gene': gene,
                'Mean Group 1': mean_group_1,
                'Mean Group 2': mean_group_2,
                'Mean Ratio': mean_ratio,
                'p-value': p_value,
                'Group 1 Values': group_1_values.tolist(),  # Store group 1 values as a list
                'Group 2 Values': group_2_values.tolist()   # Store group 2 values as a list
            })

    # Convert the list of gene statistics to a DataFrame
    stats_df = pd.DataFrame(gene_stats)

    # Filter out rows with NaN or infinite values in the p-value column
    filtered_stats_df = stats_df.replace([np.inf, -np.inf], np.nan).dropna(subset=['p-value'])

    # Apply Benjamini-Hochberg adjustment to the valid p-values
    _, pvals_corrected, _, _ = multipletests(filtered_stats_df['p-value'], method='fdr_bh')

    # Assign the adjusted p-values back to the original DataFrame
    stats_df.loc[filtered_stats_df.index, 'padj'] = pvals_corrected

    # Save the DataFrame to a CSV file
    csv_file = '/Users/kasperthorhaugechristensen/Desktop/%s.csv' % comparison_type
    stats_df.to_csv(csv_file, index=False)
    
    print(f"Gene statistics saved to {csv_file}")

    return gene_stats



def PlotTopGenes(gene_stats, top_n=10):
    # Sort the genes based on p-value and select the top N genes
    top_genes = sorted(gene_stats, key=lambda x: x['p-value'])[:top_n]

    for gene_stat in top_genes:
        gene = gene_stat['Gene']
        p_value = gene_stat['p-value']
        
        # Retrieve expression data for the gene
        values = df_PECAN.loc[df_PECAN['Gene'] == gene].iloc[0, 1:].tolist()
        pecan_samples = df_PECAN.columns[1:].tolist()

        data = {
            'Sample': [],
            'Expression': [],
            'Group': []
        }

        # Filter out samples with a '.' in CNS_at_Dx
        clin_df_filtered = clin_df[clin_df['CNS_at_Dx'] != '.']

        for sample in pecan_samples:
            # Check if the sample exists in the filtered clinical data
            match = clin_df_filtered[clin_df_filtered['RNAseq_id_D'] == sample]
            if not match.empty:
                stage = match[subanalysis_col].values[0]
                data['Sample'].append(sample)
                data['Expression'].append(values[pecan_samples.index(sample)])
                # Collapse CNS 1 and CNS 2 into 'CNS 1/2'
                if stage in ['CNS 1', 'CNS 2']:
                    data['Group'].append('CNS 1/2')
                else:
                    data['Group'].append('CNS 3')

        # Create a DataFrame from the collected data
        plot_df = pd.DataFrame(data)

        # Plotting
        plt.figure(figsize=(8, 6), dpi=200)
        sns.boxplot(x='Group', y='Expression', data=plot_df, color='white', linewidth=2, showfliers=False)
        sns.stripplot(x='Group', y='Expression', data=plot_df, color='black', size=5, alpha=0.5, jitter=True)

        # Create a title summarizing the p-value
        plt.title('Expression levels of %s (CNS 1/2 vs CNS 3): p=%.3f' % (gene, p_value), fontsize=22)
        plt.xlabel('Group', fontsize=18)
        plt.ylabel('%s (FPKM)' % gene, fontsize=18)
        plt.tick_params(axis='both', labelsize=16)

        # Save the file with 'collapsed' in the name
        file_name = 'PECAN_boxplot_collapsed_%s.png' % gene
        WriteFile(file_name)
        plt.show()
        
# Run the analysis and plot top 10 significant genes
significant_genes = AnalyzeAllGenes("CNS1_vs_CNS3")
PlotTopGenes(significant_genes, top_n=10)  # Plot top 10 genes

#%%

#%%

from tqdm import tqdm
from scipy.stats import ttest_ind
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def Grapher(gene, collapse_type=None):

    values = df_PECAN.loc[df_PECAN['Gene'] == gene].iloc[0, 1:].tolist()
    pecan_samples = df_PECAN.columns[1:].tolist()

    data = {
        'Sample': [],
        'Expression': [],
        'Group': []
    }

    clin_df_filtered = clin_df[clin_df['CNS_at_Dx'] != '.']

    for sample in pecan_samples:
        match = clin_df_filtered[clin_df_filtered['RNAseq_id_D'] == sample]
        if not match.empty:
            stage = match[subanalysis_col].values[0]

            # Collapse groups based on the specified collapse type
            if collapse_type == 'collapse_1_2' and stage in ['CNS 1', 'CNS 2']:
                stage = 'CNS 1/2'
            elif collapse_type == 'collapse_2_3' and stage in ['CNS 2', 'CNS 3']:
                stage = 'CNS 2/3'
            
            data['Sample'].append(sample)
            data['Expression'].append(values[pecan_samples.index(sample)])
            data['Group'].append(stage)

    # Create a DataFrame from the collected data
    plot_df = pd.DataFrame(data)

    # Perform pairwise comparisons between groups
    groups = plot_df['Group'].unique()
    p_values = {}
    
    for i in range(len(groups)):
        for j in range(i + 1, len(groups)):
            group1 = groups[i]
            group2 = groups[j]
            group1_values = plot_df[plot_df['Group'] == group1]['Expression']
            group2_values = plot_df[plot_df['Group'] == group2]['Expression']
            if not group1_values.empty and not group2_values.empty:
                t_stat, p_value = ttest_ind(group1_values, group2_values, equal_var=False)
                p_values[f'{group1} vs {group2}'] = p_value

    # Plotting
    plt.figure(figsize=(10, 6), dpi=200)
    sns.boxplot(x='Group', y='Expression', data=plot_df, color='white', linewidth=2, showfliers=False)
    sns.stripplot(x='Group', y='Expression', data=plot_df, color='black', size=5, alpha=0.5, jitter=True)

    # Create a title summarizing p-values
    title_text = f'Expression levels of {gene} by CNS group'
    for comparison, p_value in p_values.items():
        title_text += f'\n{comparison}: p=%.3f' % p_value

    plt.title(title_text, fontsize=22)
    plt.xlabel('Group', fontsize=18)
    plt.ylabel(f'{gene} (FPKM)', fontsize=18)
    plt.tick_params(axis='both', labelsize=16)

    # Save the file
    file_name = f'PECAN_boxplot_{gene}_{collapse_type}.png'
    WriteFile(file_name)
    plt.show()

gene_list =[
    "MVK",
    "PMVK",
    "MVD",
    "IDI1",
    "GGPS1",
    "FDPS",
    "FDFT1", #"SQS",
    "SQLE",
    "LSS", #"OSC",
    "CYP51A1",
    "DHCR24",
    "MSMO1", #"SC4MOL",
    "NSDHL",
    "HSD17B7",
    "EBP",
    "SC5D",
    "DHCR7",
    "DHCR24",
    "LDLR",
    "METTL3",
    "METTL14",
    "ALKBH5",
    "FTO",
    "YTHDF1",
    "YTHDF2",
    "HNRNPC",
    "IGF2BP2"
]

gene_list = ['NAMPT', 'NAPRT1', 'IDO1', 'DHFR']
# gene_list = ['PRPF8', 'SRRM1', 'SRRM2', 'ACIN1', 'RNPS1', 'CLK1', 'CLK2', 'CLK3', 'CLK4']

for gene in gene_list:
    try:
        Grapher(gene, collapse_type=None)  # Set collapse_type to None for separate plots for CNS 1, CNS 2, and CNS 3
        Grapher(gene, collapse_type='collapse_1_2')  # Set collapse_type to None for separate plots for CNS 1, CNS 2, and CNS 3
        Grapher(gene, collapse_type='collapse_2_3')  # Set collapse_type to None for separate plots for CNS 1, CNS 2, and CNS 3
    except:
        print(f'error for gene: {gene}')

#%%
from tqdm import tqdm
from scipy.stats import ttest_ind
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def AnalyzeAllGenes():
    gene_list = df_PECAN['Gene'].tolist()  # Get all gene names from the PECAN data
    gene_stats = []  # List to store gene statistics

    # Use tqdm to track progress over the list of genes
    for gene in tqdm(gene_list, desc="Analyzing Genes", unit="gene"):
        values = df_PECAN.loc[df_PECAN['Gene'] == gene].iloc[0, 1:].tolist()
        pecan_samples = df_PECAN.columns[1:].tolist()

        data = {
            'Sample': [],
            'Expression': [],
            'Group': []
        }

        # Filter out samples with a '.' in CNS_at_Dx and exclude CNS 2
        clin_df_filtered = clin_df[(clin_df['CNS_at_Dx'] != '.') & (clin_df['CNS_at_Dx'] != 'CNS 2')]

        for sample in pecan_samples:
            # Check if the sample exists in the filtered clinical data
            match = clin_df_filtered[clin_df_filtered['RNAseq_id_D'] == sample]
            if not match.empty:
                stage = match[subanalysis_col].values[0]
                
                # Only include CNS 1 and CNS 3
                if stage in ['CNS 1', 'CNS 3']:
                    data['Sample'].append(sample)
                    data['Expression'].append(values[pecan_samples.index(sample)])
                    data['Group'].append(stage)

        # Create a DataFrame from the collected data
        plot_df = pd.DataFrame(data)

        # Calculate p-values between CNS 1 and CNS 3
        cns_1_values = plot_df[plot_df['Group'] == 'CNS 1']['Expression']
        cns_3_values = plot_df[plot_df['Group'] == 'CNS 3']['Expression']
        
        if not cns_1_values.empty and not cns_3_values.empty:
            # Perform t-test
            t_stat, p_value = ttest_ind(cns_1_values, cns_3_values, equal_var=False)
            
            # Calculate medians for both groups
            median_cns_1 = cns_1_values.median()
            median_cns_3 = cns_3_values.median()

            # Append gene statistics to the list
            gene_stats.append({
                'Gene': gene,
                'Median CNS 1': median_cns_1,
                'Median CNS 3': median_cns_3,
                'p-value': p_value
            })

    # Convert the list of gene statistics to a DataFrame
    stats_df = pd.DataFrame(gene_stats)

    # Save the DataFrame to a CSV file
    csv_file = '/Users/kasperthorhaugechristensen/Desktop/Freya_PECAN_CNS1_vs_CNS3.csv'
    stats_df.to_csv(csv_file, index=False)
    
    print(f"Gene statistics saved to {csv_file}")

    return gene_stats


def PlotTopGenes(gene_stats, top_n=10):
    # Sort the genes based on p-value and select the top N genes
    top_genes = sorted(gene_stats, key=lambda x: x['p-value'])[:top_n]

    for gene_stat in top_genes:
        gene = gene_stat['Gene']
        p_value = gene_stat['p-value']
        
        # Retrieve expression data for the gene
        values = df_PECAN.loc[df_PECAN['Gene'] == gene].iloc[0, 1:].tolist()
        pecan_samples = df_PECAN.columns[1:].tolist()

        data = {
            'Sample': [],
            'Expression': [],
            'Group': []
        }

        # Filter out samples with a '.' in CNS_at_Dx and exclude CNS 2
        clin_df_filtered = clin_df[(clin_df['CNS_at_Dx'] != '.') & (clin_df['CNS_at_Dx'] != 'CNS 2')]

        for sample in pecan_samples:
            # Check if the sample exists in the filtered clinical data
            match = clin_df_filtered[clin_df_filtered['RNAseq_id_D'] == sample]
            if not match.empty:
                stage = match[subanalysis_col].values[0]
                
                # Only include CNS 1 and CNS 3
                if stage in ['CNS 1', 'CNS 3']:
                    data['Sample'].append(sample)
                    data['Expression'].append(values[pecan_samples.index(sample)])
                    data['Group'].append(stage)

        # Create a DataFrame from the collected data
        plot_df = pd.DataFrame(data)

        # Plotting
        plt.figure(figsize=(8, 6), dpi=200)
        sns.boxplot(x='Group', y='Expression', data=plot_df, color='white', linewidth=2, showfliers=False)
        sns.stripplot(x='Group', y='Expression', data=plot_df, color='black', size=5, alpha=0.5, jitter=True)

        # Create a title summarizing the p-value
        plt.title('Expression levels of %s (CNS 1 vs CNS 3): p=%.3f' % (gene, p_value), fontsize=22)
        plt.xlabel('Group', fontsize=18)
        plt.ylabel('%s (FPKM)' % gene, fontsize=18)
        plt.tick_params(axis='both', labelsize=16)

        # Save the file with 'CNS1_vs_CNS3' in the name
        file_name = 'PECAN_boxplot_CNS1_vs_CNS3_%s.png' % gene
        WriteFile(file_name)
        plt.show()
        
# Run the analysis and plot top 10 significant genes
significant_genes = AnalyzeAllGenes()
PlotTopGenes(significant_genes, top_n=10)  # Plot top 10 genes

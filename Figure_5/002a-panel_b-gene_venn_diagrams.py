# create venn diagrams for the overlap of differentially expressed genes, differentially methylated genes, and differentially expressed proteins

import sys, getopt, time, os, subprocess

try:
    sys.path.insert(1, "/".join(os.path.realpath(__file__).split("/")[:-2]))
except:
    sys.path.insert(1, "/".join(os.getcwd().split("/")[:-1]))

import pandas as pd 
import numpy as np

import statsmodels.stats.multitest as multitest

help_message = '''
Failed
'''

def subset_genes_by_cutoff(pvals_df, cutoff=0.1):
    res_df = []

    # iterate through columns to get genes passing cutoffs
    for col in pvals_df.columns:
        # fdr correct
        mask = multitest.fdrcorrection(pvals_df[col].to_numpy())[1] < cutoff

        print(col, np.sum(mask))

        # add description and add gene names to df
        res_df.append(pd.DataFrame([pvals_df.index[mask], np.repeat(col, np.sum(mask))], index = ['gene', 'source']).transpose())

    return pd.concat(res_df)


def subset_genes_by_direction(genes_df, logfc_df):
    # iterate through gene sources to assign direction
    for source in np.unique(genes_df['source']):
        mask = genes_df['source'] == source
        
        # are genes upregulated
        genes_of_interest = genes_df.loc[mask,'gene']

        genes_df.loc[mask,'upregulated'] = logfc_df.loc[genes_of_interest,source].to_numpy() > 0

    return genes_df




def main():
    print('XXX-XXX.py')
    t1 = time.time()

    # load dfs
    pvals_df = pd.read_csv(combined_pval_file, sep='\t', index_col=0)
    logfc_df = pd.read_csv(combined_fc_file, sep='\t', index_col=0)

    # genes that pass fdr cutoff and assign directions
    pvals_df = subset_genes_by_cutoff(pvals_df)
    pvals_df = subset_genes_by_direction(pvals_df, logfc_df)

    # save to file
    pvals_df.to_csv(figure_data_file, sep='\t', index=False)
    

    # create figure
    cline = [rscript, figure_script, '--figure_data_file', figure_data_file, '--figure_file', figure_file]

    # create comparison plots
    print(" ".join(cline))
    subprocess.run(cline)

    print(round(time.time() - t1, 2))
    print('XXX-XXX.py COMPLETE')



if __name__ == "__main__":
    # r environment
    rscript = "Rscript"

    try:
        opts, args = getopt.getopt(sys.argv[1:],"",["combined_pval_file=", "combined_fc_file=", "figure_script=", "figure_data_file=", "figure_file=", "merged_fdr_file="])
    except getopt.GetoptError:
        print (help_message)
        print(sys.argv)
        sys.exit(2)
    for opt, arg in opts:
        # file of P-values from each source
        if opt in ("--combined_pval_file"):
            combined_pval_file = str(arg)
            combined_pval_file = '~/input_data/intermediate_files/figure5_pvals.tsv'
        # file of log2 fold changes from each source
        if opt in ("--combined_fc_file"):
            combined_fc_file = str(arg)
            combined_fc_file = '~/input_data/intermediate_files/figure5_fcs.tsv'

        # script to create figure
        if opt in ("--figure_script"):
            figure_script = str(arg)
            figure_script = '~/002b-panel_b-gene_venn_diagrams.R'
        # file to save figure data
        if opt in ("--figure_data_file"):
            figure_data_file = str(arg)
            figure_data_file = '~/input_data/intermediate_files/figure5_gene_venn_diagrams.tsv'
        # file to save figure
        if opt in ("--figure_file"):
            figure_file = str(arg)
            figure_file = '~/figure5_gene_venn_diagrams.pdf'

    main()




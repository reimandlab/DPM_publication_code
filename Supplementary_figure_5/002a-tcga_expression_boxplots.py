# visualize gene expression of genes of interest in TCGA data

import sys, getopt, time, os, subprocess

try:
    sys.path.insert(1, "/".join(os.path.realpath(__file__).split("/")[:-2]))
except:
    sys.path.insert(1, "/".join(os.getcwd().split("/")[:-1]))

import pandas as pd 
import numpy as np


help_message = '''
Failed
'''


def classify_and_melt(exp_df, mut_patients):
    # labels columns as mutant or wildtype
    exp_df.columns = exp_df.columns.str.split("-").str[:3].str.join("-")
    exp_df.columns = np.where(exp_df.columns.isin(mut_patients), "mutant", "wildtype")

    # melt
    exp_df = exp_df.reset_index().melt(id_vars = "index", var_name = "mutation_status", value_name = "value").rename(columns = {'index':'gene'})

    return exp_df


def main():
    print('XXX-XXX.py')
    t1 = time.time()

    # load methylation and probe df
    exp_df = pd.read_csv(gene_expression_file, sep='\t', index_col=0)

    # subset to genes of interest
    exp_df = exp_df.loc[np.isin(exp_df.index, genes_list),:]

    # remove R132H
    exp_df = exp_df.loc[:,exp_df.columns != "TCGA-06-0221"]
    exp_df = exp_df.loc[:,exp_df.columns != "TCGA-12-0818"]

    # label and melt
    exp_df = classify_and_melt(exp_df, mut_patients)

    print(exp_df)
    print(np.unique(exp_df['mutation_status'], return_counts=True))

    # save to file
    exp_df.to_csv(figure_data_file, sep='\t', index=False)

    # create figure
    cline = [rscript, figure_script, '--figure_data_file', figure_data_file, '--figure_stats_file', figure_stats_file, '--figure_file', figure_file]

    # create comparison plots
    print(" ".join(cline))
    subprocess.run(cline)


    print(round(time.time() - t1, 2))
    print('XXX-XXX.py COMPLETE')



if __name__ == "__main__":
    rscript = "Rscript"

    # mutant patients
    mut_patients = ['TCGA-16-1460', 'TCGA-06-5417', 'TCGA-26-1442', 'TCGA-14-1456', 'TCGA-06-2570', 'TCGA-06-6701', 'TCGA-27-2521', 'TCGA-16-0849', 'TCGA-06-6389', 'TCGA-16-0850', 'TCGA-19-1788', 'TCGA-19-A6J5', 'TCGA-14-4157', 'TCGA-06-A7TL', 'TCGA-19-2629', 'TCGA-06-0129', 'TCGA-12-0827', 'TCGA-32-4208', 'TCGA-14-1821', 'TCGA-02-2483', 'TCGA-06-0128', 'TCGA-12-0818', 'TCGA-06-0221']

    # list of genes of interest
    genes_list = ['SULF1', 'FGF2', 'WNT5A', 'NRXN1']


    try:
        opts, args = getopt.getopt(sys.argv[1:],"",["gene_expression_file=", "figure_script=", "figure_data_file=", "figure_stats_file=", "figure_file="])
    except getopt.GetoptError:
        print (help_message)
        print(sys.argv)
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("--gene_expression_file"):
            gene_expression_file = str(arg)
            gene_expression_file = '~/input_data/'

        if opt in ("--figure_script"):
            figure_script = str(arg)
            figure_script = '~/002b-tcga_expression_boxplots.R'
            
        if opt in ("--figure_data_file"):
            figure_data_file = str(arg)
            figure_data_file = '~/input_data/intermediate_files/tcga_gene_expression.tsv'
        if opt in ("--figure_stats_file"):
            figure_stats_file = str(arg)
            figure_stats_file = '~/input_data/intermediate_files/tcga_gene_expression_stats.tsv'
        if opt in ("--figure_file"):
            figure_file = str(arg)
            figure_file = '~/tcga_gene_expression.pdf'


    main()




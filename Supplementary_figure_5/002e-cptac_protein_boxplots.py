# visualzie the protein expression of genes of interest in the CPTAC cohort

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
    exp_df.columns = np.where(exp_df.columns.isin(mut_patients), "mutant", "wildtype")

    print(exp_df)

    # melt
    exp_df = exp_df.reset_index().melt(id_vars = "index", var_name = "mutation_status", value_name = "value").rename(columns = {'index':'gene'})

    return exp_df


def main():
    print('XXX-XXX.py')
    t1 = time.time()

    # load methylation and probe df
    exp_df = pd.read_csv(protein_expression_file)
    exp_df.index = exp_df['Gene'].to_numpy()
    exp_df = exp_df.iloc[:,1:]

    # subset to genes of interest
    exp_df = exp_df.loc[np.isin(exp_df.index, genes_list),:]

    # remove R132H
    exp_df = exp_df.loc[:,~np.isin(exp_df.columns, ["C3N-01368","C3L-03390"])]

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
    mut_patients = ["C3L-02041", "C3N-00663", "C3L-02984", "C3L-02708", "C3L-02900","C3N-03070"]

    # list of genes of interest
    genes_list = ['SULF1', 'FGF2', 'WNT5A', 'NRXN1']


    try:
        opts, args = getopt.getopt(sys.argv[1:],"",["protein_expression_file=", "figure_script=", "figure_data_file=", "figure_stats_file=", "figure_file="])
    except getopt.GetoptError:
        print (help_message)
        print(sys.argv)
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("--protein_expression_file"):
            protein_expression_file = str(arg)
            protein_expression_file = '~/input_data/GBM-CPTAC_normalized_total_protein.csv'

        if opt in ("--figure_script"):
            figure_script = str(arg)
            figure_script = '~/002f-cptac_protein_boxplots.R'
            
        if opt in ("--figure_data_file"):
            figure_data_file = str(arg)
            figure_data_file = '~/input_data/intermediate_files/cptac_protein_boxplots_data.tsv'
        if opt in ("--figure_stats_file"):
            figure_stats_file = str(arg)
            figure_stats_file = '~/input_data/intermediate_files/cptac_protein_boxplots_stats.tsv'
        if opt in ("--figure_file"):
            figure_file = str(arg)
            figure_file = '~/cptac_protein_boxplots.pdf'

    main()




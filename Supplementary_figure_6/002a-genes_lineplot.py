# create a lineplot of genes' pvalues using the two methods

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


def main():
    print('XXX-XXX.py')
    t1 = time.time()

    # load and fdr correct, then save
    pval_df = pd.read_csv(merged_pvalues_file, sep='\t')

    print(pval_df)

    print("Directional passed:", str(np.sum(pval_df['DPM'] < cutoff)), "Directional failed:", str(np.sum(pval_df['Brown'] < cutoff) - np.sum(pval_df['DPM'] < cutoff)))


    pval_df.to_csv(figure_data_file, sep='\t')

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

    cutoff = 0.05


    try:
        opts, args = getopt.getopt(sys.argv[1:],"",["merged_pvalues_file=", "figure_script=", "figure_data_file=", "figure_file="])
    except getopt.GetoptError:
        print (help_message)
        print(sys.argv)
        sys.exit(2)
    for opt, arg in opts:
        # file with p-values for each gene from each method
        if opt in ("--merged_pvalues_file"):
            merged_pvalues_file = str(arg)
            merged_pvalues_file = '~/input_data/intermediate_files/merged_pvalues.tsv'

        # figure script
        if opt in ("--figure_script"):
            figure_script = str(arg)
            figure_script = '~/002b-genes_lineplot.R'

        # figure data file
        if opt in ("--figure_data_file"):
            figure_data_file = str(arg)
            figure_data_file = '~/input_data/intermediate_files/genes_lineplot.tsv'
        # figure file
        if opt in ("--figure_file"):
            figure_file = str(arg)
            figure_file = '~/figure_genes_lineplot.pdf'

    main()




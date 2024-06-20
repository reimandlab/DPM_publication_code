# visualize a dot-plot comparing KD and OE genes P-values and FC

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


def main():
    print('XXX-XXX.py')
    t1 = time.time()

    # load each df and merge
    pval_df = pd.read_csv(p_values_file)
    fc_df = pd.read_csv(fc_file)

    # melt dfs
    pval_df = pd.melt(pval_df, id_vars='gene', var_name='condition', value_name='pvalue')
    fc_df = pd.melt(fc_df, id_vars='gene', var_name='condition', value_name='logfc')

    # merge
    merged_df = pd.merge(pval_df, fc_df, on='gene')

    # save to file
    merged_df.to_csv(figure_data_file, sep='\t', index_col=0)

    # create figure
    cline = [rscript, figure_script, '--figure_data_file', figure_data_file, '--figure_file', figure_file]

    # create comparison plots
    print(" ".join(cline))
    subprocess.run(cline)

    print(round(time.time() - t1, 2))
    print('XXX-XXX.py COMPLETE')


if __name__ == "__main__":
    # r environment
    rscript = "~/Rscript"

    try:
        opts, args = getopt.getopt(sys.argv[1:],"",["p_values_file=", "fc_file=", "figure_script=", "figure_data_file=", "figure_file="])
    except getopt.GetoptError:
        print (help_message)
        print(sys.argv)
        sys.exit(2)
    for opt, arg in opts:
        # file of pvalues from KD and OE
        if opt in ("--p_values_file"):
            p_values_file = str(arg)
            p_values_file = '~/input_data/hoxa10_pvals.csv'
        # file of the FCs of each gene in KD and OE
        if opt in ("--fc_file"):
            fc_file = str(arg)
            fc_file = '~/input_data/hoxa10_fcs.csv'

        # script to create figure
        if opt in ("--figure_script"):
            figure_script = str(arg)
            figure_script = '~/002a-dotplots_panel_c.R'

        # file to save figure data
        if opt in ("--figure_data_file"):
            figure_data_file = str(arg)
            figure_data_file = '~/input_data/dotplots_panel_c_data.tsv'
        # file to save figure
        if opt in ("--figure_file"):
            figure_file = str(arg)
            figure_file = '~/002-dotplots_panel_c.pdf'

    main()




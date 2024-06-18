# Alec Bahcheli
# create a dot-plot / line-plot comparing the p-values of DPM and Brown's methods

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
    rscript = "~/Rscript"

    cutoff = 0.05


    try:
        opts, args = getopt.getopt(sys.argv[1:],"",["merged_pvalues_file=", "figure_script=", "figure_data_file=", "figure_file="])
    except getopt.GetoptError:
        print (help_message)
        print(sys.argv)
        sys.exit(2)
    for opt, arg in opts:
        # file of merged pvalues using DPM and Brown's methods
        if opt in ("--merged_pvalues_file"):
            merged_pvalues_file = str(arg)

        # script to create figure
        if opt in ("--figure_script"):
            figure_script = str(arg)

        # file to save figure data
        if opt in ("--figure_data_file"):
            figure_data_file = str(arg)
        # file to save figure
        if opt in ("--figure_file"):
            figure_file = str(arg)

    main()




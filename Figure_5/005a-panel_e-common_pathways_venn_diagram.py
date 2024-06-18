# Alec Bahcheli
# create a venn diagram of the common pathways between the two merging methods

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


def assign_groups_and_merge(dpm_df, brown_df):
    # assign groups
    dpm_df = dpm_df.loc[:,['term_id']]
    dpm_df.loc[:,'source'] = 'dpm'

    brown_df = brown_df.loc[:,['term_id']]
    brown_df.loc[:,'source'] = 'brown'

    return pd.concat([dpm_df, brown_df])


def main():
    print('XXX-XXX.py')
    t1 = time.time()

    # load dfs
    dpm_df = pd.read_csv(dpm_enrichment_file)
    brown_df = pd.read_csv(brown_enrichment_file)

    res_df = assign_groups_and_merge(dpm_df, brown_df)

    # save to file
    res_df.to_csv(figure_data_file, sep='\t', index=False)
    

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

    # brown_enrichment_file = "/.mounts/labs/reimandlab/private/users/abahcheli/reimand_lab_ab/small_projects_reimand/results/apw2/v2/browns/enriched_pathways.csv"
    # dpm_enrichment_file = "/.mounts/labs/reimandlab/private/users/abahcheli/reimand_lab_ab/small_projects_reimand/results/apw2/v2/dpm/enriched_pathways.csv"


    # figure_script = "/.mounts/labs/reimandlab/private/users/abahcheli/reimand_lab_ab/small_projects_reimand/bin/apw2/006-common_pathways.R"
    # figure_data_file = "/.mounts/labs/reimandlab/private/users/abahcheli/reimand_lab_ab/small_projects_reimand/results/apw2/v2/common_pathways.tsv"
    # figure_file = "/.mounts/labs/reimandlab/private/users/abahcheli/reimand_lab_ab/small_projects_reimand/results/apw2/v2/common_pathways.pdf"


    try:
        opts, args = getopt.getopt(sys.argv[1:],"",["brown_enrichment_file=", "dpm_enrichment_file=", "figure_script=", "figure_data_file=", "figure_file=", "merged_fdr_file="])
    except getopt.GetoptError:
        print (help_message)
        print(sys.argv)
        sys.exit(2)
    for opt, arg in opts:
        # enriched pathways from brown's method
        if opt in ("--brown_enrichment_file"):
            brown_enrichment_file = str(arg)
        # enriched pathways from dpm's method
        if opt in ("--dpm_enrichment_file"):
            dpm_enrichment_file = str(arg)

        # figure script
        if opt in ("--figure_script"):
            figure_script = str(arg)
        # figure data file
        if opt in ("--figure_data_file"):
            figure_data_file = str(arg)
        # figure file
        if opt in ("--figure_file"):
            figure_file = str(arg)

    main()




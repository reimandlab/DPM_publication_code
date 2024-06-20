# visualize gene promoter methlyation of genes of interest in TCGA data

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

def translate_meth_df(meth, probe_df, genes_list):
    # subset to genes of interest
    probe_df = probe_df[probe_df['genesUniq'].isin(genes_list)]

    print(f'Number of probes: {len(probe_df)} / {len(genes_list)}')

    # create translation dictionary
    probe_ids = list(map(lambda x: x.split("_")[0], probe_df['probeID']))
    probe_dict = dict(zip(probe_ids, probe_df['genesUniq']))
    
    # translate prob ids
    meth = meth.rename(index = probe_dict)

    # remove non-genes from methylation data
    # mask = np.invert(meth.index.astype('str').str.contains("^cg", regex=True))
    # meth = meth.loc[mask,:]
    meth = meth.loc[np.isin(meth.index, genes_list),:]

    print(meth)
    # get average per sample
    meth = meth.groupby(level=0).apply(lambda x: x.mean())
    print(meth)

    return meth


def classify_and_melt(meth_df, mut_patients):
    # labels columns as mutant or wildtype
    meth_df.columns = meth_df.columns.str.split("-").str[:3].str.join("-")
    meth_df.columns = np.where(meth_df.columns.isin(mut_patients), "mutant", "wildtype")

    # melt
    meth_df = meth_df.reset_index().melt(id_vars = "index", var_name = "mutation_status", value_name = "value").rename(columns = {'index':'gene'})

    return meth_df


def main():
    print('XXX-XXX.py')
    t1 = time.time()

    # load methylation and probe df
    probe_df = pd.read_csv(methyl_probes_file, sep='\t')
    meth_df = pd.read_csv(methylation_file, sep='\t', index_col=0)

    # remove R132H
    meth_df = meth_df.loc[:,meth_df.columns != "TCGA-06-0221"]
    
    # translate the methylation probes to gene promoters and subset to genes of interest
    meth_df = translate_meth_df(meth_df, probe_df, genes_list)
    
    # label and melt
    meth_df = classify_and_melt(meth_df, mut_patients)

    print(meth_df)
    print(np.unique(meth_df['mutation_status'], return_counts=True))

    # save to file
    meth_df.to_csv(figure_data_file, sep='\t', index=False)

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
        opts, args = getopt.getopt(sys.argv[1:],"",["methyl_probes_file=", "methylation_file=", "figure_script=", "figure_data_file=", "figure_stats_file=", "figure_file="])
    except getopt.GetoptError:
        print (help_message)
        print(sys.argv)
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("--methyl_probes_file"):
            methyl_probes_file = str(arg)
            methyl_probes_file = '~/input_data/README_EPICv2.hg38.manifest.gencode.v41.txt'
        if opt in ("--methylation_file"):
            methylation_file = str(arg)
            methylation_file = '~/input_data/README_tcga_gbm_methylation.txt'

        if opt in ("--figure_script"):
            figure_script = str(arg)
            figure_script = '~/002d-tcga_methylation_boxplots.R'
            
        if opt in ("--figure_data_file"):
            figure_data_file = str(arg)
            figure_data_file = '~/input_data/intermediate_files/tcga_methylation_boxplots_data.tsv'
        if opt in ("--figure_stats_file"):
            figure_stats_file = str(arg)
            figure_stats_file = '~/input_data/intermediate_files/tcga_methylation_boxplots_stats.tsv'
        if opt in ("--figure_file"):
            figure_file = str(arg)
            figure_file = '~/tcga_methylation_boxplots.pdf'


    main()




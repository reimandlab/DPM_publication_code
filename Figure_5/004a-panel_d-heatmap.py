# visualize the directional agreements and disagreement between data sources in a heatmap

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

def subset_genes_pval(scores, scores_dir, merged_pval_file, fdr_cutoff = 0.001):
    res_df = []

    # load merged fdrs
    merged_fdr = pd.read_csv(merged_pval_file, sep='\t')
    merged_fdr = merged_fdr.loc[:,['pvals_merged']]

    # fdr-correct the merged p-values
    merged_fdr.loc[:,'pvals_merged'] = multitest.fdrcorrection(merged_fdr.loc[:,'pvals_merged'].to_numpy())[1]

    # subset to genes past a fdr cutoff
    new_genes = merged_fdr.index.to_numpy()[merged_fdr['pvals_merged'] < fdr_cutoff]

    # log10 for visualization
    res_df = -np.log10(scores.loc[new_genes,:])
    res_df[res_df > 10] = 10
    res_df[res_df < -10] = -10

    # which have positive fc
    mask = scores_dir.loc[new_genes,:].to_numpy() > 0
    mask_int = mask.astype('int')

    # p-values in the other direction should be negative
    mask_int[np.invert(mask)] = -1
    res_df = res_df * mask_int
    res_df = res_df.mask(res_df == -0, 0)

    return res_df


def evaluate_directional_agreement(res_df):
    # define directions
    dir1 = ['rna', 'protein']
    dir2 = ['methylation']

    # check directional agreement
    for gene in res_df.index:
        condition1 = np.logical_and((res_df.loc[gene,dir1].to_numpy() >= 0).all(), (res_df.loc[gene,dir2].to_numpy() <= 0).all())
        condition2 = np.logical_and((res_df.loc[gene,dir1].to_numpy() <= 0).all(), (res_df.loc[gene,dir2].to_numpy() >= 0).all())
        
        # print(condition1, condition2)
        if np.logical_or(condition1, condition2):
            res_df.loc[gene,'direction'] = 'darkcyan'


    res_df = res_df.fillna('darkorange')

    return res_df



def main():
    print('XXX-XXX.py')
    t1 = time.time()

    # load pathways and scores (and directions)
    scores_dir = pd.read_csv(fc_file, sep='\t', index_col=0)
    scores = pd.read_csv(pvalues_file, sep='\t', index_col=0)

    # fill NA's with p = 1
    for df in [scores_dir, scores]:
        df = df.fillna(1)

    # subset to genes of interest for visualization, add direction to P-values (positive and negative P-values)
    res_df = subset_genes_pval(scores, scores_dir, merged_pval_file)


    # annotate direction and add cancer genes
    res_df = evaluate_directional_agreement(res_df)
    cancer_genes = pd.read_csv(cgc_file, sep='\t')['Gene Symbol'].to_numpy()

    res_df.loc[:,'cancer_gene'] = 'white'
    res_df.loc[np.isin(res_df.index, cancer_genes),'cancer_gene'] = 'firebrick'

    # cancer genes should be a different color
    mask = res_df['cancer_gene'] == 'firebrick'
    res_df.loc[mask,'significant'] = res_df.index[mask]


    # save to file
    res_df.to_csv(figure_data_file, sep='\t')

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
        opts, args = getopt.getopt(sys.argv[1:],"",["pvalues_file=", "fc_file=", "cgc_file=", "merged_pval_file=", "figure_script=", "figure_data_file=", "figure_file="])
    except getopt.GetoptError:
        print (help_message)
        print(sys.argv)
        sys.exit(2)
    for opt, arg in opts:
        # p-values from each data source
        if opt in ("--pvalues_file"):
            pvalues_file = str(arg)
            pvalues_file = '~/input_data/intermediate_files/figure5_pvals.tsv'
        # FC from each data source
        if opt in ("--fc_file"):
            fc_file = str(arg)
            fc_file = '~/input_data/intermediate_files/figure5_fcs.tsv'

        # cancer gene census (CGC) file
        if opt in ("--cgc_file"):
            cgc_file = str(arg)
            cgc_file = '~/input_data/cgc_v98_29062023.tsv'

        # direcional merged p-values
        if opt in ("--merged_pval_file"):
            merged_pval_file = str(arg)
            merged_pval_file = '~/input_data/intermediate_files/figure5_merged_pvalues.tsv'

        # figure script
        if opt in ("--figure_script"):
            figure_script = str(arg)
            figure_script = '~/004b-panel_d-heatmap.R'
        # figure data file
        if opt in ("--figure_data_file"):
            figure_data_file = str(arg)
            figure_data_file = '~/input_data/intermediate_files/figure5_heatmap.tsv'
        # figure file
        if opt in ("--figure_file"):
            figure_file = str(arg)
            figure_file = '~/figure5_heatmap.pdf'

    main()




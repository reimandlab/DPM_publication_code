# perform differential gene expression analysis on the R132G mutant vs. non-mutant patients for protein-coding genes

import sys, getopt, time, os, math

try:
    sys.path.insert(1, "/".join(os.path.realpath(__file__).split("/")[:-2]))
except:
    sys.path.insert(1, "/".join(os.getcwd().split("/")[:-1]))

import pandas as pd 
import numpy as np

import scipy.stats as stats
import statsmodels.stats.multitest as multitest
from concurrent.futures import ProcessPoolExecutor


help_message = '''
Failed
'''

def protein_coding_genes_arr(protein_coding_file):
    # load and subset to protein-coding
    pc = pd.read_csv(protein_coding_file, sep='\t')
    mask = pc['locus_group'].str.contains("protein", regex=True)
    
    protein_coding_genes = list(pc['symbol'].to_numpy()[mask])
    
    # add IG genes
    protein_coding_genes.extend(["IG_C_gene", "IG_D_gene", "IG_J_gene", "IG_LV_gene", "IG_V_gene", "TR_C_gene", "TR_J_gene", "TR_V_gene", "TR_D_gene"])
    
    return np.array(protein_coding_genes).astype("<U64")


def calculate_diff_exp(args):
    gene_list, mut_df, wt_df = args 
    
    # create a list to record p-values and fc
    res_list = []

    for gene in gene_list:
        # define the values from each group
        mut_arr = mut_df.loc[gene,:].to_numpy()
        wt_arr = wt_df.loc[gene,:].to_numpy()

        # calculate fold change of means of each group, and p-value
        fc = np.nanmean(mut_arr) / np.nanmean(wt_arr)
        p = stats.mannwhitneyu(mut_arr, wt_arr)[1]

        res_list.append([gene, p, fc])

    return np.array(res_list).astype('str')


def main():
    print('XXX-XXX.py')
    t1 = time.time()

    # load expression df
    exp_df = pd.read_csv(gene_expression_file, sep='\t', index_col=0)
    exp_df.columns = list(map(lambda x: "-".join(x.split("-")[:3]), exp_df.columns))
    exp_df = exp_df.loc[np.invert(exp_df.index.duplicated()),:]

    # remove R132G
    exp_df = exp_df.loc[:,exp_df.columns != "TCGA-06-0221"]
    exp_df = exp_df.loc[:,exp_df.columns != "TCGA-12-0818"]

    # subset to protein-coding genes
    protein_coding_arr = protein_coding_genes_arr(protein_coding_file)
    exp_df = exp_df.loc[np.isin(exp_df.index, protein_coding_arr),:]
    
    # subset to mutant and wt dfs
    mask = np.isin(exp_df.columns, mut_patients)
    mut_df = exp_df.loc[:,mask]
    wt_df = exp_df.loc[:,np.invert(mask)]

    print(mut_df.shape)
    print(np.array(mut_patients)[np.isin(mut_patients, exp_df.columns.astype('str'))])

    # separate genes into groups of approximately equal sizes
    genes_arr = wt_df.index.to_numpy(dtype='<U64')
    sublist_size = math.ceil(len(genes_arr) / threads)

    separated_lists = [genes_arr[i:i+sublist_size] for i in range(0, len(genes_arr), sublist_size)]


    # genes and dataframes
    args = [(sep_list, mut_df, wt_df) for sep_list in separated_lists]

    # multi-thread calculate p-values and fc
    with ProcessPoolExecutor(max_workers=threads) as executor:
        res = np.concatenate(list(executor.map(calculate_diff_exp, args)))


    # create df and multi-test correct
    res_df = pd.DataFrame(res, columns=['gene', 'p_value', 'fold_change'])
    res_df.loc[:,'fdr'] = multitest.fdrcorrection(res_df['p_value'].astype('float'))[1]

    # save to file
    res_df.to_csv(degs_file, sep='\t', index=False)


    print(round(time.time() - t1, 2))
    print('XXX-XXX.py COMPLETE')


if __name__ == "__main__":
    # patients with mutations according to GDC
    mut_patients = ['TCGA-16-1460', 'TCGA-06-5417', 'TCGA-26-1442', 'TCGA-14-1456', 'TCGA-06-2570', 'TCGA-06-6701', 'TCGA-27-2521', 'TCGA-16-0849', 'TCGA-06-6389', 'TCGA-16-0850', 'TCGA-19-1788', 'TCGA-19-A6J5', 'TCGA-14-4157', 'TCGA-06-A7TL', 'TCGA-19-2629', 'TCGA-06-0129', 'TCGA-12-0827', 'TCGA-32-4208', 'TCGA-14-1821', 'TCGA-02-2483', 'TCGA-06-0128']

    # threads
    threads = 1

    try:
        opts, args = getopt.getopt(sys.argv[1:],"",["gene_expression_file=", "protein_coding_file=", "degs_file=", "threads="])
    except getopt.GetoptError:
        print (help_message)
        print(sys.argv)
        sys.exit(2)
    for opt, arg in opts:
        # tcga GBM expression file (counts not TPM)
        if opt in ("--gene_expression_file"):
            gene_expression_file = str(arg)
            gene_expression_file = '~/input_data/tcga_gbm_counts.tsv'
        # HGNC-annotated file to identify protein-coding genes
        if opt in ("--protein_coding_file"):
            protein_coding_file = str(arg)
            protein_coding_file = '~/input_data/hgnc_annotated.tsv'

        # path to an output file
        if opt in ("--degs_file"):
            degs_file = str(arg)
            degs_file = '~/tcga_gbm_degs.tsv'
            
        if opt in ("--threads"):
            threads = int(arg)

    main()




# Alec Bahcheli
# perform differential methylation analysis on R132G mutant vs. non-mutant patients for protein-coding genes

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


def translate_meth_df(meth, probe_df):
    # create translation dictionary
    probe_ids = list(map(lambda x: x.split("_")[0], probe_df['probeID']))
    probe_dict = dict(zip(probe_ids, probe_df['genesUniq']))

    # translate prob ids
    meth = meth.rename(index = probe_dict)

    # remove non-genes from methylation data
    mask = np.invert(meth.index.astype('str').str.contains("^cg", regex=True))
    meth = meth.loc[mask,:]

    # collect uniq genes
    uniq_genes = []

    for genes in meth.index.astype('str'):
        uniq_genes.extend(genes.split(";"))

    uniq_genes = np.unique(uniq_genes)

    return meth, uniq_genes


def protein_coding_genes_arr(protein_coding_file):
    # load and subset to protein-coding
    pc = pd.read_csv(protein_coding_file, sep='\t')
    mask = pc['locus_group'].str.contains("protein", regex=True)
    
    protein_coding_genes = list(pc['symbol'].to_numpy()[mask])
    
    # add IG genes
    protein_coding_genes.extend(["IG_C_gene", "IG_D_gene", "IG_J_gene", "IG_LV_gene", "IG_V_gene", "TR_C_gene", "TR_J_gene", "TR_V_gene", "TR_D_gene"])
    
    return np.array(protein_coding_genes).astype("<U64")


def diff_meth_genes(args):
    gene_list, mut_df, wt_df = args

    res = []

    for gene in gene_list:
        # which rows have the genes
        mask = mut_df.index.astype('str').str.contains(gene)
        
        # for each patient, calculate mean methylation
        mut_arr = mut_df.loc[mask,:].apply(np.nanmean)
        wt_arr = wt_df.loc[mask,:].apply(np.nanmean)
        
        # drop nan
        mut_arr = mut_arr[~np.isnan(mut_arr)]
        wt_arr = wt_arr[~np.isnan(wt_arr)]
        
        if len(mut_arr) == 0:
            mut_arr = np.array([0.489])
        if len(wt_arr) == 0:
            wt_arr = np.array([0.489])
        
        
        fc = np.mean(mut_arr) / np.mean(wt_arr)
        p = stats.mannwhitneyu(mut_arr, wt_arr)[1]

        res.append([gene, p, fc])

    return np.array(res).astype('str') 



def main():
    print('XXX-XXX.py')
    t1 = time.time()


    # load methylation and probe df
    probe_df = pd.read_csv(probe_file, sep='\t')
    meth_df = pd.read_csv(gene_methylation_file, sep='\t', index_col=0)


    # remove R132H
    meth_df = meth_df.loc[:,meth_df.columns != "TCGA-06-0221"]
    
    # translate the methylation probes to gene promoters and create a list of all genes with measurements 
    meth_df, uniq_genes = translate_meth_df(meth_df, probe_df)

    print(meth_df.shape)
    
    # subset to protein-coding genes
    protein_coding_arr = protein_coding_genes_arr(protein_coding_file)
    meth_df = meth_df.loc[np.isin(meth_df.index, protein_coding_arr),:]

    print(meth_df.shape)

    # subset to mutant and wt dfs
    mask = np.isin(meth_df.columns, mut_patients)
    mut_df = meth_df.loc[:,mask]
    wt_df = meth_df.loc[:,np.invert(mask)]

    print(mut_df.shape)
    print(np.array(mut_patients)[np.isin(mut_patients, meth_df.columns.astype('str'))])

    
    # separate genes into groups of approximately equal sizes
    genes_arr = np.unique(wt_df.index.to_numpy(dtype='<U64'))
    sublist_size = math.ceil(len(genes_arr) / threads)

    separated_lists = [genes_arr[i:i+sublist_size] for i in range(0, len(genes_arr), sublist_size)]


    # genes and dataframes
    args = [(sep_list, mut_df, wt_df) for sep_list in separated_lists]

    # multi-thread calculate p-values and fc
    with ProcessPoolExecutor(max_workers=threads) as executor:
        res = np.concatenate(list(executor.map(diff_meth_genes, args)))

    
    # create df, multi-test correct, and save
    res = pd.DataFrame(res, columns=['gene', 'pvalue', 'FC'])

    # filter by methylation signal
    mask = np.abs(np.log2(res['FC'].astype('float'))) < min_fc
    res.loc[mask, 'pvalue'] = 1

    # save to file
    res.to_csv(degs_file, sep='\t', index=False)

    res.loc[:,'fdr'] = multitest.fdrcorrection(res['pvalue'].astype('float'))[1]

    res.to_csv(degs_file, sep='\t', index=False)


    print(round(time.time() - t1, 2))
    print('XXX-XXX.py COMPLETE')


if __name__ == "__main__":
    # R132G mutant patients
    mut_patients = ['TCGA-16-1460', 'TCGA-06-5417', 'TCGA-26-1442', 'TCGA-14-1456', 'TCGA-06-2570', 'TCGA-06-6701', 'TCGA-27-2521', 'TCGA-16-0849', 'TCGA-06-6389', 'TCGA-16-0850', 'TCGA-19-1788', 'TCGA-19-A6J5', 'TCGA-14-4157', 'TCGA-06-A7TL', 'TCGA-19-2629', 'TCGA-06-0129', 'TCGA-12-0827', 'TCGA-32-4208', 'TCGA-14-1821', 'TCGA-02-2483', 'TCGA-06-0128', 'TCGA-12-0818', 'TCGA-06-0221']

    # minimum FC cutoff
    min_fc = 0.25

    threads = 1

    try:
        opts, args = getopt.getopt(sys.argv[1:],"",["probe_file=", "gene_methylation_file=", "protein_coding_file=", "degs_file=", "threads="])
    except getopt.GetoptError:
        print (help_message)
        print(sys.argv)
        sys.exit(2)
    for opt, arg in opts:
        # EPICv2.hg38.manifest.gencode.v41.tsv describing how the probes map to gene promoters
        # epic_v2.txt
        if opt in ("--probe_file"):
            probe_file = str(arg)
        # tcga GBM illumina 450k methylation file containing the methylation data for probes
        # tcga_methylation.txt
        if opt in ("--gene_methylation_file"):
            gene_methylation_file = str(arg)
        # protein coding genes file annotated by HGNC
        # hgnc_gene_annotations.tsv
        if opt in ("--protein_coding_file"):
            protein_coding_file = str(arg)

        # output file for differentially methylated genes
        if opt in ("--degs_file"):
            degs_file = str(arg)

        if opt in ("--threads"):
            threads = int(arg)
            
    main()



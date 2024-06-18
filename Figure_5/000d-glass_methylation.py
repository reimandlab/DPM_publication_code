# Alec Bahcheli
# differential methylation analysis for glass cohort

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

def subset_to_cohort(df, clinical_df):
    # process patient ids
    patient_ids = np.array(list(map(lambda x: "-".join(str(x).split("-")[:3]), df)))

    # subset to the correct samples: GLASS cohort only with high-grade gliomas
    glass_mask = clinical_df['project'] == 'GLSS'
    barcodes = clinical_df['barcode'][glass_mask]
    barcodes = barcodes.to_numpy(dtype='<U64')

    # find common ids
    common_ids = np.intersect1d(patient_ids, barcodes)

    # number of mutants
    n_mutant = np.sum(clinical_df.loc[common_ids,'idh_codel_subtype'].str.contains('mut'))
    print(f"Number of samples: {len(common_ids)} \nNumber of IDHmut samples: {n_mutant}")

    # separate into mutant and wt
    mutant_ids = common_ids[clinical_df.loc[common_ids,'idh_status'].str.contains('IDHmut') == True]
    mut_df = df.loc[:,np.isin(patient_ids, mutant_ids)]
    wt_df = df.loc[:,np.invert(np.isin(patient_ids, mutant_ids))]

    # subset to primary samples
    sorted = mut_df.columns[mut_df.columns.str.split("-").str[3].argsort()]
    duplicated = sorted.str.split("-").str[2].duplicated()

    mut_df = mut_df.loc[:,~duplicated]

    # repeat for wt
    sorted = wt_df.columns[wt_df.columns.str.split("-").str[3].argsort()]
    duplicated = sorted.str.split("-").str[2].duplicated()

    wt_df = wt_df.loc[:,~duplicated]


    return mut_df, wt_df


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
    gene_list, mut_df, wt_df, mut_median, wt_median = args

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
        

        if len(mut_arr) == 0 or len(wt_arr) == 0:
            print(mut_arr)
            print(wt_arr)
            pass
        
        fc = np.mean(mut_arr) / np.mean(wt_arr)
        p = stats.mannwhitneyu(mut_arr, wt_arr)[1]

        res.append([gene, p, fc])

    return np.array(res).astype('str')



def main():
    print('XXX-XXX.py')
    t1 = time.time()

    # load methylation df
    meth_df = pd.read_csv(methylation_file, sep='\t', index_col=0)
    meth_df.columns = meth_df.columns.str.replace(".", "-")

    # load methylation probes
    probe_df = pd.read_csv(methyl_probes_file, sep='\t')

    # translate the methylation probes to gene promoters and create a list of all genes with measurements 
    meth_df, uniq_genes = translate_meth_df(meth_df, probe_df)
    

    # save to file 
    meth_df.to_csv(methylation_by_gene_file, sep='\t')

    
    # subset to protein-coding genes
    protein_coding_arr = protein_coding_genes_arr()
    meth_df = meth_df.loc[np.isin(meth_df.index, protein_coding_arr),:]

    print(f"Number of methylation genes: {len(uniq_genes)} \nNumber of protein-coding genes: {meth_df.index.size}")

    # load clinical df
    clinical_df = pd.read_csv(clinical_file, sep='\t')
    clinical_df.index = clinical_df['barcode'].to_numpy()

    # subset to the correct samples: GLASS cohort only with high-grade gliomas
    mut_df, wt_df = subset_to_cohort(meth_df, clinical_df)

    # calculate median of each group
    mut_median = np.nanmedian(mut_df)
    wt_median = np.nanmedian(wt_df)

    # separate genes into groups of approximately equal sizes
    genes_arr = np.unique(wt_df.index.to_numpy(dtype='<U64'))
    sublist_size = math.ceil(len(genes_arr) / threads)

    separated_lists = [genes_arr[i:i+sublist_size] for i in range(0, len(genes_arr), sublist_size)]



    # genes and dataframes
    args = [(sep_list, mut_df, wt_df, mut_median, wt_median) for sep_list in separated_lists]

    # multi-thread calculate p-values and fc
    with ProcessPoolExecutor(max_workers=threads) as executor:
        res = np.concatenate(list(executor.map(diff_meth_genes, args)))

    
    # create df, multi-test correct, and save
    res = pd.DataFrame(res, columns=['gene', 'pvalue', 'FC'])
    res.to_csv(dgea_file, sep='\t', index=False)

    res.loc[:,'fdr'] = multitest.fdrcorrection(res['pvalue'].astype('float'))[1]

    res.to_csv(dgea_file, sep='\t', index=False)

    print(res.index.size)


    print(round(time.time() - t1, 2))
    print('XXX-XXX.py COMPLETE')


if __name__ == "__main__":
    # protein_coding_file = "/.mounts/labs/reimandlab/private/users/abahcheli/human_genome/hgnc_gene_annotations.tsv"

    # threads
    threads = 1

    try:
        opts, args = getopt.getopt(sys.argv[1:],"",["methylation_file=", "clinical_file=", "protein_coding_file=", "methyl_probes_file=", "methylation_by_gene_file=", "dgea_file=", "threads="])
    except getopt.GetoptError:
        print (help_message)
        print(sys.argv)
        sys.exit(2)
    for opt, arg in opts:
        # GLASS methylation file
        # glass_methylation.txt
        if opt in ("--methylation_file"):
            methylation_file = str(arg)
        # GLASS clinical file
        # glass_clinical.txt
        if opt in ("--clinical_file"):
            clinical_file = str(arg)
        
        # protein coding genes file HGNC-annotated
        # hgnc_gene_annotations.tsv
        if opt in ("--protein_coding_file"):
            protein_coding_file = str(arg)
        # methylation probes file associating probes with gene promoters
        # epic_v2.txt
        if opt in ("--methyl_probes_file"):
            methyl_probes_file = str(arg)

        # output file with methylation per gene
        if opt in ("--methylation_by_gene_file"):
            methylation_by_gene_file = str(arg)
        # results differential methylation file    
        if opt in ("--dgea_file"):
            dgea_file = str(arg)

        if opt in ("--threads"):
            threads = int(arg)
            
    main()




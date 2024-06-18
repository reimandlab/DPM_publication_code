# Alec Bahcheli
# differential gene expression analysis (DGEA) for glass cohort

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
    df = pd.read_csv(expression_file, sep='\t', index_col=0)
    df.columns = df.columns.str.replace(".", "-")

    # subset to protein-coding genes and primary samples
    protein_coding_arr = protein_coding_genes_arr(protein_coding_file)
    df = df.loc[np.isin(df.index.to_list(), protein_coding_arr),:]

    # load clinical df
    clinical_df = pd.read_csv(clinical_file, sep='\t')
    clinical_df.index = clinical_df['barcode'].to_numpy()

    # subset to the correct samples: GLASS cohort only with high-grade gliomas
    mut_df, wt_df = subset_to_cohort(df, clinical_df)


    # separate genes into groups of approximately equal sizes
    genes_arr = wt_df.index.to_numpy(dtype='<U64')
    sublist_size = math.ceil(len(genes_arr) / threads)

    separated_lists = [genes_arr[i:i+sublist_size] for i in range(0, len(genes_arr), sublist_size)]

    # combine and save
    mut_df.columns = list(map(lambda x: x + '-mutant', mut_df.columns))
    wt_df.columns = list(map(lambda x: x + '-wildtype', wt_df.columns))
    
    combined_df = pd.concat([mut_df, wt_df], axis=1)
    combined_df.to_csv('/.mounts/labs/reimandlab/private/users/abahcheli/reimand_lab_ab/small_projects_reimand/results/apw2/2024_02_16/_figure_data/002-gene_expression.tsv', sep='\t')

    # genes and dataframes
    args = [(sep_list, mut_df, wt_df) for sep_list in separated_lists]

    # multi-thread calculate p-values and fc
    with ProcessPoolExecutor(max_workers=threads) as executor:
        res = np.concatenate(list(executor.map(calculate_diff_exp, args)))


    # create df and multi-test correct
    res_df = pd.DataFrame(res, columns=['gene', 'pvalue', 'FC'])
    res_df.loc[:,'fdr'] = multitest.fdrcorrection(res_df['pvalue'].astype('float'))[1]

    # save to file
    res_df.to_csv(dgea_file, sep='\t', index=False)


    print(round(time.time() - t1, 2))
    print('XXX-XXX.py COMPLETE')


if __name__ == "__main__":
    threads = 1

    try:
        opts, args = getopt.getopt(sys.argv[1:],"",["expression_file=", "clinical_file=", "protein_coding_file=", "dgea_file=", "threads="])
    except getopt.GetoptError:
        print (help_message)
        print(sys.argv)
        sys.exit(2)
    for opt, arg in opts:
        # glass gene expression file (TPM)
        # glass_expression.txt
        if opt in ("--expression_file"):
            expression_file = str(arg)
        # glass clinical file with subtypes / mutations
        # glass_clinical.txt
        if opt in ("--clinical_file"):
            clinical_file = str(arg)

        # HGNC-annotated file describing protein-coding genes
        # hgnc_gene_annotations.tsv
        if opt in ("--protein_coding_file"):
            protein_coding_file = str(arg)

        # results DGEA file
        if opt in ("--dgea_file"):
            dgea_file = str(arg)
        if opt in ("--threads"):
            threads = int(arg)

    main()




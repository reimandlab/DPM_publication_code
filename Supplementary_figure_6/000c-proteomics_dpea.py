# differential protein expression analysis

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


def subset_to_cohort_proteomics(df, clinical_df):
    # samples
    samples = df.columns.str.strip("' ")

    # classify samples into mut and wt
    mut_samples = clinical_df['Sample_name'].to_numpy(dtype='<U64')[clinical_df['IDH1_mut'].str.contains('R132') == True]

    # subset to mut and wt
    mut_df = df.loc[:,np.isin(samples, mut_samples)]
    wt_df = df.loc[:,~np.isin(samples, mut_samples)]

    return mut_df, wt_df

def calculate_diff_exp(args):
    gene_list, mut_df, wt_df = args 
    
    # create a list to record p-values and fc
    res_list = []

    for gene in gene_list:
        # define the values from each group
        mut_arr = mut_df.loc[gene,:].to_numpy()
        wt_arr = wt_df.loc[gene,:].to_numpy()

        # calculate fold change of means of each group, and p-value
        fc = (1000 + np.nanmean(mut_arr)) / (1000 + np.nanmean(wt_arr))
        p = stats.mannwhitneyu(mut_arr, wt_arr)[1]

        res_list.append([gene, p, fc])

    return np.array(res_list).astype('str')


def main():
    print('XXX-XXX.py')
    t1 = time.time()

    # load expression df
    df = pd.read_csv(expression_file, sep='\t')
    df.index = df['Symbol'].str.split(";").str[-1].str.strip("'")
    df = df.iloc[:,4:]
    df = df.loc[~df.index.duplicated(),:]

    # load clinical df
    clinical_df = pd.read_csv(clinical_file, sep='\t')

    # separate into mut and wt
    mut_df, wt_df = subset_to_cohort_proteomics(df, clinical_df)


    # combine and save
    mut_df.columns = list(map(lambda x: x + '-mutant', mut_df.columns))
    wt_df.columns = list(map(lambda x: x + '-wildtype', wt_df.columns))

    print(mut_df.head())
    print(f"Mutant: {mut_df.shape}\nWT: {wt_df.shape}")
    
    combined_df = pd.concat([mut_df, wt_df], axis=1)
    combined_df.to_csv(protein_expression_file, sep='\t')

    
    
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
    res_df.to_csv(dgea_file, sep='\t', index=False)


    print(round(time.time() - t1, 2))
    print('XXX-XXX.py COMPLETE')


if __name__ == "__main__":
    # threads
    threads = 1

    try:
        opts, args = getopt.getopt(sys.argv[1:],"",["expression_file=", "clinical_file=", "protein_expression_file=", "dgea_file=", "threads="])
    except getopt.GetoptError:
        print (help_message)
        print(sys.argv)
        sys.exit(2)
    for opt, arg in opts:
        # raw expression file
        if opt in ("--expression_file"):
            expression_file = str(arg)
            expression_file = '/input_data/validation_proteomics_expression.tsv'
        # clinical file with subtypes / mutations
        if opt in ("--clinical_file"):
            clinical_file = str(arg)
            clinical_file = '/input_data/validation_proteomics_clinical.tsv'

        # output processed protein expression file
        if opt in ("--protein_expression_file"):
            protein_expression_file = str(arg)
            protein_expression_file = '/input_data/intermediate_files/validation_proteomics_expression.tsv'

        # results protein expression file
        if opt in ("--dgea_file"):
            dgea_file = str(arg)
            dgea_file = '/input_data/intermediate_files/validation_proteomics_dgea.tsv'

        if opt in ("--threads"):
            threads = int(arg)

    main()




# create dotplots of the contribution of each gene to the pathway

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

def load_gmt_file(file_path):
    gmt_dict = {}
    with open(file_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            key = parts[1]
            values = parts[2:]
            gmt_dict[key] = values
    return gmt_dict


def create_dfs(gmt_data, df):
    # Create list of dataframes
    dfs = []

    # Iterate through gmt_data and subset df
    for key, genes in gmt_data.items():
        subset_df = df[df['gene'].isin(genes)]
        # If subset_df is not empty, add gmt_key to subset_df and append to dfs
        if not subset_df.empty:
            subset_df['gmt_key'] = key
            dfs.append(subset_df)

    return pd.concat(dfs)



def melt_df(combined_df):
    # copy combined_df
    combined_df2 = combined_df.copy()

    # subset to only p-values
    cols_oi = ['gene', 'gmt_key', 'adjusted_p_val']
    cols_oi.extend(combined_df2.columns[combined_df2.columns.str.contains('Pvalue')])

    # subset to only cols_oi
    combined_df2 = combined_df2[cols_oi]

    # melt the DataFrame
    combined_df2 = pd.melt(combined_df2, id_vars=['gene', 'gmt_key', 'adjusted_p_val'], var_name='sample', value_name='pvalue').copy()


    # subset to only fc
    cols_oi = ['gene', 'gmt_key', 'adjusted_p_val']
    cols_oi.extend(combined_df.columns[combined_df.columns.str.contains('FC')])

    # subset to only cols_oi
    combined_df = combined_df[cols_oi]

    # melt the DataFrame
    combined_df = pd.melt(combined_df, id_vars=['gene', 'gmt_key', 'adjusted_p_val'], var_name='sample', value_name='log2fc').copy()


    # add pvalue to combined_df
    combined_df['pvalue'] = combined_df2['pvalue']

    combined_df['sample'] = combined_df['sample'].str.split("_").str[0]


    # add direction to fc
    combined_df['fc_dir'] = np.where(combined_df['log2fc'] > 0, 'up', 'down')

    # subset to genes that pass p-value in at least one condition
    sig_genes = np.unique(combined_df['gene'][combined_df['pvalue'] < 0.05])
    combined_df = combined_df[np.isin(combined_df['gene'], sig_genes)]

    return combined_df


def correct_protein_fc(df):
    # set protein fc to 5 and -5 
    up_mask = np.logical_and(df['sample'] == 'protein', df['log2fc'] > 0)
    down_mask = np.logical_and(df['sample'] == 'protein', df['log2fc'] < 0)

    df.loc[up_mask, 'log2fc'] = 5
    df.loc[down_mask, 'log2fc'] = -5

    return df



def main():
    print('XXX-XXX.py')
    t1 = time.time()
    
    # load gmt and data files
    gmt_data = load_gmt_file(gmt_file)
    pval_df = pd.read_csv(pval_file, sep='\t', index_col=0)
    fc_df = pd.read_csv(fc_file, sep='\t', index_col=0)

    # add descriptions to columns
    pval_df.columns = list(map(lambda x: x + '_Pvalue', pval_df.columns))
    fc_df.columns = list(map(lambda x: x + '_FC', fc_df.columns))
    
    # merge dfs
    combined_df = pd.concat([pval_df, fc_df], axis=1)
    combined_df['gene'] = combined_df.index.to_numpy()

    # Create list of dataframes
    combined_df = create_dfs(gmt_data, combined_df)

    # load pathway p-values
    adjusted_p_val_df = pd.read_csv(pathway_pval_file, sep='\t').set_index('term_name')

    # Merge the adjusted_p_val to each DataFrame in dfs
    combined_df['adjusted_p_val'] = combined_df['gmt_key'].map(adjusted_p_val_df['adjusted_p_val'])

    

    # melt the DataFrame
    melted_df = melt_df(combined_df)

    # correct protein fc
    melted_df = correct_protein_fc(melted_df)

    # save to file
    melted_df.to_csv(figure_data_file, sep='\t', index=False)

    # create figure
    cline = [rscript, figure_script, '--figure_data_file', figure_data_file, '--figure_file', figure_file]

    # create comparison plots
    print(" ".join(cline))
    subprocess.run(cline)


    print(round(time.time() - t1, 2))
    print('XXX-XXX.py COMPLETE')


if __name__ == "__main__":
    rscript = "Rscript"

    try:
        opts, args = getopt.getopt(sys.argv[1:],"",["gmt_file=", "pathway_pval_file=", "pval_file=", "fc_file=", "figure_script=", "figure_data_file=", "figure_file="])
    except getopt.GetoptError:
        print (help_message)
        print(sys.argv)
        sys.exit(2)
    for opt, arg in opts:
        # original gene set matrix file
        if opt in ("--gmt_file"):
            gmt_file = str(arg)
            gmt_file = '/input_data/hsapiens.GO_BP_REACTOME.name.gmt'
        # enriched pathway file
        if opt in ("--pathway_pval_file"):
            pathway_pval_file = str(arg)
            pathway_pval_file = '/input_data/enriched_pathways.tsv'
        # file of p-values
        if opt in ("--pval_file"):
            pval_file = str(arg)
            pval_file = '/input_data/hoxa10_fdr.csv'
        # file of fold changes
        if opt in ("--fc_file"):
            fc_file = str(arg)
            fc_file = '/input_data/hoxa10_fc.csv'

        # figure script
        if opt in ("--figure_script"):
            figure_script = str(arg)
            figure_script = '/001b-gene_pathway_dotplots.R'
        # figure data file
        if opt in ("--figure_data_file"):
            figure_data_file = str(arg)
            figure_data_file = '/hoxa10_pathway_contributions.tsv'
        # figure file
        if opt in ("--figure_file"):
            figure_file = str(arg)
            figure_file = '/hoxa10_pathway_contributions.pdf'
            
    main()





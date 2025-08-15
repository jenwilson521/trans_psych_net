# written to merge and streamline the zebrafish code for
# transPsych analysis, 8-14-25 JLW

import pandas as pd
import numpy as np
import os

fpath = '~/Documents/UCLA/BIG_summer_2022/Source_Data/Figure_2b.csv'

#PRE:: df is the drugbank_vocabulary csv dataframe on page DrugsToTargets
#POST::Creates the dictionary for fast time for find_drugs function
def create_dict_drug_bank(df):
    drug_bank = {}
    for i, row in df.iterrows():
        drug_bank[row["Common name"].capitalize()] = row["DrugBank ID"]
        syn = row["Synonyms"]
        if(type(syn) != float):
          arr = syn.split(" | ")
          for drug in arr:
            drug_bank[drug.capitalize()] = row["DrugBank ID"]
    return drug_bank

#PRE:: drug_bank is the output of create_dict_drug_bank and data is the zebrafish data
#POST::Returns the drugs that can be identified and the drugs that cannot be
def find_drugs(drug_bank, data):
    in_drug_bank = []
    not_drug_bank = []
    for i, row in data.iterrows():
        #.capitalize turns PROGESTERONE -> Progesterone to fit the drugbank
        if row["chemName"].capitalize() in drug_bank:
            in_drug_bank.append([drug_bank[row["chemName"].capitalize()], row["1-corr"]])
        else:
            not_drug_bank.append(row["chemName"])
    return [in_drug_bank, not_drug_bank]

#PRE:: file in folder with sheet
#POST::Returns dataframe of the specific sheet
def read_drug2targets():
    return pd.read_excel('../data/Drugbank050120.xlsx', sheet_name='DrugsToTargets')

#PRE::
#POST::Creates a set of drugbank_id
def get_dbids_with_prot_targets():
    df = read_drug2targets()
    drugs = set()
    for i, row in df.iterrows():
        drugs.add(row['drugbank_id'])
    return drugs

def compare_dbids_w_targs(dbids_w_targs, list_drugs):
  drugs_in = []
  drugs_not = []
  for i, row in list_drugs.iterrows():
    if row["DBID"] in dbids_w_targs:
      drugs_in.append([row["DBID"], row["1-corr"]])
    else:
      drugs_not.append(row["DBID"])
  return [drugs_in, drugs_not]

print('reading drugbank vocabulary and raw data')
drug_bank = pd.read_csv ('../data/drugbank_vocabulary.csv',encoding= 'unicode_escape')
figure_2b = pd.read_csv(fpath)

# generate a dictionary of drug names and synonyms -> DBID
drug_bank_dict = create_dict_drug_bank(drug_bank)

# loop through all chicmal names and keep those which map to DBID
compare_data = find_drugs(drug_bank_dict, figure_2b)

# keep those with DBID, assign columns
print('filter chemical names for DBIDs')
df_drugs = pd.DataFrame(compare_data[0])
df_drugs.columns = ['DBID', '1-corr']
df_drugs_dropped = df_drugs.drop_duplicates(subset=['DBID']) # keeps the max phenoscore

# for those with DrugBank IDs, further restrict to those with protein targets
dbids_w_targs = get_dbids_with_prot_targets()
map_to_PathFX = compare_dbids_w_targs(dbids_w_targs, df_drugs_dropped)

# for those mapped, save their max phenoscores
print('keep drugs with protein-binding targets and save')
df_drugs_mapped = pd.DataFrame(map_to_PathFX[0])
df_drugs_mapped.columns = ['DBID', '1-corr']
df_drugs_mapped.to_csv('data_mapped.csv', index=False)

#PRE:: df is from the _merged_neighborhood_.txt
#POST::Returns Set of genes connected to the drug
def extract_genes2(df):
  genes = set()
  for i, row in df.iterrows():
    if row[0] != "" and type(row[0])!=type(1.0):
      genes.add(row[0])
    if row[1] != "" and type(row[1])!=type(1.0):
      genes.add(row[1])
  return genes

#PRE::
#POST::Iterate through the files to put in genes and value dictionary
def iter_thr_files():
  data_mapped = pd.read_csv('data_mapped.csv')
  dict_to_genes = {}
  dict_to_values = {}
  skipped_drugs = []
  #for root, dirs, files in os.walk("../results/alldrugbank"):
  for i, row in data_mapped.iterrows():
    path = "../results/alldrugbank/" + row["DBID"] + "/" + row["DBID"] + "_merged_neighborhood_.txt"
    if os.stat(path).st_size == 0:
      skipped_drugs.append(row["DBID"])
    else:
      df_temp = pd.read_csv(path,sep='\t', header=None)
      genes = extract_genes2(df_temp)
      dict_to_genes[row["DBID"]] = genes
      dict_to_values[row["DBID"]] = row["1-corr"]
  return [dict_to_genes,dict_to_values,skipped_drugs]

def get_df_drugs():
  ana_drugs = iter_thr_files()
  drugs_to_netGenes = ana_drugs[0]# a dictionary where keys are DBIDs and values are a set of network genes
  drugs_to_values = ana_drugs[1]# a dictionary where keys are DBIDs and values are either from zebra fish screen or from clinical data
  all_rows = []
  for (dbid, effect_size) in drugs_to_values.items():
    net_genes = drugs_to_netGenes[dbid]
    row_data = {'name':dbid, 'effect':effect_size}
    for g in net_genes:
      row_data[g]=1
    all_rows.append(row_data)
  df = pd.DataFrame(all_rows)
  df = df.fillna(0)   
  if len(ana_drugs[2]) == 0:
    print("All Drugs Found a Home")
  else:
    for i in range(len(ana_drugs[2])):
      print(ana_drugs[2][i])
  df.to_csv('df_zebrafish.csv', index=False)
  return df

print('create a save a dataframe for DBIDs with PathFX networks')
get_df_drugs()

### Testing
#ots = pd.read_csv('df_zebrafish.csv')  #output from this script
#oig = pd.read_csv('/Users/jenniferwilson/Documents/trans_psych_net/code/df_zebrafish.csv') #output saved in GitHub and used for everything else

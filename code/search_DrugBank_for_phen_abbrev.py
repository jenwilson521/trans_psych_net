# written to search DrugBank results
# for a list of CUI terms
# revised 7-14-23 JLW to parse essential results
# written 2-10-21 JLW

import pickle,os,csv
from collections import defaultdict
import pandas as pd
import argparse

# first gather relevant phenotypes
d2n = pickle.load(open('../rscs/Pfx050120_dbid2name.pkl','rb'))
n2db = dict([(v.lower(),k) for (k,v) in d2n.items()])
# gather arguments
parser = argparse.ArgumentParser(description='Parameters for the analysis')
parser.add_argument('-aname', action='store', dest='aname', help='Analysis Name, no spaces')
parser.add_argument('-rdir', action='store', dest='rdir', help='The path of the directory where the cui list is pickled')
parser.add_argument('-cuipkl', action='store', dest='cpkl', help='The pickled file name containing the list of CUI terms')
args = parser.parse_args()

rdir = args.rdir
keep_cuis = pickle.load(open(os.path.join(rdir,args.cpkl),'rb'))
aname = args.aname
out_excel_name = os.path.join(rdir,"all_drug_summary_"+aname+".xlsx") 
writer = pd.ExcelWriter(out_excel_name, engine='xlsxwriter')

# gather all association table files
rdir = '../results/alldrugbank/' # updated file path on macbook pro
allf = [(ssd,sd,flist) for (ssd,sd,flist) in os.walk(rdir)]
asf_dic = dict([(os.path.split(ssd)[-1],os.path.join(ssd,f)) for (ssd,sd,flist) in allf for f in flist if '_merged_neighborhood__assoc_table_.txt' in f])

print('Reading all Drug Bank')
res_df = pd.DataFrame()
for (dbid,dfile) in asf_dic.items():
	dname = d2n[dbid]
	df = pd.read_table(dfile)
	if not df.empty:
		df['DrugName'] = dname.lower()
		df["DrugBankID"] = dbid
		df_short = df[df['cui'].isin(keep_cuis)]
		if not df_short.empty:
#			print(dname)
			res_df = pd.concat([res_df,df_short])
	 
# add drug targets as FYI
print('Adding Drug Target and Description Info')
dint = pickle.load(open('../rscs/Pfx050120_dint.pkl','rb'))
dint_str = dict([(d,','.join(glist)) for (d,glist) in dint.items()])
res_df['DrugTargetProteins'] = res_df['DrugBankID'].map(dint_str)
res_df.drop(['rank','assoc in neigh', 'assoc in intom','probability','Unnamed: 8'], axis=1)
ord_cols = ['DrugName',"DrugBankID",'DrugTargetProteins','phenotype','Benjamini-Hochberg', 'cui','genes']
res_df = res_df = res_df[ord_cols]
res_df = res_df.rename(columns={"genes": "NetworkGenes",})


# add some more info from DrugBank as fyi
dbank = pd.read_excel('../data/Drugbank050120.xlsx')
dbank = dbank.rename(columns={'drugbank_id':'DrugBankID'})

res_df = res_df.merge(dbank[['DrugBankID','type','description']],on='DrugBankID')
print(res_df.head)


# look at common network mechanisms
net_genes = res_df['NetworkGenes'].to_list()
gene_counts = defaultdict(int)
for gentry in net_genes:
	glist = gentry.split(',')
	for g in glist:
		gene_counts[g]+=1

# keep the top 20 genes associated wih a drug
# pickle.dump(gene_counts,open(os.path.join(rdir,'top_net_genes_all.pkl'),'wb'))
top_mechs = sorted([x for x in gene_counts.items()],key=lambda x:x[1],reverse=True)[0:20]
top_genes = [x[0] for x in top_mechs]
keep_gene_strs = [gentry for gentry in net_genes for tg in top_genes if tg in gentry]
count_top_mechs = defaultdict(int)
for gentry in list(set(net_genes)):
	for tg in top_genes:
		if tg in gentry:
			count_top_mechs[gentry]+=1
top_mech_df = res_df[res_df['NetworkGenes'].isin(keep_gene_strs)]
top_mech_df['TopMechCount'] = top_mech_df['NetworkGenes'].map(count_top_mechs)
top_mech_df.to_excel(writer,index=False,sheet_name="TopPredicted")

top_mech_sum = pd.DataFrame(top_mechs,columns=['GeneName','Occured in X drug networks'])
top_mech_sum.to_excel(writer,index=False,sheet_name="SummOfMostCommonGenes")

res_df.to_excel(writer,index=False,sheet_name="allResults")

writer.save()
print('Results saved in '+out_excel_name)



# writen 12-18-25 JLW
# to look at coefficient values with shuffled scores

import os
import pandas as pd
import matplotlib.pyplot as plt
import statistics as stat
import seaborn as sns
import numpy as np
from sklearn.linear_model import LinearRegression, Lasso, Ridge
from collections import defaultdict

# load all network proteins
def get_net_nodes(f):
        d = [l.strip().split() for l in open(f,'r').readlines()]
        prot_a = [x[0] for x in d]
        prot_b = [x[1] for x in d]
        return set(prot_a+prot_b)

rdir = '/Users/jenniferwilson/Documents/PathFXv2/results/alldrugbank/' # updated file path on macbook pro, new user will need to update this to their file path
allf = [(ssd,sd,flist) for (ssd,sd,flist) in os.walk(rdir)]
mnh_dic = dict([(os.path.split(ssd)[-1],os.path.join(ssd,f)) for (ssd,sd,flist) in allf for f in flist if 'merged_neighborhood_.txt' in f])
asf_dic = dict([(os.path.split(ssd)[-1],os.path.join(ssd,f)) for (ssd,sd,flist) in allf for f in flist if '_merged_neighborhood__assoc_table_.txt' in f])

db_to_prots = defaultdict(set)
for (db,dbfile) in mnh_dic.items():
        net_prot = get_net_nodes(dbfile)
        db_to_prots[db] = net_prot

db_to_asf_prot = defaultdict(set)
for (db,dbfile) in asf_dic.items():
	adf = pd.read_csv(dbfile,delimiter='\t')
	for gs in adf.genes:
		split_gs = gs.split(',')
		for gname in split_gs:
			db_to_asf_prot[db].add(gname)

# load all DrugBank IDs
xlfile = '../supplemental_files/SF2_ridgeregcoeffs_GoEnrich.xlsx'
xlfile = '/Users/jenniferwilson/Documents/trans_psych_net/supplemental_files/SF2_ridgeregcoeffs_GoEnrich.xlsx'
snames = ['InSilicoDrugs','PhenScoreDrugs','EfficacyDrugs']

###########################	
# repeat in silico analysis
##########################
df = pd.read_excel(xlfile,sheet_name='InSilicoDrugs')
drugs_to_values = dict(zip(df.DBID,df.Effect))
all_rows = []
for (dbid, effect_size) in drugs_to_values.items():
    net_genes = db_to_prots[dbid]
    row_data = {'name':dbid, 'effect':effect_size}
    for g in net_genes:
        row_data[g]=1
    all_rows.append(row_data)

df = pd.DataFrame(all_rows)
df = df.fillna(0)

# analysis with real data
dbids = df['name'].to_list()
df = df.drop(['name'],axis=1)
y = df.loc[:,'effect']
X = df.drop(['effect'],axis=1)
gnames = list(X.columns)
alpha = 1.0
model = Ridge(alpha=alpha, random_state=0)
n = model.fit(X,y)
coefficients = model.coef_
reg_coef_real = dict(zip(gnames, coefficients))

# now shuffle y 100 times with random seed
reg_coeff_suf = defaultdict(list)
y_a = np.array(y)
for i in range(100):
	np.random.seed(i)
	np.random.shuffle(y_a)
	n = model.fit(X,y_a)
	shuf_coefficients = model.coef_
	for (g,s) in zip(gnames,shuf_coefficients):
		reg_coeff_suf[g].append(s)

# write new data to output
all_rows = []
for (g,s) in reg_coef_real.items():
	shuf_scores = reg_coeff_suf[g]
	ms = np.mean(shuf_scores)
	sh_med = np.median(shuf_scores)
	sh_min = min(shuf_scores)
	sh_max = max(shuf_scores)
	row_data = {'Gene Names':g,'Regression Coefficients':f"{s:.3f}",'Shuffled Mean':f"{ms:.3f}",'Shuffled Median':f"{sh_med:.3f}",'Shuffled Min':f"{sh_min:.3f}",'Shuffled Max':f"{sh_max:.3f}"}
	all_rows.append(row_data)

df_new = pd.DataFrame(all_rows)
df_new.to_csv("/Users/jenniferwilson/Documents/trans_psych_net/code/insilico_phagocyt_new.csv",index=False)
 
#############################
# recreate zebrafish analysis
############################
zpath = '/Users/jenniferwilson/Documents/trans_psych_net/code/df_zebrafish.csv'
zdf = pd.read_csv(zpath)
zdbids = zdf['name'].to_list()
zdf = zdf.drop(['name'],axis=1)
zy = zdf.loc[:,'effect']
zX = zdf.drop(['effect'],axis=1)
zgnames = list(zX.columns)

alpha = 1.0
zmodel = Ridge(alpha=alpha, random_state=0)
n = zmodel.fit(zX,zy)
zcoefficients = zmodel.coef_
zreg_coef_real = dict(zip(zgnames, zcoefficients))
# now shuffle y 100 times with random seed
zreg_coeff_suf = defaultdict(list)
zy_a = np.array(zy)
for i in range(100):
        np.random.seed(i)
        np.random.shuffle(zy_a)
        n = zmodel.fit(zX,zy_a)
        zshuf_coefficients = zmodel.coef_
        for (g,s) in zip(zgnames,zshuf_coefficients):
                zreg_coeff_suf[g].append(s)

zall_rows = []
for (g,s) in zreg_coef_real.items():
        shuf_scores = zreg_coeff_suf[g]
        ms = np.mean(shuf_scores)
        sh_med = np.median(shuf_scores)
        sh_min = min(shuf_scores)
        sh_max = max(shuf_scores)
        row_data = {'Gene Names':g,'Regression Coefficients':f"{s:.3f}",'Shuffled Mean':f"{ms:.3f}",'Shuffled Median':f"{sh_med:.3f}",'Shuffled Min':f"{sh_min:.3f}",'Shuffled Max':f"{sh_max:.3f}"}
        zall_rows.append(row_data)

df_new = pd.DataFrame(zall_rows)
df_new.to_csv("/Users/jenniferwilson/Documents/trans_psych_net/code/zebrafish_new.csv",index=False)


##############################
# recreate efficacy analysis
#############################
# use association table genes only

df = pd.read_excel(xlfile,sheet_name='EfficacyDrugs')
drugs_to_values = dict(zip(df.DBID,df['efficacy.or'].to_list()))
all_rows = []
for (dbid, effect_size) in drugs_to_values.items():
    net_genes = db_to_asf_prot[dbid]
    row_data = {'name':dbid, 'effect':effect_size}
    for g in net_genes:
        row_data[g]=1
    all_rows.append(row_data)

efficacy = pd.DataFrame(all_rows)
efficacy = efficacy.fillna(0)

# analysis with real data
x_eff = efficacy[efficacy.columns[2:326]]
y_eff = efficacy['effect']

# linear regression model
ridge = Ridge(alpha = 0.1)
ridge_mod = ridge.fit(x_eff, y_eff)
eff_coef_ridge = ridge_mod.coef_
reg_coef_real = dict(zip(efficacy.columns[2:326],eff_coef_ridge))

# now shuffle y 100 times with random seed
reg_coeff_suf = defaultdict(list)
y_a = np.array(y_eff)
for i in range(100):
        np.random.seed(i)
        np.random.shuffle(y_a)
        n = model.fit(x_eff,y_a)
        shuf_coefficients = model.coef_
        for (g,s) in zip(efficacy.columns[2:326],shuf_coefficients):
                reg_coeff_suf[g].append(s)

# write new data to output
all_rows = []
for (g,s) in reg_coef_real.items():
	shuf_scores = reg_coeff_suf[g]
	ms = np.mean(shuf_scores)
	sh_med = np.median(shuf_scores)
	sh_min = min(shuf_scores)
	sh_max = max(shuf_scores)
	row_data = {'Gene Names':g,'Regression Coefficients':f"{s:.3f}",'Shuffled Mean':f"{ms:.3f}",'Shuffled Median':f"{sh_med:.3f}",'Shuffled Min':f"{sh_min:.3f}",'Shuffled Max':f"{sh_max:.3f}"}
	all_rows.append(row_data)

df_new = pd.DataFrame(all_rows)
df_new.to_csv("/Users/jenniferwilson/Documents/trans_psych_net/code/antidepressants_new.csv",index=False)


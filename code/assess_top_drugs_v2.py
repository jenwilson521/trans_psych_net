# written 12-11-25 JLW to look at top
# drugs, their mechanisms, and patterns
# of "directionality"

import pickle,os,csv
import pandas as pd
from collections import defaultdict
from openpyxl import load_workbook

# look at favorable network proteins
cdir = '.'
f = os.path.join(cdir,'assess_top_drugs.xlsx')
df = pd.read_excel(f,sheet_name='fav_net_prots')
prot_to_avg_coeff = dict([(p,a) for (p,a) in zip(df.iloc[:,0].values,df.iloc[:,4].values)])
high_avg_prot = set(prot_to_avg_coeff.keys())

# load phen scores and clinical Odds-ratios
sdir = '../supplemental_files/'
sf2 = os.path.join(sdir,'SF2_ridgeregcoeffs_GoEnrich.xlsx')
phen_scores = pd.read_excel(sf2,sheet_name='PhenScoreDrugs')
odds_ratios = pd.read_excel(sf2,sheet_name='EfficacyDrugs')

# get all network proteins
def get_net_nodes(f):
        d = [l.strip().split() for l in open(f,'r').readlines()]
        prot_a = [x[0] for x in d]
        prot_b = [x[1] for x in d]
        return set(prot_a+prot_b)

rdir = '../PathFXv2/results/alldrugbank/' # updated file path on macbook pro, new user will need to update this to their file path
allf = [(ssd,sd,flist) for (ssd,sd,flist) in os.walk(rdir)]
mnh_dic = dict([(os.path.split(ssd)[-1],os.path.join(ssd,f)) for (ssd,sd,flist) in allf for f in flist if 'merged_neighborhood_.txt' in f])

db_to_prots = defaultdict(set)
for (db,dbfile) in mnh_dic.items():
	net_prot = get_net_nodes(dbfile)
	db_to_prots[db] = net_prot

# exploratory - first make a table with all ORs and avg coeffs
all_rows = []
for (dname,dbid,effor) in zip(odds_ratios.DrugName,odds_ratios.DBID,odds_ratios['efficacy.or'].to_list()):
	net_prot = db_to_prots[dbid]
	net_with_score =net_prot.intersection(high_avg_prot) 
	with_score = []
	if len(net_with_score)> 0:
		for np in net_with_score:
			with_score.append(np + ':' + str(prot_to_avg_coeff[np]))
		with_score_str = ' | '.join(sorted(with_score,reverse=True))
	else:
		with_score_str = ''
	row_data = {'DrugName':dname,'DBID':dbid,'efficacy.or':effor,'NetProtswithCoeffs':with_score_str}
	all_rows.append(row_data)

eff_with_scores = pd.DataFrame(all_rows)
eff_with_scores = eff_with_scores.sort_values(by='efficacy.or', ascending=False)
eff_with_scores.to_csv(os.path.join(cdir,'eff_ORs_with_net_coeffs.csv'),index=False)

# exploratory, next make a table with Phenoscores and avg coeffs
all_rows = []
for (dname,dbid,pscore) in zip(phen_scores.Drugname,phen_scores.DBID,phen_scores.PhenScore):
	net_prot = db_to_prots[dbid]
	net_with_score =net_prot.intersection(high_avg_prot) 
	with_score = []
	if len(net_with_score)> 0:
		for np in net_with_score:
			with_score.append(np + ':' + str(prot_to_avg_coeff[np]))
		with_score_str = ' | '.join(sorted(with_score,reverse=True))
	else:
		with_score_str = ''
	row_data = {'DrugName':dname,'DBID':dbid,'PhenScore':pscore,'NetProtswithCoeffs':with_score_str}
	all_rows.append(row_data)

phen_with_scores = pd.DataFrame(all_rows)
phen_with_scores = phen_with_scores.sort_values(by='PhenScore',ascending=False)
phen_with_scores.to_csv(os.path.join(cdir,'phenScores_with_net_coeffs.csv'),index=False)
	
	

# exploratory for sertraline
snet_file = "DB01104_merged_neighborhood__withDrugTargsAndPhens.txt"
snet_lines = [l.strip().split('\t') for l in open(snet_file,'r').readlines()]

drug_targets = [l[1] for l in snet_lines if l[0]=='DB01104']
complex_pro  = [l[1] for l in snet_lines if l[0] == 'complex:lAo9+jzjVOaRJSWZQvoh0G1mUDc'] + [l[0] for l in snet_lines if l[1]=='complex:lAo9+jzjVOaRJSWZQvoh0G1mUDc'] #DERL1 is connected to complex "complex:lAo9+jzjVOaRJSWZQvoh0G1mUDc", len = 37

set(drug_targets).intersection(complex_pro) # {'ALB'}

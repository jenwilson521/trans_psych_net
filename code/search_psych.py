# looks at all drugs associated with a search
# phenotype, creates some summary tables
# update for recreating psychiatry analysis, 1-10-25 JLW

import pickle,os,csv
from collections import defaultdict
import pandas as pd
import argparse

# pathway phenotypes that had sufficient connections to drug
# targets and were retained in this analysis

('Unipolar Depression', 'C0041696')
('Paranoid Schizophrenia', 'C0036349')
('Major depressive disorder', 'C1269683')
('Bipolar Disorder', 'C0005586')
('Schizophrenia', 'C0036341')

keep_cuis = ['C0036341', 'C0041696', 'C1269683', 'C0005586', 'C0036349']

aname = "psych_networks"
rdir = os.path.join("../results/",aname)
if not os.path.exists(rdir):
	os.makedirs(rdir)
cui_file =os.path.join(rdir,aname+"_keep_cuis.pkl") 
pickle.dump(keep_cuis,open(cui_file,'wb'))

cmd = "python search_DrugBank_for_phen_abbrev.pyy -aname %s -rdir %s -cuipkl %s"%(aname,rdir,aname+"_keep_cuis.pkl")
print(cmd)



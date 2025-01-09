# trans_psych_net

_code for recreating results, most will require updated file paths_
##1. Logistic and ridge regression for in silico screen, heatmap for in silico screen: "in_silico_logreg.ipynb", "in_silico_ridgereg.ipynb", "Heatmap_Insilico.ipynb"

Requirements:
- results directory "alldrugbank" containing PathFXv2 networks for all drugs in DrugBank. This can be created by running the script, "run_PathFX_allDrugBank.py". You will need to clone the PathFX repository directly to run this code: https://github.com/jenwilson521/PathFX
- a file of all drugnames and synonyms, '../drugbank_vocabulary.csv', which can be directly downloaded from DrugBank after registration
- ranked predictions from Gravina et al. PLoS Comp Bio 2022, '../sweetlead_ranked.xlsx', which can be obtained by emailing the authors.

##2. Ridge regression, heatmap image for the in vivo screen: "zebrafish_ridgereg.ipynb", "Heatmap_Zebrafish.ipynb"

Requirements:
- PathFX network + Phenoscore data, "df_zebrafish.csv"

##3. Ridge regression, heatmap image for the clinical data: "antidep_ridgereg.ipynb", "Heatmap_Antidepressants.ipynb"

Requirements:
- PathFX network + clinical OR data, "21AntidepressantData.csv"
- To create the above file, a directory of PathFv2 network results is needed, "results/antidepressants/". This is a subset of the results generated for step 1a.

##4. Plotting top coefficients for each screen: "top_ridge_coeffs.ipynb"

Requirements:
- formatted coefficient files, "insilico_phagocyt_ridgeregcoeffs.xlsx", "'efficacy_coef - efficacy_coef.csv'", "zebrafish_ridgeregcoeffs.xlsx"

##5. Plotting results of GO enrichment, "GO_Enrichment_Screen_Comparison.ipynb", "GO_Graphs_with_P=values.ipynb"

Requirements:
- GO results as generated and saved from EnrichR site: "GO clinical favorable.csv", "GO clinical unfavorable.csv", "GO insilico favorable.csv" ,"GO insilico unfavorable.csv", "GO zebrafish favorable.csv", "GO zebrafish unfavorable.csv"



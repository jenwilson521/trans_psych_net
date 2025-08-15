# trans_psych_net

_code for recreating results, most will require updated file paths_

[1] Running PathFX for all drugs, finding associations to phenotypes of interest, "search_psych.py"

Requirements:
- a copy of "search_DrugBank_for_phen_abbrev.py" moved into the PathFX/scripts/ directory after cloning
- results directory "PathFX/results/alldrugbank" containing PathFXv2 networks for all drugs in DrugBank (version 5.1.6). This can be created by running the script, "run_PathFX_allDrugBank.py" after moving it into PathFX/scripts/. You will need to clone the PathFX repository directly to run this code: https://github.com/jenwilson521/PathFX and will also need to create two objects after downloading the DrugBank data. Once access to the DrugBank database (version 5.1.6) is obtained, we use example parsing code from the link below to create other inputs to our analysis.
https://github.com/dhimmel/drugbank/blob/gh-pages/parse.ipynb and we created an excel spreadsheet using pandas with the default output of this script, and a second sheet, 'DrugsToTargets', which is a pandas table linking 'DrugBank_ID' to a comma-joined list of all unique drug targets. This second sheet is used in subsequent steps. The second data object, '../drugbank_vocabulary.csv',  can be directly downloaded from DrugBank after registration.

[2] Summarizing PathFX associations to psychiatric diseases, "figure1B.ipynb", "figures1C_1D.ipynb", "figure1E_ATC_codes.ipynb"				

Requirements:
- a formatted file linking DrugBank identifiers to ATC codes. This will be contained in the above Excel generated from the DrugBank .xml database.
- results of PathFX analysis explained above, "../supplemental_files/Supplemental File 1_ SCZ_PSCZ_MDD_RD_UD_BPDv3.xlsx"

[3] Logistic and ridge regression for in silico screen, heatmap for in silico screen: "in_silico_logreg.ipynb", "in_silico_ridgereg.ipynb", "Heatmap_Insilico.ipynb"

Requirements:
- results directory "alldrugbank" containing PathFXv2 networks for all drugs in DrugBank. This can be created by following step [1]
- a file of all drugnames and synonyms, '../drugbank_vocabulary.csv'
- ranked predictions from Gravina et al. PLoS Comp Bio 2022, '../sweetlead_ranked.xlsx', which can be obtained by emailing the authors. We provide code for parsing the raw data in '/code/in_silico_logreg.ipynb'

[4] Ridge regression, heatmap image for the in vivo screen: "zebrafish_ridgereg.ipynb", "Heatmap_Zebrafish.ipynb"

Requirements:
- PathFX network + Phenoscore data, "df_zebrafish.csv", we provide code ("read_and_make_df_zebrafish.py") for parsing the raw data file, "Figure_2b.csv", which can be downloaded from https://zenodo.org/records/3336616#.YrzrrRPMLLA 
 

[5] Ridge regression, heatmap image for the clinical data: "antidep_ridgereg.ipynb", "Heatmap_Antidepressants.ipynb"

Requirements:
- PathFX network + clinical OR data, "21AntidepressantData.csv". We manually created this .csv file and provide it in the /code/ directory.
- To create the above file, a directory of PathFv2 network results is needed, "results/antidepressants/". This is a subset of the results generated for step 1a.

[6] Plotting top coefficients for each screen: "top_ridge_coeffs.ipynb"

Requirements:
- formatted coefficient files, "insilico_phagocyt_ridgeregcoeffs.xlsx", "'efficacy_coef - efficacy_coef.csv'", "zebrafish_ridgeregcoeffs.xlsx"

[7] Plotting results of GO enrichment, "GO_Enrichment_Screen_Comparison.ipynb", "GO_Graphs_with_P=values.ipynb"

Requirements:
- GO results as generated and saved from EnrichR site: "GO clinical favorable.csv", "GO clinical unfavorable.csv", "GO insilico favorable.csv" ,"GO insilico unfavorable.csv", "GO zebrafish favorable.csv", "GO zebrafish unfavorable.csv"



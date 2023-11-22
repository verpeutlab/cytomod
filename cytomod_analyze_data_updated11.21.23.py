'''
Welcome to the CytoMod example code!

The full paper describing the method can be found at
https://www.frontiersin.org/articles/10.3389/fimmu.2019.01338/

The CytoMod folder contains a folder named data_files/data that contains
files named cytokine_data.xlsx and patient_data.xlsx, which hold
the data for the code's analysis.
You can replace these files with your own data files
while following the format instructions bellow.

# cytokine_data.xlsx: the first column is the subject IDs
(named PTID in the example file) which will be converted
to row indexes. If your dataset has no subject IDs,
change 'indexCol=0' to 'indexCol=None', in both cy_data and patients_data
initialization under the "Get data" header.
The next columns are your cytokines data. Each column should have
the raw cytokine measurment for the subject indicated in the specific row.

# patient_data.xlsx: Optional file for patient outcomes, for the associations
to outcomes analysis. The first column is the subject IDs
(named PTID in the example file). The instructions regarding
the IDs are the same as for the cytokine_data.xlsx file. This
data-frame should contain outcome variables to be analyzed in the
associations to outcomes analysis. It may also contain covariate
variables for controlling the regression models built for the
associations calculation.
*Make sure binary columns contain 0 and 1 values, or True and False values
(and cells with unknown values are left empty)

A folder named 'output' will be created by the code inside the
data_files folder. The code writes all results and figures into this folder.

Arguments to be manually defined:
* args.name_data
       Name of dataset/cohort (for writing files)
* args.name_compartment
       Name of compartment from which cytokines were extracted, e.g., serum (for writing files)
* args.cytokines
        List of cytokines to be analyzed. If None, will analyze all cytokines in the cytokine_data file.
* args.log_transform
       Boolean indicating whether to perform a log (base 10) transformation (True) or not (False).
* args.outcomes
       Optional. Names of outcome variables from the patients_data.xlsx data-frame to be analyzed.
       If list is left empty (i.e., []), will not perform the associations to outcomes analysis.
       This code supports binary or continuous outcome variables. Note - binary and continuous variables
       must be analyzed separately, i.e., you can analyze binary *or* continuous variables.
       (Binary outcome columns should only contain 0/1 or True/False values.)
* args.logistic
       Set to True if outcomes variables are binary (then, logistic regression will be used).
       Set to False if outcomes variables are continuous (then, linear regression will be used).
       *According to chosen value, see "associations" -> "output figures" for correct definition of
       colorscale_value and colorscale_labels variables, which define the colorscale for associations figures
* args.covariates
       Optional. Names of covariate variables (columns) from the patients_data.xlsx data-frame
       to be controlled for in the regression models. If list is left empty (i.e., []),
       will not control the associations to outcomes analysis with any covariate variables.
       Categorical covariates should be inserted as dummy variables.
* args.log_column_names
        List with names of covariate columns (from patient_data.xlsx) to be log-transformed,
        only if args.log_transform = True.
        If there are no columns you wish to transform, leave empty (i.e., [])
* args.max_testing_k
        Maximal number of clusters to test the gap statistic for.
* args.max_final_k
        The maximal number of clusters that can be chosen based on the
        gap statistic.
* args.recalculate_best_k
        Boolean. Set to True if you want the gap statistic (for finding the best K)
        to be recalcultaed in the current run.
        After calculation, the calculated best K will be saved to files,
        and used by the code until the next time you decide to recalculate them,
        or if the best K files are deleted. If no best K files are found, will
        recalculate best K anyway.
* args.seed
        Seed for random numbers stream set before cytomod calculations.

The code was written using the Anaconda3 Python interpreter and packages.
Recommended versions: Python 3.7.1, Pandas 0.23.4, Numpy 1.16.2
The palettable module (https://pypi.org/project/palettable/) must also be installed.
use: pip install palettable
'''

########### ------------ Imports & folders ------------- ###########
import os
import sys
import pandas as pd
sys.path.append(os.getcwd())
sys.path.append(os.path.join(os.getcwd(), 'cytomod', 'otherTools'))
import matplotlib.pyplot as plt
import cytomod
import cytomod.run_gap_statistic as gap_stat
import cytomod.assoc_to_outcome as outcome
from cytomod import plotting as cyplot
from hclusterplot import plotHColCluster
import tools
import numpy as np
import random

#%%
########### ------------------- Define manual arguments ----------------- ###########

args = tools.Object()

args.name_data = 'FLU09'
args.name_compartment = 'Plasma'

args.log_transform = False
args.max_testing_k = 8
args.max_final_k = 6 # Must be <= max_testing_k
args.recalculate_modules = False
args.outcomes = ['MDPV'] # names of outcome columns
args.logistic = True # True if outcomes are binary. False if outcomes are continuous. **## see "associations" -> "output figures" for correct definition of colorscale_value and colorscale_labels variables
args.covariates = ['Sex'] # names of regression covariates to control for
args.log_column_names = [] # or empty list: []
args.cytokines = None # if none, will take all

args.seed = 1234
#%%
########### ------------------- Define Output Folders ----------------- ###########

# Create output folders (if they don't already exist)
args.path_files = os.path.join(os.getcwd(), 'data_files')

args.paths = {'files': os.path.join(os.getcwd(), 'data_files'),
              'data': os.path.join(os.getcwd(), 'data_files', 'data'),
              'gap_statistic': os.path.join(os.getcwd(), 'data_files', 'output', 'gap_statistic'),
              'clustering': os.path.join(os.getcwd(), 'data_files', 'output', 'clustering'),
              'clustering_info': os.path.join(os.getcwd(), 'data_files', 'output', 'clustering', 'info'),
              'clustering_figures': os.path.join(os.getcwd(), 'data_files', 'output', 'clustering', 'figures'),
              'correlation_figures': os.path.join(os.getcwd(), 'data_files', 'output', 'correlations'),
              'association_figures': os.path.join(os.getcwd(), 'data_files', 'output', 'associations'),
              }

tools.create_folder(args.paths['gap_statistic'])
tools.create_folder(args.paths['clustering_info'])
tools.create_folder(args.paths['clustering_figures'])
tools.create_folder(args.paths['correlation_figures'])
tools.create_folder(args.paths['association_figures'])
#%%
########### ------------------- Assert args are valid ----------------- ###########

assert type(args.name_data) is str
assert type(args.name_compartment) is str
assert type(args.log_transform) is bool
assert type(args.logistic) is bool
assert type(args.max_testing_k) is int
assert type(args.max_final_k) is int
assert args.max_final_k <= args.max_testing_k
assert type(args.outcomes) is list
assert type(args.covariates) is list

for col_name in args.outcomes + args.covariates + args.log_column_names:
    assert type(col_name) is str
    tools.assert_column_exists_in_path(os.path.join(args.paths['data'], 'patient_data.xlsx'),
                                       col_name)
#%%
########### ------------ Get data ------------- ###########

# cytokines
cy_data = tools.read_excel(os.path.join(args.paths['data'], 'cytokine_data.xlsx'),
                           indexCol=0)
cy_data.dropna(axis='index', how='all', inplace=True)

if args.cytokines is None:
    args.cytokines = list(cy_data.columns)
else:
    cy_data = cy_data[args.cytokines]

# patients info
if args.outcomes != []:
    patient_data = tools.read_excel(os.path.join(args.paths['data'], 'patient_data.xlsx'),
                                    indexCol=0)
    patient_data.dropna(axis='index', how='all', inplace=True)

    # Check if args.logistic flag is correct for all outcomes defined in args.outcomes
    for outcome_col in args.outcomes:
        is_logistic = np.isin(patient_data[outcome_col].unique(), [0, 1]).all()  # checks if the data in outcomes column is binary (0,1 or true,false)
        if args.logistic != is_logistic: # mismatch! check which case
            if args.logistic:
                raise Exception('args.logistic defined as True. '
                                'However, outcome variable ' + outcome_col + ' seems '
                                'to be continuous and not binary. Please check and fix!')
            else:
                raise Exception('args.logistic defined as False. '
                                'However, outcome variable ' + outcome_col + ' seems '
                                'to be binary and not continuous. Please check and fix!')

# log transform cytokines and args.log_column_names
if args.log_transform:
    cy_data = np.log10(cy_data)

    if args.log_column_names != [] and args.outcomes != []:
        for col_name in args.log_column_names:
            new_col_name = 'log_' + col_name

            # log transform variable
            patient_data[new_col_name] = np.log10(patient_data[col_name])

            # replace column with new log transformed column
            if col_name in args.covariates:
                args.covariates.remove(col_name)
                args.covariates.append(new_col_name)
#%%
########### ------------ Adjust and Cluster ------------- ###########

do_recalculate = args.recalculate_modules or \
        not os.path.exists(os.path.join(args.paths['clustering'], 'cyto_mod_adj.dill'))

# Get best K. If first time or args.recalculate_modules=True - compute best K. Otherwise - read from file.
if do_recalculate:
    random.seed(args.seed)
    cyto_mod_adj = cytomod.cytomod_class(args.name_data, args.name_compartment, True, cy_data)
    cyto_mod_abs = cytomod.cytomod_class(args.name_data, args.name_compartment, False, cy_data)

    bestK = {}
    bestK['adj'] = gap_stat.getBestK(cyto_mod_adj.cyDf,
                                       max_testing_k = args.max_testing_k,
                                       max_final_k=args.max_final_k,
                                       save_fig_path=os.path.join(args.paths['gap_statistic'], 'gap_stat_adj.png'))

    bestK['abs'] = gap_stat.getBestK(cyto_mod_abs.cyDf,
                                       max_testing_k=args.max_testing_k,
                                       max_final_k=args.max_final_k,
                                       save_fig_path=os.path.join(args.paths['gap_statistic'], 'gap_stat_abs.png'))

    tools.write_DF_to_excel(os.path.join(args.paths['clustering'], 'bestK.xlsx'), bestK)

    # Cluster
    cyto_mod_adj.cluster_cytokines(K=bestK['adj'])
    cyto_mod_abs.cluster_cytokines(K=bestK['abs'])
    tools.write_to_dill(os.path.join(args.paths['clustering'], 'cyto_mod_adj.dill'), cyto_mod_adj)
    tools.write_to_dill(os.path.join(args.paths['clustering'], 'cyto_mod_abs.dill'), cyto_mod_abs)
else:
    # Get modules from storage
    bestK = tools.read_excel(os.path.join(args.paths['clustering'], 'bestK.xlsx'))
    bestK = dict(bestK['value'])
    cyto_mod_adj = tools.read_from_dill(os.path.join(args.paths['clustering'], 'cyto_mod_adj.dill'))
    cyto_mod_abs = tools.read_from_dill(os.path.join(args.paths['clustering'], 'cyto_mod_abs.dill'))

cyto_modules = {'adj': cyto_mod_adj, 'abs': cyto_mod_abs}
#%%
########### ------------ Output clustering results ------------- ###########

# Output clustering results
if do_recalculate:
    for cyto_object in cyto_modules.values():
        cyplot.plotMeanCorr(cyto_object.withMean, cyto_object.meanS.name,
                            cyList=sorted(cyto_object.cyDf.columns), figsize=(10, 6),
                            save_path=os.path.join(args.paths['correlation_figures'],
                                                   '%s_cy_mean_correlation.png' % cyto_object.name))

        plotHColCluster(cyto_object.cyDf, method='complete',
                             metric='pearson-signed', figsize=(10, 6),
                             save_path=os.path.join(args.paths['correlation_figures'],
                                                    '%s_correlation_heatmap.png' % cyto_object.name))

        cytomod.io.write_modules(cyto_object, args.paths['clustering_info'])
        cytomod.io.plot_modules(cyto_object, args.paths['clustering_figures'],
                                heatmap_figsize=(10, 6))
#%%
########### ------------ Associations ------------- ###########

# Analyze associations
if args.outcomes != []:
    df_outcomes = patient_data[args.outcomes + args.covariates].join(cy_data)

    #### Absolute outcomes calculation

    mod_outcome_abs_df = outcome.outcomeAnalysis(cyto_modules['abs'], patient_data,
                        analyzeModules=True,
                        outcomeVars=args.outcomes,
                        adjustmentVars=args.covariates,
                        logistic=args.logistic,
                        standardize=True)
    cy_outcome_abs_df = outcome.outcomeAnalysis(cyto_modules['abs'], patient_data,
                        analyzeModules=False,
                        outcomeVars=args.outcomes,
                        adjustmentVars=args.covariates,
                        logistic=args.logistic,
                        standardize=True)

    #### Adjusted outcomes calculation

    mod_outcome_adj_df = outcome.outcomeAnalysis(cyto_modules['adj'], patient_data,
                        analyzeModules=True,
                        outcomeVars=args.outcomes,
                        adjustmentVars=args.covariates,
                        logistic=args.logistic,
                        standardize=True)
    cy_outcome_adj_df = outcome.outcomeAnalysis(cyto_modules['adj'], patient_data,
                        analyzeModules=False,
                        outcomeVars=args.outcomes,
                        adjustmentVars=args.covariates,
                        logistic=args.logistic,
                        standardize=True)

#%%
    ########### ------------ Output Figures ------------- ###########

    # If args.logistic = False, these variables must be symmetric around 0.
    # For example:
    #       colorscale_values = [-0.8, -0.4, 0, 0.4, 0.8]
    #       colorscale_labels = [-0.8, -0.4, 0, 0.4, 0.8]
    # If args.logistic = True, these variables must be symmetric around 1.
    # For example:
    #       colorscale_values = [1 / 2.5, 1 / 2, 1 / 1.5, 1, 1.5, 2, 2.5]
    #       colorscale_labels = ['1/2.5', '1/2', '1/1.5', 1, 1.5, 2, 2.5]
    colorscale_values = [1 / 2.5, 1 / 2, 1 / 1.5, 1, 1.5, 2, 2.5]
    colorscale_labels = ['1/2.5', '1/2', '1/1.5', 1, 1.5, 2, 2.5]

    # Does colorscale_values range include all regression coefficients issued by outcomeAnalysis?
    outcome.check_colorscale_range(colorscale_values, args.logistic,
                                   cy_outcome_abs_df, mod_outcome_abs_df, cy_outcome_adj_df, mod_outcome_adj_df)

    #### Absolute
    outcome.plotResultSummary(cyto_modules['abs'],
                              mod_outcome_abs_df,
                              cy_outcome_abs_df,
                              args.outcomes,
                              logistic=args.logistic,
                              fdr_thresh_plot=0.2,
                              compartmentName=args.name_compartment,
                              figsize=(6,9),
                              scale_values=colorscale_values,
                              scale_labels=colorscale_labels,
                              save_fig_path=os.path.join(args.paths['association_figures'], 'associations_abs.png'))

    #### Adjusted
    outcome.plotResultSummary(cyto_modules['adj'],
                              mod_outcome_adj_df,
                              cy_outcome_adj_df,
                              args.outcomes,
                              logistic=args.logistic,
                              fdr_thresh_plot=0.2,
                              compartmentName=args.name_compartment,
                              figsize=(6,9),
                              scale_values=colorscale_values,
                              scale_labels=colorscale_labels,
                              save_fig_path=os.path.join(args.paths['association_figures'], 'associations_adj.png'))

#%%
#---Doesn't work below---
    # Adjusted modules
    outcome.printTable(mod_outcome_adj_df,
                       title=args.name_compartment + ' (Adjusted)',
                       fdr_output_limit=1, fwer_output_limit=1, pval_output_limit=1,
                       output_file_path=os.path.join(args.paths['association_figures'], 'associations_adj_mod_pvals.csv'))
 
    # Absolute modules
    outcome.printTable(mod_outcome_abs_df,
                       title=args.name_compartment + ' (Absolute)',
                       fdr_output_limit=1, fwer_output_limit=1, pval_output_limit=1,
                       output_file_path=os.path.join(args.paths['association_figures'], 'associations_abs_mod_pvals.csv'))

    # Adjusted cytokines
    outcome.printTable(cy_outcome_adj_df,
                       title=args.name_compartment + ' (Adjusted)',
                       fdr_output_limit=1, fwer_output_limit=1, pval_output_limit=1,
                       output_file_path=os.path.join(args.paths['association_figures'], 'associations_adj_cy_pvals.csv'))

    # Absolute cytokines
    outcome.printTable(cy_outcome_abs_df,
                       title=args.name_compartment + ' (Absolute)',
                       fdr_output_limit=1, fwer_output_limit=1, pval_output_limit=1,
                       output_file_path=os.path.join(args.paths['association_figures'], 'associations_abs_cy_pvals.csv'))
    
    


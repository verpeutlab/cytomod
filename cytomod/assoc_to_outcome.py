import os.path as op
import matplotlib.pyplot as plt
import scipy.stats.stats
import numpy as np
import palettable
import itertools
import seaborn as sns


sns.set(style='darkgrid', font_scale=1.5, palette='muted')

####----------------- helpers --------------####
from copy import deepcopy
import seaborn as sns
# from glm_compare import compare_lr_test
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.api as sm
from dfprint import toPDF, greek2latex
import warnings

sns.set(style='darkgrid', palette='muted', font_scale=2.0)

def GLMResults(df, outcome, predictors, adj=[], logistic=True):
    if logistic:
        family = sm.families.Binomial()
        coefFunc = np.exp
        cols = ['OR', 'LL', 'UL', 'pvalue', 'Diff', 'N']
    else:
        family = sm.families.Gaussian()
        coefFunc = lambda x: x
        cols = ['Coef', 'LL', 'UL', 'pvalue', 'Diff', 'N']

    k = len(predictors)
    assoc = np.zeros((k, 6))
    params = []
    pvalues = []
    resObj = []
    for i, predc in enumerate(predictors):
        exogVars = list(set([predc] + adj))
        tmp = df[[outcome] + exogVars].dropna()

        model = sm.GLM(endog=tmp[outcome].astype(float), exog=sm.add_constant(tmp[exogVars].astype(float)),
                       family=family)
        try:
            res = model.fit()
            assoc[i, 0] = coefFunc(res.params[predc])
            assoc[i, 3] = res.pvalues[predc]
            assoc[i, 1:3] = coefFunc(res.conf_int().loc[predc])
            assoc[i, 4] = tmp[predc].loc[tmp[outcome] == 1].mean() - tmp[predc].loc[tmp[outcome] == 0].mean()
            params.append(res.params.to_dict())
            pvalues.append(res.pvalues.to_dict())
            resObj.append(res)
        except sm.tools.sm_exceptions.PerfectSeparationError:
            assoc[i, 0] = np.nan
            assoc[i, 3] = 0
            assoc[i, 1:3] = [np.nan, np.nan]
            assoc[i, 4] = tmp[predc].loc[tmp[outcome] == 1].mean() - tmp[predc].loc[tmp[outcome] == 0].mean()
            params.append({k: np.nan for k in [predc] + adj})
            pvalues.append({k: np.nan for k in [predc] + adj})
            resObj.append(None)
            print('PerfectSeparationError: %s with %s' % (predc, outcome))
        assoc[i, 5] = tmp.shape[0]
    outDf = pd.DataFrame(assoc[:, :6], index=predictors, columns=cols)
    outDf['params'] = params
    outDf['pvalues'] = pvalues
    outDf['res'] = resObj
    return outDf

def outcomeAnalysis(cytomod_obj, patient_data,
                    analyzeModules=True,
                    outcomeVars=[],
                    adjustmentVars=[],
                    standardize=True,
                    logistic=True):
    """Do these FLU-positive clusters correlate with outcome,
    with/without adjustment for bacterial coinfection?"""

    modStr = 'Module' if analyzeModules else 'Analyte'
    resL = []
    standardizeFunc = lambda col: (col - np.nanmean(col)) / np.nanstd(col)
    for outcome in outcomeVars:
        """Logistic regression on outcome"""
        if analyzeModules:
            dataDf = cytomod_obj.modDf
        else:
            dataDf = cytomod_obj.cyDf
            if standardize:  # standardize cytokine values
                dataDf = dataDf.apply(standardizeFunc)

        predictors = dataDf.columns
        data_outcome_Df = patient_data[outcomeVars + adjustmentVars].join(dataDf).copy() # why outcomeVars and not outcome? make sure we can change this

        if standardize:
            if not logistic: # if continuous outcomes, standardize to normal distribution Z
                data_outcome_Df[[outcome]] = data_outcome_Df[[outcome]].apply(standardizeFunc)

            if adjustmentVars != []: # if there are any covariates
                for covariate in adjustmentVars: # standardize each to normal distribution Z
                    if len(data_outcome_Df[covariate].unique()) > 2: # but only if the covariate is not binary
                        data_outcome_Df[[covariate]] = data_outcome_Df[[covariate]].apply(standardizeFunc)

        tmpres = GLMResults(data_outcome_Df, outcome, predictors, adj=adjustmentVars, logistic=logistic)
        tmpres['Outcome'] = outcome
        tmpres['Compartment'] = cytomod_obj.sampleStr
        tmpres['Adjusted'] = 'Yes' if cytomod_obj.adjusted else 'No'
        tmpres['Fold-diff'] = np.exp(tmpres['Diff'])
        tmpres[modStr] = predictors
        resL.append(tmpres)

    resDf = pd.concat(resL, axis=0, ignore_index=True)
    return resDf


####### outcome #######
def mapColors2Labels(labels, setStr='MapSet', cmap=None):
    """Return pd.Series of colors based on labels"""
    if cmap is None:
        N = max(3, min(12, len(np.unique(labels))))
        cmap = palettable.colorbrewer.get_map(setStr, 'Qualitative', N).mpl_colors
        """Use B+W colormap"""
    cmapLookup = {k:col for k, col in zip(sorted(np.unique(labels)), itertools.cycle(cmap))}
    return labels.map(cmapLookup.get)

def adjust_pvals(res_df):
    res_df = deepcopy(res_df)
    res_df.loc[:, 'FWER'] = sm.stats.multipletests(res_df.pvalue.values, method='holm')[1]
    res_df.loc[:, 'FDR'] = sm.stats.multipletests(res_df.pvalue.values, method='fdr_bh')[1]
    res_df.loc[:, 'Bonferroni'] = sm.stats.multipletests(res_df.pvalue.values, method='bonferroni')[1]
    return res_df

''' For logistic=True, scale_values should be symmetric around 1 in log scale.
For example, [1/2.5, 1/2, 1/1.5, 1, 1.5, 2, 2.5]. Labels will then be
['1/2.5', '1/2', '1/1.5', 1, 1.5, 2, 2.5].
For logistic=False, scale_values should be symmetric around 0. 
For example, [-0.8, -0.4, 0, 0.4, 0.8]. 
Labels will then be [-0.8, -0.4, 0, 0.4, 0.8]. '''
def plotResultSummary(cytomod_obj,
                      mod_res_df,
                      cy_res_df,
                      outcomeVars,
                      fdr_thresh_plot=0.2,
                      compartmentName='BS',
                      showScalebar=True,
                      figsize=(6,9),
                      save_fig_path=None,
                      logistic=True,
                      scale_values=[1 / 2.5, 1 / 2, 1 / 1.5, 1, 1.5, 2, 2.5],
                      scale_labels=['1/2.5', '1/2', '1/1.5', 1, 1.5, 2, 2.5]):
    mod_res_df = mod_res_df.copy()
    cy_res_df = cy_res_df.copy()

    mod_res_df.loc[:, 'Name'] = mod_res_df['Module']
    cy_res_df.loc[:, 'Name'] = cy_res_df['Analyte']

    cy_res_df = adjust_pvals(cy_res_df)
    mod_res_df = adjust_pvals(mod_res_df)

    name2mod = lambda a: '%s%1.0f' % (compartmentName, cytomod_obj.labels[a])

    cy_res_df.loc[:, 'Module'] = cy_res_df['Analyte'].map(name2mod)

    if logistic:
        cols = ['Outcome', 'Name', 'Module', 'Fold-diff', 'OR', 'N', 'FWER', 'FDR']
    else:
        cols = ['Outcome', 'Name', 'Module', 'Fold-diff', 'Coef', 'N', 'FWER', 'FDR']

    hDf = pd.concat((mod_res_df[cols], cy_res_df[cols]), axis=0)
    hDf.loc[:, 'isAnalyte'] = (hDf['Module'] != hDf['Name'])
    order = hDf[['Module', 'Name', 'isAnalyte']].drop_duplicates().sort_values(by=['Module', 'isAnalyte', 'Name'])
    fdrH = hDf.pivot(index='Name', columns='Outcome', values='FDR').loc[order.Name, outcomeVars]
    fdrH = fdrH.fillna(1)
    fwerH = hDf.pivot(index='Name', columns='Outcome', values='FWER').loc[order.Name, outcomeVars]
    fwerH = fwerH.fillna(1)

    censorInd = fdrH.values > fdr_thresh_plot
    fdrH.values[censorInd] = 1.


    cmap = palettable.colorbrewer.diverging.PuOr_9_r.mpl_colormap

    if logistic:
        foldH = hDf.pivot(index='Name', columns='Outcome', values='Fold-diff').loc[order.Name, outcomeVars]
        foldH.values[censorInd] = 1.
        foldH = foldH.fillna(1)

        vals = np.log(foldH.values)
        pcParams = dict(vmin=np.log(scale_values[0]), vmax=np.log(scale_values[-1]), cmap=cmap)

        scaleLabel = 'Odds Ratio'
        ytl = np.array(scale_labels)
        yt = np.log(scale_values)

        plt.figure(figsize=figsize)
        figh = plt.gcf()
        plt.clf()
        axh = figh.add_subplot(plt.GridSpec(1, 1, left=0.6, bottom=0.05, right=0.95, top=0.85)[0, 0])
        axh.grid(None)
        pcolOut = plt.pcolormesh(vals, **pcParams)
        plt.yticks(())
        plt.xticks(np.arange(fdrH.shape[1]) + 0.5, fdrH.columns, size=11, rotation=90)
        axh.xaxis.set_ticks_position('top')
        plt.xlim((0, fdrH.shape[1]))
        plt.ylim((0, fdrH.shape[0]))
        axh.invert_yaxis()
        for cyi, cy in enumerate(foldH.index):
            for outi, out in enumerate(foldH.columns):
                if fwerH.loc[cy, out] < 0.0005:
                    ann = '***'
                elif fwerH.loc[cy, out] < 0.005:
                    ann = '**'
                elif fwerH.loc[cy, out] < 0.05:
                    ann = '*'
                else:
                    ann = ''
                if not ann == '':
                    plt.annotate(ann, xy=(outi + 0.5, cyi + 0.75), weight='bold', size=14, ha='center', va='center')

        """Colorbar showing module membership: Add labels, make B+W"""
        cbAxh = figh.add_subplot(plt.GridSpec(1, 1, left=0.5, bottom=0.05, right=0.59, top=0.85)[0, 0])
        cbAxh.grid(None)
        cmap = [(0.3, 0.3, 0.3),
                (0.7, 0.7, 0.7)]
        cbS = mapColors2Labels(order.set_index('Name')['Module'], cmap=cmap)
        _ = cbAxh.imshow([[x] for x in cbS.values], interpolation='nearest', aspect='auto', origin='lower')
        plt.ylim((0, fdrH.shape[0]))
        plt.yticks(np.arange(fdrH.shape[0]), fdrH.index, size=11)
        plt.xlim((0, 0.5))
        plt.ylim((-0.5, fdrH.shape[0] - 0.5))
        plt.xticks(())
        cbAxh.invert_yaxis()
    else:
        betaVals = hDf.pivot(index='Name', columns='Outcome', values='Coef').loc[order.Name, outcomeVars]  # LIEL
        betaVals.values[censorInd] = 0.
        vals = betaVals.values
        pcParams = dict(vmin=scale_values[0], vmax=scale_values[-1], cmap=cmap)
        scaleLabel = 'Beta Coefficient'
        ytl = np.array(scale_labels)
        yt = np.array(scale_values)
        plt.figure(figsize=figsize)
        figh = plt.gcf()
        plt.clf()
        axh = figh.add_subplot(plt.GridSpec(1, 1, left=0.6, bottom=0.05, right=0.95, top=0.85)[0, 0])
        axh.grid(None)
        pcolOut = plt.pcolormesh(vals, **pcParams)
        plt.yticks(())

        plt.xticks(np.arange(betaVals.shape[1]) + 0.5, betaVals.columns, size=11, rotation=90)
        axh.xaxis.set_ticks_position('top')
        plt.xlim((0, betaVals.shape[1]))
        plt.ylim((0, betaVals.shape[0]))
        axh.invert_yaxis()
        for cyi, cy in enumerate(betaVals.index):
            for outi, out in enumerate(betaVals.columns):
                if fwerH.loc[cy, out] < 0.0005:
                    ann = '***'
                elif fwerH.loc[cy, out] < 0.005:
                    ann = '**'
                elif fwerH.loc[cy, out] < 0.05:
                    ann = '*'
                else:
                    ann = ''
                if not ann == '':
                    plt.annotate(ann, xy=(outi + 0.5, cyi + 0.75), weight='bold', size=14, ha='center', va='center')

        """Colorbar showing module membership: Add labels, make B+W"""
        cbAxh = figh.add_subplot(plt.GridSpec(1, 1, left=0.5, bottom=0.05, right=0.59, top=0.85)[0, 0])
        cbAxh.grid(None)
        cmap = [(0.3, 0.3, 0.3),
                (0.7, 0.7, 0.7)]
        cbS = mapColors2Labels(order.set_index('Name')['Module'], cmap=cmap)
        _ = cbAxh.imshow([[x] for x in cbS.values], interpolation='nearest', aspect='auto', origin='lower')
        plt.ylim((0, betaVals.shape[0]))
        plt.yticks(np.arange(betaVals.shape[0]), betaVals.index, size=11)
        plt.xlim((0, 0.5))
        plt.ylim((-0.5, betaVals.shape[0] - 0.5))
        plt.xticks(())
        cbAxh.invert_yaxis()

    for lab in order['Module'].unique():
        y = np.mean(np.arange(order.shape[0])[np.nonzero(order['Module'] == lab)]) - 0.5
        plt.annotate(lab, xy=(0.25, y), ha='center', va='center', rotation=90, color='white', size=12)

    """Scale colorbar"""
    if showScalebar:
        scaleAxh = figh.add_subplot(plt.GridSpec(1, 1, left=0.1, bottom=0.87, right=0.2, top=0.98)[0, 0])
        cb = figh.colorbar(pcolOut, cax=scaleAxh, ticks=yt)
        cb.set_label(scaleLabel, size=9)
        cb.ax.set_yticklabels(ytl, fontsize=8)

    if save_fig_path is None:
        plt.show()
    else:
        plt.savefig(save_fig_path)


''' Print association pvalues table to file and/or console.
    @ output_file_path possible file types: pdf, csv. 
      *If None - no file will be saved. 
      *for pdf, will also write .tex files.
    @ fdr_output_limit / fwer_output_limit / pval_output_limit:
    function only writes modules/cytokines with values under these limits.'''
def printTable(assoc_res_df, title='', output_file_path = True, print_to_console=True,
               fdr_output_limit=1, fwer_output_limit=1, pval_output_limit=1):

    resDf = assoc_res_df.copy()

    if 'OR' in resDf.columns: # logistic regression
        value_name = 'OR'
    elif 'Coef' in resDf.columns: # linear regression
        value_name = 'Coef'
    else:
        raise Exception('printTable: assoc_res_df invalid. Must contain OR or Coef column.')

    # define columns to be written to file
    if 'Module' in resDf.columns:
        cols = ['Outcome', 'Module', value_name, 'pvalue', 'FWER', 'FDR']
    elif 'Analyte' in resDf.columns:
        cols = ['Outcome', 'Analyte', value_name, 'pvalue', 'FWER', 'FDR']
    else:
        raise Exception('printTable: assoc_res_df format error: '
                        'must contain column Module or Analyte')

    # get adjusted pvalues
    resDf['FWER'] = sm.stats.multipletests(resDf.pvalue.values, method='holm')[1]
    resDf['FDR'] = sm.stats.multipletests(resDf.pvalue.values, method='fdr_bh')[1]

    # Get relevant rows by value limits defined
    ind = (resDf['FDR'] <= fdr_output_limit) & \
          (resDf['FWER'] <= fwer_output_limit) & \
          (resDf['pvalue'] <= pval_output_limit)

    # print tables to console
    if print_to_console:
        print('\n' + title)
        print('==================================')
        float_format = lambda f: '%1.2g' % f
        print(resDf[cols].loc[ind].sort_values(by='pvalue').to_string(float_format=float_format))

    # change greek letters format for files
    if 'Analyte' in resDf.columns:
        resDf['Analyte'] = resDf['Analyte'].map(greek2latex)

    # output tables to files
    if output_file_path[-3:] == 'pdf':

        # write pdf and tex files
        toPDF(resDf[cols].loc[ind].sort_values(by='pvalue'),
              output_file_path,titStr=title, float_format='%1.3g')

    elif output_file_path[-3:] == 'csv':
        # write csv file
        resDf[cols].loc[ind].sort_values(by='pvalue').to_csv(output_file_path,
                                                             index=False,
                                                             float_format='%1.3g')

    elif output_file_path is None:
        pass  # do nothing

    else:
        raise Exception('printTable: unknown output file format.')

def check_colorscale_range(colorscale_values, logistic,
                           cy_outcome_abs_df, mod_outcome_abs_df,
                           cy_outcome_adj_df, mod_outcome_adj_df):
    coefficient_col = 'OR' if logistic else 'Coef'
    min_coefficient = min(min(cy_outcome_abs_df[coefficient_col]),
                          min(mod_outcome_abs_df[coefficient_col]),
                          min(cy_outcome_adj_df[coefficient_col]),
                          min(mod_outcome_adj_df[coefficient_col]))
    max_coefficient = max(max(cy_outcome_abs_df[coefficient_col]),
                          max(mod_outcome_abs_df[coefficient_col]),
                          max(cy_outcome_adj_df[coefficient_col]),
                          max(mod_outcome_adj_df[coefficient_col]))
    if min_coefficient < colorscale_values[0]:
        warnings.warn('Minimal regression coefficient calculated is ' + str(min_coefficient) +
                      ' , however, minimal value for associations figures colorscale defined '
                      'is ' + str(colorscale_values[0]) + ' . It is advised to change the colorscale '
                      'range (in variable colorscale_values) to include the minimal value. (warning issued in check_colorscale_range)')
    if max_coefficient > colorscale_values[-1]:
        warnings.warn('Maximal regression coefficient calculated is ' + str(max_coefficient) +
                      ' , however, maximal value for associations figures colorscale defined '
                      'is ' + str(colorscale_values[-1]) + ' . It is advised to change the colorscale '
                      'range (in variable colorscale_values) to include the maximal value. (warning issued in check_colorscale_range)')

import pandas as pd
from os.path import join as opj

from hclusterplot import plotHColCluster
import cytomod as cy
from cytomod import plotting as cyplot
import palettable
from custom_legends import colorLegend
import matplotlib.pyplot as plt
import numpy as np
import os

__all__ = ['write_modules', 'plot_modules']

def write_modules(clust_object, folder):
    """rawDf: pd.DataFrame with cytokines as columns and sample IDs along the index"""
    clust_object.rCyDf.to_csv(opj(folder, '{name}_raw_log-conc.csv'.format(name=clust_object.name)))

    """cyDf: same as rawDf, but with missing values filled (if needed), log-transformed and possibly normalized"""
    clust_object.cyDf.to_csv(opj(folder, '{name}_normalized_conc.csv'.format(name=clust_object.name)))

    """dmatDf: pd.DataFrame representation of pairwise distance matrix of cytokines (index and columns of cytokines)"""
    clust_object.dmatDf.to_csv(opj(folder, '{name}_dmat.csv'.format(name=clust_object.name)))
    
    """pwrelDf: pd.DataFrame of pairwise cluster reliability (as a distance) from a bootstrap (index and columns of cytokines)"""
    clust_object.pwrel.to_csv(opj(folder, '{name}_pwrel_dmat.csv'.format(name=clust_object.name)))

    """labels: pd.Series containing cluster labels, indexed by cytokine"""
    clust_object.labels.to_csv(opj(folder, '{name}_cluster_labels.csv'.format(name=clust_object.name)))

    """modDf: pd.DataFrame of summarized module variables"""
    clust_object.modDf.to_csv(opj(folder, '{name}_normalized_modules.csv'.format(name=clust_object.name)))

    clust_object.rModDf.to_csv(opj(folder, '{name}_raw_modules.csv'.format(name=clust_object.name)))

def plot_modules(clust_object, folder):
    '''Plot cytomod object modules information'''

    """Hierarchical clustering heatmap"""
    plt.figure(41, figsize=(15.5, 9.5))
    # colInds = plotHColCluster(ds[s].cyDf, method='complete', metric='pearson-signed', col_labels=ds[s].labels, col_dmat=ds[s].dmatDf)
    colInds = plotHColCluster(clust_object.cyDf, method='complete', metric='pearson-signed',
                              col_labels=clust_object.labels,
                              save_path=os.path.join(folder, '%s_hierchical_clust_heatmap.png' % clust_object.name))
    # plt.figure(41).savefig(os.path.join(folder, '%s_hierchical_clust_heatmap.png' % clust_object.name))

    """Heatmaps of pairwise reliability"""
    plt.figure(43, figsize=(15.5, 9.5))
    colInds = cyplot.plotHierClust(1 - clust_object.pwrel, clust_object.Z, labels=clust_object.labels,
                                   titleStr='Pairwise reliability (%s)' % clust_object.name.replace('_', ' '), vRange=(0, 1))
    plt.figure(43).savefig(os.path.join(folder, '%s_pwrel.png' % clust_object.name), dpi=300)

    """color_label_legend"""
    plt.figure(48, figsize=(3, 3))
    plt.clf()
    axh = plt.subplot(1, 1, 1)
    colorLegend(palettable.colorbrewer.qualitative.Set3_6.mpl_colors, ['%s' % s for s in clust_object.modDf.columns],
                loc=10)
    axh.spines['right'].set_color('none')
    axh.spines['left'].set_color('none')
    axh.spines['top'].set_color('none')
    axh.spines['bottom'].set_color('none')
    axh.set_xticks([])
    axh.set_yticks([])
    axh.set_facecolor('white')

    plt.figure(48).savefig(os.path.join(folder, 'color_label_legend.png'), dpi=300)

    """Plot intra-module correlation"""
    plt.figure(50, figsize=(15, 9))
    for lab in list(cy.labels2modules(clust_object.labels, clust_object.dropped).keys()):
        cyplot.plotModuleCorr(clust_object.cyDf, clust_object.labels, lab, dropped=clust_object.dropped)
        plt.figure(50).savefig(os.path.join(folder, '%s_mod_corr_%s.png' % (clust_object.name, lab)), dpi=300)

    """Cytokine embedding"""
    plt.figure(901, figsize=(13, 9.7))
    cyplot.plotModuleEmbedding(clust_object.dmatDf, clust_object.labels, method='kpca')
    if len(np.unique(clust_object.labels)) > 2:
        colors = palettable.colorbrewer.get_map('Set1', 'qualitative', len(np.unique(clust_object.labels))).mpl_colors
    else:
        colors = [(0.8941176470588236, 0.10196078431372549, 0.10980392156862745),
                   (0.21568627450980393, 0.49411764705882355, 0.7215686274509804)]
    colorLegend(colors, ['%s%1.0f' % (clust_object.sampleStr, i) for i in np.unique(clust_object.labels)],
                loc='lower left')
    plt.figure(901).savefig(os.path.join(folder, '%sembed.png' % clust_object.name), dpi=300)

def plot_modules(clust_object, folder, heatmap_figsize=(15.5, 9.5)):
    '''Plot cytomod object modules information'''

    """Hierarchical clustering heatmap"""
    plt.figure(41, figsize=heatmap_figsize)
    # colInds = plotHColCluster(ds[s].cyDf, method='complete', metric='pearson-signed', col_labels=ds[s].labels, col_dmat=ds[s].dmatDf)
    colInds = plotHColCluster(clust_object.cyDf, method='complete', metric='pearson-signed',
                              col_labels=clust_object.labels,
                              save_path=os.path.join(folder, '%s_hier.png' % clust_object.name))
    # plt.figure(41).savefig(os.path.join(folder, '%s_hier.png' % clust_object.name))

    """Heatmaps of pairwise reliability"""
    plt.figure(43, figsize=heatmap_figsize)
    colInds = cyplot.plotHierClust(1 - clust_object.pwrel, clust_object.Z, labels=clust_object.labels,
                                   titleStr='Pairwise reliability (%s)' % clust_object.name.replace('_', ' '), vRange=(0, 1))
    plt.figure(43).savefig(os.path.join(folder, '%s_reliability.png' % clust_object.name), dpi=300)

    """color_label_legend"""
    plt.figure(48, figsize=(3, 3))
    plt.clf()
    axh = plt.subplot(1, 1, 1)
    colorLegend(palettable.colorbrewer.qualitative.Set3_6.mpl_colors, ['%s' % s for s in clust_object.modDf.columns],
                loc=10)
    axh.spines['right'].set_color('none')
    axh.spines['left'].set_color('none')
    axh.spines['top'].set_color('none')
    axh.spines['bottom'].set_color('none')
    axh.set_xticks([])
    axh.set_yticks([])
    axh.set_facecolor('white')
    plt.figure(48).savefig(os.path.join(folder, '%s_color_label_legend.png' % clust_object.name), dpi=300)

    """Plot intra-module correlation"""
    plt.figure(50, figsize=(15, 9))
    for lab in list(cy.labels2modules(clust_object.labels, clust_object.dropped).keys()):
        cyplot.plotModuleCorr(clust_object.cyDf, clust_object.labels, lab, dropped=clust_object.dropped)
        plt.figure(50).savefig(os.path.join(folder, '%s_modules_correlations_%s.png' % (clust_object.name, lab)), dpi=300)

    """Cytokine embedding"""
    plt.figure(901, figsize=(13, 9.7))
    cyplot.plotModuleEmbedding(clust_object.dmatDf, clust_object.labels, method='kpca')
    if len(np.unique(clust_object.labels)) > 2:
        colors = palettable.colorbrewer.get_map('Set1', 'qualitative', len(np.unique(clust_object.labels))).mpl_colors
    else:
        colors = [(0.8941176470588236, 0.10196078431372549, 0.10980392156862745),
                   (0.21568627450980393, 0.49411764705882355, 0.7215686274509804)]
    colorLegend(colors, ['%s%1.0f' % (clust_object.sampleStr, i) for i in np.unique(clust_object.labels)],
                loc='lower left')
    plt.figure(901).savefig(os.path.join(folder, '%s_embedding.png' % clust_object.name), dpi=300)

def plot_clustering_heatmap(clust_object, folder, figsize=(15.5, 9.5)):
    """Hierarchical clustering heatmap"""
    plt.figure(41, figsize=figsize)
    # colInds = plotHColCluster(ds[s].cyDf, method='complete', metric='pearson-signed', col_labels=ds[s].labels, col_dmat=ds[s].dmatDf)
    colInds = plotHColCluster(clust_object.cyDf, method='complete', metric='pearson-signed',
                              col_labels=clust_object.labels, figsize=figsize,
                              save_path=os.path.join(folder, '%s_hierchical_clust_heatmap.png' % clust_object.name))
    # plt.figure(41).savefig(os.path.join(folder, '%s_hierchical_clust_heatmap.png' % clust_object.name))

def plot_reliability(clust_object, folder, figsize=(15.5, 9.5)):
    """Heatmaps of pairwise reliability"""
    plt.figure(43, figsize=figsize)
    colInds = cyplot.plotHierClust(1 - clust_object.pwrel, clust_object.Z, labels=clust_object.labels,
                                   titleStr='Pairwise reliability (%s)' % clust_object.name.replace('_', ' '), vRange=(0, 1))
    plt.figure(43).savefig(os.path.join(folder, '%s_reliability.png' % clust_object.name), dpi=300)

def plot_color_legend(clust_object, folder):
    """color_label_legend"""
    plt.figure(48, figsize=(3, 3))
    plt.clf()
    axh = plt.subplot(1, 1, 1)
    colorLegend(palettable.colorbrewer.qualitative.Set3_6.mpl_colors, ['%s' % s for s in clust_object.modDf.columns],
                loc=10)
    axh.spines['right'].set_color('none')
    axh.spines['left'].set_color('none')
    axh.spines['top'].set_color('none')
    axh.spines['bottom'].set_color('none')
    axh.set_xticks([])
    axh.set_yticks([])
    axh.set_facecolor('white')

    plt.figure(48).savefig(os.path.join(folder, '%s_color_label_legend.png' % clust_object.name), dpi=300)

def plot_module_correl(clust_object, folder):
    """Plot intra-module correlation"""
    i = 0
    for lab in list(cy.labels2modules(clust_object.labels, clust_object.dropped).keys()):
        plt.figure(50+i, figsize=(15, 9))
        cyplot.plotModuleCorr(clust_object.cyDf, clust_object.labels, lab, dropped=clust_object.dropped)
        plt.figure(50+i).savefig(os.path.join(folder, '%s_modules_correlations_%s.png' % (clust_object.name, lab)), dpi=300)
        i += 1

def plot_cy_embedding(clust_object, folder):
    """Cytokine embedding"""
    plt.figure(901, figsize=(13, 9.7))
    cyplot.plotModuleEmbedding(clust_object.dmatDf, clust_object.labels, method='kpca')
    if len(np.unique(clust_object.labels)) > 2:
        colors = palettable.colorbrewer.get_map('Set1', 'qualitative', len(np.unique(clust_object.labels))).mpl_colors
    else:
        colors = [(0.8941176470588236, 0.10196078431372549, 0.10980392156862745),
                   (0.21568627450980393, 0.49411764705882355, 0.7215686274509804)]
    colorLegend(colors, ['%s%1.0f' % (clust_object.sampleStr, i) for i in np.unique(clust_object.labels)],
                loc='lower left')
    plt.figure(901).savefig(os.path.join(folder, '%s_embedding.png' % clust_object.name), dpi=300)
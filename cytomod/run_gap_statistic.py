import gapstat
import cytomod as cy
import matplotlib.pyplot as plt
from functools import partial
import pandas as pd


def getBestK(cyData, b=1000, max_testing_k=12, bootstraps=200, max_final_k=6, save_fig_path=None):
    max_testing_k += 1

    def clusteringFunc(data, K, b=b):
        pwrel, labels, dropped = cy.formReliableClusters(data, cy.corrDmatFunc, partial(cy.hierClusterFunc, K=K), threshold=0, bootstraps=b)
        return labels

    lsICD, mBSICD, errBSICD, gap = gapstat.computeGapStatistic(cyData, cy.corrDmatFunc, clusteringFunc, cy.hierClusterFunc, max_testing_k, bootstraps=bootstraps, clusFuncRealGetsRawData=True)

    stat = pd.DataFrame(gap[:-1] - (gap[1:] - errBSICD[1:]), index=range(1, max_testing_k), columns=['stat'])
    bestK = extractK(stat, max_testing_k=max_testing_k, max_final_k=max_final_k)
    gapstat.plotGapStat(lsICD, mBSICD, errBSICD, gap)
    if save_fig_path is None:
        plt.show()
    else:
        plt.savefig(save_fig_path)

    return(bestK)

def extractK(stat, max_testing_k, max_final_k=6):
    final_k = 0
    for k in range(1, max_testing_k):
        if (stat.loc[k, 'stat'] > 0):
            final_k = k
            break

    if final_k < 2 or final_k > max_final_k:
        final_k = stat.loc[2:6, 'stat'].idxmax()

    return(final_k)
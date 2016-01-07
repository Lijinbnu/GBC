# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:l

"""Compute the connectivity matrix. """



import numpy as np
import nibabel as nib
from scipy import spatial as sp
from scipy import stats

class conn(object):

      """A class to compute the connectivity matrix.

    Attributes
    ----------
    metric : metric to compute the connectivity
    tm : task modulated
    cm: connnectivity matrix

    """
    def __init__(self, metric = 'pearson', tm = False):
        """
        metric: A metric to compute connectivity
        tm: indicate the conn is modulated by task or not, tm = task modulation

        """
        self.metric = metric

        self.tm = tm


    def compute(metric = 'pearson', ):
        if ~tm:
            if self.metrc == 'pearson':
                conn_mat =  stats.pearsonr(ds.tc,'correlation')

            elif self.metric == 'wavelet':
                print 'wavelet'
        else:
            cond_num = ds.cond.shape[1]
            # calculate connectivity
            for c in range(0,cond_num):
                conn_mat = weightedcorr(ds.tc,ds.cond[:,c])


        ds.cm = conn_mat


        return ds































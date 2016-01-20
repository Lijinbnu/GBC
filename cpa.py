import os

import nibabel as nib
import numpy as np
from scipy import stats


class UserDefinedException(Exception):
    """
    Exception defined by user
    """

    def __init__(self, str):
        """

        Parameters
        ----------
        str : a string to indicate the exception

        """
        Exception.__init__(self)
        self._str = str


def pearson_correlation(D, w=None):
    """

    Parameters
    ----------
    D : A 2-D array containing multiple variables and observations.
        Each row of `D` represents a variable, and each column a single
        observation of all those variables.

    w : 1-D array of observation vector weights.

    Returns
    -------
    R : 2-D array, the corrcoef matrix of the variables.


    """
    if w is None:
        R = np.corrcoef(D)
    else:
        c = np.cov(D, aweights=w)
        d = np.diag(c)
        R = c / np.sqrt(np.outer(d, d))

    return R


class DataSet(object):
    def __init__(self, ftarg_img, fnode_img, level='voxel', flabel_img=None, cond_file=None):
        """

        Parameters
        ----------
        ftarg_img : target image file, str
        fnode_img : node image file, str
        flabel_img : label image file for node
        cond_file : condition file, str
        level : level of interest, str, voxel or roi

        """
        # load target image
        targ_img = nib.load(ftarg_img)
        if len(targ_img.get_shape()) == 4:
            targ = targ_img.get_data()
            self.header = targ_img.header
        else:
            raise UserDefinedException('target image is not a 4D Nifti volume!')

        # load node image
        node_img = nib.load(fnode_img)
        if (len(node_img.get_shape()) == 3) and (node_img.get_shape() == targ_img.get_shape()[:3]):
            node = node_img.get_data()
        else:
            raise UserDefinedException('Node image and target image are not match!')

        self.level = level
        # extract tc for voxel
        if self.level == 'voxel':
            imgdim = self.header.get_data_shape()[:3]
            self.nid = node.reshape(np.prod(imgdim))
            self.tc = targ[node.astype(np.bool), :]

        # extract tc for roi
        elif self.level == 'roi':
            nid = np.unique(node)
            self.nid = nid[1:]
            tc = np.zeros((self.nid.shape[0], targ.shape[3]))
            for i in range(0, self.nid.shape[0]):
                tc[i, :] = np.mean(targ[node == i, :], axis=0)
            self.tc = tc
        else:
            self.tc = []
            raise UserDefinedException('Wrong level! it should be voxel or roi.')

        if flabel_img is None:
            self.nlabel = np.ones(self.nid.shape[0])
        else:
            label_img = nib.load(flabel_img)
            label = label_img.get_data()
            if (node_img.get_shape() == label_img.get_shape()) and (node.astype(np.bool) == label.astype(np.bool)).all():
                if self.level == 'voxel':
                    imgdim = self.header.get_data_shape()[:3]
                    self.nlabel = label.reshape(np.prod(imgdim))
                else:
                    for i in range(0, self.nid.shape[0]):
                        self.nlabel[i] = np.mean(label[node == i])
            else:
                raise UserDefinedException('Label image and Node image are not match!')

        # Read design matrix from design file
        if cond_file is None:
            self.cond = None
        else:
            cond = np.loadtxt(cond_file, skiprows=5)
            self.cond = cond[:, np.arange(0, cond.shape[1] - 6, 2)]

    def set_tc(self, tc):
        self.tc = tc

    def set_cond(self, cond):
        self.cond = cond

    def set_nid(self, nid):
        self.nid = nid


class Connectivity(object):
    def __init__(self, ds, metric='pearson', tm=False):
        """

        Parameters
        ----------
        ds : DataSet object
        metric : metric to compute the connectivity
        tm : is task modulated?

        """

        self.metric = metric
        self.tm = tm
        self.ds = ds
        self.mat = []

    def compute(self):
        """

        Parameters
        ----------
        ds : DataSet object

        Returns
        -------
        self : A connecivity object

        """
        ds = self.ds
        if self.metric == 'pearson':
            if not self.tm:
                self.mat = pearson_correlation(ds.tc)
            else:
                self.mat = np.zeros((ds.tc.shape[0], ds.tc.shape[0], ds.cond.shape[1]))
                # standardize  weights to [0,1]
                W = ds.cond
                wmax, wmin = W.max(), W.min()
                W = (W - wmin) / (wmax - wmin)
                for c in range(0, W.shape[1]):
                    self.mat[:, :, c] = pearson_correlation(ds.tc, W[:, c])

        elif self.metric == 'wavelet':
            print 'Wavelet metric does not work now, which is not implemented.'

        return  self

    def set_ds(self, ds):
        self.ds = ds


class Measure(object):
    def __init__(self, conn, metric='sum', ntype='weighted'):
        """

        Parameters
        ----------
        metric: metric to measure the connectivity pattern,str
        thr: threshold, scalar
        ntype: network type, str

        """
        self.metric = metric
        if self.metric == 'sum':
            self.cpu = np.nansum
        elif self.metric == 'std':
            self.cpu = np.nanstd
        elif self.metric == 'skewness':
            self.cpu = stats.skew
        elif self.metric == 'kurtosis':
            self.cpu = stats.kurtosis

        self.conn = conn
        self.ntype = ntype
        self.partition = []
        self.thr = []
        self.value = []

    def compute(self, thr=None, partition=None):
        """

        Parameters
        ----------
        thr : threshod to remove non-interest edge, scalar
        partition : node partition, 1-D array


        Returns
        -------
        self : A Measure object

        """

        self.thr = thr
        if self.thr is None and (self.ntype == 'binary'):
            raise UserDefinedException('Thr is necessary for binary network!')

        if self.thr is None:
            mat = self.conn.mat
        else:
            mat = self.conn.mat > thr
            if self.ntype == 'weighted':
                mat = mat * self.conn.mat

        if partition is None:
            self.partition = self.conn.ds.nlabel
        else:
            self.partition = partition

        P = np.unique(self.partition).tolist()
        for i in P[1:]:
            I = np.where(self.partition[self.partition.astype(np.bool)] == i)
            for j in P[1:]:
                J = np.where(self.partition[self.partition.astype(np.bool)] == j)
                sub_mat = mat[np.ix_(I[0], J[0])]
                self.value.append(self.cpu(sub_mat, axis=1))

        return self

    def set_conn(self, conn):
        self.conn = conn

    def save(self, outdir='.'):
        """

        Parameters
        ----------
        ds : DataSet object which the measure was based on
        outdir : dir to save the measures


        """

        P = np.unique(self.partition).tolist()
        NP = np.count_nonzero(P)
        ds = self.conn.ds
        if ds.level == 'roi':
            # convert self.value to 2D array, every coloumn correspond a seed module
            value = np.zeros((ds.nid.shape[0], NP * NP))
            for i in P:
                I = (self.partition == i)
                for j in P:
                    J = int((i - 1) * NP + (j - 1))
                    value[I, J] = self.value[J]

            np.savetxt(os.path.join(outdir, self.metric), value, fmt= '%.3f')

        elif ds.level == 'voxel':
            imgdim = ds.header.get_data_shape()[:3]
            value = np.zeros((np.prod(imgdim), NP * NP))

            # convert self.value to 4D array, every 3D volume correspond a seed module
            for i in P[1:]:
                I = (self.partition == i)
                for j in P[1:]:
                    J = int((i - 1) * NP + (j - 1))
                    value[I, J] = self.value[J]

            # save voxel-wise inter-module measure in 4D volume
            value = np.reshape(value, (imgdim[0], imgdim[1], imgdim[2], NP * NP))
            header = ds.header
            header['cal_max'] = value.max()
            header['cal_min'] = value.min()
            img = nib.Nifti1Image(value, None, header)
            nib.save(img, os.path.join(outdir, self.metric + '.nii.gz'))

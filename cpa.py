import os
import nibabel as nib
import numpy as np
from scipy import spatial as sp
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

        -------

        """
        Exception.__init__(self)
        self._str = str


def pearson_correlation(D, w = None):
    """

    Parameters
    ----------
    D : A 2-D array containing multiple variables and observations.
        Each column of `D` represents a variable, and each row a single
        observation of all those variables.

    w : 1-D array of observation vector weights.

    Returns

    R : 2-D array, the corrcoef matrix of the variables.
    -------

    """
    if w is None:
        R = np.corrcoef(D)
    else:
        c = np.cov(D, aweights=w)
        d = np.diag(c)
        R = c / np.sqrt(np.outer(d, d))

    return R


class DataSet(object):
    def __init__(self, targ_img_file, node_img_file, level='voxel', cond_file=None):
        """

        Parameters
        ----------
        targ_img_file : target image file, str
        mask_img_file : mask image file, str
        cond_file : condition file, str
        level : level of interest, str, voxel or roi

        -------

        """
        # load target image
        targ_img = nib.load(targ_img_file)
        if len(targ_img.get_shape()) != 4:
            raise UserDefinedException('target image is not a 4D Nifti volume!')
        targ = targ_img.get_data()
        self.header = targ_img.header

        # load node image
        node_img = nib.load(node_img_file)
        if len(node_img.get_shape()) != 3:
            raise UserDefinedException('node image is not a 3D Nifti volume!')
        node = node_img.get_data()

        self.level = level
        # extract tc for voxel
        if self.level == 'voxel':
            self.nid = node[node.astype(np.bool)]
            self.tc = targ[node.astype(np.bool),:]

        # extract tc for roi
        elif self.level == 'roi':
            nid = np.unique(node)
            self.nid = nid[1:]
            tc = np.zeros((self.nid.shape[0],targ.shape[3]))
            for i in range(0,self.nid.shape[0]):
                tc[i,:] = np.mean(targ[node == i,:],axis=0)
            self.tc = tc
        else:
            self.tc = []
            raise UserDefinedException('Wrong level! it should be voxel or roi.')

        # Read design matrix from design file
        if cond_file is not None:
            cond = np.loadtxt(cond_file, skiprows=5)
            self.cond  = cond[:, np.arange(0, cond.shape[1]-6, 2)]

    def set_tc(self, tc):
        self.tc = tc

    def set_cond(self, cond):
        self.cond = cond

    def set_nid(self, nid):
        self.nid = nid

    def set_header(self, header):
       self.header = header


class Connectivity(object):
    def __init__(self, metric='pearson', tm=False):
        """

        Parameters
        ----------
        ds
        metric
        tm


        Returns
        -------

        """

        self.metric = metric
        self.tm = tm
        self.mat = []

    def compute(self, ds):
        """

        Parameters
        ----------
        ds : DataSet object

        Returns
        -------

        """
        if self.metric == 'pearson':
            if not self.tm:
                self.mat = pearson_correlation(ds.tc)
            else:
                self.mat = np.zeros((ds.tc.shape[0],ds.tc.shape[0],ds.cond.shape[1]))
                W = ds.cond
                wmax, wmin = W.max(), W.min()
                W = (W-wmin) / (wmax - wmin)
                for c in range(0, W.shape[1]):
                    self.mat[:, :, c] = pearson_correlation(ds.tc, W[:, c])

        elif self.metric == 'wavelet':
            print 'wavelet metric is not implemented'


class Measure(object):
    def __init__(self, metric='sum', ntype='weighted'):
        """

        Parameters
        ----------
        metric: metric to measure the connectivity pattern,str
        thr: threshold, scalar

        ntype: network type, str

        Returns
        -------

        """
        self.metric = metric
        if self.metric == 'sum':
            self.cpu = np.nansum
        elif self.metric == 'std':
            self.cpu = np.nanstd
        elif self.metric == 'skewness':
            self.cpu = stats.skew
        elif self.metric == 'kurtosis':
            self.cpu = stats.skew

        self.ntype = ntype
        self.partition = []
        self.thr = []
        self.value = []

    def compute(self, conn, thr=None, partition=None):
        """

        Parameters
        ----------
        conn : Connectivity object
        thr : threshod to remove non-interest edge, scalar
        partition : node partition, 1-D array


        """

        self.thr = thr
        if self.thr is None and (self.ntype == 'binary'):
            raise UserDefinedException('Thr is necessary for binary network!')

        if self.thr is None:
            mat = conn.mat
        else:
            mat = conn.mat > thr

        if self.ntype == 'weighted':
            mat = mat * conn.mat

        if partition is None:
            self.partition = np.ones(conn.mat.shape[0])
        else:
            self.partition = partition

        P = np.unique(self.partition).tolist()
        for i in P:
           I = np.where(self.partition == i)
           for j in P:
               J = np.where(self.partition == j)
               sub_mat = mat[np.ix_(I[0], J[0])]
               self.value.append(self.cpu(sub_mat,axis=1))

    def save(self, ds, outdir='.'):
        """

        Parameters
        ----------
        ds : DataSet object which the measure was based on
        outdir : dir to save the measures

        -------

        """

        P = np.unique(self.partition).tolist()
        NP = len(P)
        if ds.level == 'roi':
            # convert self.value to 2D array, every coloumn correspond a seed module
            value = np.zeros((ds.nid.shape[0], NP * NP))
            for i in P:
                I = (self.partition == i)
                for j in P:
                    J = int((i-1) * NP + (j - 1))
                    value[I, J] = self.value[J]

            np.savetxt(os.path.join(outdir,self.metric),value)

        elif ds.level == 'voxel':
            imgdim = ds.header.get_data_shape()
            value = np.zeros((np.prod(imgdim), NP * NP))

            # convert self.value to 4D array, every 3D volume correspond a seed module
            for i in P:
                I = (self.partition == i)
                for j in P:
                    J =  int((i-1) * NP + (j - 1))
                    value[I, J] = self.value[J]

            # save voxel-wise inter-module measure in 4D volume
            value = np.reshape(value, (imgdim[0], imgdim[1], imgdim[2], NP*NP))
            header = ds.header
            header['cal_max'] = value.max()
            header['cal_min'] = value.min()
            img = nib.Nifti1Image(value, None, header)
            nib.save(img, os.path.join(outdir,self.metric+'.nii.gz'))

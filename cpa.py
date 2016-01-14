import os
import nibabel as nib
import numpy as np
from scipy import spatial as sp


class UserDefinedException(Exception):
    def __init__(self, str):
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
    print D
    if w is None:
        R = np.corrcoef(D)
    else:
        c = np.cov(D, aweights=w)
        d = np.diag(c)
        #print d
        R = c / np.sqrt(np.outer(d, d))

    return R


class DataSet(object):
    def __init__(self, targ_img_file, node_img_file, level='voxel', cond_file=None):
        """

        Parameters
        ----------
        targ_img_file
        node_img_file
        level
        cond_file

        Returns
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
        if self.level == 'voxel':
            self.nid = node[node.astype(np.bool)]
            self.tc = targ[node.astype(np.bool),:]

        elif self.level == 'roi':
            nid = np.unique(node)
            self.nid = nid[1:]
            tc = np.zeros(len(self.nid),targ_img.get_shape()[3])
            for i in range(0,len(self.nid)):
                tc[i,:] = np.mean(targ[node == i,:])
            self.tc = tc
        else:
            self.tc = []
            raise UserDefinedException('Wrong level! it should be voxel or roi.')

        # Read design matrix from design.mat file
        if cond_file is not None:
            cond = np.loadtxt(cond_file, skiprows=5)
            self.cond  = cond[:, np.arange(0, cond.shape[1]-6, 2)]
            print self.cond.shape

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
#
# class Measure(object):
#     def __init__(self, ds, conn, metric = 'sum', ntype = 'weighted'):
#
#         self.metric = metric
#         mat = conn.mat
#         mat[np.diag_indices_from(mat)] = np.nan
#
#         if self.metric == 'sum':
#             self.cpu = np.nansum(mat, axis=1)
#         elif self.metric == 'std':
#             self.cpu = np.nanstd(mat, axis=1)
#         elif self.metric == 'skewness':
#             masked = np.ma.masked_array(mat, np.isnan(mat))
#             self.cpu = stats.skew(masked, axis=1)
#         elif self.metric == 'kurtosis':
#             masked = np.ma.masked_array(mat, np.isnan(mat))
#             self.cpu = stats.skew(masked, axis=1)
#
#         self.ntype = ntype
#         self.thr = []
#
#     def compute(self, ds, conn, thr=None):
#
#         self.thr = thr
#         self.value = []
#         self.index = []
#
#         if self.thr is None:
#             mat = conn.mat
#         else:
#             if self.ntype == 'binary':
#                 mat = conn.mat >= threshold
#               #  mat[mat == 0] = np.nan
#             else:
#                 mat = conn.mat*(conn.mat >= threshold)
#                 mat[mat == 0] = np.nan
#
#         for i in np.unique(ds.module):
#            i_index = np.asarray(np.nonzero(ds.module == i)).T
#            for j in np.unique(ds.module):
#                j_index = np.asarray(np.nonzero(ds.module == j)).T
#                sub_mat = np.zeros((i_index.shape[0],j_index.shape[1]),dtype=float)
#                sub_mat = mat[i_index].reshape(i_index.shape[0],-1)[:,j_index].reshape(i_index.shape[0],-1)
#                if self.metric == 'sum':
#                    cpu = np.nansum(sub_mat, axis=1)
#                elif self.metric == 'std' and self.ntype == 'weighted':
#                    cpu = np.nanstd(sub_mat, axis=1)
#                elif self.metric == 'skewness' and self.ntype == 'weighted':
#                    masked = np.ma.masked_array(sub_mat, np.isnan(sub_mat))
#                    cpu = stats.skew(masked, axis=1)
#                elif self.metric == 'kurtosis' and self.ntype == 'weighted':
#                    masked = np.ma.masked_array(sub_mat, np.isnan(sub_mat))
#                    cpu = stats.skew(masked, axis=1)
#
#                self.value.append(cpu)
#                self.index.append([i,j])
#
#
# class CPA(object):
#     def __init__(self, ds, conn, meas):
#         self.ds = ds
#         self.conn = conn
#         self.meas = meas
#
#     def comp_conn(self):
#         self.conn.compute(self.ds)
#
#     def meas_conn(self, thr=None):
#         self.meas.compute(self.ds, self.conn, thr)
#
#     def save(self, outdir):
#         if self.ds.level == 'roi':
#             for i in range(0, len(self.meas.value)):
#                 index = self.meas.index[i]
#                 np.savetxt(outdir+'/'+'roi'+'-'+str(index[0])+'to'+str(index[1]), self.meas.value[i])
#
#         elif self.ds.level == 'voxel':
#             node_img = nib.load(self.ds.fnode_img)
#             node = node_img.get_data()
#             cell = np.zeros((node.shape[0],node.shape[1],node.shape[2]))
#             for i in range(0, len(self.meas.value)):
#                 index = self.meas.index[i]
#                 cell[(node == index[0]).astype(np.bool)] = self.meas.value[i]
#                 img = nib.Nifti1Image(cell, self.ds.affine)
#                 img.to_filename(os.path.join(outdir,self.meas.metric+'-'+str(index[0])+'to'+str(index[1])+'.nii.gz'))


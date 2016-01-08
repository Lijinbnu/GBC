# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:l

"""Extract time series for connectivity computation. """

import os

import nibabel as nib
import numpy as np
from scipy import stats


def load_files(in_files):
    vol_file = in_files[0]
    mask_file = in_files[1]
    vol = nib.load(vol_file)
    mask = nib.load(mask_file)

    return vol,mask



class user_defined_exception(Exception):
    """User defined exception.
    """
    def __init__(self, str):
        Exception.__init__(self)
        self._str = str



class cpa(object):
    """
    Connectivity pattern analysis(cpa) class

    Attributes
    ----------
    ds : data set to compute the connectivity
    conn : connectivity pattern
    meas: measures for the connectivity pattern
    """
    def __init__(self, targ_img_file, mask_img_file,cond_file = None):

        """
        Parameters
        ----------
        targ_img_file: target image file
        mask_img_file: mask image file
        cond_file: condtion file

        Returns
        -------

        """


        self.file.targ_img = targ_img_file
        self.file.mask_img = mask_img_file
        self.file.cond = cond_file


    def extract_ds(self, level = 'voxel'):
        """

        Parameters
        ----------
        level: node level, i.e., voxel or level

        Returns
        -------

        """


        self.ds.level = level

        # load target image
        targ_img = nib.load(self.file.targ_img)
        if len(targ_img.get_shape()) != 4:
            raise user_defined_exception('targ_img is not a Nifti image about 4D volume!')
        targ = targ_img.get_data()

        # load  mask image
        mask_img = nib.load(self.file.mask_img)
        if len(mask_img.get_shape()) != 3:
            raise user_defined_exception('mask_img is not a Nifti image about 3D volume!')
        mask = mask_img.get_data()

        self.ds.header = targ_img.header
        self.ds.coord = np.nonzero(mask)

        if self.level == 'voxel':
            self.ds.tc = targ[mask.astype(np.bool), :]
            # label for each non-zeros voxel
            self.ds.label = mask[mask.astype(np.bool), :]

        elif self.level == 'roi':
            label  = np.unique(mask)
            # label for each ROI
            self.label = label[1:]
            # compute ROI time course
            tc = np.zeros(targ_img.get_shape()[3], len(self.label))
            for i in self.label:
                tc[:,i] =  np.mean(targ[mask[mask.astype(np.bool)] == i,:])
            self.ds.tc = tc
        else:
            print 'wrong level'


        # read cond file and assigh cond array to ds.cond
        self.ds.cond = cond

    def comp_conn(self,metric = 'pearson', tm = False):
        """

        Parameters
        ----------
        metric: metric to compute connnectivity,str
        tm: task modulation, bool value

        Returns
        -------

        """
        self.conn.metric =  metric
        self.conn.tm = tm
        if not tm:
            if self.metrc == 'pearson':
                mat =  stats.pearsonr(self.ds.tc,'correlation')

            elif self.metric == 'wavelet':
                print 'wavelet'
        else:
            cond_num = self.ds.cond.shape[1]
            # calculate weighted correlation for each condition
            for c in range(0,cond_num):
                mat = weighted_corr(self.ds.tc,self.ds.cond[:,c])


        self.conn.mat = mat


    def meas_conn(self,metric = 'sum', thr = None, module = None, ntype = 'weighted'):
        """

        Parameters
        ----------
        metric: metric to measure the connectivity pattern,str
        thr: threshold, scalar
        module: module assignment, 1-D array
        type: network type, str

        Returns
        -------

        """
        self.meas.metric = metric
        if self.meas.metric == 'sum':
            comp_meas = np.nansum
        elif self.meas.metric == 'std':
            comp_meas = np.std
        elif self.meas.metric == 'skewness':
            comp_meas = stats.skew
        elif self.meas.metric == 'kurtosis':
            comp_meas= stats.kurtosis

        self.meas.thr = thr
        if thr is not None:
            mat = self.conn.mat > thr

        self.meas.ntype = type
        if ntype == 'weighted':
            mat = np.multiply(mat,self.conn.mat)

        # module ID for each node
        if module is None:
            module = np.ones(self.conn.mat.shape[0])
        self.meas.module = module;


        M = np.unique(self.meas.module)
        self.meas.value = []

        for i in M:
            I = np.asarray(self.meas.module == i).T
            for j in M:
                J = np.asarray(self.meas.module == j)
                sub_mat = mat[I].reshape(I.shape[0],-1)[:,J].reshape(I.shape[0],-1)

        self.meas.value.append(comp_meas(sub_mat))

    def save_conn(self, outdir):
        """

        Parameters
        ----------
        outdir: dir to save the connectivity measures

        Returns
        -------

        """
        if self.ds.level == 'roi':
        # save meas.value as txt

        elif self.ds.level == 'voxel':
            # reshape meas.value and save it as nii
            nib.save(self.meas,os.path.join(outdir,self.meas.metric+'.nii.gz'))



    def meas_anat(self):
        print 'anat measure'





    def set_tc(self, tc):
        self.tc = tc

    def set_cond(self,cond):
        self.cond = cond

    def set_coord(self, coord):
        self.coord = coord

    def  set_label(self, label):
        self.label = label

    def set_header(self,header):
        self.header = header

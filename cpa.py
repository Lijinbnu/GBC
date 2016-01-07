# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:l

"""Extract time series for connectivity computation. """

import nibabel as nib
import numpy as np
from scipy import spatial as sp
import numpy as np


def load_files(in_files):
    vol_file = in_files[0]
    mask_file = in_files[1]
    vol = nib.load(vol_file)
    mask = nib.load(mask_file)

    return vol,mask

def save_Nifti_files(imgs, filepath):
    i = 0
    for index in imgs:
        nib.save(index, filepath+ '/'+str(i)+'.nii.gz')
        i = i+1



class user_defined_exception(Exception):
    """User defined exception.
    """
    def __init__(self, str):
        Exception.__init__(self)
        self._str = str



class cpa(object):
    """Storage class for datasets with multiple attributes.

    A dataset consists of three pieces: data matrix, coords and cond.


    Attributes
    ----------
    tc : Time courses for all nodes(NxM array). N is the number of samples(e.g., TR),
        M is the number of nodes(e.g., roi)

    coord : Spatial coordinates for all nodes(3xM array)

    label: label ID for each node(Mx1 array)

    cond : Design info for each sample (e.g., TR), NxC array. C is the number of conditions

    """
    def __init__(self, targ_img_file, mask_img_file,cond_file = None):
        """
        targ_img: A 4D Nifti images object
        mask_img: A 3D Nifti images object
        level: indicate the level of connectivity, i.e., roi or voxel
        cond: NxC array
        """

        self.file.targ_img = targ_img_file

        self.file.mask_img = mask_img_file

        self.file.cond = cond_file


    def extract_ds(self, level = 'voxel'):
        # node level
        self.ds.evel = level

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
        self.conn.metric =  metric
        self.conn.tm = tm
        if not tm:
            if self.metrc == 'pearson':
                conn_mat =  stats.pearsonr(ds.tc,'correlation')

            elif self.metric == 'wavelet':
                print 'wavelet'
        else:
            cond_num = self.ds.cond.shape[1]
            # calculate connectivity
            for c in range(0,cond_num):
                conn_mat = weightedcorr(self.ds.tc,self.ds.cond[:,c])


        self.conn.cm = conn_mat


    def meas_conn(self,metric = ['sum'], thr = None, graph_type = 'weighted'):
        self.meas.metric = metric
        self.meas.thr = thr
        self.meas.graph_type = graph_type

        if thr is not None:
            cm = self.conn.cm > thr

        if graph_type == 'weighted':
            cm = cm.* self.conn.cm

         if metric == 'sum':
             corr_strength = np.nansum(cm, axis=1)
         elif index == 'std':
             corr_strength = np.std(cm, axis=1)
         elif index == 'skewness':
             corr_strength = stats.skew(cm, axis=1, bias=False)
         elif index == 'kurtosis':
             corr_strength = stats.kurtosis(cm, axis=1, bias=False)

        return self.conn.corr_str

    def save_conn(self, map):



    def meas_anat(self):






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

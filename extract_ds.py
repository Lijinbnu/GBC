# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:l

"""Extract time series for connectivity computation. """

import nibabel as nib
import numpy as np
from scipy import spatial as sp
import numpy as np




class user_defined_exception(Exception):
    """User defined exception.
    """
    def __init__(self, str):
        Exception.__init__(self)
        self._str = str



class dataset(object):
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
    def __init__(self, targ_img, mask_img, cond=None):
        """
        targ_img: A 4D Nifti images object
        mask_img: A 3D Nifti images object
        cond: NxC array
        """
        if len(targ_img.get_shape()) != 4:
            raise user_defined_exception('targ_img is not a Nifti image about 4D volume!')
        if len(mask_img.get_shape()) != 3:
            raise user_defined_exception('mask_img is not a Nifti image about 3D volume!')

        vol = targ_img.get_data()
        mask = mask_img.get_data()

        self.tc = vol[mask.astype(np.bool), :]

        self.label = mask[mask.astype(np.bool), :]


        self.coord = np.nonzero(mask) * targ_img.get_affine()

        self.cond = cond



    def set_tc(self, tc):
        self.tc = tc

    def set_cond(self,cond):
        self.cond = cond

    def set_coord(self, coord):
        self.coord = coord

    def  set_label(self, label):
        self.label = label

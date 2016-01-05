# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:l

"I/O for brain connectivity pattern analysis"

import nibabel as nib
import numpy as np
from scipy import spatial as sp
import numpy as np
import os
from rehoneib import reho_volneighbors
from statsmodels.stats.weightstats import DescrStatsW # Old version of this class has bug !!!

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

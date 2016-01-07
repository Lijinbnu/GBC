# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:l

"I/O for brain connectivity pattern analysis"

import cpa as cpa
import numpy as np
from scipy import spatial as sp


targ_img_file = '/nfs/j3/userhome/zhenzonglei/workingdir/bn/S0001/mem/002/func.feat/stanard_filtered_func_data.nii.gz'
mask_img_file = ''
cond_file = ''


mycpa = cpa(targ_img_file,mask_img_file,cond_file)
mycpa.extraxt_ds()
mycpa.comp_conn()
mycpa.meas_conn()
mycpa.save()

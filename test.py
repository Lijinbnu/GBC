# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:l

"Test for connectivity pattern analysis"

from cpa import *
import numpy as np
from scipy import spatial as sp


targ_img_file = '/nfs/j3/userhome/zhenzonglei/workingdir/bn/S0001/mem/002/func.feat/stanard_filtered_func_data.nii.gz'
mask_img_file = ''
cond_file = ''

outdir = ''

# define dataset
ds = DataSet(targ_img_file, mask_img_file,cond_file = None)
ds.load()

# define connectivity
conn = Connectivity('pearson', False)

# define measures
meas = Measure('sum')

# define cpa pipeline
cpa = CPA(ds,conn,meas)

# run cpa to compute conn
cpa.comp_conn()

# run cpa to measure conn
cpa.meas_conn()

# save connS
cpa.save(outdir)

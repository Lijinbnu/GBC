# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:l

"Test for connectivity pattern analysis"

from cpa import *


targ_img_file = './data/S0001_obj_004.nii.gz'
node_img_file = './data/face.nii.gz'
cond_file = './data/design.mat'


# define dataset
ds = DataSet(targ_img_file, node_img_file,'roi',cond_file)

# define connectivity
conn = Connectivity('pearson')
conn.compute(ds)


# define measures
meas = Measure('sum')
meas.compute(conn)

# outdir = ''
meas.save(ds)

# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:l

"Test for connectivity pattern analysis"

from cpa import *


targ_img_file = './data/S0001_obj_004.nii.gz'
node_img_file = './data/face.nii.gz'
label_img_file = None
cond_file = './data/design.mat'


# define dataset
ds = DataSet(ftarg_img=targ_img_file, fnode_img=node_img_file,
             flabel_img=label_img_file, level='voxel')

# define and compute connectivity
conn = Connectivity(ds, metric='pearson').compute()


# define and compute global measures
glob_meas = Measure(conn, metric='iqr').compute()
glob_meas.save()


# define and compute local measures
# local_meas = LocalMeasure(conn, radius=6, metric='sum').compute()
# local_meas.save()


# define and compute spatial measures
#spat_meas = SpatialMeasure(conn).compute()
#spat_meas.save()

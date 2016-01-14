import os
import nibabel as nib
import numpy as np
from scipy import spatial as sp
from scipy import stats
from statsmodels.stats.weightstats import DescrStatsW

class user_defined_exception(Exception):

    def __init__(self, str):
        Exception.__init__(self)
        self._str = str

# define each sub structures
#
##class file(object):
##   def __init__(self):
#        pass
#
#class ts(object):
#    def __init__(self):
#        pass
#
#class conn(object):
#    def __init__(self):
#        pass
#



class DataSet(object):
    
    def __init__(self, targ_img_file, node_img_file, cond_file=None):

        self.ftarg_img = targ_img_file
        self.fnode_img = node_img_file
        self.fcond = cond_file
        self.tc = []
        self.affine = []
        self.label = []

    def load(self, level='voxel', module_file=None):

        self.fmodule = module_file

        # load target image
        targ_img = nib.load(self.ftarg_img)
        if len(targ_img.get_shape()) != 4:
            raise user_defined_exception('targ_img is not a Nifti image about 4D volume!')
        targ = targ_img.get_data()
        self.affine = targ_img.get_affine()

        # load node image
        node_img = nib.load(self.fnode_img)
        if len(node_img.get_shape()) != 3:
            raise user_defined_exception('node_img is not a Nifti image about 3D volume!')
        node = node_img.get_data()

        self.level = level
        if self.level == 'voxel':
            self.label = node[node.astype(np.bool)]
            if self.fmodule is not None:
                module_img = nib.load(self.fmodule)
                module = module_img.get_data()
                self.module = module[module.astype(np.bool)]
            else:
                self.module = node[node.astype(np.bool)]
            
            self.tc = targ[node.astype(np.bool),:]

        elif self.level == 'roi':
            label = np.unique(node)
            self.label = label[1:]
            if self.fmodule is not None:
                self.module = np.loadtxt(self.fmodule)
            else:
                self.module = np.ones(self.label[0])

            tc = np.zeros(self.label,targ_img.get_shape()[3])
            for i in range(0,len(self.label)):
                tc[i,:] = np.mean(targ[node[node.astype(np.bool)] == i,:])
            self.tc = tc
        else:
            print 'wrong level'
       
        #Read design matrix from design.mat file
        if self.fcond is not None:
            self.cond = np.loadtxt(self.fcond, skiprows=5)

    def set_module(self, module_file):
        self.fmodule = module_file
        if self.level == 'voxel':
            module_img = nib.load(self.fmodule)
            module = module_img.get_data()
            self.module = module[module.astype(np.bool)]
        elif self.level == 'roi':
            self.module = np.loadtxt(self.fmodule)





        #self.cond = 'cond'
    
    #def set_tc(self, tc):
    #    self.tc = tc

    #def set_cond(self, cond):
    #    self.cond = cond

    #def set_label(self, label):
    #    self.label = label

    #def set_header(self, header):
    #    self.header = header
     

class Connectivity(object):
    def __init__(self, ds, metric='pearson', tm=False):
        self.metric = metric
        self.tm = tm
        self.mat = np.zeros((ds.label.shape[0], ds.label.shape[0]))

    def compute(self, ds):
        if not self.tm:
            if self.metric == 'pearson':
                corr = 1 - sp.distance.pdist(ds.tc, 'correlation')
                self.mat = sp.distance.squareform(corr)
            elif self.metric == 'wavelet':
                pass
        else:
            if ds.cond.ndim == 1:
                cond_num = 1
            else:
                cond_num = ds.cond.shape[1]
            for condition in range(0, cond_num):
                weights = cond[:,condition]
                weights = weights - np.min(weights)
                weights = weights/np.sum(weights)
                d1 = DescrStatsW(ds.tc.T, weights=weights)
                cov = np.dot(weights*d1.demeaned.T, d1.demeaned)
                corr = cov/d1.std/d1.std[:,None]
                self.mat = corr - np.eye(corr.shape[0], dtype=float)




class Measure(object):
    def __init__(self, ds, conn, metric='sum', ntype='weighted'):
        
        self.metric = metric
        mat = conn.mat
        mat[np.diag_indices_from(mat)] = np.nan
        
        if self.metric == 'sum':
            self.cpu = np.nansum(mat, axis=1)
        elif self.metric == 'std':
            self.cpu = np.nanstd(mat, axis=1)
        elif self.metric == 'skewness':
            masked = np.ma.masked_array(mat, np.isnan(mat))
            self.cpu = stats.skew(masked, axis=1)
        elif self.metric == 'kurtosis':
            masked = np.ma.masked_array(mat, np.isnan(mat))
            self.cpu = stats.skew(masked, axis=1)

        self.ntype = ntype
        self.thr = []

    def compute(self, ds, conn, thr=None):

        self.thr = thr
        self.value = []
        self.index = []

        if self.thr is None and (self.ntype == 'binary'):
            raise user_defined_exception('you should set threshold for binary image!')
        elif self.thr is None and (self.ntype == 'weighted'):
            mat = conn.mat
        else:
            if self.ntype == 'binary':
                mat = conn.mat >= self.thr
              #  mat[mat == 0] = np.nan
            else:
                mat = conn.mat*(conn.mat >= self.thr)
                mat[mat == 0] = np.nan
        
        for i in np.unique(ds.module):
           i_index = np.asarray(np.nonzero(ds.module == i)).T
           for j in np.unique(ds.module):
               j_index = np.asarray(np.nonzero(ds.module == j)).T
               sub_mat = np.zeros((i_index.shape[0],j_index.shape[1]),dtype=float)
               sub_mat = mat[i_index].reshape(i_index.shape[0],-1)[:,j_index].reshape(i_index.shape[0],-1)
               if self.metric == 'sum':
                   cpu = np.nansum(sub_mat, axis=1)
               elif self.metric == 'std' and self.ntype == 'weighted':
                   cpu = np.nanstd(sub_mat, axis=1)
               elif self.metric == 'skewness' and self.ntype == 'weighted':
                   masked = np.ma.masked_array(sub_mat, np.isnan(sub_mat))
                   cpu = stats.skew(masked, axis=1)
               elif self.metric == 'kurtosis' and self.ntype == 'weighted':
                   masked = np.ma.masked_array(sub_mat, np.isnan(sub_mat))
                   cpu = stats.skew(masked, axis=1)

               self.value.append(cpu)
               self.index.append([i,j])


class CPA(object):
    def __init__(self, ds, conn, meas):
        self.ds = ds
        self.conn = conn
        self.meas = meas

    def comp_conn(self):
        self.conn.compute(self.ds)

    def meas_conn(self, thr=None):
        self.meas.compute(self.ds, self.conn, thr)

    def save(self, outdir):
        if self.ds.level == 'roi':
            for i in range(0, len(self.meas.value)):
                index = self.meas.index[i]
                np.savetxt(outdir+'/'+'roi'+'-'+str(index[0])+'to'+str(index[1]), self.meas.value[i])

        elif self.ds.level == 'voxel':
            if self.ds.fmodule is not None:
                module_img = nib.load(self.ds.fmodule)
                module = module_img.get_data()
            else:
                module_img = nib.load(self.ds.fnode_img)
                module = module_img.get_data()
            cell = np.zeros((module.shape[0],module.shape[1],module.shape[2]))
            for i in range(0, len(self.meas.value)):
                index = self.meas.index[i]
                cell[(module == index[0]).astype(np.bool)] = self.meas.value[i]
                img = nib.Nifti1Image(cell, self.ds.affine)
                img.to_filename(os.path.join(outdir,self.meas.metric+'-'+str(index[0])+'to'+str(index[1])+'.nii.gz'))


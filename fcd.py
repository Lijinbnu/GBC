# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:l

__author__ = 'zhouguangfu,wangxu'

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


## weighted type 1
#def weighted(ts,weights):
#    """Get weighted ts.
#
#    Contributions
#    -------------
#        Author:wx
#        Editor:
#    parameters
#    ----------
#    ts: input time seriers, N(voxels) x M(timepoints).
#    weights: the weights of each timepoint, usually the row is design matrix, one row for one condition, 
#             1 x M(timepoints).
#
#    """
#    demeants = ts - np.mean(ts)
#    weightedts = weights * demeants
#    return weightedts


def weightedcorr(data1,data2,weights):
    """Get weighted correlation coefficients of data1 and data2.

    Contributions
    -------------
        Author:wx
        Editor:
    parameters
    ----------
    data1, data2: array_like, input time seriers, N(voxels) x M(timepoints).
    weights: the weights of each timepoint, usually the row is design matrix, one row for one condition, 
             1 x M(timepoints).

    """
    s = data1.shape
    sB = data2.shape

    if len(s) != 2:
        raise ValueError('data1 must be a 2-dimensional array.');
    if len(sB) != 2:
        raise ValueError('data1 must be a 2-dimensional array.');
    if s[1] != sB[1]:
        raise ValueError('data1 and data2 must have the same number of columns (i.e. number of time points.)')

    data = np.concatenate((data1.T, data2.T), axis=1)
    d1 = DescrStatsW(data, weights=weights)
    cov = np.dot(weights*d1.demeaned.T,d1.demeaned)
    weightedcorr = cov/d1.std/d1.std[:,None]
    wcorrdata1_data2 = weightedcorr[0:s[0],:][:,s[0]:]

    return wcorrdata1_data2


class user_defined_exception(Exception):
    """User defined exception.

    Contributions
    -------------
        Author:zhouguangfu
        Editor:

    """
    def __init__(self, str):
        Exception.__init__(self)
        self._str = str

def voxfcd(vol,mask,thr, type='single',cond=None,radius='global',shape='fast_cube',tail='right', graphtype='binary'):
    """
    :param vol: struct of Nifti images about 4D volume
    :param mask: struct of Nifti images about about 3D mask volume
    :param type: 'single' or 'multiple', default is 'single'
    :param thr: a vector of different threshold values
    :param cond: a design matrix for different conditions, default is None
    :param radius: only be used when calculate lfcd, default is 'global' which means gfcd
    :param shape: if we compute lfcd, the shape of the local region, default is 'fast_cube'
    :param tail: 'both' , 'left' or 'right', default is 'both'
    :param graphtype: 'binary' or 'strength', default is 'binary'
    :return: a list of Nifti images
    """
    pfcd = []
    nfcd = []
    fcdstrength = []

    try:
        if len(vol.get_shape()) != 4:
            raise user_defined_exception('Vol is not a Nifti image about 4D volume!')
        if len(mask.get_shape()) != 3:
            raise user_defined_exception('Mask is not a Nifti image about 3D volume!')
        if (type != 'single') and (type != 'multiple'):
            raise user_defined_exception("Type error!The value of type must be 'single' or 'multiple'")
        if thr.shape[0] == 0:
            raise user_defined_exception('The thr is empty!')
        vollen = vol.get_shape()[3]
        if cond == None:
            print 'Unweighted correlation now!'
        elif cond.shape[0] != vollen:
            print vollen,cond.shape[0]
            raise user_defined_exception('The dimension of the condition(design matrix) is mismatch!')
        if radius != 'global':
            if (radius <= 0) or (radius > vollen):
                raise user_defined_exception('Radius is too small or too big!')
        if shape == 'fast_cube':
            if radius not in [6,18,26]:
                raise user_defined_exception('The radius \''+radius+'\' not in the fast_cube list!')
        if shape not in [None,'fast_cube','cube','sphere']:
            raise user_defined_exception("Shape error! The value of shape must be 'fast_cube' or 'cube' or 'sphere'! ")
        if (tail != 'left') and (tail!= 'right') and (tail!='both'):
            raise user_defined_exception("Tail error! The value of tail must be 'both' or 'left' or 'right!' ")
            print 'Tail error!'
        if (graphtype != 'binary') and (graphtype!= 'strength'):
            raise user_defined_exception("Graphtype error! The valuse of tagraphtype  must be 'bianry' or 'strength'")
            print 'Graphtype error!'
    except user_defined_exception,str:
        print 'user_defined_exception : ',str._str
        import sys
        sys.exit(-1)

    vol_data = vol.get_data()
    mask_data = mask.get_data()
    affine = vol.get_affine()
    mask_num = len(np.unique(mask_data)) - 1
    #cond_num = len(np.unique(cond))
    if cond == None:
        cond_num = 1
    else:
        if len(cond.shape) == 1:
            cond_num = 1
        else:
            cond_num = cond.shape[1]
    thr_num = len(thr)
    print 'Mask_num is  ', mask_num

    if radius == 'global':
        print 'global fcd----------------------------------------------'
        if type =='single' or (type=='multiple' and mask_num == 1):
            maskts = vol_data[mask_data.astype(np.bool), :]
            #no matter the type is 'single' or 'multiple', the result is the same.
            for condition in range(0,cond_num):
                if cond != None:
                    if cond_num == 1:
                        weights = cond
                    else:
                        weights = cond[:,condition]
                    # We first transform the weights with min of zero, then we norm the 
                    # weights to have the sum of one. 
                    weights = weights - np.min(weights)
                    weights = weights/np.sum(weights)
                    d1 = DescrStatsW(maskts.T, weights=weights)
                    # Compute the weighted covariance of the data
                    cov = np.dot(weights*d1.demeaned.T,d1.demeaned)
                    corr = cov/d1.std/d1.std[:,None]
                    # Transform dignonal value of 1 to 0.
                    corr = corr - np.eye(corr.shape[0],dtype=float)
                else:
                    corr = 1 - sp.distance.pdist(maskts, 'correlation')
                    # Transform the vector to 2d array, diagonal value is 0.
                    corr = sp.distance.squareform(corr)
                    
                # Default: save raw strength map no matter 'binary' or 'strength' graphtype.
                # Use nansum: since the corr has 'nan' value, why?
                corr_strength = np.nansum(corr, axis=1)
                mask_copy_data_strength = mask_data.copy()
                mask_copy_data_strength[mask_copy_data_strength.astype(np.bool)] = corr_strength
                fcdstrength.append(nib.Nifti1Image(mask_copy_data_strength, affine))
                
                # Save binary or strength map after different threshold.  
                for threshold in thr:
                    if graphtype == 'binary':
                        corr_p = np.sum(corr >= threshold, axis=1)
                        corr_n = np.sum(corr <= -threshold, axis=1)
                    else:
                        corr_p = np.nansum(corr*(corr >= threshold), axis=1)
                        corr_n = np.nansum(corr*(corr <= -threshold), axis=1)
                    if tail == 'both':
                        mask_copy_data_p = mask_data.copy()
                        mask_copy_data_p[mask_copy_data_p.astype(np.bool)] = corr_p
                        pfcd.append(nib.Nifti1Image(mask_copy_data_p, affine))
                        mask_copy_data_n = mask_data.copy()
                        mask_copy_data_n[mask_copy_data_n.astype(np.bool)] = corr_n
                        nfcd.append(nib.Nifti1Image(mask_copy_data_n, affine))
                    elif tail == 'left':
                        mask_copy_data = mask_data.copy()
                        mask_copy_data[mask_copy_data.astype(np.bool)] = corr_n
                        nfcd.append(nib.Nifti1Image(mask_copy_data, affine))
                    else:
                        mask_copy_data = mask_data.copy()
                        mask_copy_data[mask_copy_data.astype(np.bool)] = corr_p
                        pfcd.append(nib.Nifti1Image(mask_copy_data, affine))

        elif type=='multiple':
            print 'Multi mask sub mode'
            mask_subs = np.unique(mask_data)
            #np.delete(mask_subs,0)
            mask_subs = mask_subs[1:]
            pcells = np.zeros((mask_data.shape[0],mask_data.shape[1],mask_data.shape[2],mask_num*mask_num))
            ncells = np.zeros((mask_data.shape[0],mask_data.shape[1],mask_data.shape[2],mask_num*mask_num))
            scells = np.zeros((mask_data.shape[0],mask_data.shape[1],mask_data.shape[2],mask_num*mask_num))
            maskts = vol_data[mask_data.astype(np.bool), :]                
            for condition in range(0,cond_num):
                if cond != None:
                    if cond_num == 1:
                        weights = cond
                    else:
                        weights = cond[:,condition]
                    weights = weights - np.min(weights)
                    weights = weights/np.sum(weights)
                    d1 = DescrStatsW(maskts.T, weights=weights)
                    cov = np.dot(weights*d1.demeaned.T,d1.demeaned)
                    corr = cov/d1.std/d1.std[:,None]
                    # Transform dignonal value of 1 to 0.
                    corr = corr - np.eye(corr.shape[0],dtype=float)
                else:
                    corr = 1 - sp.distance.pdist(maskts, 'correlation')
                    corr = sp.distance.squareform(corr)

                # Default: save raw strength map no matter 'binary' or 'strength' 
                # graphtype.
                strength_index = 0
                for i in mask_subs:
                    i_index = np.asarray(np.nonzero(mask_data[mask_data.astype(np.bool)] == i)).T
                    for j in mask_subs:
                        j_index = np.asarray(np.nonzero(mask_data[mask_data.astype(np.bool)] == j))
                        result = np.zeros((i_index.shape[0],j_index.shape[1]),dtype=float)
                        result = corr[i_index].reshape(i_index.shape[0],-1)[:,j_index].reshape(i_index.shape[0],-1)
                        corr_strength = np.sum(result, axis=1)
                        scells[(mask_data == i).astype(np.bool),strength_index]= corr_strength
                        strength_index += 1
                fcdstrength.append(nib.Nifti1Image(scells, affine))
                        
                for threshold in thr:
                    index = 0
                    for i in mask_subs:
                        i_index = np.asarray(np.nonzero(mask_data[mask_data.astype(np.bool)] == i)).T
                        for j in mask_subs:
                            j_index = np.asarray(np.nonzero(mask_data[mask_data.astype(np.bool)] == j))
                            result = np.zeros((i_index.shape[0],j_index.shape[1]),dtype=float)
                            result = corr[i_index].reshape(i_index.shape[0],-1)[:,j_index].reshape(i_index.shape[0],-1)
                            if graphtype == 'binary':
                                corr_p = np.sum(result >= threshold, axis=1)
                                corr_n = np.sum(result <= -threshold, axis=1)
                            else:
                                corr_p = np.nansum(result*(result >= threshold), axis=1)
                                corr_n = np.nansum(result*(result <= -threshold), axis=1)
                                """pcells contains the value which is bigger than positive threshold, each value is a 3D volume data. 
                                ncells contains the value which is smaller than negative threshold, each value is a 3D volume data.
                                The rule of storing the values is just one 3D volume by 3D volume per threshold per condition. For 
                                example, if the mask volume contains just two masks named A and B per threshold per condition, the 
                                arrangement of the data will be AA,AB,BA,BB."""
                            pcells[(mask_data == i).astype(np.bool),index]= corr_p
                            ncells[(mask_data == i).astype(np.bool),index] = corr_n
                            index += 1
 
                    if tail == 'both':
                        pfcd.append(nib.Nifti1Image(pcells, affine))
                        nfcd.append(nib.Nifti1Image(ncells, affine))
                    elif tail == 'left':
                        nfcd.append(nib.Nifti1Image(ncells, affine))
                    else:
                        pfcd.append(nib.Nifti1Image(pcells, affine))

    else:
        print 'local fcd--------------------------------------------------'
        timepoints = vol_data.shape[3]
        dims = vol_data.shape[0:3]
        # We calculate the local FCD for each voxel
        vol_data = np.reshape(vol_data,[dims[0]*dims[1]*dims[2],timepoints])
        maskhdr = mask.get_header()
        res = maskhdr['pixdim'][1:4]
        # neighbor.ncut_volneighbors(imgdata, ndim(2 or 3), size (radius of
        # the sphere), shape='sphere')
        neib = reho_volneighbors(mask_data, 3, radius, res, shape)
        neib = neib.compute_offsets()
               
        # We calculate the local FCD for each voxel
        mask_copy_data = mask_data.copy()
        mask_copy_data_s = mask_data.copy()
        for condition in range(0,cond_num):
            for element in neib:
                targts = vol_data[element[0], :]
                allts = vol_data[element[1], :]
                # correlation matrix, each row is a seed voxel
                if cond != None:
                    if cond_num == 1:
                        weights = cond
                    else:
                        weights = cond[:,condition]
                    weights = weights - np.min(weights)
                    weights = weights/np.sum(weights)
                    corr = weightedcorr(targts,allts,weights)
                else:
                    corr = 1 - sp.distance.cdist(targts, allts, 'correlation')
                corr_strength = np.nansum(corr, axis=1)
                # Exculde the self correlation by minus 1
                mask_copy_data_s[np.unravel_index(element[0], dims)] = corr_strength - 1
                fcdstrength.append(nib.Nifti1Image(mask_copy_data_s, affine))

            for threshold in thr:
                for element in neib:
                    targts = vol_data[element[0], :]
                    allts = vol_data[element[1], :]
                    # correlation matrix, each row is a seed voxel
                    if cond != None:
                        if cond_num == 1:
                            weights = cond
                        else:
                            weights = cond[:,condition]
                        weights = weights - np.min(weights)
                        weights = weights/np.sum(weights)
                        corr = weightedcorr(targts,allts,weights)
                    else:
                        corr = 1 - sp.distance.cdist(targts, allts, 'correlation')

                    if graphtype == 'binary':
                        corr_p = np.sum(corr >= threshold, axis=1) - 1
                        corr_n = np.sum(corr <= -threshold, axis=1)
                    else:
                        corr_p = np.sum(corr*(corr >= threshold), axis=1) - 1
                        corr_n = np.sum(corr*(corr <= -threshold), axis=1)
                    mask_data[np.unravel_index(element[0], dims)] = corr_p
                    mask_copy_data[np.unravel_index(element[0], dims)] = corr_n

                if tail == 'both':
                    pfcd.append(nib.Nifti1Image(mask_data, affine))
                    nfcd.append(nib.Nifti1Image(mask_copy_data, affine))
                elif tail == 'left':
                    nfcd.append(nib.Nifti1Image(mask_copy_data, affine))
                else:
                    pfcd.append(nib.Nifti1Image(mask_data, affine))

    return pfcd, nfcd, fcdstrength



#import time
#print time.clock()

#fsessid = open('sessid')
#subject_list  = [line.strip() for line in fsessid]
#for subject_id in subject_list:
   # vol,mask= load_files(['/nfs/j3/userhome/wangxu/workingdir/nitk/nitk-pipeline/trunk/nkpi/rfmri/sess/fcdtestdata/'+subject_id+'/obj/002/filtered_func_data_std.nii.gz','/nfs/j3/userhome/wangxu/workingdir/nitk/nitk-pipeline/trunk/nkpi/rfmri/sess/fcdtestdata/'+subject_id+'/obj/002/face_obj_label.nii.gz'])
#    vol,mask= load_files(['/nfs/j3/userhome/wangxu/workingdir/nitk/nitk-pipeline/trunk/nkpi/rfmri/sess/fcdtestdata/'+subject_id+'/rest/002/bp0.01_0.1confrm.nii.gz','/nfs/j3/userhome/wangxu/workingdir/nitk/nitk-pipeline/trunk/nkpi/rfmri/sess/fcdtestdata/'+subject_id+'/rest/002/faceact_invlin.nii.gz'])

    #Read design matrix from design.mat file
#    dm = np.loadtxt('fcdtestdata/'+subject_id+'/obj/002/design.mat',skiprows=5)
    #Choose first, third and fifth column as input weights
#    dm = dm[:,[0,2,4]]
#    print subject_id
#    thr = np.zeros((1,))
 #   thr[:] = 0.6
#    radius = 10
#    pfcd, nfcd, fcdstrength = voxfcd(vol,mask,thr,type='single',cond=None,radius='global',shape=None,tail='both',graphtype='binary')

#    print time.clock()

#    datadir = '/nfs/j3/userhome/wangxu/workingdir/nitk/nitk-pipeline/trunk/nkpi/rfmri/sess/fcdtestdata/'
#    if pfcd != []:
#        cur = datadir+'/'+subject_id+'/pfcd'
#        if not os.path.exists(cur):
#            os.mkdir(cur)
#        save_Nifti_files(pfcd,cur)
#    if nfcd != []:
#        cur = datadir+'/'+subject_id+'/nfcd'
#        if not os.path.exists(cur):
#            os.mkdir(cur)
#        save_Nifti_files(nfcd,cur)
#    if fcdstrength != []:
#        cur = datadir+'/'+subject_id+'/fcdthrength'
#        if not os.path.exists(cur):
#            os.mkdir(cur)
#        save_Nifti_files(fcdstrength,cur)

#print time.clock()
#print 'program finished!'




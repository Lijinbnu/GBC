# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:l

"""Compute the connectivity matrix. """



import numpy as np
import nibabel as nib
from scipy import spatial as sp
from scipy import stats
from statsmodels.stats.weightstats import DescrStatsW

def load_files(in_files):
    vol_file = in_files[0]
    mask_file = in_files[1]
    vol = nib.load(vol_file)
    mask = nib.load(mask_file)

    return vol,mask



def connectivity(maskts, method='pearson',cond=None):

    if cond == None:
        cond_num = 1
    else:
        cond_num = cond.shape[1]
    # calculate connectivity
    for condition in range(0,cond_num):
        if cond != None:
            if cond_num == 1:
                weights = cond
            else:
                weigths = cond[:,condition]
            # we first transform the weights with min of zero, then we
            # transform the weights to have the sum of one.
            weights = weights - np.min(weights)
            weights = weights/np.sum(weights)
            d1 = DescrStatsW(maskts.T, weights=weights)
            # Compute the weighted covariance of the data
            cov = np.dot(weights*d1.demeaned.T,d1.demeaned)
            corr = cov/d1.std/d1.std[:,None]
            # Transform dignonal value of 1 to 0.
            corr = corr - np.eye(corr.shape[0],dtype=float)
        else:
            if method == 'pearson':
                corr = 1 - sp.distance.pdist(maskts,'correlation')
                # Transform dignonal value of 1 to 0. 
                corr = sp.distance.squareform(corr)
            elif method == 'wavelet':
                .....

    return corr


def distance(coord):
    # Compute pairwise eucilidean distance between selected voxels
    distan = 1 - sp.distance.pdist(maskts,'euclidean')
    return distan

def metrics(datatype = 'corr', mask, thr, index='sum', tail='right', graphtype='binary'):
    
    mask_data = mask.get_data()
    mask_num = len(np.unique(mask_data)) - 1
    
    if datatype == 'corr':
        # get the diagnal offf
        corr_digoff = corr.reshape(corr.shape[1]**2,)
        corr_digoff = corr_digoff[np.where(corr_digoff != 0)].reshape(corr_digoff.shape[1],corr_digoff.shape[1]-1)
        
        if mask_num == 1:
            print 'Single network mode'
            if index == 'sum':
                corr_strength = np.nansum(corr_digoff, axis=1)
            elif index == 'std':
                corr_strength = np.std(corr_digoff, axis=1)
            elif index == 'skewness':
                corr_strength = stats.skew(corr_digoff, axis=1, bias=False)
            elif index == 'kurtosis':
                corr_strength = stats.kurtosis(corr_digoff, axis=1, bias=False)
            mask_copy_data_strength = mask_data.copy()
            mask_copy_data_strength[mask_copy_data_strength.astype(np.bool)] = corr_strength
            fcdstrength.append(nib.Nifti1Image(mask_copy_data_strength, affine))
            
            for threshold in thr:
                if graphtype == 'binary':
                    # binary only produce sum metrics
                    corr_p = np.sum(corr >= threshold, axis=1)
                    corr_n = np.sum(corr <= -threshold, axis=1)
                else:
                    if index == 'sum':
                        corr_p = np.nansum(corr*(corr >= threshold), axis=1)
                        corr_n = np.nansum(corr*(corr <= -threshold), axis=1)
                    elif index  == 'std':
                        corr_p = np.std(corr*(corr >= threshold), axis=1)
                        corr_n = np.std(corr*(corr <= -threshold), axis=1)
                    elif index == 'skewness':
                        corr_p = stats.skew(corr*(corr >= threshold), axis=1)
                        corr_n = stats.skew(corr*(corr <= -threshold), axis=1)
                    elif index == 'kurtosis':
                        corr_p = stats.kurtosis(corr*(corr >= threshold), axis=1)
                        corr_n = stats.kurtosis(corr*(corr <= -threshold), axis=1)
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
        else:
            print 'Multiple network mode'
            mask_subs = np.unique(mask_data)
            mask_subs = mask_subs[1:]
            pcells = np.zeros((mask_data.shape[0],mask_data.shape[1],mask_data.shape[2],mask_num*mask_num))
            ncells = np.zeros((mask_data.shape[0],mask_data.shape[1],mask_data.shape[2],mask_num*mask_num))
            scells = np.zeros((mask_data.shape[0],mask_data.shape[1],mask_data.shape[2],mask_num*mask_num))
            strength_index = 0
            for i in mask_subs:
                i_index = np.asarray(np.nonzero(mask_data[mask_data.astype(np.bool)] == i)).T
                for j in mask_subs:
                    j_index = np.asarray(np.nonzero(mask_data[mask_data.astype(np.bool)] == j))
                    result = np.zeros((i_index.shape[0],j_index.shape[1]),dtype=float)
                    result = corr[i_index].reshape(i_index.shape[0],-1)[:,j_index].reshape(i_index.shape[0],-1)
                    # remove diagnal within a network
                    if i == j:
                         result = result.reshape(result.shape[1]**2,)
                         result = result[np.where(result != 0)].reshape(result.shape[1],result.shape[1]-1)
                    # calculate metrics
                    if index == 'sum':
                        corr_strength = np.nansum(result, axis=1)
                    elif index == 'std':
                        corr_strength = np.std(result, axis=1)
                    elif index == 'skewness':
                        corr_strength = stats.skew(result, axis=1, bias=False)
                    elif index == 'kurtosis':
                        corr_strength = stats.kurtosis(result, axis=1, bias=False)

                    scells[(mask_data == i).astype(np.bool),strength_index] = corr_strength
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
                                    if index == 'sum':
                                        corr_p = np.nansum(result*(result >= threshold), axis=1)
                                        corr_n = np.nansum(result*(result <= -threshold), axis=1)
                                    elif index  == 'std':
                                        corr_p = np.std(result*(result >= threshold), axis=1)
                                        corr_n = np.std(result*(result <= -threshold), axis=1)
                                    elif index == 'skewness':
                                        corr_p = stats.skew(result*(result >= threshold), axis=1)
                                        corr_n = stats.skew(result*(result <= -threshold), axis=1)
                                    elif index == 'kurtosis':
                                        corr_p = stats.kurtosis(result*(result >= threshold), axis=1)
                                        corr_n = stats.kurtosis(result*(result <= -threshold), axis=1)

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
    
    elif datatype == 'distance':



























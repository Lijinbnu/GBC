#!/usr/bin/env python
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:l

import os                                    # system functions
import argparse
import nipype.interfaces.io as nio           # Data i/o
import nipype.interfaces.fsl as fsl          # fsl
import nipype.interfaces.utility as util     # utility
import nipype.pipeline.engine as pe          # pypeline engine
import time

# Setup any package specific configuration. The output file format
# for FSL routines is being set to uncompressed NIFTI

print fsl.Info.version()
fsl.FSLCommand.set_default_output_type('NIFTI_GZ')

def subjrlf(subject_id, data_dir, fsd, rlf, fstem, mask):
    """
    Get input file information.
    """
    import os
    frlf = open(os.path.join(data_dir,subject_id,fsd,rlf))
    run_list = [line.strip() for line in frlf]
    info = dict(prefcd=[[subject_id,fsd,run_list,fstem]],
                masks=[[subject_id,fsd,run_list,mask]],)
    return info

def fcdf(in_files,shape,size,threshold):
    """
    Compute voxel-wise fcd in the cube with edge of size
    """
    import os
    import nibabel as nib
    import numpy as np
    #from nkpi.rfmri.sess import load_files,save_Nifti_files,voxfcd 
    from nkpi.rfmri.sess.fcd import load_files,save_Nifti_files,voxfcd
     
    if not isinstance(in_files, list):
        # For single run's processing
        in_files = [in_files]

    for in_file in in_files:
        mask = in_file[0]
        func = in_file[1]
        vol,mask= load_files([func,mask])
        thr = threshold
        pfcd, nfcd, fcdstrength = voxfcd(vol,mask,thr,type='single',cond=None,radius='global',shape=None,tail='both',graphtype='strength')
        fcddir = os.path.dirname(in_file[0])
        if pfcd != []:
            cur = fcddir+'/pfcdrsts'
            if not os.path.exists(cur):
                os.mkdir(cur)
            save_Nifti_files(pfcd,cur)
        if nfcd != []:
            cur = fcddir+'/nfcdrsts'
            if not os.path.exists(cur):
                os.mkdir(cur)
            save_Nifti_files(nfcd,cur)
        if fcdstrength != []:
            cur = fcddir+'/strengthrsts'
            if not os.path.exists(cur):
                os.mkdir(cur)
            save_Nifti_files(fcdstrength,cur)

    return cur

def main():
    """
    usage: lfcd-sess [-h] (-datadir datadir | -datadirf datadir-file)
                     (-sess sessid | -sessf sessid-file) -fsd func-subdir -rlf rlf
                     -fstem func-file [-radius RADIUS]
                     [-plugin {Linear,Multiproc,IPython}] [-debug] [-v]

    Do ReHo analysis.

    Parameters
    ----------
      -h, --help            show this help message and exit
      -datadir datadir      Source directory contains data file
      -datadirf datadir-file
                            File contains the source data directory
      -sess sessid          Input the sessid
      -sessf sessid-file    Input the sessid file
      -fsd func-subdir      Functional sub directory, e.g. bold
      -rlf rlf              Run list file
      -fstem func-file      The file name(suffix) of the functional image
      -shape {fast_cube,cube,sphere}
                            The shape of the neighbour voxels
      -size size/radius     The radius of the sphere to compute Kendal W, voxel
                            based
      -plugin {Linear,Multiproc,IPython}
                            The name of the plugin, the available plugins allow
                            local and distributed execution of workflows, default
                            is IPython
      -debug                Debug mode, save mediate results in present dir
      -v, --version         show program's version number and exit

    Examples
    --------
    Specify the sphere with the radius of 3 as the neighbours of the voxel
    and compute the lfcd value for all the voxels in the sphere:
    lfcd-sess -datadirf sesspar -sess S0001 -fsd rest -rlf rfMRI.rlf 
              -fstem confrm -shape sphere -size 3 

    """
    parser = argparse.ArgumentParser(prog='vfcd-sess',
                                     prefix_chars='-',
                                     description='Do voxel-fcd analysis.')
    group1 = parser.add_mutually_exclusive_group(required=True)
    group1.add_argument('-datadir',
                        help='Source directory contains data file',
                        metavar='datadir',
                        dest='datadir')
    group1.add_argument('-datadirf',
                        help='File contains the source data directory',
                        metavar='datadir-file',
                        dest='datadirf')
    group2 = parser.add_mutually_exclusive_group(required=True)
    group2.add_argument('-sess',
                        help='Input the sessid',
                        metavar='sessid',
                        dest='sess')
    group2.add_argument('-sessf',
                        help='Input the sessid file',
                        metavar='sessid-file',
                        dest='sessf')
    parser.add_argument('-fsd',
                        help='Functional sub directory, e.g. bold',
                        dest='fsd',
                        metavar='func-subdir',
                        required=True)
    parser.add_argument('-rlf',
                        help='Run list file',
                        dest='rlf',
                        metavar='rlf',
                        required=True)
    parser.add_argument('-fstem',
                        help='The file name(suffix) of the functional image',
                        dest='fstem',
                        metavar='func-file',
                        required=True)
    parser.add_argument('-mask',
                        help='Seed mask file',
                        dest='mask',
                        required=True)
    parser.add_argument('-shape',
                        help='The shape of the neighbour voxels',
                        dest='shape',
                        choices=['fast_cube','cube','sphere'],
                        default='fast_cube')
    parser.add_argument('-size',
                        help='The radius of the sphere to compute Kendal W, voxel based',
                        dest='size',
                        metavar='size/radius',
                        type=int,
                        default='26')
    parser.add_argument('-thr',
                        help='Correlation coefficient threshold',
                        dest='thr',
                        nargs='*',
                        type=float,
                        default=[0.6])
    parser.add_argument('-plugin',
                        help='The name of the plugin, the available plugins '
                              'allow local and distributed execution of '
                              'workflows, default is IPython',
                        dest='plugin',
                        default = 'IPython',
                        choices=['Linear','Multiproc','IPython'])
    parser.add_argument('-debug',
                        help='Debug mode, save mediate results in present dir',
                        dest='debug',
                        default = False,
                        action='store_true')
    parser.add_argument('-v','--version',
                        action='version',
                        version='%(prog)s 0.1')

    args = parser.parse_args()

    # Parallel computation exec config

    pluginName = args.plugin
    # Specify the location of the data

    fsessid = args.sessf
    sessid = args.sess
    if fsessid:
        fsessid = open(fsessid)
        subject_list  = [line.strip() for line in fsessid]
    elif sessid:
        subject_list = [sessid]

    datadir = args.datadir
    datadirf = args.datadirf
    if datadir:
        data_dir = datadir
    elif datadirf:
        datadirf = open(datadirf)
        data_dir = datadirf.readline().strip()

    if args.debug:
        targetdir = './'
    elif not args.debug:
        targetdir = ''

    fsd = args.fsd
    rlf = args.rlf
    fstem = args.fstem
    mask = args.mask
    shape = args.shape
    size = args.size
    import numpy as np
    threshold = np.array(args.thr)
    

    # Set up complete workflow

    vfcd = pe.Workflow(name='vfcd')

    infosource = pe.Node(interface=util.IdentityInterface(fields=['subject_id']),
                        name="infosource")
    infosource.iterables = ('subject_id', subject_list)

    datasource = pe.Node(interface=nio.DataGrabber(infields=['subject_id'],
                                                   outfields=['prefcd','masks']),
                         name = 'datasource')
    datasource.inputs.base_directory = data_dir
    datasource.inputs.template = '%s/%s/%s/%s.nii.gz'
    datasource.inputs.sort_filelist = False

    mergenode = pe.Node(interface=util.Merge(2, axis='hstack'),
                        name='merge')

    voxfcd = pe.Node(interface=util.Function(input_names=['in_files','shape','size','threshold'],
                                             output_names=['out_file'],
                                             function=fcdf),
                     name='voxfcd')
    voxfcd.inputs.shape = shape
    voxfcd.inputs.size = size
    voxfcd.inputs.threshold = threshold

    vfcd.base_dir = os.path.abspath(targetdir)
    vfcd.connect([(infosource, datasource, [('subject_id', 'subject_id'),
                                            (('subject_id',subjrlf,data_dir,fsd,rlf,fstem,mask),'template_args')]),
                  (datasource, mergenode, [('masks', 'in1')]),
                  (datasource, mergenode,[('prefcd', 'in2')]),
                  (mergenode, voxfcd,[('out','in_files')]),
                ])

    vfcd.run(plugin=pluginName)

if __name__ == '__main__':
    main()








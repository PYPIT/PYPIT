#!/usr/bin/env python
#
# See top-level LICENSE file for Copyright information
#
# -*- coding: utf-8 -*-
"""
This script forces an extraction of undetected objects at the
expected location from the slitmask design. Run above the Science/ folder.
"""
import os
from glob import glob

import numpy as np

from IPython import embed

from astropy.io import fits
from astropy.table import Table
from astropy.stats import sigma_clipped_stats

from pypeit import msgs
from pypeit import slittrace
from pypeit import specobjs
from pypeit import spec2dobj
from pypeit import calibrations
from pypeit.images import pypeitimage
from pypeit import reduce
from pypeit import io


from configobj import ConfigObj
from pypeit.par.util import parse_pypeit_file
from pypeit.par import PypeItPar
from pypeit.spectrographs.util import load_spectrograph

from pypeit.display import display
from pypeit.core.parse import get_dnum
from pypeit.images.imagebitmask import ImageBitMask
from pypeit import masterframe



def parse_args(options=None, return_parser=False):
    import argparse
    parser = argparse.ArgumentParser(description='forces an extraction of undetected objects at the '
                                                 'expected location from the slitmask design '
                                                 '(currently only for keck_deimos).  Run above the Science/ folder',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('pypeit_file', type=str, help='PypeIt reduction file (must have .pypeit extension)')
    parser.add_argument('--science_dir', type = str, default=None, help = 'Directory where spec2d and spec1d files are. '
                                                                  'Default Science/')
    # parser.add_argument('--list', default=False, help='List the extensions only?',
    #                     action='store_true')
    # parser.add_argument('--det', default=1, type=int, help='Detector number')
    # parser.add_argument('--showmask', default=False, help='Overplot masked pixels',
    #                     action='store_true')
    # parser.add_argument('--removetrace', default=False, help="Do not overplot traces in the skysub, "
    #                                                          "sky_resid and resid channels",
    #                     action = "store_true")
    # parser.add_argument('--embed', default=False, help='Upon completion embed in ipython shell',
    #                     action='store_true')
    # parser.add_argument('--ignore_extract_mask', default=False, help='Ignore the extraction mask',
    #                     action='store_true')
    # parser.add_argument("--sensfunc", type=str, default=None, help="Pass in a sensfunc to display the sky-subtracted image with a flux calibration")
    # parser.add_argument('--channels', type=str, help='Only show a subset of the channels (0-indexed), e.g. 1,3')


    if return_parser:
        return parser

    return parser.parse_args() if options is None else parser.parse_args(options)


def read_pypeitfile(pypeit_file):
    cfg_lines, data_files, frametype, usrdata, setups = parse_pypeit_file(pypeit_file, runtime=False)
    cfg = ConfigObj(cfg_lines)
    spectrograph_name = cfg['rdx']['spectrograph']
    spectrograph = load_spectrograph(spectrograph_name)
    config_specific_file = None
    for idx, row in enumerate(usrdata):
        if ('science' in row['frametype']) or ('standard' in row['frametype']):
            config_specific_file = data_files[idx]
    spectrograph_cfg_lines = spectrograph.config_specific_par(config_specific_file).to_config()
    par = PypeItPar.from_cfg_lines(cfg_lines=spectrograph_cfg_lines, merge_with=cfg_lines)
    return par



def main(args):

    if not os.path.isfile(args.pypeit_file):
        msgs.error('PypeIt file not found')

    # Set the Science directory
    if args.science_dir is None:
        science_dir = './Science/'
    else:
        science_dir = args.science_dir

    # Check that the Science directory exists
    if not os.path.exists(science_dir):
        msgs.error('Science directory not found. Provide it with --science_dir')
    else:
        # grab th spec2d files
        spec2d_files = glob(os.path.join(science_dir, 'spec2d*.fits'))

    # find the corresponding spec1d files
    spec1d_files = []
    for spec2dfile in spec2d_files:
        spec1dfile = spec2dfile.replace('spec2d', 'spec1d')
        if not os.path.exists(spec1dfile):
            spec2d_files.remove(spec2dfile)
            msgs.warn('No spec1d file found for {spec2d_file}')
            continue
        spec1d_files.append(spec1dfile)

    # if no files return
    if (len(spec2d_files) == 0) or (len(spec1d_files) == 0):
        msgs.warn('No spec2d or spec1d file found')
        return

    # parameters
    par = read_pypeitfile(args.pypeit_file)
    spectrograph = load_spectrograph(par['rdx']['spectrograph'])
    caliBrate = calibrations.Calibrations(None, par['calibrations'], spectrograph, None)
    embed()
    for s in range(len(spec2d_files)):
        # Load spec2d
        allspec2D = spec2dobj.AllSpec2DObj.from_fits(spec2d_files[s], chk_version=True)
        allspec1D = specobjs.SpecObjs.from_fitsfile(spec1d_files[s], chk_version=True)
        for det in allspec2D.detectors:
            spec2Dobj = allspec2D[det]
            spec1Dobjs = allspec1D[allspec1D.DET == det]
            slits = spec2Dobj.slits
            slits_left, slits_right, _ = slits.select_edges()
            slits.maskdef_designtab = Table(fits.getdata(spec2d_files[s],
                                                         extname='DET{:02d}-MASKDEF_DESIGNTAB'.format(det)))
            new_spec1Dobjs = slits.mask_add_missing_obj(spec1Dobjs, par['reduce']['findobj']['find_fwhm'],
                                                        slits.maskdef_offset, slits_left, slits_right)

            sciImage = pypeitimage.PypeItImage(image=spec2Dobj.sciimg, ivar=spec2Dobj.ivarraw, rn2img=None, bpm=spec2Dobj.bpmmask,
                         crmask=None, fullmask=None, detector=spec2Dobj.detector, spat_flexure=spec2Dobj.sci_spat_flexure, PYP_SPEC=par['rdx']['spectrograph'], imgbitm=spec2Dobj.imgbitm)
            sciImage.build_mask(slitmask=slits.slit_img())
            redux=reduce.Reduce.get_instance(sciImage, spectrograph, par, caliBrate, 'science', ir_redux=False, find_negative=False, det=det, show=False)

            # # Do we have any positive objects to proceed with?
            # if len(new_spec1Dobjs) > 0:
            #     # Global sky subtraction second pass. Uses skymask from object finding
            #     if (std_redux or par['reduce']['extraction']['skip_optimal'] or
            #             par['reduce']['findobj']['skip_second_find'] or usersky):
            #         redux.global_sky = redux.initial_sky.copy()
            #     else:
            #         redux.global_sky = redux.global_skysub(skymask=self.skymask, show=self.reduce_show)
            #
            #     # Apply a global flexure correction to each slit
            #     # provided it's not a standard star
            #     if par['flexure']['spec_method'] != 'skip' and not std_redux:
            #         redux.spec_flexure_correct(mode='global')
            #
            #     # Extract + Return
            #     skymodel, objmodel, ivarmodel, outmask, sobjs \
            #         = redux.extract(redux.global_sky, new_spec1Dobjs)
            #     if redux.find_negative:
            #         self.sobjs.make_neg_pos() if return_negative else self.sobjs.purge_neg()
            # else:  # No objects, pass back what we have
            #     # Apply a global flexure correction to each slit
            #     # provided it's not a standard star
            #     if par['flexure']['spec_method'] != 'skip' and not std_redux:
            #         redux.spec_flexure_correct(mode='global')
            #     #Could have negative objects but no positive objects so purge them
            #     if self.find_negative:
            #         self.sobjs_obj.make_neg_pos() if return_negative else self.sobjs_obj.purge_neg()
            #     self.skymodel = self.initial_sky
            #     self.objmodel = np.zeros_like(self.sciImg.image)
            #     # Set to sciivar. Could create a model but what is the point?
            #     self.ivarmodel = np.copy(self.sciImg.ivar)
            #     # Set to the initial mask in case no objects were found
            #     self.outmask = self.sciImg.fullmask
            #     # empty specobjs object from object finding
            #     self.sobjs = self.sobjs_obj
            #
            # # If a global spectral flexure has been applied to all slits, store this correction as metadata in each specobj
            # if self.par['flexure']['spec_method'] != 'skip' and not self.std_redux:
            #     for iobj in range(self.sobjs.nobj):
            #         islit = self.slits.spatid_to_zero(self.sobjs[iobj].SLITID)
            #         self.sobjs[iobj].update_flex_shift(self.slitshift[islit], flex_type='global')
            #
            # # Correct for local spectral flexure
            # if self.sobjs.nobj == 0:
            #     msgs.warn('No objects to extract!')
            # elif self.par['flexure']['spec_method'] not in ['skip', 'slitcen'] and not self.std_redux:
            #     # Apply a refined estimate of the flexure to objects, and then apply reference frame correction to objects
            #     self.spec_flexure_correct(mode='local', sobjs=self.sobjs)
            #
            # # Apply a reference frame correction to each object and the waveimg
            # self.refframe_correct(ra, dec, obstime, sobjs=self.sobjs)
            #
            # # Update the mask
            # reduce_masked = np.where(np.invert(self.reduce_bpm_init) & self.reduce_bpm)[0]
            # if len(reduce_masked) > 0:
            #     self.slits.mask[reduce_masked] = self.slits.bitmask.turn_on(
            #         self.slits.mask[reduce_masked], 'BADREDUCE')


def entry_point():
    main(parse_args())


if __name__ == '__main__':
    entry_point()

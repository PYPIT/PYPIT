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
from pypeit import io

from pypeit.display import display
from pypeit.core.parse import get_dnum
from pypeit.images.imagebitmask import ImageBitMask
from pypeit import masterframe
from pypeit import spec2dobj


def parse_args(options=None, return_parser=False):
    import argparse
    parser = argparse.ArgumentParser(description='forces an extraction of undetected objects at the '
                                                 'expected location from the slitmask design '
                                                 '(currently only for keck_deimos).  Run above the Science/ folder',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

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


def main(args):

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
            new_spec1Dobjs = slits.mask_add_missing_obj(spec1Dobjs, 5.0, slits.maskdef_offset, slits_left, slits_right)




def entry_point():
    main(parse_args())


if __name__ == '__main__':
    entry_point()

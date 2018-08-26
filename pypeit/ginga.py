""" Module for ginga routines.  Mainly for debugging
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

try:
    basestring
except NameError:
    basestring = str

import os
import numpy as np

# CANNOT LOAD DEBUGGER AS THIS MODULE IS CALLED BY ARDEBUG
#from pypeit import ardebug as debugger
import pdb as debugger

from ginga.util import grc

from astropy.io import fits

from pypeit import msgs

# TODO: There needs to be a way to call show_image() without importing
# arlris, requires a code refactor
# from pypeit import arlris

def connect_to_ginga(host='localhost', port=9000, raise_err=False):
    """ Connect to an active RC Ginga
    Parameters
    ----------
    host : str, optional
    port : int, optional

    Returns
    -------
    viewer : RemoteClient
      connectoin to Ginga

    """
    # Start
    viewer = grc.RemoteClient(host, port)
    # Test
    ginga = viewer.shell()
    try:
        tmp = ginga.get_current_workspace()
    except:
        if raise_err:
            raise ValueError
        else:
            msgs.warn("Problem connecting to Ginga.  Launch an RC Ginga viewer: ginga --module=RC   then continue.")
            debugger.set_trace()
    # Return
    return viewer


def show_image(inp, chname='Image', wcs_img=None, **kwargs):
    """ Displays input image in Ginga viewer
    Supersedes method in xastropy

    Parameters
    ----------
    inp : str or ndarray (2D)
      If str, assumes the image is written to disk
    wcs_img : str, optional
      If included, use this in WCS.  Mainly to show wavelength array

    Returns
    -------

    """
    if isinstance(inp, basestring):
        if '.fits' in inp:
            if 'raw_lris' in kwargs.keys():
#                img, head, _ = arlris.read_lris(inp)
                raise NotImplementedError('ginga.show_image() cannot yet show lris images.')
            else:
                hdu = fits.open(inp)
                try:
                    exten = kwargs['exten']
                except KeyError:
                    exten = 0
                img = hdu[exten].data
    else:
        img = inp

    viewer = connect_to_ginga()
    ch = viewer.channel(chname)
    # Header
    header = {}
    header['NAXIS1'] = img.shape[1]
    header['NAXIS2'] = img.shape[0]
    if wcs_img is not None:
        header['WCS-XIMG'] = wcs_img
        #header['WCS-XIMG'] = '/home/xavier/REDUX/Keck/LRIS/2017mar20/lris_red_setup_C/MF_lris_red/MasterWave_C_02_aa.fits'
    # Giddy up
    ch.load_np(chname, img, 'fits', header)
    return viewer, ch


def show_slits(viewer, ch, lord_in, rord_in, slit_ids = None, rotate=False, pstep=1):
    """ Overplot slits on image in Ginga
    Parameters
    ----------
    viewer
    ch
    lord_in : ndarray
    rord_in : ndarray
    slit_ids : list of int
    rotate : bool, optional
      Allow for a rotated image
    pstep : int
      Show every pstep point of the edges
    """

    # This allows the input lord and rord to either be (nspec, nslit) arrays or a single
    # vectors of size (nspec)
    if lord_in.ndim == 2:
        nslit = lord_in.shape[1]
        lordloc = lord_in
        rordloc = rord_in
    else:
        nslit = 1
        lordloc = lord_in.reshape(lord_in.size,1)
        rordloc = rord_in.reshape(rord_in.size,1)

    if slit_ids is None:
        slit_ids = [str(slit) for slit in np.arange(nslit)]

    # Deal with case that slit_ids is input as a scalar
    if hasattr(slit_ids,"__len__") == False:
        slit_ids = [slit_ids]

    # Canvas
    canvas = viewer.canvas(ch._chname)
    canvas.clear()
    # y-axis
    y = (np.arange(lordloc.shape[0])).tolist()
    #ohf = lordloc.shape[0] // 2
    tthrd = int(2*lordloc.shape[0]/3.)
    # Loop on slits
    for slit in range(lordloc.shape[1]):
        # Edges
        for kk,item in enumerate([lordloc, rordloc]):
            if kk == 0:
                clr = str('green')
            else:
                clr = str('red')
            if rotate:
                points = list(zip(y[::pstep],item[::pstep,slit].tolist()))
            else:
                points = list(zip(item[::pstep,slit].tolist(),y[::pstep]))
            canvas.add(str('path'), points, color=clr)
        # Text -- Should use the 'real' name
        if rotate:
            xt, yt = float(y[tthrd]), float(lordloc[tthrd,slit])
        else:
            xt, yt = float(lordloc[tthrd,slit]), float(y[tthrd])
        canvas.add(str('text'), xt, yt, str('S{:}'.format(slit_ids[slit])), color=str('red'),
                   fontsize=20.)

def show_trace(viewer, ch, trace, trc_name = 'Trace', color='blue', clear=False,
               rotate=False, pstep=1):
    """
    rotate : bool, optional
      Allow for a rotated image
    pstep : int
      Show every pstep point of the edges
    """
    # Canvas
    canvas = viewer.canvas(ch._chname)
    if clear:
        canvas.clear()
    # Show
    y = (np.arange(trace.size)[::pstep]).tolist()
    xy = [trace[::pstep].tolist(), y]
    if rotate:
        xy[0], xy[1] = xy[1], xy[0]
    points = list(zip(xy[0], xy[1]))
    canvas.add(str('path'), points, color=str(color))
    # Text
    ohf = trace.size // (2*pstep)
    xyt = [float(trace[ohf]), float(y[ohf])]
    if rotate:
        xyt[0], xyt[1] = xyt[1], xyt[0]
    canvas.add(str('text'), xyt[0], xyt[1], trc_name, rot_deg=90., color=str(color), fontsize=17.)


def clear_canvas(cname):
    viewer = connect_to_ginga()
    ch = viewer.channel(cname)
    canvas = viewer.canvas(ch._chname)
    canvas.clear()


def chk_arc_tilts(msarc, trcdict, sedges=None, yoff=0., xoff=0., all_green=False, pstep=10,
                  cname='ArcTilts'):
    """  Display arc image and overlay the arcline tilt measurements
    Parameters
    ----------
    msarc : ndarray
    trcdict : dict
      Contains trace info
    sedges : tuple
      Arrays of the slit
    xoff : float, optional
      In case Ginga has an index offset.  It appears not to
    yoff : float, optional


    Returns
    -------

    """
    # Connect
    viewer = connect_to_ginga()
    ch = viewer.channel(cname)
    canvas = viewer.canvas(ch._chname)
    # Show image, clear canvas [in case this is a repeat]
    name='image'
    ch.load_np(name, msarc, 'fits', {})
    canvas.clear()
    # Show a trace
    ntrc = len(trcdict['xtfit'])
    for idx in range(ntrc):
        if trcdict['xtfit'][idx] is None:
            continue
        x = trcdict['xtfit'][idx] + xoff
        y = trcdict['ytfit'][idx] + yoff  # FOR IMAGING (Ginga offsets this value by 1 internally)
        gdy = y > 0.
        if np.sum(gdy) > 0:
            points = list(zip(x[gdy][::pstep].tolist(),y[gdy][::pstep].tolist()))
            if trcdict['aduse'][idx]:
                clr = 'green'  # Good line
            else:
                if not all_green:
                    clr = 'red'  # Bad line
                else:
                    clr = 'green'
            canvas.add('path', points, color=clr)
    #msgs.info("Check the Ginga viewer")
    # Show slit edges
    if sedges is not None:
        y = (np.arange(msarc.shape[0]) + yoff).tolist()
        # Left
        for edge in [0,1]:
            points = list(zip(sedges[edge][::50].tolist(),y[::50]))
            canvas.add('path', points, color='cyan')
    # ALTERNATE for comparing methods
    if 'save_yt' in trcdict.keys():
        x = trcdict['xtfit'][2]
        y = trcdict['save_yt']
        gdy = y > 0.
        points = list(zip(x[gdy].tolist(),y[gdy].tolist()))
        canvas.add('path', points, color='blue')
    #debugger.set_trace()


.. highlight:: rest

*********
Keck LRIS
*********


Overview
========

This file summarizes several instrument specific
settings that are related to the Keck/LRIS spectrogrape.


Longslit
========

If reducing data with a longslit, we recommend
that you specify that only a single slit is
desired, i.e.::

    trace slits number 1

See :ref:`trace-slit-number` for further details.

.. _LRISb:

Taking Calibrations for LRISb
=============================


Default Settings
++++++++++++++++

Here are the deviations from the default settings
for LRISb::

    settings trace dispersion direction 0
    settings trace slits tilts method spca
    settings trace slits tilts params 1,1,1
    settings trace slits pca params [3,2,1,0]
    settings trace slits sigdetect 30.0        # Good for Twilight flats; faint dome flats might fail miserably..

The last setting is fine for a relatively bright frame
taken on the twilight sky,
but we suspect a faint dome flat on the blue side will require
a lower sigdetect (and is likely to be very challenging overall).

Internal flats, meanwhile, may be too bright
and need to be tested.


Pixel Flat
++++++++++

It is recommend to correct for pixel-to-pixel variations using a slitless
flat.  If you did not take such calibration frames or cannot process them,
you may wish to use an archival.  If so, copy the file into your MasterFrame
folder (should be named MF_lris_blue and you may need to create it yourself)
and set the following in the `reduce-block` (out of date!) of the PypeIt file::


    reduce flatfield useframe MF_lris_blue/PypeIt_LRISb_pixflat_B600_2x2_17sep2009.fits.gz


Taking Calibrations for LRISr
=============================

Default Settings
++++++++++++++++

Here are the deviations from the default settings
for LRISr::

    settings trace slits sigdetect 50.0   # Good for relatively bright dome flats
    settings trace slits pca params [3,2,1,0]

Known issues
++++++++++++

Multi-slit
----------

The code may identify a 'ghost' slit in empty detector real
estate if your mask does not fill most of the field.  Be prepared
to ignore it.

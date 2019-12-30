.. highlight:: rest

.. _pypeit_file:

=====================
PypeIt Reduction File
=====================

Overview
========

The primary file which informs the PypeIt data
reduction pipeline is referred to as the PypeIt
reduction file and it has a .pypeit extension.  This
can be generated from PypeIt scripts (*recommended*)
or by hand if you are sufficiently familiar with the code.

This document provides guidance on generating and modifying
the file.

We *recommend* that you generate a unique PypeIt file for each
instrument setup (modulo detectors) or for each mask.
It is possible that you will need to modify the settings for
different gratings, etc.  It will also enable you to more
easily customize the associated calibration files to process.

Types
=====

For reference, we distinguish between several types of PypeIt
files.

Instrument PypeIt file
----------------------

For each instrument being reduced in a working folder,
the top-level PypeIt file is referred to as an *instrument*
PypeIt file.  It is intended to be used to generate the
instrument :doc:`setup` file and custom PypeIt files for the
full reductions.

The standard naming for the instrument PypeIt file is::

    instrument_date.pypeit
    e.g., lris_blue_2016-Nov-23.pypeit

Custom PypeIt file
------------------

When one performs the full reduction on a set of files for
a given setup,
the *custom* PypeIt file is used.  We refer to it as custom
because it may be significantly customized for the specifc
instrument configuration and/or target.

While it is possible for a custom PypeIt files to be used
on more than one setup grouping, it is not recommended.

A typical naming scheme is by setups, e.g.::

    lris_blue_setup_A.pypeit

although specifying by instrument configuration::

    kast_blue_600_4310_d55.pypeit

or target::

    kast_blue_3C273.pypeit

may be preferable.

.. _pypeit_setup_pypeit_files:

pypeit_setup
============

By default, the pypeit_setup script will generate a set of
custom .pypeit files, one per instrument configuration.  These
will have names like::

    lris_blue_setup_A.pypeit

This is the default because we expect that most users wish to
reduce at one time the full set of exposures taken
with the same instrument configuration.
Of course, one can create other custom .pypeit files.


By Example
==========

For reference, there are
existing PypeIt files in `PypeIt development suite
<https://github.com/pypeit/PypeIt-development-suite>`_.
The PypeIt development suite is recommended for download
(see :doc:`installing`), and the relevant PypeIt files are located
in::

    PypeIt-development-suite/pypeit_files/

You should be able to find one that matches your instrument.

.. _pypfile_by_line:

Line by line
============

This section describes the various sections of a .pypeit file.
In principle, you can use the following description to build a .pypeit
file from scratch.  This is **not** recommended.
The following documentation is mainly for guiding
modifications to an existing PypeIt file.

Naming
------

Create a .pypeit file. Name it anything you want, but for example,
it's useful to have: the instrument name, the grating or grism used,
the dichroic, etc. For example, we could call our PypeIt file
'lris_blue_long_600_4000_d560.pypeit', for our data was collected
on LRIS's blue arm, in long slit mode, using the 600/4000 grism
and d560 dichroic.

You can make any comments in your PypeIt file with a
pound sign::

    # This is a comment line

We *recommend* you separate the main blocks of the .pypeit file
with comments.

.. _run_block:

Run block
---------

The first thing to include are changes to the
default settings related to running PypeIt.
The only one required is to set the name of the
spectrograph::

    run spectograph name_of_your_spectrograph

We do recommend including several others, and the
.pypeit files made by the `pypeit_pypfiles` (out of date!) script
includes most of the following.
Here are ones that one typically sets::

    # Change the default settings
    run ncpus 1                     # number of CPUs to use; can also negative integers,
                                    so -1 means all but one CPU
    run spectrograph lris_blue      # the spectrograph (+arm, if necessary) this set of data is from;
                                    see README for list of available instruments
    output verbosity 2                   # level of screen output; 0 = no output, 1 = low level of output;
                                    2 = output everything
    output overwrite True              # overwrite any existing output files?
    output sorted lris_blue_long_600_4000_d560     # name of output files

.. _reduce-block:

Reduce block
------------

bias
~~~~

If you have no bias frames and/or wish to subtract the bias with
the overscan region, then set the following::

    bias useframe overscan


Setup block
-----------

If a Setup is defined here, the value (e.g. "A" or "D") will be
used instead of starting from the default "A" value.  But *only*
if there is a single Setup in the PypeIt file.

.. _data_block:

Data block
----------

By Files
========

This is the recommended approach when performing the
full run (as opposed to :ref:`pypeit-setup`).

By Path Only
============

Next, tell PypeIt where your raw data lives!
One specifies the full path and may use wild cards
to include a set of files.  If the data are compressed,
include that extension.  Multiple entries are allowed

Here is an example::

    # Read in the data
    data read
     /Users/path/to/your/raw/data/*.fits
    data end

If you wish to skip individual files, you can specify these
without the complete path, e.g.::

    skip LB.20160406.17832.fits

These will be ignored as if they didn't exist.

.. _spect_block:

Spect block
-----------

Then, give PypeIt some information about your raw data. For
example, PypeIt only accepts calibration files if they were
created within a time window of the science frame of interest.
You can set your own time window here. PypeIt also requires a
certain number of each type of calibration file to be matched
with the science frame, and here you can set what you want the
minimum to be::

    spect read
     #fits calwin 1000.     # calibration window; default window is 12 hrs;
                            here it is changed to 1000. hrs
     pixelflat number 1       # number of pixel flats needed for data reduction
     bias number 3          # number of bias frames; note that in this case,
                            PypeIt will combine the 3 biases into a master bias
     arc number 1           # number of arcs
     trace number 1         # number of trace frames
    spect end


In addition to the basic calibration settings above, you
may wish to redefine the frametype of a given file.
Here are some examples::

    spect read
     set bias     b150910_2036.fits.gz
     set bias     b150910_2037.fits.gz
     set bias     b150910_2038.fits.gz
     set pixelflat  b150910_2051.fits.gz
     set trace    b150910_2051.fits.gz
     set standard b150910_2083.fits.gz
    spect end


Whole enchilada
---------------
With that, the most basic PypeIt file looks something like this::

    # Change the default settings
    run ncpus 1
    run spectrograph lris_blue
    output verbosity 2
    output overwrite True
    output sorted lris_blue_long_600_4000_d560

    # Read in the data
    data read
     /Users/path/to/your/raw/data/*.fits
    data end

    spect read
     #fits calwin 1000.

     pixelflat number 1
     bias number 3
     arc number 1
     trace number 1
    spect end

You can now run PypeIt with this .pypeit settings file! See how in
:doc:`running`.



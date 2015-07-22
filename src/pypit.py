import os
import sys
import getopt
from signal import SIGINT, signal as sigsignal
from warnings import resetwarnings, simplefilter
from time import time
import traceback
import numpy as np
# Import PYPIT routines
import armsgs as msgs
import arload
import arsave
import arcomb
import arproc
import arsort
import artrace
import arutils

last_updated = "Last updated 20th July 2015"

def usage(prognm):
    print "\n#####################################################################"
    print msgs.pypitheader()
    print "##  -----------------------------------------------------------------"
    print "##  Options: (default values in brackets)"
    print "##   -c or --cpus      : (all) Number of cpu cores to use"
    print "##   -h or --help      : Print this message"
    print "##   -v or --verbose   : (2) Level of verbosity (0-2)"
    print "##  -----------------------------------------------------------------"
    print "##  %s" % last_updated
    print "#####################################################################\n"
    sys.exit()

class ClassMain:

    def __init__(self, argflag, quick=False):
        """
        argflag :: A list of arguments and flags for the reduction
        quick   :: Results in a quick reduction, if a quicker (but less accurate) reduction exists
        ---------------------------------------------------

        """

        #############################
        # Set some universal parameters
        self._argflag = argflag   # Arguments and Flags
        self._transpose = False   # Determine if the frames need to be transposed
        # Arrays to store the name for the frames that have already been combined
        self._done_bias, self._name_bias = [], []
        self._done_flat, self._name_flat = [], []
        self._done_arcs, self._name_arcs = [], []
        #############################

        # First send all signals to messages to be dealt
        # with (i.e. someone hits ctrl+c)
        sigsignal(SIGINT, msgs.signal_handler)

        # Ignore all warnings given by python
        resetwarnings()
        simplefilter("ignore")

        # Record the starting time
        self._tstart=time()

        # Load the Input file
        self._parlines, self._datlines, self._spclines = arload.load_input(self)

        # Determine the type of data that is being reduced
        msgs.work("TO BE DONE")

        # If a quick reduction has been requested, make sure the requested pipeline
        # is the quick implementation (if it exists), otherwise run the standard pipeline.
        if quick:
            # Change to a "quick" settings file
            msgs.work("TO BE DONE")

        # Load the Spectrograph settings
        self._spect = arload.load_spect(self)

        # Load any changes to the spectrograph settings
        self._spect = arload.load_spect(self, lines=self._spclines)

        # Load the important information from the fits headers
        self._fitsdict = arload.load_headers(self)

        # Load the list of standard stars
        self._standardStars = arload.load_standards(self)

        # Reduce the data!
        msgs.work("Send the data away to a definition of the type of reduction needed")
        status = 0
        if quick:
            msgs.work("define what is needed here")
        else:
            success = self.ARMLSD()
        if status==0:
            msgs.info("Reduction complete")

    ###################################
    # Reduction procedures
    ###################################

    def BadPixelMask(self, sc):
        if self._argflag['reduce']['badpix']:
            msgs.info("Preparing a bad pixel mask")
            # Get all of the bias frames for this science frame
            ind = self._spect['bias']['index'][sc]
            # Load the Bias frames
            frames = arload.load_frames(self, ind, frametype='bias', trim=False)
            tbpix = arcomb.comb_frames(frames, spect=self._spect, frametype='bias', **self._argflag['bias']['comb'])
            self._bpix = arproc.badpix(self,tbpix)
            del tbpix
        else:
            msgs.info("Not preparing a bad pixel mask")
        return

    def GetDispersionDirection(self, ind):
        if self._argflag['trace']['disp']['direction'] is None:
            self._dispaxis = artrace.dispdir(self._msarc, dispwin=self._argflag['trace']['disp']['window'], mode=0)
        elif self._argflag['trace']['disp']['direction'] in [0,1]:
            self._dispaxis = int(self._argflag['trace']['disp']['direction'])
        else:
            msgs.error("The argument for the dispersion direction (trace+disp+direction)"+msgs.newline()+
                        "must be either:"+msgs.newline()+"  0 if the dispersion axis is predominantly along a row"+msgs.newline()+
                        "  1 if the dispersion axis is predominantly along a column")
        # Perform a check to warn the user if the longest axis is not equal to the dispersion direction
        if (self._msarc.shape[0] > self._msarc.shape[1]):
            if (self._dispaxis==1): msgs.warn("The dispersion axis is set to the shorter axis, is this correct?")
        else:
            if (self._dispaxis==0): msgs.warn("The dispersion axis is set to the shorter axis, is this correct?")

        ###############
        # Change dispersion direction and files if necessary
        # The code is programmed assuming dispaxis=0
        if (self._dispaxis==1):
            msgs.info("Transposing frames")
            # Flip the transpose switch
            self._transpose = True
            # Transpose the master bias frame
            if self._msbias_name is not None:
                self._msbias = self._msbias.T
                arsave.save_master(self, self._msbias, filename=self._msbias_name, frametype=self._argflag['reduce']['usebias'], ind=ind)
            # Transpose the master arc, and save it
            self._msarc = self._msarc.T
            arsave.save_master(self, self._msarc, filename=self._msarc_name, frametype='arc', ind=ind)
            # Transpose the bad pixel mask
            self._bpix = self._bpix.T
            # Change the user-specified (x,y) pixel sizes
            tmp = self._spect['det']['xgap']
            self._spect['det']['xgap'] = self._spect['det']['ygap']
            self._spect['det']['ygap'] = tmp
            self._spect['det']['ysize'] = 1.0/self._spect['det']['ysize']
            # Update the amplifier/data/overscan sections
            for i in range(self._spect['det']['numamplifiers']):
                # Flip the order of the sections
                self._spect['det']['ampsec{0:02d}'.format(i+1)] = self._spect['det']['ampsec{0:02d}'.format(i+1)][::-1]
                self._spect['det']['datasec{0:02d}'.format(i+1)] = self._spect['det']['datasec{0:02d}'.format(i+1)][::-1]
                self._spect['det']['oscansec{0:02d}'.format(i+1)] = self._spect['det']['oscansec{0:02d}'.format(i+1)][::-1]
            # Change the user-specified (x,y) pixel sizes
            msgs.work("Transpose gain and readnoise frames")
            # Set the new dispersion axis
            self._dispaxis=0
        return

    def GetPixelLocations(self):
        if self._argflag['reduce']['locations'] is None:
            self._pixlocn = artrace.gen_pixloc(self, self._mstrace,gen=True)
        elif self._argflag['reduce']['locations'] in ["mstrace"]:
            self._pixlocn = arutils.gen_pixloc(self._spect,self._mstrace,gen=False)
        else:
            self._pixlocn = arload.load_master(self._argflag['run']['masterdir']+'/'+self._argflag['reduce']['locations'], frametype=None)
        return

    def MasterArc(self, sc):
        if self._argflag['reduce']['usearc'] in ['arc']:
            msgs.info("Preparing a master arc frame")
            ind = self._spect['arc']['index'][sc]
            # Check if an *identical* master arc frame has already been produced
            self.SetFoundArc(False)
            for gb in range(len(self._done_arcs)):
                if np.array_equal(ind, self._done_arcs[gb]):
                    msgs.info("An identical master arc frame already exists")
                    msarc = arload.load_master(self._name_arcs[gb], frametype='arc')
                    self._tltprefix = os.path.splitext(os.path.basename(msarc_name))[0]
                    msarc_name = self._name_arcs[gb]
                    self.SetFoundArc(True)
            if not self._foundarc:
                # Load the arc frames
                frames = arload.load_frames(self, ind, frametype='arc', msbias=self._msbias)
                if self._argflag['reduce']['arcmatch'] > 0.0:
                    sframes = arsort.match_frames(self, frames, self._argflag['reduce']['arcmatch'], frametype='arc', satlevel=self._spect['det']['saturation']*self._spect['det']['nonlinear'])
                    subframes = np.zeros((frames.shape[0], frames.shape[1], len(sframes)))
                    numarr = np.array([])
                    for i in range(len(sframes)):
                        numarr = np.append(numarr,sframes[i].shape[2])
                        msarc = arcomb.comb_frames(sframes[i], spect=self._spect, frametype='arc', **self._argflag['arc']['comb'])
                        msarc_name = "{0:s}/{1:s}/sub-msarc{2:s}_{3:03d}_{4:03d}.fits".format(os.getcwd(),self._argflag['run']['masterdir'],self._spect["det"]["suffix"],len(self._done_arcs),i)
                        # Send the data away to be saved
                        arsave.save_master(self, msarc, filename=msarc_name, frametype='arc', ind=ind)
                        subframes[:,:,i]=msarc.copy()
                    del sframes
                    # Combine all sub-frames
                    msarc = arcomb.comb_frames(subframes, spect=self._spect, frametype='arc', weights=numarr, **self._argflag['arc']['comb'])
                    del subframes
                else:
                    msarc = arcomb.comb_frames(frames, spect=self._spect, frametype='arc', **self._argflag['arc']['comb'])
                del frames
                # Derive a suitable name for the master arc frame
                msarc_name = "{0:s}/{1:s}/msarc{2:s}_{3:03d}.fits".format(os.getcwd(),self._argflag['run']['masterdir'],self._spect["det"]["suffix"],len(self._done_arcs))
                self._tltprefix = os.path.splitext(os.path.basename(msarc_name))[0]
                # Send the data away to be saved
                arsave.save_master(self, msarc, filename=msarc_name, frametype='arc', ind=ind)
                # Store the files used and the master bias name in case it can be used during the later reduction processes
                self._done_arcs.append(ind)
                self._name_arcs.append(msarc_name)
        else:
            msarc_name = self._argflag['run']['masterdir']+'/'+self._argflag['reduce']['usearc']
            self._tltprefix = os.path.splitext(os.path.basename(msarc_name))[0]
            msarc=arload.load_master(msarc_name, frametype=None)
        return msarc, msarc_name

    def MasterBias(self, sc):
        if self._argflag['reduce']['usebias'] in ['bias', 'dark']:
            msgs.info("Preparing a master {0:s} frame".format(self._argflag['reduce']['usebias']))
            if self._argflag['reduce']['usebias'] == 'bias':
                # Get all of the bias frames for this science frame
                ind = self._spect['bias']['index'][sc]
            elif self._argflag['reduce']['usebias'] == 'dark':
                # Get all of the dark frames for this science frame
                ind = self._spect['dark']['index'][sc]
            # Check if an *identical* master bias frame has already been produced
            found = False
            for gb in range(len(self._done_bias)):
                if np.array_equal(ind, self._done_bias[gb]):
                    msgs.info("An identical master {0:s} frame already exists.".format(self._argflag['reduce']['usebias']))
                    msbias = arload.load_master(self._name_bias[gb], frametype=self._argflag['reduce']['usebias'])
                    found = True
                    break
            if not found:
                # Load the Bias/Dark frames
                frames = arload.load_frames(self, ind, frametype=self._argflag['reduce']['usebias'], transpose=self._transpose)
                msbias = arcomb.comb_frames(frames, spect=self._spect, frametype=self._argflag['reduce']['usebias'], **self._argflag['bias']['comb'])
                # Derive a suitable name for the master bias frame
                msbias_name = "{0:s}/{1:s}/msbias{2:s}_{3:03d}.fits".format(os.getcwd(),self._argflag['run']['masterdir'],self._spect["det"]["suffix"],len(self._done_bias))
                # Send the data away to be saved
                arsave.save_master(self, msbias, filename=msbias_name, frametype=self._argflag['reduce']['usebias'], ind=ind)
                # Store the files used and the master bias name in case it can be used during the later reduction processes
                self._done_bias.append(ind)
                self._name_bias.append(msbias_name)
        elif self._argflag['reduce']['usebias'] == 'overscan':
            msbias = 'overscan'
            msbias_name = None
        elif self._argflag['reduce']['usebias'] == 'none':
            msgs.info("Not performing a bias/dark subtraction")
            msbias=None
            msbias_name = None
        else: # It must be the name of a file the user wishes to load
            msbias_name = os.getcwd()+self._argflag['run']['masterdir']+'/'+self._argflag['reduce']['usebias']
            msbias = arload.load_master(self._argflag['run']['masterdir']+'/'+self._argflag['reduce']['usebias'], frametype=None)
        return msbias, msbias_name

    def MasterFlatField(self, sc):
        if self._argflag['reduce']['flatfield']: # Only do it if the user wants to flat field
            ###############
            # Generate a master pixel flat frame
            if self._argflag['reduce']['useflat'] in ['pixflat', 'blzflat']:
                msgs.info("Preparing a master pixel flat frame with {0:s}".format(self._argflag['reduce']['useflat']))
                if self._argflag['reduce']['useflat'] == 'pixflat':
                    # Get all of the pixel flat frames for this science frame
                    ind = self._spect['pixflat']['index'][sc]
                elif self._argflag['reduce']['useflat'] == 'blzflat':
                    # Get all of the blzflat frames for this science frame
                    ind = self._spect['blzflat']['index'][sc]
                # Check if an *identical* master pixel flat frame has already been produced
                found = False
                for gb in range(len(self._done_flat)):
                    if np.array_equal(ind, self._done_flat[gb]):
                        msgs.info("An identical master flat frame already exists")
                        mspixflat_name = self._name_flat[gb]
                        mspixflat = arload.load_master(mspixflat_name, frametype='pixel flat')
                        found = True
                        break
                if not found:
                    # Load the frames for tracing
                    frames = arload.load_frames(self, ind, frametype='pixel flat', msbias=self._msbias, transpose=self._transpose)
                    if self._argflag['reduce']['flatmatch'] > 0.0:
                        sframes = arsort.match_frames(self, frames, self._argflag['reduce']['flatmatch'], frametype='pixel flat', satlevel=self._spect['det']['saturation']*self._spect['det']['nonlinear'])
                        subframes = np.zeros((frames.shape[0], frames.shape[1], len(sframes)))
                        numarr = np.array([])
                        for i in range(len(sframes)):
                            numarr = np.append(numarr,sframes[i].shape[2])
                            mspixflat = arcomb.comb_frames(sframes[i], spect=self._spect, frametype='pixel flat', **self._argflag['pixflat']['comb'])
                            mspixflat_name = "{0:s}/{1:s}/sub-msflat{2:s}_{3:03d}_{4:03d}.fits".format(os.getcwd(),self._argflag['run']['masterdir'],self._spect["det"]["suffix"],len(self._done_flat),i)
                            # Send the data away to be saved
                            arsave.save_master(self, mspixflat, filename=mspixflat_name, frametype='pixel flat', ind=ind)
                            subframes[:,:,i]=mspixflat.copy()
                        del sframes
                        # Combine all sub-frames
                        mspixflat = arcomb.comb_frames(subframes, spect=self._spect, frametype='pixel flat', weights=numarr, **self._argflag['pixflat']['comb'])
                        del subframes
                    else:
                        mspixflat = arcomb.comb_frames(frames, spect=self._spect, frametype='pixel flat', **self._argflag['pixflat']['comb'])
                    del frames
                    # Derive a suitable name for the master pixel flat frame
                    mspixflat_name = "{0:s}/{1:s}/msflat{2:s}_{3:03d}.fits".format(os.getcwd(),self._argflag['run']['masterdir'],self._spect["det"]["suffix"],len(self._done_flat))
                    # Send the data away to be saved
                    arsave.save_master(self, mspixflat, filename=mspixflat_name, frametype='pixel flat', ind=ind)
                    # Store the files used and the master bias name in case it can be used during the later reduction processes
                    self._done_flat.append(ind)
                    self._name_flat.append(mspixflat_name)
            else: # It must be the name of a file the user wishes to load
                mspixflat_name = self._argflag['run']['masterdir']+'/'+self._argflag['reduce']['usepixflat']
                mspixflat=arload.load_master(mspixflat_name, frametype=None)
            # Now that the combined, master flat field frame is loaded...
            # Normalize the flat field
            mspixflatnrm, msblaze = arproc.flatnorm(self, mspixflat, overpix=0, fname=os.path.basename(os.path.splitext(mspixflat_name)[0]))
        else:
            msgs.work("Pixel Flat arrays need to be generated when not flat fielding")
            msgs.bug("Blaze is currently undefined")
            mspixflat = np.ones_like(self._msarc)
            mspixflatnrm = np.ones_like(self._msarc)
            msblaze = None
            mspixflat_name = None
        return mspixflat, mspixflatnrm, msblaze, mspixflat_name

    def MasterTrace(self, sc):
        if self._argflag['reduce']['usetrace'] in ['trace', 'blzflat']:
            msgs.info("Preparing a master trace frame with {0:s}".format(self._argflag['reduce']['usetrace']))
            if self._argflag['reduce']['usetrace'] == 'trace':
                # Get all of the trace frames for this science frame
                ind = self._spect['trace']['index'][sc]
            elif self._argflag['reduce']['usetrace'] == 'blzflat':
                # Get all of the blzflat frames for this science frame
                ind = self._spect['blzflat']['index'][sc]
            # Check if an *identical* master trace frame has already been produced
            foundtrc = False
            for gb in range(len(self._done_flat)):
                if np.array_equal(ind, self._done_flat[gb]):
                    msgs.info("An identical master trace frame already exists")
                    mstrace_name = self._name_flat[gb]
                    mstrace = arload.load_master(self._name_flat[gb], frametype='trace')
                    self._trcprefix = os.path.splitext(os.path.basename(self._name_flat[gb]))[0]
                    foundtrc = True
                    break
            if not foundtrc:
                # Load the frames for tracing
                frames = arload.load_frames(self, ind, frametype='trace', msbias=self._msbias, trim=self._argflag['reduce']['trim'], transpose=self._transpose)
                if self._argflag['reduce']['flatmatch'] > 0.0:
                    sframes = arsort.match_frames(self, frames, self._argflag['reduce']['flatmatch'], frametype='trace', satlevel=self._spect['det']['saturation']*self._spect['det']['nonlinear'])
                    subframes = np.zeros((frames.shape[0], frames.shape[1], len(sframes)))
                    numarr = np.array([])
                    for i in range(len(sframes)):
                        numarr = np.append(numarr, sframes[i].shape[2])
                        mstrace = arcomb.comb_frames(sframes[i], spect=self._spect, frametype='trace', **self._argflag['trace']['comb'])
                        mstrace_name = "{0:s}/{1:s}/sub-msflat{2:s}_{3:03d}_{4:03d}.fits".format(os.getcwd(),self._argflag['run']['masterdir'],self._spect["det"]["suffix"],len(self._done_flat),i)
                        null = os.path.splitext(os.path.basename(mstrace_name))[0]
                        self._trcprefix = "{0:s}{1:s}".format(null,self._spect["det"]["suffix"])
                        # Send the data away to be saved
                        arsave.save_master(self, mstrace, filename=mstrace_name, frametype='trace', ind=ind)
                        subframes[:,:,i]=mstrace.copy()
                    del sframes
                    # Combine all sub-frames
                    mstrace = arcomb.comb_frames(subframes, spect=self._spect, frametype='trace', weights=numarr, **self._argflag['trace']['comb'])
                    del subframes
                else:
                    mstrace = arcomb.comb_frames(frames, spect=self._spect, frametype='trace', **self._argflag['trace']['comb'])
                del frames
                # Derive a suitable name for the master trace frame
                mstrace_name = "{0:s}/{1:s}/msflat{2:s}_{3:03d}.fits".format(os.getcwd(),self._argflag['run']['masterdir'],self._spect["det"]["suffix"],len(self._done_flat))
                self._trcprefix = os.path.splitext(os.path.basename(mstrace_name))[0]
                # Send the data away to be saved
                arsave.save_master(self, mstrace, filename=mstrace_name, frametype='trace', ind=ind)
                # Store the files used and the master bias name in case it can be used during the later reduction processes
                self._done_flat.append(ind)
                self._name_flat.append(mstrace_name)
        elif self._argflag['reduce']['usetrace'] == 'science':
            msgs.work("Tracing with a science frame is not implemented")
        else: # It must be the name of a file the user wishes to load
            mstrace_name = self._argflag['run']['masterdir']+'/'+self._argflag['reduce']['usetrace']
            mstrace=arload.load_master(mstrace_name, frametype=None)
            self._trcprefix = os.path.splitext(os.path.basename(mstrace_name))[0]
        return mstrace, mstrace_name

    def Setup(self):
        # Sort the data
        msgs.bug("Files and folders should not be deleted -- there should be an option to overwrite files automatically if they already exist, or choose to rename them if necessary")
        self._filesort = arsort.sort_data(self)
        # Write out the details of the sorted files
        if self._argflag['out']['sorted'] is not None: arsort.sort_write(self)
        # Match Science frames to calibration frames
        arsort.match_science(self)
        # If the user is only debugging, then exit now
        if self._argflag['run']['calcheck']:
            msgs.info("Calibration check complete. Change the 'calcheck' flag to continue with data reduction")
            sys.exit()
        # Make directory structure for different objects
        self._sci_targs = arsort.make_dirs(self)
        return

    # Bunch of Setters
    def SetFoundArc(self, bool): self._foundarc = bool

    ###################################
    # Reduction pipelines
    ###################################

    def ARMED(self):
        """
        Automatic Reduction & Modelling of Echelle Data
        """
        success = False
        # Insert series of reduction steps here
        return success

    def ARMLSD(self):
        """
        Automatic Reduction & Modelling of Long Slit Data
        """
        success = False
        # Sort data and match calibrations to science frames
        self.Setup()
        sci = self._filesort['science']
        numsci = np.size(sci)
        if numsci == 0:
            msgs.bug("What to do if no science frames are input? The calibrations should still be processed.")
            msgs.work("Maybe assume that each non-identical arc is a science frame, but force no science extraction")
            msgs.error("No science frames are input!!")
            numsci=1
            self._argflag['run']['preponly'] = True # Prepare the calibration files, but don't do the science frame reduction
        # Reduce each science frame entirely before moving to the next one
        for sc in range(numsci):
            ###############
            # First set the index for the science frame
            scidx = self._spect['science']['index'][sc]
            sciext_name_p, sciext_name_e = os.path.splitext(self._fitsdict['filename'][scidx[0]])
            prefix = "{0:s}{1:s}".format(sciext_name_p,self._spect["det"]["suffix"])
            ###############
            # Generate master bias frame
            self._msbias, self._msbias_name = self.MasterBias(sc)
            ###############
            # Generate a bad pixel mask
            self.BadPixelMask(sc)
            ###############
            # Estimate gain and readout noise for the amplifiers
            msgs.work("Estimate Gain and Readout noise...")
            ###############
            # Generate a master arc frame
            self._msarc, self._msarc_name = self.MasterArc(sc)
            ###############
            # Determine the dispersion direction (and transpose if necessary) only on the first pass through
            self.GetDispersionDirection(self._spect['arc']['index'][sc])
            ###############
            # Generate a master trace frame
            self._mstrace, self._mstracename = self.MasterTrace(sc)
            ###############
            # Generate an array that provides the physical pixel locations on the detector
            self.GetPixelLocations()
            ###############
            # Determine the edges of the spectrum
            # Take the edges of the (trimmed) data to be the order edges
            self._lordloc = self._pixlocn[:,0,1].reshape((self._pixlocn.shape[0],1))
            self._rordloc = self._pixlocn[:,-1,1].reshape((self._pixlocn.shape[0],1))
            # Convert physical trace into a pixel trace
            msgs.info("Converting physical trace locations to nearest pixel")
            self._pixcen  = artrace.phys_to_pix(0.5*(self._lordloc+self._rordloc), self._pixlocn, self._dispaxis, 1-self._dispaxis)
            self._pixwid  = (self._rordloc-self._lordloc).mean(0).astype(np.int)
            self._lordpix = artrace.phys_to_pix(self._lordloc, self._pixlocn, self._dispaxis, 1-self._dispaxis)
            self._rordpix = artrace.phys_to_pix(self._rordloc, self._pixlocn, self._dispaxis, 1-self._dispaxis)
            ###############
            # Prepare the pixel flat field frame
            self._mspixflat, self._mspixflatnrm, self._msblaze, self._mspixflat_name = self.MasterFlatField(sc)
            ###############
            # Derive the spectral tilt
            if self._foundarc:
                try:
                    # Load the order locations from file
                    self._tilts, self._satmask = arload.load_tilts(self._msarc_name)
                except:
                    # If this fails, rederive the order tilts
                    self._tilts, self._satmask = artrace.model_tilt(self, self._msarc, prefix=prefix, trcprefix=self._trcprefix, tltprefix=self._tltprefix) # NOTE: self._tilts = tan(tilt angle), where "tilt angle" is the angle between (1) the line representing constant wavelength and (2) the column of pixels that is most closely parallel with the spatial direction of the slit.
                    # Save the traces
                    arsave.save_tilts(self, self._msarc_name)
            else:
                # First time encountered this set of arc frames --> derive the order tilts
                self._tilts, self._satmask = artrace.model_tilt(self, self._msarc, prefix=prefix, trcprefix=self._trcprefix, tltprefix=self._tltprefix) # NOTE: self._tilts = tan(tilt angle), where "tilt angle" is the angle between (1) the line representing constant wavelength and (2) the column of pixels that is most closely parallel with the spatial direction of the slit.
                # Save the tilts
                arsave.save_tilts(self, self._msarc_name)
                msgs.error("OK?")

        # Insert remaining reduction steps here
        return success

if __name__ == "__main__":
    prognm = sys.argv[0]
    debug = True
    quick = False

    # Load options from command line
    try:
        opt,arg=getopt.getopt(sys.argv[1:],'hqc:v:', ['help',
                                                 'quick'
                                                ])
    except getopt.GetoptError, err:
        msgs.error(err.msg)
        usage(prognm)
    for o,a in opt:
        if   o in ('-h', '--help')      : usage(argflag)
        elif o in ('-q', '--quick')     : quick = True
#		elif o in ('-c', '--cpus')      : argflag['run']['ncpus']     = a
#		elif o in ('-v', '--verbose')   : argflag['out']['verbose']   = int(a)

    if debug:
        argflag = arload.optarg(sys.argv, last_updated)
        ClassMain(argflag, quick=quick)
    else:
        try:
            argflag = arload.optarg(sys.argv, last_updated)
            ClassMain(argflag, quick=quick)
        except Exception:
            # There is a bug in the code, print the file and line number of the error.
            et, ev, tb = sys.exc_info()
            while tb:
                co = tb.tb_frame.f_code
                filename = str(co.co_filename)
                line_no =  str(traceback.tb_lineno(tb))
                tb = tb.tb_next
            filename=filename.split('/')[-1]
            msgs.bug("There appears to be a bug on Line "+line_no+" of "+filename+" with error:"+msgs.newline()+str(ev)+msgs.newline()+"---> please contact the author")

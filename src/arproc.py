import numpy as np
import armsgs as msgs
from matplotlib import pyplot as plt
import arcyextract
import arcyutils
import arcyproc
import artrace
import arutils
import arplot

def background_subtraction(slf, sciframe, errframe, k=3, crsigma=20.0, maskval=-999999.9, nsample=1):
    """
    Perform a background subtraction on the science frame by
    fitting a b-spline to the background.

    This routine will (probably) work poorly if the order traces overlap (of course)
    """
    retframe = np.zeros_like(sciframe)
    norders = slf._lordloc.shape[1]
    # Look at the end corners of the detector to get detector size in the dispersion direction
    xstr = slf._pixlocn[0,0,slf._dispaxis]-slf._pixlocn[0,0,slf._dispaxis+2]/2.0
    xfin = slf._pixlocn[-1,-1,slf._dispaxis]+slf._pixlocn[-1,-1,slf._dispaxis+2]/2.0
    if slf._dispaxis == 0:
        xint = slf._pixlocn[:,0,0]
    else:
        xint = slf._pixlocn[0,:,0]
    # Find which pixels are within the order edges
    msgs.info("Identifying pixels within each order")
    ordpix = arcyutils.order_pixels(slf._pixlocn, slf._lordloc, slf._rordloc, slf._dispaxis)
    allordpix = ordpix.copy()
    msgs.info("Applying bad pixel mask")
    ordpix *= (1-slf._bpix.astype(np.int))
    whord = np.where(ordpix != 0)
    msgs.info("Masking cosmic ray hits")
    crr_id = arcyutils.crreject(sciframe/np.median(sciframe[whord]), slf._dispaxis)
    cruse = np.abs(crr_id/sciframe)[whord]
    medcr = np.median(cruse)
    madcr = 1.4826*np.median(np.abs(cruse-medcr))
    whcrt = np.where(cruse>medcr+crsigma*madcr)
    whcrr = (whord[0][whcrt],whord[1][whcrt])
    msgs.info("Identified {0:d} pixels affected by cosmic rays within orders in the science frame".format(whord[0].size))
    if whcrr[0].size != 0: ordpix[whcrr] = 0
#	temp = sciframe.copy()
#	temp[whcrr] = 0.0
#	arutils.ds9plot(temp.astype(np.float))
    msgs.info("Rectifying the orders to estimate the background locations")
    msgs.work("Multiprocess this step to make it faster")
    badorders = np.zeros(norders)
    ordpixnew = np.zeros_like(ordpix)
    for o in range(norders):
        # Rectify this order
        recframe = arcyextract.rectify(sciframe, ordpix, slf._pixcen[:,o], slf._lordpix[:,o], slf._rordpix[:,o], slf._pixwid[o], maskval, slf._dispaxis)
        recerror = arcyextract.rectify(errframe, ordpix, slf._pixcen[:,o], slf._lordpix[:,o], slf._rordpix[:,o], slf._pixwid[o], maskval, slf._dispaxis)
        #recmask = np.ones(recframe.shape, dtype=np.int)
        #wmsk = np.where(recframe==maskval)
        #recmask[wmsk] = 0
        #arutils.ds9plot(recframe.astype(np.float))
        #arutils.ds9plot(recerror.astype(np.float))
        #arutils.ds9plot((recframe/recerror).astype(np.float))
        # Create a mask where there is significant flux from the object
        #flux = arcyextract.maskedaverage_order(recframe, np.ones_like(recframe), maskval)
        #plt.plot(np.arange(flux.size),flux,'k-',drawstyle='steps')
        # At least three pixels in a given row need to be detected at 2 sigma
        rowd = np.where( ((recerror!=0.0) & (recframe/recerror > 2.0)).astype(np.int).sum(axis=1) >= 3 )
        w = np.ix_(rowd[0],np.arange(recframe.shape[1]))
        #arutils.ds9plot(recframe.astype(np.float))
        #objprof = arcyextract.maskedaverage_order(recframe[w], recerror[w]**2, maskval)
        recframesmth = arcyproc.smooth_gaussmask(recframe[w], maskval, 4.0)
        #arutils.ds9plot(recframe[w].astype(np.float))
        #arutils.ds9plot(recframesmth.astype(np.float))
        objprofm = arcyextract.maskedmedian_order(recframesmth, maskval)
        if (len(objprofm) < slf._pixwid[o]/3) or (len(objprofm) <= 5):
            badorders[o] = 1
            continue
        #plt.plot(np.arange(objprof.size),objprof,'r-',drawstyle='steps')
        #plt.plot(np.arange(objprofm.size),objprofm,'k-',drawstyle='steps')
        # Reject the flux from the 3 highest S/N pixels (originally required to be included)
        xarray = np.arange(objprofm.size)
        profmask = np.zeros(objprofm.size,dtype=np.int)
        w = np.argsort(objprofm)[-3:]
        profmask[w] = 1
        # Exclude the end points
        profmask[0] = 1
        profmask[-1] = 1
        profmask, coeff = arutils.robust_polyfit(xarray, objprofm, 0, maxone=True, sigma=2.0, function="polynomial", initialmask=profmask, forceimask=True)
        #bgfit = arutils.func_val(coeff,xarray,"polynomial")
        #plt.plot(xarray,bgfit,'r-')
        w = np.where(profmask==0)
        #plt.plot(xarray[w],bgfit[w],'go')
        #plt.axis([0,30,0,np.max(objprofm)])
        #plt.show()
        #plt.clf()
        bgloc = np.zeros_like(recframe)
        bgloc[:,w] = 1.0
        # Undo the rectification
        #arutils.ds9plot(bgloc)
        unrecmask = arcyextract.rectify_undo(bgloc, slf._pixcen[:,o], slf._lordpix[:,o], slf._rordpix[:,o], slf._pixwid[o], maskval, sciframe.shape[0], sciframe.shape[1], slf._dispaxis)
        #arutils.ds9plot(unrecmask)
        # Apply the mask to master mask
        ordpixnew += unrecmask.copy()

    # Update ordpix
    msgs.info("Masking object to determine background level")
    ordpix *= ordpixnew

    msgs.work("Plot/save background locations")
    msgs.work("Deal with bad orders")
    arutils.ds9plot(ordpix.astype(np.float))

    msgs.info("Fitting and reconstructing background")
    # Create the over-sampled array of points in the dispersion direction (detector)
    ncoeff, k = sciframe.shape[slf._dispaxis], 1
    xmod = np.linspace(xstr,xfin,sciframe.shape[slf._dispaxis]*nsample)
    ycen = 0.5*(slf._lordloc + slf._rordloc)
    """
    The b-spline algorithm takes *far* too long to compute for the entire order simultaneously.
    An iteration procedure is now needed to step along the order and fit each position as you step along the order.

    Speed this up by selecting only considering pixels inside the given lordpix/rordpix of each order, rather than grabbing all xpix and ypix
    """
    bgmod = np.zeros_like(sciframe)
    polyorder, repeat = 9, 1
    for o in range(norders):
        #if o < 3 or o > norders-5: continue
        xpix, ypix = np.where(ordpix==o+1)
        print "Preparing", o+1
        xbarr, ybarr = arcyutils.prepare_bsplfit(sciframe, slf._pixlocn, slf._tilts[:,o], xmod, ycen[:,o], xpix, ypix, slf._dispaxis)
        xapix, yapix = np.where(allordpix==o+1)
        xball = arcyproc.prepare_bgmodel(sciframe, slf._pixlocn, slf._tilts[:,o], xmod, ycen[:,o], xapix, yapix, slf._dispaxis)
        ebarr = np.ones_like(xbarr)
        print "Fitting", o+1, xbarr.size
        argsrt = np.argsort(xbarr,kind='mergesort')
        polypoints = 3*slf._pixwid[o]
        fitfunc = arcyutils.polyfit_scan(xbarr[argsrt], ybarr[argsrt], ebarr, maskval, polyorder, polypoints, repeat)
        fitfunc_model = np.interp(xball, xbarr[argsrt], fitfunc)
        bgmod += arcyproc.background_model(fitfunc_model, xapix, yapix, sciframe.shape[0], sciframe.shape[1])
#		np.save("bspl/xbarr_ord{0:d}".format(o+1),xbarr)
#		np.save("bspl/ybarr_ord{0:d}".format(o+1),ybarr)
#		np.save("bspl/ebarr_ord{0:d}".format(o+1),ebarr)
#		print min(np.min(xbarr),xstr), max(np.max(xbarr),xfin), ncoeff, k, slf._dispaxis
#		np.save("bspl/pixlocn_ord{0:d}".format(o+1),slf._pixlocn)
#		np.save("bspl/tilts_ord{0:d}".format(o+1),slf._tilts[:,o])
#		np.save("bspl/xmod_ord{0:d}".format(o+1),xmod)
#		np.save("bspl/ycen_ord{0:d}".format(o+1),ycen[:,o])
#		np.save("bspl/xpix_ord{0:d}".format(o+1),xpix)
#		np.save("bspl/ypix_ord{0:d}".format(o+1),ypix)
#		print "saved all!"
#		#mod_yarr = cybspline.bspline_fit(xmod, xbarr, ybarr, ebarr, min(np.min(xbarr),xstr), max(np.max(xbarr),xfin), ncoeff, k)
#		bgmod += arcyutils.bspline_fitmod(xbarr, ybarr, ebarr, min(np.min(xbarr),xstr), max(np.max(xbarr),xfin), ncoeff, k, slf._pixlocn, slf._tilts[:,o], xmod, ycen[:,o], xpix, ypix, slf._dispaxis)

    arutils.ds9plot(bgmod.astype(np.float))
    arutils.ds9plot((sciframe-bgmod).astype(np.float))

    #exsci, exerr = arcyextract.extract_weighted(frame, error, badpixmask, pixcen[:,o], piycen[:,o], slf._pixlocn, ordtilt[:,o], ordxcen[:,o], ordycen[:,o], ordwid[:,o], ordwnum[o], ordlen[o], ordnum[o], argf_interpnum, slf._dispaxis)

# 	ordcen = 0.5*(slf._lordloc + slf._rordloc)
# 	ordwid = np.ceil(np.median(np.abs(slf._lordloc - slf._rordloc),axis=0)).astype(np.int)/2
# 	cord_pix = artrace.phys_to_pix(ordcen, slf._pixlocn, slf._dispaxis, 1-slf._dispaxis)
# 	print cord_pix
# 	print ordwid
# 	test_rect = arcyextract.rectify_fast(sciframe, cord_pix, ordwid, -999999.9, slf._dispaxis)
# 	arutils.ds9plot(test_rect)
# 	for o in range(norders):
# 		pass










    assert(False)
    # Mask out ordpix pixels where there is target flux
    ordpix = None
    # Prepare and fit the sky background pixels in every order
    msgs.work("Multiprocess this step to make it faster")
    skybg = np.zeros_like(sciframe)
    for o in range(norders):
        xpix, ypix = np.where(ordpix==1+o)
        msgs.info("Preparing sky pixels in order {0:d}/{1:d} for a b-spline fit".format(o+1,norders))
        xbarr, ybarr = cybspline.prepare_bsplfit(arc, pixmap, tilts, xmod, ycen, xpix, ypix, dispaxis)
        msgs.info("Performing b-spline fir to oversampled sky background in order {0:d}/{1:d}".format(o+1,norders))
        ncoeff, k = flt.shape[dispaxis], 1
        #mod_yarr = cybspline.bspline_fit(xmod, xbarr, ybarr, ebarr, min(np.min(xbarr),xstr), max(np.max(xbarr),xfin), ncoeff, k)
        skybg += cybspline.bspline_fitmod(xbarr, ybarr, ebarr, min(np.min(xbarr),xstr), max(np.max(xbarr),xfin), ncoeff, k, pixmap, tilts, xmod, ycen, xpix, ypix, slf._dispaxis)

    # Subtract the background
    msgs.info("Subtracting the sky background from the science frame")
    return sciframe-skybg, skybg

def badpix(slf,frame,sigdev=10.0):
    """
    frame is a master bias frame
    sigdev is the number of standard deviations away from the median that a pixel needs to be in order to be classified as a bad pixel
    """
    bpix = np.zeros_like(frame)
    nbad = 0
    for i in range (slf._spect['det']['numamplifiers']):
        datasec = "datasec{0:02d}".format(i+1)
        x0, x1, y0, y1 = slf._spect['det'][datasec][0][0], slf._spect['det'][datasec][0][1], slf._spect['det'][datasec][1][0], slf._spect['det'][datasec][1][1]
        xv=np.arange(x0,x1)
        yv=np.arange(y0,y1)
        # Construct and array with the rows and columns to be extracted
        w = np.ix_(xv,yv)
        tframe = frame[w]
        temp = np.abs(np.median(tframe)-tframe)
        sigval = max(np.median(temp)*1.4826,1.4826)
        ws = np.where(temp > sigdev*sigval)
        subfr = np.zeros_like(tframe)
        subfr[ws] = 1.0
        bpix[w] = subfr
    del subfr, tframe, temp
    # Finally, trim the bad pixel frame
    bpix=trim(slf,bpix)
    msgs.info("Identified {0:d} bad pixels".format(int(np.sum(bpix))))
    return bpix

def variance_frame(slf, sciframe, idx):
    # Dark Current noise
    dnoise = slf._spect['det']['darkcurr'] * float(slf._fitsdict["exptime"][idx])/3600.0
    # The effective read noise
    rnoise = slf._spect['det']['ronoise']**2 + (0.5*slf._spect['det']['gain'])**2
    return np.abs(sciframe) + rnoise + dnoise

def error_frame_postext(slf, sciframe, idx):
    # Dark Current noise
    dnoise = slf._spect['det']['darkcurr'] * float(slf._fitsdict["exptime"][idx])/3600.0
    # The effective read noise
    rnoise = slf._spect['det']['ronoise']**2 + (0.5*slf._spect['det']['gain'])**2
    errframe = np.zeros_like(sciframe)
    w = np.where(sciframe != -999999.9)
    errframe[w] = np.sqrt(sciframe[w] + rnoise + dnoise)
    w = np.where(sciframe == -999999.9)
    errframe[w] = 999999.9
    return errframe

def flatfield(slf, sciframe, flatframe, snframe=None):
    retframe = np.zeros_like(sciframe)
    w = np.where(flatframe > 0.0)
    retframe[w] = sciframe[w]/flatframe[w]
    if w[0].size != flatframe.size:
        w = np.where(flatframe <= 0.0)
        slf._bpix[w] = 1.0
    if snframe is None:
        return retframe
    else:
        errframe = np.zeros_like(sciframe)
        wnz = np.where(snframe>0.0)
        errframe[wnz] = retframe[wnz]/snframe[wnz]
        return retframe, errframe

def flatnorm(slf, msflat, maskval=-999999.9, overpix=6, fname=""):
    """
    Normalize the flat-field frame
    overpix/2 = the number of pixels to extend beyond each side of the order trace
    """
    msgs.info("Normalizing the master flat field frame")
    polyorder, polypoints, repeat = 4, 128, 1
    polyord_blz = 2 # This probably doesn't need to be a parameter that can be set by the user
    norders = slf._lordloc.shape[1]
    # Look at the end corners of the detector to get detector size in the dispersion direction
    xstr = slf._pixlocn[0,0,slf._dispaxis]-slf._pixlocn[0,0,slf._dispaxis+2]/2.0
    xfin = slf._pixlocn[-1,-1,slf._dispaxis]+slf._pixlocn[-1,-1,slf._dispaxis+2]/2.0
    if slf._dispaxis == 0:
        xint = slf._pixlocn[:,0,0]
    else:
        xint = slf._pixlocn[0,:,0]
    # Find which pixels are within the order edges
    msgs.info("Identifying pixels within each order")
    ordpix = arcyutils.order_pixels(slf._pixlocn, slf._lordloc, slf._rordloc, slf._dispaxis)
    msgs.info("Applying bad pixel mask")
    ordpix *= (1-slf._bpix.astype(np.int))
    msgs.info("Rectifying the orders to estimate the background locations")
    badorders = np.zeros(norders)
    msnormflat = maskval*np.ones_like(msflat)
    msblaze = maskval*np.ones((msflat.shape[slf._dispaxis],norders))
    msgs.work("Must consider different amplifiers when normalizing and determining the blaze function")
    msgs.work("Multiprocess this step to make it faster")
    flat_ext1d = maskval*np.ones((msflat.shape[slf._dispaxis],norders))
    for o in range(norders):
        # Rectify this order
        recframe = arcyextract.rectify(msflat, ordpix, slf._pixcen[:,o], slf._lordpix[:,o], slf._rordpix[:,o], slf._pixwid[o]+overpix, maskval, slf._dispaxis)
        # Take the median along the spatial dimension
        flatmed = np.median(recframe,axis=1)
        # Perform a polynomial fitting scheme to determine the blaze profile
        xarray = np.arange(flatmed.size,dtype=np.float)
        weight = flatmed.copy()
        blazet = arcyutils.polyfit_scan(xarray, flatmed.copy(), weight, maskval, polyorder, polypoints, repeat)
        # Remove the masked endpoints
        outx, outy, outm, lox, hix = arcyproc.remove_maskedends(xarray, flatmed, blazet, maskval)
        # Inspect the end points and extrapolate from the best fitting end pixels
        derv = (outm[1:]-outm[:-1])/(outx[1:]-outx[:-1])
        dervx = 0.5*(outx[1:]+outx[:-1])
        derv2 = (derv[1:]-derv[:-1])/(dervx[1:]-dervx[:-1])
        medv = np.median(derv2)
        madv = 1.4826*np.median(np.abs(derv2-medv))
        blaze = arcyproc.blaze_fitends(outx, outy, outm, derv2-medv, madv, polyord_blz, polypoints)
        #plt.plot(xarray,flatmed,'k-',drawstyle='steps')
        #plt.plot(xarray, blaze, 'r-')
        #plt.show()
        #np.savetxt("check_blaze_ord{0:d}.txt".format(o),np.transpose((xarray,flatmed)))
        # Divide the flat by the fitted flat profile
        finalblaze = np.ones(recframe.shape[0])
        finalblaze[lox:hix] = blaze.copy()
        blazenrm = finalblaze.reshape((finalblaze.size,1)).repeat(recframe.shape[1],axis=1)
        recframe /= blazenrm
        # Store the blaze for this order
        msblaze[lox:hix,o] = blaze.copy()
        flat_ext1d[:,o] = flatmed.copy()
        # Sort the normalized frames along the dispersion direction
        recsort = np.sort(recframe,axis=0)
        # Find the mean value, but only consider the "innermost" 50 per cent of pixels (i.e. the pixels closest to 1.0)
        recmean = arcyproc.scale_blaze(recsort, maskval)
        #rows = np.arange(recsort.shape[0]/4,(3*recsort.shape[0])/4,dtype=np.int)
        #w = np.ix_(rows,np.arange(recframe.shape[1]))
        #recmean = np.mean(recsort[w],axis=0)
        for i in range(recmean.size):
            recframe[i,:] /= recmean[i]
        # Undo the rectification
        normflat_unrec = arcyextract.rectify_undo(recframe, slf._pixcen[:,o], slf._lordpix[:,o], slf._rordpix[:,o], slf._pixwid[o], maskval, msflat.shape[0], msflat.shape[1], slf._dispaxis)
        # Apply the normalized flatfield for this order to the master normalized frame
        msnormflat = arcyproc.combine_nrmflat(msnormflat, normflat_unrec, slf._pixcen[:,o], slf._lordpix[:,o], slf._rordpix[:,o], slf._pixwid[o]+overpix, maskval, slf._dispaxis)
    #arutils.ds9plot(msnormflat.astype(np.float))
    #arutils.ds9plot(msblaze.astype(np.float))
    # Send the blaze away to be plotted and saved
    msgs.work("Perform a 2D PCA analysis on echelle blaze fits?")
    arplot.plot_orderfits(slf, msblaze, flat_ext1d, plotsdir=slf._argflag['run']['plotsdir'], prefix=fname+"_blaze")
    return msnormflat, msblaze

def get_wscale(slf):
    """
    This routine calculates the wavelength array based on the sampling size (in km/s) of each pixel.
    It conveniently assumes a standard reference wavelength of 911.75348 A
    """
    lam0 = 911.75348
    step = 1.0 + slf._argflag['reduce']['pixelsize']/299792.458
    # Determine the number of pixels from lam0 that need to be taken to reach the minimum wavelength of the spectrum
    msgs.work("No orders should be masked -- remove this code when the auto wavelength ID routine is fixed, and properly extrapolates.")
    w = np.where(slf._waveids!=-999999.9)
    nmin = int(np.log10(np.min(slf._waveids[w])/lam0)/np.log10(step) )
    nmax = int(1.0 + np.log10(np.max(slf._waveids[w])/lam0)/np.log10(step) ) # 1.0+ is to round up
    wave = np.min(slf._waveids[w]) * (step**np.arange(1+nmax-nmin))
    msgs.info("Extracted wavelength range will be: {0:.5f} - {1:.5f}".format(wave.min(),wave.max()))
    msgs.info("Total number of spectral pixels in the extracted spectrum will be: {0:d}".format(1+nmax-nmin))
    return wave

def sn_frame(slf, sciframe, idx):
    # Dark Current noise
    dnoise = slf._spect['det']['darkcurr'] * float(slf._fitsdict["exptime"][idx])/3600.0
    # The effective read noise
    rnoise = np.sqrt(slf._spect['det']['ronoise']**2 + (0.5*slf._spect['det']['gain'])**2)
    errframe = np.abs(sciframe) + rnoise + dnoise
    # If there are negative pixels, mask them as bad pixels
    w = np.where(errframe <= 0.0)
    if w[0].size != 0:
        msgs.warn("The error frame is negative for {0:d} pixels".format(w[0].size)+msgs.newline()+"Are you sure the bias frame is correct?")
        msgs.info("Masking these {0:d} pixels".format(w[0].size))
        errframe[w]  = 0.0
        slf._bpix[w] = 1.0
    w = np.where(errframe > 0.0)
    snframe = np.zeros_like(sciframe)
    snframe[w] = sciframe[w]/np.sqrt(errframe[w])
    return snframe

def sub_overscan(slf,file):
    for i in range (slf._spect['det']['numamplifiers']):
        oscansec = "oscansec{0:02d}".format(i+1)
        x0, x1, y0, y1 = slf._spect['det'][oscansec][0][0], slf._spect['det'][oscansec][0][1], slf._spect['det'][oscansec][1][0], slf._spect['det'][oscansec][1][1]
        xos=np.arange(x0,x1)
        yos=np.arange(y0,y1)
        w = np.ix_(xos,yos)
        try:
            oscan = file[w]
        except IndexError:
            raise IndexError('PROBABLY NEED TO TRANSPOSE')
        if oscan.shape[0] > oscan.shape[1]:  # JXP :: Looks a bit risky
            osfit = np.median(oscan,axis=1)
        else:
            osfit = np.median(oscan,axis=0)
        # Fit the overscan region
        c=np.polyfit(np.arange(osfit.size),osfit,slf._argflag['reduce']['oscanfit'])
        ossub = np.polyval(c,np.arange(osfit.size))#.reshape(osfit.size,1)
        #plt.plot(np.arange(osfit.size),osfit,'k-')
        #plt.plot(np.arange(osfit.size),ossub,'r-')
        #plt.show()
        #plt.clf()
        ossub = ossub.reshape(osfit.size,1)
        #print osfit
        #print ossub
        # Determine the section of the chip that is read out by the amplifier
        ampsec = "ampsec{0:02d}".format(i+1)
        x0, x1, y0, y1 = slf._spect['det'][ampsec][0][0], slf._spect['det'][ampsec][0][1], slf._spect['det'][ampsec][1][0], slf._spect['det'][ampsec][1][1]
        xam=np.arange(x0,x1)
        yam=np.arange(y0,y1)
        w = np.ix_(xam,yam)
        if w[0].shape[0] == ossub.shape[0]:
            file[w] -= ossub
        elif w[1].shape[1] == ossub.shape[0]:
            file[w] -= ossub.T
        else:
            msgs.error("Could not subtract bias from overscan region --"+msgs.newline()+"size of extracted regions does not match")
    del xam, yam, xos, yos, oscan
    return file



def trim(slf,file):
    for i in range (slf._spect['det']['numamplifiers']):
        datasec = "datasec{0:02d}".format(i+1)
        x0, x1, y0, y1 = slf._spect['det'][datasec][0][0], slf._spect['det'][datasec][0][1], slf._spect['det'][datasec][1][0], slf._spect['det'][datasec][1][1]
        if i==0:
            xv=np.arange(x0,x1)
            yv=np.arange(y0,y1)
        else:
            xv = np.unique(np.append(xv,np.arange(x0,x1)))
            yv = np.unique(np.append(yv,np.arange(y0,y1)))
    # Construct and array with the rows and columns to be extracted
    w = np.ix_(xv,yv)
#	if len(file.shape) == 2:
#		trimfile = file[w]
#	elif len(file.shape) == 3:
#		trimfile = np.zeros((w[0].shape[0],w[1].shape[1],file.shape[2]))
#		for f in range(file.shape[2]):
#			trimfile[:,:,f] = file[:,:,f][w]
#	else:
#		msgs.error("Cannot trim {0:d}D frame".format(int(len(file.shape))))
    return file[w]
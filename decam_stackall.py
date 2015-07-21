#!/usr/bin/env python

import pyfits, os, sys
import numpy

sys.path.append("/work/podi_prep56")

import podi_swarpstack
from podi_commandline import *
import podi_logging
import podi_sitesetup as sitesetup

wcs_headers = [
    'CRVAL1', 'CRVAL2',
    'CD1_1', 'CD1_2', 'CD2_1', 'CD2_2',
    'CUNIT1', 'CUNIT2',
    'CTYPE1', 'CTYPE2',
    ]
del_headers = ['COMAG']



def run_swarp(input_list, outputimage, combine="MEDIAN"):

    #
    #
    # Now run swarp to create the stack
    # Use the default.swarp from QR as the fallback swarp option
    #
    #

    # make sure the input is a list and not just a single filename
    if (not type(input_list) == list):
        input_list = [input_list]

    # trim the .fits from the filename (if exists), we'll add that back in later
    if (outputimage.endswith(".fits")):
        outputimage = outputimage[:-5]

    dic = {
        'swarp_default': "/work/podi_devel/.config/swarp.default",
        'img_out': "%s.fits" % (outputimage),
        'img_w_out': "%s.weight.fits" % (outputimage),
        'resample_dir': "/scratch/", 
        'inputfile': " ".join(input_list),
        'combinetype': combine.upper(),
    }
    swarp_opts = """
             -c %(swarp_default)s 
             -IMAGEOUT_NAME %(img_out)s 
             -WEIGHTOUT_NAME %(img_w_out)s
             -PIXEL_SCALE 0
             -PIXELSCALE_TYPE MEDIAN
             -COMBINE Y 
             -COMBINE_TYPE %(combinetype)s
             -CENTER_TYPE ALL
             -RESAMPLE_DIR %(resample_dir)s 
             -SUBTRACT_BACK Y
             -FSCALE_KEYWORD XXXXXXXX
             -FSCALE_DEFAULT 1.0
             -WEIGHT_TYPE NONE
             -WEIGHT_SUFFIX .weight.fits 
             -RESCALE_WEIGHTS N
             -DELETE_TMPFILES Y
             %(inputfile)s 
             """ % dic
    swarp_cmd = "%s %s" % (sitesetup.swarp_exec, swarp_opts)

    logger.info("Creating %s" % (dic['img_out']))
    logger.debug(" ".join(swarp_cmd.split()))

    try:
        ret = subprocess.Popen(swarp_cmd.split(), 
                                   stdout=subprocess.PIPE, 
                                   stderr=subprocess.PIPE)
        (swarp_stdout, swarp_stderr) = ret.communicate()

        # logger.debug("swarp stdout:\n"+swarp_stdout)
        # if (len(swarp_stderr) > 0 and ret.returncode != 0):
        #     logger.warning("swarp stderr:\n"+swarp_stderr)
        # else:
        #     logger.debug("swarp stderr:\n"+swarp_stderr)
        #print "\n".join(swarp_stderr)
        logger.info("done, swarp returned (ret-code: %d)!" % ret.returncode)
    except OSError as e:
        podi_logging.log_exception()
        print >>sys.stderr, "Execution failed:", e


    return

if __name__ == "__main__":

    options = podi_logging.setup_logging()
    logger = logging.getLogger("StackAll")


    filename = get_clean_cmdline()[1]

    print filename
    hdulist = pyfits.open(filename)


    #
    # Gather a list of all <GRPID>s and corresponding <EXTNAME>s
    #
    grpids = {}
    print "Gathering exposure groupings"
    for ext in hdulist[1:]:
        grpid = ext.header['COGRPID']
        extname = ext.name

        if not grpid in grpids:
            grpids[grpid] = []

        grpids[grpid].append(extname)
        
    print "Found %d groups" % (len(grpids))


    # Now go through each GRPID, and extract the relevant extensions to a separate file
    print "Extracting and stacking individual objects"
    for idx, grpid in enumerate(grpids):
        print grpids[grpid]
        
        out_mef = "single_object__%05d.sidereal.fits" % (grpid)
        out_mef_nonsid = "single_object__%05d.nonsid.fits" % (grpid)

        out_hdulist = [hdulist[0]]
        for single_epoch_extname in grpids[grpid]:
            out_hdulist.append(hdulist[single_epoch_extname])
            
        # convert list of extensions into MEF
        out_hdulist = pyfits.HDUList(out_hdulist)

        # Now make every image self-contained by copying and inserting 
        # all WCS-relevant headers from the primary header
        for ext in out_hdulist[1:]:
            for hdrname in wcs_headers:
                if hdrname in out_hdulist[0].header:
                    ext.header[hdrname] = out_hdulist[0].header[hdrname]

        clobberfile(out_mef)
        out_hdulist.writeto(out_mef, clobber=True, output_verify='silentfix')

        run_swarp(input_list=out_mef, 
                  outputimage="single_object__%05d.stack" % (grpid),
                  combine='average'
        )


        #
        # Now modify the FITS file to yield the non-sidereal stack
        # This requires to change the WCS headers to make the asteroid appear 
        # at the same sky-position in all individual frames
        # Since CRVALi are already identifical for all images, make the CRPIXi
        # match as well
        #
        
        # Read the CRPIXi headers from extension 1 (the first actual image ext)
        crpix_keys = ['CRPIX1', 'CRPIX2']
        crpix_ref = numpy.zeros((2))
        for idx, key in enumerate(crpix_keys):
            crpix_ref[idx] = out_hdulist[1].header[key]
        
        # and insert into all other ImageHDUs
        for ext in out_hdulist[2:]:
            for idx, key in enumerate(crpix_keys):
                ext.header[key] = crpix_ref[idx]

        # Save the modified HDUList as a sidereal stack
        clobberfile(out_mef_nonsid)
        out_hdulist.writeto(out_mef_nonsid, clobber=True)

        run_swarp(input_list=out_mef_nonsid, 
                  outputimage="single_object__%05d.nonsidstack" % (grpid),
                  combine='median'
        )


        break

        if (idx > 5):
            break

    podi_logging.shutdown_logging(options)

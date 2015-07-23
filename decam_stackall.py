#!/usr/bin/env python

import pyfits, os, sys
import numpy

qr_dir = "/work/podi_devel/"
sys.path.insert(0, qr_dir)
sys.path.insert(0, "/work/podi_devel/test/")
print sys.path

import podi_swarpstack
from podi_commandline import *
from podi_definitions import *
import podi_logging
import podi_sitesetup as sitesetup
from astLib import astWCS

from profile import *
from meanprofile import compute_mean_profile
from compare_profiles import *

wcs_headers = [
    'CRVAL1', 'CRVAL2',
    'CD1_1', 'CD1_2', 'CD2_1', 'CD2_2',
    'CUNIT1', 'CUNIT2',
    'CTYPE1', 'CTYPE2',
    ]
del_headers = ['COMAG']



def run_swarp(input_list, outputimage, combine="MEDIAN"):

    logger = logging.getLogger("Swarp")

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
        'swarp_default': "%s/.config/swarp.default" % (qr_dir),
        'img_out': "%s.fits" % (outputimage),
        'img_w_out': "%s.weight.fits" % (outputimage),
        'resample_dir': sitesetup.scratch_dir,
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


def do_work(hdulist, grpids, grpid):

    logger = logging.getLogger("Worker(%05d)" % (grpid))

    if (True):
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

        sidereal_stack = "single_object__%05d.siderealstack" % (grpid)
        run_swarp(input_list=out_mef, 
                  outputimage="single_object__%05d.stack" % (grpid),
                  combine='average'
        )
        run_swarp(input_list=out_mef, 
                  outputimage=sidereal_stack,
                  combine='median'
        )


        #
        # Now modify the FITS file to yield the non-sidereal stack
        # This requires to change the WCS headers to make the asteroid appear 
        # at the same sky-position in all individual frames
        # Since CRVALi are already identifical for all images, make the CRPIXi
        # match as well
        #
        
        # Read the CRPIXi headers from extension 1 (the first actual image ext)
        #       print ast_radec
        crpix_keys = ['CRPIX1', 'CRPIX2']
        crpix_ref = numpy.zeros((2))
        for idx, key in enumerate(crpix_keys):
            crpix_ref[idx] = out_hdulist[1].header[key]
        
        # and insert into all other ImageHDUs
        for ext in out_hdulist[2:]:
            # Work out Ra/dec of the asteroid in the i-th frame
            for idx, key in enumerate(crpix_keys):
                ext.header[key] = crpix_ref[idx] 
            # CHANGE SIGN HERE ----
            ext.header['CRPIX1'] -= out_hdulist[1].header['COX'] - ext.header['COX']
            ext.header['CRPIX2'] -= out_hdulist[1].header['COY'] - ext.header['COY']
            # CHANGE SIGN HERE ----

        # Save the modified HDUList as a sidereal stack
        clobberfile(out_mef_nonsid)
        out_hdulist.writeto(out_mef_nonsid, clobber=True)

        nonsidereal_stack_filename = "single_object__%05d.nonsidstack" % (grpid)
        run_swarp(input_list=out_mef_nonsid, 
                  outputimage=nonsidereal_stack_filename,
                  combine='median'
        )


        #
        # Now run SourceExtractor on the sidereal median stack to get a catalog 
        # of stars we can use to model the PSF 
        #
        sidereal_cat_filename = "single_object__%05d.sidereal.cat" % (grpid)
        sex_cmd = """
            %(sex)s -c %(qr_base)s/.config/wcsfix.sex 
            -PARAMETERS_NAME %(qr_base)s/.config/wcsfix.sexparam
            -CATALOG_NAME %(catfile)s
            -WEIGHT_TYPE MAP_WEIGHT
            -WEIGHT_IMAGE %(siderealstack)s.weight.fits
            %(siderealstack)s.fits """ % {
                'sex': sitesetup.sextractor,
                'qr_base': qr_dir,
                'catfile': sidereal_cat_filename,
                'siderealstack': "single_object__%05d.siderealstack" % (grpid),
            }
        print " ".join(sex_cmd.split())
        os.system(" ".join(sex_cmd.split()))


        #
        # Load the catalog
        #
        catalog = numpy.loadtxt(sidereal_cat_filename)
        print "catalog shape:", catalog.shape
        if (catalog.shape[0] == 0):
            logger.critical("no stars found!")
            return
        elif (catalog.ndim == 1):
            logger.warning("Only one source found!")
            catalog = catalog.reshape((1, -1))

        # Now exclude all sources we consider not useful
        elongated = catalog[:, SXcolumn['elongation']] > 1.5
        bad_photometry = catalog[:, SXcolumn['mag_aper_2.0']] > 50.
        flagged = (catalog[:, SXcolumn['flags']].astype(numpy.int8) & 0b00001111) > 0
        faint = catalog[:, SXcolumn['flux_max']] < 250.
        # excludes 
        # - near neighbors, 
        # - deblended sources, 
        # - saturated, 
        # - truncated (next to image boundary)
        bad_source = elongated | bad_photometry | flagged | faint
        catalog = catalog[~bad_source]

        numpy.savetxt(sys.stdout, catalog[:,:4], "%.5f")

        # Now compute mean profile for all stars combined
        stars = compute_mean_profile(filename=sidereal_stack+".fits",
                                     fxy_list=catalog[:, 2:4], 
                                     # we need x/y coords here 
                                     mode='radial',
                                     save_tmp=False,
        )
        print stars.shape
        all_profiles_filename = "single_object__%05d.sidereal.prof" % (grpid)
        numpy.savetxt(all_profiles_filename, stars)

        # Add a mirrored profile to make sure we get a symmetric profile and 
        # to avoid potentialissues at radius 0
        stars_mirrored = numpy.append((stars * numpy.array([-1.,1.]))[::-1],
                                      stars,
                                      axis=0)
        numpy.savetxt("mirrored", stars_mirrored)

        #
        # Fit the profile with a spline - use radii up to 4.8 arcsec
        #
        #base_points = numpy.linspace(-4.8,4.8,200)
        base_points = numpy.linspace(-3., 3., 90)
        psf_fit = scipy.interpolate.LSQUnivariateSpline(
            stars_mirrored[:,0], stars_mirrored[:,1],
            t=base_points,
            k=3)

        # save some data for plotting/debugging
        numpy.savetxt("fit.dump", numpy.append(base_points.reshape((-1,1)),
                                               psf_fit(base_points).reshape((-1,1)),
                                               axis=1))
        
        #
        # Now we have the full reference star PSF profile, so we can extract 
        # the profile of the actual moving target and compare them
        #
        # For now assume spherical symmetry as well for the asteroid
        #

        # Run source extractor on the median-filtered nonsidereal stack to get 
        # the position of the moving target
        nonsidereal_cat_filename = "single_object__%05d.nonsidereal.cat" % (grpid)
        sex_cmd = """
            %(sex)s -c %(qr_base)s/.config/wcsfix.sex 
            -PARAMETERS_NAME %(qr_base)s/.config/wcsfix.sexparam
            -CATALOG_NAME %(catfile)s
            -WEIGHT_TYPE MAP_WEIGHT
            -WEIGHT_IMAGE %(nonsiderealstack)s.weight.fits
            %(nonsiderealstack)s.fits """ % {
                'sex': sitesetup.sextractor,
                'qr_base': qr_dir,
                'catfile': nonsidereal_cat_filename,
                'nonsiderealstack': nonsidereal_stack_filename
            }
        print " ".join(sex_cmd.split())
        os.system(" ".join(sex_cmd.split()))

        # As above, load the catalog
        nonsid_catalog = numpy.loadtxt(nonsidereal_cat_filename)
        print "catalog shape:", nonsid_catalog.shape
        if (nonsid_catalog.shape[0] == 0):
            logger.critical("no stars found in non-sidereal catalog!")
            return
        elif (nonsid_catalog.ndim == 1):
            logger.warning("Only one source found in non-sidereal catalog!")
            nonsid_catalog = nonsid_catalog.reshape((1, -1))

        use_old_selection = False
        if (use_old_selection):
            elongated = nonsid_catalog[:, SXcolumn['elongation']] > 1.5
            bad_photometry = nonsid_catalog[:, SXcolumn['mag_aper_2.0']] > 50.
            flagged = nonsid_catalog[:, SXcolumn['flags']].astype(numpy.int8) > 0
            bad = elongated | bad_photometry | flagged
            _nonsid_catalog = nonsid_catalog[~bad]
            print _nonsid_catalog

            if (_nonsid_catalog.size > 0):
                # We still have some valid sources left
                nonsid_catalog = _nonsid_catalog
            else:
                # seelction above left no sources behind, so skip catalo pruning
                pass


            # pick the brightest source - that's very likely the asteroid
            brightest = numpy.argmin(nonsid_catalog[:, SXcolumn['mag_aper_2.0']])
            asteroid_coords = nonsid_catalog[brightest]
            print asteroid_coords[:4]
        else:
            # This is the new, preferred way:
            # Select the most central object

            # get size of output frame
            nsstack_hdu = pyfits.open(nonsidereal_stack_filename+".fits")
            center = numpy.array(nsstack_hdu[0].data.shape)/2
            
            d_center = nonsid_catalog[:, SXcolumn['x']:SXcolumn['y']+1] - center
            r_center = numpy.hypot(d_center[:,0], d_center[:,1])
            print r_center

            most_central = numpy.argmin(r_center)
            asteroid_coords = nonsid_catalog[most_central]


        #
        # Now we can also extract the profile for the asteroid
        #
        ns_hdu = pyfits.open(nonsidereal_stack_filename+".fits")
        data = ns_hdu[0].data
        arcsec_per_pixel = ns_hdu[0].header['CD2_2'] * 3600
   
        r, f, nf = get_profile(data=data.T, 
                         center_x=asteroid_coords[SXcolumn['x']], 
                         center_y=asteroid_coords[SXcolumn['y']],
                         mx=0, my=0, width=5.0/arcsec_per_pixel,
                         mode='radial',
                         normalize=False,
                         )
        r *= arcsec_per_pixel
        max_intensity = nf
        asteroid = numpy.append(r.reshape((-1,1)), f.reshape((-1,1)), axis=1)

        # asteroid = compute_mean_profile(filename=nonsidereal_stack_filename+".fits",
        #                                 fxy_list=asteroid_coords[2:4].reshape((1,-1)), 
        #                                 # we need x/y coords here 
        #                                 mode='radial',
        #                                 save_tmp=False,
        #                                 )


        #
        # Next we compute the synthetic PSF at the coordinates of the source
        #
        asteroid_synth = psf_fit(asteroid[:,0]) * max_intensity
        #
        # Make a nice output file that contains the data, sf-fit, and difference
        #
        output_buffer = numpy.empty((asteroid.shape[0], 4))
        output_buffer[:,0:2] = asteroid[:]
        output_buffer[:,2] = asteroid_synth[:]
        output_buffer[:,3] = asteroid[:,1] - asteroid_synth[:]
        psf_sub_filename = "single_object__%05d.psfsub.dat" % (grpid)
        numpy.savetxt(psf_sub_filename,
                      output_buffer,
                      header="""\
Column 01: Radius in arcsec
Column 02: Normalized flux
Column 03: Average PSF profile, normalized
Column 04: Asteroid flux - PSF flux
-------------------""",
                  )

        numpy.savetxt("asteroid.dat", 
                      numpy.append(asteroid,
                                   asteroid_synth.reshape((-1,1)), axis=1))


        #
        # Make a second plot showing the data, reference psf, and difference
        #
        plot_filename = "single_object__%05d.psfdiff.png" % (grpid)
        fig2 = matplotlib.pyplot.figure()
        ax2 = fig2.add_subplot(111)
        ax2.set_title("Asteroid @ x=%.1f y=%.1f" % (asteroid_coords[SXcolumn['x']],
                                                    asteroid_coords[SXcolumn['y']]))
        plot_x = numpy.linspace(numpy.min(stars_mirrored[:,0]), numpy.max(stars_mirrored[:,0]), 1000)
        scale_match = 1. * max_intensity
        difference = asteroid[:,1] - asteroid_synth
        ax2.scatter(asteroid[:,0], asteroid[:,1], marker=".")
        ax2.scatter(asteroid[:,0], difference-0.2*nf, marker="x")
        ax2.plot(plot_x, psf_fit(plot_x)*scale_match, "r-", linewidth=3)
        ax2.set_xlim((0, numpy.max(base_points)))
        ax2.set_ylim((-0.4*max_intensity, 1.2*max_intensity))
        ax2.axhline(color='#80ff80', linewidth=3)
        ax2.axhline(y=-0.2*max_intensity, color='#8080ff', linewidth=3)
        ax2.set_xlabel("Radius [arcsec]")
        ax2.set_ylabel("Intensity [flux-normalized counts]")
        #fig2.show()
        fig2.savefig(plot_filename)
        #matplotlib.pyplot.show()



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
        try:
            do_work(hdulist, grpids, grpid)
        except:
            podi_logging.log_exception()
            logger.error("There was a problem processing %d" % (grpid))
            pass

        break
        if (idx > 5):
            break

    podi_logging.shutdown_logging(options)

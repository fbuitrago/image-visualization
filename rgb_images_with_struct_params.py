from astropy import wcs
from astropy.io import fits
import aplpy
import numpy as np
import matplotlib.pyplot as plt
import pdb
from astropy.table import Table


#---------------------------------constants------------------------------------
fig = plt.figure(figsize=(20,20))


#-----------------------------read catalogs------------------------------------
sample       = "../galaxy_selection/galaxy_selection_in_the_i_band.txt"
tt = Table.read(sample, names=("cataid","zz","logmass","re"), format='ascii.commented_header')
cataid  = tt['cataid']
zz      = np.array(tt['zz'],dtype=float)
logmass = np.array(tt["logmass"],dtype=float)
re      = tt["re"]

#catalogs g band
tt_1g = Table.read("../analysis_filter_g/gama_kids_g_single_sersic_fits.cat",\
     names = ["gal_id","zz","mag","mag_error","mag_exist","re_pix","re_pix_error","re_exist","re_kpc","re_kpc_err","n","n_error","n_exist",\
              "ar","ar_error","ar_exist","pa","pa_error","pa_exist","index","flag_obj_detected"],\
     format="ascii.commented_header")

tt_2g = Table.read("../analysis_filter_g/gama_kids_g_bulge_disk_decompositions.cat",\
     names = ["gal_id","zz","mag1","mag_error1","mag1_exist","re_pix1",\
              "re_pix_error1","re1_exist","re_kpc1","re_kpc_error1","n1","n_error1","n1_exist",\
              "ar1","ar_error1","ar1_exist","pa1","pa_error1","pa1_exist","mag2",\
              "mag_error2","mag2_exist","re_pix2","re_pix_error2","re2_exist","re_kpc2","re_kpc_error2",\
              "n2","n_error2","n2_exist","ar2","ar_error2","ar2_exist","pa2",\
              "pa_error2","pa2_exist","index","flag_obj_detected"],\
     format="ascii.commented_header")

#catalogs r band
tt_1r = Table.read("../analysis_filter_r/gama_kids_r_single_sersic_fits.cat",\
     names = ["gal_id","zz","mag","mag_error","mag_exist","re_pix","re_pix_error","re_exist","re_kpc","re_kpc_err","n","n_error","n_exist",\
              "ar","ar_error","ar_exist","pa","pa_error","pa_exist","index","flag_obj_detected"],\
     format="ascii.commented_header")

tt_2r = Table.read("../analysis_filter_r/gama_kids_r_bulge_disk_decompositions.cat",\
     names = ["gal_id","zz","mag1","mag_error1","mag1_exist","re_pix1",\
              "re_pix_error1","re1_exist","re_kpc1","re_kpc_error1","n1","n_error1","n1_exist",\
              "ar1","ar_error1","ar1_exist","pa1","pa_error1","pa1_exist","mag2",\
              "mag_error2","mag2_exist","re_pix2","re_pix_error2","re2_exist","re_kpc2","re_kpc_error2",\
              "n2","n_error2","n2_exist","ar2","ar_error2","ar2_exist","pa2",\
              "pa_error2","pa2_exist","index","flag_obj_detected"],\
     format="ascii.commented_header")

#catalogs i band
tt_1i = Table.read("../analysis_filter_i/gama_kids_i_single_sersic_fits.cat",\
     names = ["gal_id","zz","mag","mag_error","mag_exist","re_pix","re_pix_error","re_exist","re_kpc","re_kpc_err","n","n_error","n_exist",\
              "ar","ar_error","ar_exist","pa","pa_error","pa_exist","index","flag_obj_detected"],\
     format="ascii.commented_header")

tt_2i = Table.read("../analysis_filter_i/gama_kids_i_bulge_disk_decompositions.cat",\
     names = ["gal_id","zz","mag1","mag_error1","mag1_exist","re_pix1",\
              "re_pix_error1","re1_exist","re_kpc1","re_kpc_error1","n1","n_error1","n1_exist",\
              "ar1","ar_error1","ar1_exist","pa1","pa_error1","pa1_exist","mag2",\
              "mag_error2","mag2_exist","re_pix2","re_pix_error2","re2_exist","re_kpc2","re_kpc_error2",\
              "n2","n_error2","n2_exist","ar2","ar_error2","ar2_exist","pa2",\
              "pa_error2","pa2_exist","index","flag_obj_detected"],\
     format="ascii.commented_header")

for ii in range(len(tt_1g["gal_id"])):
    gal = str(tt_1g["gal_id"][ii])

    #-----reading the fits image-----
    # img = fits.open(images_iband_path+str(gal)+"_cut.fits")
    img = fits.open(str(gal)+"_cube_2d.fits")

    subfig = aplpy.FITSFigure(img, figure=fig)

    subfig.show_rgb(str(gal)+"_rgb.png")
    # subfig.show_colorscale(vmin=-0.001,vmax=10,vmid=-0.011,stretch='log',cmap="cubehelix")

    # #-----scale bar-----
    # subfig.add_scalebar(1./3600.)
    # subfig.scalebar.set_color('red')
    # subfig.scalebar.set_label('1 arcsec')

    pos_in_first_cat = cataid == tt_1g["gal_id"][ii]

    #-----galaxy data-----
    subfig.add_label(0.94,0.9,str(gal),relative=True,color='white',size=24)
    #subfig.add_label(0.9,0.85,"r_e = {0:.2f}".format(tt_1g["re_kpc"][ii]),relative=True,color="white",size=24)
    subfig.add_label(0.88,0.85,r"r$_{e,g}$ = "+str(round(float(tt_1g["re_kpc"][ii]),2))+" kpc",relative=True,color="white",size=24)
    subfig.add_label(0.88,0.80,r"r$_{e,r}$ = "+str(round(float(tt_1r["re_kpc"][ii]),2))+" kpc",relative=True,color="white",size=24)
    subfig.add_label(0.88,0.75,r"r$_{e,i}$ = "+str(round(float(tt_1i["re_kpc"][ii]),2))+" kpc",relative=True,color="white",size=24)
    subfig.add_label(0.88,0.70,r"z = "       +str(round(float(zz     [pos_in_first_cat][0]),2)),relative=True,color="white",size=24) #the [0] is because it comes from a filter
    subfig.add_label(0.88,0.65,r"log mass = "+str(round(float(logmass[pos_in_first_cat][0]),2)),relative=True,color="white",size=24) #the [0] is because it comes from a filter
    # subfig.add_label(0.9,0.85,r"$\rm z =$"+str(zz),relative=True,color='white')
    # subfig.add_label(0.86,0.8,r"$\rm \log_{10} M_{\star} =$"+str(mass),relative=True,color='white')
    # subfig.add_label(0.13,0.15,r"$\rm mag=$ "+str(mag),relative=True,color='grey')
    # subfig.add_label(0.11,0.1,r"$\rm r_e =$"+str(nn),relative=True,color='grey')
    # subfig.add_label(0.1,0.05,r"$\rm n =$"+str(nn),relative=True,color='grey')

    #-----axis-----
    subfig.axis_labels.hide()
    subfig.tick_labels.hide()
    subfig.ticks.hide()

    plt.savefig(gal+"_rgb_with_params.pdf", format='pdf', dpi=100)
    plt.clf()

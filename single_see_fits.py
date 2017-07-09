from astropy import wcs
from astropy.io import fits
import aplpy
import numpy as np
import matplotlib.pyplot as plt
import pdb
from astropy.table import Table
import subprocess
import os

def getting_filename_from_pattern_running_ls(gal,band,type_of_file):
    #I get the filename if matching the pattern I am looking for
    result = subprocess.check_output("ls "+str(gal)+"*_"+band+"_*"+type_of_file+"_sb.fits",shell=True)
    #making readable the result string
    filename = str(result[:-1],'utf-8')

    return(filename)

def plotting_options():
    #-----axis-----
    subfig.axis_labels.hide()
    subfig.tick_labels.hide()
    subfig.ticks.hide()

def finding_coo_brightest_SB(image):
    #if the brightest pixel was the central, this function will only be return np.unravel_index(np.nanargmin( image ),image.shape)
    x_size,y_size = image.shape
    #looking in the inner part of the stamp
    offset_x = x_size/50
    offset_y = y_size/50
    submatrix = image[round((x_size/2)-offset_x):round((x_size/2)+offset_x),round((y_size/2)-offset_y):round((y_size/2)+offset_y)]
    indices_min_value = np.unravel_index(np.nanargmin(submatrix),submatrix.shape)
    filter_pos = image == submatrix[indices_min_value[0],indices_min_value[1]] #I assume that there is only a single pixel with the brightest value
    pos = np.where(filter_pos)
    return(pos[0][0],pos[1][0])
    #np.nanargmin->gets the minimum element of the image without taking into account the NaN values
    #np.unravel_index->gets the two dimensional position from a coordinate in a flatten vector

def printing_effective_radius(image,pix_scale,re_arcsec,ar,pa):
    offset_x_by_eye = (-1.*pix_scale)/3600.
    offset_y_by_eye = (1.*pix_scale)/3600.

    yy,xx=finding_coo_brightest_SB(image)
    #ellipses
    ra, dec = subfig.pixel2world(xx,yy)
    aa = re_arcsec/3600.
    bb = (re_arcsec*ar)/3600.
    pa = pa-90.
    #     pa=   catalog['pa_sex'][(ii*8)+(jj*4)]-90.
    subfig.show_ellipses(ra+offset_x_by_eye, dec+offset_y_by_eye, 2.*aa, 2.*bb, pa, edgecolor='#FFA500', linewidth=1)

def get_table_file(band,type_of_analysis):
    path = "../analysis_filter_"+band+"/"
    if type_of_analysis == "single":
        file_to_read = "gama_kids_"+band+"_single_sersic_fits.cat"
    else:
        file_to_read = "gama_kids_"+band+"_bulge_disk_decompositions.cat"
    return(path+file_to_read)

def get_struct_param_best_fit(gal,band,type_of_analysis):
    stars = np.array(["psf_for_"+str(gal)+"_"+band,"star_in_185122_kids_"+band],dtype=str) #same order as they appear in autofit
    masks = np.array(["normal_mask","central_mask"]                            ,dtype=str) #same order as they appear in autofit

    table_to_read = get_table_file(band,type_of_analysis)

    if type_of_analysis == "single":
        tt = Table.read(table_to_read, \
             names = ["gal_id","zz","mag","mag_error","mag_exist","re_pix","re_pix_error","re_exist","re_kpc","re_kpc_err","n","n_error","n_exist",\
                      "ar","ar_error","ar_exist","pa","pa_error","pa_exist","index","flag_obj_detected"],\
             format="ascii.commented_header")
        filter_gal = tt["gal_id"] == gal
        index = tt["index"][filter_gal]
        index = np.float(index[0]) #in order to be able to do np.isnan afterwards
        if np.isnan(index): #if the there is no best fit, return nothing
            return()
        else:
            index = int(index)
            re_pix = tt["re_pix"][filter_gal]; re_pix = np.float(re_pix[0])
            re_kpc = tt["re_kpc"][filter_gal]; re_kpc = np.float(re_kpc[0])
            nn     = tt["n"]     [filter_gal]; nn     = np.float(nn[0])
            ar     = tt["ar"]    [filter_gal]; ar     = np.float(ar[0])
            pa     = tt["pa"]    [filter_gal]; pa     = np.float(pa[0])
            return(re_pix,re_kpc,nn,ar,pa)
    elif type_of_analysis == "bulge_disk":
        tt = Table.read(table_to_read, \
             names = ["gal_id","zz","mag1","mag_error1","mag1_exist","re_pix1",\
                      "re_pix_error1","re1_exist","re_kpc1","re_kpc_error1","n1","n_error1","n1_exist",\
                      "ar1","ar_error1","ar1_exist","pa1","pa_error1","pa1_exist","mag2",\
                      "mag_error2","mag2_exist","re_pix2","re_pix_error2","re2_exist","re_kpc2","re_kpc_error2",\
                      "n2","n_error2","n2_exist","ar2","ar_error2","ar2_exist","pa2",\
                      "pa_error2","pa2_exist","index","flag_obj_detected"],\
             format="ascii.commented_header")
        filter_gal = tt["gal_id"] == gal
        index = tt["index"][filter_gal]
        index = np.float(index[0]) #in order to be able to do np.isnan afterwards
        if np.isnan(index): #if the there is no best fit, return nothing
            return()
        else:
            index = int(index)
            return()

def plotting_data_original(gal,band,filename,zz,logmass):
    subfig.add_label(0.40,0.95,"{0}"            .format(str(filename[:-17])),relative=True,color="black")
    subfig.add_label(0.80,0.85,"{0}"            .format(str(gal))           ,relative=True,color="black")
    subfig.add_label(0.80,0.80,"Filter: {0}"    .format(str(band))          ,relative=True,color="black")
    subfig.add_label(0.80,0.75,"z = {0:4.2f}"   .format(zz)                 ,relative=True,color="black")
    subfig.add_label(0.80,0.70,"mass = {0:6.2e}".format(10.**logmass)         ,relative=True,color="black")
    
def plotting_data_model(re_kpc,re_pix,nn):
    subfig.add_label(0.80,0.90,"re_kpc = {0:4.2f}".format(re_kpc),relative=True,color='black')
    subfig.add_label(0.80,0.85,"re_pix = {0:4.2f}".format(re_pix),relative=True,color='black')
    subfig.add_label(0.80,0.80,"n      = {0:4.2f}".format(nn)    ,relative=True,color='black')

def plotting_colorbar():
    subfig.add_colorbar()
    subfig.colorbar.show(box=[0.4,0.05,0.2,0.025],box_orientation='horizontal',axis_label_text=r"Surf. brightness / mag arcsec$^{2}$") #box=[xmin, ymin, dx, dy]

def plotting_scalebar():
    subfig.add_scalebar(1./3600.)
    subfig.scalebar.set_color('black')
    subfig.scalebar.set_label('1 arcsec')

#it needs the pixel scale in the header
#def plotting_seeing(psf_fwhm):
#    psf_fwhm = psf_fwhm / 3600.
#    subfig.add_beam(major=psf_fwhm,minor=psf_fwhm,angle=0.,corner="bottom left")
#def plotting_seeing(psf_fwhm):
#    ra, dec = subfig.pixel2world(20.,20.)
#    subfig.show_circles(ra, dec, psf_fwhm, edgecolor='#FFA500')


#CONSTANTS
catalog_path = "/mnt/disk1/fb/gama_compacts/galaxy_selection/"
catalog      = "galaxy_selection_in_the_i_band.txt"
bands = np.array(["g","r","i"])
max_seeing_values = np.array([1.,0.7,1.]) #arcsec
pix_scale = 0.21 #arcsec/pix in KiDS
empty_images = np.array(["196099","377359","372766"])
#plotting
min_SB = 27.
max_SB = 18.
plotting_pos = [[0.00,0.66,0.33,0.33],
                [0.33,0.66,0.33,0.33],
                [0.66,0.66,0.33,0.33],
                [0.00,0.33,0.33,0.33],
                [0.33,0.33,0.33,0.33],
                [0.66,0.33,0.33,0.33],
                [0.00,0.00,0.33,0.33],
                [0.33,0.00,0.33,0.33],
                [0.66,0.00,0.33,0.33],
               ]
# plotting_pos= [[0.00,0.75,0.25,0.25],
#                [0.25,0.75,0.25,0.25],
#                [0.50,0.75,0.25,0.25],
#                [0.75,0.75,0.25,0.25],
#                [0.00,0.50,0.25,0.25],
#                [0.25,0.50,0.25,0.25],
#                [0.50,0.50,0.25,0.25],
#                [0.75,0.50,0.25,0.25],
#                [0.00,0.25,0.25,0.25],
#                [0.25,0.25,0.25,0.25],
#                [0.50,0.25,0.25,0.25],
#                [0.75,0.25,0.25,0.25],
#                [0.00,0.00,0.25,0.25],
#                [0.25,0.00,0.25,0.25],
#                [0.50,0.00,0.25,0.25],
#                [0.75,0.00,0.25,0.25],
#               ]

#reading the catalog
tt = Table.read(catalog_path+catalog, names=("cataid","zz","logmass","re"), format='ascii.commented_header')
cataid  = tt['cataid']
zz      = tt['zz']
logmass = tt["logmass"]
re_cat  = tt["re"]

fig = plt.figure( figsize = (15,15) )

for cont_gal,gal in enumerate(cataid):
    cont = 0
    print(gal)
    if str(gal) in empty_images:
        continue
    for cont_band,band in enumerate(bands):
        #printing the original image
        filename = getting_filename_from_pattern_running_ls(gal,band,"original")
        if os.path.isfile(filename):
            img = fits.open(filename)
            header_original = img[0].header
            subfig = aplpy.FITSFigure(img, figure = fig, subplot=plotting_pos[cont], origin='lower')
            re_pix,re_kpc,nn,ar,pa = get_struct_param_best_fit(gal,band,"single")
            plotting_options()
            subfig.show_colorscale(vmin=max_SB,vmax=min_SB,cmap='cubehelix')
            plotting_data_original(gal,band,filename,zz[cont_gal],logmass[cont_gal])
            plotting_scalebar()
            #if cont_band == 0: plotting_seeing(max_seeing_values[cont_band])
        cont = cont+1
        #printing the model image
        filename = getting_filename_from_pattern_running_ls(gal,band,"model")
        if os.path.isfile(filename):
            img = fits.open(filename)
            img[0].header = header_original #in order to be able to print the ellipse with the effective radius
            subfig = aplpy.FITSFigure(img, figure = fig, subplot=plotting_pos[cont], origin='lower')
            printing_effective_radius(img[0].data,pix_scale,re_pix*pix_scale,ar,pa) #only in the original image, as it is the only one with "good" header
            plotting_options()
            subfig.show_colorscale(vmin=max_SB,vmax=min_SB,cmap='cubehelix')
            plotting_data_model(re_kpc,re_pix,nn)
            #if cont_band == 0: plotting_seeing(max_seeing_values[cont_band])
        cont = cont+1
        #printing the residual image
        filename = getting_filename_from_pattern_running_ls(gal,band,"residual")
        if os.path.isfile(filename):
            img = fits.open(filename)
            subfig = aplpy.FITSFigure(img, figure = fig, subplot=plotting_pos[cont], origin='lower')
            plotting_options()
            subfig.show_colorscale(vmin=max_SB,vmax=min_SB,cmap='cubehelix')
            #if cont_band == 0: plotting_seeing(max_seeing_values[cont_band])

        #=============================
        #only for the last image, create the colorbar
        if band == bands[-1]:
            plotting_colorbar()
        #the colorbar changes the subimage layout, I need to print it again
        if band == bands[-1]:
            #printing the residual image
            filename = getting_filename_from_pattern_running_ls(gal,band,"residual")
            if os.path.isfile(filename):
                img = fits.open(filename)
                subfig = aplpy.FITSFigure(img, figure = fig, subplot=plotting_pos[cont], origin='lower')
                #print(plotting_pos[cont],band)
                plotting_options()
                subfig.show_colorscale(vmin=max_SB,vmax=min_SB,cmap='cubehelix')
        #=============================

        cont = cont+1

    plt.savefig(str(gal)+"_galfit_fits_single.pdf", format='pdf', dpi=100)
    plt.clf()
    
import aplpy
import numpy as np
from astropy.io import fits
from astropy.table import Table
import matplotlib.pyplot as plt
from matplotlib import gridspec
import os


def create_sb_image(path,filename,extension,zp,pix_scale):
    img = fits.open(path+filename)
    image  = img[extension].data
    header = img[extension].header
    sb = -2.5*np.log10(image)+zp+5.*np.log10(pix_scale)

    if extension == 1:
        img_type = "original"
    elif extension == 2:
        img_type = "model"
    elif extension == 3:
        img_type = "residual"

    fits.writeto(filename[:-5]+"_"+img_type+"_sb.fits",sb,header,clobber=True)

def get_table_file(band,type_of_analysis):
    path = "../analysis_filter_"+band+"/"
    if type_of_analysis == "single":
        file_to_read = "gama_kids_"+band+"_single_sersic_fits.cat"
    else:
        file_to_read = "gama_kids_"+band+"_bulge_disk_decompositions.cat"
    return(path+file_to_read)

def get_filename(gal,band):
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
        if np.isnan(index): #if the there is no best fit, return nothing
            pass
        else:
            index = int(index)
            star_number = index / 2; star_number = int(star_number)
            mask_number = index % 2; mask_number = int(mask_number)
            ini_cond_number = 0    ; ini_cond_number = int(ini_cond_number)
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
        if np.isnan(index): #if the there is no best fit, return nothing
            pass
        else:
            index = int(index)
        star_number = index / 4          ; star_number = int(star_number)
        mask_number     = (index / 2) / 2; mask_number = int(mask_number)
        ini_cond_number = (index / 2) % 2; ini_cond_number = int(ini_cond_number)

    if np.isnan(index): #if the there is no best fit, return nothing
        return("")
    else: 
        return( str(gal)+"_"+stars[star_number]+"_"+masks[mask_number]+"_"+str(ini_cond_number)+".fits" )

def get_path(gal,band,type_of_analysis,number_cpus):
    base_path = "/mnt/disk1/fb/gama_compacts/"

    #band
    if band == "g":
        path = base_path+"analysis_filter_g/"
    elif band == "r":
        path = base_path+"analysis_filter_r/"
    elif band == "i":
        path = base_path+"analysis_filter_i/"
    elif band == "Z":
        path = base_path+"viking_z_band_images/"

    #type of analysis
    if type_of_analysis == "single":
        path = path + "single_sersic_fits/"
    elif type_of_analysis == "bulge_disk":
        path = path + "bulge_disk_decompositions/"

    for ii in range(1,number_cpus+1):
        total_path = path+"cpu"+str(ii)+"/"+str(gal)+"/"
        if os.path.exists(total_path):
            print(path,ii,gal)
            return(total_path)

def select_conditions(band,camera,survey,special):

    zp_kids = 0.
    pix_scale_kids = 0.21
    zp_viking_z = 30.
    pix_scale_viking_z = 0.339

    if   band == "g":
        zp = zp_kids
        pix_scale = pix_scale_kids
    elif band == "r":
        zp = zp_kids
        pix_scale = pix_scale_kids
    elif band == "i":
        zp = zp_kids
        pix_scale = pix_scale_kids
    elif band == "Z":
        zp = zp_viking_z
        pix_scale = pix_scale_viking_z

    return( (zp,pix_scale,None) )


#CONSTANTS
catalog_path = "/mnt/disk1/fb/gama_compacts/galaxy_selection/"
catalog      = "galaxy_selection_in_the_i_band.txt"
bands = ["g","r","i"]
type_of_analysis = "single"
number_cpus = 3
#see also the constants in get_filename


#reading the catalog
tt = Table.read(catalog_path+catalog, names=("cataid","zz","logmass","re"), format='ascii.commented_header')
cataid  = tt['cataid']
zz      = tt['zz']
logmass = tt["logmass"]
re_cat  = tt["re"]


for cont_gal,gal in enumerate(cataid):
    for cont_band,band in enumerate(bands):
        path = get_path(gal,band,type_of_analysis,number_cpus)

        zp,pix_scale,nothing = select_conditions(band,None,None,None)

        filename = get_filename(gal,band)
        print(filename)

        if os.path.isfile(path+filename):
            create_sb_image(path,filename,1,zp,pix_scale)
            create_sb_image(path,filename,2,zp,pix_scale)
            create_sb_image(path,filename,3,zp,pix_scale)

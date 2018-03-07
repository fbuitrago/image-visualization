import aplpy
import numpy as np
from astropy.io import fits
from astropy.table import Table
import os
import aplpy
import matplotlib.pyplot as plt
from matplotlib import rc,rcParams
from astropy.cosmology import FlatLambdaCDM
import pdb


def create_sb_image(gal,path,band,zp,pix_scale):
    img = fits.open(path+str(gal)+"-kids-"+band+".fits")
    image  = img[0].data
    header = img[0].header
    sb = -2.5*np.log10(image)+zp+5.*np.log10(pix_scale)
    fits.writeto(path+str(gal)+"_kids_"+band+"_sb.fits",sb,header,clobber=True)
    
plotting_pos=[[0./6.,5./6.,1./6.,1./6.],
              [1./6.,5./6.,1./6.,1./6.],
              [2./6.,5./6.,1./6.,1./6.],
              [3./6.,5./6.,1./6.,1./6.],
              [4./6.,5./6.,1./6.,1./6.],
              [5./6.,5./6.,1./6.,1./6.],
              [0./6.,4./6.,1./6.,1./6.],
              [1./6.,4./6.,1./6.,1./6.],
              [2./6.,4./6.,1./6.,1./6.],
              [3./6.,4./6.,1./6.,1./6.],
              [4./6.,4./6.,1./6.,1./6.],
              [5./6.,4./6.,1./6.,1./6.],
              [0./6.,3./6.,1./6.,1./6.],
              [1./6.,3./6.,1./6.,1./6.],
              [2./6.,3./6.,1./6.,1./6.],
              [3./6.,3./6.,1./6.,1./6.],
              [4./6.,3./6.,1./6.,1./6.],
              [5./6.,3./6.,1./6.,1./6.],
              [0./6.,2./6.,1./6.,1./6.],
              [1./6.,2./6.,1./6.,1./6.],
              [2./6.,2./6.,1./6.,1./6.],
              [3./6.,2./6.,1./6.,1./6.],
              [4./6.,2./6.,1./6.,1./6.],
              [5./6.,2./6.,1./6.,1./6.],
              [0./6.,1./6.,1./6.,1./6.],
              [1./6.,1./6.,1./6.,1./6.],
              [2./6.,1./6.,1./6.,1./6.],
              [3./6.,1./6.,1./6.,1./6.],
              [4./6.,1./6.,1./6.,1./6.],
              [5./6.,1./6.,1./6.,1./6.],
              [0./6.,0./6.,1./6.,1./6.],
              [1./6.,0./6.,1./6.,1./6.],
              [2./6.,0./6.,1./6.,1./6.],
              [3./6.,0./6.,1./6.,1./6.],
              [4./6.,0./6.,1./6.,1./6.],
              [5./6.,0./6.,1./6.,1./6.]
              ]
max_g = 18.
max_r = 18.
max_i = 18.
min_g = 27.
min_r = 27.
min_i = 27.

#paths
images_gband_path = "/mnt/disk1/fb/gama_compacts/galaxies_all_selections_lee_18jul17/analysis_filter_g/galaxy_images/"
images_rband_path = "/mnt/disk1/fb/gama_compacts/galaxies_all_selections_lee_18jul17/analysis_filter_r/galaxy_images/"
images_iband_path = "/mnt/disk1/fb/gama_compacts/galaxies_all_selections_lee_18jul17/analysis_filter_i/galaxy_images/"
catalogs_path = "./"
#constants
min_value = -0.0001
sample       = catalogs_path+"list_good_gals_with_radii.cat"
zp = 30.
pix_scale = 0.21

rc('axes', linewidth=2)
rc('font', weight='bold')

#reading the catalog
#tt = Table.read(sample, names=("cataid","re_kpc_g","good_analysis_g","re_kpc_r","good_analysis_r","re_kpc_i","good_analysis_i","re_kpc_z","good_analysis_z"), format='ascii.commented_header')
#cataid  = tt['cataid']
#re_kpc_g      = tt['re_kpc_g']
#re_kpc_r      = tt['re_kpc_r']
#re_kpc_i      = tt['re_kpc_i']
#re_kpc_z      = tt['re_kpc_z']
#===
#cataid = np.loadtxt(sample, dtype = int)
#cataid = cataid.astype(str) 
#===
tt = Table.read("../galaxy_selection/galaxies_for_Ignacio_with_ra_dec_z_mass_selection35objects.cat",
                     names = ["cataid","ra","dec","z_tonry","logmstar"],
                     format="ascii.commented_header")
cataid    = tt['cataid']
zz        = tt["z_tonry"]
logmstar  = tt["logmstar"]

#I order the galaxies according to their redshift
zz_ordered_indices = np.argsort(zz)
cataid   = cataid  [zz_ordered_indices]
zz       = zz      [zz_ordered_indices]
logmstar = logmstar[zz_ordered_indices]

#moving distances from kpc to arcsec
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

fig = plt.figure(figsize=(20,20))

for cont_gal,gal in enumerate(cataid):

    print(gal)

    create_sb_image(gal,images_gband_path,"g",zp,pix_scale)
    create_sb_image(gal,images_rband_path,"r",zp,pix_scale)
    create_sb_image(gal,images_iband_path,"i",zp,pix_scale)

    #aplpy.make_rgb_cube([images_gband_path+str(gal)+'_kids_g.fits', images_rband_path+str(gal)+'_kids_r.fits', images_iband_path+str(gal)+'_kids_i.fits'], str(gal)+'_cube.fits')
    aplpy.make_rgb_cube([images_iband_path+str(gal)+'_kids_i_sb.fits', images_rband_path+str(gal)+'_kids_r_sb.fits', images_gband_path+str(gal)+'_kids_g_sb.fits'], str(gal)+'_cube.fits')

    img_g = fits.open(images_gband_path+str(gal)+'_kids_g_sb.fits')
    img_r = fits.open(images_rband_path+str(gal)+'_kids_r_sb.fits')
    img_i = fits.open(images_iband_path+str(gal)+'_kids_i_sb.fits')

    aplpy.make_rgb_image(str(gal)+'_cube.fits','./color_images/'+str(gal)+'_rgb.png',stretch_r='linear',stretch_g='linear',stretch_b='linear',vmin_r=min_i,vmax_r=max_i,vmin_g=min_r,vmax_g=max_r,vmin_b=min_g,vmax_b=max_g)

    os.remove(images_gband_path+str(gal)+'_kids_g_sb.fits')
    os.remove(images_rband_path+str(gal)+'_kids_r_sb.fits')
    os.remove(images_iband_path+str(gal)+'_kids_i_sb.fits')
    os.remove(str(gal)+'_cube.fits')
    
    img = fits.open(str(gal)+"_cube_2d.fits")
    subfig = aplpy.FITSFigure(img, figure=fig, subplot=plotting_pos[cont_gal])
    subfig.set_axis_labels_size(20)
    subfig.set_axis_labels_weight("bold")
    subfig.set_tick_labels_size(20)
    subfig.add_scalebar(10./3600)
    kpc_per_arcsec = 1./cosmo.arcsec_per_kpc_proper(zz[cont_gal])
    
    #axes===================
    #subfig.set_axis_labels_size(20)
    #subfig.set_axis_labels_weight("bold")
    #subfig.set_tick_labels_size(20)
    subfig.axis_labels.hide()
    subfig.tick_labels.hide()
    subfig.ticks.hide()
    ##========================
    
    subfig.scalebar.set_label('10 arcsec = {0:4.2f} kpc'.format(10*kpc_per_arcsec.value))
    subfig.scalebar.set_linewidth(2)
    subfig.scalebar.set_font(size=10, weight='bold')
    subfig.scalebar.set_color('white')
    subfig.show_rgb("./color_images/"+str(gal)+"_rgb.png")
    
    subfig.add_label(0.25,0.90,str(gal)                                 ,relative=True,color='white',weight="bold",size=14)
    subfig.add_label(0.25,0.85,"z = {0:4.2f}".format(zz      [cont_gal]),relative=True,color='white',weight="bold",size=10)
    subfig.add_label(0.42,0.80,r"log$_{10}$(M$_{\rm stellar}$/M$_{\odot}$) = "+"{0:5.2f}".format(logmstar[cont_gal]),relative=True,color='white',weight="bold",size=10)    
    
    os.remove(str(gal)+"_cube_2d.fits")
    os.remove('./color_images/'+str(gal)+'_rgb.png')
    
plt.savefig("./color_images/mosaic_rgb.pdf", format='pdf', dpi=100)
plt.clf()
plt.close(fig)


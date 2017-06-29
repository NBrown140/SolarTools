import os

import arcpy, SolarTools

############# Functions ################

def get_file_name(file_path):
    """Extract the file name without the file extension from a file path."""
    basename = os.path.basename(file_path)
    return os.path.splitext(basename)[0]
  
  
#########################################

# Define directories
# Clip polygons list
clip_list = [path_to_poly_1, path_to_poly_2]
# DSM path
path_dsm = r'path/to/dsm.tif'
dsm_name = get_file_name(path_dsm)
# Output directory
path_out_dir = r'path/to/outdir'


# Loop through all clipping polygons
for clip_poly in clip_list:
  # Execute clip
  print('Executing clip...')
  clip_name = get_file_name(clip_poly)
  path_dsm_clip = os.path.join(path_out_dir,dsm_name+'_'+clip_name+'.tif')
  arcpy.Clip_management (path_dsm, clip_poly, path_dsm_clip, "ClippingGeometry")
  # Calculate solar radiation
  print('Calculating Solar Radiation for: {path_dsm_clip}'.format())
  outGlobalRad = arcpy.sa.AreaSolarRadiation(path_dsm_clip, TimeWholeYear({2017}), "INTERVAL")
  # Aply SolarTools- Calibrate solar potential
  

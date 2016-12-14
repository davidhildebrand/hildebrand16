from ij import IJ
from ij.io import FileSaver
import os
from os import path
import string

source_path = "Correspondence_ZBrain/ZBtoEM/ZBdata_warped/"
reformat_path = "Correspondence_ZBrain/ZBtoEM/ZBdata_imseq/"

for filename in os.listdir(source_path):
  if filename.endswith(".tif"):
    print "Processing file:", filename

    ## Open image from file
    imp = IJ.openImage(os.path.join(source_path, filename))
    if imp is None:
      print "Could not open image from file:", filename
      continue
    ## Set appropriate resolution (matching BigWarp settings here)
    IJ.run(imp, "Properties...", "channels=1 slices=811 frames=1 unit=um pixel_width=0.6 pixel_height=0.6 voxel_depth=1.2")
    ## Make sure full 16-bit or 8-bit range is used
    IJ.setMinAndMax(imp, 0, 65535) # ZBrain/ZBB data stacks
    imp.setDefault16bitRange(16)
    #IJ.setMinAndMax(imp, 0, 255); # ZBrain masks/labels
    #imp.setDefault16bitRange(8);
    ## Convert to 8-bit (optional, if desired)
    IJ.run(imp, "8-bit", "")
    ## If converted, make sure full 8-bit range is used before saving
    IJ.setMinAndMax(imp, 0, 255);
    imp.setDefault16bitRange(8);
    
    fs = FileSaver(imp)
    if not path.exists(reformat_path):
      os.makedirs(reformat_path)
    if path.exists(reformat_path) and path.isdir(reformat_path):
      dirpath = path.join(reformat_path, string.replace(filename, '_Warped.tif', ''))
      if not path.exists(dirpath):
        os.makedirs(dirpath)
      IJ.run(imp, "Image Sequence... ", "name= format=PNG start=1 digits=3 save=[" + dirpath + "/000.png]");
    else:
      print "Folder does not exist", reformat_path

    imp.close()
  else:  
    print "Ignoring non-TIFF", filename
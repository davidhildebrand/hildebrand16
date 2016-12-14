from ij import IJ
from ij.io import FileSaver
import os
from os import path
import string

source_path = "Correspondence_ZBrain/ZB/data/"
reformat_path = "Correspondence_ZBrain/ZB/data_ready/"

for filename in os.listdir(source_path):
  if filename.endswith(".tif"):
    print "Processing file:", filename

    ## Open image from file
    imp = IJ.openImage(os.path.join(source_path, filename))
    if imp is None:
      print "Could not open image from file:", filename
      continue
    ## Make sure full 16-bit or 8-bit range is used
    IJ.setMinAndMax(imp, 0, 65535) # ZBrain/ZBB data stacks
    imp.setDefault16bitRange(16)
    #IJ.setMinAndMax(imp, 0, 255); # ZBrain masks/labels
    #imp.setDefault16bitRange(8);
    ## Convert to 8-bit (optional, if desired)
    #IJ.run(imp, "8-bit", "")
    ## Set appropriate resolution (matching original here)
    IJ.run(imp, "Properties...", "channels=1 slices=138 frames=1 unit=um pixel_width=0.798 pixel_height=0.798 voxel_depth=2.0")
    
    fs = FileSaver(imp)
    if path.exists(reformat_path) and path.isdir(reformat_path):
      filepath = path.join(reformat_path, string.replace(filename, '_BeforeWarp', '_ReadyWarp'))
      if path.exists(filepath):
        print "File exists, NOT overwriting."
      elif fs.saveAsTiff(filepath):
        print "File saved successfully: ", filepath
    else:
      print "Folder does not exist", reformat_path

    imp.close()
  else:  
    print "Ignoring non-TIFF", filename
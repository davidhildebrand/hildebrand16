#!/usr/bin/python
__author__ = 'David Grant Colburn Hildebrand'

import argparse
import math
import numpy
import os
from PIL import Image


def directory_check(path):
    if not os.path.isdir(path):
        err_msg = "path is not a directory (%s)"
        raise argparse.ArgumentTypeError(err_msg)
    return path

def file_check(path):
    if not os.path.isfile(path):
        err_msg = "path is not a file (%s)"
        raise argparse.ArgumentTypeError(err_msg)
    return path

# parse command line options
parser = argparse.ArgumentParser()
parser.add_argument(
    '-s', '--source', type=file_check, required=True,
    help="Path to a source single image. [required]")
parser.add_argument(
    '-d', '--dest', type=directory_check, required=True,
    help="Path to a destination directory where PNG tiles for CATMAID" \
         "tile source format 4 (section/scale/row_col.png) will be" \
         "output. [required]")
parser.add_argument(
    '-z', '--sect', type=int, required=True,
    help="Index associated with the section for which >0-scale pyramid" \
         "tiles are to be generated. [required]")
parser.add_argument(
    '-t', '--size', type=int, required=False, default=1024,
    help="Tile size to be used. [optional, default 1024]")
opts = parser.parse_args()


trans_intens = 255


source_path = opts.source
dest_path = opts.dest
sect = opts.sect
tile_size = opts.size

ext = '.png'

# open tile
try:
    sect_image = Image.open(source_path)
    if sect_image.mode != 'L':
        print "WARNING: source image is not mode L (" + sect_path + ")"
except IOError:
    print "ERROR: could not open source image file (" + source_path + ")"
    exit()

# detect image size
w = sect_image.size[0]
h = sect_image.size[1]
#print "DEBUG: source image is " + str(w) + " px x " + str(h) + " px"

# calculate number of tile rows and columns based on size
#print "DEBUG: tile size is " + str(tile_size)
colnum = int(math.ceil(float(w) / tile_size))
rownum = int(math.ceil(float(h) / tile_size))
#print "DEBUG: tiles " + str(colnum) + " cols x " + str(rownum) + " rows"

# set 0-scale, pyramids can be generated later with a different script
scale = 0

# determine number of rows/cols
rows = range(rownum)
cols = range(colnum)

save_path = os.path.join(dest_path, str(sect), str(scale))
if not os.path.exists(save_path):
    os.makedirs(save_path)

tile_info = sect_image.info
if not 'transparency' in tile_info.keys():
    # set the tRNS flag to intensity value trans_intens in the image info
    tile_info['transparency'] = trans_intens
    # clip intensity histogram to make room for tRNS
    sect_arr = numpy.array(sect_image)
    sect_arr = numpy.clip(sect_arr, 0, trans_intens - 1)
    sect_data = Image.fromarray(sect_arr)
elif tile_info['transparency'] != trans_intens:
    # set the tRNS flag to intensity value trans_intens in the image info
    tile_info['transparency'] = trans_intens
    # TODO: scale the rest of the intensities rather than clipping
    sect_arr = numpy.array(sect_image)
    sect_arr = numpy.clip(sect_arr, 0, trans_intens - 1)
    sect_data = Image.fromarray(sect_arr)
    print "WARNING: old and new tRNS values not the same, clipping performed"

tile_mode = 'L'
#print "DEBUG: padded image size " + str(colnum*tile_size) + " x " + str(rownum*tile_size)
data = Image.new(tile_mode, (colnum * tile_size, rownum * tile_size))
data.paste(trans_intens)
data.paste(sect_data, (0, 0, w, h))
#data.load()
#data.save(os.path.join(save_path, "data.png"))

for row in rows:
    top = row * tile_size
    for col in cols:
        left = col * tile_size

        tile_path = os.path.join(save_path, str(row) + "_" + str(col) + ext)

        # allocate tile
        tile = Image.new(tile_mode, (tile_size, tile_size))
        #print "DEBUG: left="+str(left)+" top="+str(top)+" right="+str(left+tile_size)+" bot="+str(top+tile_size)
        #print "DEBUG: sect_mode="+sect_image.mode+" data_mode="+data.mode+" tile_mode="+tile_mode
        # copy data from region of padded section image
        data_crop = data.crop((left, top, left + tile_size, top + tile_size))
        tile.paste(data_crop)

        # save tile with info including transparency flag
        tile.save(tile_path, **tile_info)

        print "tile row " + str(row) + " col " + str(col) + " saved with tRNS " + \
              str(trans_intens)

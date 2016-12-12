#!/usr/bin/python
__author__ = 'David Grant Colburn Hildebrand'

import argparse
import math
import os
from PIL import Image


def directory(path):
    if not os.path.isdir(path):
        err_msg = "path is not a directory (%s)"
        raise argparse.ArgumentTypeError(err_msg)
    return path


# parse command line options
parser = argparse.ArgumentParser()
parser.add_argument(
    '-s', '--source', type=directory, required=True,
    help="Path to a source directory containing a subdirectory corresponding" \
         "to the section number that contains 0-scale PNG tiles in CATMAID" \
         "tile source format 4 (section/scale/row_col.png). [required]")
parser.add_argument(
    '-d', '--dest', type=directory, required=True,
    help="Path to a destination directory. [required]")
parser.add_argument(
    '-z', '--sect', type=int, required=True,
    help="Index associated with the section for which >0-scale pyramid" \
         "tiles are to be generated. [required]")
opts = parser.parse_args()

source_path = opts.source
dest_path = opts.dest
sect = opts.sect
sect_path = os.path.join(source_path, str(sect))
s0_path = os.path.join(sect_path, '0')
ext = '.png'

# automatically determine tile information
#   and montage width (cols) and height (rows)
s0_tiles = [ f for f in os.listdir(s0_path) if os.path.isfile(os.path.join(s0_path, f)) \
                                            and os.path.splitext(f)[1].lower() == ext ]

# automatically detect tile size
s0_tile = Image.open(os.path.join(s0_path, s0_tiles[0]))
if s0_tile.size[0] == s0_tile.size[1]:
    tile_size = s0_tile.size[0]
else:
    print "ERROR: 0-scale tiles are not square"
    exit()
s0_tile_mode = s0_tile.mode
s0_tile_info = s0_tile.info
if s0_tile_mode != 'RGBA':
    print "ERROR: source tiles are not RGBA mode"
    exit()

s0_rows = max([ int(os.path.splitext(f)[0].split('_')[0]) for f in s0_tiles ]) + 1
s0_cols = max([ int(os.path.splitext(f)[0].split('_')[1]) for f in s0_tiles ]) + 1
s0_height = s0_rows * tile_size 
s0_width = s0_cols * tile_size

if len(s0_tiles) != (s0_rows * s0_cols):
    print "ERROR: expected number of 0-scale tiles not found" \
          "       ... expected " + str(s0_rows * s0_cols) + " but found " + len(s0_tiles)
    exit()

print "scale 0"
print "    path " + s0_path

# # create empty tile in case tiles do not exist or cannot be opened
# empty_tile = Image.new(tile_mode, (tile_size, tile_size))
# if 'transparency' in tile_info.keys():
#     empty_tile.paste(tile_info['transparency'], (0, 0, tile_size, tile_size))

# set 0-scale parameters
scale = 0

# determine number of rows/cols
rows = range(int(s0_height / tile_size))
cols = range(int(s0_width / tile_size))

save_path = os.path.join(dest_path, str(sect), str(scale))
if not os.path.exists(save_path):
    os.makedirs(save_path)

print "    rows " + str(len(rows)) + " cols " + str(len(cols))
for row in rows:
    for col in cols:
        # open tile
        try:
            tile_path = os.path.join(sect_path, str(scale), \
                                     str(row) + "_" + str(col) + ext)
            tile = Image.open(tile_path)
            tile.load() # required for tile.split()
        except IOError:
            print "ERROR: could not open source tile (" + tile_path + ")"
            exit()

        new_tile = Image.new("RGB", tile.size, (250, 250, 250))
        new_tile.paste(tile, mask=tile.split()[-1]) # -1 is the alpha channel
        new_tile = new_tile.convert('L')

        # save re-composed sect_image
        new_tile.save(os.path.join(save_path, str(row) + "_" + str(col) + ext))
        print "       row " + str(row) + " col " + str(col) + " saved"


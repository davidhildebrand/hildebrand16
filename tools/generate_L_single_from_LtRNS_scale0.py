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
tile_mode = s0_tile.mode
tile_info = s0_tile.info
if tile_mode != 'L':
    print "ERROR: only luminance mode without transparency is supported"
    exit()

s0_rows = max([ int(os.path.splitext(f)[0].split('_')[0]) for f in s0_tiles ]) + 1
s0_cols = max([ int(os.path.splitext(f)[0].split('_')[1]) for f in s0_tiles ]) + 1
s0_height = s0_rows * tile_size 
s0_width = s0_cols * tile_size

if len(s0_tiles) != (s0_rows * s0_cols):
    print "ERROR: expected number of 0-scale tiles not found" \
          "       ... expected " + str(s0_rows * s0_cols) + " but found " + len(s0_tiles)
    exit()

# set 0-scale parameters
scale = 0
top = 0
left = 0

# determine number of rows/cols
rows = range(int(s0_height / tile_size))
cols = range(int(s0_width / tile_size))

sect_image = Image.new(tile_mode, (s0_width, s0_height))

print "sect " + str(sect)
print "    path " + s0_path
print "    w" + str(s0_width) + "px x h" + str(s0_height) + "px"
print "    rows " + str(len(rows)) + " cols " + str(len(cols))
for row in rows:
    top = row * tile_size
    for col in cols:
        left = col * tile_size

        # open tile
        try:
            tile_path = os.path.join(sect_path, str(scale), \
                                     str(row) + "_" + str(col) + ext)
            tile = Image.open(tile_path)
            if 'transparency' in tile.info.keys():
                print "WARNING: tRNS transparency is not preserved"
        except IOError:
            print "ERROR: could not open tile (" + tile_path + ")"
            exit()

        # copy tile into section image
        sect_image.paste(tile, (left, top))

# save re-composed sect_image
save_path = os.path.join(dest_path, '{0:05d}'.format(sect) + ext)

# preserve tile info
sect_image_info = tile_info
if 'transparency' in sect_image_info.keys():
    sect_image_info.pop('transparency', None)

sect_image.save(save_path, **sect_image_info)


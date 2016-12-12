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
    '-p', '--path', type=directory, required=True,
    help="Path to a directory containing a subdirectory corresponding" \
         "to the section number that contains 0-scale PNG tiles in CATMAID" \
         "tile source format 4 (section/scale/row_col.png). [required]")
parser.add_argument(
    '-z', '--sect', type=int, required=True,
    help="Index associated with the section for which >0-scale pyramid" \
         "tiles are to be generated. [required]")
parser.add_argument(
    '-v', '--verbose', action='store_true',
    help="Verbose output. [optional]" )
opts = parser.parse_args()

path = opts.path
sect = opts.sect
verb = opts.verbose
sect_path = os.path.join(path, str(sect))
s0_path = os.path.join(sect_path, '0')
ext = '.png'


background = 250


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
    print "ERROR: support currently exists for luminance mode"
    exit()
if 'transparency' in tile_info.keys():
    print "WARNING: tRNS transparency found for this tile set and ignored"

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
print "    w" + str(s0_width) + "px x h" + str(s0_height) + "px"

# allocate empty tile in case tiles do not exist or cannot be opened
empty_tile = Image.new(tile_mode, (tile_size, tile_size))
empty_tile.paste(background, (0, 0, tile_size, tile_size))

# set initial 1-scale parameters
scale = 1
height = math.ceil(s0_height / 2)
width = math.ceil(s0_width / 2)

while height >= math.ceil(tile_size / 2) and width >= math.ceil(tile_size / 2): 
    print "  scale " + str(scale)
    scale_path = os.path.join(sect_path, str(scale))
    if not os.path.isdir(scale_path):
        print "    path not found"
        print "      created directory " + scale_path
        os.makedirs(scale_path, 0755)

    print "    path " + scale_path
    print "    w" + str(width) + "px x h" + str(height) + "px"

    # initialize for this scale
    top = 0
    left = 0
    # determine number of rows/cols for this scale
    rows = range(int(math.ceil(height / tile_size)))
    cols = range(int(math.ceil(width / tile_size)))
    print "    rows " + str(len(rows)) + " cols " + str(len(cols))
    for row in rows:
        top = row * tile_size
        for col in cols:
            left = col * tile_size

            # open tiles for previous scale
            try:
                tl_path = os.path.join(sect_path, str(scale - 1), \
                                       str(2 * row) + "_" + str(2 * col) + ext)
                if verb:
                    print "tl_path" + tl_path
                tl = Image.open(tl_path)
            except IOError:
                print "        WARNING: could not open (" + \
                      os.path.join(str(scale - 1), \
                                   str(2 * row) + "_" + str(2 * col) + ext) + \
                      ")... replaced with empty tile"
                if verb:
                    print "tl empty mode " + tile_mode + " size " + tile_size
                tl = Image.new(tile_mode, (tile_size, tile_size))
                tl.paste(empty_tile)

            try:
                tr_path = os.path.join(sect_path, str(scale - 1), \
                                       str(2 * row) + "_" + str((2 * col) + 1) + ext)
                tr = Image.open(tr_path)
            except IOError:
                print "        WARNING: could not open (" + \
                      os.path.join(str(scale - 1), \
                                   str(2 * row) + "_" + str((2 * col) + 1) + ext) + \
                      ")... replaced with empty tile"
                tr = Image.new(tile_mode, (tile_size, tile_size))
                tr.paste(empty_tile)

            try:
                bl_path = os.path.join(sect_path, str(scale - 1), \
                                       str((2 * row) + 1) + "_" + str(2 * col) + ext)
                bl = Image.open(bl_path)
            except IOError:
                print "        WARNING: could not open (" + \
                      os.path.join(str(scale - 1), \
                                   str((2 * row) + 1) + "_" + str(2 * col) + ext) + \
                      ")... replaced with empty tile"
                bl = Image.new(tile_mode, (tile_size, tile_size))
                bl.paste(empty_tile)

            try:
                br_path = os.path.join(sect_path, str(scale - 1), \
                                       str((2 * row) + 1) + "_" + str((2 * col) + 1) + ext)
                br = Image.open(br_path)
            except IOError:
                print "        WARNING: could not open (" + \
                      os.path.join(str(scale - 1), \
                                   str((2 * row) + 1) + "_" + str((2 * col) + 1) + ext) + \
                      ")... replaced with empty tile"
                br = Image.new(tile_mode, (tile_size, tile_size))
                br.paste(empty_tile)

            # allocate 2x2 tile array that will be downsampled
            tile = Image.new(tile_mode, (2 * tile_size, 2 * tile_size))

            # copy opened tiles into 2x2 tile array
            tile.paste(tl, (0, 0))
            tile.paste(tr, (tile_size, 0))
            tile.paste(bl, (0, tile_size))
            tile.paste(br, (tile_size, tile_size))

            # downsample 2x2 tile array to single tile
            #  okay to use Image.ANTIALIAS or Image.BICUBIC,
            #  without tRNS value
            #  ... which could leads to unwanted transparency if present
            tile = tile.resize((tile_size, tile_size), Image.ANTIALIAS)
            
            # preserve tile info except tRNS value
            if 'transparency' in tile_info.keys():
                tile_info.pop('transparency', None)
            tile.info = tile_info

            # save new tile
            save_path = os.path.join(sect_path, str(scale), \
                                     str(row) + "_" + str(col) + ext)
            tile.save(save_path, **tile_info)
            print "      row " + str(row) + " col " + str(col) + \
                  " (top " + str(top) + " left " + str(left) + ") saved"

    # set up next iteration
    scale += 1
    height = math.ceil(height / 2)
    width = math.ceil(width / 2)


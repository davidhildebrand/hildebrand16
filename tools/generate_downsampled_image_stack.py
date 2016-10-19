#!/usr/bin/env python
"""
This script requires four flags to be called when run:
    -s /path/to/source/directory
    -d /path/to/output/directory
    -i input resolution ex. 56.4 56.4 60
    -o output resolution ex. 1200 1200 1200

    The script also has an optional flag:
    -m toggle on mean averaging of images (default is median)

    This script expects a directory of images that are zero padded to 5 digits
"""

import argparse
import glob
import imghdr
import numpy
import os
from PIL import Image
import scipy.misc
import scipy.ndimage
import sys
import cv2


def create_dict_of_image_metadata(source_path):
    """This function takes a list of files (best done with glob) as input and
       will iterate through all the files, creating a dictonary that stores
       the path, image name, extension, and fullpath for each file."""
    print 'Retrieving image info.  This may take some time.'
    # Create empty dictionary
    imageinfo = {}
    global shape
    shape = None
    # Iterate through all image files in directory
    for fil in os.listdir(source_path):
        # Create pathname
        path = source_path + fil
        # Ensure file is an image
        if imghdr.what(path):
            # Extract root and extension from file name
            root, ext = os.path.splitext(fil)
            # Extract numberical slice value from file
            #slc = int(float(root.split('_')[0]) / (300./56.4))
            slc = int(float(root.split('_')[0]) / (300./56.4))
            # Add path, image, extension, and fullpath to the dictionary
            imageinfo[int(slc)] = {'path': source_path,
                                   'image': fil,
                                   'extension': ext,
                                   'fullpath': path}
            if shape is None:
                # Read the shape of the image
                shape = scipy.misc.imread(path).shape
    return imageinfo

def build_array_from_images(images, minim, maxim, max_slice_buff):
    """This function takes a given buffer minimum and maximum and creates an
       array of the images that fall within that range."""
    # Calculate buffer size from min and max
    slice_buff = maxim - minim
    print "Slice Buffer: {}".format(slice_buff)
    # Check if the buffer is greater than the maximum buffer value
    if (slice_buff > max_slice_buff):
        raise Exception('distance between {} and {} too great'
                        'for buffer of size {} images.'.format(minim,
                                                               maxim,
                                                               max_slice_buff))
    # Create empty list to collect all available images within a buffer range
    available_slices = []
    # Iterate through all images within the buffer range
    for slc in range(minim, (maxim)):
        # Check if this image exists
        if slc in images.keys():
            # If so, append it to the available_slices list
            available_slices.append(slc)
    print ("Total available slices "
           "within buffer: {}".format(len(available_slices)))
    # Create an empty array that has the shape of the images and the depth
    # of the number of available slices found
    stack = numpy.empty([shape[0], shape[1], len(available_slices)])
    # Iterate through all the images in the available slices list
    for i, slc in enumerate(available_slices):
        # Double check that the image exists
        if slc in images.keys():
            # Output to screen the path for each image
            print images[slc]['fullpath']
            # Add this image to the stack
            stack[:, :, (i)] = scipy.misc.imread(images[slc]['fullpath'])
    # Create a variable that holds the depth of the image stack
    depth = stack[2]
    print "array built!"
    print "stack shape: {}".format(stack.shape)
    return stack, depth

def average_and_resize_images(stack, secperpix, method):
    """This function takes a stack of images (generated bfrom the
       build_array_from_images function), and will resize the array to the
       desired size, as dictated by the inres and outres parameters. If -m
       is passed through as a command line argument, this function will
       average the stack using the mean. The defualt is the median."""
    # Create an array with the scaled shape
    scaledshape = (scale * numpy.array(stack.shape))
    # Average the stack with the mean method if called, OR
    if method == 'mean':
        mms = numpy.mean(stack, axis=2)
    # Average the stack with the median method. DEFAULT
    if method == 'median':
        mms = numpy.median(stack, axis=2)
    ss = tuple(map(int, numpy.ceil(scaledshape[:2])))
    # return the resized image stack
    return cv2.resize(mms, (ss[1], ss[0]))[:,:,numpy.newaxis]

def save_images(img, start, zscale):
    """This function takes the resized images an saves them to a dir"""
    # Check if directory exists, if not, create it
    if not os.path.exists(dest_path):
        os.makedirs(dest_path)
    # Iterate through each image in the stack, save it.
    for zcoord in range(len(img[0, 0])):
        fn = dest_path + '{}_downsampled.png'.format(
                 str(int(start + (float(zcoord) * 1. / zscale))).zfill(5) ) 
        if numpy.std(img[:, :, zcoord]) < 10:
            print "WARNING: Low standard deviation on {}".format(fn)
        scipy.misc.imsave(fn, img[:, :, zcoord])

def run_downsampling_protocol(images, scale, secsize, slice_buff, method):
    # Find the minimum and maximum slice indices from the image dictionary
    minim = min(images.keys())
    maxim = max(images.keys())
    # Iterate through all images in dictionary
    for i in range(minim, maxim, slice_buff):
        # Create current buffer
        j = i + slice_buff
        print "Processing slices {} to {}...".format(i, j)
        # Create stack, and obtain depth of stack
        stack, depth = build_array_from_images(images, i, j, slice_buff)
        # Downsample the image stack with selected method (Mean, median)
        downsampled_stack = average_and_resize_images(stack, depth, method)
        # Save the images to the provided directory
        save_images(downsampled_stack, i, scale[2])
        print "Downsampled images saved."
        print "Moving onto next set..."

def directory(path):
    if not os.path.isdir(path):
        err_msg = "path is not a directory (%s)"
        raise argparse.ArgumentTypeError(err_msg)
    return path


# parse command line options
parser = argparse.ArgumentParser()
parser.add_argument(
   '-s', '--source', type=directory, required=True,
   help="Path to a source directory containing image files to be downsampled."
        "[required]")
parser.add_argument(
   '-d', '--dest', type=directory, required=True,
   help="Path to a destination directory. [required]")
parser.add_argument(
   '-o', '--outres', nargs=3, required=True, action='append',
   help="Desired resolution of output image. [required]")
parser.add_argument(
   '-i', '--inres', nargs=3, required=True, action='append',
   help="Desired resolution of input image. [required]")
parser.add_argument(
   '-m', '--mean', required=False, action="store_true",
   help="Use this toggle to use Mean downsampling. [not required]")
opts = parser.parse_args()
# Set up source, output, and both resolutions from system arguments
source_path = opts.source
dest_path = opts.dest
out_res = opts.outres
in_res = opts.inres
if opts.mean:
    method = 'mean'
if not opts.mean:
    method = 'median'

if __name__ == "__main__":
    outres = numpy.array([float(out_res[0][0]), float(out_res[0][1]),
                         float(out_res[0][2])])
    inres = numpy.array([float(in_res[0][0]), float(in_res[0][1]),
                        float(in_res[0][2])])
    print "Output resolution: {}".format(outres)
    print "Input resolution: {}".format(inres)
    print "Averaging method is {}".format(method)
    img_files = glob.glob(source_path)
    imgs = create_dict_of_image_metadata(img_files[0])
    scale = inres / outres
    secsize = int(numpy.ceil(1 / scale[2]))
    default_slice_buff = int(numpy.ceil(1/scale[2]))
    run_downsampling_protocol(imgs, scale, secsize, default_slice_buff, method)

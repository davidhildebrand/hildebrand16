#!/usr/bin/env python
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


inres = numpy.array([56.4, 56.4, 60.])
outres = numpy.array([300., 300., 300.])
# outres = numpy.array([1200., 1200., 1200.])
scale = inres / outres

secsize = int(numpy.ceil(1 / scale[2]))
default_slice_buff = int(numpy.ceil(1/scale[2]))


def getimageinfo(source_path):
    print 'Retrieving image info.  This may take some time.'
    imageinfo = {}
    global shape
    shape = None
    for fil in os.listdir(source_path):
        path = source_path + fil
        if imghdr.what(path):
            root, ext = os.path.splitext(fil)
            slc = slicefromroot(root)
            imageinfo[slc] = {'path': source_path,
                              'image': fil,
                              'extension': ext,
                              'fullpath': path}
            if shape is None:
                shape = scipy.misc.imread(path).shape
    return imageinfo


def slicefromroot(root):
    # return int(root)
    return int(root.rstrip('T'))


def buildarray(images, minim, maxim, max_slice_buff):
    slice_buff = maxim - minim
    print "Slice Buffer: {}".format(slice_buff)
    if (slice_buff > max_slice_buff):
        raise Exception('distance between {} and {} too great'
                        'for buffer of size {} images.'.format(minim,
                                                               maxim,
                                                               max_slice_buff))
    available_slices = []
    for slc in range(minim, (maxim)):
        if slc in images.keys():
            available_slices.append(slc)
    print ("Total available slices "
           "within buffer: {}".format(len(available_slices)))
    stack = numpy.empty([shape[0], shape[1], len(available_slices)])
    for i, slc in enumerate(available_slices):
        if slc in images.keys():
            print images[slc]['fullpath']
            stack[:, :, (i)] = scipy.misc.imread(images[slc]
                                                 ['fullpath'])
    depth = stack[2]
    print "array built!"
    print "stack shape: {}".format(stack.shape)
    return stack, depth


def save_images(img, start, zscale=scale[2]):
    if not os.path.exists(dest_path):
        os.makedirs(dest_path)
    for zcoord in range(len(img[0, 0])):
        fn = dest_path + '{}T_down.PGM'.format(str(int(start + (float(zcoord) *
                                                1. / zscale))).zfill(5))
        if numpy.std(img[:, :, zcoord]) < 10:
            print "WARNING!  Low standard deviation on {}".format(fn)
        scipy.misc.imsave(fn, img[:, :, zcoord])


def intrazpix(stack, secperpix, mode):
    scaledshape = (scale * numpy.array(stack.shape))
    if mode == 'mean':
        mms = numpy.mean(stack, axis=2)
    if mode == 'median':
        mms = numpy.mean(stack, axis=2)
    ss = tuple(map(int, numpy.ceil(scaledshape[:2])))
    return cv2.resize(mms, ss)[:,:,numpy.newaxis]


def process_stack(images, scale, secsize, slice_buff, mode,
                  onlyintra_zpix=True):
    minim = min(images.keys())
    maxim = max(images.keys())

    for i in range(minim, maxim, slice_buff):
        j = i + slice_buff
        print "processing slices {} to {}.".format(i, j)
        stack, depth = buildarray(images, i, j, slice_buff)
        if onlyintra_zpix:
            stack = intrazpix(stack, depth, mode)
        save_images(stack, i, zscale=scale[2])
        print "Downsampled Imgaes Saved!"
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
   '-o', '--outres', nargs=3, required=False, action='append',
   help="Desired resolution of output image. [not required]")
parser.add_argument(
   '-i', '--inres', nargs=3, required=False, action='append',
   help="Desired resolution of input image. [not required]")
parser.add_argument(
   '-m', '--mean', required=False, action="store_true",
   help="Use this toggle to use Mean downsampling")
opts = parser.parse_args()

source_path = opts.source
dest_path = opts.dest
out_res = opts.outres
in_res = opts.inres
if opts.mean:
    mode = 'mean'
if not opts.mean:
    mode = 'median'

if __name__ == "__main__":
    if out_res:
        outres = [float(out_res[0][0]), float(out_res[0][1]),
                  float(out_res[0][2])]
    if in_res:
        inres = [float(in_res[0][0]), float(in_res[0][1]),
                 float(in_res[0][2])]
    print "Output Resolution: {}".format(outres)
    print "Input Resolution: {}".format(inres)
    print "Mode is set to {}".format(mode)
    img_files = glob.glob(source_path)
    imgs = getimageinfo(img_files[0])
    scale = inres / outres
    secsize = int(numpy.ceil(1 / scale[2]))
    default_slice_buff = int(numpy.ceil(1/scale[2]))
    process_stack(imgs, scale, secsize, default_slice_buff, mode)

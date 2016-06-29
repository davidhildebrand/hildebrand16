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


inres = numpy.array([56.3, 56.3, 60.])
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


def buildarray(images, minim, maxim, max_slice_buff=default_slice_buff):
    slice_buff = maxim - minim
    if (slice_buff > max_slice_buff):
        raise Exception('distance between {} and {} too great'
                        'for buffer of size {} images.'.format(minim,
                                                               maxim,
                                                               max_slice_buff))
    stack = numpy.empty([shape[0], shape[1], slice_buff]) * numpy.nan
    for slc in range(minim, (maxim)):
        if slc in images.keys():
            print images[slc]['fullpath']
            stack[:, :, (slc - minim)] = scipy.misc.imread(images[slc]
                                                           ['fullpath'])
    print "array built!"
    print "stack shape: {}".format(stack.shape)
    return stack


def save_images(img, start, zscale=scale[2]):
    if not os.path.exists(dest_path):
        os.makedirs(dest_path)
    for zcoord in range(len(img[0, 0])):
        fn = dest_path + '{}T_down.PGM'.format(str(int(start + (float(zcoord) *
                                                1. / zscale))).zfill(5))
        if numpy.std(img[:, :, zcoord]) < 10:
            print numpy.std(img[:, :, zcoord])
            print "WARNING!  Low standard deviation on {}".format(fn)
        scipy.misc.imsave(fn, img[:, :, zcoord])


def intrazpix(stack, secperpix):
    scaledshape = (scale * numpy.array(stack.shape))
    '''
    stack = numpy.ma.masked_invalid(stack)
    stack.data = numpy.nan_to_num(stack.data)
    newstack = scipy.ndimage.zoom(stack.data,
                                  [scale[0], scale[1], 1.],
                                  order=1)
    newmask = numpy.resize(stack.mask, newstack.shape)
    del stack
    stack = numpy.ma.array(newstack, mask=newmask)
    del newstack
    del newmask
    '''
    for j, subarr in enumerate(numpy.array_split(stack,
                                                 stack.shape[2] / secperpix,
                                                 axis=2)):
        print "averaging array of shape {}".format(subarr.shape)
        substack = numpy.ma.masked_invalid(subarr)
        print "Substack Created!"
        subzs = numpy.nan_to_num(substack.data)
        print 'Subzs Created!'
        newstack = scipy.ndimage.zoom(subzs,
                                      [scale[0], scale[1], 1.],
                                      order=1)
        print "New Stack Created!"
        del subzs
        marray = numpy.ma.array(newstack,
                                mask=numpy.resize(substack.mask,
                                                  newstack.shape))
        meanarr = numpy.ma.mean(marray, axis=2)
        print "Mean Array Created!"
        try:
            newarr[:, :, j] = meanarr
        except:
            newshape = list((marray.shape))
            newshape[2] = int(numpy.ceil(scaledshape[2]))
            newarr = numpy.empty(newshape)
            newarr[:, :, j] = meanarr
        print "Averaged with stdev {}".format(numpy.std(meanarr))
        print "New Array with shape {}".format(newarr.shape)
    return newarr


def process_stack(images, slice_buff=default_slice_buff, onlyintra_zpix=True):
    minim = min(images.keys())
    maxim = max(images.keys())

    for i in range(minim, maxim, slice_buff):
        j = i + slice_buff
        print "processing slices {} to {}.".format(i, j)
        stack = buildarray(images, i, j)
        if onlyintra_zpix:
            stack = intrazpix(stack, secsize)
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
   '-o', '--outres', type=int, required=False,
   help="Desired resolution of output image. [not required]")
parser.add_argument(
   '-i', '--inres', type=int, required=False,
   help="Desired resolution of input image. [not required]")
opts = parser.parse_args()

source_path = opts.source
dest_path = opts.dest
out_res = opts.outres
in_res = opts.inres


if __name__ == "__main__":
    if out_res:
        outres = numpy.array([out_res, out_res, out_res])
    if in_res:
        inres = numpy.array([in_res, in_res, in_res])
    print "Output Resolution: {}".format(outres)
    img_files = glob.glob(source_path)
    imgs = getimageinfo(img_files[0])
    process_stack(imgs)

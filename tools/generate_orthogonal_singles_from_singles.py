#!/usr/bin/env python
'''
simple script using numpy slicing to build orthogonal views of image stacks
'''
import numpy
import scipy.misc
import scipy.ndimage
import imghdr
import sys
import os

#import os
#import shutil
#
#for i in os.listdir(dimdir):
#    shutil.move(i, '{}_160515.tif'.format(int(float(i.split('_')[0]))/(300./56.4)))


dimdir = '/zflode/130201zf142/Serving/160515_SWiFT_60nmpx_singles_300iso_from56p4_med_cleanfinal/'
doutdir = '/zflode/130201zf142/Serving/160515_SWiFT_60nmpx_singles_300iso_from56p4_med_cleanfinal_reslices/xz/'

# new normal axis -- 0 is xz, 1 is yz
newnormal = 0
buffershapenormal = 28


def getimageinfo(imdir):
    print 'Retrieving image info.  This may take some time.'

    #def slicefromroot(root):
    #    return int(root.split('_')[0])
    def slicefromroot(root):
        return int(float(root.split('_')[0]) / (300./56.4))
        # return int(root.rstrip('T'))

    imageinfo = {}
    global shape  # FIXME globals are scary
    shape = None
    for fil in os.listdir(imdir):
        path = imdir + fil
        if imghdr.what(path):
            root, ext = os.path.splitext(fil)
            slc = slicefromroot(root)
            imageinfo[slc] = {'path': imdir,
                              'image': fil,
                              'extension': ext,
                              'fullpath': path,
                              'shape': scipy.misc.imread(path).shape}
            if shape is None:
                shape = scipy.misc.imread(path).shape
    return imageinfo

def buildarray(images, minax, maxax, norm=newnormal):
    imnormax = {0: lambda img: img[minax:maxax, :],
                1: lambda img: img[:, minax:maxax]}
    stackdims = buffershape
    stackdims[norm] = maxax - minax
    stack = 250 * numpy.ones(stackdims)
    flatim = stack[:, :, 0]
    for i, slc in enumerate(range(min(images.keys()), max(images.keys()))):
        if slc in images.keys():
            flatim = imnormax[norm](scipy.misc.imread(images[slc]['fullpath']))
        stack[:, :, i] = flatim
    return stack

def save_reslices(imstack, start=0, norm=None, outdir='./'):
    # TODO check if these are correct orientations for catmaid coordinates
    resliceori = {0: lambda xystack, i: xystack[i, :, :],
                  1: lambda xystack, i: xystack[:, i, :]}
    if norm in resliceori.keys():
        for i in range(imstack.shape[norm]):
            img = resliceori[norm](imstack, i)
            outfn = os.path.join(doutdir, '{}.png'.format(str(int(start + i)).zfill(5)))
            if not os.path.isdir(os.path.dirname(outfn)):
                os.makedirs(os.path.dirname(outfn))
            scipy.misc.imsave(outfn, img)
    else:
        raise Exception('Invalid new reslice normal {}.  '
                        'Please choose 0, 1'.format(norm))

def process_stack(images, buffshape=None, norm=None, outdir='./', multiprocess=False):
    if buffshape is None:
        buffshape = shape
    minim = min(images.keys())
    maxim = max(images.keys())

    for i in range(0, shape[norm], buffshape[norm]):
        print "building range {}-{}".format(i, min(i + buffshape[norm], shape[norm]))
        imstack = buildarray(images, i, min(i+buffshape[norm], shape[norm]))
        print "built!"
        save_reslices(imstack, start=i, norm=norm, outdir=outdir)


if __name__ == "__main__":
    imdir = dimdir
    outdir = doutdir

    imgs = getimageinfo(imdir)
    buffershape = [shape[0], shape[1], (max(imgs.keys()) - min(imgs.keys()))]
    buffershape[newnormal] = buffershapenormal
    process_stack(imgs, norm=newnormal, buffshape=buffershape, outdir=outdir)

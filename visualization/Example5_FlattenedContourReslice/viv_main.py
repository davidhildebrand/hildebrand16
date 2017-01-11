def main():
    import numpy, math
    volume = VIVALDI_WRITE('volume', load_data_3d('data/Zebrafish.dat'))

    alignPoints = VIVALDI_WRITE('alignPoints', None)
    colorPoints = VIVALDI_WRITE('colorPoints', None)
    alignPoints, colorPoints = VIVALDI_WRITE('alignPoints, colorPoints', readSkeleton())
    size = VIVALDI_WRITE('size', len(alignPoints))



    result = VIVALDI_WRITE('result', run_function(return_name='result', func_name='cutSkeleton', args=[volume, 'x', 'y', alignPoints, size, colorPoints, ], arg_names=['volume', 'x', 'y', 'alignPoints', 'size', 'colorPoints'], execid=[], work_range={'y':(0,500),'x':(0,size)}, merge_func='composite', merge_order='front-to-back', dtype_dict={'result':'uchar_volume'}))
    result = VIVALDI_WRITE('result', VIVALDI_GATHER(result))
    save_image(result,'result.png')


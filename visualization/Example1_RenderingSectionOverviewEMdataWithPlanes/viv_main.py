def main():
    volume = VIVALDI_WRITE('volume', load_data_3d('data/david_1G.dat'))
    mask = VIVALDI_WRITE('mask', load_data_3d('data/mask.dat'))
    hori = VIVALDI_WRITE('hori', load_data_2d('data/hori_filled.png'))

    colors = VIVALDI_WRITE('colors', numpy.fromstring(open('data/1.tf', 'rb').read(), dtype=numpy.uint8).reshape(256, 4))





    LoadIdentity()
    Rotate(-145, 0, 1, 0)
    Rotate(-18, 1, 0, 0)
    Rotate(135, 0, 0, 1)
    Translate(-980/2.0, -1740/2.0, -844/2.0)


    result = VIVALDI_WRITE('result', run_function(return_name='result', func_name='render', args=[mask, hori, volume, 'x', 'y', colors, ], arg_names=['mask', 'hori', 'volume', 'x', 'y', 'colors'], execid=[], work_range={'y':(-1536,1536),'x':(-1536,1536)}, dtype_dict={'volume':'uchar_volume','mask':'uchar_volume','hori':'uchar_volume'}))

    result = VIVALDI_WRITE('result', VIVALDI_GATHER(result))
    save_image(result, 'result.tif')

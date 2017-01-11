def main():
    volume = VIVALDI_WRITE('volume', load_data_3d('data/new_990_860_1600.dat'))

    xy_mask = VIVALDI_WRITE('xy_mask', load_data_2d('data/XY_guide.png'))
    xz136_mask = VIVALDI_WRITE('xz136_mask', load_data_2d('data/top.png'))
    xz332_mask = VIVALDI_WRITE('xz332_mask', load_data_2d('data/bottom.png'))
    green_mask = VIVALDI_WRITE('green_mask', load_data_2d('data/601_mask1.png'))

    import numpy
    data2 = VIVALDI_WRITE('data2', open('data/new_points1_out.txt').readlines())
    data1 = VIVALDI_WRITE('data1', [[int(float(val)) for val in elem.split(' ')] for elem in data2])
    data1 = VIVALDI_WRITE('data1', data1 + [[0, data1[-1][1], data1[-1][2]]])
    data = VIVALDI_WRITE('data', numpy.array(data1, dtype=numpy.float32))

    LoadIdentity()
    Rotate(150, 0, 1, 0)
    Rotate(-8, 1, 0, 0)
    Translate(-990/2, -860/2, -1600/2)


    tr_data = VIVALDI_WRITE('tr_data', open("data/1.tf","rb").read())
    transf = VIVALDI_WRITE('transf', numpy.fromstring(tr_data, dtype=numpy.uint8).reshape(256,4).astype(numpy.float32))

    data_size = VIVALDI_WRITE('data_size', len(data))

    result = VIVALDI_WRITE('result', run_function(return_name='result', func_name='render', args=[volume, 'x', 'y', xy_mask, xz136_mask, xz332_mask, data, data_size, transf, green_mask, ], arg_names=['volume', 'x', 'y', 'xy_mask', 'xz136_mask', 'xz332_mask', 'data', 'data_size', 'transf', 'green_mask'], execid=[], work_range={'y':(-1536,1536),'x':(-1536,1536)}))

    result = VIVALDI_WRITE('result', VIVALDI_GATHER(result))
    save_image(result)

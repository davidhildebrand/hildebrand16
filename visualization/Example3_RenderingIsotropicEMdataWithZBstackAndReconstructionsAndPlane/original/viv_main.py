def main():
    import sys
    if len(sys.argv) < 3:
        print_bold("USAGE : vivaldi %s COLOR [START=90] [END=91]"%sys.argv[1])
        print_bold("THIS SOURCE : Rendering New aligned data with skeletons and stack for SpinalBackfills and a plane")
        print_bold("Default Source is used [Blue & Purple Skeletons & Topview]")


    volume = VIVALDI_WRITE('volume', load_data_3d("data/161210_All.dat"))
    colorMapName = VIVALDI_WRITE('colorMapName', "data/161201_all_color")

    volume_mask  = VIVALDI_WRITE('volume_mask', load_data_3d("data/Mask.dat"))

    tr_data = VIVALDI_WRITE('tr_data', open("data/Zebrafish.tf","rb").read())
    transf = VIVALDI_WRITE('transf', numpy.fromstring(tr_data, dtype=numpy.uint8).reshape(256,4).astype(numpy.float32))


    colors = VIVALDI_WRITE('colors', open(colorMapName).readlines())
    colors = VIVALDI_WRITE('colors', [[int(elem.split(' ')[0]), int(elem.split(' ')[1]), int(elem.split(' ')[2])] for elem in colors])
    colors = VIVALDI_WRITE('colors', numpy.array(colors, dtype=numpy.uint8))
    colors = VIVALDI_WRITE('colors', colors[:5])
    len_col = VIVALDI_WRITE('len_col', len(colors))

    Rot_Start = VIVALDI_WRITE('Rot_Start', 90 if len(sys.argv) < 4 else int(sys.argv[3]))
    Rot_End   = VIVALDI_WRITE('Rot_End', 91 if len(sys.argv) < 5 else int(sys.argv[4]))


    Plane_normal = VIVALDI_WRITE('Plane_normal', numpy.array([[-0.999854064595206, 0.016858099721706, -0.002766583853512]], dtype=numpy.float32))
    Plane_loc    = VIVALDI_WRITE('Plane_loc', numpy.array([[ 268761.67377/600 , 135687.563833/600, 404620.713987/600 ]], dtype=numpy.float32))


    input_col = VIVALDI_WRITE('input_col', numpy.array([[171, 15, 15, 255]],dtype=numpy.float32))

    start = VIVALDI_WRITE('start', 3700)
    end = VIVALDI_WRITE('end', 65535)



    volume_brain = VIVALDI_WRITE('volume_brain', load_data_3d('data/SpinalBackfills.dat'))
    tr_data1 = VIVALDI_WRITE('tr_data1', open("data/SpinalBackfills_23zf.tf","rb").read())

    transf1  = VIVALDI_WRITE('transf1', numpy.fromstring(tr_data1, dtype=numpy.uint8).reshape(256,4).astype(numpy.float32))


    from datetime import datetime
    today = VIVALDI_WRITE('today', "%02d%02d%02d_"%(int(str(datetime.now().year)[2:]),(datetime.now().month),(datetime.now().day)))


    import os
    if not os.path.exists('result/%srotate'%today):
        os.system('mkdir -p result/%srotate'%today)

    for elem in range(Rot_Start,Rot_End):
        val = VIVALDI_WRITE('val', 1.0)

        LoadIdentity()
        Rotate(elem+0.01, 1, 0, 0)
        Rotate(90, 0, 1, 0)
        Translate(-963/2, -866/2, -1621/2)

        result = VIVALDI_WRITE('result', run_function(return_name='result', func_name='render', args=[volume, volume_mask, volume_brain, 'x', 'y', colors, transf, transf1, start, end, input_col, len_col, Plane_normal, Plane_loc, ], arg_names=['volume', 'volume_mask', 'volume_brain', 'x', 'y', 'colors', 'transf', 'transf1', 'start', 'end', 'input_col', 'len_col', 'Plane_normal', 'Plane_loc'], execid=[], work_range={'y':(-800,1000),'x':(-1024,1024)}))

        result = VIVALDI_WRITE('result', VIVALDI_GATHER(result))
        save_image(result, "%srotate/%s%03d.png"%(today,today,elem))


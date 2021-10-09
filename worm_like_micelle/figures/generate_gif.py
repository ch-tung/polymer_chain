from IPython import display

import glob
import imageio

anim_file = './ring/ring.gif'

with imageio.get_writer(anim_file, mode='I') as writer:
    filenames = glob.glob('./ring/ring_test_*.png')
    filenames = sorted(filenames)
    for filename in filenames:
        image = imageio.imread(filename)
        writer.append_data(image)
    image = imageio.imread(filename)
    writer.append_data(image)

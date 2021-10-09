from IPython import display

import glob
import imageio

anim_file = './chain/chain.gif'

with imageio.get_writer(anim_file, mode='I') as writer:
    filenames = glob.glob('./chain/chain_test_*.png')
    filenames = sorted(filenames)
    for filename in filenames:
        image = imageio.imread(filename)
        writer.append_data(image)
    image = imageio.imread(filename)
    writer.append_data(image)

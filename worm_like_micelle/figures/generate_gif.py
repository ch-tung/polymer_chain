from IPython import display

import glob
import imageio

anim_file = './WLM_a2_DP100.gif'

with imageio.get_writer(anim_file, mode='I') as writer:
    filenames = glob.glob('./WLM_a2_DP100_*.png')
    filenames = sorted(filenames)
    for filename in filenames:
        image = imageio.imread(filename)
        writer.append_data(image)
    image = imageio.imread(filename)
    writer.append_data(image)

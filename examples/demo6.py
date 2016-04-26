import os
import pylab
import sys
import time

BUILD_PATH = os.environ.get('BUILD_PATH') + '/xGadgetron'
SRC_PATH = os.environ.get('SRC_PATH') + '/xGadgetron/pGadgetron'

sys.path.append(BUILD_PATH)
sys.path.append(SRC_PATH)

import pGadgetron
import pGadgets

def main():
    # acquisitions will be read from this HDF file
    file = str(input('raw data file (with apostrophys in Python2.*): '))
    input_data = pGadgetron.MR_Acquisitions(file)
##    input_data = pGadgetron.MR_Acquisitions('opismrmrd.h5')

    # define gadgets
    gadget1 = pGadgets.RemoveROOversamplingGadget()
    gadget2 = pGadgets.SimpleReconstructionGadget()
    gadget3 = pGadgets.ExtractGadget()

    # set gadgets parameters
    gadget2.set_property('trigger_dimension', 'repetition')
    gadget2.set_property('split_slices', 'true')

    # create reconstruction object
    recon = pGadgetron.ImagesReconstructor()

    # build gadgets chain
    recon.add_gadget('g1', gadget1)
    recon.add_gadget('g2', gadget2)
    recon.add_gadget('g3', gadget3)

    # connect to input data
    recon.set_input(input_data)
    # perform reconstruction
    recon.process()
    
    # get reconstructed images
    images = recon.get_output()

    # plot reconstructed images

    nz = images.number()
    print('%d images' % nz)

    print('Please enter z-coordinate of the slice to view it')
    print('(a value outside the range [0 : %d] will stop this loop)'%(nz - 1))
    while True:
        s = str(input('z-coordinate: '))
        if len(s) < 1:
            break
        z = int(s)
        if z < 0 or z >= nz:
            break
        data = images.image_as_array(z)
        pylab.figure(z)
        pylab.imshow(data[0,0,:,:])
        print('delete the plot window to continue...')
        pylab.show()

    # write images to a new group in 'output6.h5'
    # named after the current date and time
    time_str = time.asctime()
    images.write('output6.h5', time_str)

try:
    main()
except pGadgetron.error as err:
    # display error information
    print ('Gadgetron exception occured:\n', err.value)


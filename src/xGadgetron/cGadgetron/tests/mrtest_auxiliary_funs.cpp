#include "mrtest_auxiliary_funs.h"


#include <vector>
#include <cmath>
#include <fstream>
#include <iostream>

#include <ismrmrd/ismrmrd.h>

void sirf::preprocess_acquisition_data(MRAcquisitionData& ad)
{
    std::cout << "Processing Acquisition Data" << std::endl;

    sirf::AcquisitionsProcessor preprocessing_chain;

    auto sptr_noise_gadget = std::make_shared<Gadget>(NoiseAdjustGadget());
    auto sptr_ro_overs_gadget = std::make_shared<Gadget>(RemoveROOversamplingGadget());
    auto sptr_asymmecho_gadget = std::make_shared<Gadget>(AsymmetricEchoAdjustROGadget());

    preprocessing_chain.add_gadget("dummy1", sptr_noise_gadget);
    preprocessing_chain.add_gadget("dummy2", sptr_asymmecho_gadget);
    preprocessing_chain.add_gadget("dummy3", sptr_ro_overs_gadget);

    preprocessing_chain.process(ad);
    auto sptr_preproc_ad =preprocessing_chain.get_output();

    ISMRMRD::Acquisition acq;
    for(int i=0; i<sptr_preproc_ad->number(); ++i)
    {
        sptr_preproc_ad->get_acquisition(i, acq);
        ad.set_acquisition(i, acq);
    }
    ad.set_acquisitions_info( sptr_preproc_ad->acquisitions_info());

}

void sirf::write_cfimage_to_raw(std::string const fname_prefix, CFImage& img)
{


    std::stringstream fname_out;
    fname_out << fname_prefix;
    fname_out << "_" << img.getMatrixSizeX();
    fname_out << "x" << img.getMatrixSizeY();
    fname_out << "x" << img.getMatrixSizeZ() * img.getNumberOfChannels();
    fname_out << ".raw";

    std::cout << "Writing " << fname_out.str() << std::endl;

    std::vector<float> dat;

    for(size_t i=0; i<img.getNumberOfDataElements(); ++i)
        dat.push_back( std::abs( *(img.getDataPtr() + i )));


    std::ofstream myfile( fname_out.str(), std::ios::out | std::ios::binary);
    myfile.write((char*)(&dat[0]), sizeof(float)*dat.size());

}


void sirf::write_cfimage_to_raw(std::string const fname_prefix, ImageWrap& iw)
{
    void* vptr_img = iw.ptr_image();
    CFImage* ptr_img = static_cast<CFImage*>(vptr_img);

    sirf::write_cfimage_to_raw(fname_prefix, *ptr_img);
}

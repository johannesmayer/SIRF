#include "mrtest_auxiliary_funs.h"


void sirf::preprocess_acquisition_data(MRAcquisitionData& ad, MRAcquisitionData& preproc_ad){

    sirf::AcquisitionsProcessor preprocessing_chain;

    auto sptr_noise_gadget = std::make_shared<Gadget>(NoiseAdjustGadget());
    auto sptr_ro_overs_gadget = std::make_shared<Gadget>(RemoveROOversamplingGadget());
    auto sptr_asymmecho_gadget = std::make_shared<Gadget>(AsymmetricEchoAdjustROGadget());
    auto sptr_acquisition_finish_gadget = std::make_shared<Gadget>(AcquisitionFinishGadget());

    preprocessing_chain.add_gadget("", sptr_noise_gadget);
    preprocessing_chain.add_gadget("", sptr_asymmecho_gadget);
    preprocessing_chain.add_gadget("", sptr_ro_overs_gadget);
    preprocessing_chain.add_gadget("", sptr_acquisition_finish_gadget);

    preprocessing_chain.process(mr_rawdata);
    preproc_ad = *(preprocessing_chain.get_output());

}





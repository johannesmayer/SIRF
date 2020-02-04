#include "sirf/Gadgetron/encoding.h"



#include <sstream>

#include "sirf/iUtilities/LocalisedException.h"


using namespace sirf;
using namespace ISMRMRD;

void sirf::aTrajectoryPreparation::update_acquisitions_info(MRAcquisitionData& mr_acq)
{

    IsmrmrdHeader hdr = mr_acq.acquisitions_info().get_IsmrmrdHeader();

    if(hdr.encoding.size() > 0)
        throw LocalisedException("Currrently only files with one encoding are supported", __FILE__, __LINE__);

    hdr.encoding[0].trajectory = this->traj_type_;

    this->kspace_encoding_ = hdr.encoding[0];

    std::stringstream info_stream;
    serialize(hdr, info_stream);
    mr_acq.set_acquisitions_info(info_stream.str());

}

void sirf::CartesianTrajectoryPrep::set_trajectory(MRAcquisitionData& mr_acq)
{
    update_acquisitions_info(mr_acq); // do nothing for cartesian trajectories
}


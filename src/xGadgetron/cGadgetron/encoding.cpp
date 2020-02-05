#include "sirf/Gadgetron/encoding.h"



#include <sstream>
#include <math.h>

#include "sirf/iUtilities/LocalisedException.h"


using namespace sirf;
using namespace ISMRMRD;

#define SIRF_GOLDEN_ANGLE M_PI*(3-sqrt(5))


void sirf::aTrajectoryPreparation::update_acquisitions_info(MRAcquisitionData& mr_acq)
{

    IsmrmrdHeader hdr = mr_acq.acquisitions_info().get_IsmrmrdHeader();


    if(hdr.encoding.size() != 1)
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


void sirf::GRPETrajectoryPrep::set_trajectory(MRAcquisitionData& mr_acq)
{
    update_acquisitions_info(mr_acq);

    for(size_t ia=0; ia<mr_acq.number(); ++ia)
    {
        Acquisition acq;
        mr_acq.get_acquisition(ia, acq);
        this->set_acquisition_trajectory(acq);
        mr_acq.set_acquisition(ia, acq);
    }
}

void sirf::GRPETrajectoryPrep::set_acquisition_trajectory(Acquisition& acq)
{
    acq.resize(acq.number_of_samples(),acq.active_channels(), this->traj_dim_);
    std::vector<float> acq_traj = this->calculate_trajectory(acq);
    acq.setTraj(&acq_traj[0]);
}

std::vector<float> sirf::GRPETrajectoryPrep::calculate_trajectory(Acquisition& acq)
{
    const ISMRMRD::Limit rad_lims = this->kspace_encoding_.encodingLimits.kspace_encoding_step_1.get();
    const ISMRMRD::Limit ang_lims = this->kspace_encoding_.encodingLimits.kspace_encoding_step_2.get();

    const ISMRMRD::EncodingCounters idx = acq.idx();

    float const pe_angle = SIRF_GOLDEN_ANGLE * idx.kspace_encode_step_2;

    size_t const num_diff_shifts = this->rad_shift_.size();
    float rad_shift = float( this->rad_shift_.at(this->circ_mod(idx.kspace_encode_step_2 - ang_lims.center,num_diff_shifts))) / float(num_diff_shifts);

    float pe_radius = idx.kspace_encode_step_1 - rad_lims.center;
    pe_radius = (pe_radius==0) ? pe_radius : pe_radius+rad_shift;

    float const traj_norm = 2*std::max<float>(( rad_lims.center - rad_lims.minimum + 0), (rad_lims.maximum - rad_lims.center + (num_diff_shifts-1)/num_diff_shifts));
    pe_radius /= traj_norm;

    std::vector<float> traj;

    for(size_t i_sample=0; i_sample<acq.number_of_samples();++i_sample)
    {
        traj.push_back(0); //dummy for RPE as the readout is cartesian
        traj.push_back(pe_radius * cos( pe_angle ));
        traj.push_back(pe_radius * sin( pe_angle ));
    }

    return traj;
}

#ifndef ENCODING_H
#define ENCODING_H

#include <vector>

#include <ismrmrd/xml.h>
#include "sirf/Gadgetron/gadgetron_data_containers.h"

class aTrajectoryPreparation
{
public:
    aTrajectoryPreparation(){};
    virtual void set_trajectory(sirf::MRAcquisitionData& mr_acq)=0;
protected:

    ISMRMRD::EncodingLimits enc_lims_;
    ISMRMRD::TrajectoryType traj_type_;
    int traj_dim_;
    virtual std::vector<float> compute_trajectory(ISMRMRD::Acquisition& acq);

};























#endif // ENCODING_H

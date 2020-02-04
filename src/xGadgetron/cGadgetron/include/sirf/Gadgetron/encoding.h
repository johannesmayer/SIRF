#ifndef ENCODING_H
#define ENCODING_H

#include <vector>

#include <ismrmrd/xml.h>
#include "sirf/Gadgetron/gadgetron_data_containers.h"

namespace sirf{

class aTrajectoryPreparation{

public:
    aTrajectoryPreparation(){}
    virtual void set_trajectory(sirf::MRAcquisitionData& mr_acq)=0;
protected:

    void update_acquisitions_info(sirf::MRAcquisitionData& mr_acq);

    ISMRMRD::Encoding kspace_encoding_;
    ISMRMRD::TrajectoryType traj_type_;
    int traj_dim_;
    virtual std::vector<float> compute_trajectory(ISMRMRD::Acquisition& acq)=0;

};



class CartesianTrajectoryPrep : public aTrajectoryPreparation{

public:
    CartesianTrajectoryPrep(): aTrajectoryPreparation() {
        traj_type_ = ISMRMRD::TrajectoryType::CARTESIAN;
        traj_dim_ = 0;
    }
    virtual void set_trajectory(sirf::MRAcquisitionData& mr_acq);
protected:
    virtual std::vector<float> compute_trajectory(ISMRMRD::Acquisition& acq){return std::vector<float>{};}

};




} // end namespace
#endif // ENCODING_H

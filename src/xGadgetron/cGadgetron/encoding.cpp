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




void sirf::Cartesian3DFourierEncoding::forward(CFImage* ptr_img, MRAcquisitionData& ac)
{
    CFImage& img = *ptr_img;

    unsigned int nx = img.getMatrixSizeX();
    unsigned int ny = img.getMatrixSizeY();
    unsigned int nz = img.getMatrixSizeZ();
    unsigned int nc = img.getNumberOfChannels();

    unsigned int readout = nx; //assumes the acquisitions have oversampling removed

    std::vector<size_t> dims;
    dims.push_back(readout);
    dims.push_back(ny);
    dims.push_back(nz);
    dims.push_back(nc);

    ISMRMRD::NDArray<complex_float_t> ci(dims);
    memset(ci.getDataPtr(), 0, ci.getDataSize());

    for (unsigned int c = 0; c < nc; c++) {
        for (unsigned int z = 0; z < nz; z++) {
            for (unsigned int y = 0; y < ny; y++) {
                for (unsigned int x = 0; x < nx; x++) {
                    ci(x, y, z, c) = (complex_float_t)img(x, y, z, c);
                }
            }
        }
    }

    fft3c(ci);

    memset((void*)acq.getDataPtr(), 0, acq.getDataSize());

    ISMRMRD::Acquisition acq;

    for(size_t i =0; i<ac.items(); ++i)
    {
        ac.get_acquisition(i, acq);
        acq.resize(nx, nc, 3);
        ky = acq.idx().kspace_encode_step_1;
        kz = acq.idx().kspace_encode_step_2;
        for (unsigned int c = 0; c < nc; c++) {
            for (unsigned int s = 0; s < nx; s++) {
                acq.data(xout, c) = ci(s, ky, kz, c);
            }
        }
    }

    ac.append_acquisition(acq);
}

void sirf::Cartesian3DFourierEncoding::backward(CFImage* ptr_img, MRAcquisitionData& ac)
{

    std::string par;
    ISMRMRD::IsmrmrdHeader header;
    par = ac.acquisitions_info();
    ISMRMRD::deserialize(par.c_str(), header);
    ISMRMRD::Encoding e = header.encoding[0];
    ISMRMRD::Acquisition acq;
    for (unsigned int i = 0; i < ac.number(); i++) {
        ac.get_acquisition(i, acq);
        if (!TO_BE_IGNORED(acq))
            break;
    }

    unsigned int nx = e.reconSpace.matrixSize.x;
    unsigned int ny = e.reconSpace.matrixSize.y;
    unsigned int nz = e.reconSpace.matrixSize.z;
//    unsigned int ny = e.encodedSpace.matrixSize.y;
//    unsigned int nz = e.encodedSpace.matrixSize.z;
    unsigned int nc = acq.active_channels();
    unsigned int readout = acq.number_of_samples();

    std::vector<size_t> dims;
    dims.push_back(readout);
    dims.push_back(ny);
    dims.push_back(nz);
    dims.push_back(nc);

    ISMRMRD::NDArray<complex_float_t> ci(dims);
    memset(ci.getDataPtr(), 0, ci.getDataSize());
    const int NUMVAL = 4;
    typedef std::array<int, NUMVAL> tuple;
    tuple t_first;
    bool first = true;
    unsigned int& a = off;
    for (; a < ac.number(); a++) {
        ac.get_acquisition(a, acq);
        if (TO_BE_IGNORED(acq))
            continue;
        tuple t;
        t[0] = acq.idx().repetition;
        t[1] = acq.idx().phase;
        t[2] = acq.idx().contrast;
        t[3] = acq.idx().slice;
        if (first) {
            t_first = t;
            first = false;
            std::cout << "new slice: ";
            for (int i = 0; i < NUMVAL; i++)
                std::cout << t[i] << ' ';
        }
        else if (t != t_first &&
            !acq.isFlagSet(ISMRMRD::ISMRMRD_ACQ_LAST_IN_MEASUREMENT))
            break;
        int yy = acq.idx().kspace_encode_step_1;
        int zz = acq.idx().kspace_encode_step_2;
        for (unsigned int c = 0; c < nc; c++) {
            for (unsigned int s = 0; s < readout; s++) {
                ci(s, yy, zz, c) = acq.data(s, c);
            }
        }
    }
    std::cout << "done\n";
    /*
    int y = 0;
    for (;;){
        ac.get_acquisition(off + y, acq);
        if (acq.isFlagSet(ISMRMRD::ISMRMRD_ACQ_FIRST_IN_SLICE))
            break;
        y++;
    }
    for (;;) {
        ac.get_acquisition(off + y, acq);
        int yy = acq.idx().kspace_encode_step_1;
        int zz = acq.idx().kspace_encode_step_2;
        for (unsigned int c = 0; c < nc; c++) {
            for (unsigned int s = 0; s < readout; s++) {
                ci(s, yy, zz, c) = acq.data(s, c);
            }
        }
        y++;
        if (acq.isFlagSet(ISMRMRD::ISMRMRD_ACQ_LAST_IN_SLICE))
            break;
    }
    off += y;
    */

    ifft3c(ci);

    T* ptr = im.getDataPtr();
    T s;
    memset(ptr, 0, im.getDataSize());
    long long int i = 0;
    for (unsigned int c = 0; c < nc; c++) {
        i = 0;
        for (unsigned int z = 0; z < nz; z++) {
            for (unsigned int y = 0; y < ny; y++) {
                for (unsigned int x = 0; x < nx; x++, i++) {
                    uint16_t xout = x + (readout - nx) / 2;
                    complex_float_t zi = ci(xout, y, z, c);
                    complex_float_t zc = csm(x, y, z, c);
                    xGadgetronUtilities::convert_complex(std::conj(zc) * zi, s);
                    ptr[i] += s;
                }
            }
        }
    }
}




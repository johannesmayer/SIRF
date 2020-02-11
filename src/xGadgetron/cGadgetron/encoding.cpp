#include "sirf/Gadgetron/encoding.h"



#include <sstream>
#include <math.h>

#include <cassert>
#define assertm(exp, msg) assert(((void)msg, exp)

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


void sirf::FourierEncoding::match_img_header_to_acquisition(CFImage& img, const ISMRMRD::Acquisition& acq)
{

    auto acq_hdr = acq.getHead();
    auto idx = acq_hdr.idx;

    img.setAverage(idx.average);
    img.setSlice(idx.slice);
    img.setContrast(idx.contrast);
    img.setPhase(idx.phase);
    img.setRepetition(idx.repetition);
    img.setSet(idx.set);

}



void sirf::Cartesian3DFourierEncoding::forward(CFImage* ptr_img, MRAcquisitionData& ac)
{
//    CFImage& img = *ptr_img;

//    unsigned int nx = img.getMatrixSizeX();
//    unsigned int ny = img.getMatrixSizeY();
//    unsigned int nz = img.getMatrixSizeZ();
//    unsigned int nc = img.getNumberOfChannels();

//    unsigned int readout = nx; //assumes the acquisitions have oversampling removed

//    std::vector<size_t> dims;
//    dims.push_back(readout);
//    dims.push_back(ny);
//    dims.push_back(nz);
//    dims.push_back(nc);

//    ISMRMRD::NDArray<complex_float_t> ci(dims);
//    memset(ci.getDataPtr(), 0, ci.getDataSize());

//    for (unsigned int c = 0; c < nc; c++) {
//        for (unsigned int z = 0; z < nz; z++) {
//            for (unsigned int y = 0; y < ny; y++) {
//                for (unsigned int x = 0; x < nx; x++) {
//                    ci(x, y, z, c) = (complex_float_t)img(x, y, z, c);
//                }
//            }
//        }
//    }

//    fft3c(ci);

//    memset((void*)acq.getDataPtr(), 0, acq.getDataSize());

//    ISMRMRD::Acquisition acq;

//    for(size_t i =0; i<ac.items(); ++i)
//    {
//        ac.get_acquisition(i, acq);
//        acq.resize(nx, nc, 3);
//        ky = acq.idx().kspace_encode_step_1;
//        kz = acq.idx().kspace_encode_step_2;
//        for (unsigned int c = 0; c < nc; c++) {
//            for (unsigned int s = 0; s < nx; s++) {
//                acq.data(xout, c) = ci(s, ky, kz, c);
//            }
//        }
//    }

//    ac.append_acquisition(acq);
}

void sirf::Cartesian3DFourierEncoding::backward(CFImage* ptr_img, MRAcquisitionData& ac)
{

    std::string par;
    ISMRMRD::IsmrmrdHeader header;
    par = ac.acquisitions_info();
    ISMRMRD::deserialize(par.c_str(), header);

    if( header.encoding.size() > 1)
        LocalisedException("Currently only one encoding is supported per rawdata file.", __FUNCTION__, __LINE__);

    ISMRMRD::Encoding e = header.encoding[0];

    ISMRMRD::Acquisition acq;

    ac.get_acquisition(0, acq);

    unsigned int nx = e.encodedSpace.matrixSize.x;
    unsigned int ny = e.encodedSpace.matrixSize.y;
    unsigned int nz = e.encodedSpace.matrixSize.z;
    unsigned int nc = acq.active_channels();

    unsigned int readout = acq.number_of_samples();

    if( readout > nx)
        LocalisedException("The readout is larger than the encoded space. Possibly need to remove oversampling in r.o. direction first.", __FUNCTION__, __LINE__);

    std::vector<size_t> dims, dims_dcf;

    dims.push_back(nx);
    dims.push_back(ny); dims_dcf.push_back(ny);
    dims.push_back(nz); dims_dcf.push_back(nz);
    dims.push_back(nc);

    ISMRMRD::Limit ky_lim = e.encodingLimits.kspace_encoding_step_1.get();
    ISMRMRD::Limit kz_lim = e.encodingLimits.kspace_encoding_step_2.get();

    ISMRMRD::NDArray<complex_float_t> ci(dims);
    ISMRMRD::NDArray<complex_float_t> dcf(dims);

    memset(ci.getDataPtr(), 0, ci.getDataSize());
    memset(dcf.getDataPtr(), 1, dcf.getDataSize());

    for (int a=0; a < ac.number(); a++) {
        ac.get_acquisition(a, acq);
        int yy = ny/2 - ky_lim.center + acq.idx().kspace_encode_step_1 ;
        int zz = nz/2 - kz_lim.center + acq.idx().kspace_encode_step_2;
        int ro_offset = nx/2 - acq.center_sample();
        for (unsigned int c = 0; c < nc; c++) {
            for (unsigned int s = 0; s < readout; s++) {
                ci(ro_offset + s, yy, zz, c) += acq.data(s, c);
                dcf(yy, zz) += 1;
            }
        }
    }

    // correct for double acquired PE points
    for(unsigned int c=0; c<nc; ++c)
    for(unsigned int z=0; z<nz; ++z)
    for(unsigned int y=0; y<ny; ++y)
    for(unsigned int x=0; x<nx; ++x)
        ci(x,y,z,c) /= dcf(x,y);


    // now if image and kspace have different dimension then you need to interpolate or pad with zeros here

    ifft3c(ci);
    CFImage img = *ptr_img;

    img.resize(nx, ny, nz, nc);
    memcpy(img.begin(), ci.begin(), ci.getDataSize());


    // set the header correctly of the image

}




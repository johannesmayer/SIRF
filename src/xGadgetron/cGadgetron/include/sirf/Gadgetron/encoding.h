#ifndef ENCODING_H
#define ENCODING_H

#include <vector>

#include <ismrmrd/xml.h>
#include "sirf/Gadgetron/gadgetron_data_containers.h"

#include "sirf/iUtilities/LocalisedException.h"

#include <gadgetron/hoNDArray.h>
#include <gadgetron/vector_td.h>
#include <gadgetron/vector_td_utilities.h>

#include <gadgetron/hoNDFFT.h>
#include <gadgetron/hoNFFT.h>

/*!
\file
\ingroup Fourier Encoding
\brief Specification file for preparing MRAcquisitionData for Fourier encoding.

\author Johannes Mayer
\author CCP PETMR
*/

/*!
\ingroup Fourier Encoding
\brief Abstract class for trajectory preparation

*/
namespace sirf{

class aTrajectoryPreparation{

public:
    aTrajectoryPreparation(){}
    virtual void set_trajectory(sirf::MRAcquisitionData& mr_acq)=0;



protected:

    void update_acquisitions_info(sirf::MRAcquisitionData& mr_acq);

    ISMRMRD::Encoding kspace_encoding_;
    ISMRMRD::TrajectoryType traj_type_;

    uint16_t traj_dim_;
    virtual void set_acquisition_trajectory(ISMRMRD::Acquisition& acq)=0;
    virtual std::vector<float> calculate_trajectory(ISMRMRD::Acquisition& acq)=0;
};

/*!
\ingroup Fourier Encoding
\brief Cartesian trajectory preparation class

*/

class CartesianTrajectoryPrep : public aTrajectoryPreparation{

public:
    CartesianTrajectoryPrep(): aTrajectoryPreparation() {
        traj_type_ = ISMRMRD::TrajectoryType::CARTESIAN;
        traj_dim_ = 0;
    }

    virtual void set_trajectory(sirf::MRAcquisitionData& mr_acq);

protected:
    virtual void set_acquisition_trajectory(ISMRMRD::Acquisition& acq){}
    virtual std::vector<float> calculate_trajectory(ISMRMRD::Acquisition& acq){return std::vector<float>{};}
};

/*!
\ingroup Fourier Encoding
\brief Golden Radial Phase Encoding interleaved trajectory preparation class.

*/



class GRPETrajectoryPrep : public aTrajectoryPreparation {

public:
    GRPETrajectoryPrep(): aTrajectoryPreparation() {
        traj_type_ = ISMRMRD::TrajectoryType::OTHER;
        traj_dim_ = 3;
    }

    virtual void set_trajectory(sirf::MRAcquisitionData& mr_acq);

protected:
    virtual void set_acquisition_trajectory(ISMRMRD::Acquisition& acq);
    virtual std::vector<float> calculate_trajectory(ISMRMRD::Acquisition& acq);
    std::vector< uint16_t > const rad_shift_ = {0, 2, 1, 3}; //this is bit-reversed {0 1 2 3}
    uint16_t circ_mod(uint16_t const a, uint16_t const b){ return (((a%b) + b ) % b);}
};




class FourierEncoding
{
public:
    FourierEncoding(){}

    virtual void forward(MRAcquisitionData& ac, const CFImage* ptr_img)=0;
    virtual void backward(CFImage* ptr_img, const MRAcquisitionData& ac)=0;

    void match_img_header_to_acquisition(CFImage& img, const ISMRMRD::Acquisition& acq);
};

class CartesianFourierEncoding : public FourierEncoding
{
public:
    CartesianFourierEncoding() : FourierEncoding() {}

    virtual void forward(MRAcquisitionData& ac, const CFImage* ptr_img);
    virtual void backward(CFImage* ptr_img, const MRAcquisitionData& ac);

};

typedef Gadgetron::hoNDArray<Gadgetron::floatd2> SirfTrajectoryType2D;

class RPEFourierEncoding : public FourierEncoding
{
public:
    RPEFourierEncoding(): FourierEncoding() {}

    virtual void forward(MRAcquisitionData& ac, const CFImage* ptr_img);
    virtual void backward(CFImage* ptr_img, const MRAcquisitionData& ac);

    SirfTrajectoryType2D get_trajectory(const MRAcquisitionData& ac) const;

};


using namespace Gadgetron;

class Gridder_2D
{

public:

    Gridder_2D(std::vector<size_t> img_dims_output, const SirfTrajectoryType2D &traj) : nufft_operator_(from_std_vector<size_t, 2>(img_dims_output), (float)this->oversampling_factor_, (float)this->kernel_size_)
    {
        setup_nufft(img_dims_output, traj);
    }

    void setup_nufft(std::vector<size_t> img_dims_output, const SirfTrajectoryType2D &traj)
    {
        if( img_dims_output.size() != 2)
            throw LocalisedException("The image dimensions of the output should be of size 2." , __FILE__, __LINE__);

        traj.get_dimensions(this->trajdims_);

        this->output_dims_ = img_dims_output;

//        this->nufft_operator_ = hoNFFT_plan<float,2>(from_std_vector<size_t, 2>(img_dims_output), (float)this->oversampling_factor_, (float)this->kernel_size_);
        this->nufft_operator_.preprocess(traj);
    }

    void fft()
    {
        throw LocalisedException("Forward gridding not implemented yet." , __FILE__, __LINE__);
    }

    void ifft(const hoNDArray<std::complex<float> >& kdata)
    {

        auto sptr_unit_dcw = std::make_shared<Gadgetron::hoNDArray<float> >( this->trajdims_);
        sptr_unit_dcw ->fill(1.f);

        Gadgetron::hoNDArray<std::complex<float> > result(this->output_dims_);
        result.fill(std::complex<float>(0.f, 0.f));

        this->nufft_operator_.compute(kdata, result, sptr_unit_dcw.get(), Gadgetron::NFFT_comp_mode::BACKWARDS_NC2C);

    }

protected:
    static const size_t oversampling_factor_ = 2;
    static size_t const kernel_size_ = 2;

    std::vector<size_t> trajdims_;
    std::vector<size_t> output_dims_;

    Gadgetron::hoNFFT_plan<float, 2> nufft_operator_;
};

} // namespace sirf
#endif // ENCODING_H
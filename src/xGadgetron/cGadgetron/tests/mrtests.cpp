#include <iostream>
#include <cstdlib>


#include "sirf/Gadgetron/encoding.h"
#include "sirf/Gadgetron/gadgetron_data_containers.h"


using namespace sirf;

bool test_TrajectoryPreparation_constructors( void )
{
    try
    {
        std::cout << "Running test " << __FUNCTION__ << std::endl;

        sirf::CartesianTrajectoryPrep cart_tp;
        sirf::GRPETrajectoryPrep rpe_tp;

        return true;

    }
    catch( std::runtime_error const &e)
    {
        std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
        std::cout << e.what() << std::endl;
        throw e;
    }
}


bool test_GRPETrajectoryPrep_set_trajectory( void )
{
    try
    {
        std::cout << "Running test " << __FUNCTION__ << std::endl;
        sirf::GRPETrajectoryPrep rpe_tp;

        std::string const fname_input = "/media/sf_CCPPETMR/TestData/Input/xGadgetron/cGadgetron/ResPhantom_CV_nav_GRPE_itl_ismrmrd.h5";
        std::string const fname_output = "/media/sf_CCPPETMR/TestData/Output/xGadgetron/cGadgetron/res_phant_cv_rpe_itl_gc_ismrmrd.h5";

        sirf::AcquisitionsVector mr_dat;
        mr_dat.read(fname_input);
        rpe_tp.set_trajectory(mr_dat);
        mr_dat.write(fname_output);

        return true;

    }
    catch( std::runtime_error const &e)
    {
        std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
        std::cout << e.what() << std::endl;
        throw e;
    }
}

bool test_apply_combine_coil_sensitivities( void )
{
    try
    {
        std::cout << "Running test " << __FUNCTION__ << std::endl;

        size_t const Nx=64;
        size_t const Ny=64;
        size_t const Nz=64;
        size_t const Nc=4;

        CFImage unit_img  = CFImage(Nx,Ny,Nz,1);

        for(size_t i=0;i<unit_img.getNumberOfDataElements();++i)
        {
            *(unit_img.begin()+i) = complex_float_t(1.f, 0.f);
        }
        sirf::GadgetronImagesVector giv, miv;
        giv.append(unit_img);

        // generate artifical coil maps
        size_t const num_dat_pts = Nx*Ny*Nz*Nc;
        std::vector<float> csm_real(num_dat_pts);
        std::vector<float> csm_imag(num_dat_pts);

        for(size_t nc=0; nc<Nc; nc++)
        for(size_t nz=0; nz<Nz; nz++)
        for(size_t ny=0; ny<Ny; ny++)
        for(size_t nx=0; nx<Nx; nx++)
        {
            size_t access_idx = (((nc*Nz + nz)*Ny + ny)*Nx + nx);
            csm_real.at(access_idx) = nx;
            csm_imag.at(access_idx) = nc;
        }

        CoilSensitivitiesAsImages csm;
        csm.append_csm(Nx, Ny, Nz, Nc, &csm_real[0], &csm_imag[0]);

        CFImage csm_img = csm.get_csm_as_CFImage();

        csm.apply_coil_sensitivities(miv, giv);

        bool store_dcm=true;
        if(store_dcm)
        {
            std::stringstream fname_output;
            fname_output << "/media/sf_CCPPETMR/TestData/Output/xGadgetron/cGadgetron/output_" << __FUNCTION__ <<"_application" ;

            fname_output<<".h5";
            miv.write(fname_output.str());
        }

        csm.combine_coils(giv, miv);

        if(store_dcm)
        {
            std::stringstream fname_output;
            fname_output << "/media/sf_CCPPETMR/TestData/Output/xGadgetron/cGadgetron/output_" << __FUNCTION__ <<"_combination" ;

            fname_output<<".h5";
            giv.write(fname_output.str());
        }

        return true;

    }
    catch( std::runtime_error const &e)
    {
        std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
        std::cout << e.what() << std::endl;
        throw e;
    }
}

bool test_get_kspace_order(void)
{
    try
    {
        std::cout << "Running test " << __FUNCTION__ << std::endl;

        std::string const fpath_input = "/media/sf_CCPPETMR/TestData/Input/xGadgetron/cGadgetron/";
        std::string fname_input = fpath_input + "CV_SR_64Cube_1Echo_10Dyn.h5";

        sirf::AcquisitionsVector av;
        av.read(fname_input);

        auto kspace_sorting = av.get_kspace_order();

        fname_input = fpath_input + "CV_SR_128Cube_1Echo_3Dyn.h5";

        sirf::AcquisitionsVector av_contrast;
        av_contrast.read(fname_input);

        auto kspace_sorting_contrast = av_contrast.get_kspace_order();

        return true;

    }
    catch( std::runtime_error const &e)
    {
        std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
        std::cout << e.what() << std::endl;
        throw e;
    }
}

bool test_get_subset()
{
    try
    {
        std::cout << "Running test " << __FUNCTION__ << std::endl;

        std::string const fpath_input = "/media/sf_CCPPETMR/TestData/Input/xGadgetron/cGadgetron/";
        std::string fname_input = fpath_input + "CV_SR_64Cube_1Echo_10Dyn.h5";

        sirf::AcquisitionsVector av;
        av.read(fname_input);

        std::vector<int> subset_idx;
        for(int i=0; i<av.number()/10; ++i)
            subset_idx.push_back(i);

        sirf::AcquisitionsVector subset;
        av.get_subset(subset, subset_idx);

        return true;

    }
    catch( std::runtime_error const &e)
    {
        std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
        std::cout << e.what() << std::endl;
        throw e;
    }
}

int main ()
{
	try{

//        test_TrajectoryPreparation_constructors();
//        test_GRPETrajectoryPrep_set_trajectory();
//        test_apply_combine_coil_sensitivities();
//        test_get_kspace_order();
//        test_get_subset();
        return 0;
	}
    catch(const std::exception &error) {
        std::cerr << "\nHere's the error:\n\t" << error.what() << "\n\n";
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}


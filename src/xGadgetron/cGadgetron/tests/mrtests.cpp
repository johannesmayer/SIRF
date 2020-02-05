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

        std::string const fname_input = "/media/sf_CCPPETMR/TestData/Input/xGadgetron/cGadgetron/20180716-113005,ResPhant,CV_nav_rpe_itl_golden,2090,26_ismrmrd.h5";
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



int main ()
{
	try{

        test_TrajectoryPreparation_constructors();
        test_GRPETrajectoryPrep_set_trajectory();

        return 0;
	}
	catch(...)
	{

	}
}

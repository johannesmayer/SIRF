#include <iostream>
#include <cstdlib>


#include "sirf/Gadgetron/encoding.h"

using namespace sirf;

bool test_TrajectoryPreparation_constructors( void )
{
    try
    {
        std::cout << "Running test " << __FUNCTION__ << std::endl;
        sirf::CartesianTrajectoryPrep();
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
        return 0;
	}
	catch(...)
	{

	}
}

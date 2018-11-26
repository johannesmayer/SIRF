/* ================================================

Author: Johannes Mayer
Date: 2018.04.03
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */



#pragma once 


#include <string>
#include <ismrmrd/ismrmrd.h>
#include <vector>

#include "auxiliary_testing_functions.h"

namespace test_aux_test_funs
{
	
bool test_get_serialized_ismrmrd_header( void );
bool test_get_mock_acquisition_vector( void );
bool test_get_mock_csm( void );
bool test_get_mock_gaussian_csm( void );
bool test_get_mock_coildata_as_cfimage( void );
bool test_get_mock_ismrmrd_image_with_cube( void );
bool test_get_mock_contrast_generator( void );
bool test_get_mock_pet_contrast_generator( void );
bool test_get_mock_sawtooth_signal( void );

}
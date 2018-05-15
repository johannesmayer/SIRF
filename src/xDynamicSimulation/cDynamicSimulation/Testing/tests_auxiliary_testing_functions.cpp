/* ================================================

Author: Johannes Mayer
Date: 2018.04.03
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */
#include "tests_auxiliary_testing_functions.h"



#include "sstream"

#include <ismrmrd/ismrmrd.h>
#include "gadgetron_data_containers.h" 

#include "auxiliary_input_output.h"






bool test_aux_test_funs::test_get_serialized_ismrmrd_header( void )
{
	try
	{
		std::string serialized_hdr = aux_test::get_serialized_mock_ismrmrd_header();
		std::cout << serialized_hdr << std::endl;
	}
	catch( std::runtime_error const &e)
	{
		std::cout << "Exception caught in " <<__FUNCTION__ <<" .!" <<std::endl;
		std::cout << e.what() << std::endl;
		throw e;
	}
	return true;
}



bool test_aux_test_funs::test_get_mock_acquisition_vector( void )
{

	ISMRMRD::IsmrmrdHeader hdr = aux_test::get_mock_ismrmrd_header();
	AcquisitionsVector acq_vec = aux_test::get_mock_acquisition_vector( hdr );

	unsigned int num_acquis = acq_vec.number();
	std::cout<< epiph(num_acquis) << std::endl;


	ISMRMRD::Acquisition acq;
	int const check_aqu_num[3] = {0, 10, 100};

	for( int i=0; i<3; i++)
	{	
		std::cout << epiph( check_aqu_num[i] ) << std::endl;
		acq_vec.get_acquisition(check_aqu_num[i], acq);

		uint16_t const available_channels = acq.available_channels();
		
	}

	return true;
}
bool test_aux_test_funs::test_get_mock_csm( void )

{

try
	{
		ISMRMRD::NDArray<complex_float_t> csm = aux_test::get_mock_csm();
		
		size_t const * dim = csm.getDims();
		

		ISMRMRD::Image<float> csm_abs(dim[0], dim[1], dim[2], dim[3]);
		for( int i=0; i<csm.getNumberOfElements(); i++)			
			*(csm_abs.begin() + i) = std::abs( *(csm.begin() + i) );

		std::stringstream name_stream;
		name_stream << "/media/sf_SharedFiles/test_mock_csm_" << dim[0] << "x" <<  dim[1] << "x" << dim[2]*dim[3];
		data_io::write_raw<float>(name_stream.str(), csm_abs.begin(), csm_abs.getNumberOfDataElements());

	}
	catch( std::runtime_error const &e)
	{
		std::cout << "Exception caught in " <<__FUNCTION__ <<" .!" <<std::endl;
		std::cout << e.what() << std::endl;
		throw e;
	}
	return true;

}

bool test_aux_test_funs::test_get_mock_coildata_as_cfimage( void )
{

try
	{
		CoilDataAsCFImage csm = aux_test::get_mock_coildata_as_cfimage();
		ISMRMRD::Image<complex_float_t> csm_image = csm.image();

		int dim[4];
		csm.get_dim(dim);

		ISMRMRD::Image<float> csm_abs((size_t)dim[0], (size_t)dim[1], (size_t)dim[2], (size_t)dim[3]);
		for( int i=0; i<csm_image.getNumberOfDataElements(); i++)			
			*(csm_abs.begin() + i) = std::abs( *(csm_image.begin() + i) );


		std::stringstream name_stream;
		name_stream << "/media/sf_SharedFiles/test_mock_coildata_asscfimage_" << dim[0] << "x" <<  dim[1] << "x" << dim[2]*dim[3];

		data_io::write_raw<float>(name_stream.str(), csm_abs.begin(), csm_abs.getNumberOfDataElements());

	}
	catch( std::runtime_error const &e)
	{
		std::cout << "Exception caught in " <<__FUNCTION__ <<" .!" <<std::endl;
		std::cout << e.what() << std::endl;
		throw e;
	}
	return true;

}


bool test_aux_test_funs::test_get_mock_ismrmrd_image_with_cube( void )
{

try
	{
		ISMRMRD::Image< complex_float_t > img = aux_test::get_mock_ismrmrd_image_with_cube();

		std::vector< size_t > dim;
		dim.push_back(img.getMatrixSizeX ());
		dim.push_back(img.getMatrixSizeY ());
		dim.push_back(img.getMatrixSizeZ ());

		ISMRMRD::NDArray<float> dat( dim );
		for( int i=0; i<img.getNumberOfDataElements(); i++)
		{
			*(dat.begin() + i)  = std::abs( *(img.begin() + i ) );
		}


		std::stringstream name_stream;
		name_stream << "/media/sf_SharedFiles/test_get_mock_ismrmrd_image_with_cube_" << dim[0] << "x" <<  dim[1] << "x" << dim[2];

		data_io::write_raw<float>(name_stream.str(), dat.begin(), dat.getNumberOfElements());

	}
	catch( std::runtime_error const &e)
	{
		std::cout << "Exception caught in " <<__FUNCTION__ <<" .!" <<std::endl;
		std::cout << e.what() << std::endl;
		throw e;
	}
	return true;

}

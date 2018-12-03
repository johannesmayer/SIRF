/* ================================================

Author: Johannes Mayer
Date: 2018.11.13
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */

#pragma once 


// for easier logging 
#define epiph(x) #x << " = " << x


#define SHARED_FOLDER_PATH "/media/sf_SharedFolder/CCPPETMR/"
#define ANALYZE_OUTPUT_TESTPATH SHARED_FOLDER_PATH "analyze_test_output"

#define USE_128_CUBE_INPUT

#ifdef USE_64_CUBE_INPUT

	// #define ISMRMRD_H5_TEST_PATH  SHARED_FOLDER_PATH "h5_source_files/CV_SR_64Cube_1Echo_10Dyn.h5"
	#define ISMRMRD_H5_TEST_PATH  SHARED_FOLDER_PATH "h5_source_files/CV_nav_cart_64Cube_1Echo.h5"
	// #define ISMRMRD_H5_TEST_PATH  SHARED_FOLDER_PATH "h5_source_files/CV_nav_cart_64Cube_3Echo.h5" 
	#define H5_XCAT_PHANTOM_PATH  SHARED_FOLDER_PATH "h5_phantom_input/xcat_phantom_64_cubed.h5"
	// #define H5_XCAT_PHANTOM_PATH  SHARED_FOLDER_PATH "h5_phantom_input/xcat_phantom_32x32x64.h5"
	#define DISPLACEMENT_FIELD_PATH SHARED_FOLDER_PATH "temp_folder_motion_dyn_0/motion_field_0.hdr"

	#define PET_TEMPLATE_CONTRAST_IMAGE_DATA_PATH SHARED_FOLDER_PATH ""
	#define PET_TEMPLATE_ACQUISITION_IMAGE_DATA_PATH SHARED_FOLDER_PATH ""
	#define PET_TEMPLATE_ACQUISITION_DATA_PATH SHARED_FOLDER_PATH ""

#elif defined(USE_128_CUBE_INPUT)

	#define H5_XCAT_PHANTOM_PATH  SHARED_FOLDER_PATH "h5_phantom_input/xcat_phantom_128_cubed.h5"

	#define ISMRMRD_H5_TEST_PATH  SHARED_FOLDER_PATH "h5_source_files/CV_SR_128Cube_1Echo_10Dyn.h5"
	// #define ISMRMRD_H5_TEST_PATH  SHARED_FOLDER_PATH "h5_source_files/CV_nav_cart_128Cube_1Echo.h5"   
	// #define ISMRMRD_H5_TEST_PATH  SHARED_FOLDER_PATH "h5_source_files/CV_nav_128_rpe_itl_golden.h5"  
	// #define ISMRMRD_H5_TEST_PATH  SHARED_FOLDER_PATH "h5_source_files/CV_nav_128_rpe_sfl_gc_usos8.h5" 


	#define DISPLACEMENT_FIELD_PATH SHARED_FOLDER_PATH ""
	
	#define PET_TEMPLATE_CONTRAST_IMAGE_DATA_PATH SHARED_FOLDER_PATH "pet_source_files/template_image_input_contgen.hv"
	#define PET_TEMPLATE_ACQUISITION_IMAGE_DATA_PATH SHARED_FOLDER_PATH "pet_source_files/template_image_input_acquisition.hv"
	#define PET_TEMPLATE_ACQUISITION_DATA_PATH SHARED_FOLDER_PATH "pet_source_files/template_acquisition_input.hs"

#elif defined(USE_192_CUBE_INPUT)
	// #define ISMRMRD_H5_TEST_PATH  SHARED_FOLDER_PATH "h5_source_files/20180720-105656,SimulationDummy,CV_nav_192Cube_3Echo,2528,141_ismrmrd.h5" 
	// #define H5_XCAT_PHANTOM_PATH  SHARED_FOLDER_PATH "h5_phantom_input/.h5"
	#define DISPLACEMENT_FIELD_PATH SHARED_FOLDER_PATH ""

#elif defined(USE_208_CUBE_INPUT)
	// #define H5_XCAT_PHANTOM_PATH  SHARED_FOLDER_PATH "h5_phantom_input/xcat_phantom_128x128x208.h5"
	#define H5_XCAT_PHANTOM_PATH  SHARED_FOLDER_PATH "h5_phantom_input/xcat_phantom_208_cubed.h5"

	#define ISMRMRD_H5_TEST_PATH  SHARED_FOLDER_PATH "ISMRMSimInput/MR/meas_MID00241_FID69145_Tho_T1_fast_ismrmrd.h5" 
	#define DISPLACEMENT_FIELD_PATH SHARED_FOLDER_PATH ""

	#define PET_TEMPLATE_CONTRAST_IMAGE_DATA_PATH SHARED_FOLDER_PATH "ISMRMSimInput/PET/template_image_input_contgen.hv"
	#define PET_TEMPLATE_ACQUISITION_IMAGE_DATA_PATH SHARED_FOLDER_PATH "ISMRMSimInput/PET/template_image_input_acquisition.hv"
	#define PET_TEMPLATE_ACQUISITION_DATA_PATH SHARED_FOLDER_PATH "ISMRMSimInput/PET/template_span11.hs"

#endif

#define TIME_POINTS_CARDIAC_PATH SHARED_FOLDER_PATH "ISMRMSimInput/card_time"
#define CARDIAC_SIGNAL_PATH SHARED_FOLDER_PATH "ISMRMSimInput/card_signal"

#define TIME_POINTS_RESP_PATH SHARED_FOLDER_PATH "ISMRMSimInput/resp_time"
#define RESP_SIGNAL_PATH SHARED_FOLDER_PATH "ISMRMSimInput/resp_signal"


#define XML_TEST_PATH SHARED_FOLDER_PATH "XMLTestData/test_TissueParameters_XML.xml" 
#define XML_XCAT_PATH SHARED_FOLDER_PATH "XMLTestData/XCAT_TissueParameters_XML.xml" 

// #define H5_PHANTOM_TEST_PATH  SHARED_FOLDER_PATH "h5_testfile_cube_size3.h5"
#define H5_PHANTOM_TEST_PATH  SHARED_FOLDER_PATH "h5_phantom_input/xcat_phantom_128_cubed.h5"


#define ACQU_FILE_NAME  SHARED_FOLDER_PATH "acquisitions_file_fwd_test.h5"

#define FILENAME_MR_RPE_SIM SHARED_FOLDER_PATH "testoutput_mr_rpe_simulation.h5"
#define FILENAME_MR_CONTRAST_DYNSIM  SHARED_FOLDER_PATH "testoutput_mr_dynamic_contrast_simulation.h5"
#define FILENAME_MR_MOTION_DYNSIM SHARED_FOLDER_PATH "testoutput_mr_dynamic_motion_simulation.h5"
#define FILENAME_MR_MOTION_CONTRAST_DYNSIM SHARED_FOLDER_PATH "testoutput_mr_dynamic_motion_contrast_simulation.h5"


#define FILENAME_STATICSIM_PET SHARED_FOLDER_PATH "testoutput_pet_static_simulation.hs"


#define FILENAME_MR_DEFORM_TEST SHARED_FOLDER_PATH "output_deforming_tests/deformed_img"



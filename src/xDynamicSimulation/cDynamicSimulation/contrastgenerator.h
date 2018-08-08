/* ================================================

Author: Johannes Mayer
Date: 2018.03.28
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */

#pragma once

#include <string>
#include <memory>

#include <ismrmrd/ismrmrd.h>
#include <ismrmrd/xml.h>

#include "tissueparameters.h"
#include "tissuelabelmapper.h"


#include "stir_types.h"
#include "stir_data_containers.h"

// base class for contrast generators. inherit for different modalities.
// Reading the header is the same for each modality (hopefully!!!).
// otherwise the header reading can be just inherited and implemented
// as a virutal method

using ISMRMRD::NDArray;
using ISMRMRD::IsmrmrdHeader;




#define CASE_MAP_PET_CONTRAST 0
#define CASE_MAP_PET_ATTENUATION 1


typedef std::vector<float> SeqParamType;

class AbstractContrastGenerator {

public:


	AbstractContrastGenerator(LabelArray tissue_labels, std::string const filename_tissue_parameter_xml);
	

	// pure virtual since formats are very diff for pet and mri and ct
	virtual void map_contrast()=0;

	

protected:

	TissueLabelMapper tlm_;

};



class MRContrastGenerator : public AbstractContrastGenerator {

public:
	
	MRContrastGenerator (LabelArray tissue_labels, std::string const filename_tissue_parameter_xml);

	void set_rawdata_header(IsmrmrdHeader hdr);
	void map_contrast();

	std::vector< ISMRMRD::Image< complex_float_t> > get_contrast_filled_volumes();

	void match_output_dims_to_headerinfo( void );

private:

	std::vector< ISMRMRD::Image< complex_float_t> > contrast_filled_volumes_;
	IsmrmrdHeader hdr_;

};


std::vector < complex_float_t > map_flash_contrast( TissueParameter const * const ptr_to_tiss_par, 
													ISMRMRD::IsmrmrdHeader * ptr_to_header);





class PETContrastGenerator : public AbstractContrastGenerator {

public:

	PETContrastGenerator ( LabelArray tissue_labels, std::string const filename_tissue_parameter_xml );


	void set_template_image_from_file ( std::string const filename_header_with_ext ) 
	{
	 	 this->template_pet_image_data_ = sirf::PETImageData(filename_header_with_ext);
	 	 this->template_img_is_set_ = true;
	 	 
	};

	std::vector< int > get_dimensions( void );
	std::vector< float > get_voxel_sizes( void );

	std::vector< sirf::PETImageData > get_contrast_filled_volumes();

	void map_contrast();
	void map_attenuation();

private:

	bool template_img_is_set_ = false;

	std::vector< int > img_template_dims_;
	std::vector< float > img_template_spacing_; 

	std::vector < sirf::PETImageData > contrast_filled_volumes_;
	void map_tissueparams_member(int const case_map);
	
	sirf::PETImageData template_pet_image_data_;

};

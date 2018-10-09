/* ================================================

Author: Johannes Mayer
Date: 2018.03.27
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */

#pragma once


#include <stdio.h>
#include <stdlib.h>

#include <sstream>
#include <string>
#include <vector>
#include <ismrmrd/ismrmrd.h>
#include <map>

#include <utility>
#include <memory>

#include "tissueparameters.h"


typedef unsigned int LabelType;
typedef std::vector< std::shared_ptr<TissueParameter> > TissueVector;
typedef ISMRMRD::NDArray<LabelType> LabelArray;

class TissueLabelMapper{

public:
	TissueLabelMapper();
	TissueLabelMapper(LabelArray const label_array, std::string const xml_path);

	inline TissueVector get_segmentation_tissues (void)
	{
		return this->segmentation_tissues_;
	};

	std::string get_filepath_tissue_parameter_xml( void );
	const size_t* get_segmentation_dimensions( void );

	LabelArray get_segmentation_labels( void );
	TissueParameterList get_tissue_parameter_list( void );
	
	void map_labels_to_tissue_from_xml( void );

	void replace_petmr_tissue_parameters( const LabelType&  label, const TissueParameter& tiss);

	// NEEDS TO STAY PUBLIC!
	void assign_tissues_to_labels( void );
	

// private:

	std::string filepath_tissue_parameter_xml_;

	TissueParameterList tissue_parameter_list_;
	
	LabelArray segmentation_labels_;
	TissueVector segmentation_tissues_;
	
};


// public methods for the class
TissueVector assign_tissue_parameters_to_labels( TissueParameterList &tiss_list, LabelArray label_list );				
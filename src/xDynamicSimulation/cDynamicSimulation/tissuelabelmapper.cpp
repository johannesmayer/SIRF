/* ================================================

Author: Johannes Mayer
Date: 2018.03.27
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */


#include "tissuelabelmapper.h"



TissueLabelMapper::TissueLabelMapper() {}

TissueLabelMapper::TissueLabelMapper(LabelArray const segmentation_labels, std::string const filepath_tissue_parameter_xml)
{
	this->filepath_tissue_parameter_xml_ = filepath_tissue_parameter_xml;
	this->segmentation_labels_ = segmentation_labels;

}


std::string TissueLabelMapper::get_filepath_tissue_parameter_xml()
{
	return this->filepath_tissue_parameter_xml_;
}


LabelArray TissueLabelMapper::get_segmentation_labels( void )
{
	return this->segmentation_labels_;
}

void TissueLabelMapper::map_labels_to_tissue_from_xml( void )
{
	this->tissue_parameter_list_ = read_TissueParameters_from_xml(filepath_tissue_parameter_xml_);
	this->segmentation_tissues_ = assign_tissue_parameters_to_labels( tissue_parameter_list_, segmentation_labels_);
}


const size_t* TissueLabelMapper::get_segmentation_dimensions( void )
{
	return this->segmentation_labels_.getDims();		
}

TissueVector assign_tissue_parameters_to_labels( TissueParameterList &tiss_list, LabelArray label_volume )
{

	size_t num_tissue_params = tiss_list.size();

	std::map <int, TissueParameter* >  lut;

	for(int i =0; i<num_tissue_params; i++)
	{
		lut.insert(std::make_pair( tiss_list[i].label_, &tiss_list[i]));	//map label to pointer
	}

	size_t const num_voxels = label_volume.getNumberOfElements();

	TissueVector tissue_volume;
	tissue_volume.resize(num_voxels);

	for( int i_vox =0; i_vox<num_voxels; i_vox++)
	{
		auto key_value_pair = lut.find( *(label_volume.begin() +i_vox) );
		if( key_value_pair != lut.end())
		{	
			tissue_volume[i_vox] = key_value_pair->second;
		}
		else
		{	
			std::stringstream msg;
			msg << "The label " <<  *(label_volume.begin() +i_vox) << " in your label volume does not appear in the label list.";
			throw std::runtime_error(msg.str());
		}

	}

	return tissue_volume;
}


/* ================================================

Author: Johannes Mayer
Date: 2018.08.01
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */


#include <sstream>
#include <stdexcept>
#include <deque>
#include <algorithm>

#include <boost/filesystem.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/system/error_code.hpp>


#include <ismrmrd/ismrmrd.h>

#include "dynamics.h"
#include "auxiliary_input_output.h"

using namespace sirf;


bool is_in_bin( SignalAxisType const signal, SignalBin const bin)
{

	auto bin_min = std::get<0>(bin);
	auto bin_max = std::get<2>(bin);

	if( bin_min < bin_max )
		return (signal >= bin_min && signal <= bin_max);
	else if ( bin_min > bin_max )
		return (signal >= bin_min || signal <= bin_max);
	else
		return false;
}


 
AcquisitionsVector intersect_mr_acquisition_data(AcquisitionsVector one_dat, AcquisitionsVector other_dat)
{

	bool one_dat_is_smaller = ( one_dat.items() >= other_dat.items() );

	typedef std::vector<uint32_t> CounterBox;

	CounterBox one_counters, other_counters;


	for( size_t i=0; i<one_dat.items(); i++)
	{
		ISMRMRD::Acquisition acq;

		one_dat.get_acquisition(i, acq);
		one_counters.push_back(acq.getHead().scan_counter);
	}

	for( size_t i=0; i<other_dat.items(); i++)
	{
		ISMRMRD::Acquisition acq;

		other_dat.get_acquisition(i, acq);
		other_counters.push_back(acq.getHead().scan_counter);
	}
	
	std::sort(one_counters.begin(), one_counters.end() );
	std::sort(other_counters.begin(), other_counters.end() );

	CounterBox intersected_counters(one_counters.size() + other_counters.size()); 

	CounterBox::iterator it;


	it = std::set_intersection( one_counters.begin(), one_counters.end(),
							    other_counters.begin(), other_counters.end(), 
							    intersected_counters.begin() );

	intersected_counters.resize( it - intersected_counters.begin() );
	

	MRDataType intersection;
	intersection.copy_acquisitions_info(one_dat);

	MRDataType& smaller_data_container = one_dat_is_smaller ? one_dat : other_dat;

	for( size_t i=0; i<smaller_data_container.items(); i++)
	{
		ISMRMRD::Acquisition acq;

		smaller_data_container.get_acquisition(i, acq);
		uint32_t acquis_counter = acq.getHead().scan_counter;
		if(std::find(intersected_counters.begin(), intersected_counters.end(), acquis_counter) != intersected_counters.end()) 
		{
			intersection.append_acquisition(acq);
    	} 
	}

	return intersection;

}


aDynamic::aDynamic(int const num_simul_states) : num_simul_states_(num_simul_states)
{
	set_bins( num_simul_states );
}

void aDynamic::set_bins(int const num_bins)
{
	this->signal_bins_.clear();

	for(int i_state=0; i_state<num_bins; i_state++)
	{	
		SignalBin bin;

		std::get<0>(bin) = SignalAxisType(i_state)/SignalAxisType(num_bins) - 1.f/(2*num_bins);
		std::get<1>(bin) = SignalAxisType(i_state)/SignalAxisType(num_bins);
		std::get<2>(bin) = SignalAxisType(i_state)/SignalAxisType(num_bins) + 1.f/(2*num_bins);
		
		if( std::get<0>(bin) < 0 )
			std::get<0>(bin) = ( 1 + std::get<0>(bin) );

		if( std::get<1>(bin) < 0 )
			std::get<1>(bin) = ( 1 + std::get<1>(bin) );

		if( std::get<2>(bin) < 0 )
			std::get<2>(bin) = ( 1 + std::get<2>(bin) );

		this->signal_bins_.push_back( bin );
	}
}


void aDynamic::set_num_simul_states(int const num_states)
{
	this->num_simul_states_ = num_states;
	set_bins( num_states );

}

void aDynamic::set_dyn_signal(SignalContainer signal) 
{
	this->dyn_signal_ = signal;
}


SignalAxisType aDynamic::linear_interpolate_signal(TimeAxisType time_point)
{

	size_t const num_sig_points = this->dyn_signal_.size();
	
	size_t first_bigger_thant_time_point=-1;

	for( size_t i=0; i<num_sig_points; i++)
	{
		if(this->dyn_signal_[i].first > time_point)
		{
			first_bigger_thant_time_point = i;
			break;
		}
	}

	SignalAxisType interpol_signal;

	if( first_bigger_thant_time_point == 0)
		interpol_signal = this->dyn_signal_[0].second;
	else if( first_bigger_thant_time_point == -1 )
		interpol_signal = this->dyn_signal_[num_sig_points-1].second;
	else
	{
		interpol_signal = dyn_signal_[first_bigger_thant_time_point-1].second + 
							(time_point - dyn_signal_[first_bigger_thant_time_point-1].first )
						   *(dyn_signal_[first_bigger_thant_time_point].second - dyn_signal_[first_bigger_thant_time_point-1].second)
						   /(dyn_signal_[first_bigger_thant_time_point].first  - dyn_signal_[first_bigger_thant_time_point-1].first );
	}

	return interpol_signal;

}



void aDynamic::bin_mr_acquisitions( AcquisitionsVector all_acquisitions )
{

	if(this->dyn_signal_.size() == 0)
		throw std::runtime_error( "Please set a signal first. Otherwise you cannot bin your data, you dummy!" );


	std::deque< size_t > relevant_acq_numbers;
	std::deque< size_t > acq_not_binned;


	for( size_t i=0; i<all_acquisitions.items(); i++)
		relevant_acq_numbers.push_back( i );

	for( int i_bin=0; i_bin<this->signal_bins_.size(); i_bin++)
	{

		auto bin = this->signal_bins_[i_bin];
	

		AcquisitionsVector curr_acq_vector;
		curr_acq_vector.copy_acquisitions_info(all_acquisitions);
		
		ISMRMRD::Acquisition acq;
		acq_not_binned.clear();

		while( relevant_acq_numbers.size() > 0 )	
		{
			auto curr_pos = relevant_acq_numbers[0];
			relevant_acq_numbers.pop_front();	
			
			all_acquisitions.get_acquisition(curr_pos, acq);
			
			auto acq_hdr = acq.getHead();
			
			TimeAxisType acq_time = (TimeAxisType)acq_hdr.acquisition_time_stamp;
			
			SignalAxisType signal_of_acq = this->linear_interpolate_signal( acq_time );
			if( is_in_bin(signal_of_acq, bin) )
			{
				curr_acq_vector.append_acquisition(acq);
			}
			else
			{
				acq_not_binned.push_back(curr_pos);
			}
			
		}
	
		relevant_acq_numbers.swap(acq_not_binned);
		this->binned_mr_acquisitions_.push_back( curr_acq_vector );
		
	}
}


void MotionDynamic::set_bins( int const num_bins )
{
	
	aDynamic::set_bins(num_bins);

	// this->signal_bins_.clear();

	// for(int i_state=0; i_state<num_bins; i_state++)
	// {	
	// 	SignalBin bin;

	// 	std::get<0>(bin) = SignalAxisType(i_state)/SignalAxisType(num_bins) - 1.f/(2*num_bins);
	// 	std::get<1>(bin) = SignalAxisType(i_state)/SignalAxisType(num_bins);
	// 	std::get<2>(bin) = SignalAxisType(i_state)/SignalAxisType(num_bins) + 1.f/(2*num_bins);
		
	// 	if( std::get<0>(bin) < 0 )
	// 		std::get<0>(bin) = ( 1 + std::get<0>(bin) );

	// 	if( std::get<1>(bin) < 0 )
	// 		std::get<1>(bin) = ( 1 + std::get<1>(bin) );

	// 	if( std::get<2>(bin) < 0 )
	// 		std::get<2>(bin) = ( 1 + std::get<2>(bin) );

	// 	this->signal_bins_.push_back( bin );
	// }
}


ContrastDynamic::ContrastDynamic(int const num_simul_states) : aDynamic()
{ 
	this->num_simul_states_ =num_simul_states;
	this->set_bins(num_simul_states_);
}

void ContrastDynamic::set_bins( int const num_bins )
{
	this->signal_bins_.clear();

	for(int i_state=0; i_state<num_bins; i_state++)
	{	
		SignalBin bin;

		std::get<0>(bin) = SignalAxisType(i_state)/SignalAxisType(num_bins);
		std::get<1>(bin) = SignalAxisType(i_state)/SignalAxisType(num_bins) + 1.f/(2*num_bins);
		std::get<2>(bin) = SignalAxisType(i_state)/SignalAxisType(num_bins) + 1.f/(num_bins);
	
		this->signal_bins_.push_back( bin );
	}
}

void ContrastDynamic::set_parameter_extremes(TissueParameter tiss_at_0, TissueParameter tiss_at_1)
{
	this->tissue_parameter_extremes_.first = tiss_at_0;
	this->tissue_parameter_extremes_.second = tiss_at_1;
}


TissueParameterList ContrastDynamic::get_interpolated_tissue_params(SignalAxisType const signal)
{
	TissueParameterList tiss_list;

	for(size_t i=0; i< this->list_cont_var_labels_.size(); i++)
	{
		TissueParameter curr_par = ((1-signal) * this->tissue_parameter_extremes_.first + signal * this->tissue_parameter_extremes_.second);
		curr_par.name_ = "";	// name info is lost unfortunately, but better than the wrong information
		curr_par.label_ = this->list_cont_var_labels_[i];

		tiss_list.push_back( curr_par);
	}

	return tiss_list;
}


int MotionDynamic::num_total_motion_dynamics_ = 0;

MotionDynamic::MotionDynamic():aDynamic()
{
	this->which_motion_dynamic_am_i_ = num_total_motion_dynamics_;
	this->num_total_motion_dynamics_ += 1;
	
	this->temp_folder_name_ = setup_tmp_folder_name();

}

MotionDynamic::MotionDynamic(int const num_simul_states) : aDynamic(num_simul_states)
{
	this->which_motion_dynamic_am_i_ = num_total_motion_dynamics_;
	this->num_total_motion_dynamics_ += 1;

	this->temp_folder_name_ = setup_tmp_folder_name();

}


MotionDynamic::~MotionDynamic()
{ 
	if( this->destroy_upon_deletion_)
		this->delete_temp_folder();

	this->num_total_motion_dynamics_ -= 1; 
}


int MotionDynamic::get_which_motion_dynamic_am_i(){ return this->which_motion_dynamic_am_i_; }
int MotionDynamic::get_num_total_motion_dynamics(){ return this->num_total_motion_dynamics_; }

std::string MotionDynamic::setup_tmp_folder_name()
{
	std::string const current_folder_prefix = "temp_folder_motion_dyn_";
	std::stringstream tmp_stream;
	tmp_stream << this->temp_folder_prefix_ << current_folder_prefix << this->which_motion_dynamic_am_i_;
	return tmp_stream.str();

}

std::string MotionDynamic::get_temp_folder_name()
{
	return this->temp_folder_name_;
}

bool MotionDynamic::make_temp_folder()
{
	try
	{
		std::cout << "Generating temporary folder " << this->temp_folder_name_ << std::endl;
		boost::filesystem::path dir_to_make(this->temp_folder_name_.c_str());
		bool folder_creation_worked = boost::filesystem::create_directories(dir_to_make);
		return folder_creation_worked;

	}
	catch(boost::system::error_code& e)
	{
		std::cout << e.message() << std::endl;
		throw e;	
	}
}

bool MotionDynamic::delete_temp_folder()
{
	try
	{
		boost::filesystem::path dir_to_del( this->temp_folder_name_.c_str() );
		
		if( boost::filesystem::exists(dir_to_del) )
		{
			std::cout << "Deleting temporary folder " << this->temp_folder_name_ << std::endl;

			bool folder_deletion_worked = boost::filesystem::remove_all(dir_to_del);
			return folder_deletion_worked;
		}
		else
		{
			std::cout << "Folder " << this->temp_folder_name_ << " does not exist. Deletion omitted" << std::endl;
			return false;
		}

	}
	catch(boost::system::error_code& e)
	{
		std::cout << e.message() << std::endl;
		throw e;	
	}
;
}


void MotionDynamic::set_displacment_fields( ISMRMRD::NDArray< DataTypeMotionFields >& motion_fields)
{
	using namespace ISMRMRD;

	const size_t* dimensions = motion_fields.getDims();

	std::string const output_name_mvf = "/media/sf_SharedFolder/CCPPETMR/motionfields_in_memory_10x3x64x64x64";

	data_io::write_raw(output_name_mvf, motion_fields.begin(), motion_fields.getNumberOfElements() );

	for( int i=0; i<7; i++)
		std::cout << "dim_" << i << "=" << dimensions[i] << std::endl;

	size_t const Nt = dimensions[0];
	size_t const Nv = dimensions[1];
	size_t const Nz = dimensions[2];
	size_t const Ny = dimensions[3];
	size_t const Nx = dimensions[4];

	for(size_t nt=0; nt<Nt; nt++)
	{
		
		Image<DataTypeMotionFields> img(dimensions[4],dimensions[3], dimensions[2], dimensions[1]);
 		
 		for(uint16_t  nv= 0; nv<Nv ; nv++)
		for(uint16_t  nz= 0; nz<Nz ; nz++)
		for(uint16_t  ny= 0; ny<Ny ; ny++)
		for(uint16_t  nx= 0; nx<Nx ; nx++)
		{
			// size_t const lin_index = ((((Nt-1 -nt)*Nv + Nv-1 -nv)*Nz + Nz-1 - nz)*Ny + Ny-1 - ny)*Nx + Nx-1 - nx;
			size_t const lin_index = (((nt*Nv + nv)*Nz + nz)*Ny + ny)*Nx + nx;
			img(nx,ny,nz,nv) = 	  *(motion_fields.begin() + lin_index);
		}
		this->displacment_fields_.push_back(img);
	}
}

void MotionDynamic::write_temp_displacements_fields()
{
	bool const temp_folder_creation_successful = this->make_temp_folder();

	if( temp_folder_creation_successful )
	{
		for(int i=0; i<this->displacment_fields_.size(); i++)
		{

			// std::cout << "Writing MVF #: " << i <<  std::endl;
			std::stringstream temp_filename_mvf;
			temp_filename_mvf << this->get_temp_folder_name() << this->temp_mvf_prefix_ << i;

			data_io::write_ISMRMRD_Image_to_Analyze<DataTypeMotionFields> (temp_filename_mvf.str(), this->displacment_fields_[i]);
			this->temp_mvf_filenames_.push_back(temp_filename_mvf.str());
		}
	}
	else
		throw std::runtime_error("The parent directory generation failed. Give a path to which thou hast access rights. Or maybe the directory already exists. This is dangerous. Then you should definitely choose a different temporary folder name.");
}

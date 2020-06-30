/*
SyneRBI Synergistic Image Reconstruction Framework (SIRF)
Copyright 2019 - 2020 Rutherford Appleton Laboratory STFC
Copyright 2019 - 2020 University College London

This is software developed for the Collaborative Computational
Project in Synergistic Reconstruction for Biomedical Imaging (formerly CCP PETMR)
(http://www.ccpsynerbi.ac.uk/).

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

*/

/*!
\file
\ingroup Gadgetron Extensions
\brief MR related C++ tests

\author Johannes Mayer
\author SyneRBI
*/

#include <iostream>
#include <cstdlib>

#include "sirf/Gadgetron/gadgetron_data_containers.h"
#include "sirf/Gadgetron/gadgetron_x.h"

#include "mrtest_auxiliary_funs.h"

using namespace sirf;

bool test_get_kspace_order(const std::string& fname_input)
{
    try
    {
        std::cout << "Running test " << __FUNCTION__ << std::endl;

        sirf::AcquisitionsVector av_slice;
        av_slice.read(fname_input);
        av_slice.sort();

        auto kspace_sorting_slice = av_slice.get_kspace_order();

        return true;

    }
    catch( std::runtime_error const &e)
    {
        std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
}

bool test_get_subset(const std::string& fname_input)
{
    try
    {
        std::cout << "Running test " << __FUNCTION__ << std::endl;

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
        throw;
    }
}

bool test_CoilSensitivitiesVector_calculate(const std::string& fname_input)
{
    int const num_iter = 3;

    sirf::GadgetronImagesVector ic;

    for(int i=0; i<num_iter; ++i)
    {
        CFImage img;

        void* vptr_img = new CFImage(img);// god help me I don't trust this!
        ImageWrap iw(ISMRMRD::ISMRMRD_DataTypes::ISMRMRD_CXFLOAT, vptr_img);

        ic.append(iw);
    }

    return ic.number() == num_iter;
}

bool test_memory_safety_preprocessing(const std::string& fname_input)
{
    // this test is for valgrind use only to see if there is memory trouble when connecting to Gadgetron
    try
    {
        std::cout << "Running test " << __FUNCTION__ << std::endl;

     sirf::AcquisitionsVector mr_rawdata;
     mr_rawdata.read(fname_input);

     preprocess_acquisition_data(mr_rawdata);

        CoilSensitivitiesVector csv;
        csv.set_csm_smoothness(50);
        csv.calculate(av);

}

bool test_compute_coilmaps(const std::string& fname_input)
{
    try
    {
        std::cout << "Running test " << __FUNCTION__ << std::endl;

        sirf::AcquisitionsVector mr_rawdata;
        mr_rawdata.read(fname_input);

        preprocess_acquisition_data(mr_rawdata);
        mr_rawdata.sort();

        sirf::CoilSensitivitiesAsImages csm;

        size_t const ks = 7;
        size_t const kz = 5;
        size_t const power = 7;
        csm.set_csm_gadget_params(ks,kz,power);

        csm.compute(mr_rawdata);

        for(int i=0; i<csm.items(); ++i)
        {
            std::stringstream fname_output;
            fname_output << "output_" << __FUNCTION__ << "_csm_ks_" << ks << "_kz_" << kz << "power_" << power << "_"  <<  i;
            CFImage csm_img = csm.get_csm_as_CFImage(i);
            write_cfimage_to_raw(fname_output.str(), csm_img);
        }

        return true;
    }
    catch( std::runtime_error const &e)
    {
        std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
}

bool test_bwd(const std::string& fname_input)
{
    try
    {
       std::cout << "Running test " << __FUNCTION__ << std::endl;

        sirf::AcquisitionsVector mr_rawdata;
        mr_rawdata.read(fname_input);

        preprocess_acquisition_data(mr_rawdata);
        mr_rawdata.sort();

        sirf::GadgetronImagesVector img_vec;
        sirf::MRAcquisitionModel acquis_model;

        sirf::CoilSensitivitiesAsImages csm;
        csm.compute(mr_rawdata);

        auto sptr_encoder = std::make_shared<sirf::CartesianFourierEncoding>(sirf::CartesianFourierEncoding());
        acquis_model.set_encoder(sptr_encoder);

        acquis_model.bwd(img_vec, csm, mr_rawdata);

        for(int i=0; i<img_vec.items(); ++i)
        {
            std::stringstream fname_output;
            fname_output << "output_" << __FUNCTION__ << "_image_" << i;
            write_cfimage_to_raw(fname_output.str(), img_vec.image_wrap(i));
        }

        return true;

    }
    catch( std::runtime_error const &e)
    {
        std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
}

int main (int argc, char* argv[])
{
	try{

        std::string SIRF_PATH;
        if (argc==1)
            SIRF_PATH = getenv("SIRF_PATH");
        else
            SIRF_PATH = argv[1];

        std::string simul_data_path = SIRF_PATH + "/data/examples/MR/simulated_MR_2D_cartesian.h5";
        std::string real_data_path = SIRF_PATH + "/data/examples/MR/CV_2D_Stack_144.h5";

        std::vector<bool> test_results;
        //test_results.push_back(test_GRPETrajectoryPrep_set_trajectory(simul_data_path));
        //test_results.push_back(test_apply_combine_coil_sensitivities());
        //test_results.push_back(test_get_kspace_order(simul_data_path));
        //test_results.push_back(test_get_subset(simul_data_path));
        //test_results.push_back(test_append_image_wrap());
        //test_results.push_back(test_memory_safety_preprocessing(simul_data_path));
        test_results.push_back(test_compute_coilmaps(simul_data_path));
//        test_results.push_back(test_bwd(simul_data_path));

        bool all_tests_successful = std::accumulate(std::begin(test_results), std::end(test_results), true, std::multiplies<bool>());

        if(all_tests_successful)
            return EXIT_SUCCESS;
        else
        {
            for(int i=0; i<test_results.size(); ++i)
            {
                std::cout << "Test Result #" << i <<" = " << test_results[i] << std::endl;
            }
            return EXIT_FAILURE;
        }
	}
    catch(const std::exception &error) {
        std::cerr << "\nHere's the error:\n\t" << error.what() << "\n\n";
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}


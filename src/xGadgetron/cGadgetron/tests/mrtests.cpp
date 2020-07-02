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
\ingroup Gadgetron Tests

\author Evgueni Ovtchinnikov
\author Richard Brown
\author SyneRBI
*/

#include <iostream>
#include <cstdlib>
#include <vector>
#include <numeric>

#include "mrtest_auxiliary_funs.h"

#include "sirf/Gadgetron/encoding.h"
#include "sirf/Gadgetron/gadgetron_data_containers.h"
#include "sirf/Gadgetron/gadgetron_x.h"

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
        throw;
    }
}


bool test_GRPETrajectoryPrep_set_trajectory(const std::string& fname_input)
{
    try
    {
        std::cout << "Running test " << __FUNCTION__ << std::endl;
        sirf::GRPETrajectoryPrep rpe_tp;

        std::string const fname_output = "res_phant_cv_rpe_itl_gc_ismrmrd.h5";

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
        throw;
    }
}

bool test_apply_combine_coil_sensitivities( void )
{
    try
    {
        std::cout << "Running test " << __FUNCTION__ << std::endl;

        size_t const Nx=64;
        size_t const Ny=64;
        size_t const Nz=64;
        size_t const Nc=4;

        CFImage unit_img  = CFImage(Nx,Ny,Nz,1);

        for(size_t i=0;i<unit_img.getNumberOfDataElements();++i)
        {
            *(unit_img.begin()+i) = complex_float_t(1.f, 0.f);
        }
        sirf::GadgetronImagesVector giv, miv;
        giv.append(unit_img);

        // generate artifical coil maps
        size_t const num_dat_pts = Nx*Ny*Nz*Nc;
        std::vector<float> csm_real(num_dat_pts);
        std::vector<float> csm_imag(num_dat_pts);

        for(size_t nc=0; nc<Nc; nc++)
        for(size_t nz=0; nz<Nz; nz++)
        for(size_t ny=0; ny<Ny; ny++)
        for(size_t nx=0; nx<Nx; nx++)
        {
            size_t access_idx = (((nc*Nz + nz)*Ny + ny)*Nx + nx);
            csm_real.at(access_idx) = nx;
            csm_imag.at(access_idx) = nc;
        }

//        CoilSensitivitiesAsImages csm;
//        csm.append_csm(Nx, Ny, Nz, Nc, &csm_real[0], &csm_imag[0]);

//        csm.apply_coil_sensitivities(miv, giv);
//        csm.combine_coils(giv, miv);

        throw std::runtime_error("Please write teh last part of the test again");

        return true;

    }
    catch( std::runtime_error const &e)
    {
        std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
}

bool test_get_kspace_order(const std::string& fname_input)
{
    try
    {
        std::cout << "Running test " << __FUNCTION__ << std::endl;

        sirf::AcquisitionsVector av;
        av.read(fname_input);
        av.sort();

        auto kspace_sorting = av.get_kspace_order();

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

bool test_append_image_wrap()
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

     return true;
    }
    catch( std::runtime_error const &e)
    {
        std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }

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

        sirf::CoilSensitivitiesVector csm;

        size_t const ks = 7;
        size_t const kz = 5;
        size_t const power = 7;
        csm.set_csm_gadget_params(ks,kz,power);

        csm.calculate_csm(mr_rawdata);

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


bool test_CoilSensitivitiesVector_calculate(const std::string& fname_input)
{
    try
    {
        std::cout << "Running test " << __FUNCTION__ << std::endl;

        sirf::AcquisitionsVector av;
        av.read(fname_input);

        sirf::preprocess_acquisition_data(av);

        CoilSensitivitiesVector csv;
        csv.set_csm_smoothness(50);
        csv.calculate_csm(av);

        std::cout << "We have " << csv.items() << " coilmaps" << std::endl;

        for(int i=0; i<csv.items(); ++i)
        {
            gadgetron::shared_ptr<ImageWrap> sptr_iw = csv.sptr_image_wrap(i);

            std::stringstream fname_out;
            fname_out << "output_" << __FUNCTION__ << "_" << i;

            sirf::write_cfimage_to_raw(fname_out.str(), *sptr_iw);
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

bool test_CoilSensitivitiesVector_get_csm_as_cfimage(const std::string& fname_input)
{
    try
    {
        std::cout << "Running test " << __FUNCTION__ << std::endl;

        sirf::AcquisitionsVector av;
        av.read(fname_input);

        sirf::preprocess_acquisition_data(av);

        CoilSensitivitiesVector csv;
        csv.calculate_csm(av);

        std::cout << "We have " << csv.items() << " coilmaps" << std::endl;

        for(int i=0; i<csv.items(); ++i)
        {
            CFImage img = csv.get_csm_as_CFImage(i);

            std::stringstream fname_out;
            fname_out << "output_" << __FUNCTION__ << "_" << i;

            sirf::write_cfimage_to_raw(fname_out.str(), img);
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

        sirf::CoilSensitivitiesVector csm;
        csm.calculate_csm(mr_rawdata);

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
//        test_results.push_back(test_GRPETrajectoryPrep_set_trajectory(simul_data_path));
//        test_results.push_back(test_apply_combine_coil_sensitivities());
//        test_results.push_back(test_get_kspace_order(simul_data_path));
//        test_results.push_back(test_get_subset(simul_data_path));
//        test_results.push_back(test_append_image_wrap());
//        test_results.push_back(test_memory_safety_preprocessing(simul_data_path));
        test_results.push_back(test_compute_coilmaps(simul_data_path));
        test_results.push_back(test_CoilSensitivitiesVector_calculate(simul_data_path));
        test_results.push_back(test_CoilSensitivitiesVector_get_csm_as_cfimage(simul_data_path));

        test_results.push_back(test_bwd(simul_data_path));

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
}

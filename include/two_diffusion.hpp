#ifndef _DIFFUSION_H_
#define _DIFFUSION_H_

#include <H5Cpp.h>
#include<vector>
#include <sys/stat.h>
#include "TextParser.h"
#include<string>
#include<sstream>
#include<fstream>
#include<iostream>
#include<omp.h>
#include<set>


class twodimensinal_diffusion{
    private:
        double diffusion_coefficient;
        int numOfBoundaryNode;
        std::string material_judge;
        int numofomp;
        std::string gauss_setting, outputDir;
        std::string node_file, element_file, boundary_file, phi_file;
        std::vector<std::vector<double>> D;
        std::vector<std::vector<double>> mass;
        std::vector<double> boundary_value;
        std::vector<std::vector<double>> gauss;
        std::set<int> boundary_node_judge;
        std::string output_h5_name;
    public:
        double dt;
        int time, output_interval, numOfNode, numOfElm;
        double coupling_coefficient_vi;
        double coupling_coefficient_vc;
        double coupling_coefficient_ci;
        std::vector<double> mass_centralization;
        std::vector<int> boundary_node;
        std::vector<int> update_point;
        std::vector<double> C;
        std::vector<double> phi, phi_v;
        std::vector<std::vector<double>> node;
        std::vector<std::vector<int>> element;
        TextParser tp;
        twodimensinal_diffusion(std::string mat)
        {
            material_judge = mat;
        }
        void export_vtu(const std::string &file, std::string judge, std::vector<double> output_value);
        void input_info(std::string input_file);
        void input_phi();
        int CountNumbersOfTextLines(std::string &filePath);
        void read_geometry();
        void boundary_initialize();
        void calc_dxdr(int ic, std::vector<std::vector<double>> node, std::vector<std::vector<int>> element, std::vector<std::vector<double>> &dxdr, std::vector<std::vector<double>> dNdr);
        void calc_dNdx(std::vector<std::vector<double>> &dNdx, std::vector<std::vector<double>> dNdr, std::vector<std::vector<double>> drdx);
        void calc_inverse_matrix_2x2(std::vector<std::vector<double>> dxdr, std::vector<std::vector<double>> &drdx);
        void gauss_point_setting();
        void matrix_initialize();
        void calc_matrix();
        void boundary_setting(double time, std::vector<double> Q_cv, std::vector<double> Q_iv);
        void time_step(std::vector<double> Q1, std::vector<double> Q2, double time);

        void dump(int ic);
        void hdf5_dump(int ic);
        void exportHDF5_double_1D(H5::H5File &file, const std::string &dataName, std::vector<double> i_data, int i_dim);
        double access_c(int ic);
        void transform_point_data_to_cell_data(std::vector<double> &element_C, std::vector<double> C);
        void transform_point_data_to_cell_data_phi(std::vector<double> &phiC, std::vector<double> C);
};
#endif

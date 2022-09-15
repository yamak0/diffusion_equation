#ifndef _DIFFUSION_H_
#define _DIFFUSION_H_

#include<vector>
#include <sys/stat.h>
#include "TextParser.h"
#include<string>
#include<sstream>
#include<fstream>
#include<iostream>
#include<omp.h>

class twodimensinal_diffusion{
    private:
        double dt;
        double diffusion_coefficient;
        int numOfElm, numOfBoundaryNode;
        std::string material_judge;
        std::string gauss_setting, outputDir;
        std::string node_file, element_file, boundary_file, phi_file, phi_visualize_file;
        std::vector<double> C;
        std::vector<std::vector<double>> D;
        std::vector<double> mass_centralization;
        std::vector<int> boundary_node;
        std::vector<double> boundary_value;
        std::vector<std::vector<double>> gauss;
    public:
        int time, output_interval, numOfNode;
        double coupling_coefficient;
        double coupling_coefficient_v;
        double coupling_coefficient_s;
        std::vector<double> phi, phi_v;
        std::vector<std::vector<double>> node;
        std::vector<std::vector<int>> element;
        TextParser tp;
        twodimensinal_diffusion(std::string mat)
        {
            material_judge = mat;
        }
        void export_vtu(const std::string &file);
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
        void boundary_setting();
        void time_step(std::vector<double> diff);
        void dump(int ic);
        double access_c(int ic);
};
#endif

#ifndef _DIFFUSION_H_
#define _DIFFUSION_H_

#include<vector>
#include <sys/stat.h>
#include "TextParser.h"
#include<string>
#include<sstream>
#include<fstream>
#include<iostream>

class twodimensinal_diffusion{
    private:
        double dt;
        int time, output_interval;
        double diffusion_coefficient;
        int numOfNode, numOfElm, numOfBoundaryNode;
        std::string gauss_setting, outputDir;
        std::string node_file, element_file, boundary_file;
        std::vector<std::vector<double>> node;
        std::vector<std::vector<int>> element;
        std::vector<double> C;
        std::vector<std::vector<double>> D;
        std::vector<double> mass_centralization;
        std::vector<int> boundary_node;
        std::vector<double> boundary_value;
        std::vector<std::vector<double>> gauss;
    public:
        TextParser tp;
        void export_vtu(const std::string &file);
        void input_info(std::string input_file);
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
        void time_step();
};
#endif

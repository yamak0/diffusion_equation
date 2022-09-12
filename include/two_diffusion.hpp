#ifndef _DIFFUSION_H_
#define _DIFFUSION_H_

#include<vector>

class twodimensinal_diffusion{
    private:
        double dt;
        double k;
        std::vector<std::vector<double>> node;
        std::vector<std::vector<int>> element;

    public:
        void calc_dxdr(int ic, std::vector<std::vector<double>> node, std::vector<std::vector<int>> element, std::vector<std::vector<double>> &dxdr, std::vector<std::vector<double>> dNdr);
        void calc_inverse_matrix_2x2(std::vector<std::vector<double>> dxdr, std::vector<std::vector<double>> &drdx);
};
#endif

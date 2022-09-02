#ifndef _DIFFUSION_H_
#define _DIFFUSION_H_

#include<vector>
#include<iostream>
#include<fstream>
#include<string>
#include <sys/stat.h>

class onedimensinal_diffusion{
    private:
        double element_length;
        double dt;
        double diffusion_coefficient;
        std::vector<std::vector<int>> element;
        std::vector<double> phi;
        std::vector<double> x;
        std::vector<std::vector<double>> mass;
        std::vector<std::vector<double>> K;
        std::vector<int> b_n;
        std::vector<double> i_b_n;
    
    public:
        void calc_mass_matrix();
        void calc_K_matrix();
        void set_boundary_node(int node_number, double value);
        void time_step();
        void dump(int step, std::string name);

        onedimensinal_diffusion(int numofelem, double dx, double time_step, double D)
        {
            element_length=dx;
            dt=time_step;
            diffusion_coefficient = D;
            element.resize(numofelem);
            phi.resize(numofelem+1);
            x.resize(numofelem+1);
            K.resize(numofelem+1);
            mass.resize(numofelem+1);

            for(int i=0; i<K.size(); i++){
                K[i].resize(numofelem+1);
                for(int j=0; j<K[i].size(); j++){
                    K[i][j]=0.0;
                }
            }

            for(int i=0; i<mass.size(); i++){
                mass[i].resize(numofelem+1);
                for(int j=0; j<mass[i].size(); j++){
                    mass[i][j]=0.0;
                }
            }
            
            for(int i=0; i<element.size(); i++){
                element[i].resize(2);
                element[i][0]=i;
                element[i][1]=i+1;
            }
            
            for(int i=0; i<x.size(); i++){
                x[i]=dx*i;
            }
        }
};

#endif
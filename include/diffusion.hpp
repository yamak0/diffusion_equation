#ifndef _DIFFUSION_H_
#define _DIFFUSION_H_

#include<vector>

class onedimensinal_diffusion{
    private:
        std::vector<std::vector<int>> element;
        std::vector<double> phi;
        std::vector<double> x;

    public:
        onedimensinal_diffusion(int numofelem, double dx)
        {
            element.resize(numofelem);
            phi.resize(numofelem+1);
            x.resize(numofelem+1);
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
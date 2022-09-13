#include<iostream>
#include<vector>
#include<fstream>
#include<string>
#include<cmath>
#include"two_diffusion.hpp"
#include"shapefunction.hpp"
using namespace std;

void set_field(vector<vector<double>> &node, vector<vector<int>> &element, int n, double dx)
{
    node.resize((n+1)*(n+1));
    for(int i=0; i<node.size(); i++){
        node[i].resize(2);
    }
    for(int i=0; i<n+1; i++){
        for(int j=0; j<n+1; j++){
            node[i*(n+1)+j][0]= j*dx;
            node[i*(n+1)+j][1]= i*dx;
        }
    }
    element.resize(n*n);
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            element[i*n+j].push_back(i*(n+1)+j);
            element[i*n+j].push_back(i*(n+1)+j+1);
            element[i*n+j].push_back(i*(n+1)+j+(n+2));
            element[i*n+j].push_back(i*(n+1)+j+(n+1));   
        }
    }

    ofstream ofs("node.dat");
    for(int i=0; i<node.size(); i++){
      ofs << node[i][0] << " " << node[i][1] << " " << node[i][2] << endl;
    }
    ofs.close();

    ofs.open("element.dat");
    for(int i=0; i<element.size(); i++){
      for(int j=0; j<element[i].size(); j++){
        ofs << element[i][j] << " ";
      }
      ofs << endl;
    }
    ofs.close();

    ofs.open("boundary.dat");
    for(int i=0; i<node.size(); i++){
      if(fabs(node[i][1]-0.0)<0.0001){
        ofs << i << " " << 1.0 << endl;
      }
    }
    ofs.close();
}

static void C2D4_N(std::vector<double> &N,const double &g1,const double &g2)
{
  N[0] = 2.5e-1 * (1e+0-g1) * (1e+0-g2);
  N[1] = 2.5e-1 * (1e+0+g1) * (1e+0-g2);
  N[2] = 2.5e-1 * (1e+0+g1) * (1e+0+g2);
  N[3] = 2.5e-1 * (1e+0-g1) * (1e+0+g2);
}

static void C2D4_dNdr(std::vector<std::vector<double>> &dNdr,const double &g1,const double &g2)
{
  dNdr[0][0] = -2.5e-1 * (1e+0-g2);
  dNdr[0][1] = -2.5e-1 * (1e+0-g1);
  dNdr[1][0] =  2.5e-1 * (1e+0-g2);
  dNdr[1][1] = -2.5e-1 * (1e+0+g1);
  dNdr[2][0] =  2.5e-1 * (1e+0+g2);
  dNdr[2][1] =  2.5e-1 * (1e+0+g1);
  dNdr[3][0] = -2.5e-1 * (1e+0+g2);
  dNdr[3][1] =  2.5e-1 * (1e+0-g1);
}

void calc_dxdr(int ic, vector<vector<double>> node, vector<vector<int>> element, vector<vector<double>> &dxdr, vector<vector<double>> dNdr)
{
  for(int k=0;k<2;k++){
    for(int l=0;l<2;l++){
      dxdr[k][l] = 0e0;
      for(int p=0;p<4;p++){
        dxdr[k][l] += dNdr[p][l] * node[element[ic][p]][k];
      }
    }
  }
}

void calc_inverse_matrix_2x2(vector<vector<double>> dxdr, vector<vector<double>> &drdx)
{
  double det = dxdr[0][0]*dxdr[1][1]-dxdr[0][1]*dxdr[1][0];
  drdx[0][0] = 1.0/det*dxdr[1][1];
  drdx[1][1] = 1.0/det*dxdr[0][0];
  drdx[0][1] = -1.0/det*dxdr[1][0];
  drdx[1][0] = -1.0/det*dxdr[0][1];
}

void calc_dNdx(vector<vector<double>> &dNdx, vector<vector<double>> dNdr, vector<vector<double>> drdx)
{
  for(int k=0; k<4; k++){
    for(int l=0; l<2; l++){
      dNdx[k][l] = 0.0;
      for(int p=0; p<2; p++){
        dNdx[k][l] += dNdr[k][p]*drdx[p][l];
      }
    }
  }
}

int main(int argc,char *argv[])
{
  //input argument
  if(argc!=2){
    printf("Invalid input. Please set tp file\n");
    return 1;
  }
  std::string input_file = argv[1];
  twodimensinal_diffusion Fluid;

  Fluid.input_info(input_file);
  Fluid.gauss_point_setting();
  Fluid.matrix_initialize();
  Fluid.calc_matrix();
  Fluid.time_step();
}
#include"two_diffusion.hpp"
#include"shapefunction.hpp"
#include <time.h>
using namespace std;

void twodimensinal_diffusion::boundary_initialize()
{
    string str,tmp;
    numOfBoundaryNode=CountNumbersOfTextLines(boundary_file);
    boundary_node.resize(numOfBoundaryNode);
    boundary_value.resize(numOfBoundaryNode);
    
    ifstream ifs(boundary_file);
    for(int i=0; i<numOfBoundaryNode; i++){
        getline(ifs,str);
        istringstream stream(str);
        for(int j=0;j<2;j++){
            getline(stream,tmp,' ');
            if(j==0) boundary_node[i] = stoi(tmp);
            if(j==1) boundary_value[i] = stod(tmp);
        }
    }
    ifs.close();   
}


void twodimensinal_diffusion::gauss_point_setting()
{
    if(gauss_setting=="square"){
        gauss.resize(4);
        for(int i=0; i<gauss.size(); i++){
            gauss[i].resize(2);
        }
    }

    gauss[0][0] = -0.577350296189626; gauss[0][1] = -0.577350296189626;
    gauss[1][0] = -0.577350296189626; gauss[1][1] = 0.577350296189626; 
    gauss[2][0] = 0.577350296189626; gauss[2][1] = -0.577350296189626; 
    gauss[3][0] = 0.577350296189626; gauss[3][1] = 0.577350296189626; 
}

void twodimensinal_diffusion::matrix_initialize()
{
    D.resize(node.size());
    for(int i=0; i<node.size(); i++){
        D[i].resize(node.size());
    }

    mass_centralization.resize(node.size());
}

void twodimensinal_diffusion::calc_dxdr(int ic, std::vector<std::vector<double>> node, std::vector<std::vector<int>> element, std::vector<std::vector<double>> &dxdr, std::vector<std::vector<double>> dNdr)
{
  int k,l,p;
  for(k=0;k<2;k++){
    for(l=0;l<2;l++){
      dxdr[k][l] = 0e0;
      for(p=0;p<4;p++){
        dxdr[k][l] += dNdr[p][l] * node[element[ic][p]][k];
      }
    }
  }
}

void twodimensinal_diffusion::calc_inverse_matrix_2x2(std::vector<std::vector<double>> dxdr, std::vector<std::vector<double>> &drdx)
{
  double det = dxdr[0][0]*dxdr[1][1]-dxdr[0][1]*dxdr[1][0];
  drdx[0][0] = 1.0/det*dxdr[1][1];
  drdx[1][1] = 1.0/det*dxdr[0][0];
  drdx[0][1] = -1.0/det*dxdr[1][0];
  drdx[1][0] = -1.0/det*dxdr[0][1];
}

void twodimensinal_diffusion::calc_dNdx(vector<vector<double>> &dNdx, vector<vector<double>> dNdr, vector<vector<double>> drdx)
{
  int k,l,p;
  for(k=0; k<4; k++){
    for(l=0; l<2; l++){
      dNdx[k][l] = 0.0;
      for(p=0; p<2; p++){
        dNdx[k][l] += dNdr[k][p]*drdx[p][l];
      }
    }
  }
}

void twodimensinal_diffusion::calc_matrix()
{
  int i,j,k,l,p;
  for(i=0; i<element.size(); i++){
    cout << i << endl;
    vector<double> N(4); 
    vector<vector<double>> dNdr(4, vector<double>(2));
    double volume=0.0;
    vector<vector<double>> element_D(4, vector<double>(4,0.0));
    vector<vector<double>> element_mass(4, vector<double>(4,0.0));
    vector<double> element_mass_centralization(4,0.0);
    
    for(j=0; j<4; j++){
      ShapeFunction2D::C2D4_N(N,gauss[j][0],gauss[j][1]);
      ShapeFunction2D::C2D4_dNdr(dNdr,gauss[j][0],gauss[j][1]);
      vector<vector<double>> dxdr(2, vector<double>(2)), drdx(2, vector<double>(2));
      calc_dxdr(i, node, element, dxdr, dNdr);
      calc_inverse_matrix_2x2(dxdr, drdx);
      double detJ = dxdr[0][0] * dxdr[1][1]  - dxdr[1][0] * dxdr[0][1];
      vector<vector<double>> dNdx(4, vector<double>(2, 0.0));
      calc_dNdx(dNdx, dNdr, drdx);
      
      for(k=0; k<4; k++){
        for(l=0; l<4; l++){
          for(p=0; p<2; p++){
            element_D[k][l] += dNdx[k][p]*dNdx[l][p];
          }
        }
      }
      
      for(k=0; k<4; k++){
        for(l=0; l<4; l++){
          element_mass[k][l] = N[k] * N[l];
        }
      }
      
      for(k=0; k<4; k++){
        for(l=0; l<4; l++){
          element_mass_centralization[k] += element_mass[k][l];
        }
      }
      
      for(k=0; k<4; k++){
        for(l=0; l<4; l++){
          D[element[i][k]][element[i][l]] += diffusion_coefficient * element_D[k][l] * phi[i] * detJ;
        }
        mass_centralization[element[i][k]] += element_mass_centralization[k] * detJ;
      }
    }
  }
}

void twodimensinal_diffusion::boundary_setting()
{
  for(int i=0; i<numOfBoundaryNode; i++){
    C[boundary_node[i]] = boundary_value[i]*phi_node[boundary_node[i]];
  }
  phi_node.resize(numOfNode);
  for(int i=0; i<boundary_node.size(); i++){
    phi_node[boundary_node[i]] = 1.0;
  }
}

void twodimensinal_diffusion::time_step(vector<double> diff)
{
  vector<double> R(numOfNode);
  vector<double> DC(numOfNode,0.0);
  vector<double> DcR(numOfNode,0.0);
  vector<double> MDcR(numOfNode,0.0);
  vector<double> MDC(node.size(),0.0);
    
  int i,j;
  #pragma omp parallel for private(j)
  for(i=0; i<numOfNode; i++){
    for(j=0; j<node.size(); j++){
        DC[i] += D[i][j] * (C[j]);
    }
  }
  #pragma omp parallel for
  for(i=0; i<numOfNode; i++){
    R[i] = mass_centralization[i]*diff[i];
    DcR[i] = DC[i]-R[i];
    MDcR[i] = 1.0/mass_centralization[i]*DcR[i];
    C[i] = C[i] - dt * MDcR[i];
  }  
  
  if(material_judge=="V") boundary_setting();
}

double twodimensinal_diffusion::access_c(int ic)
{
  return C[ic];
}

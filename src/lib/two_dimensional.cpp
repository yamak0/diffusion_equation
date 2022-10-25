#include"two_diffusion.hpp"
#include"shapefunction.hpp"
using namespace std;
using namespace H5;



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
            if(j==0){
              boundary_node[i] = stoi(tmp);
              boundary_node_judge.insert(stoi(tmp));
            }
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

    mass.resize(node.size());
    for(int i=0; i<node.size(); i++){
        mass[i].resize(node.size());
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
  #pragma omp parallel for
  for(int i=0; i<element.size(); i++){
    vector<double> N(4); 
    vector<vector<double>> dNdr(4, vector<double>(2));
    double volume=0.0;
    vector<vector<double>> element_D(4, vector<double>(4,0.0));
    vector<vector<double>> element_mass(4, vector<double>(4,0.0));
    vector<double> element_mass_centralization(4,0.0);
    for(int j=0; j<4; j++){
      ShapeFunction2D::C2D4_N(N,gauss[j][0],gauss[j][1]);
      ShapeFunction2D::C2D4_dNdr(dNdr,gauss[j][0],gauss[j][1]);
      vector<vector<double>> dxdr(2, vector<double>(2)), drdx(2, vector<double>(2));
      calc_dxdr(i, node, element, dxdr, dNdr);
      calc_inverse_matrix_2x2(dxdr, drdx);
      double detJ = dxdr[0][0] * dxdr[1][1]  - dxdr[1][0] * dxdr[0][1];
      vector<vector<double>> dNdx(4, vector<double>(2, 0.0));
      calc_dNdx(dNdx, dNdr, drdx);
      for(int k=0; k<4; k++){
        for(int l=0; l<4; l++){
          for(int p=0; p<2; p++){
            element_D[k][l] += dNdx[k][p]*dNdx[l][p];
          }
        }
      }
      
      for(int k=0; k<4; k++){
        for(int l=0; l<4; l++){
          element_mass[k][l] = N[k] * N[l];
        }
      }
      
      for(int k=0; k<4; k++){
        for(int l=0; l<4; l++){
          D[element[i][k]][element[i][l]] += diffusion_coefficient * element_D[k][l] * phi[i] * detJ;
          mass[element[i][k]][element[i][l]] += element_mass[k][l] * phi[i] * detJ;
        }
      }
    }
  }

  for(int i=0; i<mass.size(); i++){
    mass_centralization[i]=0.0;
    for(int j=0; j<mass[i].size(); j++){
      mass_centralization[i] += mass[i][j];
    }
  }
}

void twodimensinal_diffusion::boundary_setting(double time_t, vector<double> Q_cv, vector<double> Q_iv)
{
  if(material_judge=="F" || material_judge== "S"){
    for(int i=0; i<numOfBoundaryNode; i++){
      C[boundary_node[i]] = boundary_value[i];
    }
  }
  else{
    vector<pair<double, int>> cell_to_point_cv(numOfNode);
    vector<pair<double, int>> cell_to_point_iv(numOfNode);
    for(int i=0; i<numOfElm; i++){
      for(int j=0; j<element[i].size(); j++){
        cell_to_point_cv[element[i][j]].first += Q_cv[i];
        cell_to_point_iv[element[i][j]].first += Q_iv[i];
        cell_to_point_cv[element[i][j]].second+=1;
        cell_to_point_iv[element[i][j]].second+=1;
      }
    }
    for(int i=0; i<numOfNode; i++){
      cell_to_point_cv[i].first /= cell_to_point_cv[i].second;
      cell_to_point_iv[i].first /= cell_to_point_iv[i].second;
    }
    double sum_1=0.0;
    double sum_2=0.0;
    for(int i=0; i<numOfBoundaryNode; i++){
      double param = 0.020, sigma=0.8, mu=-20.0;
      double time_tmp = ((time_t-mu)-(time*dt-mu)/2.0)*param;
      C[boundary_node[i]] = (exp(-pow(time_tmp,2.0)/(2.0*pow(sigma,2.0)))-dt*cell_to_point_cv[boundary_node[i]].first-dt*cell_to_point_iv[boundary_node[i]].first)*boundary_value[i];
    }
  }
}

void twodimensinal_diffusion::time_step(vector<double> Q1, vector<double> Q2, double time_t)
{
  vector<double> DC(numOfNode,0.0);
  vector<double> DcR(numOfNode,0.0);
  vector<double> MDcR(numOfNode,0.0);
  vector<double> MDC(node.size(),0.0);

  if(material_judge=="F"){
    vector<pair<double, int>> cell_to_point_cv(numOfNode);
    vector<pair<double, int>> cell_to_point_ci(numOfNode);
    for(int i=0; i<numOfElm; i++){
      for(int j=0; j<element[i].size(); j++){
        cell_to_point_cv[element[i][j]].first += Q1[i];
        cell_to_point_ci[element[i][j]].first += Q2[i];
        cell_to_point_cv[element[i][j]].second+=1;
        cell_to_point_ci[element[i][j]].second+=1;
      }
    }
    for(int i=0; i<numOfNode; i++){
      cell_to_point_cv[i].first /= double(cell_to_point_cv[i].second);
      cell_to_point_ci[i].first /= double(cell_to_point_ci[i].second);
    }
    #pragma omp parallel for
    for(int i=0; i<numOfNode; i++){
      for(int j=0; j<numOfNode; j++){
          DC[i] += D[i][j] * (C[j]);
      }
    }
    #pragma omp parallel for
    for(int i=0; i<numOfNode; i++){
      DcR[i] = DC[i]-(cell_to_point_cv[i].first+cell_to_point_ci[i].first);
      if(boundary_node_judge.find(i)==boundary_node_judge.end()){
        MDcR[i] = 1.0/mass_centralization[i]*DcR[i];
        C[i] = C[i] - dt * MDcR[i];
      }
    }  
    boundary_setting(time_t, Q1, Q2);
  }
  else if(material_judge=="S"){
    vector<pair<double, int>> cell_to_point_iv(numOfNode);
    vector<pair<double, int>> cell_to_point_ci(numOfNode);
    for(int i=0; i<numOfElm; i++){
      for(int j=0; j<element[i].size(); j++){
        cell_to_point_iv[element[i][j]].first += Q1[i];
        cell_to_point_ci[element[i][j]].first += Q2[i];
        cell_to_point_iv[element[i][j]].second+=1;
        cell_to_point_ci[element[i][j]].second+=1;
      }
    }
    for(int i=0; i<numOfNode; i++){
      cell_to_point_iv[i].first /= (double)cell_to_point_iv[i].second;
      cell_to_point_ci[i].first /= (double)cell_to_point_ci[i].second;
    }
    #pragma omp parallel for
    for(int i=0; i<numOfNode; i++){
      for(int j=0; j<numOfNode; j++){
          DC[i] += D[i][j] * (C[j]);
      }
    }
    #pragma omp parallel for
    for(int i=0; i<numOfNode; i++){
      DcR[i] = DC[i]-(cell_to_point_iv[i].first+cell_to_point_ci[i].first);
      if(boundary_node_judge.find(i)==boundary_node_judge.end()){
        MDcR[i] = 1.0/mass_centralization[i]*DcR[i];
        C[i] = C[i] - dt * MDcR[i];
      }
    }  
    boundary_setting(time_t, Q1, Q2);
  }
}

double twodimensinal_diffusion::access_c(int ic)
{
  return C[ic];
}

void twodimensinal_diffusion::transform_point_data_to_cell_data(std::vector<double> &element_C, std::vector<double> C)
{
  for(int i=0; i<numOfElm; i++){
    double tmp_C=0.0;
    for(int j=0; j<element[i].size(); j++){
      tmp_C+=C[element[i][j]];
    }
    tmp_C/=element[i].size();
    element_C[i]=tmp_C;
  } 
}

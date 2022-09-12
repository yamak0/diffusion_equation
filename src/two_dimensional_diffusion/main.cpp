#include<iostream>
#include<vector>
#include<fstream>
#include<string>

using namespace std;

void export_vtu(const std::string &file, vector<vector<double>> x, vector<vector<int>> element, vector<double> p)
{
  FILE *fp;
  if ((fp = fopen(file.c_str(), "w")) == NULL)
  {
    cout << file << " open error" << endl;
    exit(1);
  }

  fprintf(fp, "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\">\n");
  fprintf(fp, "<UnstructuredGrid>\n");
  fprintf(fp, "<Piece NumberOfPoints= \"%d\" NumberOfCells= \"%d\" >\n", x.size(), element.size());
  fprintf(fp, "<Points>\n");
  int offset = 0;
  fprintf(fp, "<DataArray type=\"Float64\" Name=\"Position\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%d\"/>\n",offset);
  offset += sizeof(int) + sizeof(double) * x.size() * 3;
  fprintf(fp, "</Points>\n");

  fprintf(fp, "<Cells>\n");
  fprintf(fp, "<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n");
  for (int i = 0; i < element.size(); i++){
    for (int j = 0; j < element[i].size(); j++) fprintf(fp, "%d ", element[i][j]);
    fprintf(fp, "\n");
  }
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n");
  int num = 0;
  for (int i = 0; i < element.size(); i++)
  {
    num += element[i].size();
    fprintf(fp, "%d\n", num);
  }
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
  for (int i = 0; i < element.size(); i++) fprintf(fp, "%d\n", 9);
    
    //fprintf(fp, "%d\n", element[i].meshType);

  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "</Cells>\n");

  fprintf(fp, "<PointData>\n");
  //fprintf(fp, "<DataArray type=\"Float64\" Name=\"velocity[m/s]\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%d\"/>\n",offset);
  //offset += sizeof(int) + sizeof(double) * numOfNode * 3;

  fprintf(fp, "<DataArray type=\"Float64\" Name=\"pressure[Pa]\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%d\"/>\n",offset);
  offset += sizeof(int) + sizeof(double) * x.size();

  //fprintf(fp, "<DataArray type=\"Float64\" Name=\"wall_share_stress[Pa]\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%d\"/>\n",offset);
  //offset += sizeof(int) + sizeof(double) * numOfNode * 3;
  fprintf(fp, "</PointData>\n");

  fprintf(fp, "<CellData>\n");
  fprintf(fp, "</CellData>\n");
  fprintf(fp, "</Piece>\n");
  fprintf(fp, "</UnstructuredGrid>\n");
  fprintf(fp, "<AppendedData encoding=\"raw\">\n");
  fprintf(fp, "_");
  fclose(fp);

  fstream ofs;
  ofs.open(file.c_str(), ios::out | ios::app | ios_base::binary);
  double *data_d = new double[x.size()*3];
  num = 0;
  int size=0;
  for (int ic = 0; ic < x.size(); ic++){
    //for(int j=0;j<2;j++){
      data_d[num] = x[ic][0];
      num++;
      data_d[num] = x[ic][1];
      num++;
      data_d[num] = 0.0;
      num++;
    //}
  }
  size=sizeof(double)*x.size()*3;
  ofs.write((char *)&size, sizeof(size));
  ofs.write((char *)data_d, size);

  //num=0;
  //for (int ic = 0; ic < numOfNode; ic++){
  //    data_d[num]   = u(ic);
  //    data_d[num+1] = v(ic);
  //    data_d[num+2] = w(ic);
  //    num=num+3;
  //}
  //size=sizeof(double)*numOfNode*3;
  //ofs.write((char *)&size, sizeof(size));
  //ofs.write((char *)data_d, size);
//
  num=0;
  for (int ic = 0; ic < x.size(); ic++){
      data_d[num]   = p[ic];
      num++;
  }
  size=sizeof(double)*x.size();
  ofs.write((char *)&size, sizeof(size));
  ofs.write((char *)data_d, size);
//
  //num=0;
  //for (int ic = 0; ic < numOfNode; ic++){
  //    data_d[num]   = wall_share_stress_u(ic);
  //    data_d[num+1] = wall_share_stress_v(ic);
  //    data_d[num+2] = wall_share_stress_w(ic);
  //    num=num+3;
  //}
  //size=sizeof(double)*numOfNode*3;
  //ofs.write((char *)&size, sizeof(size));
  //ofs.write((char *)data_d, size);

  delete data_d;

  ofs.close();

  if ((fp = fopen(file.c_str(), "a")) == NULL)
  {
    cout << file << " open error" << endl;
    exit(1);
  }
  fprintf(fp, "\n</AppendedData>\n");
  fprintf(fp, "</VTKFile>\n");
  fclose(fp);
}

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

int main()
{
  vector<vector<int>> element;
  vector<vector<double>> node;
  double dt = 0.0001;
  double diffusion_coefficient = 0.1;
  set_field(node, element, 30, 0.1);
  string filename = "test_grid.vtu";
  vector<vector<double>> gauss(4, vector<double>(2));
  gauss[0][0] = -0.577350296189626; gauss[0][1] = -0.577350296189626;
  gauss[1][0] = -0.577350296189626; gauss[1][1] = 0.577350296189626; 
  gauss[2][0] = 0.577350296189626; gauss[2][1] = -0.577350296189626; 
  gauss[3][0] = 0.577350296189626; gauss[3][1] = 0.577350296189626; 
  vector<vector<double>> D(node.size(), vector<double>(node.size()));
  vector<vector<double>> mass(node.size(), vector<double>(node.size()));

  vector<double> mass_centralization(node.size());
  
  for(int i=0; i<element.size(); i++){
    vector<double> N(4); 
    vector<vector<double>> dNdr(4, vector<double>(2));
    double volume=0.0;
    vector<vector<double>> element_D(4, vector<double>(4,0.0));
    vector<vector<double>> element_mass(4, vector<double>(4,0.0));
    vector<double> element_mass_centralization(4,0.0);
    for(int j=0; j<4; j++){
      C2D4_N(N,gauss[j][0],gauss[j][1]);
      C2D4_dNdr(dNdr,gauss[j][0],gauss[j][1]);
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
          D[element[i][k]][element[i][l]] += diffusion_coefficient * element_D[k][l] * detJ;
          mass[element[i][k]][element[i][l]] += element_mass[k][l] * detJ;
        }
      }
    }
  } 

  for(int i=0; i<node.size(); i++){
    mass_centralization[i] = 0.0;
    for(int j=0; j<node.size(); j++){
      mass_centralization[i] += mass[i][j];
    }
  }

  vector<double> C(node.size(),0.0);
  for(int i=0; i<31; i++){
    C[i] = 1.0;
  }

  for(int ic=0; ic<10000; ic++){
    cout << ic << endl;
    vector<double> DC(node.size(),0.0);
    for(int i=0; i<node.size(); i++){
      for(int j=0; j<node.size(); j++){
        DC[i] += D[i][j] * C[j];
      }
    }

    vector<double> MDC(node.size(),0.0);
    for(int i=0; i<node.size(); i++){
      MDC[i] = 1.0/mass_centralization[i]*DC[i];
    }

    for(int i=0; i<node.size(); i++){
      C[i] = C[i] - dt * MDC[i];
    } 
    
    for(int i=0; i<31; i++) C[i] = 1.0;
    
    if(ic%100==0){
      string filename = "test/test_" + to_string(ic/100) + ".vtu";
      export_vtu(filename, node, element,C);
    }
  }
}
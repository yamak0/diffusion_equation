#include<iostream>
#include<fstream>
#include<vector>
#include<string>
#include<sstream>
#include<cmath>
using namespace std;

void export_vtu(const std::string &file, vector<vector<int>> element, vector<vector<double>> node, vector<double> sdf)
{
    FILE *fp;
  if ((fp = fopen(file.c_str(), "w")) == NULL)
  {
    cout << file << " open error" << endl;
    exit(1);
  }

  fprintf(fp, "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\">\n");
  fprintf(fp, "<UnstructuredGrid>\n");
  fprintf(fp, "<Piece NumberOfPoints= \"%d\" NumberOfCells= \"%d\" >\n", node.size(), element.size());
  fprintf(fp, "<Points>\n");
  int offset = 0;
  fprintf(fp, "<DataArray type=\"Float64\" Name=\"Position\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%d\"/>\n",offset);
  offset += sizeof(int) + sizeof(double) * node.size() * 3;
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
    
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "</Cells>\n");

  fprintf(fp, "<PointData>\n");
  //fprintf(fp, "<DataArray type=\"Float64\" Name=\"pressure[Pa]\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%d\"/>\n",offset);
  //offset += sizeof(int) + sizeof(double) * node.size();
 
  fprintf(fp, "</PointData>\n");

  fprintf(fp, "<CellData>\n");
  
  fprintf(fp, "<DataArray type=\"Float64\" Name=\"pressure[Pa]\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%d\"/>\n",offset);
  offset += sizeof(int) + sizeof(double) * element.size();

  fprintf(fp, "</CellData>\n");
  fprintf(fp, "</Piece>\n");
  fprintf(fp, "</UnstructuredGrid>\n");
  fprintf(fp, "<AppendedData encoding=\"raw\">\n");
  fprintf(fp, "_");
  fclose(fp);

  fstream ofs;
  ofs.open(file.c_str(), ios::out | ios::app | ios_base::binary);
  double *data_d = new double[node.size()*3];
  num = 0;
  int size=0;
  for (int ic = 0; ic < node.size(); ic++){
    data_d[num] = node[ic][0];
    num++;
    data_d[num] = node[ic][1];
    num++;
    data_d[num] = 0.0;
    num++;
  }
  size=sizeof(double)*node.size()*3;
  ofs.write((char *)&size, sizeof(size));
  ofs.write((char *)data_d, size);

  //num=0;
  //for (int ic = 0; ic < node.size(); ic++){
  //    data_d[num]   = sdf[ic];
  //    num++;
  //}
  //size=sizeof(double)*node.size();
  //ofs.write((char *)&size, sizeof(size));
  //ofs.write((char *)data_d, size);

  num=0;
  for (int ic = 0; ic < element.size(); ic++){
      data_d[num]   = sdf[ic];
      num++;
  }
  size=sizeof(double)*element.size();
  ofs.write((char *)&size, sizeof(size));
  ofs.write((char *)data_d, size);


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

int CountNumbersOfTextLines(string &filePath)
{
  long i = 0;
  ifstream ifs( filePath );
  if(!ifs){
    cout << "can't file open!" << endl;
    exit(1);
  }
  if( ifs ){
    string line;
    while( true ){
      getline( ifs, line );
      i++;
      if( ifs.eof() )
        break;
    }
  }
  return i-1;
}

int main()
{
    ifstream ifs("test.csv");
    string file = "test.csv";
    string str,tmp;
    vector<vector<double>> x;
    vector<double> sdf;
    int line = CountNumbersOfTextLines(file);
    for(int i=0; i<line; i++){
        getline(ifs,str);
        if(i==0) continue;
        istringstream stream(str);
        vector<double> tmp_x;
        for(int j=0;j<4;j++){
            getline(stream,tmp,',');
            if(j==0 || j==1 || j==2) tmp_x.push_back(stod(tmp));
            else sdf.push_back(stod(tmp)); 
        }
        x.push_back(tmp_x);
    }
    ifs.close();

    for(int i=0; i<x.size(); i++){
        x[i][2] = 0.0;
    }

    vector<vector<int>> element;
    for(int i=0; i<127; i++){
        for(int j=0; j<127; j++){
            vector<int> tmp_element;
            tmp_element.push_back(i*128+j);
            tmp_element.push_back((i+1)*128+j);
            tmp_element.push_back((i+1)*128+j+1);
            tmp_element.push_back(i*128+j+1);
            
            element.push_back(tmp_element);
            for(int k=0; k<tmp_element.size(); k++){
                cout << tmp_element[k] << " ";
            }
            //if(j==10) exit(1);
            cout << endl;
        }
    }

    ifs.open("O17_cell_value.csv");
    string sdf_file = "O17_cell_value.csv";
    int sdf_line = CountNumbersOfTextLines(sdf_file);
    cout << sdf_line << endl;
    vector<double> cell_sdf;
    for(int i=0; i<sdf_line; i++){
        getline(ifs,str);
        if(i==0) continue;
        cell_sdf.push_back(stod(str));
    }
    for(int i=0; i<cell_sdf.size(); i++){
      if(cell_sdf[i]<0.0) cell_sdf[i] = 0.0;
    }

    double muximum = 0.0;
    for(int i=0; i<cell_sdf.size(); i++){
      muximum = max(muximum, cell_sdf[i]);
    }

    for(int i=0; i<cell_sdf.size(); i++){
      cell_sdf[i] = cell_sdf[i] / muximum;
    }
    //double minimum = 10000;
    //for(int i=0; i<cell_sdf.size(); i++){
    //    if(cell_sdf[i]<minimum){
    //        minimum=cell_sdf[i];
    //    }
    //}

    //for(int i=0; i<cell_sdf.size(); i++){
    //    cell_sdf[i] += fabs(minimum);
    //    if(cell_sdf[i]>32) cell_sdf[i] = -100;
    //}

    //double muximum = 0.0;

    //for(int i=0; i<cell_sdf.size(); i++){
    //    if(cell_sdf[i]>muximum){
    //        muximum = cell_sdf[i];
    //    }
    //}

    //for(int i=0; i<cell_sdf.size(); i++){
    //    cell_sdf[i] = cell_sdf[i] /muximum;
    //}
    ofstream ofs("node.dat");
    for(int i=0; i<x.size(); i++){
      ofs << x[i][0] << " " << x[i][1] << " " << x[i][2] << endl;
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

    ofs.open("vessel_phi.dat");
    for(int i=0; i<cell_sdf.size(); i++){
      ofs << cell_sdf[i] << endl;
    }
    ofs.close();

    string solid_node_file ="solid_phi_node.csv";
    ifs.open("solid_phi_node.csv");
    vector<int> solid_phi_node;
    int solid_node_line = CountNumbersOfTextLines(solid_node_file);
    for(int i=0; i<solid_node_line; i++){
      getline(ifs,str);
      if(i==0) continue;
      solid_phi_node.push_back(stoi(str));
    }
    ifs.close();
    vector<double> solid_phi(cell_sdf.size());
    for(int i=0; i<cell_sdf.size(); i++){
      solid_phi[i] = 1.0-cell_sdf[i];
      if(fabs(solid_phi[i]-1.0)<0.001) solid_phi[i] = 0.0;
    }
    for(int i=0; i<solid_phi_node.size(); i++){
      solid_phi[solid_phi_node[i]] = 1.0;
    }

    ofs.open("ISF_phi.dat");
    for(int i=0; i<cell_sdf.size(); i++){
      ofs << solid_phi[i] << endl;
    }
    ofs.close();

    export_vtu("17_cell_test.vtu",element,x,solid_phi);
}
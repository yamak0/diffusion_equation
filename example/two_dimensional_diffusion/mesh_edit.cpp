#include<iostream>
#include<fstream>
#include<vector>
#include<string>
#include<sstream>
#include<cmath>
#include<set>
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
  //out_region.csv
  //test.csv
  //O17_cell_value.csv
  //vessel_cell_data.csv
    ifstream ifs("input/test.csv");
    string file = "input/test.csv";
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
        }
    }

    ifs.open("input/O17_cell_value.csv");
    string fluid_file = "input/O17_cell_value.csv";
    int fluid_line = CountNumbersOfTextLines(fluid_file);
    vector<double> fluid_phi(element.size());
    for(int i=0; i<fluid_line; i++){
        getline(ifs,str);
        if(i==0) continue;
        fluid_phi[i-1]=stod(str);
    }
    for(int i=0; i<fluid_phi.size(); i++){
      if(fluid_phi[i]<0.0) fluid_phi[i] = 0.0;
    }

    double muximum = 0.0;
    for(int i=0; i<fluid_phi.size(); i++){
      muximum = max(muximum, fluid_phi[i]);
    }

    for(int i=0; i<fluid_phi.size(); i++){
      fluid_phi[i] = fluid_phi[i] / muximum;
    }
    
    ifs.close();

    //脳の外側の領域を決めている
    vector<int> out_region_cell;
    ifs.open("input/out_region.csv");
    string out_region_file = "input/out_region.csv";
    int out_region_line = CountNumbersOfTextLines(out_region_file);
    for(int i=0; i<out_region_line; i++){
        getline(ifs,str);
        if(i==0) continue;
        if(stod(str)<0.0001) out_region_cell.push_back(i-1);
    }
    ifs.close();

    

    vector<double> vessel_phi(element.size());
    ifs.open("input/vessel_cell_data.csv");
    string vessel_file = "input/vessel_cell_data.csv";
    int vessel_line = CountNumbersOfTextLines(vessel_file);
    for(int i=0; i<vessel_line; i++){
        getline(ifs,str);
        if(i==0) continue;
        vessel_phi[i-1]=stod(str);
    }
    ifs.close();

    vector<double> solid_phi(element.size());
    
    for(int i=0; i<element.size(); i++){
      solid_phi[i] = (1.0-vessel_phi[i])*(1.0-fluid_phi[i]);
      fluid_phi[i] = 1.0-vessel_phi[i]-solid_phi[i];
    }

    for(int i=0; i<out_region_cell.size(); i++){
      solid_phi[out_region_cell[i]] = 0.0;
      fluid_phi[out_region_cell[i]] = 0.0;
      vessel_phi[out_region_cell[i]] = 0.0;
    }

    //for(int i=0; i<element.size(); i++){
    //  if(vessel_phi[i] + solid_phi[i] + fluid_phi[i]>0.001) cout << vessel_phi[i] + solid_phi[i] + fluid_phi[i] << endl;
    //}

    //export_vtu("solid_phi.vtu",element,x,solid_phi);
    //export_vtu("fluid_phi.vtu",element,x,fluid_phi);
    //export_vtu("vessel_phi.vtu",element,x,vessel_phi);

    //output_calculation_file

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

    ofs.open("fluid_phi.dat");
    for(int i=0; i<fluid_phi.size(); i++){
      ofs << fluid_phi[i] << endl;
    }
    ofs.close();

    ofs.open("solid_phi.dat");
    for(int i=0; i<solid_phi.size(); i++){
      ofs << solid_phi[i] << endl;
    }
    ofs.close();

    ofs.open("vessel_phi.dat");
    for(int i=0; i<vessel_phi.size(); i++){
      ofs << vessel_phi[i] << endl;
    }
    ofs.close();

    ofs.open("boundary_solid.dat");
    set<int> solid_boundary_node;
    for(int i=0; i<solid_phi.size(); i++){
      if(solid_phi[i]<0.01){
        for(int j=0; j<element[i].size(); j++){
          solid_boundary_node.insert(element[i][j]);
        }
      }
    }
    for(auto itr=solid_boundary_node.begin(); itr!=solid_boundary_node.end(); itr++){
      ofs << *(itr) << " " << 0.0 << endl;
    }
    ofs.close();

    ofs.open("boundary_fluid.dat");
    set<int> fluid_boundary_node;
    for(int i=0; i<fluid_phi.size(); i++){
      if(fluid_phi[i]<0.01){
        for(int j=0; j<element[i].size(); j++){
          fluid_boundary_node.insert(element[i][j]);
        }
      }
    }
    for(auto itr=fluid_boundary_node.begin(); itr!=fluid_boundary_node.end(); itr++){
      ofs << *(itr) << " " << 0.0 << endl;
    }
    ofs.close();

    ofs.open("boundary_vessel.dat");
    set<int> vessel_boundary_zero_node;
    for(int i=0; i<vessel_phi.size(); i++){
      if(vessel_phi[i]<0.01){
        for(int j=0; j<element[i].size(); j++){
          vessel_boundary_zero_node.insert(element[i][j]);
        }
      }
    }
    set<int> vessel_boundary_source_node;
    for(int i=0; i<vessel_phi.size(); i++){
      if(vessel_phi[i]>0.0){
        for(int j=0; j<element[i].size(); j++){
          vessel_boundary_source_node.insert(element[i][j]);
        }
      }
    }
    for(auto itr=vessel_boundary_zero_node.begin(); itr!=vessel_boundary_zero_node.end(); itr++){
      ofs << *(itr) << " " << 0.0 << endl;
    }
    for(auto itr=vessel_boundary_source_node.begin(); itr!=vessel_boundary_source_node.end(); itr++){
      ofs << *(itr) << " " << 1.0 << endl;
    }
    ofs.close();

    //nodeのphi
    vector<double> vessel_node_phi(x.size());
    ifs.open("input/vessel_node_phi.csv");
    string vessel_node_file = "input/vessel_node_phi.csv";
    int vessel_node_line = CountNumbersOfTextLines(vessel_node_file);
    for(int i=0; i<vessel_node_line; i++){
        getline(ifs,str);
        if(i==0) continue;
        vessel_node_phi[i-1]=stod(str);
    }
    ifs.close();

    vector<double> fluid_node_phi(x.size());
    ifs.open("input/fluid_node_phi.csv");
    string fluid_node_file = "input/fluid_node_phi.csv";
    int fluid_node_line = CountNumbersOfTextLines(fluid_node_file);
    for(int i=0; i<fluid_node_line; i++){
        getline(ifs,str);
        if(i==0) continue;
        fluid_node_phi[i-1]=stod(str);
    }
    ifs.close();

    vector<double> solid_node_phi(x.size());
    ifs.open("input/solid_node_phi.csv");
    string solid_node_file = "input/solid_node_phi.csv";
    int solid_node_line = CountNumbersOfTextLines(solid_node_file);
    for(int i=0; i<solid_node_line; i++){
        getline(ifs,str);
        if(i==0) continue;
        solid_node_phi[i-1]=stod(str);
    }
    ifs.close();

    ofs.open("fluid_node_phi.dat");
    for(int i=0; i<fluid_node_phi.size(); i++){
      ofs << fluid_node_phi[i] << endl;
    }
    ofs.close();

    ofs.open("solid_node_phi.dat");
    for(int i=0; i<solid_node_phi.size(); i++){
      ofs << solid_node_phi[i] << endl;
    }
    ofs.close();

    ofs.open("vessel_node_phi.dat");
    for(int i=0; i<vessel_node_phi.size(); i++){
      ofs << vessel_node_phi[i] << endl;
    }
    ofs.close();
}
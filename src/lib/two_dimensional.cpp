#include"two_diffusion.hpp"
#include"shapefunction.hpp"
using namespace std;

void twodimensinal_diffusion::input_info(std::string input_file)
{
    TextParser tp;

    int ierror;
    if ((ierror = tp.read(input_file)) != TP_NO_ERROR){
        printf("\tError at reading '%s' file\n", input_file.c_str());
        exit(1);
    }

    string str,base_label,label;
    base_label = "/Geometry";
    label = base_label + "/node_file";
    if ( !tp.getInspectedValue(label,node_file)){
      cout << label << " is not set" << endl;
      exit(0);
    }

    label = base_label + "/element_file";
    if ( !tp.getInspectedValue(label,element_file)){
      cout << label << " is not set" << endl;
      exit(0);
    }

    label = base_label + "/gauss_setting";
    if ( !tp.getInspectedValue(label,gauss_setting)){
      cout << label << " is not set" << endl;
      exit(0);
    }
    if(material_judge=="F"){
        base_label = "/Fluid";
        label = base_label + "/dt";
        if ( !tp.getInspectedValue(label,dt)){
          cout << label << " is not set" << endl;
          exit(0);
        }

        label = base_label + "/diffusion_coefficient";
        if ( !tp.getInspectedValue(label,diffusion_coefficient)){
          cout << label << " is not set" << endl;
          exit(0);
        }

        label = base_label + "/time_step";
        if ( !tp.getInspectedValue(label,time)){
          cout << label << " is not set" << endl;
          exit(0);
        }
        label = base_label + "/output_interval";
        if ( !tp.getInspectedValue(label,output_interval)){
          cout << label << " is not set" << endl;
          exit(0);
        }
        label = base_label + "/outputDir";
        if ( !tp.getInspectedValue(label,outputDir)){
          cout << label << " is not set" << endl;
          exit(0);
        }

        label = base_label + "/boundary_file";
        if ( !tp.getInspectedValue(label,boundary_file)){
            cout << label << " is not set" << endl;
            exit(0);
        }
        label = base_label + "/couplint_coefficient";
        if ( !tp.getInspectedValue(label,coupling_coefficient)){
            cout << label << " is not set" << endl;
            exit(0);
        }
        label = base_label + "/phi_file";
        if ( !tp.getInspectedValue(label,phi_file)){
            cout << label << " is not set" << endl;
            exit(0);
        }
    }
    if(material_judge=="S"){
        base_label = "/Solid";
        label = base_label + "/dt";
        if ( !tp.getInspectedValue(label,dt)){
          cout << label << " is not set" << endl;
          exit(0);
        }

        label = base_label + "/diffusion_coefficient";
        if ( !tp.getInspectedValue(label,diffusion_coefficient)){
          cout << label << " is not set" << endl;
          exit(0);
        }

        label = base_label + "/time_step";
        if ( !tp.getInspectedValue(label,time)){
          cout << label << " is not set" << endl;
          exit(0);
        }
        label = base_label + "/output_interval";
        if ( !tp.getInspectedValue(label,output_interval)){
          cout << label << " is not set" << endl;
          exit(0);
        }
        label = base_label + "/outputDir";
        if ( !tp.getInspectedValue(label,outputDir)){
          cout << label << " is not set" << endl;
          exit(0);
        }
        label = base_label + "/boundary_file";
        if ( !tp.getInspectedValue(label,boundary_file)){
            cout << label << " is not set" << endl;
            exit(0);
        }
        label = base_label + "/couplint_coefficient";
        if ( !tp.getInspectedValue(label,coupling_coefficient)){
            cout << label << " is not set" << endl;
            exit(0);
        }
        label = base_label + "/phi_file";
        if ( !tp.getInspectedValue(label,phi_file)){
            cout << label << " is not set" << endl;
            exit(0);
        }
    }

    if(material_judge=="V"){
        base_label = "/Vessel";
        label = base_label + "/dt";
        if ( !tp.getInspectedValue(label,dt)){
          cout << label << " is not set" << endl;
          exit(0);
        }

        label = base_label + "/diffusion_coefficient";
        if ( !tp.getInspectedValue(label,diffusion_coefficient)){
          cout << label << " is not set" << endl;
          exit(0);
        }

        label = base_label + "/time_step";
        if ( !tp.getInspectedValue(label,time)){
          cout << label << " is not set" << endl;
          exit(0);
        }
        label = base_label + "/output_interval";
        if ( !tp.getInspectedValue(label,output_interval)){
          cout << label << " is not set" << endl;
          exit(0);
        }
        label = base_label + "/outputDir";
        if ( !tp.getInspectedValue(label,outputDir)){
          cout << label << " is not set" << endl;
          exit(0);
        }
        label = base_label + "/boundary_file";
        if ( !tp.getInspectedValue(label,boundary_file)){
            cout << label << " is not set" << endl;
            exit(0);
        }
        label = base_label + "/couplint_coefficient";
        if ( !tp.getInspectedValue(label,coupling_coefficient)){
            cout << label << " is not set" << endl;
            exit(0);
        }
        label = base_label + "/phi_file";
        if ( !tp.getInspectedValue(label,phi_file)){
            cout << label << " is not set" << endl;
            exit(0);
        }
    }

    mkdir(outputDir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
 
    read_geometry();
    input_phi();
    if(material_judge=="V") boundary_initialize();
    export_vtu("test.vtu");
}

int twodimensinal_diffusion::CountNumbersOfTextLines(string &filePath )
{
  long i = 0;
  ifstream ifs( filePath );
  if(!ifs){
    cout << "can't file open!" << endl;
    cout << filePath << endl;
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

void twodimensinal_diffusion::read_geometry()
{
    string str,tmp;
    numOfNode=CountNumbersOfTextLines(node_file);
    node.resize(numOfNode);
    C.resize(numOfNode);
    for(int i=0; i<C.size(); i++){
        C[i] = 0.0;
    }
    for(int i=0; i<numOfNode; i++){
        node[i].resize(3);
    }
    ifstream ifs(node_file);
    for(int i=0; i<numOfNode; i++){
        getline(ifs,str);
        istringstream stream(str);
        for(int j=0;j<3;j++){
            getline(stream,tmp,' ');
            node[i][j] = stof(tmp);
        }
    }
    ifs.close();

    numOfElm=CountNumbersOfTextLines(element_file);
    element.resize(numOfElm);
    for(int i=0; i<element.size(); i++){
        element[i].resize(4);
    }
    ifs.open(element_file);
    for(int i=0; i<numOfElm; i++){
        getline(ifs,str);
        istringstream stream(str);
        for(int j=0;j<4;j++){
            getline(stream,tmp,' ');
            element[i][j] = stoi(tmp);
        }
    }
    ifs.close();
}

void twodimensinal_diffusion::input_phi()
{
  string str,tmp;
  ifstream ifs(phi_file);
  phi.resize(numOfNode);
  for(int i=0; i<numOfNode; i++){
    getline(ifs,str);
    phi[i] = stod(str);
  }
}


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

void twodimensinal_diffusion::export_vtu(const std::string &file)
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

  fprintf(fp, "<DataArray type=\"Float64\" Name=\"pressure[Pa]\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%d\"/>\n",offset);
  offset += sizeof(int) + sizeof(double) * node.size();

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

  num=0;
  for (int ic = 0; ic < node.size(); ic++){
      data_d[num]   = C[ic];
      num++;
  }
  size=sizeof(double)*node.size();
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
  for(int k=0;k<2;k++){
    for(int l=0;l<2;l++){
      dxdr[k][l] = 0e0;
      for(int p=0;p<4;p++){
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
  for(int k=0; k<4; k++){
    for(int l=0; l<2; l++){
      dNdx[k][l] = 0.0;
      for(int p=0; p<2; p++){
        dNdx[k][l] += dNdr[k][p]*drdx[p][l];
      }
    }
  }
}

void twodimensinal_diffusion::calc_matrix()
{
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
                element_mass_centralization[k] += element_mass[k][l];
              }
            }

            for(int k=0; k<4; k++){
              for(int l=0; l<4; l++){
                D[element[i][k]][element[i][l]] += diffusion_coefficient * element_D[k][l] * detJ;
                mass_centralization[element[i][k]] += element_mass_centralization[k] * detJ;
              }
            }
        }
    }
}

void twodimensinal_diffusion::boundary_setting()
{
    for(int i=0; i<numOfBoundaryNode; i++){
        C[boundary_node[i]] = boundary_value[i];
    }
}


void twodimensinal_diffusion::time_step(vector<double> diff)
{
    vector<double> R(numOfNode);
    vector<double> DC(numOfNode,0.0);
    vector<double> DcR(numOfNode, 0.0);
    vector<double> MDcR(numOfNode, 0.0);
    for(int i=0; i<node.size(); i++){
        for(int j=0; j<node.size(); j++){
            DC[i] += D[i][j] * (phi[j]*C[j]);
        }
    }

    for(int i=0; i<numOfNode; i++){
        R[i] = mass_centralization[i]*diff[i];
        DcR[i] = DC[i]-R[i];
    }

    
    vector<double> MDC(node.size(),0.0);
    for(int i=0; i<node.size(); i++){
        MDcR[i] = 1.0/mass_centralization[i]*DcR[i];
    }
    
    for(int i=0; i<node.size(); i++){
        C[i] = C[i] - dt * MDcR[i];
    } 
    
    if(material_judge=="V") boundary_setting();
}

void twodimensinal_diffusion::dump(int ic)
{
    if(material_judge=="S"){
        string filename = outputDir + "/solid_" + to_string(ic) + ".vtu";
        export_vtu(filename);
    }
    if(material_judge=="F"){
        string filename = outputDir + "/fluid_" + to_string(ic) + ".vtu";
        export_vtu(filename);
    }
    if(material_judge=="V"){
        string filename = outputDir + "/vessel_" + to_string(ic) + ".vtu";
        export_vtu(filename);
    }
}

double twodimensinal_diffusion::access_c(int ic)
{
    return C[ic];
}

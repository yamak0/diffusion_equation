#include"two_diffusion.hpp"
#include"shapefunction.hpp"
using namespace std;

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
  ifs.close();

  ifs.open(phi_visualize_file);
  phi_v.resize(numOfNode);
  for(int i=0; i<numOfNode; i++){
    getline(ifs,str);
    phi_v[i] = stod(str);
  }
  ifs.close();
}

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
        label = base_label + "/couplint_coefficient_v";
        if ( !tp.getInspectedValue(label,coupling_coefficient_v)){
            cout << label << " is not set" << endl;
            exit(0);
        }
        label = base_label + "/couplint_coefficient_s";
        if ( !tp.getInspectedValue(label,coupling_coefficient_s)){
            cout << label << " is not set" << endl;
            exit(0);
        }
        label = base_label + "/phi_file";
        if ( !tp.getInspectedValue(label,phi_file)){
            cout << label << " is not set" << endl;
            exit(0);
        }
        label = base_label + "/phi_visualize_file";
        if ( !tp.getInspectedValue(label,phi_visualize_file)){
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
        label = base_label + "/phi_visualize_file";
        if ( !tp.getInspectedValue(label,phi_visualize_file)){
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
        label = base_label + "/phi_visualize_file";
        if ( !tp.getInspectedValue(label,phi_visualize_file)){
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
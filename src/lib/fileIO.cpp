#include"two_diffusion.hpp"
#include"shapefunction.hpp"
using namespace H5;
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
    if(!ifs){
      cout << "can't open " + node_file << endl;
    }
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
    if(!ifs){
      cout << "can't open " << element_file << endl;
    }
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
  if(!ifs){
    cout << "can't open " << phi_file << endl;
  }
  phi.resize(numOfElm);
  for(int i=0; i<numOfElm; i++){
    getline(ifs,str);
    phi[i] = stod(str);
  }
  ifs.close();

  if(material_judge=="F") export_vtu("fluid_phi.vtu","CELL",phi);
  if(material_judge=="S") export_vtu("solid_phi.vtu","CELL",phi);
  if(material_judge=="V") export_vtu("vessel_phi.vtu","CELL",phi);
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
    label = base_label + "/omp_num_threads";
    if ( !tp.getInspectedValue(label,numofomp)){
      cout << label << " is not set" << endl;
      exit(0);
    }
    omp_set_num_threads(numofomp);

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
        label = base_label + "/couplint_coefficient_vc";
        if ( !tp.getInspectedValue(label,coupling_coefficient_vc)){
            cout << label << " is not set" << endl;
            exit(0);
        }
        label = base_label + "/couplint_coefficient_ci";
        if ( !tp.getInspectedValue(label,coupling_coefficient_ci)){
            cout << label << " is not set" << endl;
            exit(0);
        }
        label = base_label + "/phi_file";
        if ( !tp.getInspectedValue(label,phi_file)){
            cout << label << " is not set" << endl;
            exit(0);
        }
        label = base_label + "/output_h5_name";
        if ( !tp.getInspectedValue(label,output_h5_name)){
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
        label = base_label + "/couplint_coefficient_vi";
        if ( !tp.getInspectedValue(label,coupling_coefficient_vi)){
            cout << label << " is not set" << endl;
            exit(0);
        }
        label = base_label + "/phi_file";
        if ( !tp.getInspectedValue(label,phi_file)){
            cout << label << " is not set" << endl;
            exit(0);
        }
        label = base_label + "/output_h5_name";
        if ( !tp.getInspectedValue(label,output_h5_name)){
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
        label = base_label + "/phi_file";
        if ( !tp.getInspectedValue(label,phi_file)){
            cout << label << " is not set" << endl;
            exit(0);
        }
      
        label = base_label + "/output_h5_name";
        if ( !tp.getInspectedValue(label,output_h5_name)){
            cout << label << " is not set" << endl;
            exit(0);
        }
    }

    mkdir(outputDir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
    cout << "Read geometry" << endl;
    read_geometry();
    cout << "Read phi" << endl;
    input_phi();
    cout << "boundary initialize" << endl;
    boundary_initialize();
}

void twodimensinal_diffusion::export_vtu(const std::string &file, std::string judge, vector<double> output_value)
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
  if(judge=="point"){
    fprintf(fp, "<DataArray type=\"Float64\" Name=\"pressure[Pa]\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%d\"/>\n",offset);
    offset += sizeof(int) + sizeof(double) * node.size();
  }
  fprintf(fp, "</PointData>\n");

  fprintf(fp, "<CellData>\n");
  if(judge=="CELL"){
    fprintf(fp, "<DataArray type=\"Float64\" Name=\"pressure[Pa]\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%d\"/>\n",offset);
    offset += sizeof(int) + sizeof(double) * element.size();
  }
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
  if(judge=="point"){
    for (int ic = 0; ic < node.size(); ic++){
        data_d[num]   = output_value[ic];
        num++;
    }
    size=sizeof(double)*node.size();
    ofs.write((char *)&size, sizeof(size));
    ofs.write((char *)data_d, size);
  }
  if(judge=="CELL"){
    for (int ic = 0; ic < element.size(); ic++){
        data_d[num]   = output_value[ic];
        num++;
    }
    size=sizeof(double)*element.size();
    ofs.write((char *)&size, sizeof(size));
    ofs.write((char *)data_d, size);
  }
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

void twodimensinal_diffusion::transform_point_data_to_cell_data_phi(vector<double> &phi_C, vector<double> C)
{
  for(int i=0; i<numOfElm; i++){
    double tmp_C=0.0;
    for(int j=0; j<element[i].size(); j++){
      tmp_C+=C[element[i][j]];
    }
    tmp_C/=element[i].size();
    phi_C[i]=tmp_C*phi[i];
  }
}

void twodimensinal_diffusion::dump(int ic)
{
    if(material_judge=="S"){
      vector<double> phiC(numOfElm);
      transform_point_data_to_cell_data_phi(phiC, C);
      string filename = outputDir + "/solid_" + to_string(ic) + ".vtu";
      export_vtu(filename, "CELL", phiC);
    }
    if(material_judge=="F"){
      vector<double> phiC(numOfElm);
      transform_point_data_to_cell_data_phi(phiC, C);
      string filename = outputDir + "/fluid_" + to_string(ic) + ".vtu";
      export_vtu(filename, "CELL", phiC);
    }
    if(material_judge=="V"){
      vector<double> phiC(numOfElm);
      transform_point_data_to_cell_data_phi(phiC, C);
      string filename = outputDir + "/vessel_" + to_string(ic) + ".vtu";
      export_vtu(filename, "CELL", phiC);
    }
}

void twodimensinal_diffusion::exportHDF5_double_1D(H5::H5File &file, const std::string &dataName, vector<double> i_data, int i_dim) 
{
  H5std_string DATASET_NAME(dataName.c_str());

  hsize_t dim[1] = {i_dim}; // dataset dimensions
  H5::DataSpace dataspace(1, dim);

  double *data;
  data = new double[i_data.size()];

  for (int i = 0; i < i_data.size(); i++) {
    data[i] = i_data[i];
  }

  H5::IntType datatype(H5::PredType::NATIVE_DOUBLE);
  datatype.setOrder(H5T_ORDER_LE);
  H5::DataSet dataset = file.createDataSet(DATASET_NAME, datatype, dataspace);
  dataset.write(&data[0], H5::PredType::NATIVE_DOUBLE);
  delete[] data;
}

void twodimensinal_diffusion::hdf5_dump(int ic)
{
  H5std_string FILE_NAME(output_h5_name.c_str());
  if (ic == 0) {
    H5File file(FILE_NAME, H5F_ACC_TRUNC);
    std::string dataName;
    std::string Gr = "/"+to_string(ic);
    file.createGroup(Gr.c_str());
    Group group = file.openGroup(Gr.c_str());
    dataName = Gr + "/O17_concentration";
    vector<double> phiC(numOfElm);
    transform_point_data_to_cell_data(phiC, C);
    exportHDF5_double_1D(file, dataName, phiC, phiC.size());
  }
  else {
    H5File file(FILE_NAME, H5F_ACC_RDWR);
    std::string dataName;
    std::string Gr = "/" + to_string(ic);
    file.createGroup(Gr.c_str());
    Group group = file.openGroup(Gr.c_str());
    dataName = Gr + "/O17_concentration";
    vector<double> phiC(numOfElm);
    transform_point_data_to_cell_data(phiC, C);
    exportHDF5_double_1D(file, dataName, phiC, phiC.size());
  }
}
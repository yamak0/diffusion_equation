#include"two_diffusion.hpp"
#include"shapefunction.hpp"
#include <sys/stat.h>
#include <time.h>

using namespace H5;
using namespace std;

void export_vtu(const std::string &file, vector<vector<int>> element, vector<vector<double>> node, vector<double> C)
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

  num=0;
  for (int ic = 0; ic < element.size(); ic++){
      data_d[num]   = C[ic];
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
void exportHDF5_double_1D(H5::H5File &file, const std::string &dataName, vector<double> i_data, int i_dim) 
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

void hdf5_dump(string output_h5_name, int ic, vector<double> C_sum)
{
  H5std_string FILE_NAME(output_h5_name.c_str());
  if (ic == 0) {
    H5File file(FILE_NAME, H5F_ACC_TRUNC);
    std::string dataName;
    std::string Gr = "/"+to_string(ic);
    file.createGroup(Gr.c_str());
    Group group = file.openGroup(Gr.c_str());
    dataName = Gr + "/O17_concentration";
    exportHDF5_double_1D(file, dataName, C_sum, C_sum.size());
  }
  else {
    H5File file(FILE_NAME, H5F_ACC_RDWR);
    std::string dataName;
    std::string Gr = "/" + to_string(ic);
    file.createGroup(Gr.c_str());
    Group group = file.openGroup(Gr.c_str());
    dataName = Gr + "/O17_concentration";
    exportHDF5_double_1D(file, dataName, C_sum, C_sum.size());
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
  twodimensinal_diffusion Solid("S"), Vessel("V"), Fluid("F");
  vector<double> C_sum;
  cout << "input infomation" << endl;
  Fluid.input_info(input_file);
  Solid.input_info(input_file);
  Vessel.input_info(input_file);
  
  cout << "gauss_point_setting" << endl;
  Fluid.gauss_point_setting();
  Solid.gauss_point_setting();
  Vessel.gauss_point_setting();

  cout << "matrix initialize" << endl;
  Fluid.matrix_initialize();
  Solid.matrix_initialize();
  Vessel.matrix_initialize();

  cout << "calc matrix" << endl;
  Fluid.calc_matrix();
  Solid.calc_matrix();
  Vessel.calc_matrix();

  cout << "set boundary" << endl;
  vector<double> Q_vc(Vessel.numOfElm), Q_cv(Vessel.numOfElm), Q_vi(Vessel.numOfElm), Q_iv(Vessel.numOfElm), Q_ci(Vessel.numOfElm), Q_ic(Vessel.numOfElm);
  Vessel.boundary_setting(0.0, Q_cv, Q_iv);
  C_sum.resize(Vessel.numOfElm);

  cout << "main loop" << endl;
  for(int i=0; i<Vessel.time; i++){
    vector<double> element_C_vessel(Vessel.numOfElm), element_C_Fluid(Fluid.numOfElm), element_C_Solid(Solid.numOfElm);
    Vessel.transform_point_data_to_cell_data(element_C_vessel, Vessel.C);
    Fluid.transform_point_data_to_cell_data(element_C_Fluid, Fluid.C);
    Solid.transform_point_data_to_cell_data(element_C_Solid, Solid.C);
    for(int j=0; j<Vessel.numOfElm; j++){
      Q_vc[j] = Fluid.coupling_coefficient_vc*sqrt(Vessel.phi[j])*(element_C_vessel[j]-element_C_Fluid[j]);
      Q_vi[j] = Solid.coupling_coefficient_vi*sqrt(Vessel.phi[j])*(element_C_vessel[j]-element_C_Solid[j]);
      Q_cv[j] = Fluid.coupling_coefficient_vc*sqrt(Vessel.phi[j])*(element_C_Fluid[j]-element_C_vessel[j]);
      Q_ci[j] = Fluid.coupling_coefficient_ci*sqrt(Fluid.phi[j])*(element_C_Fluid[j]-element_C_Solid[j]);
      Q_iv[j] = Solid.coupling_coefficient_vi*sqrt(Vessel.phi[j])*(element_C_Solid[j]-element_C_vessel[j]);
      Q_ic[j] = Fluid.coupling_coefficient_ci*sqrt(Fluid.phi[j])*(element_C_Solid[j]-element_C_Fluid[j]);
    }
    Vessel.boundary_setting(Vessel.dt*i, Q_vc, Q_vi);
    Fluid.time_step(Q_cv, Q_ci, Vessel.dt*i);
    Solid.time_step(Q_iv, Q_ic, Vessel.dt*i);
    
    if(i%Fluid.output_interval==0){
      Fluid.dump(i/Fluid.output_interval);
      Fluid.hdf5_dump(i/Fluid.output_interval);
    }
    if(i%Solid.output_interval==0){
      Solid.dump(i/Solid.output_interval);
      Solid.hdf5_dump(i/Solid.output_interval);
    }
    if(i%Vessel.output_interval==0){
      Vessel.dump(i/Vessel.output_interval);
      Vessel.hdf5_dump(i/Vessel.output_interval);
    }
    if(i%Vessel.output_interval==0){
      vector<double> Vessel_phiC(Vessel.numOfElm),CSF_phiC(Fluid.numOfElm),ISF_phiC(Solid.numOfElm);
      Vessel.transform_point_data_to_cell_data_phi(Vessel_phiC, Vessel.C);
      Fluid.transform_point_data_to_cell_data_phi(CSF_phiC, Fluid.C);
      Solid.transform_point_data_to_cell_data_phi(ISF_phiC, Solid.C);
      for(int j=0; j<Vessel.numOfElm; j++){
        C_sum[j] = Vessel_phiC[j] + CSF_phiC[j] + ISF_phiC[j];
      }
      string dir = "out_C";
      mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
      string filename = dir + "/test_" + to_string(i/Vessel.output_interval) + ".vtu";
      cout << filename << endl;
      hdf5_dump("sum_O17.h5",i/Vessel.output_interval,C_sum);
      export_vtu(filename,Vessel.element,Vessel.node,C_sum);
    }
  }
  vector<double> evaluation_phi(Vessel.numOfElm);
  ifstream ifs("evaluation_C.dat");
  if(!ifs){
    cout << "can't open evaluation_C.dat" << endl;
  }
  string str;
  for(int i=0; i<Vessel.numOfElm; i++){
    getline(ifs,str);
    evaluation_phi[i] = stod(str);
  }
  ifs.close();
  double evaluation=0.0;
  for(int i=0; i<Vessel.numOfElm; i++){
    evaluation += sqrt(pow(evaluation_phi[i]-C_sum[i],2.0));
  }
  cout << evaluation << endl;
}
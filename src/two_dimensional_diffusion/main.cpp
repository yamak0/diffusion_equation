#include"two_diffusion.hpp"
#include"shapefunction.hpp"
#include <sys/stat.h>

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
  fprintf(fp, "<DataArray type=\"Float64\" Name=\"pressure[Pa]\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%d\"/>\n",offset);
  offset += sizeof(int) + sizeof(double) * element.size();

  

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

  Fluid.matrix_initialize();
  cout << "matrix initialize" << endl;
  Solid.matrix_initialize();
  Vessel.matrix_initialize();

  Fluid.calc_matrix();
  cout << "calc matrix" << endl;
  Solid.calc_matrix();
  //Vessel.calc_matrix();
  cout << "set boundary" << endl;
  Vessel.boundary_setting();
  C_sum.resize(Vessel.numOfNode);
  int j,ic;

  for(int i=0; i<Vessel.boundary_node.size(); i++){
    Vessel.phi_node[Vessel.boundary_node[i]]=0.7;
    Fluid.phi_node[Vessel.boundary_node[i]]=0.2;
    Solid.phi_node[Vessel.boundary_node[i]]=0.1;
  }
  
  //export_vtu("sdf_test_v.vtu",Vessel.element,Vessel.node,element_C);
  //exit(1);
  cout << "main loop" << endl;
  for(int i=0; i<Vessel.time; i++){
    cout << i << endl;
    vector<double> fluid_vessel_diff(Vessel.numOfNode,0.0), vessel_fluid_diff(Vessel.numOfNode); 
    vector<double> fluid_solid_diff(Solid.numOfNode), solid_fluid_diff(Solid.numOfNode,0.0); 
    //vector<double> fluid_source(Fluid.numOfNode);
    //#pragma omp parallel for
    for(j=0; j<Vessel.numOfNode; j++){
      //fluid_vessel_diff[j]=-Vessel.coupling_coefficient*(Vessel.phi_node[j]*Vessel.access_c(j)-Fluid.phi_node[j]*Solid.access_c(j));
      vessel_fluid_diff[j]=-Fluid.coupling_coefficient_v*Fluid.phi_node[j]*(Fluid.access_c(j)-Vessel.access_c(j));
      fluid_solid_diff[j]=-Solid.coupling_coefficient*Solid.phi_node[j]*(Solid.access_c(j)-Fluid.access_c(j));
      //solid_fluid_diff[j]=-Fluid.coupling_coefficient_s*(Fluid.phi[j]*Fluid.access_c(j)-Solid.phi[j]*Solid.access_c(j));
      //fluid_source[j] = vessel_fluid_diff[j] + solid_fluid_diff[j]; 
    }
    //Vessel.time_step(fluid_vessel_diff);
    Fluid.time_step(vessel_fluid_diff);
    Solid.time_step(fluid_solid_diff);

    if(i%Fluid.output_interval==0){
      Fluid.dump(i/Fluid.output_interval);
    }
    if(i%Solid.output_interval==0){
      Solid.dump(i/Solid.output_interval);
    }
    if(i%Vessel.output_interval==0){
      Vessel.dump(i/Vessel.output_interval);
    }
    if(i%Vessel.output_interval==0){
    //  #pragma omp parallel for
      for(ic=0; ic<C_sum.size(); ic++){
        C_sum[ic] = Solid.phi_node[ic]*Solid.access_c(ic)+Fluid.phi_node[ic]*Fluid.access_c(ic)+Vessel.phi_node[ic]*Vessel.access_c(ic);
      }
      string dir = "out_C";
      mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
      string filename = dir + "/test_" + to_string(i/Vessel.output_interval) + ".vtu";
      export_vtu(filename,Vessel.element,Vessel.node,C_sum);
    }
  }
}
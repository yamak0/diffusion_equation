#include"two_diffusion.hpp"
#include"shapefunction.hpp"

using namespace std;

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

    ofstream ofs("node.dat");
    for(int i=0; i<node.size(); i++){
      ofs << node[i][0] << " " << node[i][1] << " " << node[i][2] << endl;
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

    ofs.open("boundary.dat");
    for(int i=0; i<node.size(); i++){
      if(fabs(node[i][1]-0.0)<0.0001){
        ofs << i << " " << 1.0 << endl;
      }
    }
    ofs.close();
}

int main(int argc,char *argv[])
{
  //input argument
  if(argc!=2){
    printf("Invalid input. Please set tp file\n");
    return 1;
  }
  std::string input_file = argv[1];
  twodimensinal_diffusion Fluid("F"), Solid("S"), Vessel("V");

  Fluid.input_info(input_file);
  Solid.input_info(input_file);
  Vessel.input_info(input_file);

  Fluid.gauss_point_setting();
  Solid.gauss_point_setting();
  Vessel.gauss_point_setting();

  Fluid.matrix_initialize();
  Solid.matrix_initialize();
  Vessel.matrix_initialize();

  Fluid.calc_matrix();
  Solid.calc_matrix();
  Vessel.calc_matrix();

  Vessel.boundary_setting();

  for(int i=0; i<Fluid.time; i++){
    cout << i << endl;
    vector<double> fluid_vessel_diff(Fluid.numOfNode), vessel_fluid_diff(Fluid.numOfNode); 
    vector<double> fluid_solid_diff(Solid.numOfNode), solid_fluid_diff(Solid.numOfNode); 
    vector<double> fluid_source(Fluid.numOfNode);
    
    for(int j=0; j<Fluid.numOfNode; j++){
      fluid_vessel_diff[j]=-Vessel.coupling_coefficient*(Vessel.access_c(j)-Fluid.access_c(j));
      vessel_fluid_diff[j]=-Fluid.coupling_coefficient*(Fluid.access_c(j)-Vessel.access_c(j));
      fluid_solid_diff[j]=-Solid.coupling_coefficient*(Solid.access_c(j)-Fluid.access_c(j));
      solid_fluid_diff[j]=-Fluid.coupling_coefficient*(Fluid.access_c(j)-Solid.access_c(j));
      fluid_source[j] = vessel_fluid_diff[j] + solid_fluid_diff[j]; 
    }
    Vessel.time_step(fluid_vessel_diff);
    Fluid.time_step(fluid_source);
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
  }
}
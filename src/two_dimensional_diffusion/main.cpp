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
  twodimensinal_diffusion Fluid;

  Fluid.input_info(input_file);
  Fluid.gauss_point_setting();
  Fluid.matrix_initialize();
  Fluid.calc_matrix();
  Fluid.time_step();
}
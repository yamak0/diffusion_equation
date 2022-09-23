#include"one_diffusion.hpp"
#include<cmath>
#include<vector>
using namespace std;

int main(int argc,char *argv[])
{
  cout << "check" << endl;
    //input argument
    if(argc!=2){
      printf("Invalid input. Please set tp file\n");
      return 1;
    }

    std::string input_file = argv[1];

    TextParser tp;
    int div, itr, output_t;
    double coupling_f, coupling_s;
    double phi_s, phi_f;
    string outputDir;
    
    cout << "create problem" << endl;
    onedimensinal_diffusion Fluid("F"), Solid("S");
    Fluid.element_phi.resize(10);
    Solid.element_phi.resize(10);
    Fluid.element_phi[0] = 0.9; Solid.element_phi[0] = 0.1;
    Fluid.element_phi[1] = 0.8; Solid.element_phi[1] = 0.2; 
    Fluid.element_phi[2] = 0.7; Solid.element_phi[2] = 0.3; 
    Fluid.element_phi[3] = 0.6; Solid.element_phi[3] = 0.4; 
    Fluid.element_phi[4] = 0.5; Solid.element_phi[4] = 0.5; 
    Fluid.element_phi[5] = 0.4; Solid.element_phi[5] = 0.6; 
    Fluid.element_phi[6] = 0.3; Solid.element_phi[6] = 0.7; 
    Fluid.element_phi[7] = 0.2; Solid.element_phi[7] = 0.8; 
    Fluid.element_phi[8] = 0.1; Solid.element_phi[8] = 0.9; 
    Fluid.element_phi[9] = 0.1; Solid.element_phi[9] = 0.9; 

    vector<double> fluid_node_phi(11);
    vector<double> solid_node_phi(11);

    for(int i=1; i<=9; i++){
      fluid_node_phi[i] = (Fluid.element_phi[i-1] + Fluid.element_phi[i])/2.0;
      solid_node_phi[i] = (Solid.element_phi[i-1] + Solid.element_phi[i])/2.0;
    }
    fluid_node_phi[0] =0.9; fluid_node_phi[10] = 0.1;
    fluid_node_phi[0] =0.1; solid_node_phi[10] = 0.9;

    int ierror;
    if ((ierror = tp.read(input_file)) != TP_NO_ERROR) {
     printf("\tError at reading '%s' file\n", input_file.c_str());
      return 1;
    }
    string str,base_label,label,disp;
    base_label = "/Coupling";
    label = base_label + "/numOfelement";
    if ( !tp.getInspectedValue(label,div)){
      cout << label << " is not set" << endl;
      exit(0);
    }
    label = base_label + "/iteration";
    if ( !tp.getInspectedValue(label,itr)){
      cout << label << " is not set" << endl;
      exit(0);
    }
    label = base_label + "/coupling_f";
    if ( !tp.getInspectedValue(label,coupling_f)){
      cout << label << " is not set" << endl;
      exit(0);
    }
    label = base_label + "/coupling_s";
    if ( !tp.getInspectedValue(label,coupling_s)){
      cout << label << " is not set" << endl;
      exit(0);
    }
    label = base_label + "/output_t";
    if ( !tp.getInspectedValue(label,output_t)){
      cout << label << " is not set" << endl;
      exit(0);
    }
    label = base_label + "/outputDir";
    if ( !tp.getInspectedValue(label,outputDir)){
      cout << label << " is not set" << endl;
      exit(0);
    }
    vector<double> c(div+1);
    //vector<double> phi(div+1);
//
    //for(int i=0; i<phi.size(); i++){
    //  if(i==70) phi[i] =0.1;
    //  else if(i==69 || i==71) phi[i]=0.2;
    //  else if(i==68 || i==72) phi[i]=0.3;
    //  else phi[i] =0.9;
    //  //phi[i] = 1.0/(1.0+exp(-((double(i)-50.0)/10.0)));
    //}

    if ((ierror = Fluid.tp.read(input_file)) != TP_NO_ERROR) {
     printf("\tError at reading '%s' file\n", input_file.c_str());
      return 1;
    }
    if ((ierror = Solid.tp.read(input_file)) != TP_NO_ERROR) {
     printf("\tError at reading '%s' file\n", input_file.c_str());
      return 1;
    }
    cout << "fluid input" << endl;
    Fluid.input_parameter();
    Fluid.initialize();
    cout << "solid input" << endl;
    Solid.input_parameter();
    Solid.initialize();

    //Fluid.setting_phi(0.2);
    //Solid.setting_phi(0.8);
    cout << "fluid initialize" << endl;
    Fluid.set_initial_boundary_node();
    cout << "fluid calc matrix" << endl;
    Fluid.calc_mass_matrix();
    Fluid.calc_K_matrix();
    Solid.set_initial_boundary_node();
    cout << "solid calc matrix" << endl;
    Solid.calc_mass_matrix();
    Solid.calc_K_matrix();
    cout << "main loop" << endl;
    for(int i=0; i<itr; i++){
        cout << i << endl;
        vector<double> fluid_diff(div), solid_diff(div);
        for(int j=0; j<fluid_diff.size(); j++){
            fluid_diff[j]=-coupling_f*(fluid_node_phi[j]*Fluid.access_c(j)-solid_node_phi[j]*Solid.access_c(j));
            solid_diff[j]=-coupling_s*(solid_node_phi[j]*Solid.access_c(j)-fluid_node_phi[j]*Fluid.access_c(j));
            if(i==itr-1) cout << j << " " << fluid_diff[j] << " " << solid_diff[j] << endl;
        }
        Fluid.time_step(fluid_diff, i);
        Solid.time_step(solid_diff, i);
        for(int j=0; j<solid_node_phi.size(); j++){
            c[j] = (1.0-solid_node_phi[j])*Fluid.access_c(j) + solid_node_phi[j]*Solid.access_c(j);
        }
        if(i%output_t==0){
            Fluid.dump(i/output_t);
            Solid.dump(i/output_t);
            mkdir(outputDir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
            string filename = outputDir+"/"+to_string(i/output_t) + ".dat";
            ofstream ofs(filename);
            for(int j=0; j<c.size(); j++){
                ofs << c[j] << endl;
            }
            ofs.close();
        }
    }
}
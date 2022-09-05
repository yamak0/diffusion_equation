#include"diffusion.hpp"
#include<cmath>
using namespace std;

int main(int argc,char *argv[])
{
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
    
    onedimensinal_diffusion Fluid("F"), Solid("S");

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
    vector<double> phi(div+1);

    for(int i=0; i<phi.size(); i++){
      phi[i] = 1.0/(1.0+exp(-((double(i)-50.0)/10.0)));
    }
    
    if ((ierror = Fluid.tp.read(input_file)) != TP_NO_ERROR) {
     printf("\tError at reading '%s' file\n", input_file.c_str());
      return 1;
    }
    if ((ierror = Solid.tp.read(input_file)) != TP_NO_ERROR) {
     printf("\tError at reading '%s' file\n", input_file.c_str());
      return 1;
    }
    Fluid.input_parameter();
    Fluid.initialize();
    Solid.input_parameter();
    Solid.initialize();

    //Fluid.setting_phi(0.2);
    //Solid.setting_phi(0.8);

    Fluid.set_initial_boundary_node();
    Fluid.calc_mass_matrix();
    Fluid.calc_K_matrix();
    Solid.set_initial_boundary_node();
    Solid.calc_mass_matrix();
    Solid.calc_K_matrix();

    for(int i=0; i<itr; i++){
        cout << i << endl;
        vector<double> fluid_diff(div), solid_diff(div);
        for(int j=0; j<fluid_diff.size(); j++){
            fluid_diff[j]=-coupling_f*(Fluid.access_c(j)-Solid.access_c(j));
            solid_diff[j]=-coupling_s*(Solid.access_c(j)-Fluid.access_c(j));
            //cout << j << " " << fluid_diff[j] << " " << solid_diff[j] << endl;
        }
        Fluid.time_step(fluid_diff);
        Solid.time_step(solid_diff);
        for(int j=0; j<phi.size(); j++){
            c[j] = (1.0-phi[j])*Fluid.access_c(j) + phi[j]*Solid.access_c(j);
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
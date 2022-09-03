#include"diffusion.hpp"

using namespace std;

void calc_diff(vector<double> &fluid_diff, vector<double> &solid_diff, onedimensinal_diffusion Fluid, onedimensinal_diffusion Solid)
{
    for(int i=0; i<fluid_diff.size(); i++){
        fluid_diff[i]=-0.02*(Fluid.access_c(i)-Solid.access_c(i));
        solid_diff[i]=-0.02*(Solid.access_c(i)-Fluid.access_c(i));
        cout << fluid_diff[i] << " " << solid_diff[i] << endl;
    }
    cout << endl;
}

int main()
{
    int div =100;
    double c(div);
    double phi_s, phi_f;
    phi_s = 0.8;
    phi_f = 0.2;
    onedimensinal_diffusion Fluid(div, 0.1, 0.1, 0.4), Solid(div, 0.1, 0.01, 0.01);
    Fluid.setting_phi(0.2);
    Solid.setting_phi(0.8);

    Fluid.set_initial_boundary_node(0,1.0);
    Fluid.calc_mass_matrix();
    Fluid.calc_K_matrix();
    Solid.set_initial_boundary_node(50,1.0);
    Solid.calc_mass_matrix();
    Solid.calc_K_matrix();
    for(int i=0; i<1000; i++){
        vector<double> Fluid_diff(div), Solid_diff(div);
        calc_diff(Fluid_diff, Solid_diff, Fluid, Solid);
        //cout << "fluid itr: " << i << endl;
        //cout << "solid itr: " << i << endl;
        Fluid.time_step(Fluid_diff);
        Solid.time_step(Solid_diff);
        if(i%10==0){
            Fluid.dump(i/10,"fluid");
            Solid.dump(i/10,"solid");
            //cout << Fluid.return_sum() << " " << Solid.return_sum() << endl;
        }
        //exit(1);
    }
    //exit(1);
}
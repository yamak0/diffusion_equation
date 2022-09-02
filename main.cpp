#include"diffusion.hpp"

using namespace std;

int main()
{
    onedimensinal_diffusion Fluid(100, 0.1, 0.01, 0.1), Solid(100, 0.1, 0.01, 0.1);
    Fluid.calc_mass_matrix();
    Fluid.calc_K_matrix();
    Fluid.set_boundary_node(0,1.0);
    for(int i=0; i<10000; i++){
        cout << "fluid itr: " << i << endl;
        Fluid.time_step();
        if(i%100==0) Fluid.dump(i/100,"fluid");
    }

    Solid.calc_mass_matrix();
    Solid.calc_K_matrix();
    Solid.set_boundary_node(50,1.0);
    for(int i=0; i<10000; i++){
        cout << "solid itr: " << i << endl;
        Solid.time_step();
        if(i%100==0) Solid.dump(i/100,"solid");
    }
}
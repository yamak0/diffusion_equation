#include"diffusion.hpp"

using namespace std;

int main()
{
    onedimensinal_diffusion Fluid(10, 0.1, 0.01, 0.1), Solid(10, 0.1, 0.01, 0.1);
    Fluid.calc_mass_matrix();
    Fluid.calc_K_matrix();
    Fluid.set_boundary_node(0,1.0);
    for(int i=0; i<100; i++){
        Fluid.time_step();
        Fluid.dump(i,"fluid");
    }
}
#include"diffusion.hpp"

using namespace std;

void onedimensinal_diffusion::calc_mass_matrix()
{
    for(int i=0; i<element.size(); i++){
        mass[element[i][0]][element[i][0]]+=1.0/3.0*element_length;
        mass[element[i][0]][element[i][1]]+=1.0/6.0*element_length;
        mass[element[i][1]][element[i][0]]+=1.0/6.0*element_length;
        mass[element[i][1]][element[i][1]]+=1.0/3.0*element_length;
    }

    for(int i=0; i<mass.size(); i++){
        double sum =0.0;
        for(int j=0; j<mass[i].size(); j++){
            sum += mass[i][j];
        }
        for(int j=0; j<mass[i].size(); j++){
            if(i==j) mass[i][j] = 1.0/sum;
            else mass[i][j]=0.0;
        }
    }
}

void onedimensinal_diffusion::calc_K_matrix()
{
    for(int i=0; i<element.size(); i++){
        K[element[i][0]][element[i][0]] += diffusion_coefficient*1.0/element_length;
        K[element[i][0]][element[i][1]] += -diffusion_coefficient*1.0/element_length;
        K[element[i][1]][element[i][0]] += -diffusion_coefficient*1.0/element_length;
        K[element[i][1]][element[i][1]] += diffusion_coefficient*1.0/element_length;
    }
}

void onedimensinal_diffusion::set_boundary_node(int node_number, double value)
{
    b_n.push_back(node_number);
    i_b_n.push_back(value);
}

void onedimensinal_diffusion::time_step()
{
    phi[b_n[0]]=i_b_n[0];
    vector<double> Dphi(phi.size(),0.0);
    vector<double> MDphi(phi.size(),0.0);
    for(int i=0; i<K.size(); i++){
        for(int j=0; j<phi.size(); j++){
            Dphi[i] += K[i][j]*phi[j];
        }
    }
    for(int i=0; i<mass.size(); i++){
        for(int j=0; j<Dphi.size(); j++){
            MDphi[i] += mass[i][j]*Dphi[j];
        }
    }
    for(int i=0; i<phi.size(); i++){
        phi[i]=phi[i]-dt*MDphi[i];
    }

    for(int i=0; i<phi.size(); i++){
        cout << phi[i] << " ";
    }
    cout << endl;
}

void onedimensinal_diffusion::dump(int step, std::string name)
{
    string filename = name + "_" + to_string(step) + ".dat";
    ofstream ofs(filename);
    for(int i=0; i<phi.size(); i++){
        ofs << phi[i] << endl;
    }
    ofs.close();
}
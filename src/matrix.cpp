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

    //for(int i=0; i<mass.size(); i++){
    //    double sum =0.0;
    //    for(int j=0; j<mass[i].size(); j++){
    //        sum += mass[i][j];
    //    }
    //    for(int j=0; j<mass[i].size(); j++){
    //        if(i==j) mass[i][j] = 1.0/sum;
    //        else mass[i][j]=0.0;
    //    }
    //}
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

void onedimensinal_diffusion::set_initial_boundary_node(int node_number, double value)
{
    b_n.push_back(node_number);
    i_b_n.push_back(value);
    c[b_n[0]]=i_b_n[0];
}

void onedimensinal_diffusion::boundary_setting(std::vector<double> boundary)
{
    c[b_n[0]]=i_b_n[0];
}

void onedimensinal_diffusion::time_step(std::vector<double> boundary)
{
    //boundary_setting(boundary);
    vector<double> Dc(c.size(),0.0);
    vector<double> MDc(c.size(),0.0);
    for(int i=0; i<K.size(); i++){
        for(int j=0; j<c.size(); j++){
            Dc[i] += K[i][j]*c[j];
        }
    }
    for(int i=0; i<mass.size(); i++){
        for(int j=0; j<Dc.size(); j++){
            MDc[i] += mass[i][j]*Dc[j];
        }
    }
    for(int i=0; i<c.size(); i++){
        c[i]=c[i]-dt*MDc[i];
    }
    for(int i=0; i<boundary.size(); i++){
        c[i] += boundary[i];
        //cout << c[i] << endl;
    }
    //boundary_setting(boundary);
}

void onedimensinal_diffusion::dump(int step, std::string name)
{
    string result_folder = name;
    mkdir(result_folder.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);

    string filename = name+"/"+to_string(step) + ".dat";
    ofstream ofs(filename);
    for(int i=0; i<c.size(); i++){
        ofs << c[i] << endl;
    }
    ofs.close();
}

double onedimensinal_diffusion::access_c(int i)
{
    return c[i];
}

void onedimensinal_diffusion::setting_phi(double p)
{
    for(int i=0; i<phi.size(); i++) phi[i] = p;
}

double onedimensinal_diffusion::return_sum()
{
    double sum =0.0;
    for(int i=0; i<c.size(); i++){
        sum += c[i];
    }
    return sum;
}

#include<iostream>
#include<vector>
#include<fstream>
#include<string>

using namespace std;

int CountNumbersOfTextLines(string &filePath)
{
  long i = 0;
  ifstream ifs( filePath );
  if(!ifs){
    cout << "can't file open!" << endl;
    exit(1);
  }
  if( ifs ){
    string line;
    while( true ){
      getline( ifs, line );
      i++;
      if( ifs.eof() )
        break;
    }
  }
  return i-1;
}

int main()
{
    //string point_sdf = "point_sdf.csv";
    //string chorid_point = "chroid_point.csv";
//
    //int sdf_line = CountNumbersOfTextLines(point_sdf);
    //int choroid_line = CountNumbersOfTextLines(chorid_point);
//
    //vector<double> point_sdf_fluid;
    //vector<double> point_sdf_solid;
    //vector<double> point_sdf_vessel;
//
    //vector<double> read_sdf;
    //vector<int> choroid_point_id;
    //string str;
    //ifstream ifs(point_sdf);
    //for(int i=0; i<sdf_line; i++){
    //    getline(ifs,str);
    //    if(i==0) continue;
    //    read_sdf.push_back(stod(str));
    //}   
    //ifs.close();
//
    //ifs.open(chorid_point);
    //for(int i=0; i<choroid_line; i++){
    //    getline(ifs,str);
    //    if(i==0) continue;
    //    choroid_point_id.push_back(stoi(str));
    //}   
    //ifs.close();
//
    //for(int i=0; i<read_sdf.size(); i++){
    //    if(read_sdf[i]<-3.0){
    //        point_sdf_fluid.push_back(0.0);
    //        point_sdf_solid.push_back(0.0);
    //        point_sdf_vessel.push_back(0.0);
    //        continue;
    //    }
    //    else{
    //        point_sdf_solid.push_back(1.0);
    //        point_sdf_fluid.push_back(1.0);
    //        point_sdf_vessel.push_back(1.0);
    //    }
    //    //if(read_sdf[i]<0.0){
    //    //    point_sdf_solid.push_back(0.99);
    //    //    point_sdf_fluid.push_back(0.0);
    //    //    point_sdf_vessel.push_back(0.01);
    //    //}
    //    //else{
    //    //    point_sdf_solid.push_back(1.0-(read_sdf[i]));
    //    //    point_sdf_fluid.push_back(read_sdf[i]-0.05);
    //    //    point_sdf_vessel.push_back(0.05);
    //    //}
    //}
//
    //for(int i=0; i<choroid_point_id.size(); i++){
    //    point_sdf_fluid[choroid_point_id[i]] = 0.2;
    //    point_sdf_solid[choroid_point_id[i]] = 0.0;
    //    point_sdf_vessel[choroid_point_id[i]] = 0.8;
    //}
//
    //ofstream ofs("fluid_sdf.dat");
    //for(int i=0; i<point_sdf_fluid.size(); i++){
    //    ofs << point_sdf_fluid[i] << endl;
    //}
    //ofs.close();
//
    //ofs.open("solid_sdf.dat");
    //for(int i=0; i<point_sdf_solid.size(); i++){
    //    ofs << point_sdf_solid[i] << endl;
    //}
    //ofs.close();
//
    //ofs.open("vessel_sdf.dat");
    //for(int i=0; i<point_sdf_vessel.size(); i++){
    //    ofs << point_sdf_vessel[i] << endl;
    //}
    //ofs.close();
//
    vector<int> source_node;
    string str;

    ifstream ifs("source_node.csv");
    int count =0 ;
    while(getline(ifs,str)){
        count++;
        if(count==1) continue;
        source_node.push_back(stoi(str));
    }

    ofstream ofs("boundary_vessel.dat");
    for(int i=0; i<source_node.size(); i++){
        if(source_node[i]>=250) ofs << i << " " << 1.0 << endl;
    }
    ofs.close();


}
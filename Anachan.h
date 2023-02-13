#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <format>

#ifndef Anachan_h
#define Anachan_h

using namespace std;

class Anachan
{
    enum {
        WD = 17,
        HT = 41,
    };
public:
    double gain[HT][WD];
    double resolution[HT][WD];
    int    ananumber;
    string gain_file_name;
    string resolution_file_name;
    string ana_dir_name;

    //graph
    TGraphErrors *gr_cnc1, *gr_cnc2, *gr_lg000, *gr_lg115, *gr_lg230, *gr_lg345;

    Anachan(int id=0);
    ~Anachan();

    virtual void Load(int id);
    virtual int  Save_gain(int id);
    virtual void Save_gain_All();
    virtual int  Save_reso(int id, int cut=0);
    virtual void Save_reso_All(int cut=0);
    virtual void Ana_gain(int id);
    virtual double Get_Gain_Value(int id, int runno0, int runno1);
    virtual void Ana_resolution(int id);

private:
    int cnc_run0[5] = {27, 30, 32, 35, 37};
    int cnc_run1[5] = {98, 99, 100, 101, 102};
    int lg000_run[9]; //=
    int lg115_run[9];
    int lg230_run[6];
    int lg345_run[6];
};

#endif


#ifdef Anachan_cxx
Anachan::Anachan(int id=0)
{   //initiarise
    for(int i=0; i<HT; i++){
        for(int j=0; j<WD; j++){
            gain[i][j] = -1.;
            resolution[i][j] = -1.;
        }
    }
    
    ananumber = id;
    
    ana_dir_name = "./anachan/";

    Load(ananumber);
}

Anachan::~Anachan()
{
}

void Anachan::Load(int id){
    //return;
    for(int i=0; i<HT; i++){
        for(int j=0; j<WD; j++){
            gain[i][j] = -1.;
            resolution[i][j] = -1.;
        }
    }
    ananumber = id;

    if(id==0){
        gain_file_name = "./gain.txt";
        resolution_file_name = "./resolution.txt";
    } else {
        gain_file_name = "./gain" + to_string(id) + ".txt";
        resolution_file_name = "./resolution" + to_string(id) + ".txt";
    }

    ifstream gain_file(gain_file_name.c_str());
    if(gain_file){
        for(int i=0; i<HT; i++){
            for(int j=0; j<WD; j++){
                gain_file >> gain[i][j];
            }
        }
        
        cout << "Load " << gain_file_name << endl;
    } else {
        cout << "Cant Load " << gain_file_name << endl;
    }
    gain_file.close();

    ifstream resolution_file(resolution_file_name.c_str());
    if(resolution_file){
        for(int i=0; i<HT; i++){
            for(int j=0; j<WD; j++){
                resolution_file >> resolution[i][j];
            }
        }

        cout << "Load " << resolution_file_name << endl;
    } else {
        cout << "Cant Load " << resolution_file_name << endl;
    }
    resolution_file.close();

}

#endif // #ifdef Anachan_cxx
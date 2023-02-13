#define Anachan_cxx
#include "Anachan.h"
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <format>



int Anachan::Save_gain(int id){
    //return;
    string ofilename, pmfilename;
    double runno = (double)id + 100000000;
    char runfname[16];
    sprintf(runfname,"%9.0lf",runno);
    pmfilename = "./ana/" + to_string(id) + "/run" + string(runfname) + ".ped_mip";
    ofilename = "./gain.txt";

    double pedestal[16], ped_error[16], mip[16];
    ifstream pmfile(pmfilename.c_str());
    if(pmfile){
        for(int i=0; i<16; i++){
            pmfile >> pedestal[i] >> ped_error[i] >> mip[i];
        }
    } else {
        //cout << "no file" << endl;
        return -1;
    }
    pmfile.close();

    ofstream ofile(ofilename.c_str(),ios::app);
    if(ofile){
        ofile << id << " ";
        for(int i=0; i<16; i++){
            if(pedestal[i]!=0. && mip[i]!=0){
                ofile << (mip[i] - pedestal[i]) << " ";
            } else {
                ofile << -1 << " ";
            }
        }
        ofile << endl;
    }
    ofile.close();
    return 0;
}

void Anachan::Save_gain_All(){
    //return;
    int max = 102;
    int gain;
    for(int i=0; i<max; i++){
        gain = Save_gain(i+1);
        if(gain==0){
            cout << "comprete " << (i+1) << endl;
        } else {
            cout << "file " << (i+1) << " does not exist" << endl;
        }
    }
}

int Anachan::Save_reso(int id, int cut=0){
    string ofilename, pmfilename;
    double runno = (double)id + 100000000;
    char runfname[16];
    sprintf(runfname,"%9.0lf",runno);
    pmfilename = "./ana/" + to_string(id) + "/run" + string(runfname) + ".reso";
    ofilename = "./resolution.txt";

    double num[16], reso[16];
    ifstream pmfile(pmfilename.c_str());
    if(pmfile){
        for(int i=0; i<16; i++){
            pmfile >> num[i] >> reso[i];
        }
    } else {
        //cout << "no file" << endl;
        return -1;
    }
    pmfile.close();

    ofstream ofile(ofilename.c_str(),ios::app);
    if(ofile){
        ofile << (id + cut*1000) << " ";
        for(int i=0; i<16; i++){
            ofile << reso[i] << " ";
        }
        ofile << endl;
    }
    ofile.close();
    return 0;
}

void Anachan::Save_reso_All(int cut=0){
    //return;
    int max = 102;
    int reso;
    for(int i=0; i<max; i++){
        reso = Save_reso(i+1,cut);
        if(reso==0){
            cout << "comprete " << (i+1) << endl;
        } else {
            cout << "file " << (i+1) << " does not exist" << endl;
        }
    }
}

void Anachan::Ana_gain(int id = 1){
    //return;
    //cnc run  98-102(36-40)   or   27-37()
    int work[5];
    if      (id==0) for(int i=0; i<5; i++){work[i] = cnc_run0[i];}
    else if (id==1) for(int i=0; i<5; i++){work[i] = cnc_run1[i];}

    int start = -1200;
    int interval = 300;
    const int NPOINT = 9;
    double X[NPOINT], X_ERR[NPOINT], GAIN1[NPOINT], GAIN1_ERR[NPOINT], GAIN2[NPOINT], GAIN2_ERR[NPOINT];
    double gain1_r_at0, gain1_l_at0, gain2_r_at0, gain2_l_at0;
    for(int i=0; i<NPOINT; i++){
        X[i] = start + interval*i;
        X_ERR[i] = 1;
        GAIN1[i] = 1.;   //initiarise
        GAIN2[i] = 1.;   //initiarise
        GAIN1_ERR[i] = 0.;
        GAIN2_ERR[i] = 0.;
    }

    int runno;
    for(int i=0; i<5; i++){
        if(1){
            runno = 36 + i;
            if(i==0){
                gain1_l_at0 = gain[runno][9];
                gain1_r_at0 = gain[runno][10];
                gain2_l_at0 = gain[runno][11];
                gain2_r_at0 = gain[runno][12];
                GAIN1[4] = 1.;
                GAIN2[4] = 1.;
            } else {
                GAIN1[4-i] = gain[runno][9]/gain1_l_at0;
                GAIN1[4+i] = gain[runno][10]/gain1_r_at0;
                GAIN2[4-i] = gain[runno][11]/gain2_l_at0;
                GAIN2[4+i] = gain[runno][12]/gain2_r_at0;
            }
        }
    }

    double start2 = -600;
    double interval2 = 300;
    int NPOINT2 = 5;
    double X2[NPOINT2], X2_ERR[NPOINT2];
    double GAIN000[NPOINT2], GAIN000_ERR[NPOINT2], GAIN115[NPOINT2], GAIN230[NPOINT2], GAIN345[NPOINT2];
    double at0[4][2];
    int runno2[4];
    for(int i=0; i<NPOINT2; i++){
        X2[i] = start2 + interval2*i;
        X2_ERR[i] = 1;
        GAIN000_ERR[i] = 0.;
    }
    for(int i=0; i<3; i++){
        if(1){
            runno2[0] = 6 + i;
            runno2[1] = 16 + i;
            runno2[2] = 31 + i;
            runno2[3] = 25 + i;
            if(i==0){
                GAIN000[2] = 1.;
                GAIN115[2] = 1.;
                GAIN230[2] = 1.;
                GAIN345[2] = 1.;
                for(int j=0; j<4; j++){
                    at0[j][0] = gain[runno2[j]][9];
                    at0[j][1] = gain[runno2[j]][10];
                }
            } else {
                GAIN000[2-i] = gain[runno2[0]][9]/at0[0][0];
                GAIN000[2+i] = gain[runno2[0]][10]/at0[0][1];
                GAIN115[2-i] = gain[runno2[1]][9]/at0[1][0];
                GAIN115[2+i] = gain[runno2[1]][10]/at0[1][1];
                GAIN230[2-i] = gain[runno2[2]][9]/at0[2][0];
                GAIN230[2+i] = gain[runno2[2]][10]/at0[2][1];
                GAIN345[2-i] = gain[runno2[3]][9]/at0[3][0];
                GAIN345[2+i] = gain[runno2[3]][10]/at0[3][1];
            }
        }
    }

    TLegend *leg = new TLegend(0.1,0.75,0.2,0.95);
    gr_cnc1 = new TGraphErrors(NPOINT,X,GAIN1,X_ERR,GAIN1_ERR);
    gr_cnc2 = new TGraphErrors(NPOINT,X,GAIN1,X_ERR,GAIN1_ERR);
    gr_lg000 = new TGraphErrors(NPOINT2,X2,GAIN000,X2_ERR,GAIN000_ERR);
    gr_lg115 = new TGraphErrors(NPOINT2,X2,GAIN115,X2_ERR,GAIN000_ERR);
    gr_lg230 = new TGraphErrors(NPOINT2,X2,GAIN230,X2_ERR,GAIN000_ERR);
    gr_lg345 = new TGraphErrors(NPOINT2,X2,GAIN345,X2_ERR,GAIN000_ERR);

    leg->AddEntry(gr_cnc1,"cnc1","l");
    leg->AddEntry(gr_cnc2,"cnc2","l");
    leg->AddEntry(gr_lg000,"lg0","l");
    leg->AddEntry(gr_lg115,"lg115","l");
    leg->AddEntry(gr_lg230,"lg230","l");
    leg->AddEntry(gr_lg345,"lg345","l");

    int num = 2;
    TCanvas *c1 = new TCanvas("c1","c1");
    //gPad->SetLogy();
    gr_cnc1->SetMaximum(1.5);
    gr_cnc1->SetMinimum(0.8);
    /* gr_cnc1->GetXAxis()->SetTitle("position(mm)");
    gr_cnc1->GetYAxis()->SetTitle("gain(1/center)"); */
    gr_cnc1->SetTitle("gain vs position;position(mm);gain(1/center)");

    gr_cnc1->SetMarkerStyle(2);
    gr_cnc1->SetMarkerSize(1);
    gr_cnc1->SetMarkerColor(num);
    gr_cnc1->SetLineColor(num);
    gr_cnc1->Draw("APL");
    num++;
    gr_cnc2->SetMarkerStyle(2);
    gr_cnc2->SetMarkerSize(1);
    gr_cnc2->SetMarkerColor(num);
    gr_cnc2->SetLineColor(num);
    gr_cnc2->Draw("PL");
    num++;
    gr_lg000->SetMarkerStyle(2);
    gr_lg000->SetMarkerSize(1);
    gr_lg000->SetMarkerColor(num);
    gr_lg000->SetLineColor(num);
    gr_lg000->Draw("PL");
    num++;
    gr_lg115->SetMarkerStyle(2);
    gr_lg115->SetMarkerSize(1);
    gr_lg115->SetMarkerColor(num);
    gr_lg115->SetLineColor(num);
    gr_lg115->Draw("PL");
    num++;
    gr_lg230->SetMarkerStyle(2);
    gr_lg230->SetMarkerSize(1);
    gr_lg230->SetMarkerColor(num);
    gr_lg230->SetLineColor(num);
    gr_lg230->Draw("PL");
    num++;
    gr_lg345->SetMarkerStyle(2);
    gr_lg345->SetMarkerSize(1);
    gr_lg345->SetMarkerColor(num);
    gr_lg345->SetLineColor(num);
    gr_lg345->Draw("PL");

    leg->Draw();

    string ofname;
    ofname = "./pdf/gain_vs_position.pdf";
    c1->Print(ofname.c_str());

}

double Anachan::Get_Gain_Value(int id, int runno0, int runno1){
    return 0;
}

void Anachan::Ana_resolution(int id=1){
    //return;
    //cnc run  98-102(36-40)   or   27-37()
    int work[5];
    if      (id==0) for(int i=0; i<5; i++){work[i] = cnc_run0[i];}
    else if (id==1) for(int i=0; i<5; i++){work[i] = cnc_run1[i];}

    int start = -1200;
    int interval = 300;
    const int NPOINT = 9;
    double X[NPOINT], X_ERR[NPOINT], resolution1[NPOINT], resolution1_ERR[NPOINT], resolution2[NPOINT], resolution2_ERR[NPOINT];
    double resolution1_r_at0, resolution1_l_at0, resolution2_r_at0, resolution2_l_at0;
    for(int i=0; i<NPOINT; i++){
        X[i] = start + interval*i;
        X_ERR[i] = 1;
        resolution1[i] = 1.;   //initiarise
        resolution2[i] = 1.;   //initiarise
        resolution1_ERR[i] = 0.;
        resolution2_ERR[i] = 0.;
    }

    int runno;
    for(int i=0; i<5; i++){
        if(1){
            runno = 36 + i;
            if(i==0){
                resolution1_l_at0 = resolution[runno][9];
                resolution1_r_at0 = resolution[runno][10];
                resolution2_l_at0 = resolution[runno][11];
                resolution2_r_at0 = resolution[runno][12];
                resolution1[4] = 1.;
                resolution2[4] = 1.;
            } else {
                resolution1[4-i] = resolution[runno][9]/resolution1_l_at0;
                resolution1[4+i] = resolution[runno][10]/resolution1_r_at0;
                resolution2[4-i] = resolution[runno][11]/resolution2_l_at0;
                resolution2[4+i] = resolution[runno][12]/resolution2_r_at0;
            }
        }
    }

    double start2 = -600;
    double interval2 = 300;
    int NPOINT2 = 5;
    double X2[NPOINT2], X2_ERR[NPOINT2];
    double resolution000[NPOINT2], resolution000_ERR[NPOINT2], resolution115[NPOINT2], resolution230[NPOINT2], resolution345[NPOINT2];
    double at0[4][2];
    int runno2[4];
    for(int i=0; i<NPOINT2; i++){
        X2[i] = start2 + interval2*i;
        X2_ERR[i] = 1;
        resolution000_ERR[i] = 0.;
    }
    for(int i=0; i<3; i++){
        if(1){
            runno2[0] = 6 + i;
            runno2[1] = 16 + i;
            runno2[2] = 31 + i;
            runno2[3] = 25 + i;
            if(i==0){
                resolution000[2] = 1.;
                resolution115[2] = 1.;
                resolution230[2] = 1.;
                resolution345[2] = 1.;
                for(int j=0; j<4; j++){
                    at0[j][0] = resolution[runno2[j]][9];
                    at0[j][1] = resolution[runno2[j]][10];
                }
            } else {
                resolution000[2-i] = resolution[runno2[0]][9]/at0[0][0];
                resolution000[2+i] = resolution[runno2[0]][10]/at0[0][1];
                resolution115[2-i] = resolution[runno2[1]][9]/at0[1][0];
                resolution115[2+i] = resolution[runno2[1]][10]/at0[1][1];
                resolution230[2-i] = resolution[runno2[2]][9]/at0[2][0];
                resolution230[2+i] = resolution[runno2[2]][10]/at0[2][1];
                resolution345[2-i] = resolution[runno2[3]][9]/at0[3][0];
                resolution345[2+i] = resolution[runno2[3]][10]/at0[3][1];
            }
        }
    }

    TLegend *leg = new TLegend(0.1,0.75,0.2,0.95);
    gr_cnc1 = new TGraphErrors(NPOINT,X,resolution1,X_ERR,resolution1_ERR);
    gr_cnc2 = new TGraphErrors(NPOINT,X,resolution1,X_ERR,resolution1_ERR);
    gr_lg000 = new TGraphErrors(NPOINT2,X2,resolution000,X2_ERR,resolution000_ERR);
    gr_lg115 = new TGraphErrors(NPOINT2,X2,resolution115,X2_ERR,resolution000_ERR);
    gr_lg230 = new TGraphErrors(NPOINT2,X2,resolution230,X2_ERR,resolution000_ERR);
    gr_lg345 = new TGraphErrors(NPOINT2,X2,resolution345,X2_ERR,resolution000_ERR);

    leg->AddEntry(gr_cnc1,"cnc1","l");
    leg->AddEntry(gr_cnc2,"cnc2","l");
    leg->AddEntry(gr_lg000,"lg0","l");
    leg->AddEntry(gr_lg115,"lg115","l");
    leg->AddEntry(gr_lg230,"lg230","l");
    leg->AddEntry(gr_lg345,"lg345","l");

    int num = 2;
    TCanvas *c1 = new TCanvas("c1","c1");
    //gPad->SetLogy();
    gr_cnc1->SetMaximum(1.5);
    gr_cnc1->SetMinimum(0.8);
    /* gr_cnc1->GetXAxis()->SetTitle("position(mm)");
    gr_cnc1->GetYAxis()->SetTitle("resolution(1/center)"); */
    gr_cnc1->SetTitle("resolution vs position;position(mm);resolution(1/center)");

    gr_cnc1->SetMarkerStyle(2);
    gr_cnc1->SetMarkerSize(1);
    gr_cnc1->SetMarkerColor(num);
    gr_cnc1->SetLineColor(num);
    gr_cnc1->Draw("APL");
    num++;
    gr_cnc2->SetMarkerStyle(2);
    gr_cnc2->SetMarkerSize(1);
    gr_cnc2->SetMarkerColor(num);
    gr_cnc2->SetLineColor(num);
    gr_cnc2->Draw("PL");
    num++;
    gr_lg000->SetMarkerStyle(2);
    gr_lg000->SetMarkerSize(1);
    gr_lg000->SetMarkerColor(num);
    gr_lg000->SetLineColor(num);
    gr_lg000->Draw("PL");
    num++;
    gr_lg115->SetMarkerStyle(2);
    gr_lg115->SetMarkerSize(1);
    gr_lg115->SetMarkerColor(num);
    gr_lg115->SetLineColor(num);
    gr_lg115->Draw("PL");
    num++;
    gr_lg230->SetMarkerStyle(2);
    gr_lg230->SetMarkerSize(1);
    gr_lg230->SetMarkerColor(num);
    gr_lg230->SetLineColor(num);
    gr_lg230->Draw("PL");
    num++;
    gr_lg345->SetMarkerStyle(2);
    gr_lg345->SetMarkerSize(1);
    gr_lg345->SetMarkerColor(num);
    gr_lg345->SetLineColor(num);
    gr_lg345->Draw("PL");

    leg->Draw();

    string ofname;
    ofname = "./pdf/resolution_vs_position.pdf";
    c1->Print(ofname.c_str());
}
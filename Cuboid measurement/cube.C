#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <fstream>
#include <TString.h>
#include "TMath.h"

#endif

using namespace TMath;
using namespace std;

// Declare functions
void histogram(TH1D*, const TString, TCanvas*, const TString, const TString, const TString);

// Initialize histograms
TH1D *histCubeX = new TH1D("histCubeX", "histCubeX", 25, 8, 10);
TH1D *histCubeY = new TH1D("histCubeY", "histCubeY", 25, 9, 12);
TH1D *histCubeZ = new TH1D("histCubeZ", "histCubeZ", 25, 35, 39);
TH1D *histCubeV = new TH1D("histCubeV", "histCubeV", 25, 3100, 3800);

/*
 * MAIN FUNCTION
 */

void cube(TString inputFile = "cubeData.txt"){
    
    cout << "\n\nStarting process...\n\n";
    
    ifstream ifs(inputFile); if(!ifs.is_open()){cout << "Error. File " << inputFile << " not found. Exiting...\n"; return;}

    TString id;
    Double_t cubeNumber, cubeX, cubeY, cubeZ, cubeV, cubeXErr, cubeYErr, cubeZErr, cubeVErr;
    
    while(ifs >> id >> cubeNumber >> cubeX >> cubeXErr >> cubeY >> cubeYErr >> cubeZ >> cubeZErr >> cubeV >> cubeVErr){
        
        if(id.Contains("#")) continue;
        
        histCubeX->Fill(cubeX);
        histCubeY->Fill(cubeY);
        histCubeZ->Fill(cubeZ);
        histCubeV->Fill(cubeV);
        
    }
    
    TCanvas *c1 = new TCanvas("Histogram", "Histogram", 1000, 600);

    //gStyle->SetOptStat(kFALSE);

    histCubeX->Fit("gaus");
    histCubeY->Fit("gaus");
    histCubeZ->Fit("gaus");
    histCubeV->Fit("gaus");

    histogram(histCubeX, "", c1, "X measurements", "Count", "cubeX");
    histogram(histCubeY, "", c1, "Y measurements", "Count", "cubeY");
    histogram(histCubeZ, "", c1, "Z measurements", "Count", "cubeZ");
    histogram(histCubeV, "", c1, "V measurements", "Count", "cubeV");
    
}

/*
 * FUNCTION FOR SAVING ONE HISTOGRAM
 */

void histogram(TH1D *histo, const TString histName, TCanvas *can, const TString xTitle, const TString yTitle, const TString name){

    histo->SetLineWidth(3);
    histo->Draw("E1");
    // add axis labels
    histo->GetXaxis()->SetTitle(xTitle);
    histo->GetYaxis()->SetTitle(yTitle);
    histo->SetTitle(histName); // title on top

    can->SaveAs(name + ".jpg");
}

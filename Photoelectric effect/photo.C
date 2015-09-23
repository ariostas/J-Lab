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
#include <TGraph.h>

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

void photo(TString inputFile = "photoData365.txt"){
    
    cout << "\n\nStarting process...\n\n";
    
    ifstream ifs(inputFile); if(!ifs.is_open()){cout << "Error. File " << inputFile << " not found. Exiting...\n"; return;}

    Double_t voltage[100], current[100], tempVoltage, tempCurrent;
    Int_t x=0;
    
    while(ifs >> tempVoltage >> tempCurrent){
        
        voltage[x] = tempVoltage;
        current[x] = tempCurrent;
        
        x++;
    }
    
    TCanvas *c1 = new TCanvas("Histogram", "Histogram", 1000, 600);

    //gStyle->SetOptStat(kFALSE);

    TGraph *graph = new TGraph(x, voltage, current);

    graph->SetLineWidth(3);
    graph->Draw();
    graph->GetXaxis()->SetTitle("Voltage (V)");
    graph->GetYaxis()->SetTitle("Current (pA)");
    graph->SetTitle("Photo electric effect");
    
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

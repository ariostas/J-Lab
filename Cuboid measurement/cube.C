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
#include <TPaveStats.h>
#include <TText.h>
#include <TLatex.h>

#endif

using namespace TMath;
using namespace std;

// Declare functions
void histogram(TH1D*, const TString, TCanvas*, const TString, const TString, const TString);

// Initialize histograms
TH1D *histCubeX = new TH1D("histCubeX", "histCubeX", 15, 8, 10);
TH1D *histCubeY = new TH1D("histCubeY", "histCubeY", 15, 9, 12);
TH1D *histCubeZ = new TH1D("histCubeZ", "histCubeZ", 15, 34.5, 39);
TH1D *histCubeV = new TH1D("histCubeV", "histCubeV", 15, 2880, 4131);

/*
 * MAIN FUNCTION
 */

void cube(TString inputFile = "cubeData.txt"){

    TH1::StatOverflows(kTRUE);
    
    cout << "\n\nStarting process...\n\n";
    
    ifstream ifs(inputFile); if(!ifs.is_open()){cout << "Error. File " << inputFile << " not found. Exiting...\n"; return;}

    TString id;
    Double_t cubeNumber, cubeX, cubeY, cubeZ, cubeV, cubeXErr, cubeYErr, cubeZErr, cubeVErr;

    Double_t mean = 0, nEntries = 0;
    
    while(ifs >> id >> cubeNumber >> cubeX >> cubeXErr >> cubeY >> cubeYErr >> cubeZ >> cubeZErr >> cubeV >> cubeVErr){
        
        if(id.Contains("#")) continue;
        
        histCubeX->Fill(cubeX);
        histCubeY->Fill(cubeY);
        histCubeZ->Fill(cubeZ);
        histCubeV->Fill(cubeV);

        mean += cubeX;
        nEntries++;
    }

    mean = mean/nEntries;

    Double_t sumSquared =0;

    ifs.clear();
    ifs.seekg(0, std::ios::beg);

    while(ifs >> id >> cubeNumber >> cubeX >> cubeXErr >> cubeY >> cubeYErr >> cubeZ >> cubeZErr >> cubeV >> cubeVErr){
        
        sumSquared += (cubeX-mean)*(cubeX-mean);
    }

    sumSquared = sumSquared/(nEntries-1);

    cout << TMath::Sqrt(sumSquared) << endl;
    
    TCanvas *c1 = new TCanvas("Histogram", "Histogram", 1600, 900);

    gStyle->SetOptStat(2210);
    //gStyle->SetOptStat(kFALSE);
    gStyle->SetOptFit(1111);

    histCubeX->Fit("gaus", "M");
    histCubeY->Fit("gaus", "M");
    histCubeZ->Fit("gaus", "M");
    histCubeV->Fit("gaus", "M");

    histogram(histCubeX, "", c1, "Length of X side (mm)", "Count", "cubeX");
    histogram(histCubeY, "", c1, "Length of Y side (mm)", "Count", "cubeY");
    histogram(histCubeZ, "", c1, "Length of Z side (mm)", "Count", "cubeZ");
    histogram(histCubeV, "", c1, "Volume of cuboid (mm^{3})", "Count", "cubeV");
    
}

/*
 * FUNCTION FOR SAVING ONE HISTOGRAM
 */

void histogram(TH1D *histo, const TString histName, TCanvas *can, const TString xTitle, const TString yTitle, const TString name){

    histo->SetLineWidth(3);
    histo->Draw("E1");

    TF1 *f1 = new TF1(histName,TString::Format("%4.3f*TMath::Gaus(x,%4.3f,%4.3f,1)",
        histo->Integral()*histo->GetBinWidth(1), histo->GetMean(), histo->GetStdDev()),
        histo->GetXaxis()->GetXmin(), histo->GetXaxis()->GetXmax()); 
    f1->SetLineWidth(3);
    f1->SetLineColor(kRed);
    //f1->Draw("same");

    TPaveText *pt = new TPaveText(0.605,0.675,0.885,0.875, "brNDC");
    // pt->AddText(TString::Format("Entries = %3i", histo->Integral()));
    pt->AddText("Total entries = 103");
    pt->AddText(TString::Format("Mean = %4.1f #pm %4.1f_{stat} #pm 0.5_{syst} mm", histo->GetMean(), Max(histo->GetMeanError(),0.1)));
    pt->AddText(TString::Format("Std dev = %4.2f mm", histo->GetStdDev()));
    //pt->Draw();

    // add axis labels
    histo->GetXaxis()->SetTitle(xTitle);
    histo->GetXaxis()->SetTitleSize(0.045);
    histo->GetXaxis()->SetTitleOffset(1.05);
    histo->GetXaxis()->SetLabelOffset(0.010);
    histo->GetXaxis()->SetLabelSize(0.05);
    histo->GetYaxis()->SetTitle(yTitle);
    histo->GetYaxis()->SetTitleSize(0.045);
    histo->GetYaxis()->SetTitleOffset(0.70);
    histo->GetYaxis()->SetLabelSize(0.05);
    histo->SetTitle(histName); // title on top
    can->Update();

    can->SaveAs(name + ".png");
}
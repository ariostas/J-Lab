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
TH1D *histCubeY = new TH1D("histCubeY", "histCubeY", 12, 8.5, 12.5);
TH1D *histCubeZ = new TH1D("histCubeZ", "histCubeZ", 15, 34.5, 39);
TH1D *histCubeV = new TH1D("histCubeV", "histCubeV", 14, 2.8, 4.2);

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
        histCubeV->Fill(cubeV/1000.0);

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

    c1->Divide(2,2);

    c1->cd(1);
    histogram(histCubeX, "", c1, "Length of X side (mm)", "Count", "cubeX");
    c1->cd(2);
    histogram(histCubeY, "", c1, "Length of Y side (mm)", "Count", "cubeY");
    c1->cd(3);
    histogram(histCubeZ, "", c1, "Length of Z side (mm)", "Count", "cubeZ");
    c1->cd(4);
    histogram(histCubeV, "", c1, "Volume of cuboid (cm^{3})", "Count", "cubeV");
    
}

/*
 * FUNCTION FOR SAVING ONE HISTOGRAM
 */

void histogram(TH1D *histoData, const TString histName, TCanvas *can, const TString xTitle, const TString yTitle, const TString name){

    histoData->SetLineWidth(2);
    histoData->SetMarkerStyle(20);
    histoData->SetMarkerSize(1.5);
    histoData->SetMarkerColor(kRed);
    histoData->SetLineColor(kBlack);
    TF1 *f = new TF1(histName, "[0]*TMath::Gaus(x,[1],[2])");
    f->SetParameter(0, histoData->GetMaximum());
    f->SetParameter(1, histoData->GetMean());
    f->SetParameter(2, histoData->GetStdDev());
    f->SetParName(0, "Constant");
    f->SetParName(1, "Mean");
    f->SetParName(2, "Sigma");
    f->SetLineColor(kBlue);
    f->SetLineWidth(4);
    histoData->Fit(f, "M");

    gStyle->SetLegendBorderSize(0);
    TLegend *leg = new TLegend(0.65,0.6,0.885,0.875);
    leg->SetTextSize(0.055);
    leg->AddEntry(histoData, "Data","lep");
    leg->AddEntry(f, "Gaussian fit","l");
    
    histoData->Draw("E1");
    leg->Draw("same");

    // add axis labels
    histoData->GetXaxis()->SetTitle(xTitle);
    histoData->GetXaxis()->CenterTitle();
    histoData->GetXaxis()->SetTitleSize(0.055);
    histoData->GetXaxis()->SetTitleOffset(0.87);
    histoData->GetXaxis()->SetLabelOffset(0.010);
    histoData->GetXaxis()->SetLabelSize(0.05);
    histoData->GetYaxis()->SetTitle(yTitle);
    histoData->GetYaxis()->CenterTitle();
    histoData->GetYaxis()->SetTitleSize(0.055);
    histoData->GetYaxis()->SetTitleOffset(0.70);
    histoData->GetYaxis()->SetLabelSize(0.05);
    gStyle->SetTitleSize(0.08, "t");
    histoData->SetTitle(histName);
    can->Update();

    can->SaveAs(name + ".png");
    can->SaveAs(name + ".pdf");
}

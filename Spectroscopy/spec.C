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
#include <TGraph.h>

#endif

using namespace TMath;
using namespace std;

// Declare functions
void histogram(TH1D*, const TString, TCanvas*, const TString, const TString, const TString);
void histogram(TH1D*, TH1D*, const TString, TCanvas*, const TString, const TString, const TString);
TH1D* readData(TString);
void calibration();

/*
 * MAIN FUNCTION
 */

void spec(){

    TH1::StatOverflows(kTRUE);
    
    cout << "\n\nStarting process...\n\n";
    
    TH1D *merc3650 = readData("mercury_3650_Oct13");
    TH1D *merc4046 = readData("mercury_4046_Oct13");
    TH1D *merc4358 = readData("mercury_4358_Oct13");
    TH1D *merc5460 = readData("mercury_5460_Oct13");
    TH1D *merc5769 = readData("mercury_5769_Oct13");
    TH1D *merc5790 = readData("mercury_5790_Oct13");

    TH1D *hydr6562 = readData("hydrogen_6562_Oct14_goodlamp");
    TH1D *hydr4861 = readData("hydrogen_4861_Oct14_goodlamp");
    TH1D *hydr4340 = readData("hydrogen_4340_Oct14_goodlamp");
    TH1D *hydr4101 = readData("hydrogen_4101_Oct14_goodlamp");

    TH1D *deut6562 = readData("deuterium_6562_Oct14_goodlamp");
    TH1D *deut4861 = readData("deuterium_4861_Oct14_goodlamp");
    TH1D *deut4340 = readData("deuterium_4340_Oct14_goodlamp");
    TH1D *deut4101 = readData("deuterium_4101_Oct14_goodlamp");

    TCanvas *can = new TCanvas("Histogram", "Histogram", 1600, 900);

    //gStyle->SetOptStat(2210);
    //gStyle->SetOptFit(1111);
    gStyle->SetOptStat(kFALSE);
    gStyle->SetOptFit(kFALSE);

    histogram(merc3650, "Mercury 3650", can, "Wavelength", "CPS", "Plot_Mercury_3650");
    histogram(merc4046, "Mercury 4046", can, "Wavelength", "CPS", "Plot_Mercury_4046");
    histogram(merc4358, "Mercury 4358", can, "Wavelength", "CPS", "Plot_Mercury_4358");
    histogram(merc5460, "Mercury 5460", can, "Wavelength", "CPS", "Plot_Mercury_5460");
    histogram(merc5769, "Mercury 5769", can, "Wavelength", "CPS", "Plot_Mercury_5769");
    histogram(merc5790, "Mercury 5790", can, "Wavelength", "CPS", "Plot_Mercury_5790");

    histogram(hydr6562, "Hydrogen 6562", can, "Wavelength", "CPS", "Plot_Hydrogen_6562");
    histogram(hydr4861, "Hydrogen 4861", can, "Wavelength", "CPS", "Plot_Hydrogen_4861");
    histogram(hydr4340, "Hydrogen 4340", can, "Wavelength", "CPS", "Plot_Hydrogen_4340");
    histogram(hydr4101, "Hydrogen 4101", can, "Wavelength", "CPS", "Plot_Hydrogen_4101");

    histogram(deut6562, "Deuterium 6562", can, "Wavelength", "CPS", "Plot_Deuterium_6562");
    histogram(deut4861, "Deuterium 4861", can, "Wavelength", "CPS", "Plot_Deuterium_4861");
    histogram(deut4340, "Deuterium 4340", can, "Wavelength", "CPS", "Plot_Deuterium_4340");
    histogram(deut4101, "Deuterium 4101", can, "Wavelength", "CPS", "Plot_Deuterium_4101");

    gStyle->SetOptStat(kFALSE);
    gStyle->SetOptFit(kFALSE);

    histogram(hydr6562, deut6562, "Hydrogen vs Deuterium 6562", can, "Wavelength", "CPS", "Plot_Hydr_vs_Deut_6562");
    histogram(hydr4861, deut4861, "Hydrogen vs Deuterium 4861", can, "Wavelength", "CPS", "Plot_Hydr_vs_Deut_4861");
    histogram(hydr4340, deut4340, "Hydrogen vs Deuterium 4340", can, "Wavelength", "CPS", "Plot_Hydr_vs_Deut_4340");
    histogram(hydr4101, deut4101, "Hydrogen vs Deuterium 4101", can, "Wavelength", "CPS", "Plot_Hydr_vs_Deut_4101");

    gStyle->SetOptStat(2210);
    gStyle->SetOptFit(1111);

    calibration();
}

TH1D* readData(TString name){

    TString fileName = name + ".txt";

    Double_t histoMin = 99999, histoMax = 0;

    ifstream ifs(fileName); if(!ifs.is_open()){cout << "Error. File " << fileName << " not found. Exiting...\n"; return NULL;}

    Double_t entry, wavelength, photons;
    
    while(ifs >> entry >> wavelength >> photons){
        
        if(wavelength < histoMin) histoMin = wavelength;
        if(wavelength > histoMax) histoMax = wavelength;

    }

    ifs.close();

    TH1D *histo = new TH1D(name, name, Nint(entry), histoMin, histoMax);


    ifs.open(fileName); if(!ifs.is_open()){cout << "Error. File " << fileName << " not found. Exiting...\n"; return NULL;}
    
    while(ifs >> entry >> wavelength >> photons){
        
        for(Int_t x = 0; x< photons; x++){
            histo->Fill(wavelength);
        }

    }

    ifs.close();

    return histo;

}

void calibration(){

    TH1::StatOverflows(kTRUE);

    Double_t wavelength[10], marker[10];
    // From maximum values
    wavelength[0] = 3131.548; marker[0] = 3134.08;
    wavelength[1] = 3131.839; marker[1] = 3134.36;
    wavelength[2] = 3650.153; marker[2] = 3649.61;
    wavelength[3] = 4046.563; marker[3] = 4043.30;
    wavelength[4] = 4358.328; marker[4] = 4352.77;
    wavelength[5] = 5460.735; marker[5] = 5445.42;
    wavelength[6] = 5769.598; marker[6] = 5751.17;
    wavelength[7] = 5790.663; marker[7] = 5771.98;

    // From fits
    // wavelength[0] = 3650.153; marker[0] = 3649.60;
    // wavelength[1] = 4046.563; marker[1] = 4043.29;
    // wavelength[2] = 4358.328; marker[2] = 4352.77;
    // wavelength[3] = 5460.735; marker[3] = 5445.41;
    // wavelength[4] = 5769.598; marker[4] = 5751.18;
    // wavelength[5] = 5790.663; marker[5] = 5771.98;
    
    TGraph *graph = new TGraph(8, marker, wavelength);
    graph->Fit("pol2", "FM");
    
    TCanvas *c1 = new TCanvas("Calibration", "Calibration", 1600, 900);

    gStyle->SetOptStat(2210);
    gStyle->SetOptFit(1111);

    graph->Draw("AC*"); 
    c1->SaveAs("Plot_Calibration.png");
    
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
    TF1 *f = new TF1(histName, "[0]+[1]*TMath::Voigt(x-[2],[3],[4])");
    f->SetParameter(0, 200);
    f->SetParameter(1, histoData->GetMaximum()/10.);
    f->SetParameter(2, histoData->GetMean());
    f->SetParameter(3, .05);
    f->SetParameter(4, .05);
    f->SetParName(0, "Offset");
    f->SetParName(1, "Normalization");
    f->SetParName(2, "Mean");
    f->SetParName(3, "Sigma");
    f->SetParName(4, "gamma");
    f->SetLineColor(kBlue);
    f->SetLineWidth(4);
    histoData->Fit(f, "ME");

    gStyle->SetLegendBorderSize(0);
    TLegend *leg = new TLegend(0.15,0.6,0.385,0.875);
    leg->SetTextSize(0.055);
    leg->AddEntry(histoData, "Data","lep");
    leg->AddEntry(f, "Voigt fit","l");
    
    histoData->Draw("E1");
    leg->Draw("same");

    // add axis labels
    histoData->GetXaxis()->SetTitle(xTitle);
    histoData->GetXaxis()->CenterTitle();
    histoData->GetXaxis()->SetTitleSize(0.055);
    histoData->GetXaxis()->SetTitleOffset(0.90);
    histoData->GetXaxis()->SetLabelOffset(0.010);
    histoData->GetXaxis()->SetLabelSize(0.05);
    histoData->GetYaxis()->SetTitle(yTitle);
    histoData->GetYaxis()->CenterTitle();
    histoData->GetYaxis()->SetTitleSize(0.055);
    histoData->GetYaxis()->SetTitleOffset(1.0);
    histoData->GetYaxis()->SetLabelSize(0.05);
    gStyle->SetTitleSize(0.08, "t");
    histoData->SetTitle(histName);
    can->Update();

    can->SaveAs(name + ".png");
}

void histogram(TH1D *histoData1, TH1D *histoData2, const TString histName, TCanvas *can, const TString xTitle, const TString yTitle, const TString name){

    histoData1->Sumw2();
    histoData1->Scale(10.0/histoData1->Integral());
    histoData1->SetLineWidth(2);
    histoData1->SetMarkerStyle(20);
    histoData1->SetMarkerSize(1.5);
    histoData1->SetMarkerColor(kRed);
    histoData1->SetLineColor(kBlack);
    // TF1 *f1 = new TF1(histName + "1", "[0]+[1]*TMath::Voigt(x-[2],[3],[4])");
    // f1->SetParameter(0, 200);
    // f1->SetParameter(1, histoData1->GetMaximum()/10.);
    // f1->SetParameter(2, histoData1->GetMean());
    // f1->SetParameter(3, .05);
    // f1->SetParameter(4, .05);
    // f1->SetParName(0, "Offset");
    // f1->SetParName(1, "Normalization");
    // f1->SetParName(2, "Mean");
    // f1->SetParName(3, "Sigma");
    // f1->SetParName(4, "gamma");
    // f1->SetLineColor(kBlue);
    // f1->SetLineWidth(4);
    // histoData1->Fit(f1, "ME");
    histoData1->Fit("gaus", "0");

    histoData2->Sumw2();
    histoData2->Scale(10.0/histoData2->Integral());
    histoData2->SetLineWidth(2);
    histoData2->SetMarkerStyle(21);
    histoData2->SetMarkerSize(1.5);
    histoData2->SetMarkerColor(kBlue);
    histoData2->SetLineColor(kBlack);
    // TF1 *f2 = new TF1(histName + "1", "[0]+[1]*TMath::Voigt(x-[2],[3],[4])");
    // f2->SetParameter(0, 200);
    // f2->SetParameter(1, histoData1->GetMaximum()/10.);
    // f2->SetParameter(2, histoData1->GetMean());
    // f2->SetParameter(3, .05);
    // f2->SetParameter(4, .05);
    // f2->SetParName(0, "Offset");
    // f2->SetParName(1, "Normalization");
    // f2->SetParName(2, "Mean");
    // f2->SetParName(3, "Sigma");
    // f2->SetParName(4, "gamma");
    // f2->SetLineColor(kBlue);
    // f2->SetLineWidth(4);
    // histoData2->Fit(f2, "ME");
    histoData2->Fit("gaus", "0");

    gStyle->SetLegendBorderSize(0);
    TLegend *leg = new TLegend(0.15,0.6,0.385,0.875);
    leg->SetTextSize(0.055);
    leg->AddEntry(histoData1, "Hydrogen","lep");
    leg->AddEntry(histoData2, "Deuterium","lep");
    // leg->AddEntry(f, "Voigt fit","l");

    Double_t max = histoData1->GetMaximum();
    if(histoData2->GetMaximum() > max) max = histoData2->GetMaximum();
    max*=1.2;
    histoData1->SetMaximum(max);
    histoData1->SetMinimum(0.);
    histoData2->SetMaximum(max);
    histoData2->SetMinimum(0.);
    
    histoData1->Draw("E1");
    histoData2->Draw("same E1");
    leg->Draw("same");

    // add axis labels
    histoData1->GetXaxis()->SetTitle(xTitle);
    histoData1->GetXaxis()->CenterTitle();
    histoData1->GetXaxis()->SetTitleSize(0.055);
    histoData1->GetXaxis()->SetTitleOffset(0.90);
    histoData1->GetXaxis()->SetLabelOffset(0.010);
    histoData1->GetXaxis()->SetLabelSize(0.05);
    histoData1->GetYaxis()->SetTitle(yTitle);
    histoData1->GetYaxis()->CenterTitle();
    histoData1->GetYaxis()->SetTitleSize(0.055);
    histoData1->GetYaxis()->SetTitleOffset(1.0);
    histoData1->GetYaxis()->SetLabelSize(0.05);
    gStyle->SetTitleSize(0.08, "t");
    histoData1->SetTitle(histName);
    can->Update();

    can->SaveAs(name + ".png");
}

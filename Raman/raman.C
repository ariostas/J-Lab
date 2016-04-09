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
#include <TGraphErrors.h>
#include "TSpectrum.h"

#endif

using namespace TMath;
using namespace std;

// Declare functions
void histogram(TH1D *histoData, const TString histName, TCanvas *can, const TString xTitle, const TString yTitle, const TString name, const TString type);
void histogram(TH1D*, TH1D*, const TString, TCanvas*, const TString, const TString, const TString);
TH1D* readData(TString, TString, Double_t, Double_t);
void findPeaks(TH1D*);

/*
 * MAIN FUNCTION
 */

void raman(){

    TH1::StatOverflows(kTRUE);
    
    cout << "\n\nStarting process...\n\n";

    TH1D *oxygen1 = readData("oxygen_0psi_25cm_mar15", "Oxygen", 25, 180.2);
    TH1D *oxygen2 = readData("oxygen_15psi_25cm_mar15", "Oxygen", 25, 0);
    TH1D *oxygen3 = readData("oxygen_30psi_25cm_mar15", "Oxygen", 25, 0);
    TH1D *oxygen4 = readData("oxygen_0psi_25cm_mar24", "Oxygen", 25, 0);

    TH1D *nitrogen1 = readData("nitrogen_15psi_25cm_mar17", "Nitrogen", 25, 0);
    TH1D *nitrogen2 = readData("nitrogen_30psi_5cm_mar17", "Nitrogen", 5, 0);
    TH1D *nitrogen3 = readData("test_117", "Nitrogen", 5, 0);

    TH1D *carbondio = readData("co2_0psi_5cm_mar24", "CO2", 5, 0);

    findPeaks(oxygen1);
    findPeaks(oxygen2);
    findPeaks(oxygen3);
    findPeaks(oxygen4);

    findPeaks(nitrogen1);
    findPeaks(nitrogen2);
    findPeaks(nitrogen3);

    findPeaks(carbondio);

    TCanvas *can = new TCanvas("Canvas", "Canvas", 1600, 900);
    histogram(oxygen1, "Oxygen 1 atm", can, "#Delta k (cm^{-1})", "Photon count", "oxygen1", "line");
    histogram(oxygen2, "Oxygen 2 atm", can, "#Delta k (cm^{-1})", "Photon count", "oxygen2", "line");
    histogram(oxygen3, "Oxygen 3 atm", can, "#Delta k (cm^{-1})", "Photon count", "oxygen3", "line");
    histogram(oxygen4, "Oxygen 1 atm (fixed slit)", can, "#Delta k (cm^{-1})", "Photon count", "oxygen4", "line");

    histogram(nitrogen1, "Nitrogen 1 atm", can, "#Delta k (cm^{-1})", "Photon count", "nitrogen1", "line");
    histogram(nitrogen2, "Nitrogen 2 atm", can, "#Delta k (cm^{-1})", "Photon count", "nitrogen2", "line");
    histogram(nitrogen3, "Nitrogen 1 atm (small slits)", can, "#Delta k (cm^{-1})", "Photon count", "nitrogen3", "line");

    histogram(carbondio, "CO2 1 atm", can, "#Delta k (cm^{-1})", "Photon count", "carbondio", "line");

}

TH1D* readData(TString dataset, TString name, Double_t speed, Double_t centralLine){

    TString fileName = "./Data/" + dataset + ".txt";

    Int_t nBins = 0;
    Double_t offset = 0;

    ifstream ifs(fileName); if(!ifs.is_open()){cout << "Error. File " << fileName << " not found. Exiting...\n"; return NULL;}

    Double_t entry, intensity;
    
    bool offsetDone = false;
    while(ifs >> entry >> intensity){

        if(!offsetDone){offset = entry; offsetDone = true;}
        nBins++;

    }

    ifs.close();

    Double_t minHist, maxHist;

    if(centralLine == 0){

        minHist = 0;
        maxHist = (nBins-0.5)*speed/120.;

    }
    else{

        minHist = -centralLine-0.5*speed/120.;
        maxHist = -centralLine+(nBins-0.5)*speed/120.;
    }

    TH1D *histo = new TH1D(dataset, dataset, nBins, minHist, maxHist);


    ifs.open(fileName); if(!ifs.is_open()){cout << "Error. File " << fileName << " not found. Exiting...\n"; return NULL;}
    
    while(ifs >> entry >> intensity){
        
        for(Int_t x = 0; x< intensity*1829; x++){
            histo->Fill(-centralLine+(nBins/2.-entry+offset)*speed/60);
        }

    }

    ifs.close();

    // TF1 *f;
    // if(peak2 == 0.){
    //     f = new TF1("fit", "[0]+[1]*TMath::Voigt(x-[2],[3],[4])");
    //     f->SetParameter(0, histo->GetMinimum());
    //     f->SetParameter(1, histo->GetMaximum()/25.);
    //     f->SetParameter(2, peak1 + shift);
    //     f->SetParameter(3, .02);
    //     f->SetParameter(4, .02);
    //     f->SetParName(0, "Offset");
    //     f->SetParName(1, "Normalization");
    //     f->SetParName(2, "Mean");
    //     f->SetParName(3, "Sigma");
    //     f->SetParName(4, "Gamma");
    // }
    // else{
    //     f = new TF1("fit", "[0]+[1]*TMath::Voigt(x-[2],[3],[4])+[5]*TMath::Voigt(x-[6],[7],[8])");
    //     f->SetParameter(0, histo->GetMinimum());
    //     f->SetParameter(1, histo->GetMaximum()/20.);
    //     f->SetParameter(2, peak1+shift-0.1);
    //     f->SetParameter(3, .02);
    //     f->SetParameter(4, .02);
    //     f->SetParameter(5, histo->GetMaximum()/20.);
    //     f->SetParameter(6, peak2+shift+0.1);
    //     f->SetParameter(7, .02);
    //     f->SetParameter(8, .02);
    //     f->SetParName(0, "Offset");
    //     f->SetParName(1, "Normalization 1");
    //     f->SetParName(2, "Mean 1");
    //     f->SetParName(3, "Sigma 1");
    //     f->SetParName(4, "Gamma 1");
    //     f->SetParName(5, "Normalization 2");
    //     f->SetParName(6, "Mean 2");
    //     f->SetParName(7, "Sigma 2");
    //     f->SetParName(8, "Gamma 2");
    // }

    // f->SetLineColor(kBlue);
    // f->SetLineWidth(4);
    // cout << "Fitting " << name << endl;
    // histo->Fit(f, "ME");

    return histo;

}

void findPeaks(TH1D* hist){

    TSpectrum *s = new TSpectrum(50);
    Int_t nfound = s->Search(hist,10,"",0.01);
    TH1 *hb = s->Background(hist,20,"same");
    // if (hb) c1->Update();

    TF1 *fline = new TF1("fline","pol1",0,1000);
    hist->Fit(fline,"qn");

    Int_t npeaks = 0;
    Float_t *xpeaks = s->GetPositionX();
    for (Int_t p=0;p<nfound;p++) {
        Double_t xp = xpeaks[p];
        Int_t bin = hist->GetXaxis()->FindBin(xp);
        Double_t yp = hist->GetBinContent(bin);
        if (yp-TMath::Sqrt(yp) < fline->Eval(xp)) continue;
        npeaks++;
        printf("%s: peak %i at position %4.2f\n",hist->GetTitle(), p, xp);
    }
    cout << endl;

}

/*
 * FUNCTION FOR SAVING ONE HISTOGRAM
 */

void histogram(TH1D *histoData, const TString histName, TCanvas *can, const TString xTitle, const TString yTitle, const TString name, const TString type){

    gPad->SetLogy( (type.Contains("log") ? 1 : 0) );
    gStyle->SetOptFit(kFALSE);
    gStyle->SetOptStat(kFALSE);
    gStyle->SetStatFormat("6.2g");
    gStyle->SetFitFormat("4.3f");

    if(!histoData) return;

    // TF1 *f = histoData->GetFunction("4");

    if(type.Contains("dot")){
        histoData->SetLineWidth(1);
        histoData->SetMarkerStyle(20);
        histoData->SetMarkerSize(0.6);
        histoData->SetMarkerColor(kRed);
        histoData->SetLineColor(kBlack);
        // f->SetLineColor(kBlue);
        // f->SetLineWidth(4);
    }

    else if(type.Contains("line")){
        histoData->SetLineWidth(1);
        histoData->SetLineColor(kBlue);
        // f->SetLineColor(kRed);
        // f->SetLineWidth(2);
    }
    else return;

    gStyle->SetLegendBorderSize(0);
    TLegend *leg = new TLegend(0.755,0.675,0.985,0.875);
    leg->SetTextSize(0.055);
    leg->AddEntry(histoData, "Data", (type.Contains("dot") ? "lep" : "l"));
    // leg->AddEntry(f, "Voigtian fits","l");
    
    histoData->Draw((type.Contains("dot") ? "E1" : ""));
    // leg->Draw("same");

    // add axis labels
    histoData->GetXaxis()->SetTitle(xTitle);
    histoData->GetXaxis()->CenterTitle();
    histoData->GetXaxis()->SetTitleSize(0.055);
    histoData->GetXaxis()->SetTitleOffset(0.84);
    histoData->GetXaxis()->SetLabelOffset(-0.001);
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

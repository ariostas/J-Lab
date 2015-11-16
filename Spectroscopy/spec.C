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

#endif

using namespace TMath;
using namespace std;

// Declare functions
void histogram(TH1D*, const TString, TCanvas*, const TString, const TString, const TString);
void histogram(TH1D*, TH1D*, const TString, TCanvas*, const TString, const TString, const TString);
TH1D* readData(TString, TString, TString, Double_t, Double_t, Double_t, Double_t);
void calibration();
Double_t getWavelength(Double_t);
Double_t getWavelengthError(Double_t, Double_t);
void computeWavelengths();

// Declare histogram vectors
vector<TH1D*> hCalibration, hSodium, hHydrogen, hDeuterium;

// Declare calibration functions
TF1 *calF=0, *calSigF=0;

/*
 * MAIN FUNCTION
 */

void spec(){

    TH1::StatOverflows(kTRUE);
    
    cout << "\n\nStarting process...\n\n";

    calibration();

    cout << getWavelength(5000) << " +- " << getWavelengthError(5000, 0) << endl;
    cout << getWavelengthError(6161,0) << endl;
    cout << getWavelengthError(6154,0) << endl;
    cout << getWavelengthError(5896,0) << endl;
    cout << getWavelengthError(5890,0) << endl;
    cout << getWavelengthError(5688,0) << endl;
    cout << getWavelengthError(5682,0) << endl;
    cout << getWavelengthError(4752,0) << endl;
    cout << getWavelengthError(4748,0) << endl;
    cout << getWavelengthError(4498,0) << endl;
    cout << getWavelengthError(4494,0) << endl;

    TString dataset, wavelength;
    Double_t peak1, peak2, intTime;

    ifstream ifs("calibrationData.txt"); if(!ifs.is_open()){cout << "Error. File " << "calibrationData.txt" << " not found. Exiting...\n"; return;}
    
    while(ifs >> dataset >> wavelength >> peak1 >> peak2 >> intTime){
        
        if(dataset.Contains("#")) continue;

        //hCalibration.push_back(readData(dataset, "Mercury", wavelength, peak1, peak2, intTime, 0));
        //hCalibration.push_back(readData(dataset, "Mercury", wavelength, peak1, peak2, intTime,
        //    getWavelength(hCalibration.at(hCalibration.size()-1)->GetFunction("fit")->GetParameter(2)) - hCalibration.at(hCalibration.size()-1)->GetFunction("fit")->GetParameter(2)));
                
    }

    ifs.close();
    ifs.open("hydrogenSplittingData.txt"); if(!ifs.is_open()){cout << "Error. File " << "hydrogenData.txt" << " not found. Exiting...\n"; return;}
    
    while(ifs >> dataset >> wavelength >> peak1 >> peak2 >> intTime){
        
        if(dataset.Contains("#")) continue;

        //hHydrogen.push_back(readData(dataset, "Hydrogen", wavelength, peak1, peak2, intTime, 0));
        //hHydrogen.push_back(readData(dataset, "Deuterium", wavelength, peak1, peak2, intTime,
        //    getWavelength(hHydrogen.at(hHydrogen.size()-1)->GetFunction("fit")->GetParameter(2)) - hHydrogen.at(hHydrogen.size()-1)->GetFunction("fit")->GetParameter(2)));
                
    }

    ifs.close();
    ifs.open("deuteriumData.txt"); if(!ifs.is_open()){cout << "Error. File " << "deuteriumData.txt" << " not found. Exiting...\n"; return;}
    
    while(ifs >> dataset >> wavelength >> peak1 >> peak2 >> intTime){
        
        if(dataset.Contains("#")) continue;

        //hDeuterium.push_back(readData(dataset, "Deuterium", wavelength, peak1, peak2, intTime, true));
                
    }

    ifs.close();
    ifs.open("sodiumData.txt"); if(!ifs.is_open()){cout << "Error. File " << "sodiumData.txt" << " not found. Exiting...\n"; return;}
    
    while(ifs >> dataset >> wavelength >> peak1 >> peak2 >> intTime){
        
        if(dataset.Contains("#")) continue;

        hSodium.push_back(readData(dataset, "Sodium", wavelength, peak1, peak2, intTime, 0));
        hSodium.push_back(readData(dataset, "Sodium", wavelength, peak1, peak2, intTime,
            getWavelength(hSodium.at(hSodium.size()-1)->GetFunction("fit")->GetParameter(2)) - hSodium.at(hSodium.size()-1)->GetFunction("fit")->GetParameter(2)));
                
    }

    ifs.close();
    

    TCanvas *can = new TCanvas("Histogram", "Histogram", 1600, 900);

    gStyle->SetOptStat(10);
    //gStyle->SetOptFit(1111);
    //gStyle->SetOptStat(kFALSE);
    gStyle->SetOptFit(kFALSE);


    for(UInt_t x = 0; x < hCalibration.size(); x++){
        //histogram(hCalibration.at(x), hCalibration.at(x)->GetTitle(), can, "", "", TString::Format("Calibration%1i", x));
    }
    for(UInt_t x = 0; x < hHydrogen.size(); x++){
        //histogram(hHydrogen.at(x), hHydrogen.at(x)->GetTitle(), can, "", "", TString::Format("Hydrogen%1i", x));
    }
    for(UInt_t x = 0; x < hDeuterium.size(); x++){
        //histogram(hDeuterium.at(x), hDeuterium.at(x)->GetTitle(), can, "", "", TString::Format("Deuterium%1i", x));
    }
    for(UInt_t x = 0; x < hSodium.size(); x++){
        histogram(hSodium.at(x), hSodium.at(x)->GetTitle(), can, "", "", TString::Format("Sodium%1i", x));
    }

     // TH1D *hDLines = new TH1D("DLines", "DLines", 2000, hSodium.at(1)->GetXaxis()->GetXmin()-1., hSodium.at(3)->GetXaxis()->GetXmax()+1.);
     // for(Int_t x =1; x <= hSodium.at(1)->GetNbinsX(); x++){
     //     for(Int_t i = 0; i < hSodium.at(1)->GetBinContent(x); i++){
     //         hDLines->Fill(hSodium.at(1)->GetBinCenter(x));
     //     }
     // }
     // for(Int_t x =1; x <= hSodium.at(3)->GetNbinsX(); x++){
     //     for(Int_t i = 0; i < hSodium.at(3)->GetBinContent(x); i++){
     //         hDLines->Fill(hSodium.at(3)->GetBinCenter(x));
     //     }
     // }

    // // TF1 *f = new TF1("Sum of histos", "[0]+[1]*TMath::Voigt(x-[2],[3],[4])+[5]*TMath::Gaus(x,[6],[7])+[8]*TMath::Gaus(x,[9],[10])+[11]*TMath::Gaus(x,[12],[13])");
    // // f->SetParameter(0, 0);
    // // f->SetParameter(1, 70);
    // // f->SetParameter(2, 5771.95);
    // // f->SetParameter(3, .02);
    // // f->SetParameter(4, .03);
    // // f->SetParameter(5, 60);
    // // f->SetParameter(6, 5771.75);
    // // f->SetParameter(7, .005);
    // // f->SetParameter(8, 30);
    // // f->SetParameter(9, 5772.05);
    // // f->SetParameter(10, .005);
    // // f->SetParameter(11, 30);
    // // f->SetParameter(12, 5772.25);
    // // f->SetParameter(13, .005);
    // // f->SetParName(0, "Offset");
    // // f->SetParName(1, "Normalization");
    // // f->SetParName(2, "Mean");
    // // f->SetParName(3, "Sigma");
    // // f->SetParName(4, "gamma");

    // hDLines->SetLineWidth(2);
    // hDLines->SetMarkerStyle(20);
    // hDLines->SetMarkerSize(1.5);
    // hDLines->SetMarkerColor(kRed);
    // hDLines->SetLineColor(kBlack);
    // hDLines->GetXaxis()->SetTitle("Wavelength (#AA)");
    // hDLines->GetXaxis()->CenterTitle();
    // hDLines->GetXaxis()->SetTitleSize(0.045);
    // hDLines->GetXaxis()->SetTitleOffset(1.00);
    // hDLines->GetXaxis()->SetLabelOffset(0.010);
    // hDLines->GetXaxis()->SetLabelSize(0.045);
    // hDLines->GetYaxis()->SetTitle("Photon count");
    // hDLines->GetYaxis()->CenterTitle();
    // hDLines->GetYaxis()->SetTitleSize(0.050);
    // hDLines->GetYaxis()->SetTitleOffset(0.95);
    // hDLines->GetYaxis()->SetLabelSize(0.045);
    // gStyle->SetTitleSize(0.070, "t");
    // hDLines->SetTitle("Sodium D-lines 5890.57 / 5896.52 #AA");

    // f->SetLineColor(kBlue);
    // f->SetLineWidth(4);
    // cout << "Fitting " << "Sum of histos" << endl;
    // histSum->Fit(f, "ME");

    //histogram(hDLines, hDLines->GetTitle(), can, "", "", "histSum");

    computeWavelengths();

}

TH1D* readData(TString dataset, TString name, TString wavel, Double_t peak1, Double_t peak2, Double_t intTime, Double_t shift){

    TString fileName = "./Data/" + dataset + ".txt";

    Double_t histoMin = 99999, histoMax = 0;

    ifstream ifs(fileName); if(!ifs.is_open()){cout << "Error. File " << fileName << " not found. Exiting...\n"; return NULL;}

    Double_t entry, wavelength, photons, stepsize;
    
    while(ifs >> entry >> wavelength >> photons){

        if(entry == 1) stepsize = wavelength;
        else if(entry == 2) stepsize = wavelength - stepsize;
        
        if(wavelength < histoMin) histoMin = wavelength;
        if(wavelength > histoMax) histoMax = wavelength;

    }

    ifs.close();

    TH1D *histo = new TH1D(name + wavel + TString::Format(" %1.2f", shift), name + wavel + TString::Format(" %1.2f", shift), Nint(entry), histoMin + shift - stepsize/2.0, histoMax + shift + stepsize/2.0);


    ifs.open(fileName); if(!ifs.is_open()){cout << "Error. File " << fileName << " not found. Exiting...\n"; return NULL;}
    
    while(ifs >> entry >> wavelength >> photons){
        
        for(Int_t x = 0; x< photons/(1000.0/intTime); x++){
            histo->Fill(wavelength + shift);
        }

    }

    ifs.close();

    TF1 *f;
    if(peak2 == 0.){
        f = new TF1("fit", "[0]+[1]*TMath::Voigt(x-[2],[3],[4])");
        f->SetParameter(0, histo->GetMinimum());
        f->SetParameter(1, histo->GetMaximum()/25.);
        f->SetParameter(2, peak1 + shift);
        f->SetParameter(3, .02);
        f->SetParameter(4, .02);
        f->SetParName(0, "Offset");
        f->SetParName(1, "Normalization");
        f->SetParName(2, "Mean");
        f->SetParName(3, "Sigma");
        f->SetParName(4, "Gamma");
    }
    else{
        f = new TF1("fit", "[0]+[1]*TMath::Voigt(x-[2],[3],[4])+[5]*TMath::Voigt(x-[6],[7],[8])");
        f->SetParameter(0, histo->GetMinimum());
        f->SetParameter(1, histo->GetMaximum()/20.);
        f->SetParameter(2, peak1+shift-0.1);
        f->SetParameter(3, .02);
        f->SetParameter(4, .02);
        f->SetParameter(5, histo->GetMaximum()/20.);
        f->SetParameter(6, peak2+shift+0.1);
        f->SetParameter(7, .02);
        f->SetParameter(8, .02);
        f->SetParName(0, "Offset");
        f->SetParName(1, "Normalization 1");
        f->SetParName(2, "Mean 1");
        f->SetParName(3, "Sigma 1");
        f->SetParName(4, "Gamma 1");
        f->SetParName(5, "Normalization 2");
        f->SetParName(6, "Mean 2");
        f->SetParName(7, "Sigma 2");
        f->SetParName(8, "Gamma 2");
    }

    histo->SetLineWidth(2);
    histo->SetMarkerStyle(20);
    histo->SetMarkerSize(1.5);
    histo->SetMarkerColor(kRed);
    histo->SetLineColor(kBlack);
    histo->GetXaxis()->SetTitle("Wavelength (#AA)");
    histo->GetXaxis()->CenterTitle();
    histo->GetXaxis()->SetTitleSize(0.045);
    histo->GetXaxis()->SetTitleOffset(1.00);
    histo->GetXaxis()->SetLabelOffset(0.010);
    histo->GetXaxis()->SetLabelSize(0.045);
    histo->GetYaxis()->SetTitle("Photon count");
    histo->GetYaxis()->CenterTitle();
    histo->GetYaxis()->SetTitleSize(0.050);
    histo->GetYaxis()->SetTitleOffset(0.95);
    histo->GetYaxis()->SetLabelSize(0.045);
    gStyle->SetTitleSize(0.070, "t");
    if(peak2 == 0.) histo->SetTitle(name + TString::Format(" %4.2f #AA", f->GetParameter(2)));
    else histo->SetTitle(name + TString::Format(" %4.2f / %4.2f #AA", f->GetParameter(2), f->GetParameter(6)));
    histo->SetTitle("Sodium 5890.522 / 5890.567 #AA");
    cout << shift << endl;

    f->SetLineColor(kBlue);
    f->SetLineWidth(4);
    cout << "Fitting " << name << endl;
    histo->Fit(f, "ME0");

    return histo;

}

void calibration(){

    Double_t wavelength[10], marker[10], wavelengthError[10], markerError[10];

    wavelength[0] = 3131.548; marker[0] = 3133.78;
    wavelength[1] = 3131.839; marker[1] = 3134.07;
    wavelength[2] = 3650.153; marker[2] = 3649.38;
    wavelength[3] = 4046.563; marker[3] = 4043.15;
    wavelength[4] = 4358.328; marker[4] = 4352.62;
    wavelength[5] = 5460.735; marker[5] = 5445.24;
    wavelength[6] = 5769.598; marker[6] = 5751.11;
    wavelength[7] = 5790.663; marker[7] = 5771.92;

    // wavelengthError[0] = 0.0068; markerError[0] = 5.80457e-2;
    // wavelengthError[1] = 0.0068; markerError[1] = 4.18374e-2;
    // wavelengthError[2] = 0.0068; markerError[2] = 3.65332e-2;
    // wavelengthError[3] = 0.0068; markerError[3] = 4.64119e-2;
    // wavelengthError[4] = 0.0068; markerError[4] = 2.41694e-2;
    // wavelengthError[5] = 0.0068; markerError[5] = 2.45160e-2;
    // wavelengthError[6] = 0.0068; markerError[6] = 3.25336e-2;
    // wavelengthError[7] = 0.0068; markerError[7] = 3.44601e-2;

    // wavelengthError[0] = 0.0068; markerError[0] = 5.72027e-2;
    // wavelengthError[1] = 0.0068; markerError[1] = 4.14216e-2;
    // wavelengthError[2] = 0.0068; markerError[2] = 3.68554e-2;
    // wavelengthError[3] = 0.0068; markerError[3] = 5.18892e-2;
    // wavelengthError[4] = 0.0068; markerError[4] = 2.84812e-2;
    // wavelengthError[5] = 0.0068; markerError[5] = 2.64957e-2;
    // wavelengthError[6] = 0.0068; markerError[6] = 3.22631e-2;
    // wavelengthError[7] = 0.0068; markerError[7] = 4.06117e-2;

    wavelengthError[0] = 0.0068; markerError[0] = 5.72027e-2+0.000945915;
    wavelengthError[1] = 0.0068; markerError[1] = 4.14216e-2+0.000622423;
    wavelengthError[2] = 0.0068; markerError[2] = 3.68554e-2+0.000184086;
    wavelengthError[3] = 0.0068; markerError[3] = 5.18892e-2+0.000365154;
    wavelengthError[4] = 0.0068; markerError[4] = 2.84812e-2+0.000200005;
    wavelengthError[5] = 0.0068; markerError[5] = 2.64957e-2+0.000222945;
    wavelengthError[6] = 0.0068; markerError[6] = 3.22631e-2+0.000696721;
    wavelengthError[7] = 0.0068; markerError[7] = 4.06117e-2+0.016829400;

    // From fits
    // wavelength[0] = 3650.153; marker[0] = 3649.60;
    // wavelength[1] = 4046.563; marker[1] = 4043.29;
    // wavelength[2] = 4358.328; marker[2] = 4352.77;
    // wavelength[3] = 5460.735; marker[3] = 5445.41;
    // wavelength[4] = 5769.598; marker[4] = 5751.18;
    // wavelength[5] = 5790.663; marker[5] = 5771.98;
    
    TGraph *graph = new TGraph(8, marker, wavelength);
    TGraphErrors *graphErrors = new TGraphErrors(8, marker, wavelength, markerError, wavelengthError);
    TF1 *f = new TF1("CalibFit", "[0]+[1]*x+[2]*x*x+[3]*x*x*x");
    f->SetParameter(0, 1.98759);
    f->SetParameter(1, 1.00606);
    f->SetParameter(2, -2.20401e-06);
    f->SetParameter(3, 9.30838e-11);
    f->SetParName(0, "a");
    f->SetParName(1, "b");
    f->SetParName(2, "c");
    f->SetParName(3, "d");
    //graph->Fit(f, "FME");
    graphErrors->Fit(f, "FE");
    
    TCanvas *c1 = new TCanvas("Calibration", "Calibration", 1600, 900);

    gStyle->SetOptStat(2210);
    gStyle->SetOptFit(1111);

    graphErrors->SetMarkerSize(2.0);
    graphErrors->SetTitle("Mercury calibration");
    graphErrors->SetMarkerColor(4);
    graphErrors->SetMarkerStyle(21);
    graphErrors->GetXaxis()->SetTitle("Measured wavelength");
    graphErrors->GetXaxis()->CenterTitle();
    graphErrors->GetXaxis()->SetTitleSize(0.045);
    graphErrors->GetXaxis()->SetTitleOffset(1.00);
    graphErrors->GetXaxis()->SetLabelOffset(0.010);
    graphErrors->GetXaxis()->SetLabelSize(0.045);
    graphErrors->GetYaxis()->SetTitle("True wavelength");
    graphErrors->GetYaxis()->CenterTitle();
    graphErrors->GetYaxis()->SetTitleSize(0.050);
    graphErrors->GetYaxis()->SetTitleOffset(0.95);
    graphErrors->GetYaxis()->SetLabelSize(0.045);
    gStyle->SetTitleSize(0.070, "t");
    graphErrors->Draw("ALP");

    c1->SaveAs("Plot_Calibration.png");

    calF = new TF1("Calibration", "[0]+[1]*x+[2]*x*x+[3]*x*x*x");
    calSigF = new TF1("Calibration plus sigma", "[0]+[1]*x+[2]*x*x+[3]*x*x*x");

    // calF->FixParameter(0, f->GetParameter(0));
    // calF->FixParameter(1, f->GetParameter(1));
    // calF->FixParameter(2, f->GetParameter(2));
    // calF->FixParameter(3, f->GetParameter(3));

    // calSigF->FixParameter(0, f->GetParameter(0) + f->GetParError(0));
    // calSigF->FixParameter(1, f->GetParameter(1) + f->GetParError(1));
    // calSigF->FixParameter(2, f->GetParameter(2) + f->GetParError(2));
    // calSigF->FixParameter(3, f->GetParameter(3) + f->GetParError(3));

    calF->FixParameter(0, -1.84462);
    calF->FixParameter(1, 9.93825e-01);
    calF->FixParameter(2, 2.22390e-06);
    calF->FixParameter(3, -9.29431e-11);

    calSigF->FixParameter(0, -1.84462 + 8.25177e-02);
    calSigF->FixParameter(1, 9.93825e-01 + 2.84972e-05);
    calSigF->FixParameter(2, 2.22390e-06 + 7.45643e-09);
    calSigF->FixParameter(3, -9.29431e-11 + 1.01904e-12);

    cout << "Par error " << f->GetParError(2) << endl;
    cout << "error 5000 " << getWavelengthError(5000, 0);
    
}

Double_t getWavelength(Double_t marker){

    if(!calF){
        cout << "Error: calibration function not set." << endl;
        return 0;
    }
    else{
        return calSigF->Eval(marker);
    }
}

Double_t getWavelengthError(Double_t marker, Double_t sigma){

    if(!calF || !calSigF){
        cout << "Error: calibration functions not set." << endl;
        return 0;
    }
    else{
        return calSigF->Eval(marker + sigma) - calF->Eval(marker);
    }
}

void computeWavelengths(){

    Double_t echarge = 1.60217662e-19, hbar = 1.0545718e-34, h = 6.62607004e-34, c = 299792458, emass = 9.10938291e-31, pi = 3.14159265359, ep0 = 8.854187817e-12, ev = 1.602176565e-19;
    Double_t n = 3.0, z = 11.0, g = 2.0023193043617;
    Double_t bradius = 4.0*pi*ep0*hbar*hbar/(emass*echarge*echarge);
    Double_t alpha = echarge*echarge/(4.0*pi*ep0*hbar*c);


    //Double_t e3 = 11.0*11.0*echarge*echarge/(2.0*bradius*3.0*3.0);
    Double_t e3 = -11.0*11.0*emass*echarge*echarge*echarge*echarge/(32.0*pi*pi*ep0*ep0*hbar*hbar*3.0*3.0);
    //Double_t e3 = -z*z*alpha*alpha*emass*c*c/(n*n);

    //Double_t emax1 = e3*(1+11.0*11.0*alpha*alpha/(3.0)*(1.0/(1.0)-3.0/(4.0*3.0)));
    //Double_t emax1 = e3+e3*e3/(2.0*emass*c*c)*(3.0-4.0*n/(1.5))-n*e3*e3*g/(2.0*emass*c*c*1.5);
    //Double_t emin1 = e3+2.0*n*e3*e3/(emass*c*c);

    // Double_t emax2 = e3+e3*e3/(2.0*emass*c*c)*(3.0-4.0*n/(1.5))+n*e3*e3*g/(2.0*emass*c*c*3.0);
    // Double_t emin2 = e3+2.0*n*e3*e3/(emass*c*c);


    Double_t emax1 = e3-e3*e3/(2.0*emass*c*c)*(4.0*n/(1.5)-3.0)-e3*e3*n/(emass*c*c)*2.0/3.0;
    Double_t emin1 = e3-e3*e3/(2.0*emass*c*c)*(4.0*n/(1.5)-3.0)+2.0*n*e3*e3/(emass*c*c);

    Double_t emax2 = e3-e3*e3/(2.0*emass*c*c)*(4.0*n/(1.5)-3.0)+e3*e3*n/(emass*c*c)*1.0/3.0;
    Double_t emin2 = e3-e3*e3/(2.0*emass*c*c)*(4.0*n/(1.5)-3.0)+2.0*n*e3*e3/(emass*c*c);

    Double_t wavelength1 = c*h/(emax1-emin1), wavelength2 = c*h/(emax2-emin2);

    Double_t diff = 1.0/wavelength2 - 1.0/wavelength1;

    cout << "E1 is " << emax1/ev << " ev" << endl;
    cout << "E2 is " << emin1/ev << " ev" << endl;
    cout << "Wavelength is " << wavelength1*1.0e10 << " Angstroms." << endl;
    cout << "Wavelength is " << wavelength2*1.0e10 << " Angstroms." << endl;
    cout << "Diff Wavelength is " << wavelength1*1.0e10 - wavelength2*1.0e10 << " Angstroms." << endl;
    cout << "wavenumber diff is " << diff*100.0 << " /cm" << endl;
}

/*
 * FUNCTION FOR SAVING ONE HISTOGRAM
 */

void histogram(TH1D *histoData, const TString histName, TCanvas *can, const TString xTitle, const TString yTitle, const TString name){

    // histoData->SetLineWidth(2);
    // histoData->SetMarkerStyle(20);
    // histoData->SetMarkerSize(1.5);
    // histoData->SetMarkerColor(kRed);
    // histoData->SetLineColor(kBlack);
    // TF1 *f = new TF1(histName, "[0]+[1]*TMath::Voigt(x-[2],[3],[4])");
    // f->SetParameter(0, 200);
    // f->SetParameter(1, histoData->GetMaximum()/10.);
    // f->SetParameter(2, histoData->GetMean());
    // f->SetParameter(3, .05);
    // f->SetParameter(4, .05);
    // f->SetParName(0, "Offset");
    // f->SetParName(1, "Normalization");
    // f->SetParName(2, "Mean");
    // f->SetParName(3, "Sigma");
    // f->SetParName(4, "gamma");
    // f->SetLineColor(kBlue);
    // f->SetLineWidth(4);
    // histoData->Fit(f, "ME");
    TF1 *f = histoData->GetFunction("fit");

    gStyle->SetLegendBorderSize(0);
    TLegend *leg = new TLegend(0.655,0.675,0.885,0.875);
    leg->SetTextSize(0.055);
    leg->AddEntry(histoData, "Data","lep");
    //leg->AddEntry(f, "Double Voigtian fit","l");
    
    histoData->Draw("E1");
    leg->Draw("same");

    // // add axis labels
    // histoData->GetXaxis()->SetTitle(xTitle);
    // histoData->GetXaxis()->CenterTitle();
    // histoData->GetXaxis()->SetTitleSize(0.055);
    // histoData->GetXaxis()->SetTitleOffset(0.90);
    // histoData->GetXaxis()->SetLabelOffset(0.010);
    // histoData->GetXaxis()->SetLabelSize(0.05);
    // histoData->GetYaxis()->SetTitle(yTitle);
    // histoData->GetYaxis()->CenterTitle();
    // histoData->GetYaxis()->SetTitleSize(0.055);
    // histoData->GetYaxis()->SetTitleOffset(1.0);
    // histoData->GetYaxis()->SetLabelSize(0.05);
    // gStyle->SetTitleSize(0.08, "t");
    // histoData->SetTitle(histName);
    // can->Update();

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

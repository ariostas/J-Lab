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
#include <TF2.h>
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
#include "TTimer.h"

#endif

using namespace TMath;
using namespace std;

// Declare functions
void histogram(TH1D*, const TString, TCanvas*, const TString, const TString, const TString, const TString);
void histogram(TH1D*, TH1D*, const TString, TCanvas*, const TString, const TString, const TString, const TString);
void graph(TGraph*, const TString, TCanvas*, const TString, const TString, const TString);
void plot(TF1*, const TString, TCanvas*, const TString, const TString, const TString, const TString type = "");
TH1D* readFile(TString, TF1*);
TGraph* calibration();
TH1D* produceMonteCarlo(TF1*, TString);
TH1D* readMonteCarlo(TF1*, TString);

UInt_t h = 0;
TF1 *fcal;

/*
 * MAIN FUNCTION
 */

void muons(){

    TH1::StatOverflows(kTRUE);
    gStyle->SetOptFit(1111);
    gStyle->SetOptStat(10);
    // TCanvas *can = new TCanvas("Canvas", "Canvas", 1600, 900);

    TGraph *gcal = calibration();

    TH1D *lifetime = readFile("lifetime_data2", fcal);

    // TF1 *fLife = new TF1("lifetime", "[0]+[1]*TMath::Exp(-x/[2])");
    // fLife->SetParameter(0, 0);
    // fLife->SetParameter(1, lifetime->GetMaximum()*Exp(1));
    // fLife->SetParameter(2, 2.197);
    // fLife->SetParName(0, "Constant");
    // fLife->SetParName(1, "Normalization");
    // fLife->SetParName(2, "Mean lifetime");
    // lifetime->Fit(fLife, "ME");

    TF1 *fLife = new TF1("lifetime", "[0]+[1]*(0.44*TMath::Exp(-x/[2])+0.56*TMath::Exp(-x/[3]))");
    fLife->SetParameter(0, 0);
    fLife->SetParameter(1, lifetime->GetMaximum()*Exp(1));
    fLife->FixParameter(2, 2.023);
    fLife->SetParameter(3, 2.197);
    lifetime->Fit(fLife, "ME");

    // histogram(lifetime, "Data", can, "Time (#mus)", "Counts", "lifetime_dot", "dot");
    // histogram(lifetime, "Data", can, "Time (#mus)", "Counts", "lifetime_line", "line");
    // histogram(lifetime, "Data", can, "Time (#mus)", "Counts", "lifetime_dot_log", "dot log");
    // histogram(lifetime, "Data", can, "Time (#mus)", "Counts", "lifetime_line_log", "line log");
    //graph(gcal, "speed", can, "times", "distances", "calibration");

    cout << "The mean life of muon is " << fLife->GetParameter(2) << " +- " << fLife->GetParError(2) << " microseconds" << endl;

    Double_t gf = 1.16637e-11, pi = Pi(), hbar = 6.582119e-22, muT = 2.166, muTError = 0.035, muMass, muMassError;

    muMass = Power(192.*pi*pi*pi*hbar/(gf*gf*muT*1e-6), 1./5.);

    muMassError = muTError*muMass/(5*muT);

    cout << "The mass of muon is " << muMass << " +- " << muMassError << " MeV" << endl;

    Double_t totalCounts = lifetime->Integral(), totalTime = 239660, detectionRate = totalCounts/totalTime;

    cout << "The stopping rate of muons was " << detectionRate << " Hz" << endl;

    // TH1D *mcHist = produceMonteCarlo(fcal, "MonteCarloData");
    // histogram(mcHist, "meas 1", can, "Time (#mus)", "Counts", "lifetimeMC");

    // TH1D *mcHist = readMonteCarlo(fcal, "MonteCarloData1");
    // histogram(mcHist, "Monte Carlo", can, "Time (#mus)", "Counts", "lifetimeMC_dot", "dot");
    // histogram(mcHist, "Monte Carlo", can, "Time (#mus)", "Counts", "lifetimeMC_line", "line");
    // histogram(mcHist, "Monte Carlo", can, "Time (#mus)", "Counts", "lifetimeMC_dot_log", "dot log");
    // histogram(mcHist, "Monte Carlo", can, "Time (#mus)", "Counts", "lifetimeMC_line_log", "line log");

    // TF1 *lognormal = new TF1("lognormal", "[0]*TMath::LogNormal(TMath::Abs([1]*x-[2]), [3])", 0, 2048);
    // lognormal->SetParameter(0, 20);
    // lognormal->SetParameter(1, 1./700.);
    // lognormal->SetParameter(2, -10/25.);
    // lognormal->SetParameter(3, 0.20);

    // can->Clear();
    // gPad->SetLogy(0);
    // lognormal->Draw();

}

// TGraphErrors* readFile(TString inputfile, TF1* cal){

//     ifstream ifs("./Data/" + inputfile + ".txt"); if(!ifs.is_open()){cout << "Error. File " << inputfile << " not found. Exiting...\n"; return NULL;}

//     Double_t bins[2048], counts[2048], binError[2048], countError[2048];
//     Int_t bin = 0, count = 0;
    
//     while(ifs >> count){
        
//         bins[bin] = cal->Eval(bin);
//         binError[bin] = (cal->Eval(bin+1)-cal->Eval(bin-1))/2.;
//         counts[bin] = count;
//         countError[bin] = Sqrt(count);
//         bin++;
                
//     }

//     ifs.close();

//     TGraphErrors *p = new TGraphErrors(2048, bins, counts, binError, countError);

//     //TF1 *f = new TF1(inputfile, "[0]*TMath::Landau(x,[1],[2])+[3]");
//     // TF1 *f = new TF1(inputfile, "[0]*TMath::Gaus(x,[1],[2])+[3]");
//     // f->SetParameter(0, p->GetMaximum());
//     // f->SetParameter(1, p->GetMean());
//     // f->SetParameter(2, p->GetStdDev());
//     // f->SetParameter(3, p->GetMinimum());
//     // f->SetParLimits(0, 4, 30);
//     // f->SetParLimits(1, 300, 600);
//     // f->SetParLimits(2, 20, 50);
//     // f->SetParLimits(3, -1, 1);
//     //p->Fit(f, "ME");

//     return p;

// }

TH1D* readFile(TString inputfile, TF1* cal){

    ifstream ifs("./Data/" + inputfile + ".txt"); if(!ifs.is_open()){cout << "Error. File " << inputfile << " not found. Exiting...\n"; return NULL;}

    TH1D *p = new TH1D(inputfile, inputfile, 1024, cal->Eval(-0.5), cal->Eval(2047.5));

    Int_t bin = 0, count = 0;
    
    while(ifs >> count){
        
        if(bin < 122){bin++; continue;}
        //if(cal->Eval(bin) > 15) continue;
        for(Int_t x = 0; x < count; x++){
            p->Fill(cal->Eval(bin));
        }
        bin++;
                
    }

    ifs.close();

    //TF1 *f = new TF1(inputfile, "[0]*TMath::Landau(x,[1],[2])+[3]");
    TF1 *f = new TF1(inputfile, "[0]*TMath::Gaus(x,[1],[2])+[3]");
    f->SetParameter(0, p->GetMaximum());
    f->SetParameter(1, p->GetMean());
    f->SetParameter(2, p->GetStdDev());
    f->SetParameter(3, p->GetMinimum());
    f->SetParLimits(0, 4, 30);
    f->SetParLimits(1, 300, 600);
    f->SetParLimits(2, 20, 50);
    f->SetParLimits(3, -1, 1);
    //p->Fit(f, "ME");

    return p;

}

TGraph* calibration(){

    fcal = new TF1("calibration", "[0]+[1]*x");

    Double_t calTimes[20], calBins[20];

    calBins[0] = 131; calTimes[0] = 1.28;
    calBins[1] = 262; calTimes[1] = 2.56;
    calBins[2] = 393; calTimes[2] = 3.84;
    calBins[3] = 524; calTimes[3] = 5.12;
    calBins[4] = 655; calTimes[4] = 6.40;
    calBins[5] = 786; calTimes[5] = 7.68;
    calBins[6] = 917; calTimes[6] = 8.96;
    calBins[7] = 1048; calTimes[7] = 10.24;
    calBins[8] = 1179; calTimes[8] = 11.52;
    calBins[9] = 1310; calTimes[9] = 12.8;
    calBins[10] = 1441; calTimes[10] = 14.08;
    calBins[11] = 1572; calTimes[11] = 15.36;
    calBins[12] = 1703; calTimes[12] = 16.64;
    calBins[13] = 1834; calTimes[13] = 17.92;
    calBins[14] = 1965; calTimes[14] = 19.20;

    TGraph *graphCal = new TGraph(15, calBins, calTimes);
    graphCal->Fit(fcal, "ME");

    return graphCal;

}

TH1D* produceMonteCarlo(TF1* cal, TString outputFile){

    ofstream outFile;
    outFile.open(outputFile + ".txt");
    if(!outFile.is_open()) return NULL;

    TH1D *histo = new TH1D("monte carlo", "monte carlo", 256, cal->Eval(-0.5), cal->Eval(2047.5));

    Double_t expectedRate = 0.0207, integrationTime = 239660, diameterConsidered = 500, diameterDetector = 27.94, height = 30.48;
    Double_t stoppingRate = 2, meanLifePos = 2.197, meanLifeNeg = 2.043, muonMass = 105.6583715, timeOffset = cal->Eval(122), minDecayEnergy = 10;
    expectedRate *= Pi()*diameterConsidered*diameterConsidered/4.;

    TF1 *constant = new TF1("constant", "1", 0, 10000);
    TF1 *phiDist = new TF1("phi", "TMath::Cos(x)*TMath::Cos(x)", 0, Pi()/2.);
    TF1 *radiusDist = new TF1("radius", "TMath::Abs(x-[0])", 0, diameterConsidered);
    radiusDist->FixParameter(0, diameterConsidered/2.);
    TF1 *decayTimePosDist = new TF1("decay time", "TMath::Exp(-x/[0])", 0, 25);
    decayTimePosDist->FixParameter(0, meanLifePos);
    TF1 *decayTimeNegDist = new TF1("decay time", "TMath::Exp(-x/[0])", 0, 25);
    decayTimeNegDist->FixParameter(0, meanLifeNeg);
    TF1 *decayEnergyDist = new TF1("decay energy", "3*x*x-2*x*x*x", 0, 1);

    TF1 *energyDist = new TF1("energy", "[0]*TMath::LogNormal(TMath::Abs([1]*x-[2]), [3])", 10, 10000);
    energyDist->SetParameter(0, 1);
    energyDist->SetParameter(1, 1./8000.);
    energyDist->SetParameter(2, 1./500.);
    energyDist->SetParameter(3, 1.5);

    Double_t energies[10], fluxes[10];

    energies[0] = 10; fluxes[0] = 1;
    energies[1] = 100; fluxes[1] = 4;
    energies[2] = 400; fluxes[2] = 6;
    energies[3] = 1000; fluxes[3] = 5;
    energies[4] = 10000; fluxes[4] = 0.7;

    TGraph *enDist = new TGraph(5, energies, fluxes);
    enDist->Fit(energyDist);

    // TCanvas *c2 = new TCanvas("Plots", "Plots", 1600, 900);
    // gPad->SetLogx();
    // gPad->SetLogy();
    // gStyle->SetTitleSize(0.08, "t");
    // energyDist->SetParameter(0, energyDist->GetParameter(0)*1e-6);
    // plot(phiDist, "#phi distribution", c2, "#phi angle", "dP/d#phi", "phiDist");
    // plot(energyDist, "Energy distribution", c2, "Energy (MeV)", "Muon flux (cm^{-2}MeV^{-1}s^{-1})", "enDist");
    // plot(decayTimePosDist, "Decay time distribution", c2, "Time (#mus)", "dP/dt", "timeDist", "log");
    // plot(decayEnergyDist, "Decay energy distribution", c2, "Energy (#times100 MeV)", "dP/dE", "decEnDist");

    Double_t nCounts = 0;

    for(Int_t x = 0; x < expectedRate*integrationTime; x++){

        if(x%100000 == 0) cout << "Time is " << x/expectedRate << " s" << endl;

        Double_t muonType = constant->GetRandom(0, 1);
        Double_t theta = constant->GetRandom(0, Pi()/4.);
        Double_t distance = radiusDist->GetRandom(0, diameterConsidered);
        Double_t energy = energyDist->GetRandom(10, 10000);
        Double_t phi = phiDist->GetRandom(0, Pi()/2.);
        Double_t decayTime = (muonType <= 0.56 ? decayTimePosDist->GetRandom(0, 25) : decayTimeNegDist->GetRandom(0, 25));
        Double_t decayEnergy = decayEnergyDist->GetRandom(0, 1)*muonMass/2.;

        if((distance <= diameterConsidered/2. && diameterConsidered/2. - distance <= diameterDetector/2.) || 
            (distance >= diameterConsidered/2. && distance - diameterConsidered/2.<= diameterDetector/2.)){

            Double_t d, pathPlane, pathLength;

            distance = distance - (diameterConsidered/2. - diameterDetector/2.);

            if(distance >= diameterDetector/2.){
                d = distance - diameterDetector/2.;
                pathPlane = -d*Cos(theta)+Sqrt(d*d*Cos(theta)*Cos(theta)+diameterDetector*diameterDetector/4.-d*d);
            }
            else{
                d = diameterDetector/2. - distance;
                pathPlane = d*Cos(theta)+Sqrt(d*d*Cos(theta)*Cos(theta)+diameterDetector*diameterDetector/4.-d*d);
            }

            Double_t depth;
            if(phi != 0) depth = pathPlane/Tan(phi);

            if(phi == 0) pathLength = height;
            else if(depth <= height) pathLength = pathPlane/Sin(phi);
            else{
                pathLength = pathPlane/Sin(phi) - (pathPlane/Tan(phi)-height)/Cos(phi);
            }

            Double_t remainingEnergy = energy - pathLength*stoppingRate;

            if(remainingEnergy <= 0 && decayEnergy >= minDecayEnergy){
                nCounts++;
                histo->Fill(decayTime + timeOffset);
                outFile << muonType << " " << theta << " " << distance << " " << energy << " " << phi << " " << decayTime << " " << decayEnergy << endl;
            }

        }

        else if(distance >= diameterConsidered/2. && distance - diameterConsidered/2. >= diameterDetector/2.) continue;

        else{

            if(theta >= ASin((diameterDetector/2.)/(diameterConsidered/2. - distance))) continue;
            if(phi == 0) continue;

            Double_t entryDistance = diameterConsidered/2. - distance;
            Double_t pathConsidered = entryDistance*Cos(theta)-Sqrt(entryDistance*entryDistance*Cos(theta)*Cos(theta)+diameterDetector*diameterDetector/4.-entryDistance*entryDistance);

            Double_t effectiveHeight =  height - pathConsidered/Tan(phi);

            if(effectiveHeight <= 0) continue;

            Double_t effectiveTheta = Pi() - ACos((pathConsidered*pathConsidered+diameterDetector*diameterDetector/4.-entryDistance*entryDistance)/(2.*pathConsidered*diameterDetector/2.));

            if(effectiveTheta < 0 || effectiveTheta > Pi()/2.) cout << "Error in theta. Check code. " << effectiveTheta << endl;

            Double_t d, pathPlane, pathLength;

            d = diameterDetector/2.;
            pathPlane = d*Cos(effectiveTheta)+Sqrt(d*d*Cos(effectiveTheta)*Cos(effectiveTheta)+diameterDetector*diameterDetector/4.-d*d);

            Double_t depth;
            if(phi != 0) depth = pathPlane/Tan(phi);

            if(phi == 0) pathLength = effectiveHeight;
            else if(depth <= effectiveHeight) pathLength = pathPlane/Sin(phi);
            else{
                pathLength = pathPlane/Sin(phi) - (pathPlane/Tan(phi)-effectiveHeight)/Cos(phi);
            }

            Double_t remainingEnergy = energy - pathLength*stoppingRate;

            if(remainingEnergy <= 0 && decayEnergy >= minDecayEnergy){
                nCounts++;
                histo->Fill(decayTime + timeOffset);
                outFile << muonType << " " << theta << " " << distance << " " << energy << " " << phi << " " << decayTime << " " << decayEnergy << endl;
            }
        }

    }

    outFile.close();

    cout << "Frequency from Monte Carlo is " << nCounts/integrationTime << " +- " << Sqrt(nCounts)/integrationTime << " Hz" << endl;

    TF1 *f = new TF1("lifetime", "[0]+[1]*TMath::Exp(-x/[2])");
    f->SetParameter(0, 0);
    f->SetParameter(1, histo->GetMaximum()*Exp(1));
    f->SetParameter(2, 2.25);
    histo->Fit(f, "ME");

    return histo;

}

TH1D* readMonteCarlo(TF1* cal, TString inputfile){

    TH1D *histo = new TH1D(inputfile, inputfile, 256, cal->Eval(-0.5), cal->Eval(2047.5));

    ifstream ifs("./Data/" + inputfile + ".txt"); if(!ifs.is_open()){cout << "Error. File " << inputfile << " not found. Exiting...\n"; return NULL;}

    Double_t timeOffset = cal->Eval(119.5);
    Double_t muonType, theta, distance, energy, phi, decayTime, decayEnergy;
    
    while(ifs >> muonType >> theta >> distance >> energy >> phi >> decayTime >> decayEnergy){
        
        histo->Fill(decayTime + timeOffset);
        
    }
    cout << cal->Eval(122) << "  " << (cal->Eval(2047.5) - cal->Eval(-0.5))/1024. << endl;

    TF1 *f = new TF1("lifetime", "[0]+[1]*TMath::Exp(-x/[2])");
    f->SetParameter(0, 0);
    f->SetParameter(1, histo->GetMaximum()*Exp(1));
    f->SetParameter(2, 2.19);
    f->SetParName(0, "Constant");
    f->SetParName(1, "Normalization");
    f->SetParName(2, "Mean lifetime");
    histo->Fit(f, "ME");

    // TF1 *f = new TF1("lifetime", "[0]+[1]*(0.44*TMath::Exp(-x/[2])+0.56*TMath::Exp(-x/[3]))");
    // f->SetParameter(0, 0);
    // f->SetParameter(1, histo->GetMaximum()*Exp(1));
    // f->FixParameter(2, 2.043);
    // f->SetParameter(3, 2.197);
    // histo->Fit(f, "ME");

    return histo;
}


// Double_t monteCarloSpeed(){

//     Double_t expectedRate = 0.0207, integrationTime = 100000, diameter = 27.94, height = 30.48, stoppingRate = 5, meanEnergy = 50, energySpread = 2;
//     expectedRate *= Pi()*diameter*diameter/4.;

//     TF1 *constant = new TF1("constant", "1", 0, 10000);
//     TF1 *phiDist = new TF1("phi", "TMath::Cos(x)*TMath::Cos(x)", 0, Pi()/2.);
//     //TF1 *energyDist = new TF1("energy", "TMath::Gaus(x, [0], [1])", meanEnergy - 5*energySpread, meanEnergy + 5*energySpread);
//     //energyDist->FixParameter(0, meanEnergy);
//     //energyDist->FixParameter(1, energySpread);

//     Double_t nCounts = 0;

//     for(Int_t x = 0; x < expectedRate*integrationTime; x++){

//         Double_t theta = constant->GetRandom(0, Pi()/4.);
//         Double_t distance = constant->GetRandom(0, diameter);
//         //Double_t energy = energyDist->GetRandom();
//         Double_t energy = constant->GetRandom(10, 10000);
//         Double_t phi = phiDist->GetRandom(0, Pi()/2.);

//         Double_t d, pathPlane, pathLength;

//         if(distance >= diameter/2.){
//             d = distance - diameter/2.;
//             pathPlane = -d*Cos(theta)+Sqrt(d*d*Cos(theta)*Cos(theta)+diameter*diameter/4.-d*d);
//         }
//         else{
//             d = diameter/2. - distance;
//             pathPlane = d*Cos(theta)+Sqrt(d*d*Cos(theta)*Cos(theta)+diameter*diameter/4.-d*d);
//         }

//         Double_t depth;
//         if(phi != 0) depth = pathPlane/Tan(phi);

//         if(phi == 0) pathLength = height;
//         else if(depth <= height) pathLength = pathPlane/Sin(phi);
//         else{
//             pathLength = pathPlane/Sin(phi) - (pathPlane/Tan(phi)-height)/Cos(phi);
//         }

//         Double_t remainingEnergy = energy - pathLength*stoppingRate;

//         if(remainingEnergy <= 0) nCounts++;

//     }

//     return nCounts/integrationTime;


// }


void graph(TGraph *graph, const TString graphName, TCanvas *can, const TString xTitle, const TString yTitle, const TString name){

    //gStyle->SetOptStat(2210);
    gStyle->SetOptFit(1111);
    gStyle->SetOptStat(kFALSE);
    //gStyle->SetOptFit(1100);

    if(!graph){
        cout << "Error: Graph \"" << graphName << "\" not defined" << endl;
        return;
    }

    graph->SetLineWidth(2);
    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(1.5);
    graph->SetMarkerColor(kRed);
    graph->SetLineColor(kBlack);
    // TF1 *f = new TF1(graphName, "[0]+[1]*TMath::Power(TMath::BesselJ1([2]*x-[3])*TMath::BesselJ1([2]*x-[3])/(([2]*x-[3])*([2]*x-[3])),0.5)");
    TF1 *f = new TF1(graphName, "[0]+[1]/TMath::Log([4]*TMath::BesselJ1([2]*x-[3])*TMath::BesselJ1([2]*x-[3])/(([2]*x-[3])*([2]*x-[3])) + [5])");
    f->SetParameter(0, 250);
    f->SetParameter(1, 100);
    f->SetParameter(2, 0.395);
    f->SetParameter(3, 1);
    f->SetParameter(4, 10000);
    f->SetParameter(5, 100);
    f->SetParName(0, "Offset");
    f->SetParName(1, "Peak intensity");
    f->SetParName(2, "Aperture");
    f->SetParName(3, "Center value");
    f->SetParName(4, "Amp");
    f->SetParLimits(0, 100, 300);
    f->SetParLimits(1, 10, 100);
    f->SetParLimits(2, 0.1, 1);
    f->SetParLimits(3, -0.2, 2);
    f->SetParLimits(4, 50, 10000);
    f->SetParLimits(5, 1, 10000);
    f->SetLineColor(kBlue);
    f->SetLineWidth(4);
    //graph->Fit(f, "ME");

    gStyle->SetLegendBorderSize(0);
    TLegend *leg = new TLegend(0.655,0.675,0.885,0.875);
    leg->SetTextSize(0.055);
    leg->AddEntry(graph, "Data","lep");
    //leg->AddEntry(f, "Airy disk fit","l");
    
    graph->Draw("ap");
    //leg->Draw("same");

    // add axis labels
    graph->GetXaxis()->SetTitle(xTitle);
    graph->GetXaxis()->CenterTitle();
    graph->GetXaxis()->SetTitleSize(0.055);
    graph->GetXaxis()->SetTitleOffset(0.82);
    //graph->GetXaxis()->SetLabelOffset(0.010);
    graph->GetXaxis()->SetLabelSize(0.05);
    graph->GetYaxis()->SetTitle(yTitle);
    graph->GetYaxis()->CenterTitle();
    graph->GetYaxis()->SetTitleSize(0.055);
    graph->GetYaxis()->SetTitleOffset(0.9);
    graph->GetYaxis()->SetLabelSize(0.05);
    gStyle->SetTitleSize(0.08, "t");
    graph->SetTitle(graphName);
    can->Update();

    can->SaveAs(name + ".png");

    //can->Clear();
}

void plot(TF1 *function, const TString fName, TCanvas *can, const TString xTitle, const TString yTitle, const TString name, const TString type){

    if(!function){
        cout << "Error: Graph \"" << fName << "\" not defined" << endl;
        return;
    }

    function->SetLineWidth(2);
    function->SetLineColor(kBlack);
    
    //can->Clear();
    function->Draw();

    // add axis labels
    function->GetXaxis()->SetTitle(xTitle);
    function->GetXaxis()->CenterTitle();
    function->GetXaxis()->SetTitleSize(0.055);
    function->GetXaxis()->SetTitleOffset(0.82);
    //graph->GetXaxis()->SetLabelOffset(0.010);
    function->GetXaxis()->SetLabelSize(0.05);
    function->GetYaxis()->SetTitle(yTitle);
    function->GetYaxis()->CenterTitle();
    function->GetYaxis()->SetTitleSize(0.055);
    function->GetYaxis()->SetTitleOffset(0.8);
    function->GetYaxis()->SetLabelSize(0.05);
    gStyle->SetTitleSize(0.08, "t");
    function->SetTitle(fName);
    can->Update();

    can->SaveAs(name + ".png");

    //can->Clear();
}

/*
 * FUNCTION FOR SAVING ONE HISTOGRAM
 */

void histogram(TH1D *histoData, const TString histName, TCanvas *can, const TString xTitle, const TString yTitle, const TString name, const TString type){

    gPad->SetLogy( (type.Contains("log") ? 1 : 0) );

    if(!histoData) return;

    TF1 *f = histoData->GetFunction("lifetime");

    if(type.Contains("dot")){
        histoData->SetLineWidth(1);
        histoData->SetMarkerStyle(20);
        histoData->SetMarkerSize((histName.Contains("Monte") ? 1.0 : 0.9));
        histoData->SetMarkerColor(kRed);
        histoData->SetLineColor(kBlack);
        f->SetLineColor(kBlue);
        f->SetLineWidth(4);
    }

    else if(type.Contains("line")){
        histoData->SetLineWidth(1);
        histoData->SetLineColor(kBlue);
        f->SetLineColor(kRed);
        f->SetLineWidth(2);
    }
    else return;

    gStyle->SetLegendBorderSize(0);
    TLegend *leg = new TLegend(0.355,0.675,0.585,0.875);
    leg->SetTextSize(0.055);
    leg->AddEntry(histoData, (histName.Contains("Monte") ? "Monte Carlo" : "Data"), (type.Contains("dot") ? "lep" : "l"));
    leg->AddEntry(f, "Exponential fit","l");
    
    histoData->Draw((type.Contains("dot") ? "E1" : ""));
    leg->Draw("same");

    // add axis labels
    histoData->GetXaxis()->SetTitle(xTitle);
    histoData->GetXaxis()->CenterTitle();
    histoData->GetXaxis()->SetTitleSize(0.055);
    histoData->GetXaxis()->SetTitleOffset(0.84);
    histoData->GetXaxis()->SetLabelOffset(0.010);
    histoData->GetXaxis()->SetLabelSize(0.05);
    histoData->GetYaxis()->SetTitle(yTitle);
    histoData->GetYaxis()->CenterTitle();
    histoData->GetYaxis()->SetTitleSize(0.055);
    histoData->GetYaxis()->SetTitleOffset(1.0);
    histoData->GetYaxis()->SetLabelSize(0.05);
    gStyle->SetTitleSize(0.08, "t");
    histoData->SetTitle("");
    can->Update();

    can->SaveAs(name + ".png");
}

void histogram(TH1D *histoData, TH1D *histoMC, const TString histName, TCanvas *can, const TString xTitle, const TString yTitle, const TString name, const TString type){

    gPad->SetLogy( (type.Contains("log") ? 1 : 0) );

    if(!histoData || !histoMC) return;

    if(type.Contains("dot")){
        histoData->SetLineWidth(1);
        histoData->SetMarkerStyle(20);
        histoData->SetMarkerSize(0.9);
        histoData->SetMarkerColor(kRed);
        histoData->SetLineColor(kBlack);
        histoMC->SetLineWidth(1);
        histoMC->SetMarkerStyle(21);
        histoMC->SetMarkerSize(0.9);
        histoMC->SetMarkerColor(kBlue);
        histoMC->SetLineColor(kBlack);
    }

    else if(type.Contains("line")){
        histoData->SetLineWidth(1);
        histoData->SetLineColor(kBlue);
        histoMC->SetLineWidth(1);
        histoMC->SetLineColor(kRed);
    }
    else return;

    gStyle->SetLegendBorderSize(0);
    TLegend *leg = new TLegend(0.355,0.675,0.585,0.875);
    leg->SetTextSize(0.055);
    leg->AddEntry(histoData, "Data", (type == "dot" ? "lep" : "l"));
    leg->AddEntry(histoData, "Monte Carlo", (type == "dot" ? "lep" : "l"));

    if(histoData->GetMaximum() > histoMC->GetMaximum()){
        histoData->Draw((type.Contains("dot") ? "E1" : ""));
        histoMC->Draw((type.Contains("dot") ? "same E1" : "same"));

        // add axis labels
        histoData->GetXaxis()->SetTitle(xTitle);
        histoData->GetXaxis()->CenterTitle();
        histoData->GetXaxis()->SetTitleSize(0.055);
        histoData->GetXaxis()->SetTitleOffset(0.84);
        histoData->GetXaxis()->SetLabelOffset(0.010);
        histoData->GetXaxis()->SetLabelSize(0.05);
        histoData->GetYaxis()->SetTitle(yTitle);
        histoData->GetYaxis()->CenterTitle();
        histoData->GetYaxis()->SetTitleSize(0.055);
        histoData->GetYaxis()->SetTitleOffset(1.0);
        histoData->GetYaxis()->SetLabelSize(0.05);
        gStyle->SetTitleSize(0.08, "t");
        histoData->SetTitle("");
    }
    else{
        histoMC->Draw((type == "dot" ? "E1" : ""));
        histoData->Draw((type == "dot" ? "same E1" : "same"));

        // add axis labels
        histoMC->GetXaxis()->SetTitle(xTitle);
        histoMC->GetXaxis()->CenterTitle();
        histoMC->GetXaxis()->SetTitleSize(0.055);
        histoMC->GetXaxis()->SetTitleOffset(0.84);
        histoMC->GetXaxis()->SetLabelOffset(0.010);
        histoMC->GetXaxis()->SetLabelSize(0.05);
        histoMC->GetYaxis()->SetTitle(yTitle);
        histoMC->GetYaxis()->CenterTitle();
        histoMC->GetYaxis()->SetTitleSize(0.055);
        histoMC->GetYaxis()->SetTitleOffset(1.0);
        histoMC->GetYaxis()->SetLabelSize(0.05);
        gStyle->SetTitleSize(0.08, "t");
        histoMC->SetTitle("");
    }
    
    leg->Draw("same");

    can->Update();

    can->SaveAs(name + ".png");
}

//--------------------------------------------------------------------------------------------------
//
// Macro to parse .rad files and combine measurements from the same objects into single ROOT 
// TProfiles (http://root.cern.ch/root/html/TProfile.html)
//
// run over file called "inputfile.rad" with:
//
//   root (-l -q) parser.C+\(\"inputfile\"\)
//
// 1/21/14 Jay Lawhorn
//--------------------------------------------------------------------------------------------------
#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TProfile.h>
#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <TCanvas.h>
#include <TH1.h>
#include <TStyle.h>
#endif

// string parser for targets of interest
Bool_t nameMatch(string line);
Bool_t badData(TProfile*);
void histogram(TH1D *histo, const TString graphName, TCanvas *can, const TString xTitle, const TString yTitle, const TString name);

// variable for histograms
UInt_t h = 0;

TH1D *vlsr = new TH1D("vlsr", "vlsr", 1000, -50, 50);

void parser(TString infile="try2" // name of file to parse
    ) { // start main fxn

    // string variables
    char filename[50];
    ifstream ifs;
    string line;

    // variables to read in
    TString datetime;
    Float_t azimuth, elevation, az_off, ele_off, glon, glat, startF, stepF;
    Int_t mode, nStep;
    Float_t mag;

    // open file, check that file is open
    sprintf(filename, "%s.rad",infile.Data());
    ifs.open(filename); assert(ifs.is_open());

    // make new output file with same name as input file
    sprintf(filename, "%s.root",infile.Data());
    TFile *outFile = new TFile(filename, "RECREATE");
    vector<TProfile*> outVect;

    Int_t n=0;                  // line counter (debugging purposes)
    Int_t found=0;              // counter variable for asterisks
    TString sampleType;         // string to hold sample name
    Bool_t store=kFALSE;        // boolean flag for trying to parse and store next line
    Bool_t repeatSource=kFALSE; // boolean flag for making a new TProfile or appending to an old one
    Int_t nRepeat=0;            // integer used to set boolean flag for repeatSource
    Double_t tsys=350;          // t_sys variable
    TString blah1, blah2;       // placeholder variables

    while(getline(ifs,line)) { // loop over lines
        stringstream ss(line);
        stringstream ss2(line);
        n++; // increment line counter

        if (nameMatch(line)==1) { // if the line contains sample info of interest
            n=0; // reset line counter

            if(h!=0) cout << "VLSR was " << vlsr->GetMean() << " +- " << vlsr->GetStdDev() << " km/s" << endl;
            vlsr->Clear();
            vlsr->Reset();

            found = line.find("*",0); // find first asterisk in line

            if (line.find(":",0)>0) { // if there is a colon in the line
                found = line.find(":",0); // find the first colon
                found = line.find(":",found+1); // and the second - sample title is actually everything after second colon
            }

            sampleType=line.substr(found+2); // either save whole line or just what's after second colon
            sampleType.Resize(sampleType.Sizeof()-2);

            cout << sampleType.Data() << endl; // print sample name

            store=kTRUE; // set flag to store next line

            continue; // go to next line

        } // end name match statement

        if (line.find("tsys")!=string::npos) { // if the line contains "tsys"

            ss >> blah1 >> blah2 >> tsys; // read in the t_sys

            cout << "Tsys: " << tsys << " DeltaT: " << tsys/TMath::Sqrt(0.0078125e6*0.52488) << endl; // output t_sys and associated uncertainty
            // uncertainty from "Useful Technical Information" pg 15"

        } // end "tsys" statement

        if (line[0]=='*') continue; // skip lines that have leading asterisk

        if (tsys==0) continue; // if tsys isn't set, skip following lines

        if (store==kTRUE) { 

            ss >> datetime >> azimuth >> elevation >> az_off >> ele_off >> glon >> glat >> startF >> stepF >> mode >> nStep;
            ss2 >> datetime >> azimuth >> elevation >> az_off >> ele_off >> glon >> glat >> startF >> stepF >> mode >> nStep;

            if (outVect.size()==0) { // if output vector is empty

                outVect.push_back(new TProfile(sampleType.Data(), sampleType.Data(), nStep, startF, startF+nStep*stepF)); // make new TProfile of appropriate size

                TProfile *tempP = new TProfile((TString::Format("%s_%1u", sampleType.Data(), h)).Data(), sampleType.Data(), nStep, startF, startF+nStep*stepF);

                // fill temp TProfile
                for (Int_t i=0; i<nStep; i++) { ss >> mag; if ((i>10)&&(i<nStep-10)) tempP->Fill(startF+i*stepF,mag,tsys/TMath::Sqrt(0.0078125e6*0.52488));}

                // check if data point is good
                if(badData(tempP)){
                    delete tempP;
                    continue;
                }

                // fill new TProfile
                for (Int_t i=0; i<nStep; i++) { ss2 >> mag; outVect[outVect.size()-1]->Fill(startF+i*stepF,mag,tsys/TMath::Sqrt(0.0078125e6*0.52488));} 
                //outVect[outVect.size()-1]->Add(tempP);
                delete tempP;
                ss >> blah1 >> blah2;
                vlsr->Fill(atof(blah2.Data()));

            }
            else {
                nRepeat=-1;
                for (Int_t i=0; i<outVect.size(); i++) {
                    // check if sample name already exists
                    if (sampleType.CompareTo(outVect[i]->GetName()) == 0) {repeatSource=kTRUE; nRepeat=i; }

                }

                if (nRepeat==-1) { repeatSource=kFALSE; }

                if (repeatSource) { // if it already exists

                    TProfile *tempP = new TProfile((TString::Format("%s_%1u", sampleType.Data(), h)).Data(), sampleType.Data(), nStep, startF, startF+nStep*stepF);

                    // fill temp TProfile
                    for (Int_t i=0; i<nStep; i++) { ss >> mag; tempP->Fill(startF+i*stepF,mag,tsys/TMath::Sqrt(0.0078125e6*0.52488));}

                    // check if data point is good
                    if(badData(tempP)){
                        //delete tempP;
                        continue;
                    }

                    // fill extant TProfile
                    for (Int_t i=0; i<nStep; i++) { ss2 >> mag; outVect[nRepeat]->Fill(startF+i*stepF,mag, tsys/TMath::Sqrt(0.0078125e6*0.52488)); }
                    //outVect[nRepeat]->Add(tempP);
                    delete tempP;
                    ss >> blah1 >> blah2;
                    vlsr->Fill(atof(blah2.Data()));
                }
                else { // if it doesn't exist, make new TProfile and fill it
                    outVect.push_back(new TProfile(sampleType.Data(), sampleType.Data(), nStep, startF, startF+nStep*stepF));

                    TProfile *tempP = new TProfile((TString::Format("%s_%1u", sampleType.Data(), h)).Data(), sampleType.Data(), nStep, startF, startF+nStep*stepF);

                    // fill temp TProfile
                    for (Int_t i=0; i<nStep; i++) { ss >> mag; tempP->Fill(startF+i*stepF,mag,tsys/TMath::Sqrt(0.0078125e6*0.52488));}

                    // check if data point is good
                    if(badData(tempP)){
                        //delete tempP;
                        continue;
                    }

                    for (Int_t i=0; i<nStep; i++) { ss2 >> mag; outVect[outVect.size()-1]->Fill(startF+i*stepF,mag, tsys/TMath::Sqrt(0.0078125e6*0.52488));}
                    //outVect[outVect.size()-1]->Add(tempP);
                    delete tempP;
                    ss >> blah1 >> blah2;
                    vlsr->Fill(atof(blah2.Data()));
                }
            }
        }
    } // end line loop

    if(h!=0) cout << "VLSR was " << vlsr->GetMean() << " +- " << vlsr->GetStdDev() << " km/s" << endl;

    ifs.close(); // close input file

    TCanvas *can = new TCanvas("Histogram", "Histogram", 1600, 900);

    // write output vector to file
    for (UInt_t i=0; i<outVect.size(); i++) {
        outVect.at(i)->Write();
        //outVect.at(i)->Draw();
        //can->SaveAs(TString::Format("hist%1i.png", i));
        //histogram(outVect.at(i), "Solar radio emission", can, "Frequency (MHz)", "Brightness temperature (K)", TString::Format("histo_%1i", i));
    }
    cout << "Saved branches for " << outVect.size() << " objects" << endl;

    // save and close file
    //outFile->Write();
    outFile->Close();

} // end main fxn

//
// Function to recognize lines of interest based on string parsing. Look for lines containing "G" but not 
// "G " (some garbage lines in file have G at end of word instead of G00, G80, etc), "Sun", "npoint", "offset",
// "galactic", but not containing "file", "record", "freq".
//

Bool_t nameMatch(string line) {

    return  ((((line.find("G", 0)!=string::npos) && (line.find("G ",0,2)==string::npos)) || ((line.find("Sun",0)!=string::npos) && (line.find("npoint",0)==string::npos)) || (line.find("offset",0)!=string::npos) || (line.find("galactic",0)!=string::npos)) && (line.find("file",0)==string::npos) && (line.find("record",0)==string::npos) && (line.find("freq",0)==string::npos));

}

Bool_t badData(TProfile *p){

    if(!p) return kTRUE;

    h++;

    TH1D *histo = new TH1D((TString::Format("tempHisto_%1u", h)).Data(), (TString::Format("tempHisto_%1u", h)).Data(), 1000,
        p->GetBinContent(p->GetMinimumBin()), p->GetBinContent(p->GetMaximumBin()));

    for(UInt_t bin = 11; bin <= p->GetXaxis()->GetNbins()-2-10; bin++){
        histo->Fill(p->GetBinContent(bin));
    }

    Double_t stdDev = histo->GetStdDev();

    //delete histo;

    if(stdDev > 150) return kTRUE;

    return kFALSE;
}

void histogram(TH1D *histo, const TString graphName, TCanvas *can, const TString xTitle, const TString yTitle, const TString name){

    gStyle->SetOptStat(2210);
    //gStyle->SetOptFit(1111);
    //gStyle->SetOptStat(kFALSE);
    gStyle->SetOptFit(100);

    if(!histo){
        cout << "Error: Graph \"" << graphName << "\" not defined" << endl;
        return;
    }

    histo->SetLineWidth(2);
    histo->SetMarkerStyle(20);
    histo->SetMarkerSize(1.5);
    histo->SetMarkerColor(kRed);
    histo->SetLineColor(kBlack);

    // gStyle->SetLegendBorderSize(0);
    // TLegend *leg = new TLegend(0.655,0.675,0.885,0.875);
    // leg->SetTextSize(0.055);
    // leg->AddEntry(graph, "Data","lep");
    // leg->AddEntry(f, "Gaussian fit","l");
    
    histo->Draw();
    // leg->Draw("same");

    // add axis labels
    histo->GetXaxis()->SetTitle(xTitle);
    histo->GetXaxis()->CenterTitle();
    histo->GetXaxis()->SetTitleSize(0.055);
    histo->GetXaxis()->SetTitleOffset(0.90);
    histo->GetXaxis()->SetLabelOffset(0.010);
    histo->GetXaxis()->SetLabelSize(0.05);
    histo->GetYaxis()->SetTitle(yTitle);
    histo->GetYaxis()->CenterTitle();
    histo->GetYaxis()->SetTitleSize(0.055);
    histo->GetYaxis()->SetTitleOffset(1.0);
    histo->GetYaxis()->SetLabelSize(0.05);
    gStyle->SetTitleSize(0.08, "t");
    histo->SetTitle(graphName);
    can->Update();

    can->SaveAs(name + ".png");

    can->Clear();
}

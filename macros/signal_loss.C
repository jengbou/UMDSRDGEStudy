//STL headers
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cmath>

//ROOT headers
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH3.h"
#include "TH2.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TSpectrum.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TLine.h"
#include "TStyle.h"
#include "TMath.h"
#include "TGaxis.h"
#include "TExec.h"

//CMSSW headers
#include "HEDarkening.C"

#define maxHDeta 14
#define maxHDlay 19

using namespace std;

//---------------------------------------------------
//function to make plot of signal loss for given lumi
void signal_loss(double lumi, bool do_print=false){
	gStyle->SetPalette(1); //rainbow
	
	//HEDarkening object for darkening weights
	HEDarkening darkening;
	
	TH2F* hist = new TH2F("signal_loss","",18,-0.5,17.5,14,-29.5,-15.5);

	for(int j = 0; j < maxHDeta; j++){
		//fill histo (negative ieta)
		for(int i = 0; i <= maxHDlay; i++){
			int bin = hist->FindBin(i,-(j+16));
			double sigloss = darkening.degradation(lumi,j+16,i);
			hist->SetBinContent(bin,sigloss); //be careful of eta and layer numbering
			//hist->SetBinContent(i+1,14-j,darkening.degradation(lumi,j+16,i)); //be careful of eta and layer numbering
		}
	}

	//formatting
	hist->SetTitle("");
	hist->GetXaxis()->SetTitle("Layer");
	hist->GetYaxis()->SetTitle("Tower");
	stringstream sz;
	sz << "relative signal after " << lumi << " fb^{-1}";
	hist->GetZaxis()->SetTitle((sz.str()).c_str());
	hist->GetZaxis()->SetRangeUser(0.,1.);
	
	//drawing
	TCanvas* can = new TCanvas("signal_loss","signal_loss",900,600);
	TPad* pad1 = new TPad("graph","",0,0,1,1);
	pad1->SetMargin(0.1,0.175,0.125,0.05);
	pad1->SetTicks(1,1);
	pad1->Draw();
	pad1->cd();

	//more formatting
	hist->GetYaxis()->SetTitleOffset(0.75);
	hist->GetZaxis()->SetTitleOffset(0.85);
	hist->GetXaxis()->SetTitleOffset(0.85);
	hist->GetZaxis()->SetTitleSize(38/(pad1->GetWh()*pad1->GetAbsHNDC()));
	hist->GetZaxis()->SetLabelSize(32/(pad1->GetWh()*pad1->GetAbsHNDC()));
	hist->GetYaxis()->SetTitleSize(38/(pad1->GetWh()*pad1->GetAbsHNDC()));
	hist->GetYaxis()->SetLabelSize(32/(pad1->GetWh()*pad1->GetAbsHNDC()));
	hist->GetXaxis()->SetTitleSize(38/(pad1->GetWh()*pad1->GetAbsHNDC()));
	hist->GetXaxis()->SetLabelSize(32/(pad1->GetWh()*pad1->GetAbsHNDC()));
	hist->GetYaxis()->SetTickLength(12/(pad1->GetWh()*pad1->GetAbsHNDC()));
	hist->GetXaxis()->SetTickLength(12/(pad1->GetWh()*pad1->GetAbsHNDC()));
	
	//gStyle->SetPaintTextFormat(".2g"); //g is 'adaptive', uses shorter of e or f
	gStyle->SetPaintTextFormat(".2f");
	hist->Draw("COLZ TEXT");
	
	// Remove the current axis
	hist->GetYaxis()->SetLabelOffset(999);
	//hist->GetYaxis()->SetTickLength(0);
	
	// Redraw the new axis 
	gPad->Update();
	TGaxis *newaxis = new TGaxis(gPad->GetUxmin(), 
								gPad->GetUymax(),
								gPad->GetUxmin()-0.001,
								gPad->GetUymin(),
								-hist->GetYaxis()->GetXmax(),
								-hist->GetYaxis()->GetXmin(),
								510,"+");
	newaxis->SetTitleSize(38/(pad1->GetWh()*pad1->GetAbsHNDC()));
	newaxis->SetLabelSize(32/(pad1->GetWh()*pad1->GetAbsHNDC()));
	//newaxis->SetTickSize(12/(pad1->GetWh()*pad1->GetAbsHNDC()));
	newaxis->SetTickSize(0);
	newaxis->SetLabelOffset(-0.05);
	newaxis->SetLabelFont(42);
	newaxis->Draw();
	
	if(do_print){
		stringstream outname;
		outname << "signal_loss_lumi" << lumi;
		outname << ".png";
		can->Print((outname.str()).c_str(),"png");
	}
}
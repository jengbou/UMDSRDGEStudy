#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <map>
#include <TROOT.h>
#include <TMath.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TImageDump.h>
#include <TFile.h>
#include <TTree.h>
#include "Math/Interpolator.h"
#include <AnalysisOutput.C>

using TMath::Exp;
using TMath::Power;
using TMath::Log;

void readRatios(string inputFile);

string constructFileName(const char* runName, int layerNo, int ieta)
{
	TString fileName("../RUNNAME/Analysis-RUNNAME_layerLAYERNUMBER_ietaIETA.txt");
	fileName.ReplaceAll("RUNNAME", runName);
	fileName.ReplaceAll("LAYERNUMBER", TString::Itoa(layerNo, 10));
	fileName.ReplaceAll("IETA", TString::Itoa(ieta, 10));
	string fileNameString(fileName.Data());
	return  fileNameString;
}

void processDoseLine(string line, int& ieta, double& dose, double& fluence)
{
	stringstream ss(line);
	string field;
	int ctr = 0; //column index of tab-delimited row
	while(getline(ss, field, '\t'))
	{
		if(ctr==0) {stringstream temp(field); temp >> ieta;}
		else if (ctr==4) {stringstream temp(field); temp >> dose;}
		else if (ctr==5) {stringstream temp(field); temp >> fluence;}
		else {}
		
		ctr++;
	}
}

void readDoseFile(string inputFile, std::map<int, double>& doseMap, std::map<int, double>& fluenceMap)
{
	ifstream inputStream(inputFile.c_str());
	string line;
	if(inputStream.is_open())
	{
		while(getline(inputStream, line))
		{
			if(line[0]=='#') continue; //Skip comment lines
			else 
			{
				//Process the line to get ieta and dose
				int ieta;
				double dose, fluence;
				processDoseLine(line, ieta, dose, fluence);
				doseMap.insert(std::pair<int, double>(ieta, dose));
				fluenceMap.insert(std::pair<int, double>(ieta, fluence / 1e12));
			}
		}
	}
	else {cout <<"Error: Input file not found." <<endl;}
}

void dosePerTile()
{
	////////////////////////////////////////////////////////
	// Calculate dose per 100fb^-1 for each tile
	////////////////////////////////////////////////////////
	std::map<int, double> doseMapLayer1, fluenceMapLayer1;
	readDoseFile("doseLayer1.txt", doseMapLayer1, fluenceMapLayer1);
	cout << "ieta \t dose" << endl;
	for (int i=17; i<=29; i++) {
		cout << i << "\t" << doseMapLayer1[i] << endl;
	}
}

// Function to model efficiency vs attenuation length.
Double_t myfunction(const Double_t *x, Double_t *par)
{
	Double_t xx =x[0];
	Double_t f = par[1] * Exp(-par[0]/(xx));
	return f;
}

// Function to model efficiency vs 1 / attenuation length.
Double_t myfunction2(const Double_t *x, Double_t *par)
{
	Double_t xx =x[0];
	Double_t f = par[1] * Exp(-par[0]*(xx));
	return f;
}

// Function to model efficiency vs 1 / attenuation length.
Double_t EfficiencyVsMu(const Double_t *x, Double_t *par)
{
	Double_t xx =x[0];
	Double_t f = Exp(par[3]*pow(xx, 3) + par[2]*pow(xx, 2) + par[1]*(xx) + par[0]);
	return f;
}

// Inverse of EfficiencyVsMu
/*
Double_t MuVsEfficiency(const Double_t *x, Double_t *par)
{
	Double_t xx =x[0];
	Double_t f = (-par[1] - sqrt(pow(par[1], 2) - 4*par[2]*(par[0] - Log(xx)))) / (2*par[2]);
	return f;
}
*/

// Test EfficiencyVsMu and inverse
/*
void test()
{
	double par[3] = {-1.5, -46, 254};
	double x[1], y[1], z[1];
	x[0] = 0.0101;
	y[0] = EfficiencyVsMu(x, par);
	z[0] = MuVsEfficiency(y, par);
	cout << y[0] << endl;
	cout << z[0] << endl;
}
*/

// Function to model efficiency vs 1 / attenuation length.
Double_t LogEfficiencyVsMu(const Double_t *x, Double_t *par)
{
	Double_t xx = x[0];
	Double_t f = par[2] * pow(xx, 2) + par[1] * xx + par[0];
	return f;
}

// Inverse of myfunction 
Double_t invmyfunction(const Double_t *f, Double_t *par)
{
	Double_t ff = f[0];
	Double_t x = par[0] / Log( par[1]/ff );
	return x;
}

// Inverse of myfunction2 
Double_t invmyfunction2(const Double_t *f, Double_t *par)
{
	Double_t ff = f[0];
	Double_t x = - Log(ff / par[1]) / par[0];
	return x;
}

// Function to map dose in kGy to mu(Tile).
// mu(Fiber) is implicitly set in simulation. 
Double_t DoseToMu(const Double_t *x, Double_t *par)
{
	Double_t xx = x[0];
	Double_t f = par[1] * xx + par[0];
	return f;
}

// Function to map (dose, dose rate) in kGy to mu(Tile).
// mu(Fiber) is implicitly set in simulation. 
Double_t CalculateMu(const Double_t *x, Double_t *par)
{
	Double_t dose = x[0];
	Double_t doseRate = x[1];
	Double_t doseRateBase = par[2];
	Double_t A = par[3];
	Double_t f = ( 1 - A * log10(doseRate/doseRateBase) ) * par[1] * dose + par[0];
	return f;
}

void myfunc()
{
//	TF1 *f1 = new TF1("myfunc", myfunction, 0.1, 100.0, 2);
//	Double_t par[] = {1.0, 1.0};
//	f1->SetParameters(par);
//	f1->SetParNames("average track length", "constant");
//	TGraphErrors* gr = longrun6test();
//	gr->GetXaxis()->SetTitle("Absorption length (cm)");
//	gr->GetYaxis()->SetTitle("Efficiency");
//	gr->GetXaxis()->SetLimits(0.1,300);
//	gr->GetYaxis()->SetRangeUser(1e-7,1e-0);
//	gr->Draw("AP");
//	gr->Fit(f1, "", "", 0.5, 100);
//	f1->GetParameters(par);
//	cout << par[0] << "\t" << par[1] << endl;
}

void processLaserDataLine(string line, double& lumi, double* ratios)
{
	// ratio[0] through ratio[12] corresponds to ieta 17 to 29.
	// Make sure "ratios" has 13 elements.
	stringstream ss(line);
	string field;
	int ctr = 0; //column index of tab-delimited row
	while(getline(ss, field, '\t'))
	{
		if(ctr==0) {
			stringstream temp(field);
			temp >> lumi;
			//cout << lumi;
		}
		else if (true) {//ctr >= 1 && ctr <= 13) {
			stringstream temp(field); 
			double ratio;
			temp >> ratio;
			//cout << ratio << "\t";
			ratios[ctr - 1] = ratio;
		}
		else {
			cout << endl;
			break;
		}
		
		ctr++;
	}
}

void readLaserDataFile(string inputFile, vector<double>& lumis, std::map< int, vector<double> > & ratiosMap)
{
	ifstream inputStream(inputFile.c_str());
	string line;
	if(inputStream.is_open())
	{
		while(getline(inputStream, line))
		{
			if(line[0]=='#') continue; //Skip comment lines
			else 
			{
				//Process the line to get ieta and dose
				double lumi;
				double ratios[13];
				processLaserDataLine(line, lumi, ratios);
				lumis.push_back(lumi);
				for (int i=0; i < 13; i++) {
					ratiosMap[i+17].push_back(ratios[i]);
					//cout << "ratios[i] " <<ratios[i] << endl;
					//cout << "ratiosMap[i+17] " << (ratiosMap[i+17])[1] << endl;
				}
			}
		}
	}
	else {cout <<"Error: Input file not found." <<endl;}
}

void readLaserDataFile2(string inputFile, vector<double>& lumis, std::map< int, vector<double> > & ratiosMap)
{
	ifstream inputStream(inputFile.c_str());
	string line;
	if(inputStream.is_open())
	{
		while(getline(inputStream, line))
		{
			if(line[0]=='#') continue; //Skip comment lines
			else 
			{
				//Process the line to get ieta and dose
				double lumi;
				double ratios[13];
				processLaserDataLine(line, lumi, ratios);
				lumis.push_back(lumi);
				for (int i=0; i < 13; i++) {
					ratiosMap[i+17].push_back(ratios[i]);
					//cout << "ratios[i] " <<ratios[i] << endl;
					//cout << "ratiosMap[i+17] " << (ratiosMap[i+17])[1] << endl;
				}
			}
		}
	}
	else {cout <<"Error: Input file not found." <<endl;}
}

void ratiosMap(int ieta = 29) //bool doPlot = false)
{
	vector<double> lumis;
	std::map<int, vector<double> > ratiosMap;

	readLaserDataFile("../HELaserData.txt", lumis, ratiosMap);

	TCanvas* c1 = new TCanvas("c1","canvas1",800,600);
	TGraph* graph;
	graph = new TGraph(lumis.size(), &lumis[0], &((ratiosMap[ieta])[0]));
	c1->cd();
	graph->Draw("AP");
}

// Function to model attenuation length vs dose 
Double_t LambdaVsDose(const Double_t *x, Double_t *par)
{
	Double_t xx =x[0];
	Double_t f = 1 / (par[0] * xx + par[1]);
	return f;
}

// Function to model 1/attenuation length vs dose 
Double_t LambdaVsDose2(const Double_t *x, Double_t *par)
{
	Double_t xx =x[0];
	Double_t f = par[0] * (xx - 0) + par[1];
	return f;
}

void readRatios(string inputFile)
{
	TFile* file = new TFile("dataRatios.root", "recreate");
	TTree* dataRatios = new TTree("dataRatios", "ratios from laser data");
	Int_t ieta, lumiIndex;
	lumiIndex = 0;
	Double_t lumi, ratio;
	dataRatios->Branch("ieta", &ieta, "ieta/I");
	dataRatios->Branch("lumiIndex", &lumiIndex, "lumiIndex/I");
	dataRatios->Branch("lumi", &lumi, "lumi/D");
	dataRatios->Branch("ratio", &ratio, "ratio/D");

	ifstream inputStream(inputFile.c_str());
	string line;
	if(inputStream.is_open())
	{
		while(getline(inputStream, line))
		{
			if(line[0]=='#') continue; //Skip comment lines
			else {
				// Process the line to get ieta and dose.
				//double lumi;
				double ratios[13];
				processLaserDataLine(line, lumi, ratios);
				for (int i=0; i < 13; i++) {
					ieta = i + 17;
					ratio = ratios[i];

					cout << "ieta: " << ieta << endl; 
					cout << "lumiIndex: " << lumiIndex << endl; 
					cout << "lumi: " << lumi << endl; 
					cout << "ratio: " << ratio << endl << endl; 

					dataRatios->Fill();
					//cout << "ratios[i] " <<ratios[i] << endl;
				}
				lumiIndex = lumiIndex +1; 
			}
		}
	}
	else {cout <<"Error: Input file not found." <<endl;}
	dataRatios->BuildIndex("ieta", "lumiIndex");
	dataRatios->Write();
}

void ratiosPlot() //bool doPlot = false)
{
	readRatios("../HELaserData.txt");

//	TCanvas* c1 = new TCanvas("c1","canvas1",800,600);
//	TGraph* graph;
//	graph = new TGraph(lumis.size(), &lumis[0], &((ratiosMap[ieta])[0]));
//	c1->cd();
//	graph->Draw("AP");
}


// Copy of MakeTrees to enable fitting of dose factor.
void FitEffVsMu(const char* runName, bool boolTileFiber)//, const char* treeLabel) 
{
// Vector of colors.
	vector<Color_t> colors(18, kBlack);
	{
		colors[0] = kBlue;
		colors[1] = kOrange;
		colors[2] = kGreen;
		colors[3] = kRed;
		colors[4] = kCyan;
		colors[5] = kMagenta;
		// Define colors[6] through colors[11]
		for (int i = 0; i < 6; i++ ) {
			colors[i+6] = colors[i] + 1;
		}
		// Define colors[12] through colors[17]
		for (int i = 0; i < 6; i++ ) {
			colors[i+12] = colors[i] + 2;
		}
	}
	
	TCanvas* c1;
	c1 = new TCanvas("c1", "Efficiency vs Mu", 800, 600);
	//gPad->SetLogx();
	//gPad->SetLogy();
	
	// graphsAll contains Efficiency vs Mu data for all tiles from the Geant4 simulations.
	// namesAll contains names for all graphs in graphsAll.
	// graphsTypeX contains Efficiency vs Lambda data for tiles of type X from the Geant4 simulations.
	vector<TGraphErrors*> graphsLayer1, graphsLayer7, graphsAll, graphsAllLog;
	vector<string> namesAll;

	// Read from Analysis files and add to graphsTypeX, and graphsAll.
	// (Only reading type 1 tiles.)
	{
		for(int ieta=22; ieta<30; ieta++) {
			AnalysisOutput* analysis = new AnalysisOutput(constructFileName(runName, 1, ieta));
			if(boolTileFiber) {graphsLayer1.push_back(analysis->graphTile);}
			else {graphsLayer1.push_back(analysis->graphFiber);}
		}
		
		// Add all graphs from graphsLayer1 to graphsAll. Add all legend names to namesAll.
		for (unsigned int i=0; i < graphsLayer1.size(); i++) {
			graphsAll.push_back(graphsLayer1[i]);
			stringstream name;
			name << "ieta " << (i+22) << "Layer1";
			namesAll.push_back(name.str());
		}
	}

	TLegend* leg = new TLegend(0.7, 0.6, 0.9, 0.85);
	// Run through all graphs in graphsAll and modify error bars if necessary.
	// Push log version of graphsAll into graphsAllLog.
	// Set graph properties of all graphs in graphsAll and draw to c1 pad 1.
	// Add all legend names from namesAll to leg.
	for (unsigned int i=0; i < graphsAll.size(); i++) {
		//cout << graphsAll.size() << ": " << graphsAll[i] << endl;
		
		int n = graphsAll[i]->GetN();
		double* xV = graphsAll[i]->GetX();
		// Offset x-values for clarity.
		// (Disabled)
		for (int j=0; j < n; j++) {
			//xV[j] += i*0.0002;
		}
		double* yV = graphsAll[i]->GetY();
		double* exV = graphsAll[i]->GetEX();
		double* eyV = graphsAll[i]->GetEY();
	
		graphsAll[i] = new TGraphErrors(n, xV, yV, exV, eyV);
		graphsAll[i]->SetTitle("");//("Light collection efficiency vs Absorption Length");
		graphsAll[i]->GetXaxis()->SetTitle("Absorption coefficient [cm^{-1}]");
		graphsAll[i]->GetYaxis()->SetTitle("Efficiency");
		//graphsAll[i]->GetXaxis()->SetLimits(0.0, 0.05);
		graphsAll[i]->GetYaxis()->SetRangeUser(1e-5, 2e-1);
		graphsAll[i]->SetMarkerColor(colors[i]);
		graphsAll[i]->SetMarkerSize(0.5);
		graphsAll[i]->SetMarkerStyle(21);

		graphsAllLog.push_back( (TGraphErrors*) graphsAll[i]->Clone() );

		c1->cd(1);
		c1->SetLogy();
		if(i==0) {
			graphsAll[i]->Draw("AP");
		}
		else {
			graphsAll[i]->Draw("P same");
		}
	
		leg->AddEntry(graphsAll[i], namesAll[i].c_str(), "p");
	}

	string outFileName = "FitEffVsMu-";
	outFileName.append(runName);
	outFileName.append(".txt");
	ofstream outfile;
	outfile.open(outFileName.c_str());

	vector<TF1*> functionV (graphsAll.size());
	for (unsigned int i=0; i < graphsAll.size(); i++) {
		functionV[i] = new TF1("EfficiencyVsMu", EfficiencyVsMu, 0.0, 10.0, 4);
	}
	// Fit all graphs in graphsAll to EfficiencyVsMu.
	outfile << "layer\tieta\tp1\tp0" << endl;
	for (unsigned int i=0; i < graphsAll.size(); i++) {
		Double_t par[] = {1.0, 1.0, 1.0, -1.0};
		functionV[i]->SetParameters(par);
		functionV[i]->FixParameter(3, 0.0); //*-* Only fit to exp(linear function).
		functionV[i]->FixParameter(2, 0.0); //*-* Only fit to exp(linear function).
		functionV[i]->SetLineColor(colors[i]);
		functionV[i]->SetParNames("p0", "p1", "p2", "p3");
		graphsAll[i]->Fit(functionV[i], "", "", 0.01, 0.19);

		int ieta = (i % 8) + 22;
		int layer = i<8 ? 1 : 7;
		double p1 = functionV[i]->GetParameter("p1");
		double p0 = functionV[i]->GetParameter("p0");
		outfile << layer << "\t" << ieta << "\t";//*-*
		outfile << p1 << "\t" << p0 << endl;
	}

	// Make log version of graphsAllLog
	for (unsigned int i=0; i < graphsAllLog.size(); i++) {
		//cout << graphsAllLogLog.size() << ": " << graphsAllLog[i] << endl;
		
		int n = graphsAllLog[i]->GetN();
		double* xV = graphsAllLog[i]->GetX();
		double* yV = graphsAllLog[i]->GetY();
		for (int j=0; j < n; j++) {
			yV[j] = log10( yV[j] );
		}
		double* exV = graphsAllLog[i]->GetEX();
		double* eyV = graphsAllLog[i]->GetEY();

		graphsAllLog[i] = new TGraphErrors(n, xV, yV, exV, eyV);
		graphsAllLog[i]->SetTitle("");//("Light collection efficiency vs Absorption Length");
		graphsAllLog[i]->GetXaxis()->SetTitle("Absorption coefficient [cm^{-1}]");
		graphsAllLog[i]->GetYaxis()->SetTitle("Log(Efficiency)");
		//graphsAllLog[i]->GetXaxis()->SetLimits(-0.005, 0.045);
		//graphsAllLog[i]->GetYaxis()->SetRangeUser(1e-6, 1e0);
		graphsAllLog[i]->SetMarkerColor(colors[i]);
		graphsAllLog[i]->SetMarkerSize(0.5);
		graphsAllLog[i]->SetMarkerStyle(21);

		//c1->cd(1);
		//if(i==0) {
		//	graphsAllLog[i]->Draw("AP");
		//}
		//else {
		//	graphsAllLog[i]->Draw("P same");
		//}
	}

	// Draw leg to c1 pad 1.
	{
		c1->cd(1);
		leg->SetFillColor(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.03);
		leg->SetTextFont(42);
		leg->Draw();
	}

	/*
	vector<TF1*> functionLogV (graphsAllLog.size(), new TF1("LogEfficiencyVsMu", LogEfficiencyVsMu, 0.0, 0.0, 3));
	// Fit all graphs in graphsAllLog to myfunction.
	for (unsigned int i=0; i < graphsAllLog.size(); i++) {
		Double_t par[] = {1.0, 1.0, 1.0};
		functionLogV[i]->SetParameters(par);
		functionLogV[i]->SetLineColor(colors[i]);
		functionLogV[i]->SetParNames("p0", "p1", "p2");
		//*-*graphsAllLog[i]->Fit(functionLogV[i], "", "", 0.01, 0.045);
	}

	// Read laser data file.
	// ratiosMap[ieta] contains ratios for ieta at luminosity values contained in lumis.
	vector<double> lumis;
	std::map<int, vector<double> > ratiosMapData;
	readLaserDataFile("../HELaserData.txt", lumis, ratiosMapData);

	// Read dose and fluence map data files.
	std::map<int, double> doseMapLayer1;
	std::map<int, double> fluenceMapLayer1;
	readDoseFile("../doseLayer1Mars.txt", doseMapLayer1, fluenceMapLayer1); //*-*

	// Variables to read and write from Trees.
	Int_t ieta, lumiIndex;
	Double_t lumi, ratio;
	ieta = 0;
	lumiIndex = 0;
	lumi = 0.;
	ratio = 0.;

	// Prepare to read dataRatios Tree.
	TFile* dataFile = new TFile("dataRatios.root");
	if (!dataFile->IsOpen()) {cout << "dataRatios.root not found"; return;}
	TTree* dataRatios = (TTree*)dataFile->Get("dataRatios");
	dataRatios->Show(10);
	dataRatios->SetBranchAddress("ieta", &ieta);
	dataRatios->SetBranchAddress("lumiIndex", &lumiIndex);
	dataRatios->SetBranchAddress("lumi", &lumi);
	dataRatios->SetBranchAddress("ratio", &ratio);
	
	// Prepare mcRatios Tree to be written. mcRatios shares branch addresses with dataRatios.
	TString mcFileName(runName);
	mcFileName.Prepend("mcRatios-").Append("-").Append(treeLabel).Append(".root");
	TFile* file = new TFile(mcFileName.Data(), "recreate");
	if (!file->IsOpen()) {cout << "mcRatios.root not found"; return;}
	TTree* mcRatios = new TTree("mcRatios", "ratios from Geant4");
	mcRatios->Branch("ieta", &ieta, "ieta/I");
	mcRatios->Branch("lumiIndex", &lumiIndex, "lumiIndex/I");
	mcRatios->Branch("lumi", &lumi, "lumi/D");
	mcRatios->Branch("ratio", &ratio, "ratio/D");

	for (unsigned int i=0; i < graphsAll.size(); i++) {
		// Define parameters for DoseToMu and calculate undamaged efficiency.
		//double par[2] = {0.02, 0.002};
		double doseRateBase = 2.67e-3;
		double A = 0.5;
		double par[2] = {0.023, 0.001};
		double par2[4] = {0.023, 0.001, doseRateBase, A};
		double muTileUndamaged = 0.023;

		for (unsigned int j=0; j < lumis.size(); j++) {
			dataRatios->GetEntryWithIndex(i+22, j); // Get entry with arbitrary ieta and lumiIndex == j.
			//stringstream cutFormula; cutFormula << "ieta == " << i << " && lumiIndex == " << j;
			//cout << cutFormula.str() << endl;
			//dataRatios->Draw("ieta:lumi:ratio:lumiIndex", cutFormula.str().c_str(), "goff");
			//ieta = (dataRatios->GetVal(1))[0];

			//cout << "ieta: " << ieta << endl; 
			//cout << "lumiIndex: " << lumiIndex << endl; 
			//cout << "lumi: " << lumi << endl; 
			//cout << "ratio: " << ratio << endl; 

			double dose = doseFactor * (doseMapLayer1[ieta] / 100) * lumi; // Dose[kGy] per fb-1.
			double doseRate = doseMapLayer1[ieta] * 3.6e-5; // Dose rate [kGy/hr]
			double xx[2] = {dose, doseRate};
			double muTile = DoseToMu(&dose, par);
			//double muTile = CalculateMu(xx, par2);
			//double efficiency = functionV[i]->Eval(muTile);
			double efficiency = graphsAll[i]->Eval(muTile); //*-*
			if(j==0) {
				muTileUndamaged = muTile;
			}
			//double efficiencyUndamaged = functionV[i]->Eval(muTileUndamaged);
			double efficiencyUndamaged = graphsAll[i]->Eval(muTileUndamaged); //*-*
			ratio = efficiency / efficiencyUndamaged;
			double ratioFromData = (ratiosMapData[ieta])[j];

//			if(j == lumis.size()-1) {
//				cout << "ieta: " << i+22 << endl;
//				cout << "dose: " << dose << endl; 
//				cout << "muTile: " << muTile << endl; 
//				cout << "efficiency: " << efficiency << endl; 
//				cout << "efficiencyUndamaged: " << efficiencyUndamaged << endl; 
//				cout << "ratio: " << ratio << endl; 
//				cout << "ratioFromData: " << ratioFromData << endl; 
//				cout << endl;
//			}
			mcRatios->Fill();
		}
	}

	int nEntries = dataRatios->GetEntries();
	for(int i=0; i < nEntries; i++) {
//		dataRatios->GetEntry(i);
//
//		cout << "ieta: " << ieta << endl; 
//		cout << "lumiIndex: " << lumiIndex << endl; 
//		cout << "lumi: " << lumi << endl; 
//		cout << "ratio: " << ratio << endl; 
	}
	
	mcRatios->Write();
*/
}

// Copy of FitEffVsMu for rectangular tiles.
void FitEffVsMu2(const char* runName, bool boolTileFiber)//, const char* treeLabel) 
{
// Vector of colors.
	vector<Color_t> colors(18, kBlack);
	{
		colors[0] = kBlue;
		colors[1] = kOrange;
		colors[2] = kGreen;
		colors[3] = kRed;
		colors[4] = kCyan;
		colors[5] = kMagenta;
		// Define colors[6] through colors[11]
		for (int i = 0; i < 6; i++ ) {
			colors[i+6] = colors[i] + 1;
		}
		// Define colors[12] through colors[17]
		for (int i = 0; i < 6; i++ ) {
			colors[i+12] = colors[i] + 2;
		}
	}
	
	TCanvas* c1;
	c1 = new TCanvas("c1", "Efficiency vs Mu", 800, 600);
	//gPad->SetLogx();
	//gPad->SetLogy();
	
	// graphsAll contains Efficiency vs Mu data for all tiles from the Geant4 simulations.
	// namesAll contains names for all graphs in graphsAll.
	// graphsTypeX contains Efficiency vs Lambda data for tiles of type X from the Geant4 simulations.
	vector<TGraphErrors*> graphsLayer1, graphsLayer7, graphsAll, graphsAllLog;
	vector<string> namesAll;

	// Read from Analysis files and add to graphsTypeX, and graphsAll.
	// (Only reading type 1 tiles.)
	{
		{
			AnalysisOutput* analysis;

			analysis = new AnalysisOutput("../" + string(runName) + "/Analysis-" + string(runName)+ "_Dx10_Dy10.txt");
			if(boolTileFiber) {graphsLayer1.push_back(analysis->graphTile);}
			else {graphsLayer1.push_back(analysis->graphFiber);}
			analysis = new AnalysisOutput("../" + string(runName) + "/Analysis-" + string(runName)+ "_Dx20_Dy20.txt");
			if(boolTileFiber) {graphsLayer1.push_back(analysis->graphTile);}
			else {graphsLayer1.push_back(analysis->graphFiber);}
			analysis = new AnalysisOutput("../" + string(runName) + "/Analysis-" + string(runName)+ "_Dx8_Dy12.txt");
			if(boolTileFiber) {graphsLayer1.push_back(analysis->graphTile);}
			else {graphsLayer1.push_back(analysis->graphFiber);}
			analysis = new AnalysisOutput("../" + string(runName) + "/Analysis-" + string(runName)+ "_Dx8_Dy5.txt");
			if(boolTileFiber) {graphsLayer1.push_back(analysis->graphTile);}
			else {graphsLayer1.push_back(analysis->graphFiber);}
			analysis = new AnalysisOutput("../" + string(runName) + "/Analysis-" + string(runName)+ "_Dx12_Dy8.txt");
			if(boolTileFiber) {graphsLayer1.push_back(analysis->graphTile);}
			else {graphsLayer1.push_back(analysis->graphFiber);}
			analysis = new AnalysisOutput("../" + string(runName) + "/Analysis-" + string(runName)+ "_Dx5_Dy8.txt");
			if(boolTileFiber) {graphsLayer1.push_back(analysis->graphTile);}
			else {graphsLayer1.push_back(analysis->graphFiber);}
		}
		
		// Add all graphs from graphsLayer1 to graphsAll. Add all legend names to namesAll.
		for (unsigned int i=0; i < graphsLayer1.size(); i++) {
			graphsAll.push_back(graphsLayer1[i]);
			stringstream name;
			name << "file number " << i;
			namesAll.push_back(name.str());
		}
	}

	TLegend* leg = new TLegend(0.7, 0.6, 0.9, 0.85);
	// Run through all graphs in graphsAll and modify error bars if necessary.
	// Push log version of graphsAll into graphsAllLog.
	// Set graph properties of all graphs in graphsAll and draw to c1 pad 1.
	// Add all legend names from namesAll to leg.
	for (unsigned int i=0; i < graphsAll.size(); i++) {
		//cout << graphsAll.size() << ": " << graphsAll[i] << endl;
		
		int n = graphsAll[i]->GetN();
		double* xV = graphsAll[i]->GetX();
		// Offset x-values for clarity.
		// (Disabled)
		for (int j=0; j < n; j++) {
			//xV[j] += i*0.0002;
		}
		double* yV = graphsAll[i]->GetY();
		double* exV = graphsAll[i]->GetEX();
		double* eyV = graphsAll[i]->GetEY();
	
		graphsAll[i] = new TGraphErrors(n, xV, yV, exV, eyV);
		graphsAll[i]->SetTitle("");//("Light collection efficiency vs Absorption Length");
		graphsAll[i]->GetXaxis()->SetTitle("Absorption coefficient [cm^{-1}]");
		graphsAll[i]->GetYaxis()->SetTitle("Efficiency");
		//graphsAll[i]->GetXaxis()->SetLimits(0.0, 0.05);
		graphsAll[i]->GetYaxis()->SetRangeUser(1e-5, 2e-1);
		graphsAll[i]->SetMarkerColor(colors[i]);
		graphsAll[i]->SetMarkerSize(0.5);
		graphsAll[i]->SetMarkerStyle(21);

		graphsAllLog.push_back( (TGraphErrors*) graphsAll[i]->Clone() );

		c1->cd(1);
		c1->SetLogy();
		if(i==0) {
			graphsAll[i]->Draw("AP");
		}
		else {
			graphsAll[i]->Draw("P same");
		}
	
		leg->AddEntry(graphsAll[i], namesAll[i].c_str(), "p");
	}

	string outFileName = "FitEffVsMu-";
	outFileName.append(runName);
	outFileName.append(".txt");
	ofstream outfile;
	outfile.open(outFileName.c_str());

	vector<TF1*> functionV (graphsAll.size());
	for (unsigned int i=0; i < graphsAll.size(); i++) {
		functionV[i] = new TF1("EfficiencyVsMu", EfficiencyVsMu, 0.0, 10.0, 4);
	}
	// Fit all graphs in graphsAll to EfficiencyVsMu.
	outfile << "file number\tp1\tp0" << endl;
	for (unsigned int i=0; i < graphsAll.size(); i++) {
		Double_t par[] = {1.0, 1.0, 1.0, -1.0};
		functionV[i]->SetParameters(par);
		functionV[i]->FixParameter(3, 0.0); //*-* Only fit to exp(linear function).
		functionV[i]->FixParameter(2, 0.0); //*-* Only fit to exp(linear function).
		functionV[i]->SetLineColor(colors[i]);
		functionV[i]->SetParNames("p0", "p1", "p2", "p3");
		graphsAll[i]->Fit(functionV[i], "", "", 0.01, 0.19);

		double p1 = functionV[i]->GetParameter("p1");
		double p0 = functionV[i]->GetParameter("p0");
		outfile << i << "\t";//*-*
		outfile << p1 << "\t" << p0 << endl;
	}

	// Make log version of graphsAllLog
	for (unsigned int i=0; i < graphsAllLog.size(); i++) {
		//cout << graphsAllLogLog.size() << ": " << graphsAllLog[i] << endl;
		
		int n = graphsAllLog[i]->GetN();
		double* xV = graphsAllLog[i]->GetX();
		double* yV = graphsAllLog[i]->GetY();
		for (int j=0; j < n; j++) {
			yV[j] = log10( yV[j] );
		}
		double* exV = graphsAllLog[i]->GetEX();
		double* eyV = graphsAllLog[i]->GetEY();

		graphsAllLog[i] = new TGraphErrors(n, xV, yV, exV, eyV);
		graphsAllLog[i]->SetTitle("");//("Light collection efficiency vs Absorption Length");
		graphsAllLog[i]->GetXaxis()->SetTitle("Absorption coefficient [cm^{-1}]");
		graphsAllLog[i]->GetYaxis()->SetTitle("Log(Efficiency)");
		//graphsAllLog[i]->GetXaxis()->SetLimits(-0.005, 0.045);
		//graphsAllLog[i]->GetYaxis()->SetRangeUser(1e-6, 1e0);
		graphsAllLog[i]->SetMarkerColor(colors[i]);
		graphsAllLog[i]->SetMarkerSize(0.5);
		graphsAllLog[i]->SetMarkerStyle(21);

		//c1->cd(1);
		//if(i==0) {
		//	graphsAllLog[i]->Draw("AP");
		//}
		//else {
		//	graphsAllLog[i]->Draw("P same");
		//}
	}

	// Draw leg to c1 pad 1.
	{
		c1->cd(1);
		leg->SetFillColor(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.03);
		leg->SetTextFont(42);
		leg->Draw();
	}

	/*
	vector<TF1*> functionLogV (graphsAllLog.size(), new TF1("LogEfficiencyVsMu", LogEfficiencyVsMu, 0.0, 0.0, 3));
	// Fit all graphs in graphsAllLog to myfunction.
	for (unsigned int i=0; i < graphsAllLog.size(); i++) {
		Double_t par[] = {1.0, 1.0, 1.0};
		functionLogV[i]->SetParameters(par);
		functionLogV[i]->SetLineColor(colors[i]);
		functionLogV[i]->SetParNames("p0", "p1", "p2");
		//*-*graphsAllLog[i]->Fit(functionLogV[i], "", "", 0.01, 0.045);
	}

	// Read laser data file.
	// ratiosMap[ieta] contains ratios for ieta at luminosity values contained in lumis.
	vector<double> lumis;
	std::map<int, vector<double> > ratiosMapData;
	readLaserDataFile("../HELaserData.txt", lumis, ratiosMapData);

	// Read dose and fluence map data files.
	std::map<int, double> doseMapLayer1;
	std::map<int, double> fluenceMapLayer1;
	readDoseFile("../doseLayer1Mars.txt", doseMapLayer1, fluenceMapLayer1); //*-*

	// Variables to read and write from Trees.
	Int_t ieta, lumiIndex;
	Double_t lumi, ratio;
	ieta = 0;
	lumiIndex = 0;
	lumi = 0.;
	ratio = 0.;

	// Prepare to read dataRatios Tree.
	TFile* dataFile = new TFile("dataRatios.root");
	if (!dataFile->IsOpen()) {cout << "dataRatios.root not found"; return;}
	TTree* dataRatios = (TTree*)dataFile->Get("dataRatios");
	dataRatios->Show(10);
	dataRatios->SetBranchAddress("ieta", &ieta);
	dataRatios->SetBranchAddress("lumiIndex", &lumiIndex);
	dataRatios->SetBranchAddress("lumi", &lumi);
	dataRatios->SetBranchAddress("ratio", &ratio);
	
	// Prepare mcRatios Tree to be written. mcRatios shares branch addresses with dataRatios.
	TString mcFileName(runName);
	mcFileName.Prepend("mcRatios-").Append("-").Append(treeLabel).Append(".root");
	TFile* file = new TFile(mcFileName.Data(), "recreate");
	if (!file->IsOpen()) {cout << "mcRatios.root not found"; return;}
	TTree* mcRatios = new TTree("mcRatios", "ratios from Geant4");
	mcRatios->Branch("ieta", &ieta, "ieta/I");
	mcRatios->Branch("lumiIndex", &lumiIndex, "lumiIndex/I");
	mcRatios->Branch("lumi", &lumi, "lumi/D");
	mcRatios->Branch("ratio", &ratio, "ratio/D");

	for (unsigned int i=0; i < graphsAll.size(); i++) {
		// Define parameters for DoseToMu and calculate undamaged efficiency.
		//double par[2] = {0.02, 0.002};
		double doseRateBase = 2.67e-3;
		double A = 0.5;
		double par[2] = {0.023, 0.001};
		double par2[4] = {0.023, 0.001, doseRateBase, A};
		double muTileUndamaged = 0.023;

		for (unsigned int j=0; j < lumis.size(); j++) {
			dataRatios->GetEntryWithIndex(i+22, j); // Get entry with arbitrary ieta and lumiIndex == j.
			//stringstream cutFormula; cutFormula << "ieta == " << i << " && lumiIndex == " << j;
			//cout << cutFormula.str() << endl;
			//dataRatios->Draw("ieta:lumi:ratio:lumiIndex", cutFormula.str().c_str(), "goff");
			//ieta = (dataRatios->GetVal(1))[0];

			//cout << "ieta: " << ieta << endl; 
			//cout << "lumiIndex: " << lumiIndex << endl; 
			//cout << "lumi: " << lumi << endl; 
			//cout << "ratio: " << ratio << endl; 

			double dose = doseFactor * (doseMapLayer1[ieta] / 100) * lumi; // Dose[kGy] per fb-1.
			double doseRate = doseMapLayer1[ieta] * 3.6e-5; // Dose rate [kGy/hr]
			double xx[2] = {dose, doseRate};
			double muTile = DoseToMu(&dose, par);
			//double muTile = CalculateMu(xx, par2);
			//double efficiency = functionV[i]->Eval(muTile);
			double efficiency = graphsAll[i]->Eval(muTile); //*-*
			if(j==0) {
				muTileUndamaged = muTile;
			}
			//double efficiencyUndamaged = functionV[i]->Eval(muTileUndamaged);
			double efficiencyUndamaged = graphsAll[i]->Eval(muTileUndamaged); //*-*
			ratio = efficiency / efficiencyUndamaged;
			double ratioFromData = (ratiosMapData[ieta])[j];

//			if(j == lumis.size()-1) {
//				cout << "ieta: " << i+22 << endl;
//				cout << "dose: " << dose << endl; 
//				cout << "muTile: " << muTile << endl; 
//				cout << "efficiency: " << efficiency << endl; 
//				cout << "efficiencyUndamaged: " << efficiencyUndamaged << endl; 
//				cout << "ratio: " << ratio << endl; 
//				cout << "ratioFromData: " << ratioFromData << endl; 
//				cout << endl;
//			}
			mcRatios->Fill();
		}
	}

	int nEntries = dataRatios->GetEntries();
	for(int i=0; i < nEntries; i++) {
//		dataRatios->GetEntry(i);
//
//		cout << "ieta: " << ieta << endl; 
//		cout << "lumiIndex: " << lumiIndex << endl; 
//		cout << "lumi: " << lumi << endl; 
//		cout << "ratio: " << ratio << endl; 
	}
	
	mcRatios->Write();
*/
}

// Macro to produce rectangular tile fits.
void RectangularFits() 
{
	FitEffVsMu2("rectangular1", false);
	FitEffVsMu2("rectangular2", true);
}

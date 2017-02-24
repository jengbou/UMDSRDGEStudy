#include <iostream>

// Class to read Analysis output files from LYSim (run16 onwards).
// File format:
// outputfile << "#Mu_tile [cm^-1]\tMu_fiber [cm^-1]\tEfficiency" << G4endl;
class AnalysisOutput {

public:
	AnalysisOutput() {}
	// Main constructor with input file name
	AnalysisOutput(string inputFile)
	{
		cout << "AnalysisOutput: Reading new Analysis output file" << endl;
		cout << inputFile << endl;
	
		readFile(inputFile);
		fillGraphs();
	}
	~AnalysisOutput() {}

	vector<double> absorptionTileV;
	vector<double> absorptionFiberV;
	vector<double> efficiencyV;
	vector<double> errorV;

	TGraphErrors* graphTile; //Efficiency vs tile absorption coefficient
	TGraphErrors* graphFiber; //Efficiency vs fiber absorption coefficient

	void fillGraphs()
	{
		Int_t n = efficiencyV.size();
		graphTile = new TGraphErrors(n, &absorptionTileV[0], &efficiencyV[0], NULL, &errorV[0]);
		graphFiber = new TGraphErrors(n, &absorptionFiberV[0], &efficiencyV[0], NULL, &errorV[0]);
	}

	// Read analysis file and populate vectors
	void readFile(string inputFile)
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
					//Process the line to get absorption coefficients and efficiency
					double absorptionTile;
					double absorptionFiber;
					double efficiency;
					double errorEfficiency;
					processLine(line, absorptionTile, absorptionFiber, efficiency);
					absorptionTileV.push_back(absorptionTile);
					absorptionFiberV.push_back(absorptionFiber);
					efficiencyV.push_back(efficiency);
					// Binomial error from 1e5 events.
					errorEfficiency = sqrt(efficiency * (1 - efficiency) / 1e5);
					errorV.push_back(errorEfficiency);
				}
			}
		}
		else {
			cout<<"Error: Input file not found." <<endl;
		}
	}


private:

	// Process one line of the GEANT output. 
	void processLine(string line, double& absorptionTile, double& absorptionFiber, double& efficiency)
	{
		stringstream ss(line);
		string field;
		int ctr = 0; //column index of tab-delimited row
		while(getline(ss, field, '\t'))
		{
			if(ctr==0) {stringstream temp(field); temp >> absorptionTile;}
			else if (ctr==1) {stringstream temp(field); temp >> absorptionFiber;}
			else if (ctr==2) {stringstream temp(field); temp >> efficiency;}
			else {}
			
			ctr++;
		}
	}

};

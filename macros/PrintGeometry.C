
void printGeometry()
{
	// From CMSSW/Geometry/HcalTowerAlgo/src/HcalFlexiHardcodeGeometryLoader.cc.

	// Eta bounds for ieta 16 to 29.
	float etaBounds[] = {0.087*15, 0.087*16, 0.087*17, 0.087*18, 0.087*19,
	1.74, 1.83, 1.93, 2.043, 2.172,
	2.322, 2.500, 2.650, 2.868, 3.000};

	// Z-position for layers -1 to 17.
	float layerDepths[19] = {400.458, 408.718, 416.978, 425.248, 433.508, 
	441.768, 450.038, 458.298, 466.558, 474.828, 
	483.088, 491.348, 499.618, 507.878, 516.138, 
	524.398, 532.668, 540.928, 549.268};
	float angle1 = 0;
	float angle2 = 10 * atan(1)*4 / 180;

	int layerNo = 1;
	cout << "layer\tieta\tDy\tDx2" << endl;
	for (int ieta = 22; ieta< 30; ieta++)
	{
		float etaMin = etaBounds[ieta - 16];
		float etaMax = etaBounds[ieta - 16 + 1];

		float centerZ = layerDepths[layerNo - (-1) + 1];

		double thetaMin = 2 * atan(exp(-etaMax));
		double thetaMax = 2 * atan(exp(-etaMin));
		double rMin = centerZ * tan(thetaMin);
		double rMax = centerZ * tan(thetaMax);
		
		float Dy = rMax - rMin;
		float Dx2 = rMin * tan(angle2) - rMin * tan(angle1);
		
		//cout << "thetaMin set to " << thetaMin << endl;
		//cout << "thetaMax set to " << thetaMax << endl;
		//cout << "rMin set to " << rMin << endl;
		//cout << "rMax set to " << rMax << endl;
		//cout << "Dy set to " << Dy << endl;

		cout << layerNo << "\t" << ieta << "\t" << Dy << "\t" << Dx2 << endl;
	}
}

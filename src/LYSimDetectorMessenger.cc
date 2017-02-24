#include "LYSimDetectorMessenger.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

LYSimDetectorMessenger::LYSimDetectorMessenger(LYSimDetectorConstruction * Det)
 : Detector(Det)
{
	analysis = Analysis::GetInstance();

  detDir = new G4UIdirectory("/LYSim/");
  detDir->SetGuidance(" Geometry Setup ");

  UpdateCmd = new G4UIcmdWithoutParameter("/LYSim/Update",this);
  UpdateCmd->SetGuidance("Update geometry");
  UpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  UpdateCmd->SetGuidance("if you changed geometrical value(s).");
  UpdateCmd->AvailableForStates(G4State_Idle);
  
  SetFiberHoleCmd = new G4UIcmdWithABool("/LYSim/SetFiberHole", this);
  SetFiberHoleCmd->SetGuidance("Set Fiber Hole toggle true or false");
  SetFiberHoleCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  SetWrappingCmd = new G4UIcmdWithABool("/LYSim/SetWrapping", this);
  SetWrappingCmd->SetGuidance("Set Wrapping (Tyvek around scintillator) toggle true or false");
  SetWrappingCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
	
  SetFiberCmd = new G4UIcmdWithABool("/LYSim/SetFiber", this);
  SetFiberCmd->SetGuidance("Set Fiber toggle true or false");
  SetFiberCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  SetWLSCmd = new G4UIcmdWithABool("/LYSim/SetWLS", this);
  SetWLSCmd->SetGuidance("Set WLS toggle true or false");
  SetWLSCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  SetShieldingCmd = new G4UIcmdWithABool("/LYSim/SetShielding", this);
  SetShieldingCmd->SetGuidance("Set Shielding toggle true or false");
  SetShieldingCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  SetRefIndexCmd = new G4UIcmdWithADouble("/LYSim/SetRefIndex", this);
  SetRefIndexCmd->SetGuidance("Set the refractive index");
  SetRefIndexCmd->SetParameterName("RefIndex",false);
  SetRefIndexCmd->SetRange("RefIndex>=1.");
  SetRefIndexCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  SetScintThicknessCmd = new G4UIcmdWithADoubleAndUnit("/LYSim/SetScintThickness", this);
  SetScintThicknessCmd->SetGuidance("Set the scintillator thickness");
  SetScintThicknessCmd->SetParameterName("ScintThickness",false);
  SetScintThicknessCmd->SetUnitCategory("Length");
  SetScintThicknessCmd->SetDefaultUnit("mm");
  SetScintThicknessCmd->SetRange("ScintThickness>=0.");
  SetScintThicknessCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  SetScintSizeXYCmd = new G4UIcmdWithADoubleAndUnit("/LYSim/SetScintSizeXY", this);
  SetScintSizeXYCmd->SetGuidance("Set the scintillator transverse dimenstions");
  SetScintSizeXYCmd->SetParameterName("ScintSizeXY",false);
  SetScintSizeXYCmd->SetUnitCategory("Length");
  SetScintSizeXYCmd->SetDefaultUnit("mm");
  SetScintSizeXYCmd->SetRange("ScintSizeXY>=0.");
  SetScintSizeXYCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  SetScintPMTGapThicknessCmd = new G4UIcmdWithADoubleAndUnit("/LYSim/SetGapThickness", this);
  SetScintPMTGapThicknessCmd->SetGuidance("Set the thickness of the gap between the scintillator and PMT");
  SetScintPMTGapThicknessCmd->SetParameterName("GapThickness",false);
  SetScintPMTGapThicknessCmd->SetUnitCategory("Length");
  SetScintPMTGapThicknessCmd->SetDefaultUnit("mm");
  SetScintPMTGapThicknessCmd->SetRange("GapThickness>=0.");
  SetScintPMTGapThicknessCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  SetAngle1Cmd = new G4UIcmdWithADoubleAndUnit("/LYSim/SetAngle1", this);
  SetAngle1Cmd->SetGuidance("Set angle1 for tile. Conflicts with SetTileType.");
  SetAngle1Cmd->SetParameterName("angle1",false);
  SetAngle1Cmd->SetUnitCategory("Angle");
  SetAngle1Cmd->SetDefaultUnit("degree");
  SetAngle1Cmd->SetRange("angle1>=0.");
  SetAngle1Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  SetAngle2Cmd = new G4UIcmdWithADoubleAndUnit("/LYSim/SetAngle2", this);
  SetAngle2Cmd->SetGuidance("Set angle2 for tile. Conflicts with SetTileType.");
  SetAngle2Cmd->SetParameterName("angle2",false);
  SetAngle2Cmd->SetUnitCategory("Angle");
  SetAngle2Cmd->SetDefaultUnit("degree");
  SetAngle2Cmd->SetRange("angle2>=0.");
  SetAngle2Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  SetDx2Cmd = new G4UIcmdWithADoubleAndUnit("/LYSim/SetDx2", this);
  SetDx2Cmd->SetGuidance("Set Dx2 for tile");
  SetDx2Cmd->SetParameterName("Dx2",false);
  SetDx2Cmd->SetUnitCategory("Length");
  SetDx2Cmd->SetDefaultUnit("mm");
  SetDx2Cmd->SetRange("Dx2>=0.");
  SetDx2Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  SetDyCmd = new G4UIcmdWithADoubleAndUnit("/LYSim/SetDy", this);
  SetDyCmd->SetGuidance("Set Dy for tile");
  SetDyCmd->SetParameterName("Dy",false);
  SetDyCmd->SetUnitCategory("Length");
  SetDyCmd->SetDefaultUnit("mm");
  SetDyCmd->SetRange("Dy>=0.");
  SetDyCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  SetDzCmd = new G4UIcmdWithADoubleAndUnit("/LYSim/SetDz", this);
  SetDzCmd->SetGuidance("Set Dz for tile");
  SetDzCmd->SetParameterName("Dz",false);
  SetDzCmd->SetUnitCategory("Length");
  SetDzCmd->SetDefaultUnit("mm");
  SetDzCmd->SetRange("Dz>=0.");
  SetDzCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  SetIetaCmd = new G4UIcmdWithAnInteger("/LYSim/SetIeta", this);
  SetIetaCmd->SetGuidance("Set the ieta for the tile");
  SetIetaCmd->SetParameterName("ieta",false);
  SetIetaCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  SetLayerNoCmd = new G4UIcmdWithAnInteger("/LYSim/SetLayerNo", this);
  SetLayerNoCmd->SetGuidance("Set the layer number for the tile");
  SetLayerNoCmd->SetParameterName("layerNo",false);
  SetLayerNoCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
	
  SetLayerNoCmd = new G4UIcmdWithAnInteger("/LYSim/SetTileType", this);
  SetLayerNoCmd->SetGuidance("Set the tile type. Conflicts with SetAngle1 and SetAngle2.");
  SetLayerNoCmd->SetParameterName("type",false);
  SetLayerNoCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
	  
  SetTileAbsLengthCmd = new G4UIcmdWithADoubleAndUnit("/LYSim/SetTileAbsLength", this);
  SetTileAbsLengthCmd->SetGuidance("Set the light attenuation length in the tile");
  SetTileAbsLengthCmd->SetParameterName("TileAbsLength",false);
  SetTileAbsLengthCmd->SetUnitCategory("Length");
  SetTileAbsLengthCmd->SetDefaultUnit("cm");
  SetTileAbsLengthCmd->SetRange("TileAbsLength>=0.");
  SetTileAbsLengthCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  SetInducedMuTileCmd = new G4UIcmdWithADouble("/LYSim/SetInducedMuTile", this);
  SetInducedMuTileCmd->SetGuidance("Set the induced absorption coefficient (cm^-1)");
  SetInducedMuTileCmd->SetParameterName("InducedMuTile",false);
  SetInducedMuTileCmd->SetRange("InducedMuTile>=0.");
  SetInducedMuTileCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  SetInducedMuFiberCmd = new G4UIcmdWithADouble("/LYSim/SetInducedMuFiber", this);
  SetInducedMuFiberCmd->SetGuidance("Set the induced absorption coefficient (cm^-1)");
  SetInducedMuFiberCmd->SetParameterName("InducedMuFiber",false);
  SetInducedMuFiberCmd->SetRange("InducedMuFiber>=0.");
  SetInducedMuFiberCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

LYSimDetectorMessenger::~LYSimDetectorMessenger()
{
	delete detDir;
	delete UpdateCmd;
	delete SetFiberHoleCmd;
	delete SetWrappingCmd;
	delete SetFiberCmd;
	delete SetWLSCmd;
	delete SetShieldingCmd;
	delete SetRefIndexCmd;
	delete SetScintThicknessCmd;
	delete SetScintSizeXYCmd;
	delete SetScintPMTGapThicknessCmd;
	delete SetAngle1Cmd;
	delete SetAngle2Cmd;
	delete SetDx2Cmd;
	delete SetDyCmd;
	delete SetDzCmd;
	delete SetIetaCmd;
	delete SetLayerNoCmd;
	delete SetTileTypeCmd;
	delete SetTileAbsLengthCmd;
	delete SetInducedMuTileCmd;
	delete SetInducedMuFiberCmd;
}

void LYSimDetectorMessenger::SetNewValue(G4UIcommand* command,G4String val)
{
	if( command == UpdateCmd ) {
		Detector->UpdateGeometry();
	}
	else if( command == SetFiberHoleCmd ) {
		Detector->SetFiberHoleToggle(G4UIcmdWithABool::GetNewBoolValue(val));
	}
	else if( command == SetWrappingCmd ) {
		Detector->SetWrappingToggle(G4UIcmdWithABool::GetNewBoolValue(val));
	}
	else if( command == SetFiberCmd ) {
		Detector->SetFiberToggle(G4UIcmdWithABool::GetNewBoolValue(val));
	}
	else if( command == SetWLSCmd ) {
		Detector->SetWLSToggle(G4UIcmdWithABool::GetNewBoolValue(val));
	}
	else if( command == SetShieldingCmd ) {
		Detector->SetShieldingToggle(G4UIcmdWithABool::GetNewBoolValue(val));
	}
	else if( command == SetRefIndexCmd ) {
		Detector->SetRefIndex(G4UIcmdWithADouble::GetNewDoubleValue(val));
	}
	else if( command == SetScintThicknessCmd ) {
		Detector->
		SetScintThickness(G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(val));
	}
	else if( command == SetScintSizeXYCmd ) {
		Detector->
		SetScintSizeXY(G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(val));
	}
	else if( command == SetScintPMTGapThicknessCmd ) {
		Detector->
		SetScintPMTGapThickness(G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(val));
	}
	else if( command == SetAngle1Cmd ) {
		Detector->
		SetAngle1(G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(val));
	}
	else if( command == SetAngle2Cmd ) {
		Detector->
		SetAngle2(G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(val));
	}
	else if( command == SetDx2Cmd ) {
		Detector->
		SetDx2(G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(val));
	}
	else if( command == SetDyCmd ) {
		Detector->
		SetDy(G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(val));
	}
	else if( command == SetDzCmd ) {
		Detector->
		SetDz(G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(val));
	}
	else if( command == SetIetaCmd ) {
		Detector->
		SetIeta(G4UIcmdWithAnInteger::GetNewIntValue(val));
	}
	else if( command == SetLayerNoCmd ) {
		Detector->
		SetLayerNo(G4UIcmdWithAnInteger::GetNewIntValue(val));
	}
	else if( command == SetTileTypeCmd ) {
		Detector->
		SetTileType(G4UIcmdWithAnInteger::GetNewIntValue(val));
	}
	else if( command == SetTileAbsLengthCmd ) {
		Detector->
		SetTileAbsLength(G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(val));
	}
	else if( command == SetInducedMuTileCmd ) {
		G4double value = G4UIcmdWithADouble::GetNewDoubleValue(val);
		Detector->
		SetInducedMuTile(value);
		analysis->
		SetInducedMuTile(value);
	}
	else if( command == SetInducedMuFiberCmd ) {
		G4double value = G4UIcmdWithADouble::GetNewDoubleValue(val);
		Detector->
		SetInducedMuFiber(value);
		analysis->
		SetInducedMuFiber(value);
	}
}

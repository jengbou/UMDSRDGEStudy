================================================================================
The following instructions are for the lxplus cluster at CERN.


Checkout instructions
----------------------
1. Create working area
mkdir -p workarea/geant4work
cd workarea/geant4work

2. Checkout source code and copy to working area
git clone https://github.com/jengbou/UMDSRDGEStudy.git


Build and run instructions
-----------------------------
1. Create build directory

mkdir UMDSRDGEStudy-build

2a. Check out CMSSW (Only required the first time), e.g.
cmsrel CMSSW_7_6_3

2b. cd to SCRAM based area and set up runtime environment, e.g.
cd CMSSW_7_6_3/src
cmsenv

3. Compile the source code (only needed for one time)
cd ~/workarea/geant4work/UMDSRDGEStudy-build
export G4BASE=/afs/cern.ch/sw/lcg/external/geant4
source $G4BASE/10.3/x86_64-slc6-gcc49-opt/CMake-setup.sh
cmake -DWITH_GEANT4_UIVIS=ON -DGeant4_DIR=$G4BASE/10.3/x86_64-slc6-gcc49-opt/lib64/Geant4-10.3.0/ $HOME/workarea/geant4work/UMDSRDGEStudy
make

3a. Once the code is compiled, you can make a script for automatically setting up GEANT4 environment next time when you login lxplus as follows:
    3a.1: create a G4_setup.sh file:

======= G4_setup.sh =======
#!/bin/sh
G4BASE=/afs/cern.ch/sw/lcg/external/geant4/10.3/x86_64-slc6-gcc49-opt/bin
G4WORK=~/workarea/geant4work/UMDSRDGEStudy-build
cd ${G4BASE}
source geant4.sh
cd ${G4WORK}
===========================

    3a.2: add the following line to your ~/.bashrc file:
source ~/G4_setup.sh


================================================================================
Some examples:
4a. The following command run the code interactively for two different types of geometries.
    4a.1: for tile (Trapezoid):
    ./LYSim 0

    4a.2: for rod (5x1x1 cm^3):
    ./LYSim 1

4b. The following commands executes the input macro files (photontest.mac), and pipes the output to the corresponding output files (<log filename>.log)
./LYSim 1 photontest.mac <output filename> >& <log filename>.log &

Note: two files will be created: <output filename>.root and <output filename>.txt



See LYSimDetectorConstruction.cc for details on simulation setup, detector geometry and material properties

================================================================================

Description of files in LYSim project.

LYSim.cc
Main function of the LYSim program. For simulating light yield measurements.

LYSimPhysicsList
Activates required physics processes.

LYSimScintillation
Custom physics process. Inherits from G4Scintillation, with added output when scintillating.

LYSimPrimaryGeneratorAction
Initializes GeneralParticleSource. Sets random photon polarization. 
Energy and postion of generated photons are set in input file, as GPS commands,
e.g.
/gps/particle opticalphoton
/gps/pos/type Plane
/gps/pos/shape Circle
/gps/pos/centre 0. 0. 0. mm
/gps/pos/halfx 12.5 mm
/gps/pos/halfy 12.5 mm
/gps/ang/type iso
/gps/ene/type Arb
/gps/hist/type arb
/gps/hist/point 3.44e-6 0.00
/gps/hist/point 3.26e-6 0.06
/gps/hist/point 3.18e-6 0.28
/gps/hist/point 3.10e-6 0.72
/gps/hist/point 3.02e-6 1.40
/gps/hist/point 2.95e-6 2.00
/gps/hist/point 2.88e-6 2.20
/gps/hist/point 2.82e-6 2.06
/gps/hist/point 2.70e-6 1.48
/gps/hist/point 2.58e-6 0.94
/gps/hist/point 2.48e-6 0.60
/gps/hist/point 2.38e-6 0.40
/gps/hist/point 2.30e-6 0.30
/gps/hist/point 2.21e-6 0.20
/gps/hist/point 2.14e-6 0.10
/gps/hist/point 1.82e-6 0.00
/gps/hist/inter Lin

LYSimRunAction
Prepares new run and collects analysis information for a run in the Analysis code.
Calls Analysis::PrepareNewRun() and Analysis::EndOfRun() before and after each run.

LYSimEventAction
Prints event ID every 100 events.
Calls Analysis::PrepareNewEvent() and Analysis::EndOfEvent() before and after each event.

LYSimTrackingAction/LYSimTrajectory/LYSimTrajectoryPoint
Specifies track information to save. (Default should be good enough for now.)
Can also be used to control which tracks to draw and color of tracks. (Currently draws everything.)

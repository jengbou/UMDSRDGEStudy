#!/bin/bash

#
# variables from arguments string in jdl
#

RUN_DIR=$1
INPUT_NAME=$2
ANALYSIS_FILE_NAME=$3
DET_TYPE=$4

echo $RUN_DIR
echo $INPUT_NAME
echo $ANALYSIS_FILE_NAME
echo $DET_TYPE

bash

echo ""
echo "CMSSW on Condor"
echo ""

START_TIME=`/bin/date`
echo "started at $START_TIME"

CMSSWBASE=/home/jengbou/workspace/CMSSW_7_6_3/src/
#
# setup CMSSW software environment at UMD
#
export VO_CMS_SW_DIR=/sharesoft/cmssw
. $VO_CMS_SW_DIR/cmsset_default.sh
cd ${CMSSWBASE}
eval `scramv1 runtime -sh`

#G4BASE=/data/users/jengbou/workspace/UserCode/geant4.10.03-install/bin
G4BASE=/data/users/avermeer/geant4/geant4.9.6.p03-build/bin

cd ${G4BASE}
. geant4.sh

cd $RUN_DIR

./LYSim ${DET_TYPE} ${INPUT_NAME} ${ANALYSIS_FILE_NAME}

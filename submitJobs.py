#!/usr/bin/python
import os, sys
import shlex, subprocess
from datetime import datetime, date, time
#sys.path.append(os.path.abspath(os.path.curdir))

JobTime = datetime.now()
fTag = JobTime.strftime("%Y%m%d_%H%M%S")
sTag = "AlphaSource_1kevts"
dirname = "jobs/%s_%s"%(sTag,fTag)
InFile = "AlphaSource.in"#"HERadDam.in"
DetType = "1" #rod

try:
    os.makedirs(dirname)
except:
    pass

ProdTag = "UMDSRDGEStudy-00-00-02_20170331"
OutDir  = "/data/users/jengbou/G4Results"
WorkDir = "/home/jengbou/workspace/UserCode/geant4/UMDSRDGEStudy-build"

try:
    os.makedirs(OutDir)
except:
    pass

try:
    os.makedirs(OutDir+"/"+ProdTag)
except:
    pass

try:
    os.makedirs(OutDir+"/"+ProdTag+"/logs")
except:
    pass


#########################################
# make sure OutDir is the same in main.cc
#########################################
condor_script_template = """
universe = vanilla
Executable = ./LYSim.sh
+IsLocalJob = true
Should_transfer_files = NO
Requirements = TARGET.FileSystemDomain == "privnet"
Output = %(OUTDIR)s/%(MYPREFIX)s/logs/%(FILENAME)s_sce_$(cluster)_$(process).stdout
Error  = %(OUTDIR)s/%(MYPREFIX)s/logs/%(FILENAME)s_sce_$(cluster)_$(process).stderr
Log    = %(OUTDIR)s/%(MYPREFIX)s/logs/%(FILENAME)s_sce_$(cluster)_$(process).condor
Arguments = %(WORKDIR)s %(INPUT)s %(OUTDIR)s/%(MYPREFIX)s/%(FILENAME)s %(DETTYPE)s
Queue 1

"""
#########################################

kw = {}

kw["MYPREFIX"]  = ProdTag
kw["WORKDIR"]   = WorkDir
kw["OUTDIR"]    = OutDir
kw["INPUT"]     = InFile
kw["FILENAME"]  = sTag
kw["DETTYPE"]   = DetType

script_str = condor_script_template % kw
f = open("%s/condor_jobs_%s_G4Sim.jdl"%(dirname,sTag), 'w')
f.write(script_str)
f.close()

condorcmd = "condor_submit %s/condor_jobs_%s_G4Sim.jdl"%(dirname,sTag)
print 'condorcmd: ', condorcmd
print 'Executing condorcmd'
p=subprocess.Popen(condorcmd, shell=True)
p.wait()

print "Histos output dir: %s/%s"%(OutDir,ProdTag)


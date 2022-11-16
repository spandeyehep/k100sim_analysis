import os, sys
import uproot, awkward
#import ROOT
import numpy as np
from array import array
from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt
from optparse import OptionParser
import subprocess

import utils.randomize_nrCascadeSim as nrcRand
import utils.helper as helper

parser = OptionParser()
parser.add_option("-d", "--basedir", dest="basedir",action="store",
                      help="base directory containing k100sim step1 root file",default="/Users/shubhampandey/work/geant4/k100sim_anthony/sim_files/sim_2M_neutrons")
parser.add_option("-f", "--inFileName", dest="inFileName",action="store",
                      help="Only neutron capture k100sim and nrCascadeSim step2 root file",default="sim_2M_events_PuBe_nrCascade_step2.root")

(options, args) = parser.parse_args()
print(options)

base_dir = options.basedir

inFileName = base_dir+"/"+options.inFileName
inFile = uproot.open(inFileName)
if(not inFile):
	print("could not open file: %s"%(inFileName))
	sys.exit(0)
inTree = inFile["combined"]
if(not inTree):
	print("combined does not exist in file: %s"%(inFileName))
	sys.exit(0)

event = inTree['Event'].array(library="np")
COUNT = len(event)
print("] Number of events to be produced : %d"%(COUNT))

basedir_geant4 = "/Users/shubhampandey/work/geant4/k100sim_anthony/k100Sim_ROOT_OUT_nrCascade_combined_NaI_array/build"

outMacFileName = options.inFileName.strip(".root") + ".mac"
#tempName = options.inFileName
FullMacroFileName = basedir_geant4+"/macros/"+outMacFileName
#print (outMacFileName)
# print ((tempName).strip("_step2.root"))
# print (tempName.strip("_step2.root") + ".mac")

#sys.exit(0)
outMacFile = open(FullMacroFileName,"w")
outMacFile.write("/CDMS/detector/activate Shields\n")
outMacFile.write("/CDMS/detector/activate IceBox\n")
outMacFile.write("/CDMS/detector/activate NaIArray\n")
outMacFile.write("/CDMS/PuBe/doPuBeGammas false\n")
outMacFile.write("/CDMS/updateGeometry\n")
outMacFile.write("/run/beamOn %d\n"%(COUNT))
#outMacFile.write("exit\n")
outMacFile.close()
print("] mac file: %s written for %d events"%(outMacFileName,COUNT))

executable = basedir_geant4 + "/k100Sim"
if(not os.path.exists(executable)):
	print ("ERROR: k100sim executable not found at %s"%(basedir_geant4))
	print ("Exiting.")
	sys.exit(0)

cmd = executable + " -c -i " + inFileName + " " +  FullMacroFileName
os.system(cmd)

# temp_shell_file = open("temp_shell_file.sh","w")
# temp_shell_file.write(cmd)
# os.system("source temp_shell_file.sh")

#print (cmd.split(" "))
# subprocess.run([executable,"-c","-i",inFileName,FullMacroFileName])
# print ("] Running : %s"%(cmd)) 
# print (os.system(cmd))


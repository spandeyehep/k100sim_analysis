import os, sys
import uproot, awkward
#import ROOT
import numpy as np

def getNcapIndex(list_):
	for i,it in enumerate(list_):
		if(it == 1):
			return i
	return -1

def GetNcapEvents(nCap,DT):
	nevents = 0
	for i,list__ in enumerate(nCap):
		ncap_index = getNcapIndex(list__)
		if(ncap_index >= 0 and DT[i][ncap_index] == 1):
			nevents += 1
	return nevents

def generate_nrCascadeSim(events,outfilename):
	executable = "/Users/shubhampandey/work/NuclearCascade/nrCascadeSim_build/realizeCascades"
	print ("] Using executable = %s"%(executable))
	if(not os.path.exists(executable)):
		print ("] Executable not found. Please make sure file name and path is correct.")
		print ("] Exiting.")
		sys.exit(0)
	cmd = executable + " -n " + str(events) + " -o " + outfilename + " /Users/shubhampandey/work/NuclearCascade/nrCascadeSim/levelfiles/v3_natSi.txt"
	print ("\n \n Running : "+ cmd + "\n \n" )
	os.system(cmd)
	print ("] File generated: %s \n"%(outfilename))
	if(not os.path.exists(outfilename)):
		print ("something went wrong in nrCascadeSim file generation.")
		sys.exit(0)


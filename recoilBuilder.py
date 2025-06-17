import ROOT
import sys
import numpy as np
import uproot as up
import pandas as pd
import array
from utils import collisionUtils
import os

#########
##INPUT##
#########
#Check arguments
if not len(sys.argv)==3:
  print("\nERROR! Usage is:\n\n   python3 recoilBuilder.py <parsed PTRAC file name> <desired output name>\n\nExiting!\n")
  sys.exit(1)

#Open input file, load tree
inpFilename = sys.argv[1]
branchNames = ["historyNum","eventType","targetNucleus","interactionType","node","particleType","x","y","z","u","v","w","energy","time"]
inpFile = up.open(inpFilename)
tree = inpFile["ptracTree"]
df = tree.arrays(branchNames,library="pd")

#Get unique set of events
allHistoryNums = df["historyNum"][:]
uniqueHistoryNums = list(set(allHistoryNums))
uniqueHistoryNums.sort()
print("Found {0} unique histories in file".format(len(uniqueHistoryNums)))

##########
##OUTPUT##
##########
#Open output file, prep output tree
outFilename = sys.argv[2]
outFile = ROOT.TFile(outFilename,"RECREATE")
tree = ROOT.TTree("recoilTree","Events")

#Tree branch variables
historyNum = array.array('i',[0])
pdgCode = array.array('i',[0])
x = array.array('d',[0.])
y = array.array('d',[0.])
z = array.array('d',[0.])
u = array.array('d',[0.])
v = array.array('d',[0.])
w = array.array('d',[0.])
energy = array.array('d',[0.])
time = array.array('d',[0.])
interactionType = array.array('i',[0])

tree.Branch("historyNum", historyNum, "historyNum/I")
tree.Branch("interactionType",interactionType,"interactionType/I")
tree.Branch("pdgCode", pdgCode, "pdgCode/I")
tree.Branch("x", x, "x/D")
tree.Branch("y", y, "y/D")
tree.Branch("z", z, "z/D")
tree.Branch("u", u, "u/D")
tree.Branch("v", v, "v/D")
tree.Branch("w", w, "w/D")
tree.Branch("energy", energy, "energy/D")
tree.Branch("time", time, "time/D")

#############
##MAIN LOOP##
#############
numBadEvents=0
#Loop through histories.
for ihist,histNum in enumerate(uniqueHistoryNums):
  if ihist%1000==0:
    print("On history {0} of {1}".format(ihist,len(uniqueHistoryNums)))

  #Create subset of data for this throw
  history_df = df[df["historyNum"] == histNum]
  
  #Get the neutron collisions
  neutronCollisions_df = history_df[(history_df["eventType"]==4000) & (history_df["particleType"]==1)]

  #Step through neutron collisions:
  for _, row in neutronCollisions_df.iterrows():
    #Get all neutron surface crossings and collisions prior to this collisiion
    potentialParents_df = history_df[((history_df["eventType"]==3000) | (history_df["eventType"]==4000)) & 
                                     (history_df["time"] < row["time"]) & 
                                     (history_df["particleType"]==1)] 
    
    #Get all banked particles
    bankedParticles_df = history_df[(history_df["eventType"]>=2000) & (history_df["eventType"] < 3000) &
                                   (history_df["time"]==row["time"]) &
                                   (history_df["x"]==row["x"]) & (history_df["y"] == row["y"]) & (history_df["z"] == row["z"])]
    
    ions = collisionUtils.processNeutronCollision(row,potentialParents_df,bankedParticles_df,histNum)
    if len(ions)==0:
      numBadEvents+=1

    #Fill tree
    badEv=0
    for ion in ions:
      historyNum[0] = histNum
      pdgCode[0] = ion.pdgCode
      x[0] = ion.x
      y[0] = ion.y
      z[0] = ion.z
      u[0] = ion.u
      v[0] = ion.v
      w[0] = ion.w
      energy[0] = ion.KE
      time[0] = ion.time
      interactionType[0] = int(row["interactionType"])
      if energy[0]<0:
        badEv=1
      tree.Fill()
      
    if badEv==1:
      for ion in ions:
        print(histNum,interactionType[0],ion.KE)

print(f"Out of {len(uniqueHistoryNums)} events, found {numBadEvents} that we could not process")
#Write to file
outFile.cd()
tree.Write("eventTree",ROOT.TObject.kOverwrite)
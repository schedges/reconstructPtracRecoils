import ROOT
import sys
import numpy as np
import uproot as up
import pandas as pd

#######
##I/O##
#######
#Check arguments
if not len(sys.argv)>=2:
  print("\nERROR! Usage is:\n\n python3 historyPrinter.py <name> [optional history number]\n\nExiting!\n")
  sys.exit(1)
if len(sys.argv)==3:
  historyToPrint = int(sys.argv[2])
else:
  historyToPrint=0

#Open input file, load tree
inpFilename = sys.argv[1]
branchNames = ["historyNum","eventType","targetNucleus","interactionType","node","particleType","x","y","z","u","v","w","energy","weight","time"]
inpFile = up.open(inpFilename)
tree = inpFile["ptracTree"]
df = tree.arrays(branchNames,library="pd")

#################
##PRE-LOOP PREP##
#################
#Get unique set of events
if historyToPrint==0:
  allHistoryNums = df["historyNum"][:]
  uniqueHistoryNums = list(set(allHistoryNums))
  uniqueHistoryNums.sort()
  print("Found {0} unique histories in file".format(len(uniqueHistoryNums)))
else:
  uniqueHistoryNums=[historyToPrint]

#############
##MAIN LOOP##
#############
#Loop through histories.
for historyNum in uniqueHistoryNums:

  filtered_df = df[df["historyNum"] == historyNum]
  checkIfAllSurface_df = filtered_df[filtered_df["eventType"]!=3000]
  if len(checkIfAllSurface_df)>0:

    lines = filtered_df.to_dict(orient="records")

    hasIntType=1 #Set to zero here to enable filtering
    for i,line in enumerate(lines):
      if line.get("interactionType")>100:
        hasIntType=1 

    if hasIntType==1:
      print("\nHistory #{0}".format(historyNum))
      for i,line in enumerate(lines):
        print("Line={0}, eventType={1}, targetNucleus={10:05d}, intType={2:04d}, partType={3:02d}, energy={4:.4f}, time={5:.4f}, x={6:.8f}, y={7:.8f}, z={8:.8f}, weight={9:.3f}, node={11:02d}, u={12:0.3f}, v={13:.3f}, w={14:.3f}".format(
          i,line.get("eventType"),line.get("interactionType"),line.get("particleType"),line.get("energy"),line.get("time"),line.get("x"),line.get("y"),line.get("z"),line.get("weight"),line.get("targetNucleus"),line.get("node"),line.get("u"),line.get("v"),line.get("w")))
    
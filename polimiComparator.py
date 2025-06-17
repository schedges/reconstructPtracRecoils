import ROOT
import sys
import uproot as up
import matplotlib.pyplot as plt
import subprocess
import os
import subprocess

print_events = 0

ptracFilename = sys.argv[1]
polimiFilename = sys.argv[2]

ptracBranches = ["historyNum", "interactionType", "pdgCode", "x", "y", "z", "u", "v", "w", "energy", "time"]
ptracFile = up.open(ptracFilename)
ptracTree = ptracFile["eventTree"]
ptrac_df = ptracTree.arrays(ptracBranches, library="pd")

polimiBranches = ["historyNum", "particleNum", "lineNum", "projectileType", "interactionType", "targetNucleus",
                  "cellNum", "channelID", "energy", "timestamp", "xLoc", "yLoc", "zLoc", "particleWeight",
                  "generationNum", "numScatterings", "code", "energyPreScatter"]
polimiFile = up.open(polimiFilename)
polimiTree = polimiFile["polimiTree"]
polimi_df = polimiTree.arrays(polimiBranches, library="pd")
# Apply filtering: adjust values below as needed

projectile_mask = polimi_df["projectileType"] == 1

# Now apply mask
polimi_df = polimi_df[projectile_mask].reset_index(drop=True)

#################
## PRE-LOOP PREP
#################
polimiEnergies = []
ptracEnergies = []
ratios = []

allHistoryNums = polimi_df["historyNum"][:]
uniqueHistoryNums = list(set(allHistoryNums))
uniqueHistoryNums.sort()
print("Found {0} unique histories in file".format(len(uniqueHistoryNums)))

for ihist, historyNum in enumerate(uniqueHistoryNums):
    if ihist % 1000 == 0:
        print("On history {0} of {1}".format(ihist, len(uniqueHistoryNums)))

    filteredPolimi_df = polimi_df[polimi_df["historyNum"] == historyNum]
    filteredPtrac_df = ptrac_df[ptrac_df["historyNum"] == historyNum]

    totalPolimiEnergy = filteredPolimi_df["energy"].sum()
    totalPtracEnergy = filteredPtrac_df["energy"].sum()

    polimiEnergies.append(totalPolimiEnergy)
    ptracEnergies.append(totalPtracEnergy)

    if totalPolimiEnergy > 0:
        ratios.append(totalPtracEnergy / totalPolimiEnergy)

    if print_events==1:
      
      print("{0}: ratio={1:.3f}, ptrac={2:.5f}, polimi={3:.5f}".format(historyNum,ratios[-1],totalPtracEnergy,totalPolimiEnergy))
      if abs(ratios[-1] - 1) > 0.10 and totalPtracEnergy>0.002:
        try:
          subprocess.run(["python3", "historyPrinter.py", "data/ptrac/ptrac.root", str(int(historyNum))], check=True)
        except KeyboardInterrupt:
          print("\nInterrupted during subprocess — exiting.")
          sys.exit(0)
        print()
      
          
# Energy distributions
plt.figure()
plt.hist(polimiEnergies, label="polimi", bins=500, range=(0, 10), histtype='step')
plt.hist(ptracEnergies, label="ptrac", bins=500, range=(0, 10), histtype='step')
plt.xlabel("Energy (MeV)")
plt.ylabel("Counts")
plt.yscale("log")
plt.legend()
plt.title("Total Energy per History")
plt.tight_layout()

# Ratio distribution
plt.figure()
plt.hist(ratios, bins=400, range=(0, 2), histtype='step', color='black')
plt.yscale("log")
plt.xlabel("PTRAC / PoliMi Total History Energy")
plt.ylabel("Counts")
plt.tight_layout()

plt.show()

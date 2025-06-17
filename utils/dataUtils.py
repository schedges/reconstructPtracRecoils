import numpy as np
import math
import os

m_e = 0.51099895069 #MeV/c
amu_to_MeV = 931.49410242

#Loads the masses for calculating thresholds.
def loadMassTable(massTableFilename="data/mass_1.mas20.txt"):
  #Uses Chinese Phys. C 45 030002 (2021), and Chinese Phys. C 45 030003 (2021) as reported on https://www-nds.iaea.org/amdc/ to 
  #load up masses into a table Z,A,mass in microamu. We pad this table so we can index it directly from Z,A instead of having to worry about zero-indexing.
  massTable=np.zeros((300,300))
  firstLineNum=36
  for iline,line in enumerate(open(massTableFilename,"r")):
    if iline>=firstLineNum:
      if not line=="":
        line=line.strip("\n")
        line=line.replace("#","")
        lineParts=line.split()
        nLeadingSpaces = len(line) - len(line.lstrip(" "))
        if nLeadingSpaces==0:
          Z=int(lineParts[3])
          A=int(lineParts[4])
        else:
          Z=int(lineParts[2])
          A=int(lineParts[3])
        atomic_MeV=(float(lineParts[-3])*1000000+float(lineParts[-2]))*math.pow(10,-6)*amu_to_MeV
        nuclear_MeV = atomic_MeV - Z*m_e
        massTable[Z][A]=nuclear_MeV 
  return massTable

levelsFolderName = "data/levels/"
if not levelsFolderName.endswith("/"):
  levelsFolderName+="/"
  levelsFiles = [i for i in os.listdir(levelsFolderName) if i.startswith("z") and i.endswith(".dat")]
levelsFiles = []
levelsDict = {}
def createLevelsDict():
  levelsDict = {}
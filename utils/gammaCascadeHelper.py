#Loads RIPL 3 data, build up a library of level objects. We will call this and then
#randomly assign directions to the cascade gammas (for now). We'll sum up the total momentum
#contribution from these gammas, and then use that in our inelastic calculator to 
#determine the reocoil direction.
import numpy as np
import math
import os

levelsFolder = "data/levels/"
if not levelsFolder.endswith("/"):
  levelsFolder+="/"
levelsFiles = [i for i in os.listdir(levelsFolder) if i.startswith("z") and i.endswith(".dat")]

#Defines an energy level--the energy, parity, spin, #of decays, and decay scheme (list of decays)
class Level:
  def __init__(self,energy,twoJ,parity,nDecays,decayScheme):
    self.energy = energy
    self.twoJ = twoJ
    self.parity = parity
    self.nDecays = nDecays
    self.decayScheme = decayScheme #list of decay objects
           
#Decay object, has an ennergy, branching ratio, and identifier of the final state
class Decay:
  def __init__(self,energy,branchingRatio,finalStateNumber):
    self.energy=energy
    self.branchingRatio=branchingRatio
    self.finalStateNumber=finalStateNumber

class Gamma:
  def __init__(self,energy,u,v,w):
    self.energy = energy
    self.u = u
    self.v = v
    self.w = w

def findHighestLevelIndexWithDecayScheme(levels):
  nLevels = len(levels)
  for i in range(nLevels-1,0,-1):
    if len(levels[i].decayScheme)>0:
      return i
  return 0

def findHighestLevelIndexBeforeEnergy(levels,energy):
  nLevels = len(levels)
  for i in range(nLevels-1,0,-1):
    if levels[i].energy < energy:
      return i
  return 0

#NOTE: The level identifiers are NOT zero-indexed (i.e. ground state = 1) in 
#RIPL data, so we have to subtract 1 when referring to indices in our list.
def loadLevelsData(Z,A):
  levels = {}

  fname = "z{0:03d}.dat".format(Z)
  if not fname in levelsFiles:
    print(f"Error! file {fname} not found")
    return levels
  
  #Open input file
  levelDataFile = open(levelsFolder+fname,"r")
  lines = levelDataFile.readlines()

  #Get the start name of the line
  startString = "{0}{1}".format(A,atomicSymbolMap.get(Z))

  #Loop through lines
  startLineNum=-1
  for iline,line in enumerate(lines):
    #Skip blanks
    if line=="":
      continue
    #Skip until we find the first line with our start string
    elif line.lstrip().startswith(startString):
      startLineNum=iline
      break
  
  #If we didn't find this isotope, return empty list of evels
  if startLineNum==-1:

    print(f"Error in loadLevels!")
    return levels
  
  line=lines[startLineNum]
  nLevels = int(line[15:20])
  
  currentLine=startLineNum+1
  for ilev in range(0,nLevels):
    line=lines[currentLine]
    #idx = int(line[0:3])
    energy = float(line[4:14])
    twoJ = int(float(line[15:20])*2)
    parity = int(line[21:24])
    nDecays = int(line[35:38])
    decayScheme = []
    currentLine+=1 #Advance to decay scheme
    for idecay in range(0,nDecays):
      line=lines[currentLine]
      decayIndex = int(line[39:43])
      decayEnergy = float(line[44:54])
      branchingRatio = float(line[55:65])
      decayScheme.append(Decay(decayEnergy,branchingRatio,decayIndex)) #TODO we're neglecting non-gamma decay channels
      currentLine+=1
    levels[ilev] = Level(energy,twoJ,parity,nDecays,decayScheme)
  
  return levels

def generateBankedGammas(excitedStateIndex,levels):
  bankedGammas=[]
  while excitedStateIndex > 0:
    decayScheme = levels[excitedStateIndex].decayScheme
    #If we don't have a decay scheme, emit a gamma with energy to put us on the cascade
    if len(decayScheme)==0:
      maxLevelIdx = findHighestLevelIndexWithDecayScheme(levels)
      energy = levels[excitedStateIndex].energy - levels[maxLevelIdx].energy
      excitedStateIndex = maxLevelIdx
    else:
      probs=[]
      emissionEnergies=[]
      for decay in decayScheme:
        probs.append(decay.branchingRatio)
        emissionEnergies.append(decay.energy)
      probs = np.array(probs)
      probs /= probs.sum() #Force normalization to 1
      decay = np.random.choice(decayScheme,p=probs)
      energy = decay.energy
      excitedStateIndex = decay.finalStateNumber - 1 #ACCOUNT FOR INDEXING MISMATCH HERE

    #Give random direction. #TODO we might want to calculate a direction manually
    theta = np.arccos(2*np.random.uniform() - 1)
    phi = 2 * np.pi * np.random.uniform()
    u = np.sin(theta) * np.cos(phi)
    v = np.sin(theta) * np.sin(phi)
    w = np.cos(theta)
    mag = math.sqrt(u*u+v*v+w*w)
    u = u/mag
    v = v/mag
    w = w/mag
    bankedGammas.append(Gamma(energy,u,v,w))
  return bankedGammas


atomicSymbolMap = {
    1: "H",
    2: "He",
    3: "Li",
    4: "Be",
    5: "B",
    6: "C",
    7: "N",
    8: "O",
    9: "F",
    10: "Ne",
    11: "Na",
    12: "Mg",
    13: "Al",
    14: "Si",
    15: "P",
    16: "S",
    17: "Cl",
    18: "Ar",
    19: "K",
    20: "Ca",
    21: "Sc",
    22: "Ti",
    23: "V",
    24: "Cr",
    25: "Mn",
    26: "Fe",
    27: "Co",
    28: "Ni",
    29: "Cu",
    30: "Zn",
    31: "Ga",
    32: "Ge",
    33: "As",
    34: "Se",
    35: "Br",
    36: "Kr",
    37: "Rb",
    38: "Sr",
    39: "Y",
    40: "Zr",
    41: "Nb",
    42: "Mo",
    43: "Tc",
    44: "Ru",
    45: "Rh",
    46: "Pd",
    47: "Ag",
    48: "Cd",
    49: "In",
    50: "Sn",
    51: "Sb",
    52: "Te",
    53: "I",
    54: "Xe",
    55: "Cs",
    56: "Ba",
    57: "La",
    58: "Ce",
    59: "Pr",
    60: "Nd",
    61: "Pm",
    62: "Sm",
    63: "Eu",
    64: "Gd",
    65: "Tb",
    66: "Dy",
    67: "Ho",
    68: "Er",
    69: "Tm",
    70: "Yb",
    71: "Lu",
    72: "Hf",
    73: "Ta",
    74: "W",
    75: "Re",
    76: "Os",
    77: "Ir",
    78: "Pt",
    79: "Au",
    80: "Hg",
    81: "Tl",
    82: "Pb",
    83: "Bi",
    84: "Po",
    85: "At",
    86: "Rn",
    87: "Fr",
    88: "Ra",
    89: "Ac",
    90: "Th",
    91: "Pa",
    92: "U",
    93: "Np",
    94: "Pu",
    95: "Am",
    96: "Cm",
    97: "Bk",
    98: "Cf",
    99: "Es",
    100: "Fm",
    101: "Md",
    102: "No",
    103: "Lr",
    104: "Rf",
    105: "Db",
    106: "Sg",
    107: "Bh",
    108: "Hs",
    109: "Mt",
    110: "Ds",
    111: "Rg",
    112: "Cn",
    113: "Nh",
    114: "Fl",
    115: "Mc",
    116: "Lv",
    117: "Ts",
    118: "Og"
}

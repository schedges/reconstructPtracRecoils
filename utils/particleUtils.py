import math 
import numpy as np
import sys
from utils import dataUtils

massTable = dataUtils.loadMassTable()

#############################################
#Codes for converting between naming schemes#
#############################################
def calcZA_from_pdgCode(pdgCode):
  Z,A = -1,-1
  pdgCodeStr = str(pdgCode)
  if pdgCode==2112:
    Z=0
    A=1
  elif pdgCode==2212:
    Z=1
    A=1
  elif pdgCode==22:
    Z=0
    A=0
  else:
    Z = int(pdgCodeStr[2:6])
    A = int(pdgCodeStr[6:9])
  return int(Z),int(A)

def calcZA_from_polimiCode(polimiCode):
  A = int(polimiCode%1000)
  Z = int(polimiCode-A)/1000
  if Z>1 and A==0:
    print(f"Error! Invalid polimiCode {polimiCode}")
    sys.exit()
  return int(Z),int(A)

def calcPolimiCode_from_ZA(Z,A):
  return 1000*Z + A

def calcPDGCode_from_IPT(ipt):
  pdgDict = {
    1: 2112,
    2: 22,
    9: 2212,
    31: 1000010020,
    32: 1000010030,
    33: 1000020030,
    34: 1000020040
  }
  if ipt in pdgDict:
    return pdgDict[ipt]
  else:
    print(f"ERROR! Could not determine PDG code for IPT={ipt}")
    sys.exit


def calcPolimiIPT_from_pdgCode(pdg):
  iptDict = {
    2112: 1,
    22: 2,
    2212: 9, 
    1000010020: 31,
    1000010030: 32,
    1000020030: 33,
    1000020040: 34,
  }
  if pdg in iptDict:
    return iptDict[pdg]
  else:
    print(f"ERROR! Could not determine IPT from PDG code {pdg}")
    sys.exit()

def calcPDGCode_from_ZA(Z,A):
  if Z==1 and A==1:
    return 2212
  elif Z==0 and A==1:
    return 2112
  elif Z==0 and A==0:
    return 22
  else:
    return int("100{0:03d}{1:03d}0".format(Z,A))

######################
##PARTICLE INSTANCES##
######################
class Particle:
  def __init__(self,Z,A,x,y,z,u,v,w,KE,time): 
    self.Z = Z
    self.A = A
    self.mass = massTable[Z][A]
    self.pdgCode = calcPDGCode_from_ZA(Z,A)

    self.x = x
    self.y = y
    self.z = z

    self.KE = KE
    self.time = time

    self.u = u
    self.v = v
    self.w = w
    self.fourVec = None

    norm = math.sqrt(self.u**2 + self.v**2 + self.w**2)
    if norm > 0:
        self.u /= norm
        self.v /= norm
        self.w /= norm

  def create4Vec(self):
    E = self.KE + self.mass
    pMag = math.sqrt(self.KE * (self.KE + 2 * self.mass))
    px, py, pz = pMag * np.array([self.u, self.v, self.w])
    self.fourVec = np.array([E, px, py, pz])


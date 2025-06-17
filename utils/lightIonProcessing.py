import math
from utils import particleUtils
import sys
import numpy as np

def calcQValue(inputPDGs,exitPDGs):
  inputSum,exitSum = 0,0
  for pdg in inputPDGs:
    Z,A = particleUtils.calcZA_from_pdgCode(pdg)
    inputSum += particleUtils.massTable[Z][A]
  for pdg in exitPDGs:
    Z,A = particleUtils.calcZA_from_pdgCode(pdg)
    exitSum += particleUtils.massTable[Z][A]
  return exitSum - inputSum


def handleLightIons(bankedLightIons,missingLightIonPDGCodes,
                    incidentPDGCodes,expectedBankPDGCodes, #For calculating Q values
                    incidentParticles, collision, residualNucleusPDGCode, #For calculations
                    outgoingParticles, #For updating in place
                    historyNum): #for debugging

  #if banked is two and missing is zero, this is a (n,alpha t) type rxn that polimi returns correctly. Return the bank
  if len(bankedLightIons)==2 and len(missingLightIonPDGCodes)==0:
    for bankedLightIon in bankedLightIons:
      outgoingParticles.append(bankedLightIon)
  
  ###########################
  ##n, n alpha type rection##
  ###########################
  elif collision["interactionType"] < 100:
    #If have our banked ion, the energy needs to be corrected
    if len(missingLightIonPDGCodes)==0 and len(bankedLightIons)==1:
      Q = calcQValue(incidentPDGCodes,expectedBankPDGCodes)
      bankedLightIons[0].KE -= Q
      if bankedLightIons[0].KE < 0:
        bankedLightIons[0].KE += Q
      outgoingParticles.append(bankedLightIons[0])

    #Otherwise our ion is missing. The system is underconstrained--we don't know how much energy goes to the residual nucleus vs. alpha
    #Assume kinematic mass splittiing
    elif len(missingLightIonPDGCodes)==1 and len(bankedLightIons)==0:
      Q = calcQValue(incidentPDGCodes, expectedBankPDGCodes)

      # Get incoming/outgoing neutron energies
      parentNeutronCandidates = [i for i in incidentParticles if i.pdgCode==2112]
      if not len(parentNeutronCandidates)==1:
        print(f"ERROR! Should only have one parent neutron in light ion calculation for history {historyNum}")
        sys.exit()
      parentNeutron = parentNeutronCandidates[0]

      exitingNeutronCandidates = [i for i in outgoingParticles if i.pdgCode==2112]
      if not len(exitingNeutronCandidates)==1:
        print(f"ERROR! Should only have one exiting neutron in light ion calculation for history {historyNum}")
        sys.exit()
      exitingNeutron = exitingNeutronCandidates[0]

      E_available = parentNeutron.KE - exitingNeutron.KE - Q 
      if E_available<0:
        print(f"Kinematic error in light ion calculation! E_available < 0 for history {historyNum}")
        sys.exit()

      # Masses in MeV
      ejectileZ,ejectileA = particleUtils.calcZA_from_pdgCode(missingLightIonPDGCodes[0])
      residualZ,residualA = particleUtils.calcZA_from_pdgCode(residualNucleusPDGCode)
      malpha = particleUtils.massTable[ejectileZ][ejectileA]
      mresid = particleUtils.massTable[residualZ][residualA]

      # Split energy by mass ratio (non-relativistic)
      KE_alpha = E_available * mresid / (malpha + mresid)

      # Random isotropic direction
      u, v, w = np.random.normal(size=3)
      norm = np.linalg.norm([u, v, w])
      u, v, w = u / norm, v / norm, w / norm

      outgoingParticles.append(
        particleUtils.Particle(ejectileZ, ejectileA,
                               collision["x"], collision["y"], collision["z"],
                               u, v, w,
                               KE_alpha, collision["time"])
      )
      #(residual nucleus still handled externally

    else:
      print(f"New edge case discovered for history {historyNum}") #Will add better print if we see this
      sys.exit()

  ##########################
  ##n,alpha type reactions##
  ##########################
  #For these, usually we see some other banked particle, but we're assuming the ptrac MT code is truth
  elif collision["interactionType"] > 100:
    #Particle is correct, but the KE needs to be adjusted
    if len(bankedLightIons)==1 and len(missingLightIonPDGCodes)==0:
      Q = calcQValue(incidentPDGCodes,expectedBankPDGCodes)
      bankedLightIons[0].KE -= Q 
      if bankedLightIons[0].KE < 0:
        bankedLightIons[0].KE += Q 
      outgoingParticles.append(bankedLightIons[0])

    elif len(missingLightIonPDGCodes) == 1:
      # Q-value
      Q = calcQValue(incidentPDGCodes, expectedBankPDGCodes)

      # Determine PDG and (Z,A) of missing light ion and residual
      ejectilePDG = missingLightIonPDGCodes[0]
      ejectileZ, ejectileA = particleUtils.calcZA_from_pdgCode(ejectilePDG)
      residualZ, residualA = particleUtils.calcZA_from_pdgCode(residualNucleusPDGCode)

      # Target and incident neutron
      targetParticleCandidates = [i for i in incidentParticles if i.pdgCode>1000020040]
      if not len(targetParticleCandidates)==1:
        print(f"ERROR! Should only have target in light ion calculation for history {historyNum}")
        sys.exit()
      targetParticle = targetParticleCandidates[0]

      parentNeutronCandidates = [i for i in incidentParticles if i.pdgCode==2112]
      if not len(parentNeutronCandidates)==1:
        print(f"ERROR! Should only have one parent neutron in light ion calculation for history {historyNum}")
        sys.exit()
      parentNeutron = parentNeutronCandidates[0]

      #Kinematics
      m_a = parentNeutron.mass
      m_A = targetParticle.mass
      m_b = particleUtils.massTable[ejectileZ][ejectileA]
      m_B = particleUtils.massTable[residualZ][residualA]
      M = m_a + m_A

      KE_in = parentNeutron.KE  # incoming neutron energy
      uhat = np.array([parentNeutron.u,parentNeutron.v,parentNeutron.w])

      # CM velocity magnitude
      v_cm_mag = (m_a / M) * math.sqrt(2 * KE_in / m_a)

      # Ejectile velocity magnitude in CM frame
      v_b_star = math.sqrt((2 * m_B * m_b * (KE_in + Q)) / (M**2 * (m_B + m_b)))

      # Random isotropic direction in CM frame
      theta = np.arccos(2 * np.random.uniform() - 1)
      phi = 2 * np.pi * np.random.uniform()
      u_cm = np.array([
        np.sin(theta) * np.cos(phi),
        np.sin(theta) * np.sin(phi),
        np.cos(theta)
      ])

      # Boost to lab frame
      v_lab_vec = v_cm_mag * uhat + v_b_star * u_cm
      v_lab_mag = np.linalg.norm(v_lab_vec)
      if v_lab_mag == 0:
        print(f"Zero velocity magnitude in lab frame during boost for history {historyNum}")
        sys.exit()

      # Normalize direction
      u, v, w = v_lab_vec / v_lab_mag
      KE = 0.5 * m_b * v_lab_mag**2  # in MeV

      # Append to outgoingParticles
      outgoingParticles.append(
        particleUtils.Particle(ejectileZ, ejectileA,
                               collision["x"], collision["y"], collision["z"],
                               u, v, w,
                               KE, collision["time"])
      )
    elif len(missingLightIonPDGCodes)==2:
      Q = calcQValue(incidentPDGCodes, expectedBankPDGCodes)

      # Get incident neutron
      parentNeutronCandidates = [i for i in incidentParticles if i.pdgCode == 2112]
      if len(parentNeutronCandidates) != 1:
          print("ERROR! Should only have one parent neutron in light ion calculation")
          sys.exit()
      parentNeutron = parentNeutronCandidates[0]

      # Get target
      targetCandidates = [i for i in incidentParticles if i.pdgCode > 1000020040]
      if len(targetCandidates) != 1:
        print("ERROR! Should only have one target nucleus in light ion calculation")
        sys.exit()
      target = targetCandidates[0]

      # Particle identities
      pdg_b = missingLightIonPDGCodes[0]
      pdg_B = missingLightIonPDGCodes[1]
      Zb, Ab = particleUtils.calcZA_from_pdgCode(pdg_b)
      ZB, AB = particleUtils.calcZA_from_pdgCode(pdg_B)
      m_b = particleUtils.massTable[Zb][Ab]
      m_B = particleUtils.massTable[ZB][AB]

      m_a = parentNeutron.mass
      m_A = target.mass
      KE_in = parentNeutron.KE
      M = m_a + m_A

      uhat = np.array([parentNeutron.u, parentNeutron.v, parentNeutron.w])
      v_cm_mag = (m_a / M) * math.sqrt(2 * KE_in / m_a)
      v_cm_vec = v_cm_mag * uhat
      arg = (2 * m_B * m_b * (KE_in + Q)) / (M**2 * (m_B + m_b))
      if arg < 0:
        arg=0
      # CM velocity magnitude of ejectile b (and B gets recoil)
      v_b_star = math.sqrt(arg)

      # Isotropic CM emission
      theta = np.arccos(2 * np.random.uniform() - 1)
      phi = 2 * np.pi * np.random.uniform()
      u_cm = np.array([
          np.sin(theta) * np.cos(phi),
          np.sin(theta) * np.sin(phi),
          np.cos(theta)
      ])

      # Particle b
      v_lab_b = v_cm_vec + v_b_star * u_cm
      v_lab_mag_b = np.linalg.norm(v_lab_b)
      KE_b = 0.5 * m_b * v_lab_mag_b**2
      u_b, v_b, w_b = v_lab_b / v_lab_mag_b

      outgoingParticles.append(particleUtils.Particle(Zb, Ab,
          collision["x"], collision["y"], collision["z"],
          u_b, v_b, w_b, KE_b, collision["time"]))

      # Particle B (opposite in CM)
      v_lab_B = v_cm_vec - v_b_star * u_cm
      v_lab_mag_B = np.linalg.norm(v_lab_B)
      KE_B = 0.5 * m_B * v_lab_mag_B**2
      u_B, v_B, w_B = v_lab_B / v_lab_mag_B

      outgoingParticles.append(particleUtils.Particle(ZB, AB,
          collision["x"], collision["y"], collision["z"],
          u_B, v_B, w_B, KE_B, collision["time"]))

  
    else:
      print(f"Unhandled case: MT>100 but light ion not uniquely missing for history {historyNum}")
      sys.exit()
  
  else:
    print(f"New edge case discovered in light processing code for history {historyNum}")
    sys.exit()

  return outgoingParticles
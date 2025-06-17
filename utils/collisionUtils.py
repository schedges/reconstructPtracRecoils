import numpy as np
from utils import particleUtils,gammaCascadeHelper,lightIonProcessing
import sys
import math

excitationEnergy_tol = 1e-2 
levelsDict = {}

def getExpectedOutgoingParticlePDGCodes(interactionType,targetPDGCode):
  targetZ,targetA = particleUtils.calcZA_from_pdgCode(targetPDGCode)
  if interactionType==2:
    residual_Z = targetZ
    residual_A = targetA
    return [2112,particleUtils.calcPDGCode_from_ZA(residual_Z,residual_A)]
  elif interactionType==22:
    residual_Z = targetZ-2
    residual_A = targetA-4
    return [2112,1000020040,particleUtils.calcPDGCode_from_ZA(residual_Z,residual_A)]
  elif interactionType==28:
    residual_Z = targetZ-1
    residual_A = targetA-1
    return [2112,2212,particleUtils.calcPDGCode_from_ZA(residual_Z,residual_A)]
  elif interactionType==32:
    residual_Z = targetZ-1
    residual_A = targetA-2
    return [2112,1000010020,particleUtils.calcPDGCode_from_ZA(residual_Z,residual_A)]
  elif interactionType==33:
    residual_Z = targetZ-1
    residual_A = targetA-3
    return [2112,1000010030,particleUtils.calcPDGCode_from_ZA(residual_Z,residual_A)]
  elif interactionType>50 and interactionType<=91:
    residual_Z = targetZ
    residual_A = targetA
    return [2112,22,particleUtils.calcPDGCode_from_ZA(residual_Z,residual_A)]
  #102
  elif interactionType==102:
    residual_Z = targetZ
    residual_A = targetA+1
    return [22,particleUtils.calcPDGCode_from_ZA(residual_Z,residual_A)]
  elif interactionType==103:
    residual_Z = targetZ-1
    residual_A = targetA
    return [2212,particleUtils.calcPDGCode_from_ZA(residual_Z,residual_A)]
  elif interactionType==104:
    residual_Z = targetZ-1
    residual_A = targetA-1
    return [1000010020,particleUtils.calcPDGCode_from_ZA(residual_Z,residual_A)]
  elif interactionType==105:
    residual_Z = targetZ-1
    residual_A = targetA-2
    return [1000010030,particleUtils.calcPDGCode_from_ZA(residual_Z,residual_A)]
  elif interactionType==106:
    residual_Z = targetZ-2
    residual_A = targetA-2
    return [1000020030,particleUtils.calcPDGCode_from_ZA(residual_Z,residual_A)]
  elif interactionType==107:
    residual_Z = targetZ-2
    residual_A = targetA-3
    return [1000020040,particleUtils.calcPDGCode_from_ZA(residual_Z,residual_A)]
  else:
    print(f"Undefined MT code {interactionType} not defined in collisionUtils.getExpectedOutgoingParticleDefs() for history, skipping")
    return []
  


#################
##MAIN FUNCTION##
#################
def processNeutronCollision(collision,potentialParents_df,bankedParticles_df,historyNum):

  incidentParticles = []
  outgoingParticles = []

  ######################
  ##INCIDENT PARTICLES##
  ######################
  #Determine the parent--whichever line in potentialParents closest in time to collision (potential parents already preceeding events)
  delta_times = collision["time"] - potentialParents_df["time"]
  closest_idx = delta_times.idxmin()
  parent = potentialParents_df.loc[closest_idx]

  #Make incident neutron particle
  incidentParticles.append(particleUtils.Particle(0,1,
                                                  parent["x"],parent["y"],parent["z"],
                                                  parent["u"],parent["v"],parent["w"],
                                                  parent["energy"],parent["time"]))

  #Get target from collision line
  targetPolimiCode = collision["targetNucleus"]
  targetZ,targetA = particleUtils.calcZA_from_polimiCode(targetPolimiCode)
  targetPDGCode = particleUtils.calcPDGCode_from_ZA(targetZ,targetA)
  incidentParticles.append(particleUtils.Particle(targetZ,targetA,
                                                  collision["x"],collision["y"],collision["z"],
                                                  0,0,0,
                                                  0,collision["time"]))
  incidentPDGCodes = [2112,targetPDGCode]

  ######################
  ##EXITING PARTIICLES##
  ######################
  #Determine expected banked particles
  expectedBankPDGCodes = getExpectedOutgoingParticlePDGCodes(collision["interactionType"],targetPDGCode)
  bankedPDGCodes = []
  for bankedRow in bankedParticles_df.itertuples(index=False):
    ipt = bankedRow.particleType
    bankedPDGCodes.append(particleUtils.calcPDGCode_from_IPT(ipt))

  bankedLightIons=[]
  missingLightIonPDGCodes=[]

  #See if there's a residual nucleus
  residualNucleusPDGCode, residualZ, residualA = None, None, None
  for expectedParticlePDGCode in expectedBankPDGCodes:
    if expectedParticlePDGCode>1000020040: #alpha
      residualNucleusPDGCode = expectedParticlePDGCode
      residualZ,residualA = particleUtils.calcZA_from_pdgCode(residualNucleusPDGCode)

  ######################
  ##NEUTRON PROCESSING##
  ######################
  #Step through expected particles, look for the outgoing neutron
  for expectedParticlePDGCode in expectedBankPDGCodes:
    #If neutron--assign outgoing energy from collision line.
    if expectedParticlePDGCode==2112:
      outgoingParticles.append(particleUtils.Particle(0,1,
                                                  collision["x"],collision["y"],collision["z"],
                                                  collision["u"],collision["v"],collision["w"],
                                                  collision["energy"],collision["time"]))
      
  ####################
  ##GAMMA PROCESSING##
  ####################
  #Step through the expected PDG codes, looking for gammas
  for expectedParticlePDGCode in expectedBankPDGCodes:
    #If gamma--check this is correct for the supplied MT code, if not calculate our own cascade
    if expectedParticlePDGCode==22:
      #We can't have de-excitation gammas if we don't have a residual nucleus
      if not residualNucleusPDGCode is None:

        #TODO this messes up for radiative capture, we probably just need special handling for 
        #MT=102 but for now this defaults to throwing the highest level gamma possible
        expectedGammaLevel = int(collision["interactionType"]-50)
        
        #Load levels info
        ZZAAA = "{0:02d}{1:03d}".format(residualZ,residualA) 
        if ZZAAA not in levelsDict:
          levelsDict[ZZAAA] = gammaCascadeHelper.loadLevelsData(residualZ,residualA)
        
        #If specified level above max, use max
        maxLev = gammaCascadeHelper.findHighestLevelIndexWithDecayScheme(levelsDict[ZZAAA])
        if expectedGammaLevel > maxLev:
          expectedGammaLevel = maxLev #TODO bteter handling?
       
        #Check for energy conservation E_n_in > E_n_out + gamma
        outgoingNeutronCandidates = [i for i in outgoingParticles if i.pdgCode == 2112]
        if len(outgoingNeutronCandidates) != 1:
          print(f"Multi neutron emission found for hist {historyNum}, exiting")
          return []
        
        E_n_out = outgoingNeutronCandidates[0].KE
        incidentNeutronCandidates = [i for i in incidentParticles if i.pdgCode == 2112]
        if len(incidentNeutronCandidates) != 1:
          print(f"Multi neutron incident for hist {historyNum}, exiting")
          sys.exit()
        E_n_in = incidentNeutronCandidates[0].KE
        E_remaining = E_n_in - E_n_out
        if E_remaining < 0:
          print("Error! E_n_out > E_n_in!")
          sys.exit()
        if levelsDict[ZZAAA][expectedGammaLevel].energy > E_remaining:
          expectedGammaLevel = gammaCascadeHelper.findHighestLevelIndexBeforeEnergy(levelsDict[ZZAAA],E_remaining)

        #Get expected gamma energy
        expectedGammaEnergy = levelsDict[ZZAAA][expectedGammaLevel].energy

        #Get total banked gamma energy, see if it matches the expected level. 
        totalBankedGammaEnergy = bankedParticles_df.loc[bankedParticles_df["particleType"]==2, "energy"].sum()
        #print(f"Expected gamma energy {expectedGammaEnergy}, gamma energy in ptrac: {totalBankedGammaEnergy}")
        
        #If not, re-calculate
        if abs(totalBankedGammaEnergy-expectedGammaEnergy) > excitationEnergy_tol:
          #Returns gamma objects with energy, u, v, w
          bankedGammaLines = gammaCascadeHelper.generateBankedGammas(expectedGammaLevel,levelsDict[ZZAAA])

        #If within tolerance, just use banked gammas
        else:
          gamma_df = bankedParticles_df[bankedParticles_df["particleType"] == 2]
          bankedGammaLines = [gammaCascadeHelper.Gamma(row.energy, row.u, row.v, row.w) for row in gamma_df.itertuples(index=False)]

        #Write gammas
        for g in bankedGammaLines:
          outgoingParticles.append(particleUtils.Particle(0,0,
                                                      collision["x"],collision["y"],collision["z"],
                                                      g.u,g.v,g.w,
                                                      g.energy,collision["time"]))
      else:
        print("Error! No residual nucleus defined, cannot compute de-excitation gammas")
        print(f"HistoryNum = {historyNum}")

  ########################
  ##LIGHT ION PROCESSING##
  ########################
  #Step through the expected PDG codes, looking for light ions
  for expectedParticlePDGCode in expectedBankPDGCodes:
    #And this is a process that should produce a light ion
    if collision["interactionType"] in [22, 28, 32, 33, 103, 104, 105, 106, 107]:
      #This is the light ion
      if expectedParticlePDGCode>2112 and expectedParticlePDGCode<=1000020040:
        #Check if this is banked
        if expectedParticlePDGCode in bankedPDGCodes:
          lightIonZ,lightIonA = particleUtils.calcZA_from_pdgCode(expectedParticlePDGCode)
          ipt = particleUtils.calcPolimiIPT_from_pdgCode(expectedParticlePDGCode)
          row = bankedParticles_df[bankedParticles_df["particleType"] == ipt].iloc[0]  #TODO: Issue if 2 light ions of same species present
          bankedLightIons.append(particleUtils.Particle(
              lightIonZ, lightIonA,
              collision["x"], collision["y"], collision["z"],
              row["u"], row["v"], row["w"],
              row["energy"], collision["time"]))
        else:
          missingLightIonPDGCodes.append(expectedParticlePDGCode)

    #This should just be the residual nucleus
    elif expectedParticlePDGCode != 2112 and expectedParticlePDGCode != 22 and expectedParticlePDGCode != residualNucleusPDGCode:
      print(f"Error! Found banked particle {expectedParticlePDGCode} and are not processing it")
      print(f"HistoryNum = {historyNum}")

  ###########################
  #Light ion post-processing#
  ###########################
  if collision["interactionType"] in [22, 28, 32, 33, 103, 104, 105, 106, 107]:
    outgoingParticles = lightIonProcessing.handleLightIons(bankedLightIons,missingLightIonPDGCodes,
                                                          incidentPDGCodes,expectedBankPDGCodes, #For calculating Q values
                                                          incidentParticles,collision,residualNucleusPDGCode, #For kinematic calculations
                                                          outgoingParticles, #Updated in output
                                                          historyNum) #for debugging

  ############################
  ##Calculate Nuclear Recoil##
  ############################
  #If we have a residual nucleus, we need to calculate its recoil and add it to outgoing Particles
  if not residualNucleusPDGCode is None:
    for part in incidentParticles:
      part.create4Vec()
    incident_4vecs = [part.fourVec for part in incidentParticles]
    incoming_total = np.sum(incident_4vecs, axis=0)

    for part in outgoingParticles:
      part.create4Vec()
    outgoing_4vecs = [part.fourVec for part in outgoingParticles]
    outgoing_total = np.sum(outgoing_4vecs, axis=0)

    recoil_4vec = incoming_total - outgoing_total

    #Calc u,v,w,KE
    mass = particleUtils.massTable[residualZ][residualA]
    px, py, pz = recoil_4vec[1:]
    p_mag = np.linalg.norm([px, py, pz])
    u, v, w = px/p_mag, py/p_mag, pz/p_mag
    KE = recoil_4vec[0] - mass

    '''
    print(f"Ion KE is {KE}")
    print(historyNum)
    print(incident_4vecs)
    print(outgoing_4vecs)
    print()
    '''

    outgoingParticles.append(particleUtils.Particle(residualZ,residualA,
                                                    collision["x"],collision["y"],collision["z"],
                                                    u,v,w,
                                                    KE,collision["time"]))

  #Return ions produced from this collision - We don't return gammas, neutrons #TODO
  ions = [i for i in outgoingParticles if not i.pdgCode == 22 and not i.pdgCode==2112]

  for ion in ions:
    if ion.KE < -0.000001:
      return []
    
  return ions


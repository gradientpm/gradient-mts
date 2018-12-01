import os, sys, optparse
parser = optparse.OptionParser()
parser.add_option('-i','--input', help='input XML file')
parser.add_option('-p','--python', help='python file for generate names')
parser.add_option('-n', '--name', help='prefix name XML', default="scene")
parser.add_option('-o','--output', help='output file directory')
(opts, args) = parser.parse_args()

# === Get the name of the scene
SCENE_NAME = opts.name

#########################
# Behavior configuration
#########################
# --- No reducing the radius
FIXRADII = True
# --- Ratio of sample multiplication for SPPM
MULTSPPM = 1 #< Scene dependent (equal number)
# --- Enable generation two stage rendering
ACTIVE_ALL_SPLITTING = True
# --- Enable xml generation for PT based techniques
SIMPLETECHNIQUES = True
OTHERMLTTECHNIQUES = False
PATHMLTTECHNIQUES = True
TWOSTAGETECHNIQUES = True
MULTISTAGETECHNIQUES = True
# --- Generate xml for kappa reference techniques
#TODO: Not used now
GENERATEREFKAPPA = False

# --- Use BPM variation
USE_BPM_VARIATION = False
USE_VCM_VARIATION = False
# --- 
COMPUTE_REFERENCE_CONFIG = False
SPPM_AND_OTHER = True

#######
# Global config
#######
# --- Other options
bounceRoughness = "0.41"
# Use Normal to check that the MC is stuck, default is different
defaultStrategy = "Different" 
epsilon="0.0001"
# --- Number of photon per passes
photonCount = "10000000"
photonCountSPPM = str(int(photonCount)*MULTSPPM)
# ---  Number of iteration and timeout
maxPass = "20000"
timeMax = "3700" # 1h + 100 sec

# --- Path structure 
rrDepthPhoton = "1"
rrDepthPaths = "5"
maxDepth = "-1"
# --- Radius informations
initialScale = "1.0"
alpha = "0.9"
if(FIXRADII):
    initialScale = "1.0"
    alpha = "1.0"
# --- Other configuratios
samplesPT = "32"
samplesBDPT = "16"
samplesERPT = "8"
samplesMLT = "8000" # Samples of PSSMLT & Veach MLT
DumpContribution = "false"
directSamples = "-1"
dumpTimeImage = "30"
stepSnapshot = "1"
if(COMPUTE_REFERENCE_CONFIG):
    referenceMod = "true"
    photonCount = str(int(photonCount)*4)
    stepSnapshot = "10" 
    maxPass = "20000"
    

def xmlEntry(type, name, value):
    return (type,{"name":name,"value":value})

def xmlIntAttrs(type, attrs):
    return ("integrator",{"type":type, "attrs":attrs})

class VariationGenerator:
    def __init__(self):
        self.variantes = []
    
    def addVariation(self, v):
        self.variantes.append(v)
    
    def applyVariation(self, oriInt):
        """Returns all version of the integrator
        with the differents variations"""
        integrators = [oriInt]
        for i in range(len(self.variantes)):
            currVar = self.variantes[i]
            newVersionInt = []
            # There is two cases:
            #  - If there is only one options, no variation is made (no dict struct)
            #  - If there is more than one option, generate severals version
            if(len(currVar) == 1):
                for intEl in integrators:
                    intEl["attrs"].extend(currVar[0])
                    newVersionInt.append(intEl)
            else:
                for intEl in integrators:
                    for kV, vV in currVar.iteritems():
                        typeInt = intEl["type"]
                        nameInt = intEl["name"]
                        attrsInt = intEl["attrs"]
                
                        newIntegrator = {"type" : typeInt,
                               "name" : str(nameInt+"_"+kV),
                               "attrs" : attrsInt[:] + vV}
                        print("INFO: Create variation: "+str(newIntegrator["name"]))
                        newVersionInt.append(newIntegrator)
            
            #### End new generation
            integrators = newVersionInt
        return integrators
    
#################
# Create variators
#################
localImpVar = VariationGenerator()
otherVar = VariationGenerator()
sppmVar = VariationGenerator()


if(ACTIVE_ALL_SPLITTING):
    localImpVar.addVariation([[xmlEntry("string", "splittingStrat", "variance")]])

if(USE_BPM_VARIATION or USE_VCM_VARIATION):
    algoVar = {"sppm":[xmlEntry("string", "technique", "sppm")]}
    if(USE_BPM_VARIATION):
        algoVar["bpm"] = [xmlEntry("string", "technique", "bpm")]
    if(USE_VCM_VARIATION):
        algoVar["vcm"] = [xmlEntry("string", "technique", "vcm")]
        
    localImpVar.addVariation(algoVar)
    otherVar.addVariation(algoVar)
    sppmVar.addVariation(algoVar)
else:
    # Per default we standard one SPPM
    localImpVar.addVariation([[xmlEntry("string", "technique", "sppm")]])
    otherVar.addVariation([[xmlEntry("string", "technique", "sppm")]])
    sppmVar.addVariation([[xmlEntry("string", "technique", "sppm")]])


#=============
# For the 3 chain
# We need to have differente Variation
otherVar.addVariation([[xmlEntry("integer", "numberChains", "1"), 
                        xmlEntry("boolean", "useMIS", "false")]])
localImpVar.addVariation([[xmlEntry("integer", "numberChains", "3")]])
localImpVar.addVariation([[xmlEntry("boolean", "useMIS", "true")]])
    
INTEGRATORS = []

##############################################
##############################################
##############################################
############# Local Imp ######################
##############################################
##############################################
##############################################
MSPPMDefault = [xmlEntry("integer", "maxDepth", maxDepth),
              xmlEntry("integer", "rrDepth", rrDepthPhoton),
              xmlEntry("float", "alpha", alpha),
              xmlEntry("integer", "photonCount", photonCount),
              xmlEntry("integer", "stepSnapshot", stepSnapshot),
              xmlEntry("integer", "maxPass", maxPass),
              xmlEntry("float", "initialScale", initialScale),
              xmlEntry("boolean", "dumpImagePass", DumpContribution),
              xmlEntry("float", "bounceRoughness", bounceRoughness),
              xmlEntry("string", "numberStrategy", defaultStrategy)]

########
# Local Imp function 3D
# Metropolis sampling
########
LocalImpAttrs3D = MSPPMDefault[:] + [xmlEntry("string", "impFunc", "LocalImp"),
                                          xmlEntry("boolean", "use3D", "true"),
                                          xmlEntry("boolean", "clearFrist", "true")]

LocalImpIntegrator3D = {"type" : "msppm",
                   "name" : "LocalImp3D",
                   "attrs" : LocalImpAttrs3D}

# === Add classical integrator
INTEGRATORS.extend(localImpVar.applyVariation(LocalImpIntegrator3D))

########
# Local Imp function 2D
# Progressive estimation
########
LocalImpAttrs2D = MSPPMDefault[:] + [xmlEntry("string", "impFunc", "LocalImp"),
                                     xmlEntry("boolean", "use3D", "false"),
                                     xmlEntry("string", "strategy", "NumberPhoton"),
                                     xmlEntry("boolean", "multiStage", "true"),
                                     xmlEntry("float", "maxRange", "100000")]
LocalImpIntegrator2D = {"type" : "msppm",
                             "name" : "LocalImp2D",
                             "attrs" : LocalImpAttrs2D[:]}
INTEGRATORS.extend(localImpVar.applyVariation(LocalImpIntegrator2D))

########
# Inv Statistic
# Using the estimate Kappa statistic
########
if(GENERATEREFKAPPA):
    InvDensitySPPMAttrs = MSPPMDefault[:] + [xmlEntry("string", "impFunc", "InvDensity"),
                                         xmlEntry("string", "densityImg", "data/density.hdr")]
    
    InvDensitySPPMIntegrator = {"type" : "msppm",
                            "name" : "InvDensity",
                            "attrs" : InvDensitySPPMAttrs}
    
    INTEGRATORS.extend(addMiVariationIfNecessary(InvDensitySPPMIntegrator))
    
    InvDensityEpsSPPMIntegrator = {"type" : "msppm",
                            "name" : "InvDensity_eps",
                            "attrs" : InvDensitySPPMAttrs[:] + [xmlEntry("float", "epsilon", "0.0001")]}
    
    INTEGRATORS.extend(addMiVariationIfNecessary(InvDensityEpsSPPMIntegrator))

##############################################
##############################################
##############################################
############# Photon techniques###############
############# For comparison   ###############
##############################################
##############################################

########
# VSPPM
# Visiblity driven sampling
########
VSPPMAttrs = MSPPMDefault[:] + [xmlEntry("string", "impFunc", "VSPPM"),
                                xmlEntry("boolean", "removeFirstIteration", "false")]

VSPPMIntegrator = {"type" : "msppm",
                   "name" : "VSPPM",
                   "attrs" : VSPPMAttrs}
INTEGRATORS.extend(otherVar.applyVariation(VSPPMIntegrator))

########
# InverseSurface
# INverse surface driven sampling
########
InvSurfAttrs = MSPPMDefault[:] + [xmlEntry("string", "impFunc", "InvSurf"),
                                  xmlEntry("boolean", "removeFirstIteration", "false")]
InvSurfIntegrator = {"type" : "msppm",
                   "name" : "InvSurf",
                   "attrs" : InvSurfAttrs}
INTEGRATORS.extend(otherVar.applyVariation(InvSurfIntegrator))

########
# ISPPM
# Inverse density 
########
ISPPMAttrs = MSPPMDefault[:] + [xmlEntry("string", "impFunc", "ISPPM"),
                                xmlEntry("boolean", "removeFirstIteration", "false")]
ISPPMIntegrator = {"type" : "msppm",
                   "name" : "ISPPM",
                   "attrs" : ISPPMAttrs}
INTEGRATORS.extend(otherVar.applyVariation(ISPPMIntegrator))

#######
# Visual PPM
# The new paper using the importon
#######

VisualPPMAttrs = MSPPMDefault[:] + [xmlEntry("string", "impFunc", "Visual"),
                                xmlEntry("boolean", "removeFirstIteration", "false")]

VisualPPMIntegrator = {"type" : "msppm",
                       "name" : "Visual",
                       "attrs" : VisualPPMAttrs}
INTEGRATORS.extend(otherVar.applyVariation(VisualPPMIntegrator))

########
# SPPM
# Reference technique without
# any metropolis sampling
########
if(SPPM_AND_OTHER):
    SPPMAttrs = [xmlEntry("integer", "maxDepth", maxDepth),
                xmlEntry("integer", "rrDepth", rrDepthPhoton),
                xmlEntry("float", "alpha", alpha),
                xmlEntry("integer", "photonCount", photonCountSPPM),
                xmlEntry("integer", "stepSnapshot", stepSnapshot),
                xmlEntry("integer", "maxPass", maxPass),
                xmlEntry("float", "initialScale", initialScale),
                xmlEntry("boolean", "dumpImagePass", DumpContribution),
                xmlEntry("float", "bounceRoughness", bounceRoughness)]
    SPPMIntegrator = {"type" : "sppm_splat",
                      "name" : "SPPM",
                      "attrs" : SPPMAttrs}
    INTEGRATORS.extend(sppmVar.applyVariation(SPPMIntegrator))
    
    ########
    # PT bouncing
    ########
    PTBounceAttrs = [xmlEntry("integer", "maxDepth", maxDepth),
               xmlEntry("integer", "rrDepth", rrDepthPaths),
               xmlEntry("float", "bounceRoughness", bounceRoughness)]
    PTBounceAvgAttrs = [xmlEntry("integer", "maxPass", maxPass),
                  xmlEntry("integer", "maxRenderingTime", timeMax),
                  xmlIntAttrs("path_bounce", PTBounceAttrs)]
    PTBounceIntegrator = {"type" : "avg",
                    "name" : "PTBounce",
                    "attrs" : PTBounceAvgAttrs,
                    "samples" : samplesPT} # Rewritting samples counts
    INTEGRATORS.append(PTBounceIntegrator)

if(SIMPLETECHNIQUES):
    ########
    # PT
    ########
    PTAttrs = [xmlEntry("integer", "maxDepth", maxDepth),
               xmlEntry("integer", "rrDepth", rrDepthPaths)]
    PTAvgAttrs = [xmlEntry("integer", "maxPass", maxPass),
                  xmlEntry("integer", "maxRenderingTime", timeMax),
                  xmlIntAttrs("path", PTAttrs)]
    PTIntegrator = {"type" : "avg",
                    "name" : "PT",
                    "attrs" : PTAvgAttrs,
                    "samples" : samplesPT} # Rewritting samples counts
    
    INTEGRATORS.append(PTIntegrator)
    
    #######
    # BDPT
    #######
    BDPTAttrs = [xmlEntry("integer", "maxDepth", maxDepth),
               xmlEntry("integer", "rrDepth", rrDepthPaths)]
    BDPTAvgAttrs = [xmlEntry("integer", "maxPass", maxPass),
                  xmlEntry("integer", "maxRenderingTime", timeMax),
                  xmlIntAttrs("bdpt", BDPTAttrs)]
    BDPTIntegrator = {"type" : "avg",
                    "name" : "BDPT",
                    "attrs" : BDPTAvgAttrs,
                    "samples" : samplesBDPT} # Rewritting samples counts
    
    INTEGRATORS.append(BDPTIntegrator)

if(OTHERMLTTECHNIQUES):
    #######
    # ERPT
    #######
    ERPTAttrs = [xmlEntry("integer", "maxDepth", maxDepth),
               xmlEntry("integer", "rrDepth", rrDepthPaths),
               xmlEntry("integer","directSamples", directSamples)]
    ERPTAvgAttrs = [xmlEntry("integer", "maxPass", maxPass),
                  xmlEntry("integer", "maxRenderingTime", timeMax),
                  xmlIntAttrs("erpt", ERPTAttrs)]
    ERPTIntegrator = {"type" : "avg",
                    "name" : "ERPT",
                    "attrs" : ERPTAvgAttrs,
                    "samples" : samplesERPT} # Rewritting samples counts
    
    INTEGRATORS.append(ERPTIntegrator)

    MLTAttrs = [xmlEntry("integer", "maxDepth", maxDepth),
               xmlEntry("integer", "rrDepth", rrDepthPaths),
               xmlEntry("integer","directSamples", directSamples),
               xmlEntry("integer","maxTimeImgDump",dumpTimeImage),
               xmlEntry("integer","timeout", timeMax),
               xmlEntry("string", "stagedTechnique", "NoStage")] # NEED TO HAVE STAGE AT THE END
    
    MLTIntegrator = {"type" : "mlt",
                    "name" : "MLT",
                    "attrs" : MLTAttrs,
                    "samples" : samplesMLT, # Rewritting samples counts
                    "samplerType" : "independent"}# Rewritting sampler type
    
    INTEGRATORS.append(MLTIntegrator)
    
    if(TWOSTAGETECHNIQUES):
        MLTTWOIntegrator = MLTIntegrator.copy()
        MLTTWOIntegrator["attrs"] = MLTIntegrator["attrs"][:]
        MLTTWOIntegrator["attrs"][-1] =  xmlEntry("string", "stagedTechnique", "TwoStage")
        MLTTWOIntegrator["name"] += "TwoStage"
        INTEGRATORS.append(MLTTWOIntegrator)
        
if(PATHMLTTECHNIQUES):
    #######
    # ERPT Manifold
    #######
    ERPTManiAttrs = [xmlEntry("integer", "maxDepth", maxDepth),
               xmlEntry("integer", "rrDepth", rrDepthPaths),
               xmlEntry("boolean", "manifoldPerturbation", "true"),
               xmlEntry("integer","directSamples", directSamples)]
    ERPTManiAvgAttrs = [xmlEntry("integer", "maxPass", maxPass),
                  xmlEntry("integer", "maxRenderingTime", timeMax),
                  xmlIntAttrs("erpt", ERPTManiAttrs)]
    ERPTManiIntegrator = {"type" : "avg",
                    "name" : "ERPTMani",
                    "attrs" : ERPTManiAvgAttrs,
                    "samples" : samplesERPT} # Rewritting samples counts
    
    INTEGRATORS.append(ERPTManiIntegrator)
    
    #######
    # PSSMLT
    #######
    PSSMLTAttrs = [xmlEntry("integer", "maxDepth", maxDepth),
               xmlEntry("integer", "rrDepth", rrDepthPaths),
               xmlEntry("integer","directSamples", directSamples),
               xmlEntry("integer","maxTimeImgDump",dumpTimeImage),
               xmlEntry("integer","timeout", timeMax),
               xmlEntry("boolean", "recomputeNormalisation", "true"),
               xmlEntry("string", "stagedTechnique", "NoStage")] # NEED TO HAVE STAGE AT THE END
    
    PSSMLTIntegrator = {"type" : "pssmlt",
                    "name" : "PSSMLT",
                    "attrs" : PSSMLTAttrs,
                    "samples" : samplesMLT, # Rewritting samples counts
                    "samplerType" : "independent"}# Rewritting sampler type
    
    INTEGRATORS.append(PSSMLTIntegrator)
    
    #######
    # MLT Manifold
    #######
    MLTManiAttrs = [xmlEntry("integer", "maxDepth", maxDepth),
               xmlEntry("integer", "rrDepth", rrDepthPaths),
               xmlEntry("integer","directSamples", directSamples),
               xmlEntry("integer","maxTimeImgDump",dumpTimeImage),
               xmlEntry("integer","timeout", timeMax),
               xmlEntry("boolean", "manifoldPerturbation", "true"),
               xmlEntry("string", "stagedTechnique", "NoStage")] # NEED TO HAVE STAGE AT THE END
    
    MLTManiIntegrator = {"type" : "mlt",
                    "name" : "MLTMani",
                    "attrs" : MLTManiAttrs,
                    "samples" : samplesMLT, # Rewritting samples counts
                    "samplerType" : "independent"}# Rewritting sampler type
    
    INTEGRATORS.append(MLTManiIntegrator)
    
    ########
    # PT bouncing
    ########
    PTBounceAttrsMLT = [xmlEntry("integer", "maxDepth", maxDepth),
               xmlEntry("integer", "rrDepth", rrDepthPaths),
               xmlEntry("integer","directSamples", directSamples),
               xmlEntry("integer","maxTimeImgDump",dumpTimeImage),
               xmlEntry("integer","timeout", timeMax),
               xmlEntry("boolean", "recomputeNormalisation", "true"),
               xmlEntry("string", "stagedTechnique", "NoStage"),
               xmlEntry("float", "bounceRoughness", bounceRoughness),
               xmlEntry("boolean", "bidirectional", "false"),
               xmlEntry("boolean", "missingComp", "true"),
               xmlEntry("boolean", "useKelemenMutation", "false"),
               xmlEntry("boolean", "useAMCMC", "true"),
               xmlEntry("float", "desiredSeed", "1.0")]
    PTBounceIntegratorMLT = {"type" : "pssmlt",
                    "name" : "PTBounceMLT",
                    "attrs" : PTBounceAttrsMLT,
                    "samples" : samplesMLT, # Rewritting samples counts
                    "samplerType" : "independent"}
    INTEGRATORS.append(PTBounceIntegratorMLT)
    
    ####
    # The two stage configuration ????
    ####
    if(TWOSTAGETECHNIQUES):
        MLTManiTWOIntegrator = MLTManiIntegrator.copy()
        MLTManiTWOIntegrator["attrs"] = MLTManiIntegrator["attrs"][:]
        MLTManiTWOIntegrator["attrs"][-1] =  xmlEntry("string", "stagedTechnique", "TwoStage")
        MLTManiTWOIntegrator["name"] += "TwoStage"
        INTEGRATORS.append(MLTManiTWOIntegrator)
        
        PSSMLTTWOIntegrator = PSSMLTIntegrator.copy()
        PSSMLTTWOIntegrator["attrs"] = PSSMLTIntegrator["attrs"][:]
        PSSMLTTWOIntegrator["attrs"][-1] =  xmlEntry("string", "stagedTechnique", "TwoStage")
        PSSMLTTWOIntegrator["name"] += "TwoStage"
        INTEGRATORS.append(PSSMLTTWOIntegrator)
    
    if(MULTISTAGETECHNIQUES):
        MLTManiMultiIntegrator = MLTManiIntegrator.copy()
        MLTManiMultiIntegrator["attrs"] = MLTManiIntegrator["attrs"][:]
        MLTManiMultiIntegrator["attrs"][-1] =  xmlEntry("string", "stagedTechnique", "MultiStage")
        MLTManiMultiIntegrator["name"] += "MultiStage"
        INTEGRATORS.append(MLTManiMultiIntegrator)
        
        PSSMLTMultiIntegrator = PSSMLTIntegrator.copy()
        PSSMLTMultiIntegrator["attrs"] = PSSMLTIntegrator["attrs"][:]
        PSSMLTMultiIntegrator["attrs"][-1] =  xmlEntry("string", "stagedTechnique", "MultiStage")
        PSSMLTMultiIntegrator["name"] += "MultiStage"
        INTEGRATORS.append(PSSMLTMultiIntegrator)

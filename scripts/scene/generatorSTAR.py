#######
# Global config
#######
# --- Path tracing techniques
rrDepthPaths = "12"
timeMax = "3700"  # 1h + 100 sec
samplesPT = "8"
samplesBDPT = "2"
samplesVCM = "1"

samplesGPT = "8"
samplesGBDPT = "2"
samplesGVCM = "1"
samplesMULT = 4

samplesPT = str(int(samplesPT) * samplesMULT)
samplesBDPT = str(int(samplesBDPT) * samplesMULT)
samplesVCM = str(int(samplesVCM) * samplesMULT)
samplesGPT = str(int(samplesGPT) * samplesMULT)
samplesGBDPT = str(int(samplesGBDPT) * samplesMULT)
samplesGVCM = str(int(samplesGVCM) * samplesMULT)

# --- SPPM options (basic)
# * Photons config
photonCountSPPMSplat = "10000000"
photonCountGPM = str(int(photonCountSPPMSplat) // 10)  # ensure gather is the same for both SPPM and GPM
maxDepth = "12"
rrDepthPhoton = "12"

# * Gather point config
shiftThreshold = "0.001"
shiftThresholdGPM = "0.05"
initialScale = "1.0"
alpha = "0.9"
initialRadius = "0.0"

# * Other config
maxPasses = "2000"
reconstructAlpha = "0.2"
forceBlackPixels = "true"
maxManifoldIterations = "5"

# Object where all integrators will be added
INTEGRATORS = []
separateReconstruction = True
directTracing = "false"

# Reconstruction
L1_RECONS = True
L2_RECONS = True
UNI_RECONS = True
WEIGHTED_RECONS = True
NAN_CHECK = "true"
FORCE_BLACK = "false"


##############################
# Functions
##############################
def xmlEntry(type, name, value):
    return (type, {"name": name, "value": value})


def xmlIntAttrs(type, attrs):
    return ("integrator", {"type": type, "attrs": attrs})


##############################
# SPPM splatting
##############################
SPPMAttrs = [xmlEntry("integer", "maxDepth", maxDepth),
             xmlEntry("integer", "rrDepth", rrDepthPhoton),
             xmlEntry("float", "alpha", alpha),
             xmlEntry("integer", "dumpIteration", "1"),
             xmlEntry("integer", "maxPasses", maxPasses),
             xmlEntry("float", "initialScale", initialScale),
             xmlEntry("float", "shiftThreshold", shiftThreshold),
             xmlEntry("float", "shiftThresholdGPM", shiftThresholdGPM),
             xmlEntry("boolean", "directTracing", directTracing),
             xmlEntry("integer", "volumePhotonCount", "0"),
             xmlEntry("boolean", "volumeRendering", "false"),
             xmlEntry("float", "initialRadius", initialRadius)]

SPPMIntegrator = {"type": "sppm",
                  "name": "SPPM",
                  "attrs": SPPMAttrs +
                           [xmlEntry("integer", "photonCount", photonCountGPM)]}
INTEGRATORS.append(SPPMIntegrator)

##############################
# BDPT
##############################
BDPTAttrs = [xmlEntry("integer", "maxDepth", maxDepth),
             xmlEntry("integer", "rrDepth", rrDepthPaths),
             xmlEntry("boolean", "directTracing", directTracing)]
BDPTAvgAttrs = [xmlEntry("integer", "maxPasses", maxPasses),
                xmlEntry("integer", "maxRenderingTime", timeMax),
                xmlEntry("integer", "dumpIteration", "1"),
                xmlIntAttrs("bdpt", BDPTAttrs)]
BDPTIntegrator = {"type": "avg",
                  "name": "BDPT",
                  "attrs": BDPTAvgAttrs,
                  "samples": samplesBDPT}  # Rewritting samples counts
INTEGRATORS.append(BDPTIntegrator)

##############################
# PT
##############################
PTAttrs = [xmlEntry("integer", "maxDepth", maxDepth),
           xmlEntry("integer", "rrDepth", rrDepthPaths),
           xmlEntry("boolean", "directTracing", directTracing)]
PTAvgAttrs = [xmlEntry("integer", "maxPasses", maxPasses),
              xmlEntry("integer", "maxRenderingTime", timeMax),
              xmlEntry("integer", "dumpIteration", "1"),
              xmlIntAttrs("path", PTAttrs)]
PTIntegrator = {"type": "avg",
                "name": "PT",
                "attrs": PTAvgAttrs,
                "samples": samplesPT}  # Rewritting samples counts

INTEGRATORS.append(PTIntegrator)

##############################
# BCD
##############################
BCDAttrs = PTAttrs[:] + [xmlEntry("integer", "maxPasses", maxPasses),
                         xmlEntry("integer", "maxRenderingTime", timeMax),
                         xmlEntry("integer", "dumpIteration", "1")]
BCDIntegrator = {"type": "bcd",
                 "name": "BCD",
                 "attrs": BCDAttrs,
                 "samples": samplesPT}  # Rewritting samples counts

INTEGRATORS.append(BCDIntegrator)

##############################
# NFOR
##############################
NFORAttrs = PTAttrs[:] + [xmlEntry("integer", "maxPasses", maxPasses),
                          xmlEntry("integer", "maxRenderingTime", timeMax),
                          xmlEntry("integer", "dumpIteration", "1")]

NFORIntegrator = {"type": "nfor",
                  "name": "NFOR",
                  "attrs": NFORAttrs,
                  "samples": samplesPT}  # Rewritting samples counts

INTEGRATORS.append(NFORIntegrator)

##############################
# VCM
##############################
VCMAttrs = [xmlEntry("integer", "maxDepth", maxDepth),
            xmlEntry("integer", "rrDepth", rrDepthPaths),
            xmlEntry("boolean", "directTracing", directTracing),
            xmlEntry("float", "initialRadius", initialRadius)]
VCMAvgAttrs = [xmlEntry("integer", "maxPasses", maxPasses),
               xmlEntry("integer", "maxRenderingTime", timeMax),
               xmlEntry("integer", "dumpIteration", "1"),
               xmlIntAttrs("vcm", VCMAttrs)]
VCMIntegrator = {"type": "avg",
                 "name": "VCM",
                 "attrs": VCMAvgAttrs,
                 "samples": samplesVCM}  # Rewritting samples counts
INTEGRATORS.append(VCMIntegrator)

##############################
# Gradient domain
##############################
GRADIENTAttrs = [xmlEntry("float", "shiftThreshold", shiftThreshold),
                 xmlEntry("float", "reconstructAlpha", reconstructAlpha),
                 xmlEntry("boolean", "forceBlackPixels", forceBlackPixels),
                 xmlEntry("integer", "maxManifoldIterations", maxManifoldIterations),
                 xmlEntry("boolean", "directTracing", directTracing)]

GRADIENTRecons = [("", [xmlEntry("boolean", "reconstructL1", "true"),
                        xmlEntry("boolean", "reconstructL2", "true"),
                        xmlEntry("boolean", "reconstructUni", "true"),
                        xmlEntry("boolean", "reconstructWeighted", "true")])]
if separateReconstruction:
    # GRADIENTRecons = []
    if L2_RECONS:
        GRADIENTRecons += [
            ("L2", [xmlEntry("boolean", "reconstructL2", "true"),
                    xmlEntry("boolean", "reconstructL1", "false"),
                    xmlEntry("boolean", "reconstructUni", "false"),
                    xmlEntry("boolean", "reconstructWeighted", "false")])]
    if L1_RECONS:
        GRADIENTRecons += [
            ("L1", [xmlEntry("boolean", "reconstructL2", "false"),
                    xmlEntry("boolean", "reconstructL1", "true"),
                    xmlEntry("boolean", "reconstructUni", "false"),
                    xmlEntry("boolean", "reconstructWeighted", "false")])]
    if UNI_RECONS:
        GRADIENTRecons += [
            ("Uni", [xmlEntry("boolean", "reconstructL2", "false"),
                     xmlEntry("boolean", "reconstructL1", "false"),
                     xmlEntry("boolean", "reconstructUni", "true"),
                     xmlEntry("boolean", "reconstructWeighted", "false")])]

    if WEIGHTED_RECONS:
        GRADIENTRecons += [
            ("Weighted", [xmlEntry("boolean", "reconstructL2", "false"),
                 xmlEntry("boolean", "reconstructL1", "false"),
                 xmlEntry("boolean", "reconstructUni", "false"),
                 xmlEntry("boolean", "reconstructWeighted", "true")])]

# NFOR XML generation
NFORRecons = []
NFORRecons += [("L2NFOR", [xmlEntry("boolean", "reconstructL2", "true"),
                           xmlEntry("boolean", "reconstructL1", "false"),
                           xmlEntry("boolean", "reconstructUni", "false"),
                           xmlEntry("boolean", "reconstructWeighted", "false"),
                           xmlEntry("boolean", "reconstructNfor", "true")])]
NFORRecons += [("L1NFOR", [xmlEntry("boolean", "reconstructL2", "false"),
                           xmlEntry("boolean", "reconstructL1", "true"),
                           xmlEntry("boolean", "reconstructUni", "false"),
                           xmlEntry("boolean", "reconstructWeighted", "false"),
                           xmlEntry("boolean", "reconstructNfor", "true")])]

for name, reconsAttrib in NFORRecons:
    if (name != ""):
        name = "_" + name
    GRADIENT_PT_Attrs = GRADIENTAttrs[:] + reconsAttrib[:] + \
                        [xmlEntry("integer", "maxDepth", maxDepth),
                         xmlEntry("integer", "rrDepth", rrDepthPaths),
                         xmlEntry("boolean", "forceBlackPixels", FORCE_BLACK),
                         xmlEntry("integer", "dumpIteration", "1")]
    GRADIENT_PT_Integrator = {"type": "gpt",
                              "name": "GPT" + name,
                              "attrs": GRADIENT_PT_Attrs,
                              "samples": samplesGPT,  # Rewritting samples counts
                              "filmType": "multifilm"}
    INTEGRATORS.append(GRADIENT_PT_Integrator)

for name, reconsAttrib in GRADIENTRecons:
    if (name != ""):
        name = "_" + name
    ####### PT
    GRADIENT_PT_Attrs = GRADIENTAttrs[:] + reconsAttrib[:] + \
                        [xmlEntry("integer", "maxDepth", maxDepth),
                         xmlEntry("integer", "rrDepth", rrDepthPaths),
                         xmlEntry("boolean", "forceBlackPixels", FORCE_BLACK),
                         xmlEntry("integer", "dumpIteration", "1")]

    # FIXME: Add max time and passes
    GRADIENT_PT_Integrator = {"type": "gpt",
                              "name": "GPT" + name,
                              "attrs": GRADIENT_PT_Attrs,
                              "samples": samplesGPT,  # Rewritting samples counts
                              "filmType": "multifilm"}
    INTEGRATORS.append(GRADIENT_PT_Integrator)

    ####### PT Path reuse
    GRADIENT_PT_Attrs = GRADIENTAttrs[:] + reconsAttrib[:] + \
                        [xmlEntry("integer", "maxDepth", maxDepth),
                         xmlEntry("integer", "rrDepth", rrDepthPaths),
                         xmlEntry("boolean", "forceBlackPixels", FORCE_BLACK),
                         xmlEntry("integer", "dumpIteration", "1"),
                         xmlEntry("boolean", "pathReuse", "true"),
                         xmlEntry("string", "shiftmapping", "explicit")]

    # FIXME: Add max time and passes
    GRADIENT_PT_Integrator = {"type": "gpt",
                              "name": "GPT" + name + "_PathReuse",
                              "attrs": GRADIENT_PT_Attrs,
                              "samples": samplesGPT,  # Rewritting samples counts
                              "filmType": "multifilm"}
    INTEGRATORS.append(GRADIENT_PT_Integrator)

    ####### BDPT
    GRADIENT_BDPT_Attrs = GRADIENTAttrs[:] + reconsAttrib[:] + \
                          [xmlEntry("integer", "maxDepth", maxDepth),
                           xmlEntry("integer", "rrDepth", rrDepthPaths),
                           xmlEntry("boolean", "forceBlackPixels", FORCE_BLACK),
                           xmlEntry("boolean", "lightImage", "true"),
                           xmlEntry("integer", "dumpIteration", "1")]

    # FIXME: Add max time and passes
    GRADIENT_BDPT_Integrator = {"type": "gbdpt",
                                "name": "GBDPT" + name,
                                "attrs": GRADIENT_BDPT_Attrs,
                                "samples": samplesGBDPT,  # Rewritting samples counts
                                "filmType": "multifilm"}
    INTEGRATORS.append(GRADIENT_BDPT_Integrator)

    ####### VCM
    GRADIENT_VCM_Attrs = GRADIENTAttrs[:] + reconsAttrib[:] + \
                         [xmlEntry("integer", "maxDepth", maxDepth),
                          xmlEntry("integer", "rrDepth", rrDepthPaths),
                          xmlEntry("boolean", "forceBlackPixels", FORCE_BLACK),
                          xmlEntry("boolean", "lightImage", "true"),
                          xmlEntry("integer", "dumpIteration", "1"),
                          xmlEntry("float", "initialRadius", initialRadius)]

    # FIXME: Add max time and passes
    GRADIENT_VCM_Integrator = {"type": "gvcm",
                               "name": "GVCM" + name,
                               "attrs": GRADIENT_VCM_Attrs,
                               "samples": samplesGVCM,  # Rewritting samples counts
                               "filmType": "multifilm"}
    INTEGRATORS.append(GRADIENT_VCM_Integrator)

    ##############################
    # GPM
    ##############################
    GPMAttrs = GRADIENTAttrs[:] + reconsAttrib[:] + \
               [xmlEntry("float", "shiftThresholdGPM", shiftThresholdGPM),
                xmlEntry("integer", "maxDepth", maxDepth),
                xmlEntry("integer", "rrDepth", rrDepthPhoton),
                xmlEntry("float", "alpha", alpha),
                xmlEntry("integer", "photonCount", photonCountGPM),
                xmlEntry("integer", "dumpIteration", "1"),
                xmlEntry("integer", "maxPasses", maxPasses),
                xmlEntry("float", "initialScale", initialScale),
                xmlEntry("boolean", "nanCheck", NAN_CHECK),
                xmlEntry("boolean", "forceBlackPixels", FORCE_BLACK),
                xmlEntry("float", "initialRadius", initialRadius)]

    GPMIntegrator = {"type": "gvpm",
                     "name": "GPM" + name,
                     "attrs": GPMAttrs}
    INTEGRATORS.append(GPMIntegrator)

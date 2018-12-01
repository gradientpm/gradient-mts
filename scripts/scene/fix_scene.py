import optparse
import xml.etree.ElementTree as ET

def rreplace(s, old, new, count = 1):
    return (s[::-1].replace(old[::-1], new[::-1], count))[::-1]

def BSDF(sceneRoot):
    BSDFDict = {}
    for bsdfNode in sceneRoot.iter("bsdf"):
        BSDFDict[bsdfNode.attrib["id"]] = bsdfNode

    i = 1
    for bsdfNode in sceneRoot.iter("bsdf"):
        if(bsdfNode.attrib["type"] == "mixturebsdf"):
            name = bsdfNode.attrib["id"]
            nameClean = name.replace("_None", "")
            name1 = rreplace(nameClean,"-material", "_I-material")
            name2 = rreplace(nameClean,"-material", "_II-material")

            child1 = BSDFDict[name1]
            child2 = BSDFDict[name2]

            diffuseReflectance = None

            # Get diffuse reflectance
            print(i, name, "\n *", child1.attrib["type"], "\n *", child2.attrib["type"], bsdfNode.find("string").attrib["value"])
            i+= 1
            if(child1.attrib["type"] != "diffuse"):
                raise "Non Diffuse"
            else:
                diffuseReflectance = child1.findall("rgb")
                if(not diffuseReflectance):
                    diffuseReflectance = child1.findall("ref")
                    if(not diffuseReflectance):
                        raise "no Diffuse RGB/Texture"
                diffuseReflectance = diffuseReflectance[0]

            if(child2.attrib["type"] != "roughconductor" and child2.attrib["type"] != "conductor" ):
                raise "Non Conductor"

def medLight(sceneRoot, nameSmoke):
    for emitterNode in sceneRoot.iter("emitter"):
        ref = ET.SubElement(emitterNode, 'ref')
        ref.attrib["id"] = nameSmoke

def medCamera(sceneRoot, nameSmoke):
    cameraNode = sceneRoot.find("sensor")
    refCamera = ET.SubElement(cameraNode, 'ref')
    refCamera.attrib["id"] = nameSmoke


if __name__ == "__main__":
    # --- Read all params
    parser = optparse.OptionParser()
    parser.add_option('-i','--input', help='input')
    parser.add_option('-B','--bsdf', help='fix BSDF', default=False, action="store_true")
    parser.add_option('-C','--camera', help='fix camera', default=False, action="store_true")
    parser.add_option('-L','--light', help='fix light', default=False, action="store_true")
    parser.add_option('-s','--smoke', help="smoke name", default="med")
    (opts, args) = parser.parse_args()

    tree = ET.parse(opts.input)
    sceneRoot = tree.getroot()

    if(opts.bsdf):
        BSDF(sceneRoot)
    if(opts.camera):
        medCamera(sceneRoot, opts.smoke)
    if(opts.light):
        medLight(sceneRoot, opts.smoke)


    tree.write(opts.input)
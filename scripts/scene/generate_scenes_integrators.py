"""
This script transform a given MTS 4.X files into a clay scene
"""
import os, sys, optparse
import xml.etree.ElementTree as ET
from xml.dom import minidom

def prettify(elem):
    """Return a pretty-printed XML string for the Element.
    """
    rough_string = ET.tostring(elem, 'utf-8')
    reparsed = minidom.parseString(rough_string)
    return reparsed.toprettyxml(indent="  ")

def transform_integrator(inpt, out, typeInt, listAttr, samples, samplerType, filmType):
    tree = ET.parse(inpt)
    sceneRoot = tree.getroot()
    

    # Search all nodes
    for sensorNode in sceneRoot.iter("sensor"):
        # Sampler modification if need
        for samplerNode in sensorNode.iter("sampler"):
            if(samplerType != "none"): # Change type
                samplerNode.attrib["type"] = samplerType
            
            # Go to child nodes
            for samplerSubNode in samplerNode:
                if(samplerSubNode.attrib["name"] == "sampleCount"):
                    if(samples != -1): # Change samples count
                        samplerSubNode.attrib["value"] = str(samples)

        # Film modification if needed
        if(filmType != "none"): # Change type
            for filmNode in sensorNode.iter("film"):
                filmNode.attrib["type"] = filmType

    # === Find all BSDF
    for integratorNode in sceneRoot.iter("integrator"):
        #  - Remove all attributes (Clean all)
        for integratorBSDF in integratorNode:
            integratorNode.remove(integratorBSDF)
            
        # === Change integrator type
        integratorNode.attrib["type"] = typeInt
        
        # === Add new dicts
        for typeAttr, dictV in listAttr:
            if(typeAttr == "integrator"): # Create sub integrator node
                subIntegratorNode = ET.SubElement(integratorNode, typeAttr)
                subIntegratorNode.attrib["type"] = dictV["type"]
                for typeSubAttr, dictSubV in dictV["attrs"]:
                    newNode = ET.SubElement(subIntegratorNode,typeSubAttr, dictSubV)
            else:
                newNode = ET.SubElement(integratorNode,typeAttr, dictV)
        break
    
    print("DEBUG: Write out: "+out)
    prettyTree = prettify(sceneRoot)
    with open(out, "w") as f:
        for l in prettyTree.splitlines(True):
            if "<" in l: # Check the line is not empty
                f.write(l)

def sceneDependentParameter(xmlOriFile):
    """This function is responsible to rewrite the default parameter values
    for a given scene. Indeed, some parameter are scene dependant"""

    paramFile = xmlOriFile.replace(".xml", ".param")
    if not(os.path.exists(paramFile)): # No special parameters, return empty dict
        print('[INFO] No scene dependent parameters founds')
        return {}

    # Parse paramFile
    f = open(paramFile, "r")
    lines = f.readlines()
    newValues = {}

    for l in lines:
        if ":" in l:
            name, value = l.split(":")
            name = name.strip()
            value = value.strip()
            newValues[name] = value

    return newValues

def changeAttrValue(attrs, name, newValue):
    for attr in attrs:
        # If we found a nested integrator treat it as well
        if(attr[0] == "integrator"):
            changeAttrValue(attr[1]["attrs"], name, newValue)
        else:
            if(attr[1]["name"] == name): # If the parameter name correspond, change its value
                attr[1]["value"] = newValue
                print("   = Re-write "+name)

if __name__=="__main__":
    parser = optparse.OptionParser()
    parser.add_option('-i','--input', help='input XML file')
    parser.add_option('-p','--python', help='python file for generate names')
    parser.add_option('-n', '--name', help='prefix name XML', default="scene")
    parser.add_option('-o','--output', help='output file directory')
    (opts, args) = parser.parse_args()
    
    if not(opts.input and opts.output):
        parser.error("Need input/output values")
    
    if not(os.path.exists(opts.input)):
        parser.error("Unable to find: "+str(opts.input))
    
    if not(os.path.exists(opts.python)):
        parser.error("Unable to find: "+str(opts.python))

    with open(opts.python) as f:
        code = compile(f.read(), opts.python, 'exec')
        exec(code, globals(),  locals())

    # Scene dependant parameters?
    sceneParam = sceneDependentParameter(opts.input)
    print(sceneParam)
    # If we have some, try to apply changes to all integrators
    if(sceneParam):
        print("/!\\ Scene dependant parameter /!\\")
        for integrator in INTEGRATORS:
            print(" - " + integrator["name"])
            for n,v in sceneParam.items():
                print("   * " + n + " : " + v)
                changeAttrValue(integrator["attrs"], n, v)

    techniqueNames = [""]
    print("==== Techniques")
    for integrator in INTEGRATORS:
        print(" - " + integrator["name"])
        
        # Extract samples/samplerType parameter if exist
        samples = -1
        samplerType = "none"
        filmType = "none"
        if("samples" in integrator):
            samples = int(integrator["samples"])
            print("   * Rewritting samples : "+str(samples))
        if("samplerType" in integrator):
            samplerType = integrator["samplerType"]
            print("   * Rewritting sampler type : "+samplerType)
        if("filmType" in integrator):
            filmType = integrator["filmType"]
            print("   * Rewritting film type : "+filmType)

        # Call transformation function
        transform_integrator(opts.input, opts.output+os.path.sep+opts.name+"_"+integrator["name"]+".xml", 
                             integrator["type"], integrator["attrs"], samples, samplerType, filmType)
        
        # Other stuffs
        techniqueNames.append(integrator["name"])
    
    print("=== OAR Batch")
    print(" -t ".join(techniqueNames))

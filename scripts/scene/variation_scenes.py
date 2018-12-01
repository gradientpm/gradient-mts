import os, sys, optparse

# Read the options
parser = optparse.OptionParser()
parser.add_option('-i','--input', help='input root path (directory + scene name)')
parser.add_option('-r','--remove', help='remove files not used for variation', default=False, action="store_true")
parser.add_option('-c','--config', help='config file')
(opts, args) = parser.parse_args()


import glob
import xml.etree.ElementTree as ET

def xmlEntry(type, name, value):
    return (type,{"name":name,"value":value})

def generateNewXML(xmlFile, listAttr, out):
    # Get XML tree
    tree = ET.parse(xmlFile)
    sceneRoot = tree.getroot()

    # Add the new node inside integrators
    # If the flag already exist, replace the value
    for integratorNode in sceneRoot.iter("integrator"):
        for typeSubAttr, dictSubV in listAttr:
            found = False
            for nodeSameType in integratorNode.iter(typeSubAttr):
                if ("name" in nodeSameType.attrib) \
                        and (dictSubV["name"] == nodeSameType.attrib["name"]):
                    nodeSameType.attrib = dictSubV
                    found = True

            if(not found):
                newNode = ET.SubElement(integratorNode,typeSubAttr, dictSubV)
    
    # Write out the tree
    tree.write(out)

class Change:
    def __init__(self, techniques, suffix, listAttr):
        self.techniques = techniques
        self.suffix = suffix
        self.listAttr = listAttr
        
    @staticmethod
    def XML(node):
        techniques = node.attrib["techniques"].split(";")
        suffix = node.attrib["suffix"]
        
        attribs = []
        for n in node:
            attribs.append(xmlEntry(n.tag, n.attrib["name"], n.attrib["value"]))
        return Change(techniques, suffix, attribs)
    
    def apply(self, xml):
        out = xml.replace(".xml", "_"+self.suffix+".xml")
        generateNewXML(xml, self.listAttr, out)
        return (out,out.replace(opts.input+"_","").replace(".xml",""))
    
    def __str__(self):
        s = "[ techniques: "+str(self.techniques)+", suffix: "+str(self.suffix)+"]"
        return s

def loadConfig(configFile):
    # Get XML tree
    tree = ET.parse(configFile)
    sceneRoot = tree.getroot()

    # Get all selected
    print("--- Selected:")
    selectedConfig = []
    for selected in sceneRoot.iter("Selected"):
        selectedConfig.append(selected.attrib["name"])
        print("  *", selected.attrib["name"])
    
    print("--- Changes:")
    changes = []
    for change in sceneRoot.iter("Change"):
        changes.append(Change.XML(change))
        print(changes[-1])
        
    print("--- Delete")
    deleteConfig = []
    for delected in sceneRoot.iter("Deleted"):
        deleteConfig.append(delected.attrib["name"])
        print("  *", delected.attrib["name"])
    
    return (selectedConfig, changes, deleteConfig)

# Read config files
selected, changes, deletes = loadConfig(opts.config)

# List all xml files
xmls = glob.glob(opts.input+"*.xml")
xmlsLists = []
print("Base remove: "+opts.input+"_")
for xml in xmls:
    base = xml.replace(opts.input+"_","").replace(".xml","")
    xmlsLists.append((xml, base))

if(opts.remove):
    selectedXmlsLists = []
    for xml, base in xmlsLists:
        if base in selected:
            selectedXmlsLists.append((xml, base))
        else:
            print("Remove: ", xml, "(", base, " not in ", str(selected), ")")
            os.remove(xml)
    # Copy ptr
    xmlsLists = selectedXmlsLists

for change in changes:
    for xml, base in xmlsLists:
        if base in change.techniques:
            xmlsLists.append(change.apply(xml))

xmls = glob.glob(opts.input+"*.xml")
xmlsLists = []
for xml in xmls:
    base = xml.replace(opts.input+"_","").replace(".xml","")
    xmlsLists.append((xml, base))
    
for xml, base in xmlsLists:
    if base in deletes:
        print("Delete",xml)
        os.remove(xml)

# --- Final print
xmls = glob.glob(opts.input+"*.xml")
xmlsLists = []
for xml in xmls:
    base = xml.replace(opts.input+"_","").replace(".xml","")
    xmlsLists.append((xml, base))

print("=================")
print("=== Generated ===")
print("=================")

bases = [base for xml, base in xmlsLists]
bases.sort()

print("== Number of techniques: "+str(len(bases)))

for base in bases:
    print(" *",base)
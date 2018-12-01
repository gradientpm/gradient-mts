import msetools
import optparse
import os
import glob
import math

import xml.etree.ElementTree as ET
def loadRules(configFile):
    """
    From a config file, tell for each technique which output
    need to make error analysis.
    Note that one technique can have several outputs
    :param configFile:
    :return:
    """
    # Get XML tree
    tree = ET.parse(configFile)
    sceneRoot = tree.getroot()

    # Get all selected
    print("--- Selected:")
    techniques = {}
    for tech in sceneRoot.iter("technique"):
        techniques[tech.attrib["name"]] = {}
        for prefix in tech.iter("prefix"):
            techniques[tech.attrib["name"]][prefix.attrib["name"]] = prefix.attrib["out"]
    return techniques

def findMaxIteration(filename, inputDir):
    """
    This function extract the maximum iteration number
    for a given function
    :param filename: the name of the technique
    :param inputDir: the input dir where there are the iteration
    :return:
    """
    ls = os.listdir(inputDir)
    outputFilesId = []
    for curname in ls:
        if(curname.find(filename) == 0):
            intStr = curname[len(filename):-4]
            try:
                outputFilesId.append(int(intStr))
            except ValueError:
                pass # Nothing to do, just ignore

    if(len(outputFilesId) == 0):
        print("WARN: Empty rendering !!!!!!!!!!!!! "+filename)
        return 0

    return max(outputFilesId)

def findName(outputs, tech):
    for name in outputs.keys():
        if(tech.find(name) == 0):
            return name
    return ""
    
if __name__=="__main__":
    parser = optparse.OptionParser()

    # Options in/out
    parser.add_option('-i', '--input', help="input directory", default=".")
    parser.add_option('-o','--output',help='Output directory')

    # Techniques, references and step options
    parser.add_option('-r','--reference',help='Reference image')
    parser.add_option('-s','--step', default=1,help='Images step. If only 1 image is written out of 10 computed, this argument should be 10.')
    parser.add_option('-t','--technique', help='technique name', default=[], action="append")
    parser.add_option('-A','--automatic', help='automatic technique search', default=False, action="store_true")

    # This file is quite important and need to be updated if we want to make new comparison
    parser.add_option('-c','--config',help='Configure file for mapping technique and outputs', default="rules.xml")

    # Other informations to correctly compute MSE
    parser.add_option('-e','--exposure', help='image exposure', default="0")
    parser.add_option('-m','--mask', help='image exposure', default="")
    parser.add_option('-p','--percentage', help='min percentage pixels', default="1.0")

    (opts, args) = parser.parse_args()
    steps = int(opts.step)
    finalImage = opts.reference
    mult = float(math.pow(2, float(opts.exposure)))
    percentage = float(opts.percentage)
    print("Computed with exp: ",mult)

    # Automatic detection of the techniques
    if(opts.automatic):
        outputCommand = ""
        print("==== Automatic detection ... ")
        csvFiles = glob.glob(opts.input+os.path.sep+"*_time.csv")
        for techCSV in csvFiles:
            techFiltered = os.path.basename(techCSV).replace("_time.csv","")
            print(techFiltered)
            outputCommand += "-t " + techFiltered + " "
            opts.technique.append(techFiltered)

        print("=== Manual command")
        print(outputCommand)

    # Load technique output config
    outputs = loadRules(opts.config)

    for tech in opts.technique:
        print("----------------------")
        print(tech)
        print("----------------------")
        outputName = findName(outputs, tech)

        if outputName == "":
            print("[ERROR] No output mapping is given for the technique ", tech)
            print("[ERROR] No error analysis will be perform on this technique")
            continue

        # For all the output we need to analyse
        for prefix in outputs[outputName].keys():
            # --- Construct the filename
            filename = tech+"_"+prefix+"_"
            if(prefix == ""):
                filename = tech+"_"

            # --- count how many image for the particular prefix
            nbImages = findMaxIteration(filename, opts.input)
            if(nbImages == 0):
                print("[WARN] The prefix ", prefix, "for the technique", tech)
                print("[WARN] No image is found... no analysis on this particular prefix")
                continue


            # --- create the output files
            outputCSV = []
            outputBase = tech+"_"+outputs[outputName][prefix]
            if(outputs[outputName][prefix] == ""):
                outputBase = tech

            for name in msetools.CSVNames:
                outputCSV.append(opts.output+os.path.sep+outputBase+name)

            stepsLocal = steps

            # Find nb images
            print("[INFO] Compute All metric for: ", tech, "(ref: ", opts.reference, ")")
            print("Dump info inside: "+str(outputCSV))
            filename =  opts.input + os.path.sep + filename
            msetools.computeMSEAll(filename, nbImages, stepsLocal, opts.reference, percentage, outputCSV, mult, opts.mask)

            import shutil
            inTime = opts.output+os.path.sep+tech+"_time.csv"
            outTime = opts.output+os.path.sep+outputBase+"_time.csv"
            if inTime != outTime:
                shutil.copy(inTime, outTime)
    

import os
import shutil
import glob
import math
import optparse
import sys

# === Custom import
# For read csv easily
import csv_utils


def copyFile(src, dest):
    """Just a wrapper to know what is action ...
    """
    if(not os.path.exists(src)):
        print("[WARN] Impossible to find "+ src)
    else:
        print("[INFO] Copy",src,"->",dest)
        shutil.copy(src, dest)

def extractImageTime(timeSec, technique, step, inputDir):
    """Return the image name for this technique
    with the time selected images
    """
    # === Open time .csv
    # And read the times
    timeCSVFile = inputDir+os.path.sep+technique+"_time.csv"
    times = csv_utils.extractCSVNumber(timeCSVFile)
    if(len(times) < step):
        print("No enought iteration done by the technique")
        print("Or no time file")
        return []

    # === Compute cumul times
    # Indeed inside the CSV file,
    # there is only time per iteration
    cumulTimes = []
    currentTime = 0
    for t in times:
        currentTime += t
        cumulTimes.append(currentTime)
    
    # === Find split for each times
    timesSplit = []
    for tsec in timeSec:
        # For each time split we want
        # we do a linear seach inside the timing to get the closest image
        # which have the time upper that the time limits
        found = False
        lastID = step
        for i in range(step-1,len(cumulTimes),step):
            if(tsec < cumulTimes[i]):
                timesSplit.append(i+1)
                print("[INFO] Time diff for ", tsec, "is equal to", cumulTimes[i] - tsec)
                found = True
                break # pass next
            lastID = i
        # It can happens that in some case, it is not possible
        # to find a image that is upper the time limit
        # In this case, we take the last one
        if(not found):
            print("[WARN] Force find !",  technique, "(", lastID, ")")
            timesSplit.append(lastID+1)
            print("[INFO] Time diff for ", tsec, "is equal to", cumulTimes[lastID] - tsec)
            
    # Now, we know each iteration that correspond to the time split
    # Here just
    return timesSplit

def extractAndCopyImageTime(output, timeSec, tech, step, inputDir, names):
    """
    This technique
    :param output: output dir informations
    :param timeSec: the time splitting
    :param tech: the technique name
    :param step: image step use for all the techniques
    :param inputDir: where the images are
    :param names: the additional name of all the images
    :return:
    """

    # Get all the iteration related to the technique
    iterations = extractImageTime(timeSec, tech, step, inputDir)
    if(len(iterations) == 0):
        print("WARN: Somethings goes wrong for technique: "+tech+", SKIP IT")
        return # Do nothing !

    print(tech, " iterations: ", str(iterations))

    # for all the name, make a copy of it
    # in the same time, replace the iteration number
    # by the time split value
    for j in range(len(timeSec)):
        for name in names:
            nameBase = tech+"_"+name+"_"
            if(name == ""):
                # This is a special case, in this case, change the nameBase by removing double __
                nameBase = nameBase.replace("__", "_")
            copyFile(inputDir+os.path.sep+nameBase+str(iterations[j])+".hdr",
                     output+os.path.sep+nameBase+str(timeSec[j])+".hdr")
                    
if __name__=="__main__":
    parser = optparse.OptionParser()

    # Input/Output directories
    parser.add_option('-i', '--input', help="input directory", default=".")
    parser.add_option('-o','--output',help='Output directory')

    # Options for the image step and which technique we want to pack the results
    parser.add_option('-s','--step', default=1,help='Images step. If only 1 image is written out of 10 computed, this argument should be 10.')
    parser.add_option('-t','--technique', help='technique name', default=[], action="append")
    parser.add_option('-A','--automatic', help='automatic technique search', default=False, action="store_true")
    parser.add_option('-f','--time', help='time image [XX:XX:XX]', default=[], action="append")

    # Option to tell which name file we want to pack
    parser.add_option('-n','--name', help='files names (pass, gX, etc.)', default=["pass"], action="append") # Options for computing the reference and all other values
    parser.add_option('-r','--reference',help='Reference image')

    (opts, args) = parser.parse_args()

    # Read all the input provided by the user
    step = int(opts.step)
    reference = opts.reference
    
    # === Automatic technique founding
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

    # === Create output directory
    if(not os.path.exists(opts.output)):
        os.makedirs(opts.output)
        
    # --- Copy all into res directory
    # Timing informations and log files (if it's exist)
    for tech in opts.technique:
         TimeFile = tech+"_time.csv"
         copyFile(opts.input + os.path.sep +TimeFile, opts.output+os.path.sep+TimeFile)
         OutFile = tech + ".out"
         copyFile(opts.input + os.path.sep +OutFile, opts.output+os.path.sep+OutFile)

    # === Copy the ref also
    if(opts.reference):
        #copyFile(opts.input + os.path.sep +reference, opts.output+os.path.sep+"ref.hdr")
        copyFile(opts.reference, opts.output+os.path.sep+"Ref.hdr")
    
    # === Extract step images (time)
    # --- Convert time
    timeSec = []
    for t in opts.time:
        tSplit = t.split(":")
        hours = int(tSplit[0])
        minutes = int(tSplit[1])
        seconds = int(tSplit[2])
        timeSec.append(hours*3600 + minutes*60 + seconds)

    print("Time Seconds:", timeSec)
    
    
    # --- Extract image steps with techniques
    # HACK: empty prefix becomes the double quote on Windows, so work around by adding a real empty string here
    opts.name += ['']
    for tech in opts.technique:
        stepsLocal = step
        if(tech == "GBDPT_L1" or tech == "GBDPT_L2" or
           tech == "GPT_L1" or tech == "GPT_L2"):
            print("[WARN] Change the step")
            stepsLocal = 1

        extractAndCopyImageTime(opts.output, timeSec, tech, stepsLocal, opts.input, opts.name)


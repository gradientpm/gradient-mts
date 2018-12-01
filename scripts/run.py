import optparse
import os
import sys

PYTHON_NAME="python3"

# --- Read all params
parser = optparse.OptionParser()

# Input/Output options
parser.add_option('-m','--mitsuba', help='mitsuba exe', default="")
parser.add_option('-i','--input', help='input dir', default="")
parser.add_option('-s','--scenename', help='scene name')
parser.add_option('-o','--output', help='output name')

# Timing to check for different timing
parser.add_option('-f','--time', help='time image [XX:XX:XX]', default=[], action="append")

# XML Options
parser.add_option('-v','--variation', help='xml variation if required', default="")

# Running/Rendering Options
parser.add_option('-j','--jobs', help="number of threads", default="1")

# Packing Options
parser.add_option('-r','--reference', help="reference image", default="")
parser.add_option('-l','--layout', help="html page layout (xml file)", default="")
parser.add_option('-d','--dumpiter', help="dump interval for each images", default="1")

parser.add_option('-g','--xmlgen', help="enable scene generation (XML)", default=False, action="store_true")
parser.add_option('-c','--compute', help="enable scene running mts (Rendering)", default=False, action="store_true")
parser.add_option('-p','--pack', help="enable pack results (Pack, MSE)", default=False, action="store_true")
parser.add_option('-P','--packhtml', help="enable pack results (HTML)", default=False, action="store_true")
parser.add_option('-R','--avoidMETRIC', help="disable metric computation", default=False, action="store_true")
parser.add_option('-H','--hackMETRIC', help="enable hack metric", default=False, action="store_true")
parser.add_option('-C','--cluster', help="enable cluster running mode", default=False, action="store_true")

parser.add_option('-e','--exposure', help="exposure value", default="0")
parser.add_option('-O','--htmlsuffix', help="html suffix", default="0")

parser.add_option('-G','--generator', help="python file for generating XML", default="scene/generatorSTAR.py")

# Technique
parser.add_option('-t','--technique', help="technique", default="")
 
# Read input
(opts, args) = parser.parse_args()

#########################
## Error checking
if(opts.input == ""):
    print("[ERROR] Need to where the scene are (-i)")
    parser.print_help()
    sys.exit(1)

#########################
def launch(command):
    import subprocess
    print(("[DEBUG] Launch: "+" ".join(command)))
    process = subprocess.Popen(command,shell=False)
    return process.communicate()

#########################
# Step 1: Re-generate scene
#########################scenes
currentScenePath = opts.input + os.path.sep + opts.scenename + "_scene"
if(opts.xmlgen):
    command = [PYTHON_NAME, "scene/generate_scenes_integrators.py",
               "-i", currentScenePath + os.path.sep + "ori_" + opts.scenename + ".xml",
               "-p", opts.generator,
               "-n", opts.scenename,
               "-o", currentScenePath]
    launch(command)


    if(opts.variation != ""):
        # We want a variation, run it
        command = [PYTHON_NAME, "scene/variation_scenes.py",
                   "-i", currentScenePath + os.path.sep + opts.scenename,
                   "-r", # Remove
                   "-c", opts.variation]
        launch(command)

    if not (opts.compute or opts.pack or opts.packhtml):
        print("Finish :)")
        sys.exit(0)

#########################
# Step 2: Run the computation
#########################
if(opts.mitsuba == ""):
    print("[ERROR] Need to provide mitsuba binary path (-m)")
    parser.print_help()
    sys.exit(1)

if(len(opts.time) == 0):
    print("No time is provided which is needed for the rest of the algorithm")
    parser.print_help()
    sys.exit(1)

# - Convert time sec
timingsSec = []
for t in opts.time:
    tSplit = t.split(":")
    hours = int(tSplit[0])
    minutes = int(tSplit[1])
    seconds = int(tSplit[2])
    timingsSec.append(hours*3600 + minutes*60 + seconds)
maxTimeSec = max(timingsSec)

currentOutMTS = currentScenePath + os.path.sep + "out" + opts.output
currentOutRES = currentScenePath + os.path.sep + "res" + opts.output
currentOutHTML = currentScenePath + os.path.sep + "html" + opts.output + opts.htmlsuffix

if opts.technique == "":
    #techFlag = ["-t", "GPM_L2", "-t", "SPPM"]
    #techFlag = ["-t", "GPT", "-t", "GBDPT"]
    techFlag = ["-A"]
else:
    techFlag = ["-t", opts.technique]
	
if(opts.compute):
    command = [PYTHON_NAME, "run/run_batch.py",
               "-m", opts.mitsuba, # Use this MTS
               "-s", str(maxTimeSec), # Rendering time
               "-i", currentScenePath + os.path.sep + opts.scenename,
               "-o", currentScenePath + os.path.sep + "out" + opts.output,
               "-j", opts.jobs
               ]
    if(opts.cluster):
        print("[CLUSTER MODE]")
        command += ["-C"]
    command += techFlag
    launch(command)

#########################
# Step 3: Packing...
#########################
if(opts.pack):
    #TODO: Check reference
    command = [PYTHON_NAME, "results/run_pack.py",
               "-i", currentOutMTS,
               "-o", currentOutRES,
               "-r", opts.reference,
               "-s", opts.dumpiter,
               "-n", "pass",
               "-n", '""',
               "-n", "dxAbs",
               "-n", "dyAbs",
               "-n", "L1",
               "-n", "L2",
               "-n", "L2Weighted",
               "-n", "Uni",
               "-n", "Weighted",
               "-n", "throughput",
               "-n", "recons"]
    for t in opts.time:
        command += ["-f", t]
    command += techFlag
    launch(command)

    if(not opts.avoidMETRIC):
        pourcentage = "1.0"
        if(opts.hackMETRIC):
            pourcentage = "0.9999"
        command = [PYTHON_NAME, "results/run_mse.py",
                   "-i", currentOutMTS,
                   "-o", currentOutRES,
                   "-r", opts.reference,
                   "-s", opts.dumpiter,
                   "-p", pourcentage,
                   "-e", opts.exposure,
                   "-c", "results/rules.xml"]
        command += techFlag
        launch(command)

#########################
# Step 4: Packing...
#########################
if(opts.packhtml):
    #TODO: Check reference
    #TODO: Check layout
    command = [PYTHON_NAME, "results/run_html.py",
               "-j", "./results/data/js",
               "-i", currentOutRES,
               "-o", currentOutHTML,
               "-r", opts.reference,
               "-t", str(maxTimeSec), # FIXME: Time
               "-c", opts.layout,
               "-e", opts.exposure,
               "-s", opts.dumpiter,
               "-n", opts.scenename,
               "-m"]
    print(launch(command))

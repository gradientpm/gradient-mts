import optparse
import os
import glob
import sys
import shutil
"""
The aim of this script is to copy some technique results (matched with the patern)
of differents scenes to destination result dir
With the possiblity to rename the technique 

python copy_tech.py -i ../../scenes -s cboxSmoke -o June23 -d SigRes -p "GVPM_L2_a2m_distance\S*" -b old
-f : force
-c : do the copy operation
"""

parser = optparse.OptionParser()
parser.add_option('-i','--input', help='input dir', default="")
parser.add_option('-s','--scenes', help='scene name', default=[], action="append")
parser.add_option('-p','--pattern', help='pattern name', default="")
parser.add_option('-o','--origin', help='origin name results')
parser.add_option('-d','--dest', help='dest name results')
parser.add_option('-b','--back', help="back prefix", default="")
parser.add_option('-c','--copy', help="do the copy operation", default=False, action="store_true")
parser.add_option('-f','--force', help="force the copy by ignoring the collisions", default=False, action="store_true")

(opts, args) = parser.parse_args()

# Arguments checking
if(opts.origin == ""):
    print("[ERROR] need to specify an origin")
    sys.exit(1)
if(opts.dest == ""):
    print("[ERROR] need to specify an dest")
    sys.exit(1)
if(opts.origin == opts.dest):
    print("[ERROR] Origin and Dest need to be different")

import re
p = re.compile(opts.pattern)

for s in opts.scenes:
    # Directory checking
    orgPath = opts.input + os.path.sep + s + "_scene" + os.path.sep + "res" + opts.origin + os.path.sep
    destPath = opts.input + os.path.sep + s + "_scene" + os.path.sep + "res" + opts.dest + os.path.sep
    if(not os.path.exists(orgPath)):
        print("[WARN] the origin res dir need to exist: %s" % orgPath)
        continue
    if(not os.path.exists(destPath)):
        print("[WARN] the dest res dir need to exist: %s" % destPath)
        continue

    # List all techniques
    allTechName = glob.glob(orgPath+os.path.sep+"*_time.csv")
    allTechName = [n.replace(orgPath,"").replace("_time.csv", "") for n in allTechName]
    allTechName.sort()

    print("List of techniques matches: ")
    allTechNameMatches = []
    for t in allTechName:
        if(p.match(t)):
            allTechNameMatches.append(t)
            print(" *", t)

    if(len(allTechNameMatches) == 0):
        print("[INFO] Empty match, skip")
        continue

    if(opts.back != ""):
        allTechNameMatchesRename = [n+"_"+opts.back for n in allTechNameMatches]
    else:
        allTechNameMatchesRename = allTechNameMatches
    # Check that the technique is not already in place (in dest)
    allTechNameDest = glob.glob(destPath+os.path.sep+"*_time.csv")
    allTechNameDest = [n.replace(destPath,"").replace("_time.csv", "") for n in allTechNameDest]
    allTechNameDest.sort()

    collision = False
    for t in allTechNameMatchesRename:
        if t in allTechNameDest:
            print("[WARN] Collision ! ", t)
            collision = True

    if(collision and (not opts.force)):
        print("[ERROR] Collision detected, skip")
        continue

    print("[INFO] Copy...")
    filesMoved = []
    for i in range(len(allTechNameMatches)-1, -1, -1):
        t = allTechNameMatches[i]
        tRename = allTechNameMatchesRename[i]

        currFiles = glob.glob(orgPath+os.path.sep+t+"*")
        print(" * ", t)
        for f in currFiles:
            if f in filesMoved:
                continue
            else:
                filesMoved.append(f)
                newName = destPath+os.path.sep+f.replace(orgPath,"").replace(t, tRename)
                print("    - ", f, ' -> ', newName)
                if(opts.copy):
                    shutil.copyfile(f, newName)
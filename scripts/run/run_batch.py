import optparse
import glob
import sys
import os
import subprocess
import datetime


def kill(proc_pid):
    process = psutil.Process(proc_pid)
    print(process)
    for proc in process.children(recursive=True):
        proc.kill()
    process.kill()


# TODO: Use it as a parameters
CLUSTER_MODE = "todai"

if __name__ == "__main__":
    # Options
    parser = optparse.OptionParser()
    parser.add_option("-m", "--mitsuba", help="the mitsuba .sh or MTS executable")
    parser.add_option('-t', '--technique', help='technique name', default=[], action="append")
    parser.add_option('-s', '--time', help='time of running (in sec)', default=None)
    parser.add_option('-i', '--input', help='input scene name (%path_scene%_%tech%.xml)')
    parser.add_option('-o', '--output', help='output directory')
    parser.add_option('-j', '--jobs', help='nb thread per jobs', default=12)
    parser.add_option('-A', '--automatic', help='replace -t usage by finding all files', default=False,
                      action="store_true")
    parser.add_option('-C', '--cluster', help='submit the job in the cluster', default=False, action="store_true")
    (options, args) = parser.parse_args()

    # === Error handling
    if (options.input == None):
        print("\nError: Need input values\n")
        parser.print_help()
        sys.exit(1)
    if (options.output == None):
        print("\nError: Need output values\n")
        parser.print_help()
        sys.exit(1)
    if (options.time == None):
        print("\nError: Need to specify time\n")
        parser.print_help()
        sys.exit(1)

    print("Techniques input: ", options.technique)

    # === Create output directory
    if (not os.path.exists(options.output)):
        os.makedirs(options.output)

    # List all files we need to render
    # Note that if automatic mode is used. Manually input techniques are not taken into account
    if (len(options.technique) != 0 and options.automatic):
        print(""""[WARN] Automatic mode and manually technique input is detected. 
        Just ignore manual input and usefully automatic detection""")

    # Technique array will store all techniques names
    techniques = []
    if (options.automatic):
        print(" === Automatic finding ...")
        files = glob.glob(options.input + "*.xml")
        for fileXML in files:
            filename = os.path.basename(fileXML)
            filenameRoot = os.path.basename(options.input)
            tech = filename[len(filenameRoot) + 1:-4]
            print("   * Find: %s" % tech)
            techniques.append(tech)
    else:
        # If we use manual, just copy technique names
        for tech in options.technique:
            techniques.append(tech.strip())

    dateStart = datetime.datetime.now()
    dateEnd = dateStart + datetime.timedelta(seconds=int(options.time) * len(techniques))
    for tech in techniques:
        print("[INFO] End of rendering tasks: ", "{:%H:%M}".format(dateEnd))
        command = [options.mitsuba, "-p", str(options.jobs)]
        # Add output
        command += ["-o", options.output + os.path.sep + tech]
        # Add input
        command += [options.input + "_" + tech + ".xml"]
        # Remove the progress bar
        command += ["-z"]

        if (options.cluster):
            if CLUSTER_MODE == "mcgills":
                if (int(options.jobs) != 12):
                    print("With mcgill cluster, only 12 thread at a time is possible.")
                    raise

                ## Cluster which use torque
                ## target only westmere
                TMP_FILE = "/home/neodym60/scratch/cluster_tmp_command.sh"
                fileCommand = open(TMP_FILE, "w")
                fileCommand.write("#!/bin/sh\n")
                fileCommand.write("\n")

                # Number of threads
                # force to use all the thread on westmere (SW) nodes
                fileCommand.write("#PBS -l nodes=1:ppn=12:westmere\n")
                fileCommand.write("#PBS -A YOURCODE\n")

                # Redirect regular output (+error)
                fileCommand.write("#PBS -o "+options.output + os.path.sep + tech + ".out\n")
                fileCommand.write("#PBS -e "+options.output + os.path.sep + tech + ".err\n")                

                # Write the time
                minutesTask = int(options.time) // 60
                secondsTaks = int(options.time) % 60
                fileCommand.write("#PBS -l walltime=00:" + str(minutesTask) + ":" + str(secondsTaks) + "\n")

                # Add the task (mitsuba call)
                fileCommand.write("\n")
                fileCommand.write('MITSUBA_DIR="' + str(os.path.dirname(options.mitsuba)) + '"\n')
                fileCommand.write('export LD_LIBRARY_PATH="$MITSUBA_DIR:$LD_LIBRARY_PATH"\n')
                fileCommand.write('export PATH="$MITSUBA_DIR:$PATH"\n')
                fileCommand.write(" ".join(command) + "\n")

                # Run the file on the cluster
                fileCommand.close()
                clusterCommand = ["qsub", TMP_FILE]
                # process = subprocess.run(clusterCommand,shell=False,stdout=subprocess.PIPE)
                proc = subprocess.check_output(clusterCommand, shell=False)
                print(proc)
            elif CLUSTER_MODE == "todai":
                if (int(options.jobs) != 18):
                    print("With todai cluster, only 18 thread at a time is possible. (for now)")
                    raise

                ## This are the command for the hand made cluster
                ## based on slurm
                print("[CLUSTER] Submit the job")
                TMP_FILE = "/home/u00068/cluster_tmp_command.sh"
                fileCommand = open(TMP_FILE, "w")
                fileCommand.write("#!/bin/sh\n")
                fileCommand.write("\n")

                # Write time limit
                minutesTask = int(options.time) // 60
                secondsTaks = int(options.time) % 60
                fileCommand.write("#SBATCH --time=" + str(minutesTask) + ":" + str(secondsTaks) + "\n")

                # Number of threads and the target node
                fileCommand.write("#SBATCH --ntasks=1\n")
                fileCommand.write("#SBATCH --cpus-per-task=18\n")
                fileCommand.write("#SBATCH --partition=p\n")

                # Redirect the outputs
                fileCommand.write("#SBATCH -o "+options.output + os.path.sep + tech + ".out\n")
                fileCommand.write("#SBATCH -e "+options.output + os.path.sep + tech + ".err\n")                

                
                # The main command
                fileCommand.write("\n")
                fileCommand.write("srun " + " ".join(command) + "\n")

                # Run the file on the cluster
                fileCommand.close()
                clusterCommand = ["sbatch", TMP_FILE]
                process = subprocess.run(clusterCommand, shell=False, stdout=subprocess.PIPE)
            elif CLUSTER_MODE == "custom":
                ## This are the command for the hand made cluster
                ## based on slurm
                print("[CLUSTER] Submit the job")
                TMP_FILE = "/tmp/cluster_tmp_command.sh"
                fileCommand = open(TMP_FILE, "w")
                fileCommand.write("#!/bin/sh\n")
                fileCommand.write("\n")

                # Write time limit
                minutesTask = int(options.time) // 60
                secondsTaks = int(options.time) % 60
                fileCommand.write("#SBATCH --time=" + str(minutesTask) + ":" + str(secondsTaks) + "\n")

                # Number of threads
                fileCommand.write("#SBATCH --ntasks=1\n")
                fileCommand.write("#SBACTH --cpus-per-task=12\n")
                
                # The main command
                fileCommand.write("\n")
                fileCommand.write("srun " + " ".join(command) + "\n")

                # Run the file on the cluster
                fileCommand.close()
                clusterCommand = ["sbatch", TMP_FILE]
                process = subprocess.run(clusterCommand, shell=False, stdout=subprocess.PIPE)
                print(process)
            else:
                print("bad cluster configuration selection")
                raise
        else:
            print("[DEBUG] Run computation for", tech)
            # Start foo as a process
            try:
                proc = subprocess.check_output(command, shell=False, timeout=int(options.time))
            except subprocess.TimeoutExpired as e:
                outputFile = open(options.output + os.path.sep + tech + ".out", "w")
                outputFile.write(e.output.decode("utf-8"))
                outputFile.close()
                print("Process killed :)")
            except subprocess.CalledProcessError as e:
                outputFile = open(options.output + os.path.sep + tech + ".out", "w")
                outputFile.write(e.output.decode("utf-8"))
                outputFile.close()
                print("Error during the rendering :(")

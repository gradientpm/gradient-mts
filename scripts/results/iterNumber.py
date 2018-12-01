import optparse
import glob
import os

def nbIter(filename, timeLimit):
    if(timeLimit == 0.0):
        return len(open(filename, "r").readlines())
    else:
        lines = open(filename, "r").readlines()
        times = [float(v.replace(",\n", "")) for v in lines]
        cum, nb = 0, 0
        while(cum < timeLimit and nb < len(times)):
            cum += times[nb]
            nb += 1
        return nb

if __name__=="__main__":
    parser = optparse.OptionParser()
    parser.add_option('-i','--input', help="input directory", default="")
    parser.add_option('-f','--time', help="time", default="600")
    parser.add_option('-n','--name', help="name technique", default="")
    (opts, args) = parser.parse_args()

    time = float(opts.time)
    csvFiles = glob.glob(opts.input+os.path.sep+"*_time.csv")
    for techCSV in csvFiles:
        if opts.name in techCSV:
            f = open(techCSV, "r")
            print("======================")
            print("======================")
            print(techCSV)
            print("Nb iterations: ", nbIter(techCSV, time))
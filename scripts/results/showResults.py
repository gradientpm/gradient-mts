from pylab import *
import matplotlib.pyplot as plt
import optparse
import csv_utils
import os
import glob
import logging

logger = logging.getLogger(__name__)

colors = ['darkslateblue', 'darkorange', 'rosybrown', 'lightcoral', 'magenta', 'indianred', 'orangered', 'indigo', 'yellow', 'darkturquoise', 'palevioletred', 'peru', 'sienna', 'lightgreen', 'orchid', 'lawngreen', 'steelblue', 'aquamarine', 'darkviolet', 'olive', 'cadetblue', 'darkgreen', 'darkslategray', 'red', 'black', 'darksalmon', 'lime', 'grey', 'aqua', 'violet', 'lightseagreen', 'whitesmoke', 'tan', 'dodgerblue', 'darkblue', 'deeppink', 'dimgrey', 'mediumseagreen', 'midnightblue', 'darkgray', 'purple', 'slateblue', 'blue', 'green', 'mediumblue', 'lightpink', 'greenyellow', 'mediumaquamarine', 'firebrick', 'deepskyblue', 'lightsteelblue', 'olivedrab', 'blueviolet', 'lightgray', 'mediumturquoise', 'mediumspringgreen', 'coral', 'saddlebrown', 'lemonchiffon', 'cornflowerblue', 'darkkhaki', 'springgreen', 'mediumvioletred', 'darkseagreen', 'burlywood', 'royalblue', 'limegreen', 'mediumorchid', 'darkorchid', 'darkgoldenrod', 'khaki', 'darkolivegreen', 'mediumpurple', 'hotpink', 'lightskyblue', 'mediumslateblue', 'thistle', 'darkred', 'gainsboro', 'teal', 'turquoise', 'yellowgreen', 'seagreen', 'oldlace', 'silver', 'chocolate', 'crimson']

class Technique(csv_utils.Technique):
    def __init__(self, name, color, xCSVFile, yCSVFile, step, log,shift=0, column=0, lineStyle="-"):
        csv_utils.Technique.__init__(self, name, color, xCSVFile, yCSVFile, step, log,shift=shift, column=column)
        self.lineStyle = lineStyle
        
    def addPlot(self, ax, fitting = 1, fittingClamp = 0.5):        
        if(fitting != 0):
            clampIndices = int(len(self.x)*fittingClamp)
            logger.info("Clamping indices: "+str(clampIndices))
            x = self.x[clampIndices:]
            y = self.y[clampIndices:]
            fit = polyfit(x,y,fitting)
            logger.info("Fitting: "+str(fit))
            fit_fn = poly1d(fit)
            ax.plot(self.x, fit_fn(self.x), 'r--', color=self.color)
            if(fitting == 1):
                self.name += " {:.2f}".format(fit[0])
        logger.info(" === Plot "+str(self.name))
        ax.plot(self.x, self.y, color=self.color, linestyle=self.lineStyle, label=self.name)
        
    def dump(self):
        logger.info("Technique "+self.name)
        logger.info("x length: " + str(len(self.x)))
        logger.info("y length: " + str(len(self.y)))
        
    def isValid(self):
        return len(self.y) > 0
    
    def generatePairData(self):
        temp = []
        for i in range(len(self.x)):
            temp.append([self.x[i], self.y[i]])
        return temp
    
    def clampTime(self, value):
        if self.x[-1] < value:
            return
        
        i = 0
        while self.x[i] < value:
            i += 1
        logger.debug(str(i))
        self.x = self.x[:i]
        self.y = self.y[:i]
        
    def generateConstantDataX(self):
        temp = []
        prev = 0
        for i in range(len(self.x)):
            temp.append([i, self.x[i] - prev])
            prev = self.x[i]
        return temp
    
    def generateConstantDataXLog(self):
        temp = []
        prev = 0
        for i in range(len(self.x)):
            temp.append([math.log10(self.x[i]), self.x[i] - prev])
            prev = self.x[i]
        return temp
    
    def jsEntry(self):
        temp = []
        for i in range(len(self.x)):
            temp.append([self.x[i], self.y[i]])
        return "{data: " + str(temp) + ', label: "' + self.name + '"}'

def getTechniqueNames(rep):
    list = []
    logger.info("==== Automatic detection ... ")
    csvFiles = glob.glob(rep+os.path.sep+"*_time.csv")
    for techCSV in csvFiles:
        techFiltered = os.path.basename(techCSV).replace("_time.csv","")
        logger.debug(techFiltered)
        list.append(techFiltered)
    return list

linesStyles = ["-", "."]

def readAllTechniques(list, rep, step, uselog=False, basex="_time.csv", basey="_rmse.csv"):
    linesStylesIndex = 0
    colorIndex = 0
    techniques = []
    for tech in list:
        # Extract
        splitLine = tech.split(',')
        name = splitLine[0]
        timename = name
        # If another name for time, take it into account
        if len(splitLine) > 1:
            timename = splitLine[1]
        logger.info("[READ ALL]: "+str(name)+" time="+str(timename))

        stepsLocal = step
        if(timename == "GBDPT_L1" or timename == "GBDPT_L2" or
           timename == "GPT_L1" or timename == "GPT_L2"):
            stepsLocal = 1

        color = colors[colorIndex]
        colorIndex += 1
        lineStyle = linesStyles[linesStylesIndex]
        shift = 0
        column = 0

        #No shift and column selected
        # if len(splitLine) > 2:
        #     column = int(splitLine[2])
        # if len(splitLine) > 3:
        #     shift = int(splitLine[3])
            
        xCSVFile = rep+os.path.sep+timename+basex
        yCSVFile = rep+os.path.sep+name+basey
        
        stepFile = stepsLocal

        # Create & Append
        techCurr = None
        if uselog:
            techCurr = Technique(name, color, xCSVFile, yCSVFile, stepFile, True, shift=shift, column=column, lineStyle=lineStyle)
        else:
            techCurr = Technique(name, color, xCSVFile, yCSVFile, stepFile, False, shift=shift, column=column, lineStyle=lineStyle)
        
        # === Check if it's valid or not
        if(techCurr.isValid()):
            techniques.append(techCurr)
        else:
            logger.warn("Technique invalid: "+name)
    return techniques

def createJSScript(placeName, techniques, colorDict = None):
    totalText = '$(function() { $.plot("'+placeName+'",['
    # === fill the data section
    for i in range(len(techniques)):
        if techniques[i].name in colorDict:
            totalText += techniques[i].jsEntry()
            if((i+1) != len(techniques)):
                totalText += ",\n"
    totalText += "]" 
    
    if(colorDict):
        totalText += ","
        totalText += "{colors:["
        for i in range(len(techniques)):
            if techniques[i].name in colorDict:
                c = colorDict[techniques[i].name] 
                totalText += '"'+c+'", '
        totalText += "]}"
    # close plot function
    totalText += ");\n"
    # close js call
    totalText += "});"
    return(totalText)

if __name__=="__main__":
    parser = optparse.OptionParser()
    parser.add_option('-o','--output', help="output pdf figure", default=".")
    parser.add_option('-r','--rep', help="repertory of files", default=".")
    parser.add_option('-x','--xname', help="Name x axis", default = "")
    parser.add_option('-y','--yname', help="Name y axis", default = "")
    parser.add_option('-b','--basey', help="Base name file y", default="_rmse.csv")
    parser.add_option('-s','--step', default=1,help='Images step. If only 1 image is written out of 10 computed, this argument should be 10.')
    parser.add_option('-i','--input', help='Format: <name, color, column, shift>', default=[], action="append")
    parser.add_option('-n','--normalize', help='Normalize the X absisse', default=False, action="store_true")
    parser.add_option('-l','--log', help="To create log plot", default=False, action = "store_true")
    parser.add_option('-R','--reorganise', default=False, action = "store_true")
    parser.add_option('-A','--automatic', default=False, action = "store_true")
    parser.add_option('-T','--timeNormalisation', default=False, action = "store_true")
    parser.add_option('-f','--fitting', default="0")
    parser.add_option('-F','--fittingClamp', default="0.1")
    parser.add_option('-m','--maxpass', default=-1)
    parser.add_option('-J','--javascript', default="")
    
    (options, args) = parser.parse_args()
    options.maxpass = int(options.maxpass)
    # === Convert input
    step = int(options.step)
    
    fittingClamp = min(max(float(options.fittingClamp),0.0),1.0)
    fittingPoly = int(options.fitting)
    
    if(options.automatic):
        options.input.extend(getTechniqueNames(options.rep))
    
    techniques = readAllTechniques(options.input, options.rep, step, options.log, "_time.csv", options.basey)
    
    if(options.javascript != ""):
        text = createJSScript(options.javascript, techniques)
        logger.info(text)
        
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)

    if(options.timeNormalisation):
        for tech in techniques:
            for i in range(len(tech.x)):
                tech.x[i] = i
    
    if options.normalize:
        
        # === Search min X
        minX = techniques[0].x[-1]
        for tech in techniques:
            minX = min(tech.x[-1], minX)
        logger.info("[INFO] Min x axis found: "+str(minX))
        
        # === Cut down the DATA
        for tech in techniques:
            newX = []
            for x in tech.x:
                if x > minX:
                    break
                newX.append(x)
            tech.x = newX
            tech.y = tech.y[0:len(tech.x)]

    
    ax.set_ylabel(options.yname)
    ax.set_xlabel(options.xname)
    
    if(options.maxpass != -1):
        # === Cut the axis ===
        for tech in techniques:
            tech.x = tech.x[:options.maxpass]
            tech.y = tech.y[:options.maxpass]
    
    if(options.reorganise):
        conv = {}
        for i,tech in enumerate(techniques):
            conv[tech.y[-1]] = tech 
        keysSorts = list(conv.keys())
        keysSorts.sort()
        keysSorts.reverse()
        logger.debug(conv)
        
        techniques = [conv[k] for k in keysSorts]

    logger.info("============== Plot Techniques ... ")
    for tech in techniques:
        tech.dump()
        tech.addPlot(ax, fittingPoly, fittingClamp)

    logger.info("=== Technique Order:")
    for pos,tech in enumerate(techniques):
        logger.info("%i: Tech: %s | RMSE: %f" % (pos,tech.name,tech.y[-1]))

    plt.legend()
    if(options.output != ""):
        plt.savefig( options.output, format='pdf' )
    else:
        plt.show()


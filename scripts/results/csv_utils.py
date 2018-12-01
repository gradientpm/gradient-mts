import math
import logging

logger = logging.getLogger(__name__)

def extractCSVColumn(filename):
    """ Extract the number of column 
    for the CSV file ..
    """
    fileCSV = open(filename, "r")
    lines = fileCSV.readlines()
    if(len(lines) > 1):
        line = lines[1]
        return len(line.split(','))
    return 0

def isCSVHeader(filename):
    """Test if the CSV file begun with
    an # symbol"""
    try:
        fileCSV = open(filename, "r")
        lines = fileCSV.readlines()
        if(len(lines) > 0):
            if(lines[0][0] == '#'):
                logger.debug("Found header :)")
                return True
        return False
    except:
        return False

def extractCSVHeader(filename):
    """ Extract the number of column 
    for the CSV file ..
    """
    try: 
        fileCSV = open(filename, "r")
        lines = fileCSV.readlines()
        if(len(lines) > 0):
            line = lines[0][1:]
            return line.split(',')
        return []
    except:
        logger.error("CSV file header error: " + filename)
        return []


def extractCSVNumberLog(filename, mul, column=0):
    """ The CSV file format need to be, somethings like:
        <number>,
        <number>,
        ...
        eof
    """
    beginEntry = 0
    if(isCSVHeader(filename)):
        beginEntry += 1 
    numbers = []

    try:
        fileCSV = open(filename, "r")
        lines = fileCSV.readlines()
        for line in lines:
            line = line.strip(' \r\n') # Remove all unessary symbols
            if not line:               # Skip blank lines
                continue
            if(line[-1] != ","):
                raise Exception('Error format CSV handling', line)
            num = float(line.split(',')[column])
            num = math.log10( num )
            numbers.append( mul*num ) # Remove the , character

        return numbers
    except:
        logger.error("CSV file IO (logarithm) error: " + filename)
        return []

def extractCSVNumber(filename, column=0):
    """ The CSV file format need to be, somethings like:
        <number>,
        <number>,
        ...
        eof
    """
    beginEntry = 0
    if(isCSVHeader(filename)):
        beginEntry += 1 
        
    numbers = []
    try:
        fileCSV = open(filename, "r")
        lines = fileCSV.readlines()
        for line in lines:
            line = line.strip(' \r\n') # Remove all unessary symbols
            if not line:               # Skip blank lines
                continue
            if(line[-1] != ","):
                raise Exception('Error format CSV handling', line)
            numbers.append(float(line.split(',')[column])) # Remove the , character
        return numbers
    except:
        logger.error("CSV file IO error: " + filename)
        return []

class Technique:
    def __init__(self, name, color, xCSVFile, yCSVFile, step, log, shift=0, column=0):
        self.name = name
        self.color = color
        self.shift = shift
        
        if log:
            self.x = extractCSVNumber(xCSVFile)
            self.y = extractCSVNumberLog(yCSVFile, 20, column)
        else:
            self.x = extractCSVNumber(xCSVFile)
            self.y = extractCSVNumber(yCSVFile, column)
        
        self.step = step
        
        if self.shift > 0:
            self.x = self.x[0:-self.shift]
            self.y = self.y[self.shift::]
        
        # Compute the cumulate time
        currentX = 0.0
        for i in range(len(self.x)):
            currentX = currentX + self.x[i]
            if(log):
                self.x[i] = math.log10( currentX )
            else:
                self.x[i] = currentX
            
        Xconcat = []
        for i in range(0, len(self.x), self.step):
            Xconcat.append(self.x[i])
        self.x = Xconcat

        logger.debug(self.name,self.y)
        
        if((len(self.x)/self.step) != len(self.y)):
            logger.info("===============================================")
            logger.warn("The CSV Files doesn't have the same size")
            self.dump()
            logger.info("Clip to the minimum ... :(")
            sizeMin = min(len(self.x), len(self.y))
             
            self.x = self.x[:sizeMin]
            self.y = self.y[:sizeMin]
            

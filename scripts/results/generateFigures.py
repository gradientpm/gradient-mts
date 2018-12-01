# Classical python imports
import optparse
import xml.etree.ElementTree as ET
import os
import math
import logging
import sys

# For plotting informations
import matplotlib.pyplot as plt
from matplotlib import cm

# For read images
import rgbe.io
import rgbe.utils
import rgbe.fast
try:
    import Image
    import ImageDraw
except ImportError:
    from PIL import Image
    from PIL import ImageDraw

logger = logging.getLogger(__name__)


def copyPixeltoPIL(w, h, p, im):
    clamp = lambda x : 255 if math.isnan(x) else int(x*255)
    pInt = [(clamp(c[0]),clamp(c[1]),clamp(c[2])) for c in p]
    im.putdata(pInt)

def saveNPImage(imgPath,pRef,output,scale=1.0):
     
    (width,height),pixelsHDR = rgbe.io.read(imgPath)
     
    # --- Change data to compute norm
    pRef = [lum(p) for p in pRef]
    pixelsHDR = [lum(p) for p in pixelsHDR]
     
    for i in range(len(pRef)):
        pixelsHDR[i] = 0.0 if pRef[i]+pixelsHDR[i] == 0.0 else (2.0*(pixelsHDR[i] - pRef[i]))/(pixelsHDR[i]+pRef[i])

    logger.debug("NP Image: Min",min(pixelsHDR),"Max",max(pixelsHDR))
     
    im = Image.new("RGB", (width,height))
    pixConv = lambda x,v: (int(v*255),0,0) if x < 0.0 else (0,int(v*255),0) if x > 0.0 else (0,0,0)
    buf = [pixConv(p,pow(min(abs(p)/scale,1.0),1/2.2)) for p in pixelsHDR]
    im.putdata(buf)
    im.save(output, optimize=True)
    im.close()

def saveNPImageRef(imgPath,imgRefPath,output,scale=1.0):
    """
    The difference from the previous function is here
    we read the reference again from scratch
    """

    (wRef,hRef),pRef = rgbe.io.read(imgRefPath)
    return saveNPImage(imgPath, pRef, output, scale)

def saveFig(w,h,data,output,cmapData=None, minData=None, maxData=None):
    
    # --- Save image
    fig = plt.figure(figsize=((w/100)+1, (h/100)+1), dpi=100)
    cax = plt.figimage(data, vmin=minData, vmax=maxData, cmap=cmapData)
    fig.savefig(output)#, bbox_inches=extent)
    plt.close()
    
    # Load and resave image
    im = Image.open(output)
    (widthN, heightN) = im.size
    logger.info("Detected size: ",widthN,heightN, "targeted", w, h)
    im2 = im.crop(((widthN-w), 
                  (heightN-h),
                  w,h))
    im.close()
    im2.save(output)
    
    
def lum(p):
    return 0.21268*p[0] + 0.7152*p[1] + 0.0722*p[2]

def convertImage(p, h, w, inver):
    data = []
    for x in range(h):
        tmp = []
        for y in range(w):
            p1 = lum(p[x*w+y])
            if inver > 0 and p1 > 0:
                p1 = 1.0/p1
            tmp.append(p1)
        data.append(tmp)
    return data    

def readColor(t):
    tA = [int(v) for v in t.split(",")]
    return (tA[0],tA[1],tA[2])

class MetricOp:
    def __init__(self):
        self.ref = ""
        self.img = ""
        self.exposure = 0
        self.mask = ""
        self.percent = 1.0

    def readXML(self, n):
        self.ref = n.attrib["ref"]
        self.img = n.attrib["img"]
        self.exposure = int(n.attrib["exposure"])
        if "percent" in n.attrib:
            self.percent = float(n.attrib["percent"])
    
    def show(self, wk):
        (w,h),pRef = rgbe.io.read(wk + os.path.sep + self.ref)

        maskData = []
        if(self.mask != ""):
            imMask = Image.open(wk + os.path.sep + self.mask)
            imMaskData = imMask.getdata()
            maskData = [p[0] != 255 for p in imMaskData]
        mult = float(math.pow(2, float(self.exposure)))
        images = [wk + os.path.sep + self.img]
        if(self.percent == 1.0):
            errors = rgbe.fast.rmse_all_images(w,h, images, pRef, mult, maskData)[0]
        else:
            if(self.mask != ""):
                print("ERROR: MASK and percentage is incompatible.")
                sys.exit()
            errors = rgbe.fast.rmse_all_images_percentage(w,h, self.percent, images, pRef, mult)[0]
        errorsNames = ["mse", "rmse", "mseLog", "rmseLog", "tvi", "relative", "relMSE"]
        logger.info("Metric for "+self.img+" ( with "+self.ref+") MULT: "+str(mult))
        for i in range(len(errors)):
            logger.info(errorsNames[i] + " : " + str(errors[i]))
        

class BoxOp:
    def __init__(self):
        self.x = 0
        self.y = 0
        self.sX = 0
        self.sY = 0
        self.border = 0
        self.color = (0,0,0)

    def readXML(self, n):
        self.x = int(n.attrib["x"])
        self.y = int(n.attrib["y"])
        self.sX = int(n.attrib["sX"])
        self.sY = int(n.attrib["sY"])
        self.border = int(n.attrib["border"])
        if(self.border != 0):
            self.color = readColor(n.attrib["color"])
    
    def apply(self, im, w, h):
        x0 = self.x
        y0 = self.y 
        x1 = x0 + self.sX
        y1 = y0 + self.sY
        
        im2 = im.copy()
        im2 = im2.crop((x0, y0, x1, y1))
        draw = ImageDraw.Draw(im)
        draw.rectangle((x0-self.border, 
                        y0-self.border, 
                        x1+self.border-1, 
                        y1+self.border-1), fill=self.color)
        im.paste(im2, (x0, y0))
        
        return im

class CropOp(BoxOp):
    def __init__(self):
        BoxOp.__init__(self)        

    def apply(self, im, w, h):
        x0 = self.x
        y0 = self.y 
        x1 = x0 + self.sX
        y1 = y0 + self.sY
        
        sX = self.sX + self.border*2
        sY = self.sY + self.border*2
        
        im2 = im.copy()
        im2 = im2.crop((x0-self.border, 
                        y0-self.border, 
                        x1+self.border, 
                        y1+self.border))
        
        if(self.border != 0):
            im3 = im2.copy()
            im3 = im3.crop((self.border, self.border,
                            sX - self.border, sY - self.border))
            
            draw = ImageDraw.Draw(im2)
            draw.rectangle((0, 
                            0, 
                            sX, 
                            sY), fill=self.color)
            im2.paste(im3, (self.border, self.border))
            
        return im2

class ImageOp(object):
    def __init__(self):
        self.expo = 0
        self.input = ""
        self.output = ""
        self.actions = []
        self.gamma = 2.2
        
        # Other options
        self.autoscale = False
        self.ref = ""
        
        # === internal attrib
        self.width = 0
        self.height = 0
        self.im = None
        self.pixelsHDR = []
        
        
    def readXML(self, n):
        # === Read all basic informations
        self.expo = int(n.attrib["exposure"])
        self.input = n.attrib["input"]
        self.output = n.attrib["output"]
        
        if "gamma" in n.attrib:
            self.gamma = float(n.attrib["gamma"])
        
        # === Read extra options
        if "autoscale" in n.attrib:
            self.autoscale = (n.attrib["autoscale"] == "true")
        if(self.autoscale):
            self.ref = n.attrib["ref"]
        
        # === Read all actions
        for bXML in n.iter('Box'):
            b = BoxOp()
            b.readXML(bXML)
            self.actions.append(b)
        for bXML in n.iter('Crop'):
            b = CropOp()
            b.readXML(bXML)
            self.actions.append(b)
    
    def getImg(self,wk):
        return wk + os.path.sep + self.input
    
    def loadHDR(self, wk):
        # --- Load HDR
        (self.width,self.height),self.pixelsHDR = rgbe.io.read(self.getImg(wk))
        print("Reading HDR: ", self.getImg(wk))
        
        if(self.expo != 0.0 or self.gamma != 1.0):
            logger.debug("Expo:",self.expo,"Gamma:",self.gamma)
            rgbe.utils.applyExposureGamma(self.pixelsHDR, self.expo, self.gamma)
            
        if(self.autoscale):
            (wRef,hRef),pRef = rgbe.io.read(wk + os.path.sep + self.ref)
            if(self.expo != 0.0 or self.gamma != 1.0):
                rgbe.utils.applyExposureGamma(pRef, self.expo, self.gamma)
            
            # --- Change data to compute factor
            pRefLum = [lum(p) for p in pRef]
            pixelsHDRLum = [lum(p) for p in self.pixelsHDR]
            factor = sum(pRefLum) / sum(pixelsHDRLum)
            
            # Scale the data
            for i in range(len(self.pixelsHDR)):
                r,g,b = self.pixelsHDR[i]
                self.pixelsHDR[i] = (r*factor, g*factor, b*factor)
    def loadIm(self):
        logger.debug("Simple load")
        self.im = Image.new("RGB", (self.width,self.height)) 
        copyPixeltoPIL(self.width, self.height, self.pixelsHDR, self.im)
        
    def generate(self, wk):
        self.loadHDR(wk)
        self.loadIm()
        logger.debug(self.width, self.height)
        
        for action in self.actions:
            self.im = action.apply(self.im, self.width, self.height)
        
        # --- At the end, save it
        #print("Save "+self.output)
        logger.info("Save "+self.output)
        self.im.save(wk + os.path.sep + self.output)

class ImageFalseColorOp(ImageOp):
    def __init__(self):
        ImageOp.__init__(self)
        
        # Magic numbers
        self.minV = -10000.005454
        self.maxV = 10000.005454
        self.inverse = False
        self.pMax = -1
        self.pMin = -1
        self.cmap = cm.get_cmap("viridis")
        
    def readXML(self, n):
        ImageOp.readXML(self, n)
        if "min" in n.attrib:
            self.minV = float(n.attrib["min"])
        if "max" in n.attrib:
            self.maxV = float(n.attrib["max"])
        if "pMax" in n.attrib:
            self.pMax = float(n.attrib["pMax"])
        if "pMin" in n.attrib:
            self.pMin = float(n.attrib["pMin"])
        self.inverse = (n.attrib["inverse"] == "true")
        
    def loadIm(self):
        logger.debug("Complex load")
        fig = plt.figure(figsize=((self.width/100.0), (self.height/100.0)), dpi=100)
        data = convertImage(self.pixelsHDR, self.height, 
                            self.width, self.inverse)
        
        if(self.inverse):
            pRefLum = [0.0 if lum(p)==0 else 1.0/lum(p) for p in self.pixelsHDR]
        else:
            pRefLum = [lum(p) for p in self.pixelsHDR]
        pRefLum.sort()
        
        maxData = pRefLum[-1]
        minData = pRefLum[0]

        logger.info("Find the min/max: ",str(minData),str(maxData))
  
        
        if(self.pMax != -1):
            maxData = pRefLum[int(len(pRefLum)*self.pMax)]
        if(self.pMin != -1):
            minData = pRefLum[int(len(pRefLum)*self.pMin)]
              
        if self.minV != -10000.005454:
            minData = self.minV
        if self.maxV != 10000.005454:
            maxData = self.maxV
            
        logger.info("Used min/max: ",str(minData),str(maxData))
        
        # --- Save the figure
        cax = plt.figimage(data, vmin=minData, vmax=maxData, cmap=self.cmap)
        fig.savefig(wk + os.path.sep + self.output)#, bbox_inches=extent)
        
        self.im = Image.open(wk + os.path.sep + self.output)
        (widthN, heightN) = self.im.size
        self.im = self.im.crop(((widthN-self.width), 
                                (heightN-self.height),
                                self.width,
                                self.height))
        
class ImageFalseColorNBPathsOp(ImageFalseColorOp):
    def __init__(self):
        ImageFalseColorOp.__init__(self)
    
    def loadHDR(self, wk):
        ImageFalseColorOp.loadHDR(self, wk)
        
        # --- Compute total and mean
        pLum = [p[0] for p in self.pixelsHDR]
        total = sum(pLum)
        meanPix = total / (self.height*self.width)
        logger.debug("Mean:", meanPix)
        
        # --- Normalize by mean
        pLum = [p/meanPix for p in pLum]
        
        # --- Remake 3 channel
        self.pixelsHDR = [(p,p,p) for p in pLum]
    
        
class ImageFalseColorDiffOp(ImageFalseColorOp):
    def __init__(self):
        ImageFalseColorOp.__init__(self)
        self.ref = ""
        
    def readXML(self, n):
        ImageFalseColorOp.readXML(self, n)
        self.ref = n.attrib["ref"]
    
    def loadHDR(self, wk):
        ImageFalseColorOp.loadHDR(self, wk)
        
        # --- Load reference
        (wRef,hRef),pRef = rgbe.io.read(wk + os.path.sep + self.ref)
        if(self.expo != 0.0 or self.gamma != 1.0):
            rgbe.utils.applyExposureGamma(pRef, self.expo, self.gamma)
        
        # --- Change data to compute norm
        pRef = [lum(p) for p in pRef]
        self.pixelsHDR = [lum(p) for p in self.pixelsHDR]
        
        # --- Compute relative error
        for i in range(len(pRef)):
            self.pixelsHDR[i] = 0.0 if pRef[i] == 0.0 else abs(pRef[i] - self.pixelsHDR[i])/pRef[i]
            
        #maxError = max(self.pixelsHDR)
        #self.pixelsHDR = [maxError if p == 0 else p for p in self.pixelsHDR]
        
        # --- Remake 3 channel
        self.pixelsHDR = [(p,p,p) for p in self.pixelsHDR]
        

class ImagePNDiffOp(ImageFalseColorOp):
    def __init__(self):
        ImageFalseColorOp.__init__(self)
        self.ref = ""
        self.scale = 1.0
        
    def readXML(self, n):
        ImageFalseColorOp.readXML(self, n)
        self.ref = n.attrib["ref"]
        self.cmap = None
    
    def loadHDR(self, wk):
        ImageFalseColorOp.loadHDR(self, wk)
        
        (wRef,hRef),pRef = rgbe.io.read(wk + os.path.sep + self.ref)
        if(self.expo != 0.0 or self.gamma != 1.0):
            rgbe.utils.applyExposureGamma(pRef, self.expo, self.gamma)
        
        # --- Change data to compute norm
        pRef = [lum(p) for p in pRef]
        self.pixelsHDR = [lum(p) for p in self.pixelsHDR]

                
        for i in range(len(pRef)):
            self.pixelsHDR[i] = 0.0 if pRef[i]+self.pixelsHDR[i] == 0.0 else (2.0*(self.pixelsHDR[i] - pRef[i]))/(self.pixelsHDR[i]+pRef[i])

        logger.info("Min",min(self.pixelsHDR),"Max",max(self.pixelsHDR))
        self.maxV = 1.0
        self.minV = -1.0
        
        self.pixelsHDR = [(p,p,p) for p in self.pixelsHDR]
    
    def loadIm(self):
        logger.debug("Hard convert")
        self.im = Image.new("RGB", (self.width,self.height))
        pixConv = lambda x,v: (int(v*255),0,0) if x < 0.0 else (0,int(v*255),0) if x > 0.0 else (0,0,0)
        buf = [pixConv(p[0],pow(min(abs(p[0])/self.scale,1.0),1/2.2)) for p in self.pixelsHDR]
        self.im.putdata(buf)
   
def readXMLComparisons(file):
    tree = ET.parse(file)
    root = tree.getroot()

    images = []

    for imNode in root.iter('Image'):
        im = ImageOp()
        im.readXML(imNode)
        images.append(im)
    for imNode in root.iter('ImageFalseColor'):
        im = ImageFalseColorOp()
        im.readXML(imNode)
        images.append(im)
    for imNode in root.iter('ImageFalseColorDiff'):
        im = ImageFalseColorDiffOp()
        im.readXML(imNode)
        images.append(im)
    for imNode in root.iter('ImageFalseColorNbPaths'):
        im = ImageFalseColorNBPathsOp()
        im.readXML(imNode)
        images.append(im)
    for imNode in root.iter('ImageNP'):
        im = ImagePNDiffOp()
        im.readXML(imNode)
        images.append(im)
        
    displays = []
    for disNode in root.iter('DisplayMetric'):
        dis = MetricOp()
        dis.readXML(disNode)
        displays.append(dis)
        
    return images,displays

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)

    # --- Read all params
    parser = optparse.OptionParser()
    parser.add_option('-i','--input', help='input')
    (opts, args) = parser.parse_args()
    
    inXML = opts.input
    wk = os.path.dirname(inXML)
    
    print(wk)
    images,displays = readXMLComparisons(inXML)
    
    while len(images) != 0:
        images[-1].generate(wk)
        images.pop() # Remove last one
    
    for display in displays:
        display.show(wk)
    

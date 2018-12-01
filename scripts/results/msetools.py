import optparse

# Import to read rgbe images
import rgbe.io
import rgbe.fast
import rgbe.utils
try:
    import Image
except ImportError:
    from PIL import Image

# CSV based name files
CSVNames = ["_mse.csv", "_rmse.csv", "_mseLog.csv", "_rmseLog.csv", "_tvi.csv", "_relative.csv", "_relMSE.csv"]

# === Constants
def diffImageRef(currImage, pRef, mult, metricType, maskImage = "", outputImage = ""):
    (w,h),pCur = rgbe.io.read(currImage)
    
    # Load mask if necessary
    mask = None
    if(maskImage != ""):
        imgMask = Image.open(maskImage)
        mask = list(imgMask.getdata())
        
    rmsevalue, imgDiffData = rgbe.fast.rmse(w, h, pCur, pRef, mult, metricType, mask)
    
    # TODO: Fix img
    if(outputImage != ""):
        print("[INFO] Write image: ", outputImage)
        imgDiff = Image.new("RGB", (w, h))
        rgbe.utils.copyPixeltoPIL(w, h, imgDiffData, imgDiff);
        imgDiff.save(outputImage)
    
    return rmsevalue

def diffImage(currImage, finalImage, mult, metricType=0, maskImage = "", outputImage = "", all=False):
    (w,h),pRef = rgbe.io.read(finalImage)
    if(all):
        #FIXME: Make the error management
        metrics = rgbe.fast.rmse_all_images(w,h, [currImage], pRef, mult, None)
        print(metrics)
        return metrics[0]
    else:
        return diffImageRef(currImage, pRef, mult, metricType, maskImage, outputImage)

def computeMSEAll(filename, nbImages, steps, finalImage, percentage = 1.0, outputs=[], mult=1,  maskImg = ''):
    # === Open all files
    # These files will be used to log 
    # All the metric values
    files = []
    for output in outputs:
        f = open(output,'w')
        files.append(f)


    # Read reference image
    (w,h),pRef = rgbe.io.read(finalImage)

    # List all HDR images
    imagesHDR = []
    i = steps
    while i <= nbImages:
        imagesHDR.append(filename + str(i) + '.hdr')
        i += steps

    # Debug
    #print(imagesHDR)
    #print("\n")

    # Launch the computation
    if(percentage == 1.0):
        maskData = []
        if(maskImg != ""):
            imMask = Image.open(maskImg)
            imMaskData = imMask.getdata()
            maskData = [p[0] != 255 for p in imMaskData]
            print("Read mask", maskImg)
        metrics = rgbe.fast.rmse_all_images(w,h, imagesHDR, pRef, mult, maskData)
    else:
        metrics = rgbe.fast.rmse_all_images_percentage(w,h, percentage, imagesHDR, pRef, mult)

    if(len(files) != 0):
        # Write down all files
        for k in range(len(metrics)):
            for j in range(len(metrics[0])):
                files[j].write(str(metrics[k][j]) + ',\n')
    return metrics
        
if __name__ == "__main__":
    # --- Read all params
    parser = optparse.OptionParser()
    parser.add_option('-f','--file', help='Filename')
    parser.add_option('-n','--nbImages', help='Number of images computed. If the last file name ends with 1203, this argument should be 1203.', default=1)
    parser.add_option('-s','--step', default=1,help='Images step. If only 1 image is written out of 10 computed, this argument should be 10.')
    parser.add_option('-r','--reference',help='Reference image')
    parser.add_option('-o','--output',help='csv output file', default="mse.csv")
    parser.add_option('-m','--mult', help='Multiplication factor for the difference image', default=1.0)
    parser.add_option('-M','--mask',help='Mask image (Optional)', default="")
    parser.add_option('-1','--onlyone',help='To compute the MSE for a given image', action="store_true", default=False)
    parser.add_option('-p','--metric',help='metric choose (0 = MSE, 1 = RMSE, 2 = MSE_log, 3 = RMSE_log, 4 = TVI)', default=1)
    parser.add_option('-A','--all',help='To compute the MSE for a given image', action="store_true", default=False)
    
    (opts, args) = parser.parse_args()
    
    filename = opts.file
    nbImages = int(opts.nbImages)
    steps = int(opts.step)
    finalImage = opts.reference
    output = opts.output
    mult = float(opts.mult)
    maskImage = opts.mask
    metricType = int(opts.metric)
    
    #TODO: Need to reimplement it
    #Because we change the function name/signatures

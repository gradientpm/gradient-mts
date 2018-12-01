import glob
import rgbe.io
import optparse

def crop(inFile, outFile, cX, cY, cXSize, cYSize):
    (width,height),pixelsHDR = rgbe.io.read(inFile)
    pixelsHDRCrop = [(0,0,0)] * cXSize * cYSize
    for i in range(len(pixelsHDR)):
        x = i % width
        y = (i - x) // width
        if(x >= cX and x < cX+cXSize and y >= cY and y < cY+cYSize):
            k = (x - cX) + cXSize * (y - cY)
            pixelsHDRCrop[k] = pixelsHDR[i]
    rgbe.io.write_hdr(outFile, cXSize, cYSize, pixelsHDRCrop)

if __name__ == "__main__":
    # --- Read all params
    parser = optparse.OptionParser()
    parser.add_option('-i','--input', help='input')
    parser.add_option('-o','--output', help='output')
    parser.add_option('-x','--cropx', help='crop X coordinate')
    parser.add_option('-y','--cropy', help='crop Y coordinate')
    parser.add_option('-l','--cropxsize', help='crop X size')
    parser.add_option('-k','--cropysize', help='crop Y size')
    parser.add_option('-B', '--batch', help="the input file is a list of files", default=False, action="store_true")
    (opts, args) = parser.parse_args()

    cX = int(opts.cropx)
    cY = int(opts.cropy)
    cXSize = int(opts.cropxsize)
    cYSize = int(opts.cropysize)

    if(opts.batch):
        import glob
        hdrFiles = glob.glob(DIR+"*.hdr")
        idFile = 0
        for f in hdrFiles:
            crop.crop(f,f, cX, cY, cXSize, cYSize)
            print(idFile, "/", len(hdrFiles), " cropped...")
            idFile += 1
    else:
        print("Crop: ", opts.input, " => ", opts.output)
        crop(opts.input, opts.output, cX, cY, cXSize, cYSize)


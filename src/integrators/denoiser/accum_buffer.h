#pragma once

#include <mitsuba/mitsuba.h>
#include <mitsuba/core/bitmap.h>

MTS_NAMESPACE_BEGIN

/**
 * A simple wrapper to track accumulated radiance 
 * and its population and sample variance
 */
class AccumBuffer : public SerializableObject {
public:
    AccumBuffer(int numBuffers, const Vector2i& s, bool trackVariance, bool twoHalfBuffers) :
            size(s), trackVariance(trackVariance), twoHalfBuffers(twoHalfBuffers) {

        // Value buffers 
        accumValue.resize(numBuffers);
        for (int i = 0; i < numBuffers; ++i) {
            accumValue[i] = new Bitmap(Bitmap::ESpectrum, Bitmap::EFloat, size);
            accumValue[i]->clear();
        }

        if (twoHalfBuffers) {
            accumValueA.resize(numBuffers);
            accumValueB.resize(numBuffers);
            for (int i = 0; i < numBuffers; ++i) {
                accumValueA[i] = new Bitmap(Bitmap::ESpectrum, Bitmap::EFloat, size);
                accumValueA[i]->clear();
                accumValueB[i] = new Bitmap(Bitmap::ESpectrum, Bitmap::EFloat, size);
                accumValueB[i]->clear();
            }
        }

        sVariance.resize(numBuffers);       // Just invalid references
        if (trackVariance) {
            for (int i = 0; i < numBuffers; ++i) {
                sVariance[i] = new Bitmap(Bitmap::ESpectrum, Bitmap::EFloat, size);
                sVariance[i]->clear();
            }
        }

        // Count buffers
        accumCountVec.resize(numBuffers);
        for (int i = 0; i < numBuffers; ++i) {
            accumCountVec[i] = new Bitmap(Bitmap::ELuminance, Bitmap::EFloat, size);
            accumCountVec[i]->clear();
        }
    }

    /**
     * A convenient constructor to store multiple bitmaps that track accumulation and variance
     */
    AccumBuffer(std::vector<ref<Bitmap> > vecAccum, 
                std::vector<ref<Bitmap> > vecVariance, 
                std::vector<ref<Bitmap> > vecAccumA,
                std::vector<ref<Bitmap> > vecAccumB,
                int count) :
            trackVariance(true), twoHalfBuffers(true), size(vecAccum[0]->getSize())
    {   
        int numBuffers = vecAccum.size();

        // Value buffers 
        accumValue.resize(numBuffers);
        for (int i = 0; i < numBuffers; ++i) {
            accumValue[i] = new Bitmap(Bitmap::ESpectrum, Bitmap::EFloat, size);
            accumValue[i]->copyFrom(vecAccum[i]);
        }

        if (twoHalfBuffers) {
            accumValueA.resize(numBuffers);
            accumValueB.resize(numBuffers);
            for (int i = 0; i < numBuffers; ++i) {
                accumValueA[i] = new Bitmap(Bitmap::ESpectrum, Bitmap::EFloat, size);
                accumValueA[i]->copyFrom(vecAccumA[i]);
                accumValueB[i] = new Bitmap(Bitmap::ESpectrum, Bitmap::EFloat, size);
                accumValueB[i]->copyFrom(vecAccumB[i]);
            }
        }

        sVariance.resize(numBuffers);       // Just invalid references
        if (trackVariance) {
            for (int i = 0; i < numBuffers; ++i) {
                sVariance[i] = new Bitmap(Bitmap::ESpectrum, Bitmap::EFloat, size);
                sVariance[i]->copyFrom(vecVariance[i]);
            }
        }

        // Count buffers
        accumCountVec.resize(numBuffers);
        for (int i = 0; i < numBuffers; ++i) {
            accumCountVec[i] = new Bitmap(Bitmap::ELuminance, Bitmap::EFloat, size);
            accumCountVec[i]->clear();
            Float *data = (Float *)accumCountVec[i]->getFloatData();
            for (int k = 0; k < size.x * size.y; ++k) data[k] = count;
        }
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
    };
    
    const Vector2i getSize() const { return size; }

    void clear() {
        int numBuffers = accumValue.size();
        for (int i = 0; i < numBuffers; ++i) {            
            accumValue[i]->clear();
        }

        if (twoHalfBuffers) {            
            for (int i = 0; i < numBuffers; ++i) {                
                accumValueA[i]->clear();                
                accumValueB[i]->clear();
            }
        }
        
        if (trackVariance) {
            for (int i = 0; i < numBuffers; ++i) {
                sVariance[i]->clear();
            }
        }

        for (int i = 0; i < numBuffers; ++i) {
            accumCountVec[i]->clear();
        }
    }

    ref <Bitmap> getSampleVariance(int buffer) const {
        return sVariance[buffer];
    }

    ref <Bitmap> getBuffer(int buffer) const {
        return accumValue[buffer];
    }

    ref <Bitmap> getBufferA(int buffer) const {
        return accumValueA[buffer];
    }

    ref <Bitmap> getBufferB(int buffer) const {
        return accumValueB[buffer];
    }

    ref <Bitmap> getCountBitmap(int buffer) const {
        return accumCountVec[buffer];
    }

    Spectrum getValue(int idx, int buffer) const {
        Spectrum *data = (Spectrum *)accumValue[buffer]->getFloatData();
        return data[idx];
    }

    void addSample(int x, int y, const Spectrum& value, int buffer = 0) {        
        if(x >= size.x || y >= size.y) {
            return;
        }

        if (twoHalfBuffers)
            addSampleTungsten(x, y, value, buffer);

        addSampleOriginal(x, y, value, buffer);
    }

    int getBufferCount() const {
        return accumValue.size();
    }

private:
    void addSampleTungsten(int x, int y, const Spectrum& value, int buffer = 0) {
        int idx = y * size.x + x;

        const Float* sampleCount = (const Float *) accumCountVec[buffer]->getFloatData();
        Spectrum* bufferA = (Spectrum *) accumValueA[buffer]->getFloatData();
        Spectrum* bufferB = (Spectrum *) accumValueB[buffer]->getFloatData();
        
        unsigned int sampleIdx = sampleCount[idx];

        // Accumulate alternatively to each buffer depending on the index
        Spectrum *feature = (sampleIdx & 1) ? bufferB : bufferA;
        unsigned int perBufferSampleCount = sampleIdx / 2 + 1;
        feature[idx] += (value - feature[idx]) / perBufferSampleCount;        
    }

    void addSampleOriginal(int x, int y, const Spectrum& value, int buffer = 0){

        // Compute online mean and variance
        Vector2i size = accumValue[buffer]->getSize();
        int pixelID = y * size.x + x;

        Spectrum *meanPtr = (Spectrum *)accumValue[buffer]->getFloatData();        
        Float *countPtr = (Float *)accumCountVec[buffer]->getFloatData();

        Spectrum sample = value;
        Spectrum mean = meanPtr[pixelID];
        int n = countPtr[pixelID] + 1;
        
        const Spectrum delta = sample - mean;       // sample - old mean 
        mean += delta / n;                          // estimate new mean

        if (trackVariance) {
            Spectrum *variancePtr = (Spectrum *)sVariance[buffer]->getFloatData();
            Spectrum variance(0.);
            if (n > 1) {
                // Numerically robust online variance estimation using an
                // algorithm proposed by Donald Knuth (TAOCP vol.2, 3rd ed., p.232)            
                Spectrum diffSqr = variancePtr[pixelID] * (n - 2) * (n - 1);
                diffSqr += delta * (sample - mean);         // (sample - new mean) * (sample - old mean)            

                                                            // Since we do not know the population mean, 
                                                            // we use sample mean in our variance calculation, and this causes bias. 
                                                            // We can use Bessel's correction (divide by n - 1 instead of n) to get an unbiased estimate of the population variance).
                Spectrum populationVariance = diffSqr / (n - 1);

                variance = populationVariance / n;
            }
            variancePtr[pixelID] = variance;
        }

        // Save 
        meanPtr[pixelID] = mean;        
        countPtr[pixelID] = n;
    }

private:
    const Vector2i size;
    bool trackVariance;
    std::vector< ref<Bitmap> > accumValue;      // Multiple accumulation buffers
    std::vector< ref<Bitmap> > sVariance;       // Unbiased sample variance (population variance / n)
    std::vector< ref<Bitmap> > accumCountVec;   // Number of non-empty pass per pixel. In each sample, the number of samples is >= 0.

    bool twoHalfBuffers;
    std::vector< ref<Bitmap> > accumValueA;     // Half accumulation buffers 
    std::vector< ref<Bitmap> > accumValueB;     // Half accumulation buffers
};

MTS_NAMESPACE_END

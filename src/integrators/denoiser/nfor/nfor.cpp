#include "core/thread/ThreadUtils.hpp"
#include "core/Timer.hpp"

#include "Regression.hpp"
#include "NlMeans.hpp"
#include "Pixmap.hpp"

#include "core/math/MathUtil.hpp"

#include "nfor.h"

#include <vector>
#include <cstdlib>
#include <cstring>
#include <cmath>

using namespace Tungsten;

static const int OPT_VERSION  = 0;
static const int OPT_HELP     = 1;

template<typename Texel>
struct RenderBuffer
{
    std::unique_ptr<Pixmap<Texel>> buffer;
    std::unique_ptr<Pixmap<Texel>> bufferA;
    std::unique_ptr<Pixmap<Texel>> bufferB;
    std::unique_ptr<Pixmap<Texel>> bufferVariance;
};
typedef RenderBuffer<float> RenderBufferF;
typedef RenderBuffer<Vec3f> RenderBuffer3f;

Pixmap3f nforDenoiser(RenderBuffer3f image, std::vector<RenderBufferF> features)
{
    int w = image.buffer->w(), h = image.buffer->h();

    // Feature cross-prefiltering (section 5.1)
    printTimestampedLog("Prefiltering features...");
    std::vector<PixmapF> filteredFeaturesA(features.size());
    std::vector<PixmapF> filteredFeaturesB(features.size());
    SimdNlMeans featureFilter;
    for (size_t i = 0; i < features.size(); ++i) {
        featureFilter.addBuffer(filteredFeaturesA[i], *features[i].bufferA, *features[i].bufferB, *features[i].bufferVariance);
        featureFilter.addBuffer(filteredFeaturesB[i], *features[i].bufferB, *features[i].bufferA, *features[i].bufferVariance);
    }
    featureFilter.denoise(3, 5, 0.5f, 2.0f);
    features.clear();
    printTimestampedLog("Prefiltering done");

    // Main regression (section 5.2)
    std::vector<Pixmap3f> filteredColorsA;
    std::vector<Pixmap3f> filteredColorsB;
    std::vector<Pixmap3f> mses;
    for (float k : {0.5f, 1.0f}) {
        printTimestampedLog(tfm::format("Beginning regression pass %d/2", mses.size() + 1));
        // Regression pass
        printTimestampedLog("Denosing half buffer A...");
        Pixmap3f filteredColorA = collaborativeRegression(*image.bufferA, *image.bufferB, filteredFeaturesB, *image.bufferVariance, 3, 9, k);
        printTimestampedLog("Denosing half buffer B...");
        Pixmap3f filteredColorB = collaborativeRegression(*image.bufferB, *image.bufferA, filteredFeaturesA, *image.bufferVariance, 3, 9, k);

        // MSE estimation (section 5.3)
        printTimestampedLog("Estimating MSE...");
        Pixmap3f noisyMse(w, h);
        for (int i = 0; i < w*h; ++i) {
            Vec3f mseA = sqr((*image.bufferB)[i] - filteredColorA[i]) - 2.0f*(*image.bufferVariance)[i];
            Vec3f mseB = sqr((*image.bufferA)[i] - filteredColorB[i]) - 2.0f*(*image.bufferVariance)[i];
            Vec3f residualColorVariance = sqr(filteredColorB[i] - filteredColorA[i])*0.25f;

            noisyMse[i] = (mseA + mseB)*0.5f - residualColorVariance;
        }
        filteredColorsA.emplace_back(std::move(filteredColorA));
        filteredColorsB.emplace_back(std::move(filteredColorB));

        // MSE filtering
        mses.emplace_back(nlMeans(noisyMse, *image.buffer, *image.bufferVariance, 1, 9, 1.0f, 1.0f, true));
    }
    printTimestampedLog("Regression pass done");

    // Bandwidth selection (section 5.3)
    // Generate selection map
    printTimestampedLog("Generating selection maps...");
    Pixmap3f noisySelection(w, h);
    for (int i = 0; i < w*h; ++i)
        for (int j = 0; j < 3; ++j)
            noisySelection[i][j] = mses[0][i][j] < mses[1][i][j] ? 0.0f : 1.0f;
    mses.clear();
    // Filter selection map
    Pixmap3f selection = nlMeans(noisySelection, *image.buffer, *image.bufferVariance, 1, 9, 1.0f, 1.0f, true);

    // Apply selection map
    Pixmap3f resultA(w, h);
    Pixmap3f resultB(w, h);
    for (int i = 0; i < w*h; ++i) {
        resultA[i] += lerp(filteredColorsA[0][i], filteredColorsA[1][i], selection[i]);
        resultB[i] += lerp(filteredColorsB[0][i], filteredColorsB[1][i], selection[i]);
    }
    selection.reset();
    filteredColorsA.clear();
    filteredColorsB.clear();

    // Second filter pass (section 5.4)
    printTimestampedLog("Beginning second filter pass");
    printTimestampedLog("Denoising final features...");
    std::vector<PixmapF> finalFeatures;
    for (size_t i = 0; i < filteredFeaturesA.size(); ++i) {
        PixmapF combinedFeature(w, h);
        PixmapF combinedFeatureVar(w, h);

        for (int j = 0; j < w*h; ++j) {
            combinedFeature   [j] =    (filteredFeaturesA[i][j] + filteredFeaturesB[i][j])*0.5f;
            combinedFeatureVar[j] = sqr(filteredFeaturesB[i][j] - filteredFeaturesA[i][j])*0.25f;
        }
        filteredFeaturesA[i].reset();
        filteredFeaturesB[i].reset();

        finalFeatures.emplace_back(nlMeans(combinedFeature, combinedFeature, combinedFeatureVar, 3, 2, 0.5f));
    }

    Pixmap3f combinedResult(w, h);
    Pixmap3f combinedResultVar(w, h);
    for (int j = 0; j < w*h; ++j) {
        combinedResult   [j] =    (resultA[j] + resultB[j])*0.5f;
        combinedResultVar[j] = sqr(resultB[j] - resultA[j])*0.25f;
    }
    printTimestampedLog("Performing final regression...");
    return collaborativeRegression(combinedResult, combinedResult, finalFeatures, combinedResultVar, 3, 9, 1.0f);
}

/**
 * A custom NFOR that properly use reconstruction image to guide the filtering of primal image. 
 */ 
Pixmap3f nforDenoiserReconstruction(RenderBuffer3f image, std::vector<RenderBufferF> features)
{
    int w = image.buffer->w(), h = image.buffer->h();

    // Feature cross-prefiltering (section 5.1)
    
    printTimestampedLog("Prefiltering features...");
    std::vector<PixmapF> filteredFeaturesA(features.size());
    std::vector<PixmapF> filteredFeaturesB(features.size());
    SimdNlMeans featureFilter;
    for (size_t i = 0; i < features.size(); ++i) {
        featureFilter.addBuffer(filteredFeaturesA[i], *features[i].bufferA, *features[i].bufferB, *features[i].bufferVariance);
        featureFilter.addBuffer(filteredFeaturesB[i], *features[i].bufferB, *features[i].bufferA, *features[i].bufferVariance);
    }
    featureFilter.denoise(3, 5, 0.5f, 2.0f);
    features.clear();
    printTimestampedLog("Prefiltering done");
        
    // No feature prefiltering 
    /*
    std::vector<PixmapF> filteredFeaturesA(features.size());
    std::vector<PixmapF> filteredFeaturesB(features.size());
    for (size_t i = 0; i < features.size(); ++i) {
        filteredFeaturesA[i] = PixmapF(w, h, *features[i].bufferA);
        filteredFeaturesB[i] = PixmapF(w, h, *features[i].bufferB);
    }*/

    // Main regression (section 5.2)
    std::vector<Pixmap3f> filteredColorsA;
    std::vector<Pixmap3f> filteredColorsB;
    std::vector<Pixmap3f> mses;

    /*
    // As our feature is a color map, use it as guide instead of the noisy buffer
    auto featuresColorA = std::unique_ptr<Pixmap3f>(new Pixmap3f(w, h));
    auto featuresColorB = std::unique_ptr<Pixmap3f>(new Pixmap3f(w, h));
    auto featuresColorVariance = std::unique_ptr<Pixmap3f>(new Pixmap3f(w, h));
    // Assume features only has 3 channels
    for (int j = 0; j < w*h; ++j) {
        (*featuresColorA)[j] = Vec3f((*features[0].bufferA)[j], (*features[1].bufferA)[j], (*features[2].bufferA)[j]);
        (*featuresColorB)[j] = Vec3f((*features[0].bufferB)[j], (*features[1].bufferB)[j], (*features[2].bufferB)[j]);
        (*featuresColorVariance)[j] = Vec3f((*features[0].bufferVariance)[j], (*features[1].bufferVariance)[j], (*features[2].bufferVariance)[j]);
    }*/

    for (float k : {0.5f, 1.0f}) {
        printTimestampedLog(tfm::format("Beginning regression pass %d/2", mses.size() + 1));
        // Regression pass        
        printTimestampedLog("Denosing half buffer A...");
        Pixmap3f filteredColorA = collaborativeRegression(*image.bufferA, *image.bufferB, filteredFeaturesB, *image.bufferVariance, 3, 9, k);
        printTimestampedLog("Denosing half buffer B...");
        Pixmap3f filteredColorB = collaborativeRegression(*image.bufferB, *image.bufferA, filteredFeaturesA, *image.bufferVariance, 3, 9, k);
        
        /*
        printTimestampedLog("Denosing half buffer A...");
        Pixmap3f filteredColorA = collaborativeRegression(*image.bufferA, *featuresColorB, filteredFeaturesB, *featuresColorVariance, 3, 9, k);
        printTimestampedLog("Denosing half buffer B...");
        Pixmap3f filteredColorB = collaborativeRegression(*image.bufferB, *featuresColorA, filteredFeaturesA, *featuresColorVariance, 3, 9, k);
        */

        // MSE estimation (section 5.3)
        printTimestampedLog("Estimating MSE...");
        Pixmap3f noisyMse(w, h);
        for (int i = 0; i < w*h; ++i) {
            Vec3f mseA = sqr((*image.bufferB)[i] - filteredColorA[i]) - 2.0f*(*image.bufferVariance)[i];
            Vec3f mseB = sqr((*image.bufferA)[i] - filteredColorB[i]) - 2.0f*(*image.bufferVariance)[i];
            Vec3f residualColorVariance = sqr(filteredColorB[i] - filteredColorA[i])*0.25f;

            noisyMse[i] = (mseA + mseB)*0.5f - residualColorVariance;
        }
        filteredColorsA.emplace_back(std::move(filteredColorA));
        filteredColorsB.emplace_back(std::move(filteredColorB));

        // MSE filtering
        mses.emplace_back(nlMeans(noisyMse, *image.buffer, *image.bufferVariance, 1, 9, 1.0f, 1.0f, true));
    }
    printTimestampedLog("Regression pass done");

    // Bandwidth selection (section 5.3)
    // Generate selection map
    printTimestampedLog("Generating selection maps...");
    Pixmap3f noisySelection(w, h);
    for (int i = 0; i < w*h; ++i)
        for (int j = 0; j < 3; ++j)
            noisySelection[i][j] = mses[0][i][j] < mses[1][i][j] ? 0.0f : 1.0f;
    mses.clear();
    // Filter selection map
    Pixmap3f selection = nlMeans(noisySelection, *image.buffer, *image.bufferVariance, 1, 9, 1.0f, 1.0f, true);

    // Apply selection map
    Pixmap3f resultA(w, h);
    Pixmap3f resultB(w, h);
    for (int i = 0; i < w*h; ++i) {
        resultA[i] += lerp(filteredColorsA[0][i], filteredColorsA[1][i], selection[i]);
        resultB[i] += lerp(filteredColorsB[0][i], filteredColorsB[1][i], selection[i]);
    }
    selection.reset();
    filteredColorsA.clear();
    filteredColorsB.clear();

    // Second filter pass (section 5.4)
    printTimestampedLog("Beginning second filter pass");

    /*
    printTimestampedLog("Denoising final features...");
    std::vector<PixmapF> finalFeatures;
    for (size_t i = 0; i < filteredFeaturesA.size(); ++i) {
        PixmapF combinedFeature(w, h);
        PixmapF combinedFeatureVar(w, h);

        for (int j = 0; j < w*h; ++j) {
            combinedFeature   [j] =    (filteredFeaturesA[i][j] + filteredFeaturesB[i][j])*0.5f;
            combinedFeatureVar[j] = sqr(filteredFeaturesB[i][j] - filteredFeaturesA[i][j])*0.25f;
        }
        filteredFeaturesA[i].reset();
        filteredFeaturesB[i].reset();

        finalFeatures.emplace_back(nlMeans(combinedFeature, combinedFeature, combinedFeatureVar, 3, 2, 0.5f));
    }*/

    // Final features from the given buffer because reconstruction is not linear
    std::vector<PixmapF> finalFeatures(features.size());    
    for (size_t i = 0; i < features.size(); ++i) {
        // no filter 
        // finalFeatures[i] = PixmapF(w, h, *features[i].buffer);

        // with filter 
        PixmapF combinedFeature(w, h);
        PixmapF combinedFeatureVar(w, h); 
        for (int j = 0; j < w*h; ++j) {
            auto mean = (*features[i].buffer)[j];
            combinedFeature   [j] =    mean;            
            combinedFeatureVar[j] = 0.5 * ((filteredFeaturesA[i][j] - mean) * (filteredFeaturesA[i][j] - mean) +
                                           (filteredFeaturesB[i][j] - mean) * (filteredFeaturesB[i][j] - mean));
        }
        finalFeatures.emplace_back(nlMeans(combinedFeature, combinedFeature, combinedFeatureVar, 3, 2, 0.5f));
    }

    

    Pixmap3f combinedResult(w, h);
    Pixmap3f combinedResultVar(w, h);
    for (int j = 0; j < w*h; ++j) {
        combinedResult   [j] =    (resultA[j] + resultB[j])*0.5f;
        combinedResultVar[j] = sqr(resultB[j] - resultA[j])*0.25f;
    }
    printTimestampedLog("Performing final regression...");
    return collaborativeRegression(combinedResult, combinedResult, finalFeatures, combinedResultVar, 3, 9, 1.0f);
}


// Extracts a single channel of an RGB image into a separate pixmap
std::unique_ptr<PixmapF> slicePixmap(const Pixmap3f &src, int channel)
{
    int w = src.w(), h = src.h();

    auto result = std::unique_ptr<PixmapF>(new PixmapF(w, h));
    for (int j = 0; j < w*h; ++j)
        (*result)[j] = src[j][channel];

    return std::move(result);
}

/** 
 * Denoise with auxiliary feature buffer such as albedo, normal. 
 */
void Nfor::denoise(AccumBuffer* throughput, Bitmap *reconsBitmap, AuxiliaryBuffer* auxBuffer){
    
    ThreadUtils::startThreads(std::max(ThreadUtils::idealThreadCount() - 1, 1u));

    printTimestampedLog(tfm::format("Filtering with feature buffers..."));

    RenderBuffer3f image;
    std::vector<RenderBufferF> features;

    Spectrum *throughputPtr  = (Spectrum *)throughput->getBuffer(0)->getFloatData();
    Spectrum *variancePtr    = (Spectrum *)throughput->getSampleVariance(0)->getFloatData();
    Spectrum *throughputAPtr = (Spectrum *)throughput->getBufferA(0)->getFloatData();
    Spectrum *throughputBPtr = (Spectrum *)throughput->getBufferB(0)->getFloatData();

    Vector2i size = throughput->getSize();

    image.buffer = std::unique_ptr<Pixmap3f>(new Pixmap3f(size.x, size.y));
    image.bufferVariance = std::unique_ptr<Pixmap3f>(new Pixmap3f(size.x, size.y));
    image.bufferA = std::unique_ptr<Pixmap3f>(new Pixmap3f(size.x, size.y));
    image.bufferB = std::unique_ptr<Pixmap3f>(new Pixmap3f(size.x, size.y));

    for (int j = 0; j < size.x * size.y; ++j) {        
        (*image.buffer)[j] = Vec3f(
            (float)throughputPtr[j][0], 
            (float)throughputPtr[j][1], 
            (float)throughputPtr[j][2]);

        (*image.bufferVariance)[j] = Vec3f(
            (float)variancePtr[j][0],
            (float)variancePtr[j][1],
            (float)variancePtr[j][2]);

        (*image.bufferA)[j] = Vec3f(
                (float)throughputAPtr[j][0],
                (float)throughputAPtr[j][1],
                (float)throughputAPtr[j][2]);

        (*image.bufferB)[j] = Vec3f(
                (float)throughputBPtr[j][0],
                (float)throughputBPtr[j][1],
                (float)throughputBPtr[j][2]);
    }

    // Assume each auxiliary buffer has 3 channels
    std::vector<EAuxBuffer> auxBufferIds = { EAuxBuffer::AlbedoBuffer, EAuxBuffer::PositionBuffer, EAuxBuffer::ShadingNormalW  };

    features.resize(auxBufferIds.size() * 3);
    for (size_t i = 0; i < features.size(); ++i) {
        //features[i].buffer = std::unique_ptr<PixmapF>(new PixmapF(size.x, size.y));
        features[i].bufferA = std::unique_ptr<PixmapF>(new PixmapF(size.x, size.y));
        features[i].bufferB = std::unique_ptr<PixmapF>(new PixmapF(size.x, size.y));
        features[i].bufferVariance = std::unique_ptr<PixmapF>(new PixmapF(size.x, size.y));
    }

    for (int i = 0; i < auxBufferIds.size(); i++) {
        int k = auxBufferIds[i];
        for (int j = 0; j < size.x * size.y; ++j) {

            //Vector &buffer = auxBuffer->getBuffer(k).get()->getV3Data()[j];
            //(*features[i*3+0].buffer)[j] = buffer[0];
            //(*features[i*3+1].buffer)[j] = buffer[1];
            //(*features[i*3+2].buffer)[j] = buffer[2];

            Vector &bufferA = auxBuffer->getBufferA(k).get()->getV3Data()[j];
            (*features[i*3+0].bufferA)[j] = bufferA.x;
            (*features[i*3+1].bufferA)[j] = bufferA.y;
            (*features[i*3+2].bufferA)[j] = bufferA.z;

            Vector &bufferB = auxBuffer->getBufferB(k).get()->getV3Data()[j];
            (*features[i*3+0].bufferB)[j] = bufferB.x;
            (*features[i*3+1].bufferB)[j] = bufferB.y;
            (*features[i*3+2].bufferB)[j] = bufferB.z;

            Vector &variance = auxBuffer->getVariance(k).get()->getV3Data()[j];
            (*features[i*3+0].bufferVariance)[j] = variance.x;
            (*features[i*3+1].bufferVariance)[j] = variance.y;
            (*features[i*3+2].bufferVariance)[j] = variance.z;
        }
    }

    /*
    std::cout << "Prepare to call save" << std::endl;
    (*image.buffer).save("buffer.exr", true);
    image.bufferA->save("bufferA.exr", true);
    image.bufferB->save("bufferB.exr", true);
    image.bufferVariance->save("bufferVariance.exr", true);

    for (int i = 0; i < features.size(); ++i) {
        {
            std::stringstream ss; 
            ss << "featureA" << i << ".exr";
            features[i].bufferA->save(ss.str(), false);
        }
        {
            std::stringstream ss; 
            ss << "featureB" << i << ".exr";
            features[i].bufferB->save(ss.str(), false);
        }
        {
            std::stringstream ss; 
            ss << "featureVariance" << i << ".exr";
            features[i].bufferVariance->save(ss.str(), false);
        }
    }*/

    Tungsten::Timer timer;
    Pixmap3f result = nforDenoiser(std::move(image), std::move(features));
    timer.stop();
    printTimestampedLog(tfm::format("Filtering complete! Filter time: %.1fs", timer.elapsed()));    
    //result.save("result.exr", true);

    // Write output
    Spectrum *recons = (Spectrum *)reconsBitmap->getFloatData();
    for (int i = 0; i < size.x * size.y; ++i) {
        Float val[] = { std::max(Float(result[i][0]), Float(0.0)),
                        std::max(Float(result[i][1]), Float(0.0)),
                        std::max(Float(result[i][2]), Float(0.0)) };
        recons[i] = Spectrum(val);
    }
}

/**
 * Denoise with only estimated gradients 
 */
void Nfor::denoise(AccumBuffer* throughput, Bitmap *reconsBitmap, AccumBuffer* featureBuffer, std::vector<int> bufferIndices, bool customRecon) {
    
    ThreadUtils::startThreads(std::max(ThreadUtils::idealThreadCount() - 1, 1u));

    printTimestampedLog(tfm::format("Filtering with estimated gradients..."));

    RenderBuffer3f image;
    std::vector<RenderBufferF> features;
    
    Spectrum *throughputPtr  = (Spectrum *)throughput->getBuffer(0)->getFloatData();
    Spectrum *variancePtr    = (Spectrum *)throughput->getSampleVariance(0)->getFloatData();
    Spectrum *throughputAPtr = (Spectrum *)throughput->getBufferA(0)->getFloatData();
    Spectrum *throughputBPtr = (Spectrum *)throughput->getBufferB(0)->getFloatData();

    Vector2i size = throughput->getSize();

    image.buffer = std::unique_ptr<Pixmap3f>(new Pixmap3f(size.x, size.y));
    image.bufferVariance = std::unique_ptr<Pixmap3f>(new Pixmap3f(size.x, size.y));
    image.bufferA = std::unique_ptr<Pixmap3f>(new Pixmap3f(size.x, size.y));
    image.bufferB = std::unique_ptr<Pixmap3f>(new Pixmap3f(size.x, size.y));

    for (int j = 0; j < size.x * size.y; ++j) {        
        (*image.buffer)[j] = Vec3f(
            (float)throughputPtr[j][0], 
            (float)throughputPtr[j][1], 
            (float)throughputPtr[j][2]);

        (*image.bufferVariance)[j] = Vec3f(
            (float)variancePtr[j][0],
            (float)variancePtr[j][1],
            (float)variancePtr[j][2]);

        (*image.bufferA)[j] = Vec3f(
                (float)throughputAPtr[j][0],
                (float)throughputAPtr[j][1],
                (float)throughputAPtr[j][2]);

        (*image.bufferB)[j] = Vec3f(
                (float)throughputBPtr[j][0],
                (float)throughputBPtr[j][1],
                (float)throughputBPtr[j][2]);
    }

    // Assume each auxiliary buffer has 3 channels
    AccumBuffer* auxBuffer = featureBuffer;    
    std::vector<int> auxBufferIds = bufferIndices;

    features.resize(auxBufferIds.size() * 3);
    for (size_t i = 0; i < features.size(); ++i) {
        // The main buffer is not used in the filtering process of traditional NFOR though
        features[i].buffer = std::unique_ptr<PixmapF>(new PixmapF(size.x, size.y));
        features[i].bufferA = std::unique_ptr<PixmapF>(new PixmapF(size.x, size.y));
        features[i].bufferB = std::unique_ptr<PixmapF>(new PixmapF(size.x, size.y));
        features[i].bufferVariance = std::unique_ptr<PixmapF>(new PixmapF(size.x, size.y));
    }

    for (auto i = 0; i < auxBufferIds.size(); i++) {
        int k = auxBufferIds[i];
        for (auto j = 0; j < size.x * size.y; ++j) {

            Vector &buffer = auxBuffer->getBuffer(k).get()->getV3Data()[j];
            (*features[i*3+0].buffer)[j] = buffer[0];
            (*features[i*3+1].buffer)[j] = buffer[1];
            (*features[i*3+2].buffer)[j] = buffer[2];

            Vector &bufferA = auxBuffer->getBufferA(k).get()->getV3Data()[j];
            (*features[i*3+0].bufferA)[j] = bufferA.x;
            (*features[i*3+1].bufferA)[j] = bufferA.y;
            (*features[i*3+2].bufferA)[j] = bufferA.z;

            Vector &bufferB = auxBuffer->getBufferB(k).get()->getV3Data()[j];
            (*features[i*3+0].bufferB)[j] = bufferB.x;
            (*features[i*3+1].bufferB)[j] = bufferB.y;
            (*features[i*3+2].bufferB)[j] = bufferB.z;

            Vector &variance = auxBuffer->getSampleVariance(k).get()->getV3Data()[j];
            (*features[i*3+0].bufferVariance)[j] = variance.x;
            (*features[i*3+1].bufferVariance)[j] = variance.y;
            (*features[i*3+2].bufferVariance)[j] = variance.z;
        }
    }

    // Dump data to see
    /*
    std::cout << "Prepare to call save" << std::endl;
    (*image.buffer).save("buffer.exr", true);
    image.bufferA->save("bufferA.exr", true);
    image.bufferB->save("bufferB.exr", true);
    image.bufferVariance->save("bufferVariance.exr", true);

    for (int i = 0; i < features.size(); ++i) {
        {
            std::stringstream ss; 
            ss << "featureA" << i << ".exr";
            features[i].bufferA->save(ss.str(), false);
        }
        {
            std::stringstream ss; 
            ss << "featureB" << i << ".exr";
            features[i].bufferB->save(ss.str(), false);
        }
        {
            std::stringstream ss; 
            ss << "featureVariance" << i << ".exr";
            features[i].bufferVariance->save(ss.str(), false);
        }
    }*/
    

    Tungsten::Timer timer;
    Pixmap3f result; 
    if (!customRecon)   
        // traditional NFOR
        result = nforDenoiser(std::move(image), std::move(features));
    else
        // custom NFOR with feature map from gradient reconstruction 
        result = nforDenoiserReconstruction(std::move(image), std::move(features));

    timer.stop();
    printTimestampedLog(tfm::format("Filtering complete! Filter time: %.1fs", timer.elapsed()));    
    //result.save("result.exr", true);

    // Write output
    Spectrum *recons = (Spectrum *)reconsBitmap->getFloatData();
    for (int i = 0; i < size.x * size.y; ++i) {
        Float val[] = { std::max(Float(result[i][0]), Float(0.0)),
                        std::max(Float(result[i][1]), Float(0.0)),
                        std::max(Float(result[i][2]), Float(0.0)) };
        recons[i] = Spectrum(val);
    }
}

#pragma once

#include <mitsuba/core/spectrum.h>

MTS_NAMESPACE_BEGIN

// One gather points at given pixel

class PixelData {
public:
  // Position informations
  Intersection its;  //< Intersection information
  int depth;

  // BSDF gather point dependances
  Spectrum weight;  //< Importance associate
  int sampledComponent;
  float pdfComponent;

  // Gather point associated data
  Point2i pos;  // Position on the image (pixel wise)
  Spectrum fluxDirect; // Holds flux from path-tracing for given pixel
  Spectrum flux;
  Float scale; // Radius scale
  Float radius; // Radius of gather points in this list
  Float N; // Number of gathered photons (for SPPM only)
  size_t nPhotons; // Number of gathered photons without radius reduction

  // === Sync thread data - only 1 alloc per list!
  Float *tempM;
  Spectrum *tempFlux; // To avoid thread sync

  // Max thread information
  int maxThread;
  // Directly visible emission
  Spectrum emission;
  // Cumulated importance
  Float cumulImportance;

  inline void rescaleRadii(Float v) {
    flux *= v;
  }

  static int allocTempSize(int mT) {
    return mT * (sizeof(Float) + // M Statistic = Number of photons
        sizeof(Spectrum)); // Flux collected statistic
  }

  // Alloc gather point inner structures
  void allocTemp(int mT, char *&allocPtr) {
    if (maxThread != -1) // Already allocated?
      return;

    // Ptr settings
    tempFlux = (Spectrum *) allocPtr;
    allocPtr += sizeof(Spectrum) * mT;
    tempM = (Float *) allocPtr;
    allocPtr += sizeof(Float) * mT;

    maxThread = mT;
  }

  /// Reset the temp value associated to the gather point.
  inline void resetTemp() {
    memset(tempFlux, 0, sizeof(Spectrum) * maxThread);
    memset(tempM, 0, sizeof(Float) * maxThread);
  }

  inline Spectrum getFlux() const {
    Spectrum flux(0.f);
    for (int idThread = 0; idThread < maxThread; idThread++) {
      flux += tempFlux[idThread];
    }
    return flux;
  }

  inline Float getM() const {
    Float M = 0.f;
    for (int idThread = 0; idThread < maxThread; idThread++) {
      M += tempM[idThread];
    }
    return M;
  }

  inline PixelData() :
      pos(-1, -1),
      fluxDirect(0.f),
      flux(0.f),
      N(0.f),
      nPhotons(0),
      tempM(NULL),
      tempFlux(NULL),
      maxThread(-1),
      emission(0.f),
      cumulImportance(0.f) {
  }

};

MTS_NAMESPACE_END


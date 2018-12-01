#include "shift_ME.h"

// To make possible to use some statistics
#include <mitsuba/core/statistics.h>

MTS_NAMESPACE_BEGIN

StatsCounter MESuccess("GPM",
                       "Manifold Exp. sucessful shift", EPercentage);
StatsCounter MEShiftMEFailed("GPM",
                             "ShiftME failure shift", EPercentage);

bool generateShiftPathME(const Path &source, Path &proposal, size_t b, size_t c, MemoryPool &pool,
                         ManifoldPerturbation *mePerturb, const PathVertex &shiftVertex, Float radius,
                         Point gpBasePos, Point gpOffsetPos) {
  SAssert(b >= 1);
  SAssert(b < c);
  /* Allocate memory for the proposed path */
  proposal.clear();
  proposal.append(source,
                  0,
                  b);        // Note that vertex b cannot be shared, because the pdf[EImportance] will be changed.
  proposal.append(source.edge(b - 1)->clone(pool),
                  source.vertex(b)->clone(pool));
  for (size_t i = b + 1; i <= c; ++i) {
    proposal.append(pool.allocEdge());            // outgoing edge at vertex b - 1
    memset(proposal.edge(proposal.edgeCount() - 1), 0, sizeof(PathEdge));  // make sure edges are ALWAYS initialized!

    proposal.append(pool.allocVertex());          // next vertex
    memset(proposal.vertex(proposal.vertexCount() - 1), 0, sizeof(PathVertex));
    proposal.vertex(proposal.vertexCount() - 1)->sampledComponentIndex = source.vertex(i)->sampledComponentIndex;
  }

  // Assign information for c
  *proposal.vertex(c) = shiftVertex;

  // Manifold walk between b .. c
  MESuccess.incrementBase();
  bool walkSuccess = mePerturb->manifoldWalkGPM(source, proposal, 1, // towards sensor
                                                b, c, radius, gpBasePos, gpOffsetPos);
  if (!walkSuccess) {
    return false;
  }
  ++MESuccess;

  // Sanity check
  for (size_t i = b; i < c; ++i) {
    if (source.vertex(i)->sampledComponentIndex != proposal.vertex(i)->sampledComponentIndex) {
      return false;
      //SLog(EError, "Sample component index not match during ME");
    }
  }

  // NOTE: after manifold walk, only the position of each vertex is generated and some cache values such as pdf and weight of a vertex or edge on
  // the proposal path are available.
  // So here, we calculate some additional values that we need later
  for (size_t i = b; i < c; ++i) {
    Vector d = proposal.vertex(i + 1)->getPosition() - proposal.vertex(i)->getPosition();
    Float len = d.length();
    proposal.edge(i)->d = d / len;
    proposal.edge(i)->length = len;

    // assume same medium
    proposal.edge(i)->medium = source.edge(i)->medium;
  }

  // Verify and recompute all cache values in the proposal path (edge->d, for example)
  /*
  if (scene) {
    std::stringstream ss;
    proposal.verify(scene, EImportance, ss);
    SLog(EInfo, "Verify proposal path result b = %d c = %d", b, c);
    SLog(EInfo, "%s", ss.str().c_str());
  }*/

  return true;
}

bool ShiftME(ShiftRecord &sRec, const Path &source, const Path &proposal, size_t b, size_t c,
             bool isBeam) {
  SAssert(b >= 1);
  SAssert(b < c);

  sRec.throughtput = Spectrum(1.f);
  for (size_t i = 0; i < b; ++i) {
    sRec.throughtput *= proposal.vertex(i)->weight[EImportance] *
        proposal.vertex(i)->rrWeight *
        proposal.edge(i)->weight[EImportance];
  }

  sRec.pdf = 1.f;
  for (size_t i = b; i < c; ++i) {
    // For the vertex
    {
      Float proposalPdf = proposal.vertex(i)->pdf[EImportance];
      Float basePdf = source.vertex(i)->pdf[EImportance];
      if (basePdf == 0.f) {
        SLog(EWarn, "0 PDF base path ME");
        return false;
      }

      sRec.throughtput *= proposal.vertex(i)->weight[EImportance] * (proposalPdf / basePdf);
      sRec.pdf *= proposalPdf;
    }

    if ((i == c - 1) & isBeam) {
      // Need to change the pdf to the right one (PDFFailure)
      // Only for the proposal PDF
      MediumSamplingRecord mRecShift;
      Ray mRay(proposal.vertex(i)->getPosition(), proposal.edge(i)->d, 0.f, proposal.edge(i)->length, 0.f);
      source.edge(i)->medium->eval(mRay, mRecShift);

      if (source.edge(i)->pdf[EImportance] == 0.f) {
        SLog(EWarn, "0 PDF base path ME");
        return false;
      }

      sRec.throughtput *= mRecShift.transmittance;
      sRec.pdf *= mRecShift.pdfFailure; // add the prob of pass through the medium
      sRec.throughtput /= source.edge(i)->pdf[EImportance];
    } else {
      Float proposalPdfDist = proposal.edge(i)->pdf[EImportance];
      Float basePdfDist = source.edge(i)->pdf[EImportance];
      if (basePdfDist == 0.f) {
        SLog(EWarn, "0 PDF base path ME");
        return false;
      }

      sRec.throughtput *= proposal.edge(i)->weight[EImportance] * (proposalPdfDist / basePdfDist);
      sRec.pdf *= proposalPdfDist;
    }

    sRec.throughtput *= source.vertex(i)->rrWeight;
  }

  MEShiftMEFailed.incrementBase();
  if (sRec.throughtput.isZero() || sRec.pdf == 0.f) {
    ++MEShiftMEFailed;
    //SLog(EWarn, "Impossible to evaluate this path");
    return false;
  }

  return true;
}

MTS_NAMESPACE_END
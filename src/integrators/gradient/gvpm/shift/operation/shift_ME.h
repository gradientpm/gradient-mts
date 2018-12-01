#pragma once

#include "../shift_utilities.h"

#ifndef MITSUBA_SHIFT_ME_H
#define MITSUBA_SHIFT_ME_H

MTS_NAMESPACE_BEGIN

bool generateShiftPathME(const Path &source,
                         Path &proposal,
                         size_t b,
                         size_t c,
                         MemoryPool &pool,
                         ManifoldPerturbation *mePerturb,
                         const PathVertex &shiftVertex,
                         Float radius,
                         Point gpBasePos,
                         Point gpOffsetPos);
bool ShiftME(ShiftRecord &sRec, const Path &source, const Path &proposal, size_t b, size_t c,
             bool isBeam = false);

MTS_NAMESPACE_END;

#endif //MITSUBA_SHIFT_ME_H

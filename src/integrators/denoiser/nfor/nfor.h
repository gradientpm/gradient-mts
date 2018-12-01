#pragma once

#include <mitsuba/mitsuba.h>
#include <mitsuba/core/bitmap.h>
using namespace mitsuba;

#include "../accum_buffer.h"
#include "../aux_buffer.h"

class Nfor : public Object {
public: 
    void denoise(AccumBuffer* throughput, Bitmap *reconsBitmap, AccumBuffer* featureBuffer, std::vector<int> bufferIndices, bool customRecon = false);    
    void denoise(AccumBuffer* throughput, Bitmap *reconsBitmap, AuxiliaryBuffer* auxBuffer);
};

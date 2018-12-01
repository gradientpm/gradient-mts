#pragma once

#include <mitsuba/mitsuba.h>
#include "accum_buffer.h"


MTS_NAMESPACE_BEGIN

    enum EAuxBuffer {
        PositionBuffer,
        NormalBuffer,
        AlbedoBuffer,
        ShadingNormalW,
        RayW,
        AuxBufferCount
    };

struct AuxiliaryBuffer : public SerializableObject {

    AuxiliaryBuffer(const Vector2i &s, bool trackVariance, bool twoBufferVariance) : size(s) {
        accumBuffer = new AccumBuffer(EAuxBuffer::AuxBufferCount, size, trackVariance, twoBufferVariance);
    }

    /// Unserialize a serializable object
    AuxiliaryBuffer(Stream *stream, InstanceManager *manager) {
    }

    /// Serialize this object to a stream
    virtual void serialize(Stream *stream, InstanceManager *manager) const {
    }

    void clear() {
        accumBuffer.get()->clear();
    }

    const Spectrum getValue(int idx, EAuxBuffer type) const{
        return accumBuffer.get()->getValue(idx, type);
    }

    ref <Bitmap> getBuffer(int idx) const{
        return accumBuffer.get()->getBuffer(idx);
    }

    ref <Bitmap> getBufferA(int idx) const{
        return accumBuffer.get()->getBufferA(idx);
    }

    ref <Bitmap> getBufferB(int idx) const{
        return accumBuffer.get()->getBufferB(idx);
    }

    ref <Bitmap> getVariance(int idx) const{
        return accumBuffer.get()->getSampleVariance(idx);
    }

    void addSample(int x, int y, const Vector& value, EAuxBuffer type){
        Float val[] = {value[0], value[1], value[2]};
        Spectrum v(val);
        accumBuffer.get()->addSample(x, y, v, type);
    }

private:
    ref<AccumBuffer> accumBuffer;
    Vector2i size;
};

MTS_NAMESPACE_END
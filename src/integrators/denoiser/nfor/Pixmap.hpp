#ifndef PIXMAP_HPP_
#define PIXMAP_HPP_

#include "core/Memory.hpp"
#include "core/math/Vec.hpp"
#include "core/sse/SimdFloat.hpp"

#include <cstring>

#include <ImfInputFile.h>
#include <ImfOutputFile.h>
#include <ImfRgbaFile.h>
#include <ImfChannelList.h>
#include <ImfFrameBuffer.h>
#include <ImfIO.h>
#include <half.h>

namespace Tungsten {

template<typename Texel>
class Pixmap
{
    int _w, _h;
    aligned_unique_ptr<Texel> _pixels;

public:
    Pixmap() = default;

    Pixmap(int w, int h, const Texel *src = nullptr)
    : _w(w), _h(h),
      _pixels(src ? alignedAlloc<Texel>(w*h, 16) : alignedZeroAlloc<Texel>(w*h, 16))
    {
        if (src)
            std::memcpy(_pixels.get(), src, w*h*sizeof(Texel));
    }

    Pixmap(int w, int h, const Pixmap<Texel> &src)
    : _w(w), _h(h),
      _pixels( alignedAlloc<Texel>(w*h, 16) )
    {        
        std::memcpy(_pixels.get(), src._pixels.get(), w*h*sizeof(Texel));
    }
    
    bool save(std::string filename, bool rgb = true) const
    {    
        int w = _w;
        int h = _h;      
        int channels = rgb ? 3 : 1;
        
        const float* img = reinterpret_cast<const float *>(_pixels.get());
        
        try {

        Imf::Header header(w, h, 1.0f, Imath::V2f(0, 0), 1.0f, Imf::INCREASING_Y, Imf::PIZ_COMPRESSION);
        Imf::FrameBuffer frameBuffer;
        
        std::unique_ptr<half[]> data(new half[w*h*channels]);
        for (int i = 0; i < w*h*channels; ++i)
            data[i] = half(img[i]);

        const char *channelNames[] = {"R", "G", "B", "A"};
        for (int i = 0; i < channels; ++i) {
            const char *channelName = (channels == 1) ? "Y" : channelNames[i];
            header.channels().insert(channelName, Imf::Channel(Imf::HALF));
            frameBuffer.insert(channelName, Imf::Slice(Imf::HALF, reinterpret_cast<char *>(data.get() + i),
                    sizeof(half)*channels, sizeof(half)*channels*w));
        }

        Imf::OutputFile file(filename.c_str(), header);
        file.setFrameBuffer(frameBuffer);
        file.writePixels(h);

        return true;

        } catch(const std::exception &e) {
            std::cout << "OpenEXR writer failed: " << e.what() << std::endl;
            return false;
        }
    }

    void clear()
    {
        std::memset(_pixels.get(), 0, _w*_h*sizeof(Texel));
    }

    void reset()
    {
        _w = _h = 0;
        _pixels.reset();
    }

    inline Texel &operator[](int idx)       { return _pixels[idx]; }
    inline Texel  operator[](int idx) const { return _pixels[idx]; }

    inline Texel &operator()(int x, int y)       { return _pixels[x + y*_w]; }
    inline Texel  operator()(int x, int y) const { return _pixels[x + y*_w]; }

    inline Texel &operator[](Vec2i p)       { return operator()(p.x(), p.y()); }
    inline Texel  operator[](Vec2i p) const { return operator()(p.x(), p.y()); }

    int w() const { return _w; }
    int h() const { return _h; }
};

typedef Pixmap<float4> Pixmap4pf;
typedef Pixmap<Vec3f> Pixmap3f;
typedef Pixmap<float> PixmapF;

/*
template<typename Texel>
std::unique_ptr<Pixmap<Texel>> loadPixmap(const Path &path, bool rgb = true)
{
    auto texelType = rgb ? TexelConversion::REQUEST_RGB : TexelConversion::REQUEST_AVERAGE;
    int w, h;
    std::unique_ptr<float[]> pixels = ImageIO::loadHdr(path, texelType, w, h);

    if (!pixels)
        return nullptr;

    return std::unique_ptr<Pixmap<Texel>>(new Pixmap<Texel>(w, h, reinterpret_cast<Texel *>(pixels.get())));
}*/

template<typename Texel>
class PixmapIterator
{
    Pixmap<Texel> &_pixmap;
    int _idx;

public:
    PixmapIterator(Pixmap<Texel> &pixmap, uint64 idx) : _pixmap(pixmap), _idx(idx) {}

    Texel &operator*() { return _pixmap[_idx]; }

    bool operator!=(const PixmapIterator &o) const { return _idx != o._idx; }

    PixmapIterator &operator++() { _idx++; return *this; }
    PixmapIterator operator++(int) { PixmapIterator copy(*this); ++(*this); return copy; }
};

template<typename Texel>
class ConstPixmapIterator
{
    const Pixmap<Texel> &_pixmap;
    int _idx;

public:
    ConstPixmapIterator(const Pixmap<Texel> &pixmap, uint64 idx) : _pixmap(pixmap), _idx(idx) {}

    const Texel &operator*() { return _pixmap[_idx]; }

    bool operator!=(const ConstPixmapIterator &o) const { return _idx != o._idx; }

    ConstPixmapIterator &operator++() { _idx++; return *this; }
    ConstPixmapIterator operator++(int) { ConstPixmapIterator copy(*this); ++(*this); return copy; }
};

template<typename Texel>
PixmapIterator<Texel> begin(Pixmap<Texel> &pixmap)
{
    return PixmapIterator<Texel>(pixmap, 0);
}
template<typename Texel>
PixmapIterator<Texel> end(Pixmap<Texel> &pixmap)
{
    return PixmapIterator<Texel>(pixmap, pixmap.w()*pixmap.h());
}

template<typename Texel>
ConstPixmapIterator<Texel> begin(const Pixmap<Texel> &pixmap)
{
    return ConstPixmapIterator<Texel>(pixmap, 0);
}
template<typename Texel>
ConstPixmapIterator<Texel> end(const Pixmap<Texel> &pixmap)
{
    return ConstPixmapIterator<Texel>(pixmap, pixmap.w()*pixmap.h());
}

}

#endif /* PIXMAP_HPP_ */

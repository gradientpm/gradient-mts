/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2014 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#if !defined(__GDVCM_H)
#define __GDVCM_H

#include <mitsuba/mitsuba.h>
#include "../../vcm/vcm_basics.h"

#define SEPARATE_DIRECT // use separate buffer for LS*E paths

MTS_NAMESPACE_BEGIN

/* ==================================================================== */
/*                         Configuration storage                        */
/* ==================================================================== */

/**
 * \brief Stores all configuration parameters of the
 * bidirectional path tracer
 */
struct GDVCMConfiguration : VCMConfigBase {
	int blockSize;
	bool lightImage;
	bool sampleDirect;
	Vector2i cropSize;

	int extraBorder;
	int nNeighbours;

	float m_shiftThreshold;
	bool m_reconstructL1;
	bool m_reconstructL2;
	bool m_reconstructUni;
	float m_reconstructAlpha;
	Float radiusReductionAlpha;
	bool directTracing;

	bool forceBlackPixels;
	int maxManifoldIterations;

	inline GDVCMConfiguration() { }

	inline GDVCMConfiguration(Stream *stream) {
		maxDepth = stream->readInt();
		blockSize = stream->readInt();
		lightImage = stream->readBool();
		sampleDirect = stream->readBool();
		cropSize = Vector2i(stream);
		rrDepth = stream->readInt();
                initialRadius = stream->readFloat();
                radiusReductionAlpha = stream->readFloat();
		extraBorder = stream->readInt();
		nNeighbours = stream->readInt();

		m_shiftThreshold = stream->readFloat();
		m_reconstructL1 = stream->readBool();
		m_reconstructL2 = stream->readBool();
		m_reconstructAlpha = stream->readFloat();
                mergeOnly = stream->readBool();
                phExponent = stream->readFloat();
	}

	inline void serialize(Stream *stream) const {
		stream->writeInt(maxDepth);
		stream->writeInt(blockSize);
		stream->writeBool(lightImage);
		stream->writeBool(sampleDirect);
		cropSize.serialize(stream);
		stream->writeInt(rrDepth);	//possible problem with network rendering?
                stream->writeFloat(initialRadius);
                stream->writeFloat(radiusReductionAlpha);

		stream->writeInt(extraBorder);
		stream->writeInt(nNeighbours);

		stream->writeFloat(m_shiftThreshold);
		stream->writeBool(m_reconstructL1);
		stream->writeBool(m_reconstructL2);
		stream->writeFloat(m_reconstructAlpha);
                stream->writeBool(mergeOnly);
		stream->writeFloat(phExponent);
	}

	void dump() const {
		SLog(EInfo, "Gradient-Domain Bidirectional Path Tracer configuration:");
		SLog(EDebug, "   Maximum path depth          : %i", maxDepth);
		SLog(EDebug, "   Image size                  : %ix%i",
			cropSize.x, cropSize.y);
		SLog(EDebug, "   Generate light image        : %s",
			lightImage ? "yes" : "no");
		SLog(EDebug, "   Russian roulette depth      : %i", rrDepth);
		SLog(EDebug, "   Block size                  : %i", blockSize);
	}

	bool accumulateData(ref<Bitmap> buff, ref<Film> film, int bufferIdx, int target, int iter, const std::vector<Float> &weights);
	void prepareDataForSolver(float w, float* out, Float * data, int len, Float *data2 = NULL, int offset = 0);
	void setBitmapFromArray(ref<Bitmap> &bitmap, float *img);
};

MTS_NAMESPACE_END

#endif /* __GDVCM_H */

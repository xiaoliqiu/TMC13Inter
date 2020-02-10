/* The copyright in this software is being made available under the BSD
 * Licence, included below.  This software may be subject to other third
 * party and contributor rights, including patent rights, and no such
 * rights are granted under this licence.
 *
 * Copyright (c) 2017-2018, ISO/IEC
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * * Redistributions of source code must retain the above copyright
 *   notice, this list of conditions and the following disclaimer.
 *
 * * Redistributions in binary form must reproduce the above copyright
 *   notice, this list of conditions and the following disclaimer in the
 *   documentation and/or other materials provided with the distribution.
 *
 * * Neither the name of the ISO/IEC nor the names of its contributors
 *   may be used to endorse or promote products derived from this
 *   software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#pragma once

#include <functional>
#include <map>
#include <string>
#include <vector>

#include "PayloadBuffer.h"
#include "PCCMath.h"
#include "PCCPointSet.h"
#include "hls.h"

namespace pcc {

//----------------------------------------------------------------------------

struct McParam {
  int Amotion0;
  double lambda;
  int decimate;
  int sampleW;
  int globalMotionInRdo;
};

//----------------------------------------------------------------------------

struct EncoderParams {
  SequenceParameterSet sps;
  GeometryParameterSet gps;

  // Motion for octrees
  McParam motion;

  // NB: information about attributes is split between the SPS and the APS.
  //  => The SPS enumerates the attributes, the APS controls coding params.
  std::vector<AttributeParameterSet> aps;

  // Maps indexes of sps.attributeSet[*] to aps[*] in order to select an
  // APS for each enumerated attribute.
  std::vector<int> attributeSetIdxToApsIdxMap;

  // Indicates the next frame to be encoded is a random access point
  bool randomAccessPoint;

  // Period of random access points (managed by SequenceEncoder)
  int randomAccessPeriod;
};

//============================================================================

class PCCTMC3Encoder3 {
public:
  class Callbacks;

  PCCTMC3Encoder3() { init(); }
  PCCTMC3Encoder3(const PCCTMC3Encoder3&) = default;
  PCCTMC3Encoder3& operator=(const PCCTMC3Encoder3& rhs) = default;
  ~PCCTMC3Encoder3() = default;

  // Configure the encoder (SPS,GPS,APS)
  void configure(EncoderParams* params);

  void init();

#if INTER_HIERARCHICAL
  void setPOC(int pocSeq) { poc = pocSeq; }
  void setFirstFramePOC(int firstFramePOCSeq) { firstFramePOC = firstFramePOCSeq; }
  void setLastGOP(bool lastGOPSeq) { lastGOP = lastGOPSeq; }
#endif

  int compress(
    const PCCPointSet3& inputPointCloud,
    const EncoderParams& params,
    Callbacks*,
    PCCPointSet3* reconstructedCloud = nullptr);

private:
  void reconstructedPointCloud(PCCPointSet3* reconstructedCloud);

  void encodeGeometryBrick(const EncoderParams& params, PayloadBuffer* buf);

  void computeMinPositions(const PCCPointSet3& inputPointCloud);

  void quantization(const PCCPointSet3& inputPointCloud);

private:
  // todo(df): minPositions is unscaled -- which isn't quite correct.
  PCCVector3D minPositions;
  PCCBox3<uint32_t> boundingBox;
  PCCPointSet3 pointCloud;

  // Point cloud that acts as a predictor of @pointCloud's geometry
  // occupancy.
  PCCPointSet3 predPointCloud;

#if INTER_HIERARCHICAL
  PCCPointSet3 backPredPointCloud;
  PCCPointSet3 recGOPPointCloud[9];

  int poc;
  int firstFramePOC;
  bool lastGOP;
#endif

  // The active parameter sets
  const SequenceParameterSet* _sps;
  const GeometryParameterSet* _gps;

  // configured parameter sets
  std::map<int, SequenceParameterSet> _spss;
  std::map<int, GeometryParameterSet> _gpss;
  std::map<int, AttributeParameterSet> _apss;

  // the list of attributes to code and which APS to use for each
  //  ie, maps attr_idx -> apsid, for attr_idx in sps.
  std::vector<int> _attrIdxToApsIdMap;
};

//----------------------------------------------------------------------------

class PCCTMC3Encoder3::Callbacks {
public:
  virtual void onOutputBuffer(const PayloadBuffer&) = 0;
  virtual void onPostRecolour(const PCCPointSet3&) = 0;
};

//============================================================================

}  // namespace pcc

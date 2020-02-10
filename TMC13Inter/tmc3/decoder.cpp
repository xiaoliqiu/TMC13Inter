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

#include "PCCTMC3Decoder.h"

#include <cassert>
#include <string>

#include "AttributeDecoder.h"
#include "PayloadBuffer.h"
#include "PCCPointSet.h"
#include "geometry.h"
#include "io_hls.h"
#include "io_tlv.h"
#include "pcc_chrono.h"
#include "osspecific.h"

namespace pcc {

//============================================================================

void
PCCTMC3Decoder3::init()
{
  _sps = nullptr;
  _gps = nullptr;
  _spss.clear();
  _gpss.clear();
  _apss.clear();
}

//============================================================================

int
PCCTMC3Decoder3::decompress(
  const DecoderParams& params,
#if INTER_HIERARCHICAL
  const int firstFrameNumSeq,
  const int frameCountSeq,
#endif
  const PayloadBuffer* buf,
  PCCTMC3Decoder3::Callbacks* callback)
{
  if (!buf) {
    // flush decoder, output pending cloud if any
#if INTER_HIERARCHICAL
    callback->onOutputCloud(_currentPointCloud, poc);
#else
    callback->onOutputCloud(_currentPointCloud);
#endif
    return 0;
  }

  switch (buf->type) {
  case PayloadType::kSequenceParameterSet: storeSps(parseSps(*buf)); return 0;

  case PayloadType::kGeometryParameterSet: storeGps(parseGps(*buf)); return 0;

  case PayloadType::kAttributeParameterSet: storeAps(parseAps(*buf)); return 0;

  case PayloadType::kGeometryBrick:
    // call onOutputCloud when starting a new cloud
    if (_currentPointCloud.getPointCount()) {
#if INTER_HIERARCHICAL
      callback->onOutputCloud(_currentPointCloud, poc);
#else
      callback->onOutputCloud(_currentPointCloud);
#endif
      // save the decoded pointcloud in the reference buffer

#if INTER_HIERARCHICAL
      firstFrameNum = firstFrameNumSeq;
      frameCount = frameCountSeq;
      int randomaccesPeriod = 8;

      if ( !lastGOP )
      {
        int index = (poc - firstFrameNum) % randomaccesPeriod;

        if (index != 0)
        {
          recGOPPointCloud[index] = _currentPointCloud;
        }
        else if (poc == firstFrameNum)
        {
          recGOPPointCloud[0] = _currentPointCloud;
          recGOPPointCloud[8] = _currentPointCloud;
        }
        else // for other cases
        {
          recGOPPointCloud[0] = recGOPPointCloud[8];
          recGOPPointCloud[8] = _currentPointCloud;
        }
      }
      else if ( lastGOP )
      {
        recGOPPointCloud[8] = _currentPointCloud;
      }

      const int lastFrameNum = firstFrameNum + frameCount;
      frameNum++;

      if ((frameNum + randomaccesPeriod - 1 - firstFrameNum) % randomaccesPeriod == 0)
      {
        if (frameNum + randomaccesPeriod > lastFrameNum)
        {
          lastGOP = true;
        }
    }

      if (!lastGOP)
      {
        int frameIndex = (frameNum - firstFrameNum - 1) % randomaccesPeriod;

        int pocIndex = 0;
        switch (frameIndex)
        {
        case 0: pocIndex = 8; break;
        case 1: pocIndex = 4; break;
        case 2: pocIndex = 2; break;
        case 3: pocIndex = 1; break;
        case 4: pocIndex = 3; break;
        case 5: pocIndex = 6; break;
        case 6: pocIndex = 5; break;
        case 7: pocIndex = 7; break;
        default:
          break;
        }

        poc = (frameNum - firstFrameNum - 1) / randomaccesPeriod * randomaccesPeriod + firstFrameNum + pocIndex;
      }
      else if (lastGOP)
      {
        poc = frameNum;
      }
#else
      std::swap(_refPointCloud, _currentPointCloud);
#endif
    }
    return decodeGeometryBrick(*buf);

  case PayloadType::kAttributeBrick: decodeAttributeBrick(*buf); return 0;
  }

  // todo(df): error, unhandled payload type
  return 1;
}

//--------------------------------------------------------------------------

void
PCCTMC3Decoder3::storeSps(SequenceParameterSet&& sps)
{
  // todo(df): handle replacement semantics
  _spss.emplace(std::make_pair(sps.sps_seq_parameter_set_id, sps));
}

//--------------------------------------------------------------------------

void
PCCTMC3Decoder3::storeGps(GeometryParameterSet&& gps)
{
  // todo(df): handle replacement semantics
  _gpss.emplace(std::make_pair(gps.gps_geom_parameter_set_id, gps));
}

//--------------------------------------------------------------------------

void
PCCTMC3Decoder3::storeAps(AttributeParameterSet&& aps)
{
  // todo(df): handle replacement semantics
  _apss.emplace(std::make_pair(aps.aps_attr_parameter_set_id, aps));
}

//==========================================================================
// Initialise the point cloud storage and decode a single geometry brick.
//
// NB: there is an implicit assumption that there is only a single geometry
//     brick per point cloud.
// todo(??): fixme

int
PCCTMC3Decoder3::decodeGeometryBrick(const PayloadBuffer& buf)
{
  assert(buf.type == PayloadType::kGeometryBrick);
  std::cout << "positions bitstream size " << buf.size() << " B\n";

  // HACK: assume activation of the first SPS and GPS
  // todo(df): parse brick header here for propper sps & gps activation
  //  -- this is currently inconsistent between trisoup and octree
  assert(!_spss.empty());
  assert(!_gpss.empty());
  _sps = &_spss.cbegin()->second;
  _gps = &_gpss.cbegin()->second;

  // todo(df): replace with attribute mapping
  bool hasColour = std::any_of(
    _sps->attributeSets.begin(), _sps->attributeSets.end(),
    [](const AttributeDescription& desc) {
      return desc.attributeLabel == KnownAttributeLabel::kColour;
    });

  bool hasReflectance = std::any_of(
    _sps->attributeSets.begin(), _sps->attributeSets.end(),
    [](const AttributeDescription& desc) {
      return desc.attributeLabel == KnownAttributeLabel::kReflectance;
    });

  _currentPointCloud.clear();
  _currentPointCloud.addRemoveAttributes(hasColour, hasReflectance);

  pcc::chrono::Stopwatch<pcc::chrono::utime_inc_children_clock> clock_user;
  clock_user.start();

  int gbhSize;
  GeometryBrickHeader gbh = parseGbh(*_gps, buf, &gbhSize);
  minPositions.x() = gbh.geomBoxOrigin.x();
  minPositions.y() = gbh.geomBoxOrigin.y();
  minPositions.z() = gbh.geomBoxOrigin.z();

#if INTER_HIERARCHICAL
  if (gbh.random_access_point) {
    _refPointCloud.clear();
    _backRefPointCloud.clear();
  }
  else if ( !lastGOP ) {
    int pocIndex = (poc - firstFrameNum) % 8;                       // poc supposed to be parsed. it is a fixed value currently
    if (pocIndex == 4)
    {
      _refPointCloud = recGOPPointCloud[0];
      _backRefPointCloud = recGOPPointCloud[8];
    }
    else if (pocIndex == 2 || pocIndex == 6)
    {
      _refPointCloud = recGOPPointCloud[pocIndex - 2];
      _backRefPointCloud = recGOPPointCloud[pocIndex + 2];
    }
    else
    {
      _refPointCloud = recGOPPointCloud[pocIndex - 1];
      _backRefPointCloud = recGOPPointCloud[pocIndex + 1];
    }
  }
  else {
    _refPointCloud = recGOPPointCloud[8];
    _backRefPointCloud = recGOPPointCloud[8];
  }
#else
  if (gbh.random_access_point)
    _refPointCloud.clear();
#endif

  EntropyDecoder arithmeticDecoder;
  arithmeticDecoder.setBuffer(int(buf.size()) - gbhSize, buf.data() + gbhSize);
  arithmeticDecoder.start();

  if (_gps->geom_codec_type == GeometryCodecType::kOctree) {
    _currentPointCloud.resize(gbh.geom_num_points);
    decodeGeometryOctree(
#if INTER_HIERARCHICAL
      *_gps, gbh, _currentPointCloud, _refPointCloud, _backRefPointCloud, &arithmeticDecoder);
#else
      *_gps, gbh, _currentPointCloud, _refPointCloud, &arithmeticDecoder);
#endif
  }
  if (_gps->geom_codec_type == GeometryCodecType::kTriSoup) {
    decodeGeometryTrisoup(*_gps, gbh, _currentPointCloud, &arithmeticDecoder);
  }

  arithmeticDecoder.stop();
  clock_user.stop();

  auto total_user =
    std::chrono::duration_cast<std::chrono::milliseconds>(clock_user.count());
  std::cout << "positions processing time (user): "
            << total_user.count() / 1000.0 << " s\n";
  std::cout << std::endl;

  return 0;
}

//--------------------------------------------------------------------------
// Initialise the point cloud storage and decode a single geometry brick.
//
// NB: there is an implicit assumption that there is only a single geometry
//     brick per point cloud.

void
PCCTMC3Decoder3::decodeAttributeBrick(const PayloadBuffer& buf)
{
  assert(buf.type == PayloadType::kAttributeBrick);
  // todo(df): replace assertions with error handling
  assert(_sps);
  assert(_gps);

  // todo(df): validate that sps activation is not changed via the APS
  AttributeBrickHeader abh = parseAbh(buf, nullptr);
  const auto it_attr_aps = _apss.find(abh.attr_attr_parameter_set_id);

  assert(it_attr_aps != _apss.cend());
  const auto& attr_aps = it_attr_aps->second;

  assert(abh.attr_sps_attr_idx < _sps->attributeSets.size());
  const auto& attr_sps = _sps->attributeSets[abh.attr_sps_attr_idx];
  const auto& label = attr_sps.attributeLabel;

  AttributeDecoder attrDecoder;
  pcc::chrono::Stopwatch<pcc::chrono::utime_inc_children_clock> clock_user;

  clock_user.start();
  attrDecoder.decode(attr_sps, attr_aps, buf, _currentPointCloud);
  clock_user.stop();

  std::cout << label << "s bitstream size " << buf.size() << " B\n";

  auto total_user =
    std::chrono::duration_cast<std::chrono::milliseconds>(clock_user.count());
  std::cout << label
            << "s processing time (user): " << total_user.count() / 1000.0
            << " s\n";
  std::cout << std::endl;
}

//==========================================================================

void
PCCTMC3Decoder3::inverseQuantization(PCCPointSet3& pointCloud)
{
  const size_t pointCount = pointCloud.getPointCount();
  const double invScale = 1.0 / _sps->seq_source_geom_scale_factor;

  for (size_t i = 0; i < pointCount; ++i) {
    auto& point = pointCloud[i];
    for (size_t k = 0; k < 3; ++k) {
      point[k] = point[k] * invScale + minPositions[k];
    }
  }
}

//============================================================================

}  // namespace pcc

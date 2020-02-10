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

#include "PCCTMC3Encoder.h"

#include <cassert>
#include <set>

#include "AttributeEncoder.h"
#include "PCCPointSetProcessing.h"
#include "geometry.h"
#include "io_hls.h"
#include "motionWip.h"
#include "osspecific.h"
#include "pcc_chrono.h"

namespace pcc {

//============================================================================

void
PCCTMC3Encoder3::init()
{
  minPositions = 0.0;
  boundingBox.min = uint32_t(0);
  boundingBox.max = uint32_t(0);
  pointCloud.clear();
}

//----------------------------------------------------------------------------

void
PCCTMC3Encoder3::configure(EncoderParams* params)
{
  // NB: configuration updates are not supported
  // todo(df): support configuration updates

  // configure SPS:
  //  - todo(df): allow external means to force sps id
  //  - set sps id
  //  - set profile and level based on configured options
  auto sps_inserted = _spss.insert({0, params->sps});
  if (!sps_inserted.second) {
    // todo(df): error sps already set
    return;
  }

  SequenceParameterSet& sps = sps_inserted.first->second;
  sps.sps_seq_parameter_set_id = 0;
  sps.profileCompatibility.profile_compatibility_flags = 0;
  sps.level = 0;
  params->sps = sps;

  // configure GPS:
  auto gps_inserted = _gpss.insert({0, params->gps});
  GeometryParameterSet& gps = gps_inserted.first->second;
  gps.gps_seq_parameter_set_id = sps.sps_seq_parameter_set_id;
  gps.gps_geom_parameter_set_id = 0;
  // the encoder writes out the minPositions in the GBH:
  gps.geom_box_present_flag = true;

  gps.motion.motion_max_prefix_bits =
    deriveMotionMaxPrefixBits(params->gps.motion);
  gps.motion.motion_max_suffix_bits =
    deriveMotionMaxSuffixBits(params->gps.motion);

  params->gps = gps;

  // configure APS:
  int cur_apsid = -1;
  for (auto& src_aps : params->aps) {
    cur_apsid++;
    auto aps_inserted = _apss.insert({cur_apsid, src_aps});
    AttributeParameterSet& aps = aps_inserted.first->second;
    aps.aps_seq_parameter_set_id = sps.sps_seq_parameter_set_id;
    aps.aps_attr_parameter_set_id = cur_apsid;
    /*params->*/src_aps = aps;
  }

  // configure the list of attributes to code and which APS to use for each
  //  ie, maps attr_idx -> apsid, for attr_idx in sps.
  assert(
    params->attributeSetIdxToApsIdxMap.size() == sps.attributeSets.size());

  _attrIdxToApsIdMap.clear();
  _attrIdxToApsIdMap.resize(sps.attributeSets.size());
  for (int attrIdx = 0; attrIdx < sps.attributeSets.size(); attrIdx++) {
    int paramsApsIdx = params->attributeSetIdxToApsIdxMap[attrIdx];
    int paramsApsId = params->aps[paramsApsIdx].aps_attr_parameter_set_id;
    _attrIdxToApsIdMap[attrIdx] = paramsApsId;
  }
}

//----------------------------------------------------------------------------

int
PCCTMC3Encoder3::compress(
  const PCCPointSet3& inputPointCloud,
  const EncoderParams& params,
  PCCTMC3Encoder3::Callbacks* callback,
  PCCPointSet3* reconstructedCloud)
{
  // Check that configure() has been called.
  if (_spss.empty() || _gpss.empty())
    return 1;

  init();

  if (params->randomAccessPoint) {
    // ensure that there is nothing to predict the random access point with
    predPointCloud.clear();
  }

  // placeholder to "activate" the parameter sets
  _sps = &_spss[0];
  _gps = &_gpss[0];

  // write out all parameter sets prior to encoding
  callback->onOutputBuffer(write(*_sps));
  callback->onOutputBuffer(write(*_gps));
  for (const auto aps : _apss) {
    callback->onOutputBuffer(write(aps.second));
  }

  // geometry compression consists of the following stages:
  //  - prefilter/quantize geometry (non-normative)
  //  - encode geometry
  //  - recolour

  // The quantization process will determine the bounding box
  quantization(inputPointCloud);

  // geometry encoding
  if (1) {
    PayloadBuffer payload(PayloadType::kGeometryBrick);

    pcc::chrono::Stopwatch<pcc::chrono::utime_inc_children_clock> clock_user;
    clock_user.start();

    encodeGeometryBrick(params, &payload);

    clock_user.stop();

    double bpp = double(8 * payload.size()) / inputPointCloud.getPointCount();
    std::cout << "positions bitstream size " << payload.size() << " B (" << bpp
              << " bpp)\n";

    auto total_user = std::chrono::duration_cast<std::chrono::milliseconds>(
      clock_user.count());
    std::cout << "positions processing time (user): "
              << total_user.count() / 1000.0 << " s" << std::endl;

    callback->onOutputBuffer(payload);
  }

  // recolouring

  // NB: recolouring is required if points are added / removed
  bool recolourNeeded = _gps->geom_unique_points_flag
    || _gps->geom_codec_type == GeometryCodecType::kTriSoup;

  if (recolourNeeded) {
    recolour(
      inputPointCloud, _sps->seq_source_geom_scale_factor, minPositions,
      &pointCloud);
  }

  // dump recoloured point cloud
  callback->onPostRecolour(pointCloud);

  // attributeCoding

  // for each attribute
  for (int attrIdx = 0; attrIdx < _sps->attributeSets.size(); attrIdx++) {
    const auto& attr_sps = _sps->attributeSets[attrIdx];
    const auto& attr_aps = _apss[_attrIdxToApsIdMap[attrIdx]];
    const auto& label = attr_sps.attributeLabel;

    PayloadBuffer payload(PayloadType::kAttributeBrick);

    pcc::chrono::Stopwatch<pcc::chrono::utime_inc_children_clock> clock_user;
    clock_user.start();

    // todo(df): move elsewhere?
    AttributeBrickHeader abh;
    abh.attr_attr_parameter_set_id = attr_aps.aps_attr_parameter_set_id;
    abh.attr_sps_attr_idx = attrIdx;
    abh.attr_geom_brick_id = 0;
    write(abh, &payload);

    AttributeEncoder attrEncoder;
    attrEncoder.encode(attr_sps, attr_aps, pointCloud, &payload);
    clock_user.stop();

    int coded_size = int(payload.size());
    double bpp = double(8 * coded_size) / inputPointCloud.getPointCount();
    std::cout << label << "s bitstream size " << coded_size << " B (" << bpp
              << " bpp)\n";

    auto time_user = std::chrono::duration_cast<std::chrono::milliseconds>(
      clock_user.count());
    std::cout << label
              << "s processing time (user): " << time_user.count() / 1000.0
              << " s" << std::endl;

    callback->onOutputBuffer(payload);
  }

  // Save the the reconstructed (but not rescaled) point cloud as a future
  // predictor.
  predPointCloud = pointCloud;

  reconstructedPointCloud(reconstructedCloud);

  return 0;
}

//----------------------------------------------------------------------------

void
PCCTMC3Encoder3::encodeGeometryBrick(
  const EncoderParams& params, PayloadBuffer* buf)
{
  // todo(df): confirm minimum of 1 isn't needed
  uint32_t maxBB =
    std::max({1u, boundingBox.max[0], boundingBox.max[1], boundingBox.max[2]});

  // the current node dimension (log2) encompasing maxBB
  int nodeSizeLog2 = ceillog2(maxBB + 1);

  // todo(df): consider removal of trisoup_depth
  if (_gps->geom_codec_type == GeometryCodecType::kTriSoup)
    nodeSizeLog2 = _gps->trisoup_depth;

  GeometryBrickHeader gbh;
  gbh.random_access_point = !predPointCloud.getPointCount();
  gbh.geom_geom_parameter_set_id = _gps->gps_geom_parameter_set_id;
  gbh.geomBoxOrigin.x() = int(minPositions.x());
  gbh.geomBoxOrigin.y() = int(minPositions.y());
  gbh.geomBoxOrigin.z() = int(minPositions.z());
  gbh.geom_box_log2_scale = 0;
  gbh.geom_max_node_size_log2 = nodeSizeLog2;
  gbh.geom_num_points = int(pointCloud.getPointCount());
  write(*_gps, gbh, buf);

  // todo(df): remove estimate when arithmetic codec is replaced
  int maxAcBufLen = int(pointCloud.getPointCount()) * 3 * 4 + 1024;
  EntropyEncoder arithmeticEncoder(maxAcBufLen, nullptr);
  arithmeticEncoder.start();

  if (_gps->geom_codec_type == GeometryCodecType::kOctree) {
    encodeGeometryOctree(
      params, *_sps, *_gps, gbh, pointCloud, predPointCloud,
      &arithmeticEncoder);
  }
  if (_gps->geom_codec_type == GeometryCodecType::kTriSoup) {
    encodeGeometryTrisoup(params, *_sps, *_gps, gbh, pointCloud,
      &arithmeticEncoder);
  }

  uint32_t dataLen = arithmeticEncoder.stop();
  std::copy_n(arithmeticEncoder.buffer(), dataLen, std::back_inserter(*buf));
}

//----------------------------------------------------------------------------

void
PCCTMC3Encoder3::reconstructedPointCloud(PCCPointSet3* reconstructedCloud)
{
  if (reconstructedCloud == nullptr) {
    return;
  }
  const size_t pointCount = pointCloud.getPointCount();

  reconstructedCloud->addRemoveAttributes(
    pointCloud.hasColors(), pointCloud.hasReflectances());
  reconstructedCloud->resize(pointCount);

  const double minPositionQuantizationScale = 0.0000000001;
  const double invScale =
    fabs(_sps->seq_source_geom_scale_factor) > minPositionQuantizationScale
    ? 1.0 / _sps->seq_source_geom_scale_factor
    : 1.0;
  for (size_t i = 0; i < pointCount; ++i) {
    const auto quantizedPoint = pointCloud[i];
    auto& point = (*reconstructedCloud)[i];
    for (size_t k = 0; k < 3; ++k) {
      point[k] = quantizedPoint[k] * invScale + minPositions[k];
    }
    if (pointCloud.hasColors()) {
      reconstructedCloud->setColor(i, pointCloud.getColor(i));
    }
    if (pointCloud.hasReflectances()) {
      reconstructedCloud->setReflectance(i, pointCloud.getReflectance(i));
    }
  }
}

//----------------------------------------------------------------------------

void
PCCTMC3Encoder3::computeMinPositions(const PCCPointSet3& inputPointCloud)
{
  const size_t inputPointCount = inputPointCloud.getPointCount();
  minPositions = inputPointCloud[0];
  for (size_t i = 1; i < inputPointCount; ++i) {
    const auto point = inputPointCloud[i];
    for (int k = 0; k < 3; ++k) {
      if (minPositions[k] > point[k]) {
        minPositions[k] = point[k];
      }
    }
  }
}

//----------------------------------------------------------------------------

void
PCCTMC3Encoder3::quantization(const PCCPointSet3& inputPointCloud)
{
  // if sps sequence width/height/depth is set, don't auto compute bbox
  bool computeBBox = _sps->seq_bounding_box_whd == PCCVector3<int>{0};
  if (computeBBox)
    computeMinPositions(inputPointCloud);
  else {
    for (int k = 0; k < 3; k++)
      minPositions[k] = _sps->seq_bounding_box_xyz0[k];
  }

  // Clamp all points to [clampBox.min, clampBox.max] after translation
  // and quantisation.
  PCCBox3<int32_t> clampBox{{0, 0, 0}, {INT32_MAX, INT32_MAX, INT32_MAX}};
  if (!computeBBox) {
    // todo(df): this is icky (not to mention rounding issues)
    // NB: the sps seq_bounding_box_* uses unscaled co-ordinates => convert
    // NB: minus 1 to convert to max x/y/z position
    clampBox = PCCBox3<int32_t> {{0, 0, 0}, _sps->seq_bounding_box_whd};
    for (int k = 0; k < 3; k++)
      clampBox.max[k] =
        int(ceil(clampBox.max[k] * _sps->seq_source_geom_scale_factor)) - 1;
  }

  if (_gps->geom_unique_points_flag) {
    quantizePositionsUniq(
      _sps->seq_source_geom_scale_factor, -minPositions, clampBox,
      inputPointCloud, &pointCloud);
  } else {
    quantizePositions(
      _sps->seq_source_geom_scale_factor, -minPositions, clampBox,
      inputPointCloud, &pointCloud);
  }

  if (!computeBBox) {
    boundingBox.min = uint32_t(0);
    for (int k = 0; k < 3; k++)
      boundingBox.max[k] = clampBox.max[k];
    return;
  }

  const size_t pointCount = pointCloud.getPointCount();
  boundingBox.min = uint32_t(0);
  boundingBox.max = uint32_t(0);
  for (size_t i = 0; i < pointCount; ++i) {
    const PCCVector3D point = pointCloud[i];
    for (int k = 0; k < 3; ++k) {
      const uint32_t coord = uint32_t(point[k]);
      if (boundingBox.max[k] < coord) {
        boundingBox.max[k] = coord;
      }
    }
  }
}

//============================================================================

}  // namespace pcc

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

#include <vector>

#include "PCCPointSet.h"
#include "PCCTMC3Encoder.h"
#include "entropy.h"
#include "hls.h"

namespace pcc {

//============================================================================

struct PCCOctree3Node;

//============================================================================

struct PUtree {
  std::vector<bool> popul_flags;
  std::vector<bool> split_flags;
  std::vector<PCCVector3<int>> MVs;

  bool isWorld = false;
};

//============================================================================
int plus1log2shifted4(int x);

int deriveMotionMaxPrefixBits(const GeometryParameterSet::Motion& param);
int deriveMotionMaxSuffixBits(const GeometryParameterSet::Motion& param);

bool motionSearchForNode(
  const PCCPointSet3& pointCloud,
  PCCOctree3Node* node0,
  const McParam& encParam,
  const GeometryParameterSet::Motion& param,
  int nodeSizeLog2,
  EntropyEncoder* arithmeticEncoder,
  PCCPointSet3& pointPredictor,
  int* bufferPoints,
  PUtree* local_PU_tree,
  PCCPointSet3& pointPredictorVehicle,
  PCCPointSet3& pointPredictorWorld);

void encode_splitPU_MV_MC(
  PCCOctree3Node* node0,
  PUtree* local_PU_tree,
  const GeometryParameterSet::Motion& param,
  int nodeSizeLog2,
  EntropyEncoder* arithmeticEncoder,
  PCCPointSet3& pointPredictor,
  PCCPointSet3* compensatedPointCloud,
  PCCPointSet3& pointPredictorVehicle,
  PCCPointSet3& pointPredictorWorld);

void extracPUsubtree(
  const GeometryParameterSet::Motion& param,
  PUtree* local_PU_tree,
  int block_size,
  int& pos_fs,
  int& pos_fp,
  int& pos_MV,
  PUtree* destination_tree);

// motion decoder

void decode_splitPU_MV_MC(
  PCCOctree3Node* node0,
  const GeometryParameterSet::Motion& param,
  int nodeSizeLog2,
  EntropyDecoder* arithmeticDecoder,
  PCCPointSet3& pointPredictor,
  PCCPointSet3* compensatedPointCloud,
  PCCPointSet3& pointPredictorVehicle,
  PCCPointSet3& pointPredictorWorld);

// for global motion
int
distortionAB(std::vector<PCCVector3<int>>& A, std::vector<PCCVector3<int>>& B);

void quantizeGlobalMotion(double Mat_GM[4][3], int32_t Mat_GM_Q[4][3]);
void applyGlobalMotion(
  std::vector<PCCVector3<int>>& listPoints, double Mat_GM[4][3]);
void applyGlobalMotion(
  PCCPointSet3& PC,
  int32_t Mat_GM_Q[4][3],
  PCCVector3<double> vehicle_position);

double map_reference(
  std::vector<PCCVector3<int>>& pc_world_target,
  std::vector<PCCVector3<int>>& pointPredictor_centered,
  std::vector<PCCVector3<int>>& pc_world_ref);
void LMS3D(
  std::vector<PCCVector3<int>>& P1,
  std::vector<PCCVector3<int>>& P2,
  std::vector<PCCVector3<int>>& pointPredictor_centered,
  uint32_t maxBB,
  double Mat_GM[4][3]);
void SearchGlobalMotion(
  PCCPointSet3& pointCloud,
  PCCPointSet3& pointPredictor,
  PCCPointSet3& pointPredictorWorld,
  double QS,
  const GeometryParameterSet::Motion& param,
  uint32_t maxBB,
  PCCVector3D minPositions,
  int32_t Mat_GM_Q[4][3]);
void SearchGlobalMotionPerTile(
  PCCPointSet3& pointCloud,
  PCCPointSet3& pointPredictor,
  PCCPointSet3& pointPredictorWorld,
  double QS,
  const GeometryParameterSet::Motion& param,
  uint32_t maxBB,
  PCCVector3D minPositions,
  int32_t Mat_GM_Q[4][3]);

//============================================================================

}  // namespace pcc

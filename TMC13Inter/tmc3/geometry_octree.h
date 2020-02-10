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

#include <cstdint>
#include <memory>

#include "PCCMath.h"
#include "PCCPointSet.h"
#include "entropy.h"
#include "hls.h"
#include "motionWip.h"
#include "ringbuf.h"
#include "tables.h"

namespace pcc {

//============================================================================

const int MAX_NUM_DM_LEAF_POINTS = 2;

//============================================================================

struct PCCOctree3Node {
  // 3D position of the current node's origin (local x,y,z = 0).
  PCCVector3<uint32_t> pos;

  // Range of point indexes spanned by node
  uint32_t start;
  uint32_t end;

  // Range of prediction's point indexes spanned by node
  uint32_t predStart;
  uint32_t predEnd;

  // address of the current node in 3D morton order.
  uint64_t mortonIdx;

  // pattern denoting occupied neighbour nodes.
  //    32 8 (y)
  //     |/
  //  2--n--1 (x)
  //    /|
  //   4 16 (z)
  uint8_t neighPattern = 0;

  // The current node's number of siblings plus one.
  // ie, the number of child nodes present in this node's parent.
  uint8_t numSiblingsPlus1;

  // The number of mispredictions in determining the occupancy
  // map of the child nodes in this node's parent.
  int8_t numSiblingsMispredicted;

  // The occupancy map used describing the current node and its siblings.
  uint8_t siblingOccupancy;

  // Prediction Unit tree (encoder only) attached to node
  std::unique_ptr<PUtree> PU_tree;

  // prediction ranges refer to compensated reference
  bool isCompensated : 1;
  bool hasMotion : 1;
  bool isWorld : 1;
};

//---------------------------------------------------------------------------

uint8_t mapGeometryOccupancy(uint8_t occupancy, uint8_t neighPattern);
uint8_t mapGeometryOccupancyInv(uint8_t occupancy, uint8_t neighPattern);

void updateGeometryNeighState(
  bool siblingRestriction,
  const ringbuf<PCCOctree3Node>::iterator& bufEnd,
  int64_t numNodesNextLvl,
  int childSizeLog2,
  PCCOctree3Node& child,
  int childIdx,
  uint8_t neighPattern,
  uint8_t parentOccupancy);

//---------------------------------------------------------------------------
// :: octree encoder exposing internal ringbuffer

void encodeGeometryOctree(
  const EncoderParams& params,
  const SequenceParameterSet& sps,
  const GeometryParameterSet& gps,
  const GeometryBrickHeader& gbh,
  PCCPointSet3& pointCloud,
  PCCPointSet3& predPointCloud,
  EntropyEncoder* arithmeticEncoder,
  pcc::ringbuf<PCCOctree3Node>* nodesRemaining);

void decodeGeometryOctree(
  const GeometryParameterSet& gps,
  const GeometryBrickHeader& gbh,
  PCCPointSet3& pointCloud,
  PCCPointSet3& predPointCloud,
  EntropyDecoder* arithmeticDecoder,
  pcc::ringbuf<PCCOctree3Node>* nodesRemaining);

//---------------------------------------------------------------------------
// Determine if direct coding is permitted.
// If tool is enabled:
//   - Block must not be near the bottom of the tree
//   - The parent / grandparent are sparsely occupied
//   - Occupancy is not predictable

inline bool
isDirectModeEligible(
  bool featureEnabled,
  int nodeSizeLog2,
  const PCCOctree3Node& node,
  const PCCOctree3Node& child,
  bool occupancyIsPredictable)
{
  return featureEnabled && (nodeSizeLog2 >= 2) && (node.neighPattern == 0)
    && (child.numSiblingsPlus1 == 1) && (node.numSiblingsPlus1 <= 2)
    && !occupancyIsPredictable;
}

//---------------------------------------------------------------------------
// Select the neighbour pattern reduction table according to GPS config.

using NeighPattern64toX = const uint8_t[64];

inline NeighPattern64toX&
neighPattern64toR1(const GeometryParameterSet& gps)
{
  if (!gps.neighbour_context_restriction_flag)
    return kNeighPattern64to10;
  return kNeighPattern64to6;
}

//---------------------------------------------------------------------------

struct CtxModelOctreeOccupancy {
  AdaptiveBitModelFast contexts[256];
  int ctxFactorShift;

  CtxModelOctreeOccupancy(int ctxFactorShift) : ctxFactorShift(ctxFactorShift)
  {}

  AdaptiveBitModelFast& operator[](int idx)
  {
    return contexts[idx >> ctxFactorShift];
  }
};

//---------------------------------------------------------------------------
// Encapsulates the derivation of ctxIdx for occupancy coding.

class CtxMapOctreeOccupancy {
public:
  struct CtxIdxMap {
    uint8_t b0[10];
    uint8_t b1[20];
    uint8_t b2[39];
    uint8_t b3[76];
    uint8_t b4[149];
    uint8_t b5[294];
    uint8_t b6[391];
    uint8_t b7[520];
  };

  CtxMapOctreeOccupancy();

  const uint8_t* operator[](int bit) const { return b[bit]; }

  uint8_t* operator[](int bit) { return b[bit]; }

  // return *ctxIdx and update *ctxIdx according to bit
  static uint8_t evolve(bool bit, uint8_t* ctxIdx);

private:
  std::unique_ptr<CtxIdxMap> map;
  uint8_t* b[8];
};

//----------------------------------------------------------------------------

inline uint8_t
CtxMapOctreeOccupancy::evolve(bool bit, uint8_t* ctxIdx)
{
  uint8_t retval = *ctxIdx;

  if (bit)
    *ctxIdx = kCtxMapOctreeOccupancyEvolutionOn1[*ctxIdx];
  else
    *ctxIdx = kCtxMapOctreeOccupancyEvolutionOn0[*ctxIdx];

  return retval;
}

//---------------------------------------------------------------------------

}  // namespace pcc

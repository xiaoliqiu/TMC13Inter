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

namespace pcc {

// Symmetry reduction of 64 neighbour pattern to 10
extern const uint8_t kNeighPattern64to10[64];
extern const uint8_t kNeighPattern64to6[64];
extern const uint8_t kNeighPattern10to7[10];
extern const uint8_t kNeighPattern7to5[7];

// Identifies the X-Y rotation to be performed given a neighbour pattern
//  0 => none
//  1 => X  90° ACW (encoder)  [kOccMapRotateX090]
//  1 => X 270° ACW (decoder)  [kOccMapRotateX270]
//  2 => X 270° ACW, Y 180°    [kOccMapRotateX270Y180]
//  3 => X  90° ACW, Y 180°    [kOccMapRotateX090Y180]
extern const uint8_t kOccMapRotateXIdFromPattern[64];
extern const uint8_t kOccMapRotateX270[256];
extern const uint8_t LUT_rotx_mode2[256];
extern const uint8_t LUT_rotx_mode3[256];
extern const uint8_t kOccMapRotateX090[256];
extern const uint8_t kOccMapRotateX270Y180[256];
extern const uint8_t kOccMapRotateX090Y180[256];

// Identifies the Y rotation to be performed given a neighbour pattern
//  0 => none
//  1 => Y 270° ACW (encoder)  [kOccMapRotateY270]
//  1 => Y  90° ACW (decoder)  [kOccMapRotateY090]
extern const bool kOccMapRotateYIdFromPattern[64];
extern const uint8_t kOccMapRotateY090[256];
extern const uint8_t kOccMapRotateY270[256];

// Vertical mirroring permutation
extern const uint8_t kOccMapMirrorXY[256];

// Identifies the Z rotation to be performed given a partial neighbour pattern
//  0 => none
//  1 => 270° ACW (decoder)  [kOccMapRotateZ270]
//  1 => 090° ACW (encoder)  [kOccMapRotateZ090]
//  2 => 180°                [kOccMapRotateZ180]
//  3 => 090° ACW (decoder)  [kOccMapRotateZ090]
//  3 => 270° ACW (encoder)  [kOccMapRotateZ270]
extern const uint8_t kOccMapRotateZIdFromPatternXY[16];
extern const uint8_t kOccMapRotateZ270[256];
extern const uint8_t kOccMapRotateZ180[256];
extern const uint8_t kOccMapRotateZ090[256];

// todo(df): ...
extern const uint8_t kCtxMapOctreeOccupancyEvolutionOn0[256];
extern const uint8_t kCtxMapOctreeOccupancyEvolutionOn1[256];

// LUT initialisation table
extern const uint8_t kDualLutOccupancyCoderInit[10][32];

//============================================================================
// Mapping of (x,y,z) components to 3D Morton code.
extern const uint32_t kMortonCode256X[256];
extern const uint32_t kMortonCode256Y[256];
extern const uint32_t kMortonCode256Z[256];

} /* namespace pcc */

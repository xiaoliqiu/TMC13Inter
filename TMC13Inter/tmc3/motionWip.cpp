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

#include "motionWip.h"

#include <algorithm>
#include <cfloat>
#include <climits>
#include <set>
#include <vector>

#include "PCCMath.h"
#include "PCCPointSet.h"
#include "entropy.h"
#include "geometry_octree.h"

namespace pcc {

//============================================================================

struct MotionEntropy {
  AdaptiveBitModel splitPu;
#if MV_PREDICTION_RDO
  AdaptiveBitModel bPred;
#endif

  StaticBitModel mvSign;
  AdaptiveBitModel mvIsZero;
  StaticBitModel expGolombV0;
  AdaptiveBitModel expGolombV[6];

  MotionEntropy();
};

//----------------------------------------------------------------------------

class MotionEntropyEncoder : public MotionEntropy {
public:
  MotionEntropyEncoder(EntropyEncoder* arithmeticEncoder)
    : _arithmeticEncoder(arithmeticEncoder)
  {}

  void encodeSplitPu(int symbol);
#if MV_PREDICTION_RDO
  void encodePredPu(int symbol);
#endif
  void encodeVector(
    const PCCVector3<int>& mv,
    int mvPrecision,
    int boundPrefix,
    int boundSuffix);

private:
  EntropyEncoder* _arithmeticEncoder;
};

//----------------------------------------------------------------------------

class MotionEntropyDecoder : public MotionEntropy {
public:
  MotionEntropyDecoder(EntropyDecoder* arithmeticDecoder)
    : _arithmeticDecoder(arithmeticDecoder)
  {}

  bool decodeSplitPu();
#if MV_PREDICTION_RDO
  bool decodePredPu();
#endif
  void decodeVector(PCCVector3<int>* mv, int boundPrefix, int boundSuffix);

private:
  EntropyDecoder* _arithmeticDecoder;
};

//----------------------------------------------------------------------------

struct MotionEntropyEstimate {
  double hMvIsZero[2];
  double hExpGolombV[6][2];
  double hSplitPu[2];

  MotionEntropyEstimate(const MotionEntropy& codec);

  double estimateVector(
    const PCCVector3<int>& mv,
    int mvPrecision,
    int boundPrefix,
    int boundSuffix) const;
};

//============================================================================

MotionEntropy::MotionEntropy()
{
  splitPu.reset(true);
  mvIsZero.reset(true);
  for (int i = 0; i < 6; i++)
    expGolombV[i].reset(true);
}

//----------------------------------------------------------------------------

MotionEntropyEstimate::MotionEntropyEstimate(const MotionEntropy& codec)
{
  codec.splitPu.getEntropy(hSplitPu);
  codec.mvIsZero.getEntropy(hMvIsZero);
  for (int i = 0; i < 6; i++)
    codec.expGolombV[i].getEntropy(hExpGolombV[i]);
}

//----------------------------------------------------------------------------

inline void
MotionEntropyEncoder::encodeSplitPu(int symbol)
{
  _arithmeticEncoder->encode(symbol, splitPu);
}

#if MV_PREDICTION_RDO
inline void
MotionEntropyEncoder::encodePredPu(int symbol)
{
  _arithmeticEncoder->encode(symbol, bPred);
}
#endif

//----------------------------------------------------------------------------

inline bool
MotionEntropyDecoder::decodeSplitPu()
{
  return _arithmeticDecoder->decode(splitPu);
}

#if MV_PREDICTION_RDO
inline bool
MotionEntropyDecoder::decodePredPu()
{
  return _arithmeticDecoder->decode(bPred);
}
#endif

//----------------------------------------------------------------------------

inline void
MotionEntropyEncoder::encodeVector(
  const PCCVector3<int>& mv, int mvPrecision, int boundPrefix, int boundSuffix)
{
  for (int comp = 0; comp < 3; comp++) {
    int v = mv[comp];
    if (v == 0) {
      _arithmeticEncoder->encode(1, mvIsZero);
    } else {
      _arithmeticEncoder->encode(0, mvIsZero);

      _arithmeticEncoder->encode(v < 0, mvSign);
      if (v < 0)
        v = -v;
      v /= mvPrecision;
      v--;
      // expGolomb on |v|-1 with truncation
      _arithmeticEncoder->encodeExpGolomb(
        uint32_t(v), 1, expGolombV0, expGolombV, boundPrefix, boundSuffix);
    }
  }
}

//----------------------------------------------------------------------------

inline void
MotionEntropyDecoder::decodeVector(
  PCCVector3<int>* mv, int boundPrefix, int boundSuffix)
{
  for (int comp = 0; comp < 3; comp++) {
    if (_arithmeticDecoder->decode(mvIsZero)) {
      (*mv)[comp] = 0;
      continue;
    }
    bool sign = _arithmeticDecoder->decode(mvSign);
    int v = 1
      + _arithmeticDecoder->decodeExpGolomb(
          1, expGolombV0, expGolombV, boundPrefix, boundSuffix);
    if (sign)
      v = -v;
    (*mv)[comp] = v;
  }
}

//----------------------------------------------------------------------------

double
MotionEntropyEstimate::estimateVector(
  const PCCVector3<int>& mv,
  int mvPrecision,
  int boundPrefix,
  int boundSuffix) const
{
  double r = 0.;
  for (int comp = 0; comp < 3; comp++) {
    int v = mv[comp];
    if (!v) {
      r += hMvIsZero[1];
    } else {
      r += hMvIsZero[0];
      r += 1.;  // unpredictable sign

      if (v < 0)
        v = -v;
      v /= mvPrecision;
      v--;

      int k = 1;  // initially, expgolomb order
      int ctx_number = 0;
      int num_bit_prefix = 0;
      while (1) {
        if (boundPrefix && v >= static_cast<unsigned int>(1 << k)) {
          // unary part
          r += hExpGolombV[ctx_number][1];
          v = v - (1 << k);
          k++;

          ctx_number++;
          if (ctx_number > 5)
            ctx_number = 5;
          num_bit_prefix++;
          continue;
        }

        // terminated zero of unary part + suffix
        if (num_bit_prefix < boundPrefix)
          r += hExpGolombV[ctx_number][0];
        if (num_bit_prefix == boundPrefix)
          k = boundSuffix;
        r += k;
        break;
      }
    }
  }
  return r;
}

//============================================================================

int
deriveMotionMaxPrefixBits(const GeometryParameterSet::Motion& param)
{
  if (param.motion_precision == 0)
    return 0;

  int max_MV =
    std::max(0, param.motion_window_size / param.motion_precision - 1);
  if (max_MV >= 256)
    return 31;

  static int LUT_bound_prefix[256] = {
    0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7};

  return LUT_bound_prefix[max_MV];
}

int
deriveMotionMaxSuffixBits(const GeometryParameterSet::Motion& param)
{
  if (param.motion_precision == 0)
    return 0;

  int max_MV =
    std::max(0, param.motion_window_size / param.motion_precision - 1);
  if (max_MV >= 256)
    return 31;

  static int LUT_bound_suffix[256] = {
    0, 1, 0, 1, 2, 2, 0, 1, 2, 2, 3, 3, 3, 3, 0, 1, 2, 2, 3, 3, 3, 3, 4, 4,
    4, 4, 4, 4, 4, 4, 0, 1, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 0, 1, 2, 2, 3, 3, 3, 3, 4, 4,
    4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6,
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
    6, 6, 6, 6, 6, 6, 0, 1, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 0, 1};

  return LUT_bound_suffix[max_MV];
}

//============================================================================

static double
find_motion(
  const McParam& encParam,
  const GeometryParameterSet::Motion& param,
  const MotionEntropyEstimate& motionEntropy,
  std::vector<PCCVector3<int>>& Window,
  std::vector<PCCVector3<int>>& Block0,
  int x0,
  int y0,
  int z0,
  int local_size,
  PUtree* local_PU_tree,
  int* bufferPoints)
{
  if (!Window.size())
    return DBL_MAX;

  int Amotion = encParam.Amotion0;

  // ---------------------------- test no split --------------------
  // test V=0 and  dynamic planing at Block0-Window +/- window_size
  double cost_NoSplit = DBL_MAX;
  int wSize = param.motion_window_size;

  if (Window.size() > encParam.sampleW * std::max(int(Block0.size()), 16))
    wSize >>= 1;

  int max_distance = wSize << 1;
  PCCVector3<int> bestV_NoSplit = PCCVector3<int>(0, 0, 0);
  PCCVector3<int> bestV_NoSplit_2nd = PCCVector3<int>(0, 0, 0);
  bool flag_must_split = false;  // local_size > param.minPuSize && Block0.size() >= 256;

  if (!flag_must_split) {
    double d = 0.;
    int Dist = 0;
    std::vector<int> blockEnds;
    blockEnds.push_back(0);
    int blockPos = 0;
    int* pBuffer = bufferPoints;

    std::vector<PCCVector3<int>>::iterator itB = Block0.begin();
    int jumpBlock = 1 + (Block0.size() >> encParam.decimate);  // (kind of) random sampling of the original block to code
    int Nsamples = 0;
    int a0, a1, a2;
#if INTER_HIERARCHICAL
    for (int Nb = 0; Nb < int(Block0.size()); Nb += jumpBlock) {
#else
    for (int Nb = 0; Nb < int(Block0.size()); Nb += jumpBlock, itB += jumpBlock, Nsamples++) {
#endif
      PCCVector3<int> b = *itB;
      int min_d = max_distance;
      for (const auto w : Window) {
#if MV_PREDICTION_ME
        if (w[0] >= x0 && w[0] < x0 + local_size && w[1] >= y0 &&
          w[1] < y0 + local_size && w[2] >= z0 && w[2] < z0 + local_size)
        {
          a0 = std::abs(b[0] - w[0]);
          a1 = std::abs(b[1] - w[1]);
          a2 = std::abs(b[2] - w[2]);

          a0 += a1 + a2;
          if (a0 < min_d)
            min_d = a0;
        }
      }
#else
        pBuffer[0] = b[0] - w[0];
        a0 = std::abs(pBuffer[0]);
        pBuffer[1] = b[1] - w[1];
        a1 = std::abs(pBuffer[1]);
        pBuffer[2] = b[2] - w[2];
        a2 = std::abs(pBuffer[2]);
        if (a0 <= wSize && a1 <= wSize && a2 <= wSize) {
          pBuffer += 3;
          blockPos++;
        }

        a0 += a1 + a2;  // replaces int local_d
        if (a0 < min_d)
          min_d = a0;
      }
      blockEnds.push_back(blockPos);
#endif
      Dist += plus1log2shifted4(min_d);  // 1/0.0625 = 16 times log
#if INTER_HIERARCHICAL
      Nsamples++;
      if (Nsamples * jumpBlock < Block0.size())
      {
        itB += jumpBlock;
      }
#endif
    }

    d = jumpBlock * Dist * 0.0625 + encParam.lambda * motionEntropy.estimateVector(
            PCCVector3<int>(0, 0, 0), param.motion_precision,
            param.motion_max_prefix_bits, param.motion_max_suffix_bits);

    // set loop search parameters
    double best_d[2];
    best_d[0] = d;
    best_d[1] = d;
    std::vector<PCCVector3<int>> list_tested;
    list_tested.push_back(PCCVector3<int>(0, 0, 0));

    const int searchPattern[3 * 18] = {
      1,  1, 0, 1,  0,  0,  1, -1, 0,  0, 1, 0,  0, -1, 0,  -1, 1,  0,
      -1, 0, 0, -1, -1, 0,  1, 0,  1,  0, 1, 1,  0, 0,  1,  0,  -1, 1,
      -1, 0, 1, 1,  0,  -1, 0, 1,  -1, 0, 0, -1, 0, -1, -1, -1, 0,  -1};
    //const int searchPattern[3 * 6] = {   1,0,0,   0,1,0,  0,-1,0,    -1,0,0,    0,0,1,   0,0,-1};
    //const int searchPattern[3 * 26] = { 1,1,0,  1,0,0,  1,-1,0,  0,1,0,  0,-1,0,  -1,1,0,  -1,0,0,  -1,-1,0,  1,0,1,  0,1,1,  0,0,1,  0,-1,1,  -1,0,1,  1,0,-1,  0,1,-1,  0,0,-1,  0,-1,-1,  -1,0,-1,  1,1,1,  1,1,-1,  1,-1,1,  1,-1,-1,  -1,1,1,  -1,1,-1,  -1,-1,1,  -1,-1,-1 };

    // loop MV search
    bool flag_first_search = true;  // if true one seed, otherwise two
    while (Amotion >= param.motion_precision) {
      double old_b = best_d[0];

      // loop on MV best seeds
      int Nseed = flag_first_search ? 1 : 2;
      PCCVector3<int> Vs = bestV_NoSplit;
      PCCVector3<int> old_bestV_NoSplit_2nd = bestV_NoSplit_2nd;
      for (int see = 0; see < Nseed; see++) {
        if (see == 1)
          Vs = old_bestV_NoSplit_2nd;
        const int* pSearch = searchPattern;

        // loop on searchPattern
        for (int t = 0; t < 18; t++, pSearch += 3) {
          PCCVector3<int> V = Vs + PCCVector3<int>(pSearch[0] * Amotion, pSearch[1] * Amotion, pSearch[2] * Amotion);
          if ( std::abs(V[0]) > wSize || std::abs(V[1]) > wSize || std::abs(V[2]) > wSize)
            continue;  // ensure MV does not go out of the window
          if ( std::find(list_tested.begin(), list_tested.end(), V ) != list_tested.end() )
            continue;

          Dist = 0;
          int index = 0;
          pBuffer = bufferPoints;
#if MV_PREDICTION_ME
          std::vector<PCCVector3<int>>::iterator itB = Block0.begin();
          int jumpBlock = 1 + (Block0.size() >> encParam.decimate);  // (kind of) random sampling of the original block to code
          int Nsamples = 0;
          int a0, a1, a2;
          for (int Nb = 0; Nb < int(Block0.size()); Nb += jumpBlock, itB += jumpBlock, Nsamples++) {
            PCCVector3<int> b = *itB;
            int min_d = max_distance;
            for (const auto w : Window) {
              PCCVector3<int> wV = w - V;
              if (wV[0] >= x0 && wV[0] < x0 + local_size && wV[1] >= y0 &&
                wV[1] < y0 + local_size && wV[2] >= z0 && wV[2] < z0 + local_size)
              {
                a0 = std::abs(b[0] - wV[0]);
                a1 = std::abs(b[1] - wV[1]);
                a2 = std::abs(b[2] - wV[2]);

                a0 += a1 + a2;
                if (a0 < min_d)
                  min_d = a0;
              }
            }
            Dist += plus1log2shifted4(min_d);
          }
#else
          for (int Nb = 0; Nb < Nsamples; Nb++) {
            int min_d = max_distance;
            while (index < blockEnds[Nb + 1]) {
              int local_d = std::abs(pBuffer[0] + V[0]) + std::abs(pBuffer[1] + V[1]) + std::abs(pBuffer[2] + V[2]);
              if (local_d < min_d)
                min_d = local_d;
              index++;
              pBuffer += 3;
            }
            Dist += plus1log2shifted4(min_d);  // 1/0.0625 = 16 times log
          }
#endif
          d = jumpBlock * Dist * 0.0625 + encParam.lambda * motionEntropy.estimateVector(
                  V, param.motion_precision, param.motion_max_prefix_bits,
                  param.motion_max_suffix_bits);

          // keep 2 best MV
          if (d < best_d[0]) {
            best_d[1] = best_d[0];
            bestV_NoSplit_2nd = bestV_NoSplit;
            best_d[0] = d;
            bestV_NoSplit = V;
          } else if (d < best_d[1]) {
            best_d[1] = d;
            bestV_NoSplit_2nd = V;
          }

          list_tested.push_back(V);
        }  // end loop on searchPattern
      }    // end loop on seeds

      // log reduction of search range
      if (old_b == best_d[0]) {
        Amotion >>= 1;
        if (!bestV_NoSplit[0] && !bestV_NoSplit[1] && !bestV_NoSplit[2]) {
          bool flag = Amotion > param.motion_precision;
          Amotion >>= 1;
          if (flag)
            Amotion = std::max(Amotion, param.motion_precision);
        }
      }

      flag_first_search = false;
    }  // end loop MV search

    cost_NoSplit = best_d[0];
  }

  // ---------------------------- test split --------------------
  double cost_Split = DBL_MAX;
  PUtree* Split_PU_tree = new PUtree;  // local split tree

  if (local_size > param.motion_min_pu_size && Block0.size() >= 8) {
    // condition on number of points for search acceleration
    int local_size1 = local_size >> 1;

    std::vector<PCCVector3<int>> list_xyz;
    PCCVector3<int> xyz0 = PCCVector3<int>(x0, y0, z0);

    list_xyz.push_back(PCCVector3<int>(0, 0, 0) + xyz0);
    list_xyz.push_back(PCCVector3<int>(0, 0, local_size1) + xyz0);
    list_xyz.push_back(PCCVector3<int>(0, local_size1, 0) + xyz0);
    list_xyz.push_back(PCCVector3<int>(0, local_size1, local_size1) + xyz0);
    list_xyz.push_back(PCCVector3<int>(local_size1, 0, 0) + xyz0);
    list_xyz.push_back(PCCVector3<int>(local_size1, 0, local_size1) + xyz0);
    list_xyz.push_back(PCCVector3<int>(local_size1, local_size1, 0) + xyz0);
    list_xyz.push_back(
      PCCVector3<int>(local_size1, local_size1, local_size1) + xyz0);

    // loop on 8 child PU
    cost_Split = 0.;
    for (int t = 0; t < 8; t++) {
      // child PU coordinates
      int x1 = list_xyz[t][0];
      int y1 = list_xyz[t][1];
      int z1 = list_xyz[t][2];

      // block for child PU
      std::vector<PCCVector3<int>> Block1;
      int xhigh = x1 + local_size1;
      int yhigh = y1 + local_size1;
      int zhigh = z1 + local_size1;
      for (const auto& b : Block0) {
        if (
          b[0] >= x1 && b[0] < xhigh && b[1] >= y1 && b[1] < yhigh
          && b[2] >= z1 && b[2] < zhigh)
          Block1.push_back(b);
      }

      cost_Split +=
        1.0;  // the cost due to not coding the occupancy with inter pred

      if (!Block1.size()) {  // empty PU
        Split_PU_tree->popul_flags.push_back(0);
        continue;
      }
      Split_PU_tree->popul_flags.push_back(1);

      std::vector<PCCVector3<int>> Window1;
      wSize = param.motion_window_size;
      for (const auto& b : Window) {
        if ( b[0] >= x1 - wSize && b[0] < xhigh + wSize && b[1] >= y1 - wSize
          && b[1] < yhigh + wSize && b[2] >= z1 - wSize && b[2] < zhigh + wSize )
          Window1.push_back(b);
      }
      cost_Split += find_motion(
        encParam, param, motionEntropy, Window1, Block1, x1, y1, z1,
        local_size1, Split_PU_tree, bufferPoints);
    }
  }

  // ---------------------------- choose split vs no split --------------------
  if (local_size > param.motion_min_pu_size) {
    cost_NoSplit +=
      encParam.lambda * motionEntropy.hSplitPu[0];  // cost no split flag
  }
  cost_Split +=
    encParam.lambda * motionEntropy.hSplitPu[1];  // cost split flag

  if (
    local_size <= param.motion_min_pu_size
    || cost_NoSplit <= cost_Split) {  // no split
    // push non split flag, only if size>size_min
    if (local_size > param.motion_min_pu_size) {
      local_PU_tree->split_flags.push_back(0);
    }
    // push MV
    local_PU_tree->MVs.push_back(bestV_NoSplit);

    return cost_NoSplit;
  } else {
    // split
    local_PU_tree->split_flags.push_back(1);  // push split PU flag

    // append Split_PU_tree to  local_PU_tree
    for (const auto& f : Split_PU_tree->popul_flags)
      local_PU_tree->popul_flags.push_back(f);
    for (const auto& f : Split_PU_tree->split_flags)
      local_PU_tree->split_flags.push_back(f);
    for (const auto& v : Split_PU_tree->MVs)
      local_PU_tree->MVs.push_back(v);

    delete Split_PU_tree;
    return cost_Split;
  }
}

//----------------------------------------------------------------------------
void
extracPUsubtree(
  const GeometryParameterSet::Motion& param,
  PUtree* local_PU_tree,
  int block_size,
  int& pos_fs,
  int& pos_fp,
  int& pos_MV,
  PUtree* destination_tree)
{
  // copy World vs vehicle
  destination_tree->isWorld = local_PU_tree->isWorld;

  // non-split terminal case
  if (
    block_size <= param.motion_min_pu_size
    || !local_PU_tree->split_flags[pos_fs]) {
    if (block_size > param.motion_min_pu_size) {
      destination_tree->split_flags.push_back(0);
      pos_fs++;
    }

    destination_tree->MVs.push_back(local_PU_tree->MVs[pos_MV]);
    pos_MV++;
    return;
  }

  // split case
  destination_tree->split_flags.push_back(1);
  pos_fs++;

  // loop on 8 children
  for (int s = 0; s < 8; s++) {
    if (local_PU_tree->popul_flags[pos_fp]) {  // child populated
      destination_tree->popul_flags.push_back(1);
      pos_fp++;

      extracPUsubtree(
        param, local_PU_tree, block_size >> 1, pos_fs, pos_fp, pos_MV,
        destination_tree);
    } else {  // child not pouplated
      destination_tree->popul_flags.push_back(1);
      pos_fp++;
    }
  }
}

//----------------------------------------------------------------------------
static const int LUT_LOG2[64]{
  INT_MIN, 0,  16, 25, 32, 37, 41, 45, 48, 51, 53, 55, 57, 59, 61, 63,
  64,      65, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 79,
  80,      81, 81, 82, 83, 83, 84, 85, 85, 86, 86, 87, 87, 88, 88, 89,
  89,      90, 90, 91, 91, 92, 92, 93, 93, 93, 94, 94, 95, 95, 95, 96};

int
plus1log2shifted4(int x)
{
  x++;
  int result = 0;
  while (x >= 64) {
    x >>= 1;
    result += 16;
  }

  return result + LUT_LOG2[x];
}

bool
motionSearchForNode(
  const PCCPointSet3& pointCloud,
  PCCOctree3Node* node0,
  const McParam& encParam,
  const GeometryParameterSet::Motion& param,
  const int nodeSizeLog2,
  EntropyEncoder* arithmeticEncoder,
  PCCPointSet3& pointPredictor,
#if INTER_HIERARCHICAL
  PCCPointSet3& backPointPredictor,
  int& interDir,
#endif
  int* bufferPoints,
  PUtree* local_PU_tree,
  PCCPointSet3& pointPredictorVehicle,
  PCCPointSet3& pointPredictorWorld)
{
  MotionEntropyEncoder motionEncoder(arithmeticEncoder);

  // converting the point cloud from double to int beforehand
  std::vector<PCCVector3<int>> Block0;
  for (size_t i = node0->start; i < node0->end; ++i) {
    const PCCVector3D point = pointCloud[i];
    Block0.push_back(
      PCCVector3<int>(int(point[0]), int(point[1]), int(point[2])));
  }

  // find search window while converting to integer
  std::vector<PCCVector3<int>> Window;
#if INTER_HIERARCHICAL
  std::vector<PCCVector3<int>> backWindow;
#endif
  std::vector<PCCVector3<int>> Window_V;
  std::vector<PCCVector3<int>> Window_W;
  int node_size = 1 << nodeSizeLog2;

  int xlow = node0->pos[0] - param.motion_window_size;
  int xhigh = node0->pos[0] + node_size + param.motion_window_size;
  int ylow = node0->pos[1] - param.motion_window_size;
  int yhigh = node0->pos[1] + node_size + param.motion_window_size;
  int zlow = node0->pos[2] - param.motion_window_size;
  int zhigh = node0->pos[2] + node_size + param.motion_window_size;

  // window size
  if (!param.global_motion_enabled)  // no global motion
  {
    for (size_t i = 0; i < pointPredictor.getPointCount(); ++i) {
      const PCCVector3D point = pointPredictor[i];
      const PCCVector3<int> point_int =
        PCCVector3<int>(int(point[0]), int(point[1]), int(point[2]));

      if (
        point_int[0] >= xlow && point_int[0] < xhigh && point_int[1] >= ylow
        && point_int[1] < yhigh && point_int[2] >= zlow
        && point_int[2] < zhigh)
        Window.push_back(point_int);
    }

#if INTER_HIERARCHICAL
    for (size_t i = 0; i < backPointPredictor.getPointCount(); ++i) {
      const PCCVector3D point = backPointPredictor[i];
      const PCCVector3<int> point_int =
        PCCVector3<int>(int(point[0]), int(point[1]), int(point[2]));

      if (
        point_int[0] >= xlow && point_int[0] < xhigh && point_int[1] >= ylow
        && point_int[1] < yhigh && point_int[2] >= zlow
        && point_int[2] < zhigh)
        backWindow.push_back(point_int);
    }
#endif

    // if window is empty, no compensation
#if INTER_HIERARCHICAL
    if (!Window.size() && !backWindow.size()) {
#else
    if (!Window.size()) {
#endif
      return false;
    }
  } else  //global motion
  {
    // window for vehicle referential
    for (size_t i = 0; i < pointPredictorVehicle.getPointCount(); ++i) {
      const PCCVector3D point = pointPredictorVehicle[i];
      const PCCVector3<int> point_int =
        PCCVector3<int>(int(point[0]), int(point[1]), int(point[2]));

      if (
        point_int[0] >= xlow && point_int[0] < xhigh && point_int[1] >= ylow
        && point_int[1] < yhigh && point_int[2] >= zlow
        && point_int[2] < zhigh)
        Window_V.push_back(point_int);
    }

    // window for world referential
    for (size_t i = 0; i < pointPredictorWorld.getPointCount(); ++i) {
      const PCCVector3D point = pointPredictorWorld[i];
      const PCCVector3<int> point_int =
        PCCVector3<int>(int(point[0]), int(point[1]), int(point[2]));

      if (
        point_int[0] >= xlow && point_int[0] < xhigh && point_int[1] >= ylow
        && point_int[1] < yhigh && point_int[2] >= zlow
        && point_int[2] < zhigh)
        Window_W.push_back(point_int);
    }

    // if window is empty, no compensation
    if (!Window_V.size() && !Window_W.size()) {
      return false;
    }
  }

  // entropy estimates
  MotionEntropyEstimate mcEstimate(motionEncoder);

  // motion search
  //PUtree local_PU_tree;
  int x0 = node0->pos[0];
  int y0 = node0->pos[1];
  int z0 = node0->pos[2];

  if (!param.global_motion_enabled) {
    // MV search
#if INTER_HIERARCHICAL
    if ( Window.size() && backWindow.size())
    {
      std::unique_ptr<PUtree> forward_PU_tree(new PUtree);
      std::unique_ptr<int> forwardbufferPoints;
      forwardbufferPoints.reset(new int[3 * 3000 * 10000]);
      double forwardCost = find_motion(
        encParam, param, mcEstimate, Window, Block0, x0, y0, z0,
        param.motion_block_size, forward_PU_tree.get(), forwardbufferPoints.get());

      std::unique_ptr<PUtree> backward_PU_tree(new PUtree);
      std::unique_ptr<int> backwardbufferPoints;
      backwardbufferPoints.reset(new int[3 * 3000 * 10000]);
      double backwardCost = find_motion(
        encParam, param, mcEstimate, backWindow, Block0, x0, y0, z0,
        param.motion_block_size, backward_PU_tree.get(), backwardbufferPoints.get());

      if (forwardCost < backwardCost)
      {
        *local_PU_tree = *(forward_PU_tree.get());
        *bufferPoints = *(forwardbufferPoints.get());
        interDir = 0;
      }
      else
      {
        *local_PU_tree = *(backward_PU_tree.get());
        *bufferPoints = *(backwardbufferPoints.get());
        interDir = 1;
      }
    }
    else if (Window.size())
    {
      find_motion(
        encParam, param, mcEstimate, Window, Block0, x0, y0, z0,
        param.motion_block_size, local_PU_tree, bufferPoints);
      
      interDir = 0;
    }
    else if (backWindow.size())
    {
      find_motion(
        encParam, param, mcEstimate, backWindow, Block0, x0, y0, z0,
        param.motion_block_size, local_PU_tree, bufferPoints);

      interDir = 1;
    }
#else
    find_motion(
      encParam, param, mcEstimate, Window, Block0, x0, y0, z0,
      param.motion_block_size, local_PU_tree, bufferPoints);
#endif
  } else  //global motion
  {
    if (encParam.globalMotionInRdo)  // RDO determination V vs W
    {
      // MV search Vehicle
      PUtree local_PU_tree_V;
      double cost_vehicle = find_motion(
        encParam, param, mcEstimate, Window_V, Block0, x0, y0, z0,
        param.motion_block_size, &local_PU_tree_V, bufferPoints);

      // MV search World
      PUtree local_PU_tree_W;
      double cost_world = find_motion(
        encParam, param, mcEstimate, Window_W, Block0, x0, y0, z0,
        param.motion_block_size, &local_PU_tree_W, bufferPoints);

      if (Window_W.size() && cost_world < cost_vehicle) {
        local_PU_tree->popul_flags = local_PU_tree_W.popul_flags;
        local_PU_tree->split_flags = local_PU_tree_W.split_flags;
        local_PU_tree->MVs = local_PU_tree_W.MVs;
        local_PU_tree->isWorld = true;
        //std::cout << "W";
      } else {
        local_PU_tree->popul_flags = local_PU_tree_V.popul_flags;
        local_PU_tree->split_flags = local_PU_tree_V.split_flags;
        local_PU_tree->MVs = local_PU_tree_V.MVs;
        local_PU_tree->isWorld = false;
        //std::cout << "V";
      }
    } else  // fast distortion determination V vs W
    {
      int distW =
        distortionAB(Block0, Window_W);  // +distortionAB(Window_W, Block0);
      int distV =
        distortionAB(Block0, Window_V);  // +distortionAB(Window_V, Block0);
      if (Window_W.size() && distW < distV) {
        // MV search in World
        find_motion(
          encParam, param, mcEstimate, Window_W, Block0, x0, y0, z0,
          param.motion_block_size, local_PU_tree, bufferPoints);
        local_PU_tree->isWorld = true;
      } else {
        // MV search in Vehicle
        find_motion(
          encParam, param, mcEstimate, Window_V, Block0, x0, y0, z0,
          param.motion_block_size, local_PU_tree, bufferPoints);
        local_PU_tree->isWorld = false;
      }
    }

  }  // end global motion determination

  return true;
}

void
encode_splitPU_MV_MC(
  PCCOctree3Node* node0,
  PUtree* local_PU_tree,
  const GeometryParameterSet::Motion& param,
  const int nodeSizeLog2,
  EntropyEncoder* arithmeticEncoder,
  PCCPointSet3& pointPredictor,
#if INTER_HIERARCHICAL
  PCCPointSet3& backPointPredictor,
#endif
  PCCPointSet3* compensatedPointCloud,
  PCCPointSet3& pointPredictorVehicle,
  PCCPointSet3& pointPredictorWorld)
{
  int node_size = 1 << nodeSizeLog2;
  MotionEntropyEncoder motionEncoder(arithmeticEncoder);

  // --------------  non-split / terminal case  ----------------
  if (
    node_size <= param.motion_min_pu_size || !local_PU_tree->split_flags[0]) {
    if (node_size > param.motion_min_pu_size) {
      motionEncoder.encodeSplitPu(0);
    }

    // encode MV
    PCCVector3<int> MV = local_PU_tree->MVs[0];

#if MV_PREDICTION
    // find the neighbor

    PCCVector3<int> predMV = PCCVector3<int>(0, 0, 0);

    int bestDist = 10000;
    int bestIndex = 0;
    if ( g_positionX.size() != 0 )
    {
      for ( int i = 0; i < g_positionX.size(); i++ )
      {
        PCCVector3<int> posRef = PCCVector3<int>(g_positionX[i], g_positionY[i], g_positionZ[i]);
        PCCVector3<int> posCurr = PCCVector3<int>((node0->pos)[0], (node0->pos)[1], (node0->pos)[2]);
        PCCVector3<int> diffPos = posRef - posCurr;
        int dist = diffPos.getNorm();

        if ( dist <= bestDist )
        {
          bestDist = dist;
          bestIndex = i;
        }
      }

      if ( bestDist < 256 )   // to avoid too far distance
      {
        predMV = PCCVector3<int>(g_MVX[bestIndex], g_MVY[bestIndex], g_MVZ[bestIndex]);
      }
    }

    // MV prediction
    PCCVector3<int> diffMV = MV - predMV;

    g_positionX.push_back((node0->pos)[0]);
    g_positionY.push_back((node0->pos)[1]);
    g_positionZ.push_back((node0->pos)[2]);
    g_MVX.push_back(MV[0]);
    g_MVY.push_back(MV[1]);
    g_MVZ.push_back(MV[2]);

#if MV_PREDICTION_RDO
    bool bPredict = false;

    if ( predMV == PCCVector3<int>(0, 0, 0) )
    {
      bPredict = false;
    }
    else
    {
      if ( bestDist < 256 )
      {
        int maxMV = param.motion_window_size / param.motion_precision - 1;
        if ( abs(diffMV[0]) + abs(diffMV[1]) + abs(diffMV[2]) < abs(MV[0]) + abs(MV[1]) + abs(MV[2]) && abs(diffMV[0]) <= maxMV && abs(diffMV[1]) <= maxMV && abs(diffMV[2]) <= maxMV)
        {
          bPredict = true;
        }

        motionEncoder.encodePredPu(bPredict);
      }
    }
#endif

#endif

#if MV_PREDICTION_RDO
    if ( bPredict )
    {
      g_predGood++;
      motionEncoder.encodeVector(
        diffMV, param.motion_precision, param.motion_max_prefix_bits,
        param.motion_max_suffix_bits);
    }
    else
    {
      g_predBad++;
      motionEncoder.encodeVector(
        MV, param.motion_precision, param.motion_max_prefix_bits,
        param.motion_max_suffix_bits);
    }
#else
    motionEncoder.encodeVector(
#if MV_PREDICTION
      diffMV, param.motion_precision, param.motion_max_prefix_bits,
#else
      MV, param.motion_precision, param.motion_max_prefix_bits,
#endif
      param.motion_max_suffix_bits);
#endif
    PCCVector3D MVd = PCCVector3D(double(MV[0]), double(MV[1]), double(MV[2]));

    // find search window
    std::vector<PCCVector3D> Window;
    int xlow = node0->pos[0] - param.motion_window_size;
    int xhigh = node0->pos[0] + node_size + param.motion_window_size;
    int ylow = node0->pos[1] - param.motion_window_size;
    int yhigh = node0->pos[1] + node_size + param.motion_window_size;
    int zlow = node0->pos[2] - param.motion_window_size;
    int zhigh = node0->pos[2] + node_size + param.motion_window_size;

    if (!param.global_motion_enabled)  // no global motion
    {
#if INTER_HIERARCHICAL
      if ( node0->interDir == 0 )
      {
        for (size_t i = 0; i < pointPredictor.getPointCount(); ++i) {
          const PCCVector3D point = pointPredictor[i];
          const PCCVector3<int> point_int =
            PCCVector3<int>(int(point[0]), int(point[1]), int(point[2]));
          if (
            point_int[0] >= xlow && point_int[0] < xhigh && point_int[1] >= ylow
            && point_int[1] < yhigh && point_int[2] >= zlow
            && point_int[2] < zhigh)
            Window.push_back(point);
        }
      }
      else if (node0->interDir == 1)
      {
        for (size_t i = 0; i < backPointPredictor.getPointCount(); ++i) {
          const PCCVector3D point = backPointPredictor[i];
          const PCCVector3<int> point_int =
            PCCVector3<int>(int(point[0]), int(point[1]), int(point[2]));
          if (
            point_int[0] >= xlow && point_int[0] < xhigh && point_int[1] >= ylow
            && point_int[1] < yhigh && point_int[2] >= zlow
            && point_int[2] < zhigh)
            Window.push_back(point);
        }
      }
#else
      for (size_t i = 0; i < pointPredictor.getPointCount(); ++i) {
        const PCCVector3D point = pointPredictor[i];
        const PCCVector3<int> point_int =
          PCCVector3<int>(int(point[0]), int(point[1]), int(point[2]));
        if (
          point_int[0] >= xlow && point_int[0] < xhigh && point_int[1] >= ylow
          && point_int[1] < yhigh && point_int[2] >= zlow
          && point_int[2] < zhigh)
          Window.push_back(point);
      }
#endif
    } else if (local_PU_tree->isWorld)  //  global motion World
    {
      for (size_t i = 0; i < pointPredictorWorld.getPointCount(); ++i) {
        const PCCVector3D point = pointPredictorWorld[i];
        const PCCVector3<int> point_int =
          PCCVector3<int>(int(point[0]), int(point[1]), int(point[2]));
        if (
          point_int[0] >= xlow && point_int[0] < xhigh && point_int[1] >= ylow
          && point_int[1] < yhigh && point_int[2] >= zlow
          && point_int[2] < zhigh)
          Window.push_back(point);
      }
    } else  //  global motion Vehicle
    {
      for (size_t i = 0; i < pointPredictorVehicle.getPointCount(); ++i) {
        const PCCVector3D point = pointPredictorVehicle[i];
        const PCCVector3<int> point_int =
          PCCVector3<int>(int(point[0]), int(point[1]), int(point[2]));
        if (
          point_int[0] >= xlow && point_int[0] < xhigh && point_int[1] >= ylow
          && point_int[1] < yhigh && point_int[2] >= zlow
          && point_int[2] < zhigh)
          Window.push_back(point);
      }
    }

    // create the compensated points
    xlow = node0->pos[0];
    xhigh = node0->pos[0] + node_size;
    ylow = node0->pos[1];
    yhigh = node0->pos[1] + node_size;
    zlow = node0->pos[2];
    zhigh = node0->pos[2] + node_size;
    std::vector<PCCVector3D> pointPredictorMC;
    for (auto& w : Window) {
      // apply best motion
      PCCVector3D wV = w - MVd;
      if (
        wV[0] >= xlow && wV[0] < xhigh && wV[1] >= ylow && wV[1] < yhigh
        && wV[2] >= zlow && wV[2] < zhigh)
        pointPredictorMC.push_back(wV);
    }

    //and make node0 point to them
#if INTER_HIERARCHICAL
    if ( node0->interDir == 0 ) {
      node0->predStart = compensatedPointCloud->getPointCount();
    } else if ( node0->interDir == 1 ) {
      node0->backPredStart = compensatedPointCloud->getPointCount();
    }
#else
    node0->predStart = compensatedPointCloud->getPointCount();
#endif
    compensatedPointCloud->resize(
      compensatedPointCloud->getPointCount() + pointPredictorMC.size());
    size_t counter = node0->predStart;
#if INTER_HIERARCHICAL
    if ( node0->interDir == 1 ) {
      counter = node0->backPredStart;
    }
#endif
    for (const auto& p : pointPredictorMC) {
      auto& predPoint = (*compensatedPointCloud)[counter++];
      predPoint[0] = p[0];
      predPoint[1] = p[1];
      predPoint[2] = p[2];
    }
#if INTER_HIERARCHICAL
    if ( node0->interDir == 0 ) {
      node0->predEnd = compensatedPointCloud->getPointCount();
    } else if ( node0->interDir == 1 ) {
      node0->backPredEnd = compensatedPointCloud->getPointCount();
    }
#else
    node0->predEnd = compensatedPointCloud->getPointCount();
#endif
    node0->isCompensated = true;
    return;
  }

  // --------------- split case ----------------------
  motionEncoder.encodeSplitPu(1);
}

//----------------------------------------------------------------------------

void
decode_splitPU_MV_MC(
  PCCOctree3Node* node0,
  const GeometryParameterSet::Motion& param,
  const int nodeSizeLog2,
  EntropyDecoder* arithmeticDecoder,
  PCCPointSet3& pointPredictor,
#if INTER_HIERARCHICAL
  PCCPointSet3& backPointPredictor,
#endif
  PCCPointSet3* compensatedPointCloud,
  PCCPointSet3& pointPredictorVehicle,
  PCCPointSet3& pointPredictorWorld)
{
  int node_size = 1 << nodeSizeLog2;
  MotionEntropyDecoder motionDecoder(arithmeticDecoder);

  // decode split flag
  bool split = false;
  if (node_size > param.motion_min_pu_size)
    split = motionDecoder.decodeSplitPu();

  if (!split) {  // not split
    // decode MV
    PCCVector3<int> MV;

#if MV_PREDICTION
    PCCVector3<int> diffMV;
#endif

#if MV_PREDICTION
    PCCVector3<int> predMV = PCCVector3<int>(0, 0, 0);
    int bestDist = 10000;
    int bestIndex = 0;
    if (g_positionX.size() != 0)
    {
      for (int i = 0; i < g_positionX.size(); i++)
      {
        PCCVector3<int> posRef = PCCVector3<int>(g_positionX[i], g_positionY[i], g_positionZ[i]);
        PCCVector3<int> posCurr = PCCVector3<int>((node0->pos)[0], (node0->pos)[1], (node0->pos)[2]);
        PCCVector3<int> diffPos = posRef - posCurr;
        int dist = diffPos.getNorm();

        if (dist <= bestDist)
        {
          bestDist = dist;
          bestIndex = i;
        }
      }

      if (bestDist < 256)   // to avoid too far distance
      {
        predMV = PCCVector3<int>(g_MVX[bestIndex], g_MVY[bestIndex], g_MVZ[bestIndex]);
      }
    }

#if MV_PREDICTION_RDO
    bool bPred = false;
    
    if ( predMV == PCCVector3<int>(0, 0, 0) )
    {
      bPred = false;
    }
    else
    {
      if (bestDist < 256)
      {
        bPred = motionDecoder.decodePredPu();
      }
    }
#endif

#endif

#if MV_PREDICTION_RDO
    if (bPred)
    {
      motionDecoder.decodeVector(
        &diffMV, param.motion_max_prefix_bits, param.motion_max_suffix_bits);

      MV = diffMV + predMV;
    }
    else
    {
      motionDecoder.decodeVector(
        &MV, param.motion_max_prefix_bits, param.motion_max_suffix_bits);
    }
#else
    motionDecoder.decodeVector(
#if MV_PREDICTION
      &diffMV, param.motion_max_prefix_bits, param.motion_max_suffix_bits);

    MV = diffMV + predMV;
#else
      &MV, param.motion_max_prefix_bits, param.motion_max_suffix_bits);
#endif
#endif
    MV *= param.motion_precision;

#if MV_PREDICTION
    g_positionX.push_back((node0->pos)[0]);
    g_positionY.push_back((node0->pos)[1]);
    g_positionZ.push_back((node0->pos)[2]);
    g_MVX.push_back(MV[0]);
    g_MVY.push_back(MV[1]);
    g_MVZ.push_back(MV[2]);
#endif


    PCCVector3D MVd = PCCVector3D(double(MV[0]), double(MV[1]), double(MV[2]));

    // find search window
    std::vector<PCCVector3D> Window;
    int xlow = node0->pos[0] - param.motion_window_size;
    int xhigh = node0->pos[0] + node_size + param.motion_window_size;
    int ylow = node0->pos[1] - param.motion_window_size;
    int yhigh = node0->pos[1] + node_size + param.motion_window_size;
    int zlow = node0->pos[2] - param.motion_window_size;
    int zhigh = node0->pos[2] + node_size + param.motion_window_size;

    if (!param.global_motion_enabled)  // no global motion
    {
#if INTER_HIERARCHICAL
      if (node0->interDir == 0)
      {
        for (size_t i = 0; i < pointPredictor.getPointCount(); ++i) {
          const PCCVector3D point = pointPredictor[i];
          const PCCVector3<int> point_int =
            PCCVector3<int>(int(point[0]), int(point[1]), int(point[2]));
          if (
            point_int[0] >= xlow && point_int[0] < xhigh && point_int[1] >= ylow
            && point_int[1] < yhigh && point_int[2] >= zlow
            && point_int[2] < zhigh)
            Window.push_back(point);
        }
      }
      else
      {
        for (size_t i = 0; i < backPointPredictor.getPointCount(); ++i) {
          const PCCVector3D point = backPointPredictor[i];
          const PCCVector3<int> point_int =
            PCCVector3<int>(int(point[0]), int(point[1]), int(point[2]));
          if (
            point_int[0] >= xlow && point_int[0] < xhigh && point_int[1] >= ylow
            && point_int[1] < yhigh && point_int[2] >= zlow
            && point_int[2] < zhigh)
            Window.push_back(point);
        }
      }
#else
      for (size_t i = 0; i < pointPredictor.getPointCount(); ++i) {
        const PCCVector3D point = pointPredictor[i];
        const PCCVector3<int> point_int =
          PCCVector3<int>(int(point[0]), int(point[1]), int(point[2]));
        if (
          point_int[0] >= xlow && point_int[0] < xhigh && point_int[1] >= ylow
          && point_int[1] < yhigh && point_int[2] >= zlow
          && point_int[2] < zhigh)
          Window.push_back(point);
      }
#endif
    } else if (node0->isWorld)  //  global motion World
    {
      for (size_t i = 0; i < pointPredictorWorld.getPointCount(); ++i) {
        const PCCVector3D point = pointPredictorWorld[i];
        const PCCVector3<int> point_int =
          PCCVector3<int>(int(point[0]), int(point[1]), int(point[2]));
        if (
          point_int[0] >= xlow && point_int[0] < xhigh && point_int[1] >= ylow
          && point_int[1] < yhigh && point_int[2] >= zlow
          && point_int[2] < zhigh)
          Window.push_back(point);
      }
    } else  //  global motion Vehicle
    {
      for (size_t i = 0; i < pointPredictorVehicle.getPointCount(); ++i) {
        const PCCVector3D point = pointPredictorVehicle[i];
        const PCCVector3<int> point_int =
          PCCVector3<int>(int(point[0]), int(point[1]), int(point[2]));
        if (
          point_int[0] >= xlow && point_int[0] < xhigh && point_int[1] >= ylow
          && point_int[1] < yhigh && point_int[2] >= zlow
          && point_int[2] < zhigh)
          Window.push_back(point);
      }
    }

    // create the compensated points
    xlow = node0->pos[0];
    xhigh = node0->pos[0] + node_size;
    ylow = node0->pos[1];
    yhigh = node0->pos[1] + node_size;
    zlow = node0->pos[2];
    zhigh = node0->pos[2] + node_size;
    std::vector<PCCVector3D> pointPredictorMC;
    for (auto& w : Window) {
      // apply best motion
      PCCVector3D wV = w - MVd;
      if (
        wV[0] >= xlow && wV[0] < xhigh && wV[1] >= ylow && wV[1] < yhigh
        && wV[2] >= zlow && wV[2] < zhigh)
        pointPredictorMC.push_back(wV);
    }

    //and make node0 point to them
#if INTER_HIERARCHICAL
    if ( node0->interDir == 0 ) {
      node0->predStart = compensatedPointCloud->getPointCount();
    } else if ( node0->interDir == 1 ) {
      node0->backPredStart = compensatedPointCloud->getPointCount();
    }
#else
    node0->predStart = compensatedPointCloud->getPointCount();
#endif
    compensatedPointCloud->resize(
      compensatedPointCloud->getPointCount() + pointPredictorMC.size());
    size_t counter = node0->predStart;
#if INTER_HIERARCHICAL
    if ( node0->interDir == 1 ) {
      counter = node0->backPredStart;
    }
#endif
    for (const auto& p : pointPredictorMC) {
      auto& predPoint = (*compensatedPointCloud)[counter++];
      predPoint[0] = p[0];
      predPoint[1] = p[1];
      predPoint[2] = p[2];
    }
#if INTER_HIERARCHICAL
    if ( node0->interDir == 0 ) {
      node0->predEnd = compensatedPointCloud->getPointCount();
    } else if ( node0->interDir == 1 ) {
      node0->backPredEnd = compensatedPointCloud->getPointCount();
    }
#else
    node0->predEnd = compensatedPointCloud->getPointCount();
#endif
    node0->isCompensated = true;
    return;
  }

  // split; nothing to do
}

//----------------------------------------------------------------------------
// Global motion -------------------------------------
int
distortionAB(std::vector<PCCVector3<int>>& A, std::vector<PCCVector3<int>>& B)
{
  if (A.size() == 0)
    return 0;
  if (B.size() == 0)
    return 1 << 30;
  int D = 0;

  for (const auto& a : A) {
    int dmin = 1 << 30;
    for (const auto& b : B) {
      int L =
        std::abs(a[0] - b[0]) + std::abs(a[1] - b[1]) + std::abs(a[2] - b[2]);
      if (L < dmin)
        dmin = L;
    }
    D += std::log2(1 + dmin);
  }

  return D;
}

//----------------
void
quantizeGlobalMotion(double Mat_GM[4][3], int32_t Mat_GM_Q[4][3])
{
  double scale = 65536.;
  for (int l = 0; l < 4; l++)
    for (int c = 0; c < 3; c++) {
      if (l == c)
        Mat_GM_Q[l][c] = int32_t((Mat_GM[l][c] - 1.) * scale + 0.5);
      else if (l < 3)
        Mat_GM_Q[l][c] = int32_t(Mat_GM[l][c] * scale + 0.5);
      else
        Mat_GM_Q[l][c] = int32_t(Mat_GM[l][c] + 0.5);
    }
}

//----------------
void
applyGlobalMotion(
  std::vector<PCCVector3<int>>& listPoints, double Mat_GM[4][3])
{
  for (auto& b : listPoints) {
    PCCVector3<int> point;
    point[0] = int(
      Mat_GM[0][0] * b[0] + Mat_GM[1][0] * b[1] + Mat_GM[2][0] * b[2]
      + Mat_GM[3][0]);
    point[1] = int(
      Mat_GM[0][1] * b[0] + Mat_GM[1][1] * b[1] + Mat_GM[2][1] * b[2]
      + Mat_GM[3][1]);
    point[2] = int(
      Mat_GM[0][2] * b[0] + Mat_GM[1][2] * b[1] + Mat_GM[2][2] * b[2]
      + Mat_GM[3][2]);

    b = point;
  }
}

//----------------
void
applyGlobalMotion(
  PCCPointSet3& PC,
  int32_t Mat_GM_Q[4][3],
  PCCVector3<double> vehicle_position)
{
  // unquantize
  double Mat_GM[4][3];
  double scale = 65536.;
  for (int l = 0; l < 4; l++)
    for (int c = 0; c < 3; c++) {
      if (l == c)
        Mat_GM[l][c] = double(Mat_GM_Q[l][c]) / scale + 1.;
      else if (l < 3)
        Mat_GM[l][c] = double(Mat_GM_Q[l][c]) / scale;
      else
        Mat_GM[l][c] = double(Mat_GM_Q[l][c]);
    }

  // apply
  for (int n = 0; n < PC.getPointCount(); n++) {
    PCCPoint3D b = PC[n] + vehicle_position;
    PC[n][0] = int(
      Mat_GM[0][0] * b[0] + Mat_GM[1][0] * b[1] + Mat_GM[2][0] * b[2]
      + Mat_GM[3][0] - vehicle_position[0]);
    PC[n][1] = int(
      Mat_GM[0][1] * b[0] + Mat_GM[1][1] * b[1] + Mat_GM[2][1] * b[2]
      + Mat_GM[3][1] - vehicle_position[1]);
    PC[n][2] = int(
      Mat_GM[0][2] * b[0] + Mat_GM[1][2] * b[1] + Mat_GM[2][2] * b[2]
      + Mat_GM[3][2] - vehicle_position[2]);
  }
}

//----------------
double
map_reference(
  std::vector<PCCVector3<int>>& pc_world_target,
  std::vector<PCCVector3<int>>& pointPredictor_centered,
  std::vector<PCCVector3<int>>& pc_world_ref)
{
  std::vector<int> accu_m;
  int Meanm = 0;
  for (const auto& b : pc_world_target) {
    int dmin = 1 << 30;
    PCCVector3<int> closest;
    for (const auto& w : pointPredictor_centered) {
      int L =
        std::abs(w[0] - b[0]) + std::abs(w[1] - b[1]) + std::abs(w[2] - b[2]);
      if (L < dmin) {
        dmin = L;
        closest = w;
      }
    }
    pc_world_ref.push_back(closest);

    accu_m.push_back(dmin);
    Meanm += dmin;
  }
  double err = double(Meanm) / pc_world_target.size() / 3.;

  // eliminate outlayers
  int count = 0;
  std::vector<int>::iterator it_accu = accu_m.begin();
  std::vector<PCCVector3<int>>::iterator it_target = pc_world_target.begin();
  std::vector<PCCVector3<int>>::iterator it_ref = pc_world_ref.begin();
  std::vector<PCCVector3<int>>::iterator it_target_fill =
    pc_world_target.begin();
  std::vector<PCCVector3<int>>::iterator it_ref_fill = pc_world_ref.begin();
  for (; it_accu != accu_m.end(); it_accu++, it_target++, it_ref++) {
    if (*it_accu * accu_m.size() <= 2 * Meanm) {
      *it_target_fill++ = *it_target;
      *it_ref_fill++ = *it_ref;
      count++;
    }
  }

  pc_world_target.resize(count);
  pc_world_ref.resize(count);
  std::cout << "(" << pc_world_ref.size() << "/" << err << ") ";
  return err;
}

//----------------
void
LMS3D(
  std::vector<PCCVector3<int>>& P1,
  std::vector<PCCVector3<int>>& P2,
  std::vector<PCCVector3<int>>& pointPredictor_centered,
  uint32_t maxBB,
  double Mat_GM[4][3])

{
  // determine correlation matrix M in (X,Y,Z,MV_unity)
  const int MV_unity = maxBB >> 4;  //  // for better matrix conditioning
  double M[4][4] = {
    {0., 0., 0., 0.}, {0., 0., 0., 0.}, {0., 0., 0., 0.}, {0., 0., 0., 0.}};
  std::vector<PCCVector3<int>>::iterator it1 = P1.begin();
  for (; it1 != P1.end(); it1++) {
    PCCVector3<double> Pref = PCCVector3<double>(
      double((*it1)[0]), double((*it1)[1]), double((*it1)[2]));
    // X
    M[0][0] += Pref[0] * Pref[0];
    M[0][1] += Pref[0] * Pref[1];
    M[0][2] += Pref[0] * Pref[2];
    M[0][3] += Pref[0] * MV_unity;
    // Y
    M[1][1] += Pref[1] * Pref[1];
    M[1][2] += Pref[1] * Pref[2];
    M[1][3] += Pref[1] * MV_unity;
    // Z
    M[2][2] += Pref[2] * Pref[2];
    M[2][3] += Pref[2] * MV_unity;
    // 1
    M[3][3] += MV_unity * MV_unity;
  }
  M[1][0] = M[0][1];
  M[2][0] = M[0][2];
  M[2][1] = M[1][2];
  M[3][0] = M[0][3];
  M[3][1] = M[1][3];
  M[3][2] = M[2][3];

  // inverse M by Gauss pivoting
  double invM[4][4] = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};
  for (int pivot = 0; pivot < 3; pivot++)  // descente
  {
    double value_pivot = M[pivot][pivot];
    for (int l = pivot + 1; l < 4; l++) {
      double factor = -M[l][pivot] / value_pivot;
      for (int c = 0; c < 4; c++) {
        M[l][c] += M[pivot][c] * factor;
        invM[l][c] += invM[pivot][c] * factor;
      }
    }
  }

  for (int pivot = 3; pivot > 0; pivot--)  // mont�e
  {
    double value_pivot = M[pivot][pivot];
    for (int l = pivot - 1; l >= 0; l--) {
      double factor = -M[l][pivot] / value_pivot;
      for (int c = 0; c < 4; c++) {
        M[l][c] += M[pivot][c] * factor;
        invM[l][c] += invM[pivot][c] * factor;
      }
    }
  }

  for (int pivot = 0; pivot < 4; pivot++)  // normalisation
  {
    double factor = 1 / M[pivot][pivot];
    for (int c = 0; c < 4; c++) {
      invM[pivot][c] *= factor;
    }
  }

  // determine rhs matrix R
  double R[4][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
  it1 = P1.begin();
  std::vector<PCCVector3<int>>::iterator it2 = P2.begin();
  for (; it1 != P1.end(); it1++, it2++) {
    PCCVector3<double> Pref = PCCVector3<double>(
      double((*it1)[0]), double((*it1)[1]), double((*it1)[2]));
    PCCVector3<double> Ptarget = PCCVector3<double>(
      double((*it2)[0]), double((*it2)[1]), double((*it2)[2]));

    // X
    R[0][0] += Ptarget[0] * Pref[0];
    R[1][0] += Ptarget[0] * Pref[1];
    R[2][0] += Ptarget[0] * Pref[2];
    R[3][0] += Ptarget[0] * MV_unity;
    // Y
    R[0][1] += Ptarget[1] * Pref[0];
    R[1][1] += Ptarget[1] * Pref[1];
    R[2][1] += Ptarget[1] * Pref[2];
    R[3][1] += Ptarget[1] * MV_unity;
    // Z
    R[0][2] += Ptarget[2] * Pref[0];
    R[1][2] += Ptarget[2] * Pref[1];
    R[2][2] += Ptarget[2] * Pref[2];
    R[3][2] += Ptarget[2] * MV_unity;
  }

  // apply inv M to R to get the transformation matrix T
  double T[4][3];
  for (int l = 0; l < 4; l++)
    for (int c = 0; c < 3; c++)
      T[l][c] = invM[l][0] * R[0][c] + invM[l][1] * R[1][c]
        + invM[l][2] * R[2][c] + invM[l][3] * R[3][c];

  // deconditioning of 1 <-> MV_unity
  for (int c = 0; c < 3; c++)
    T[3][c] *= double(MV_unity);

  // penalization
  double lambda = 1.0;
  for (int l = 0; l < 4; l++)
    for (int c = 0; c < 3; c++)
      T[l][c] *= lambda;
  T[0][0] += 1 - lambda;
  T[1][1] += 1 - lambda;
  T[2][2] += 1 - lambda;

  // apply T to global motion matrix  Mat_GM
  double Mat_GM1[4][3];  //copy old GM matrix
  for (int l = 0; l < 4; l++)
    for (int c = 0; c < 3; c++)
      Mat_GM1[l][c] = Mat_GM[l][c];

  for (int l = 0; l < 3; l++)  // deformation part
    for (int c = 0; c < 3; c++)
      Mat_GM[l][c] = Mat_GM1[l][0] * T[0][c] + Mat_GM1[l][1] * T[1][c]
        + Mat_GM1[l][2] * T[2][c];

  for (int c = 0; c < 3; c++)  // translation part
    Mat_GM[3][c] = Mat_GM1[3][0] * T[0][c] + Mat_GM1[3][1] * T[1][c]
      + Mat_GM1[3][2] * T[2][c] + T[3][c];
}

//----------------
void
SearchGlobalMotion(
  PCCPointSet3& pointCloud,
  PCCPointSet3& pointPredictor,
  PCCPointSet3& pointPredictorWorld,
  double QS,
  const GeometryParameterSet::Motion& param,
  uint32_t maxBB,
  PCCVector3D minPositions,
  int32_t Mat_GM_Q[4][3])
{
  // ------------- first pass: find world-referential-likely LPU ----
  //std::cout << "First pass: find world-referential-likely LPU " << std::endl;

  // number of LCU
  int bsize = param.motion_block_size;
  int TotalSearch = maxBB / bsize + 1;  // x
  TotalSearch *= maxBB / bsize + 1;     // y
  TotalSearch *= maxBB / bsize + 1;     // z

  // loop on LCU
  std::vector<PCCVector3<int>> pc_likely_world;
  int th_dist = int(1000 * QS);

  for (int x0 = 0; x0 < maxBB; x0 += bsize)  //loop on x
  {
    for (int y0 = 0; y0 < maxBB; y0 += bsize)  //loop on y
    {
      // current block in xy
      std::vector<PCCVector3<int>> Block_current_xy;
      double xlow = x0;
      double xhigh = x0 + bsize;
      double ylow = y0;
      double yhigh = y0 + bsize;
      for (size_t i = 0; i < pointCloud.getPointCount(); ++i) {
        const PCCVector3D point = pointCloud[i];
        if (
          point[0] >= xlow && point[0] < xhigh && point[1] >= ylow
          && point[1] < yhigh)
          Block_current_xy.push_back(
            PCCVector3<int>(int(point[0]), int(point[1]), int(point[2])));
      }

      if (!Block_current_xy.size())
        continue;

      // reference block in xy
      xlow -= th_dist;
      xhigh += th_dist;
      ylow -= th_dist;
      yhigh += th_dist;
      std::vector<PCCVector3<int>> Block_ref_xy;
      for (size_t i = 0; i < pointPredictor.getPointCount(); ++i) {
        const PCCVector3D point = pointPredictor[i];
        if (
          point[0] >= xlow && point[0] < xhigh && point[1] >= ylow
          && point[1] < yhigh)
          Block_ref_xy.push_back(
            PCCVector3<int>(int(point[0]), int(point[1]), int(point[2])));
      }

      for (int z0 = 0; z0 < maxBB; z0 += bsize)  //loop on z
      {
        // current block
        std::vector<PCCVector3<int>> Block_current;
        int zlow = z0;
        int zhigh = z0 + bsize;
        for (const auto& b : Block_current_xy)
          if (b[2] >= zlow && b[2] < zhigh)
            Block_current.push_back(b);

        if (!Block_current.size())
          continue;

        // reference block
        zlow -= th_dist;
        zhigh += th_dist;
        std::vector<PCCVector3<int>> Block_ref;
        for (const auto& w : Block_ref_xy)
          if (w[2] >= zlow && w[2] < zhigh)
            Block_ref.push_back(w);

        if (!Block_ref.size())
          continue;

        // for world, keep only points that have a far closest neighbourg (further than the threshold)
        for (const auto& p : Block_current) {
          int dmin = 1 << 30;
          for (const auto& p2 : Block_ref) {
            int L = std::abs(p[0] - p2[0]) + std::abs(p[1] - p2[1])
              + std::abs(p[2] - p2[2]);
            if (L < dmin)
              dmin = L;
          }
          if (dmin > th_dist)
            pc_likely_world.push_back(p);
        }

      }  // end loop z
    }    // end loop y
  }      // end loop x
  //std::cout << "Likely world has " << pc_likely_world.size() << " points" << std::endl;

  // ------------- second pass: find global motion for world-referential-likely LPU ----
  //std::cout << "Second pass: find global motion for world-referential-likely LPU " << std::endl;

  // center to vehicle ( assumed at (0,0,0) )
  PCCVector3<double> minPositions_local = minPositions * QS;
  PCCVector3<int> minPositions_int = PCCVector3<int>(
    int(minPositions_local[0]), int(minPositions_local[1]),
    int(minPositions_local[2]));
  for (auto& b : pc_likely_world)
    b += minPositions_int;

  std::vector<PCCVector3<int>> pointPredictor_centered;
  for (size_t i = 0; i < pointPredictor.getPointCount(); ++i) {
    const PCCVector3D point = pointPredictor[i] + minPositions_local;
    pointPredictor_centered.push_back(
      PCCVector3<int>(int(point[0]), int(point[1]), int(point[2])));
  }
  std::vector<PCCVector3<int>> pointPredictor_centered0 =
    pointPredictor_centered;

  // global motion found iteratively using LMS
  int NLMS = 10;
  int nb_points = 2000;
  int jump = 1 + (pc_likely_world.size() / nb_points);
  double Mat_GM[4][3] = {{1, 0, 0},
                         {0, 1, 0},
                         {0, 0, 1},
                         {0, 0, 0}};  // global motion is identity as a start

  std::cout << "(points/err) = ";
  for (int i = 0; i < NLMS; i++) {
    // sample pc_likely_world
    std::vector<PCCVector3<int>> pc_world_target;
    std::vector<PCCVector3<int>>::iterator it = pc_likely_world.begin() + i;
    for (int N = i; N < pc_likely_world.size(); N += jump, it += jump)
      pc_world_target.push_back(*it);

    // map reference to pc_world
    std::vector<PCCVector3<int>> pc_world_ref;
    map_reference(pc_world_target, pointPredictor_centered, pc_world_ref);

    // Least Mean Square 3D
    LMS3D(
      pc_world_ref, pc_world_target, pointPredictor_centered, maxBB, Mat_GM);

    // apply global motion
    pointPredictor_centered = pointPredictor_centered0;
    applyGlobalMotion(pointPredictor_centered, Mat_GM);
  }
  std::cout << std::endl;

  // quantize global motion and store it
  quantizeGlobalMotion(Mat_GM, Mat_GM_Q);

  // copy pointPredictor_centered to pointPredictorWorld
  pointPredictorWorld = pointPredictor;
  applyGlobalMotion(pointPredictorWorld, Mat_GM_Q, minPositions_local);
}

//----------------
void
SearchGlobalMotionPerTile(
  PCCPointSet3& pointCloud,
  PCCPointSet3& pointPredictor,
  PCCPointSet3& pointPredictorWorld,
  double QS,
  const GeometryParameterSet::Motion& param,
  uint32_t maxBB,
  PCCVector3D minPositions,
  int32_t Mat_GM_Q[4][3])
{
  SearchGlobalMotion(
    pointCloud, pointPredictor, pointPredictorWorld, QS, param, maxBB,
    minPositions, Mat_GM_Q);
}

//============================================================================

}  // namespace pcc

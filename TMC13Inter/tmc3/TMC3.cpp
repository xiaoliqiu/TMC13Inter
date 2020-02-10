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

#include "TMC3.h"

#include <memory>

#include "PCCTMC3Encoder.h"
#include "PCCTMC3Decoder.h"
#include "constants.h"
#include "program_options_lite.h"
#include "io_tlv.h"
#include "version.h"

using namespace std;
using namespace pcc;

//============================================================================

struct Parameters {
  bool isDecoder;

  // command line parsing should adjust dist2 values according to PQS
  bool positionQuantizationScaleAdjustsDist2;

  // output mode for ply writing (binary or ascii)
  bool outputBinaryPly;

  // when true, configure the encoder as if no attributes are specified
  bool disableAttributeCoding;

  // Frame number of first file in input sequence.
  int firstFrameNum;

  // Number of frames to process.
  int frameCount;

  std::string uncompressedDataPath;
  std::string compressedStreamPath;
  std::string reconstructedDataPath;

  // Filename for saving recoloured point cloud (encoder).
  std::string postRecolorPath;

  // Filename for saving pre inverse scaled point cloud (decoder).
  std::string preInvScalePath;

  pcc::EncoderParams encoder;
  pcc::DecoderParams decoder;

  // todo(df): this should be per-attribute
  ColorTransform colorTransform;

  // todo(df): this should be per-attribute
  int reflectanceScale;
};

//----------------------------------------------------------------------------

class SequenceEncoder : public PCCTMC3Encoder3::Callbacks {
public:
  // NB: params must outlive the lifetime of the decoder.
  SequenceEncoder(Parameters* params);

  int compress(Stopwatch* clock);

protected:
  int compressOneFrame(Stopwatch* clock);

  void onOutputBuffer(const PayloadBuffer& buf) override;
  void onPostRecolour(const PCCPointSet3& cloud) override;

private:
  Parameters* params;
  PCCTMC3Encoder3 encoder;

  std::ofstream bytestreamFile;

  int frameNum;

  // When zero, the current frame is a random access point and must
  // be independently decodable.
  int framesBeforeNextRandomAccessPoint;
};

//----------------------------------------------------------------------------

class SequenceDecoder : public PCCTMC3Decoder3::Callbacks {
public:
  // NB: params must outlive the lifetime of the decoder.
  SequenceDecoder(const Parameters* params);

  int decompress(Stopwatch* clock);

protected:
  void onOutputCloud(const PCCPointSet3& decodedPointCloud) override;

private:
  const Parameters* params;
  PCCTMC3Decoder3 decoder;

  std::ofstream bytestreamFile;

  int frameNum;
  Stopwatch* clock;
};

//============================================================================

int
main(int argc, char* argv[])
{
  cout << "MPEG PCC tmc3 version " << ::pcc::version << endl;

  Parameters params;

  try {
    if (!ParseParameters(argc, argv, params))
      return 1;
  }
  catch (df::program_options_lite::ParseFailure &e) {
    std::cerr << "Error parsing option \"" << e.arg
              << "\" with argument \""<< e.val << "\"." << std::endl;
    return 1;
  }

  // Timers to count elapsed wall/user time
  pcc::chrono::Stopwatch<std::chrono::steady_clock> clock_wall;
  pcc::chrono::Stopwatch<pcc::chrono::utime_inc_children_clock> clock_user;

  clock_wall.start();

  int ret = 0;
  if (params.isDecoder) {
    ret = SequenceDecoder(&params).decompress(&clock_user);
  } else {
    ret = SequenceEncoder(&params).compress(&clock_user);
  }

  clock_wall.stop();

  using namespace std::chrono;
  auto total_wall = duration_cast<milliseconds>(clock_wall.count()).count();
  auto total_user = duration_cast<milliseconds>(clock_user.count()).count();
  std::cout << "Processing time (wall): " << total_wall / 1000.0 << " s\n";
  std::cout << "Processing time (user): " << total_user / 1000.0 << " s\n";

  return ret;
}

//---------------------------------------------------------------------------
// motion parameter derivation
// Setup the block sizes based on a predefined mode:
//  1 => Large scale sparse point clouds (eg cat3)
//  2 => Small voxelised point clouds (eg cat2)
void
deriveMotionParams(int presetMode, EncoderParams* params)
{
  auto scaleFactor = params->sps.seq_source_geom_scale_factor;
  auto& motion = params->gps.motion;

  if (presetMode == 0) {
    motion = GeometryParameterSet::Motion();
    return;
  }

  if (presetMode == 2) {
    motion.motion_block_size = 16 * 2;
    motion.motion_window_size = std::max(2, int(std::round(8 * scaleFactor)));
    motion.motion_precision = 1;
    motion.motion_min_pu_size = motion.motion_block_size >> 2;

    // search parameters
    params->motion.decimate = 7;
    params->motion.sampleW = 1024;
    params->motion.Amotion0 = std::max(1, int(std::round(2 * scaleFactor)));
    params->motion.lambda = 0.5;

    // global motion
    motion.global_motion_enabled = false;
    params->motion.globalMotionInRdo = false;

    return;
  }

  // NB: default parameters are large so that getting it wrong doesn't
  //     run for ever.
  motion.motion_block_size = std::max(64, int(std::round(4096 * scaleFactor)));
  motion.motion_window_size = int(std::round(512 * 1 * scaleFactor * 2));
  motion.motion_precision = std::max(1, int(std::round(4 * scaleFactor)));
  motion.motion_min_pu_size = motion.motion_block_size >> 2;

  // search parameters
  params->motion.decimate = 6;
  params->motion.sampleW = 4;
  params->motion.Amotion0 = motion.motion_window_size >> 2;
  params->motion.lambda = 0.5;

  // global motion
  motion.global_motion_enabled = true;
  params->motion.globalMotionInRdo = false;
}

//---------------------------------------------------------------------------
// :: Command line / config parsing helpers

template<typename T>
static std::istream&
readUInt(std::istream& in, T& val)
{
  unsigned int tmp;
  in >> tmp;
  val = T(tmp);
  return in;
}

static std::istream&
operator>>(std::istream& in, ColorTransform& val)
{
  return readUInt(in, val);
}

namespace pcc {
static std::istream&
operator>>(std::istream& in, AttributeEncoding& val)
{
  return readUInt(in, val);
}
}  // namespace pcc

namespace pcc {
static std::istream&
operator>>(std::istream& in, GeometryCodecType& val)
{
  return readUInt(in, val);
}
}  // namespace pcc

namespace pcc {
static std::ostream&
operator<<(std::ostream& out, const AttributeEncoding& val)
{
  switch (val) {
  case AttributeEncoding::kPredictingTransform: out << "0 (Pred)"; break;
  case AttributeEncoding::kRAHTransform: out << "1 (RAHT)"; break;
  case AttributeEncoding::kLiftingTransform: out << "2 (Lift)"; break;
  }
  return out;
}
}  // namespace pcc

namespace pcc {
static std::ostream&
operator<<(std::ostream& out, const GeometryCodecType& val)
{
  switch (val) {
  case GeometryCodecType::kOctree: out << "1 (Octree)"; break;
  case GeometryCodecType::kTriSoup: out << "2 (TriSoup)"; break;
  }
  return out;
}
}  // namespace pcc

namespace df {
namespace program_options_lite {
template <typename T>
struct option_detail<pcc::PCCVector3<T>> {
  static constexpr bool is_container = true;
  static constexpr bool is_fixed_size = true;
  typedef T* output_iterator;

  static void clear(pcc::PCCVector3<T>& container) { };
  static output_iterator make_output_iterator(pcc::PCCVector3<T>& container) {
    return &container[0];
  }
};
} // program_options_lite
} // df

//---------------------------------------------------------------------------
// :: Command line / config parsing

bool
ParseParameters(int argc, char* argv[], Parameters& params)
{
  namespace po = df::program_options_lite;

  struct {
    AttributeDescription desc;
    AttributeParameterSet aps;
  } params_attr;

  // mapping between attribute name and params_attr index
  std::map<std::string, int> params_attr_idx_map;

  bool print_help = false;
  int motionParamPreset = 0;

  // a helper to set the attribute
  std::function<po::OptionFunc::Func> attribute_setter =
    [&](po::Options&, const std::string& name, po::ErrorReporter) {
      // set default attribute metadata
      // todo(df): allow configuration from cli

      // todo(df): software only supports one instance of each attribute
      params_attr.desc.attr_instance_id = 0;

      // default colour vui information says "unknown"
      params_attr.desc.cicp_colour_primaries_idx = 2;
      params_attr.desc.cicp_transfer_characteristics_idx = 2;
      params_attr.desc.cicp_matrix_coefficients_idx = 2;
      params_attr.desc.cicp_video_full_range_flag = true;

      if (name == "color") {
        params_attr.desc.attr_num_dimensions = 3;
        params_attr.desc.attributeLabel = KnownAttributeLabel::kColour;
      }

      if (name == "reflectance") {
        params_attr.desc.attr_num_dimensions = 1;
        params_attr.desc.attributeLabel = KnownAttributeLabel::kReflectance;
      }

      // copy the current state of parsed attribute parameters
      //
      // NB: this does not cause the default values of attr to be restored
      // for the next attribute block.  A side-effect of this is that the
      // following is allowed leading to attribute foo having both X=1 and
      // Y=2:
      //   "--attr.X=1 --attribute foo --attr.Y=2 --attribute foo"
      //

      // NB: insert returns any existing element
      const auto& it = params_attr_idx_map.insert(
        {name, int(params_attr_idx_map.size())});

      // if insertion was successful
      if (it.second) {
        params.encoder.sps.attributeSets.push_back(params_attr.desc);
        params.encoder.aps.push_back(params_attr.aps);
        params.encoder.attributeSetIdxToApsIdxMap.push_back(it.first->second);
        return;
      }

      // update existing entry
      params.encoder.sps.attributeSets[it.first->second] = params_attr.desc;
      params.encoder.aps[it.first->second] = params_attr.aps;
    };

  /* clang-format off */
  // The definition of the program/config options, along with default values.
  //
  // NB: when updating the following tables:
  //      (a) please keep to 80-columns for easier reading at a glance,
  //      (b) do not vertically align values -- it breaks quickly
  //
  po::Options opts;
  opts.addOptions()
  ("help", print_help, false, "this help text")
  ("config,c", po::parseConfigFile, "configuration file name")

  (po::Section("General"))

  ("mode", params.isDecoder, false,
    "The encoding/decoding mode:\n"
    "  0: encode\n"
    "  1: decode")

  // i/o parameters
  ("firstFrameNum",
     params.firstFrameNum, 0,
     "Frame number for use with interpolating %d format specifiers"
     "in input/output filenames")

  ("frameCount",
     params.frameCount, 1,
     "Number of frames to encode")

  ("reconstructedDataPath",
    params.reconstructedDataPath, {},
    "The ouput reconstructed pointcloud file path (decoder only)")

  ("uncompressedDataPath",
    params.uncompressedDataPath, {},
    "The input pointcloud file path")

  ("compressedStreamPath",
    params.compressedStreamPath, {},
    "The compressed bitstream path (encoder=output, decoder=input)")

  ("postRecolorPath",
    params.postRecolorPath, {},
    "Recolored pointcloud file path (encoder only)")

  ("preInvScalePath",
    params.preInvScalePath, {},
    "Pre inverse scaled pointcloud file path (decoder only)")

  ("outputBinaryPly",
    params.outputBinaryPly, false,
    "Output ply files using binary (or otherwise ascii) format")

  // general
  // todo(df): this should be per-attribute
  ("colorTransform",
    params.colorTransform, COLOR_TRANSFORM_RGB_TO_YCBCR,
    "The colour transform to be applied:\n"
    "  0: none\n"
    "  1: RGB to YCbCr (Rec.709)")

  // todo(df): this should be per-attribute
  ("hack.reflectanceScale",
    params.reflectanceScale, 1,
    "scale factor to be applied to reflectance "
    "pre encoding / post reconstruction")

  // NB: if adding decoder options, uncomment the Decoder section marker
  // (po::Section("Decoder"))

  (po::Section("Encoder"))

  ("seq_bounding_box_xyz0",
    params.encoder.sps.seq_bounding_box_xyz0, {0},
    "seq_bounding_box_xyz0.  NB: seq_bounding_box_whd must be set for this "
    "parameter to have an effect")

  ("seq_bounding_box_whd",
    params.encoder.sps.seq_bounding_box_whd, {0},
    "seq_bounding_box_whd")

  ("positionQuantizationScale",
    params.encoder.sps.seq_source_geom_scale_factor, 1.f,
    "Scale factor to be applied to point positions during quantization process")

  ("positionQuantizationScaleAdjustsDist2",
    params.positionQuantizationScaleAdjustsDist2, false,
    "Scale dist2 values by squared positionQuantizationScale")

  ("mergeDuplicatedPoints",
    params.encoder.gps.geom_unique_points_flag, true,
    "Enables removal of duplicated points")

  ("disableAttributeCoding",
    params.disableAttributeCoding, false,
    "Ignore attribute coding configuration")

  (po::Section("Geometry"))

  ("randomAccessPeriod",
     params.encoder.randomAccessPeriod, 1,
     "Distance (in pictures) between random access points when "
     "encoding a sequence")

  // tools
  ("geometryCodec",
    params.encoder.gps.geom_codec_type, GeometryCodecType::kOctree,
    "Controls the method used to encode geometry:\n"
    "  1: octree (TMC3)\n"
    "  2: trisoup (TMC1)")

  ("bitwiseOccupancyCoding",
    params.encoder.gps.bitwise_occupancy_coding_flag, true,
    "Selects between bitwise and bytewise occupancy coding:\n"
    "  0: bytewise\n"
    "  1: bitwise")

  ("neighbourContextRestriction",
    params.encoder.gps.neighbour_context_restriction_flag, false,
    "Limit geometry octree occupancy contextualisation to sibling nodes")

  ("neighbourAvailBoundaryLog2",
    params.encoder.gps.neighbour_avail_boundary_log2, 0,
    "Defines the avaliability volume for neighbour occupancy lookups."
    " 0: unconstrained")

  ("inferredDirectCodingMode",
    params.encoder.gps.inferred_direct_coding_mode_enabled_flag, true,
    "Permits early termination of the geometry octree for isolated points")

  ("intra_pred_max_node_size_log2",
    params.encoder.gps.intra_pred_max_node_size_log2, 0,
    "octree nodesizes eligible for occupancy intra prediction")

  ("ctxOccupancyReductionFactor",
     params.encoder.gps.geom_occupancy_ctx_reduction_factor, 3,
     "Adjusts the number of contexts used in occupancy coding")

  ("motionParamPreset",
     motionParamPreset, 0,
     "Genaralised derivation of motion compensation parameters:"
     "  1: Large scale point clouds\n"
     "  2: Small voxelised point clouds")

  // (trisoup) geometry parameters
  ("triSoupDepth",  // log2(maxBB+1), where maxBB+1 is analogous to image width
    params.encoder.gps.trisoup_depth, 10,
    "Depth of voxels (reconstructed points) in trisoup geometry")

  ("triSoupLevel",
    params.encoder.gps.trisoup_triangle_level, 7,
    "Level of triangles (reconstructed surface) in trisoup geometry")

  ("triSoupIntToOrigScale",  // reciprocal of positionQuantizationScale
    params.encoder.sps.donotuse_trisoup_int_to_orig_scale, 1.f,
    "orig_coords = integer_coords * intToOrigScale")

  (po::Section("Attributes"))

  // attribute processing
  //   NB: Attribute options are special in the way they are applied (see above)
  ("attribute",
    attribute_setter,
    "Encode the given attribute (NB, must appear after the"
    "following attribute parameters)")

  ("bitdepth",
    params_attr.desc.attr_bitdepth, 8,
    "Attribute bitdepth")

  ("transformType",
    params_attr.aps.attr_encoding, AttributeEncoding::kPredictingTransform,
    "Coding method to use for attribute:\n"
    "  0: Hierarchical neighbourhood prediction\n"
    "  1: Region Adaptive Hierarchical Transform (RAHT)\n"
    "  2: Hierarichical neighbourhood prediction as lifting transform")

  ("rahtLeafDecimationDepth",
    params_attr.aps.raht_binary_level_threshold, 3,
    "Sets coefficients to zero in the bottom n levels of RAHT tree. "
    "Used for chroma-subsampling in attribute=color only.")

  ("rahtQuantizationStep",
    params_attr.aps.quant_step_size_luma, 0,
    "deprecated -- use quantizationStepsLuma")

  ("rahtDepth",
    params_attr.aps.raht_depth, 21,
    "Number of bits for morton representation of an RAHT co-ordinate"
    "component")

  ("numberOfNearestNeighborsInPrediction",
    params_attr.aps.num_pred_nearest_neighbours, 3,
    "Attribute's maximum number of nearest neighbors to be used for prediction")

  ("adaptivePredictionThreshold",
    params_attr.aps.adaptive_prediction_threshold, -1,
    "Neighbouring attribute value difference that enables choice of "
    "single|multi predictors. Applies to transformType=2 only.\n"
    "  -1: auto = 2**(bitdepth-2)")

  ("attributeSearchRange",
    params_attr.aps.search_range, 128,
    "Range for nearest neighbor search")

  ("max_num_direct_predictors",
    params_attr.aps.max_num_direct_predictors, 3,
    "Maximum number of nearest neighbour candidates used in direct"
    "attribute prediction")

  ("levelOfDetailCount",
    params_attr.aps.num_detail_levels, 1,
    "Attribute's number of levels of detail")

  ("quantizationStepLuma",
    params_attr.aps.quant_step_size_luma, 0,
    "Attribute's luma quantization step size")

  ("quantizationStepChroma",
    params_attr.aps.quant_step_size_chroma, 0,
    "Attribute's chroma quantization step size")

  ("dist2",
    params_attr.aps.dist2, {},
    "Attribute's list of squared distances, or initial value for automatic"
    "derivation")
  ;
  /* clang-format on */

  po::setDefaults(opts);
  po::ErrorReporter err;
  const list<const char*>& argv_unhandled =
    po::scanArgv(opts, argc, (const char**)argv, err);

  for (const auto arg : argv_unhandled) {
    err.warn() << "Unhandled argument ignored: " << arg << "\n";
  }

  if (argc == 1 || print_help) {
    po::doHelp(std::cout, opts, 78);
    return false;
  }

  if (int(params.encoder.gps.geom_codec_type) == 0) {
    err.error() << "Bypassed geometry coding is no longer supported\n";
  }

  ////
  // derive parameters
  if (!params.isDecoder) {
    deriveMotionParams(motionParamPreset, &params.encoder);

    if (!motionParamPreset)
      params.encoder.randomAccessPeriod = 1;
  }

  // For trisoup, ensure that positionQuantizationScale is the exact inverse of intToOrigScale.
  if (params.encoder.gps.geom_codec_type == GeometryCodecType::kTriSoup) {
    params.encoder.sps.seq_source_geom_scale_factor =
      1.0f / params.encoder.sps.donotuse_trisoup_int_to_orig_scale;
  }

  // TriSoup geometry is only enabled when trisoup depth > triangle_level.
  // NB: this happens after the scale factor fudge above, since the user
  //     believes they are configuring trisoup
  if (params.encoder.gps.geom_codec_type == GeometryCodecType::kTriSoup) {
    const auto& gps = params.encoder.gps;
    int trisoupNodeSizeLog2 = gps.trisoup_depth - gps.trisoup_triangle_level;
    if (trisoupNodeSizeLog2 <= 0) {
      err.warn() << "TriSoup disabled when depth <= triangle level\n";
      params.encoder.gps.geom_codec_type = GeometryCodecType::kOctree;
    }
  }

  // Certain coding modes are not available when trisoup is enabled.
  // Disable them, and warn if set (they may be set as defaults).
  if (params.encoder.gps.geom_codec_type == GeometryCodecType::kTriSoup) {
    if (!params.encoder.gps.geom_unique_points_flag)
      err.warn() << "TriSoup geometry does not preserve duplicated points\n";

    if (params.encoder.gps.inferred_direct_coding_mode_enabled_flag)
      err.warn() << "TriSoup geometry is incompatable with IDCM\n";

    params.encoder.gps.geom_unique_points_flag = true;
    params.encoder.gps.inferred_direct_coding_mode_enabled_flag = false;
  }

  // support disabling attribute coding (simplifies configuration)
  if (params.disableAttributeCoding) {
    params.encoder.attributeSetIdxToApsIdxMap.clear();
    params.encoder.sps.attributeSets.clear();
    params.encoder.aps.clear();
  }

  // fixup any per-attribute settings
  int numAttributeSets = params.encoder.sps.attributeSets.size();
  for (int attrIdx = 0; attrIdx < numAttributeSets; attrIdx++) {
    int apsIdx = params.encoder.attributeSetIdxToApsIdxMap[attrIdx];
    auto& attr_sps = params.encoder.sps.attributeSets[attrIdx];
    auto& attr_aps = params.encoder.aps[apsIdx];

    // Avoid wasting bits signalling chroma quant step size for reflectance
    if (attr_sps.attributeLabel == KnownAttributeLabel::kReflectance) {
      attr_aps.quant_step_size_chroma = 0;
    }

    bool isLifting =
      attr_aps.attr_encoding == AttributeEncoding::kPredictingTransform
      || attr_aps.attr_encoding == AttributeEncoding::kLiftingTransform;

    // derive the dist2 values based on an initial value
    if (isLifting && !attr_aps.dist2.empty()) {
      if (attr_aps.dist2.size() < attr_aps.num_detail_levels) {
        attr_aps.dist2.resize(attr_aps.num_detail_levels);
        const double distRatio = 4.0;
        uint64_t d2 = attr_aps.dist2[0];
        for (int i = 0; i < attr_aps.num_detail_levels; ++i) {
          attr_aps.dist2[i] = d2;
          d2 = uint64_t(std::round(distRatio * d2));
        }
      }
    }

    // In order to simplify specification of dist2 values, which are
    // depending on the scale of the coded point cloud, the following
    // adjust the dist2 values according to PQS.  The user need only
    // specify the unquantised PQS value.
    if (params.positionQuantizationScaleAdjustsDist2) {
      double pqs = params.encoder.sps.seq_source_geom_scale_factor;
      double pqs2 = pqs * pqs;
      for (auto& dist2 : attr_aps.dist2)
        dist2 = int64_t(std::round(pqs2 * dist2));
    }

    // Set default threshold based on bitdepth
    if (attr_aps.adaptive_prediction_threshold == -1) {
      attr_aps.adaptive_prediction_threshold = 1
        << (attr_sps.attr_bitdepth - 2);
    }

    if (attr_aps.attr_encoding == AttributeEncoding::kLiftingTransform) {
      attr_aps.adaptive_prediction_threshold = 0;
    }

    // For RAHT, ensure that the unused lod count = 0 (prevents mishaps)
    if (attr_aps.attr_encoding == AttributeEncoding::kRAHTransform) {
      attr_aps.num_detail_levels = 0;
      attr_aps.adaptive_prediction_threshold = 0;

      // todo(df): suggest chroma quant_step_size for raht
      attr_aps.quant_step_size_chroma = 0;
    }
  }

  // sanity checks

  if (params.encoder.gps.intra_pred_max_node_size_log2)
    if (!params.encoder.gps.neighbour_avail_boundary_log2)
      err.error() << "Geometry intra prediction requires finite"
                     "neighbour_avail_boundary_log2\n";

  for (int attrIdx = 0; attrIdx < numAttributeSets; attrIdx++) {
    int apsIdx = params.encoder.attributeSetIdxToApsIdxMap[attrIdx];
    auto& attr_sps = params.encoder.sps.attributeSets[attrIdx];
    auto& attr_aps = params.encoder.aps[apsIdx];
    auto label = attr_sps.attributeLabel;

    bool isLifting =
      attr_aps.attr_encoding == AttributeEncoding::kPredictingTransform
      || attr_aps.attr_encoding == AttributeEncoding::kLiftingTransform;

    if (attr_sps.attributeLabel == KnownAttributeLabel::kColour) {
      // todo(??): permit relaxing of the following constraint
      if (attr_sps.attr_bitdepth > 8)
        err.error() << label << ".bitdepth must be less than 9\n";
    }

    if (attr_sps.attributeLabel == KnownAttributeLabel::kReflectance) {
      if (attr_sps.attr_bitdepth > 16)
        err.error() << label << ".bitdepth must be less than 17\n";
    }

    if (isLifting) {
      int lod = attr_aps.num_detail_levels;

      if (lod > 255 || lod < 1) {
        err.error() << label
                    << ".levelOfDetailCount must be in the range [1,255]\n";
      }
      if (attr_aps.dist2.size() != lod) {
        err.error() << label << ".dist2 does not have " << lod
                    << " entries\n";
      }

      if (attr_aps.adaptive_prediction_threshold < 0) {
        err.error() << label
                    << ".adaptivePredictionThreshold must be positive\n";
      }

      if (
        attr_aps.num_pred_nearest_neighbours
        > kAttributePredictionMaxNeighbourCount) {
        err.error() << label
                    << ".numberOfNearestNeighborsInPrediction must be <= "
                    << kAttributePredictionMaxNeighbourCount << "\n";
      }
    }
  }

  // check required arguments are specified

  if (!params.isDecoder && params.uncompressedDataPath.empty())
    err.error() << "uncompressedDataPath not set\n";

  if (params.isDecoder && params.reconstructedDataPath.empty())
    err.error() << "reconstructedDataPath not set\n";

  if (params.compressedStreamPath.empty())
    err.error() << "compressedStreamPath not set\n";

  // report the current configuration (only in the absence of errors so
  // that errors/warnings are more obvious and in the same place).
  if (err.is_errored)
    return false;

  // Dump the complete derived configuration
  cout << "+ Effective configuration parameters\n";

  po::dumpCfg(cout, opts, "General", 4);
  if (params.isDecoder) {
    po::dumpCfg(cout, opts, "Decoder", 4);
  } else {
    po::dumpCfg(cout, opts, "Encoder", 4);
    po::dumpCfg(cout, opts, "Geometry", 4);

    for (int attrIdx = 0; attrIdx < numAttributeSets; attrIdx++) {
      int apsIdx = params.encoder.attributeSetIdxToApsIdxMap[attrIdx];

      // NB: when dumping the config, opts references params_attr
      params_attr.desc = params.encoder.sps.attributeSets[attrIdx];
      params_attr.aps = params.encoder.aps[apsIdx];

      cout << "    " << params_attr.desc.attributeLabel << "\n";
      po::dumpCfg(cout, opts, "Attributes", 8);
    }
  }

  cout << endl;

  return true;
}

//============================================================================

SequenceEncoder::SequenceEncoder(Parameters* params) : params(params)
{}

//----------------------------------------------------------------------------

int
SequenceEncoder::compress(Stopwatch* clock)
{
  encoder.configure(&params->encoder);

  bytestreamFile.open(params->compressedStreamPath, ios::binary);
  if (!bytestreamFile.is_open()) {
    return -1;
  }

  framesBeforeNextRandomAccessPoint = 0;
  const int lastFrameNum = params->firstFrameNum + params->frameCount;
  for (frameNum = params->firstFrameNum; frameNum < lastFrameNum; frameNum++)
  {
    params->encoder.randomAccessPoint = !framesBeforeNextRandomAccessPoint;
    if (compressOneFrame(clock))
      return -1;

    if (!framesBeforeNextRandomAccessPoint) {
      framesBeforeNextRandomAccessPoint = params->encoder.randomAccessPeriod;
    }
    framesBeforeNextRandomAccessPoint--;
  }

  std::cout << "Total bitstream size " << bytestreamFile.tellp() << " B\n";
  bytestreamFile.close();

  return 0;
}

//----------------------------------------------------------------------------

int
SequenceEncoder::compressOneFrame(Stopwatch* clock)
{
  std::string srcName{expandNum(params->uncompressedDataPath, frameNum)};
  PCCPointSet3 pointCloud;
  if (!pointCloud.read(srcName) || pointCloud.getPointCount() == 0) {
    cout << "Error: can't open input file!" << endl;
    return -1;
  }

  // Sanitise the input point cloud
  // todo(df): remove the following with generic handling of properties
  bool codeColour = std::any_of(
    params->encoder.sps.attributeSets.begin(),
    params->encoder.sps.attributeSets.end(),
    [](const AttributeDescription& desc) {
      return desc.attributeLabel == KnownAttributeLabel::kColour;
    });
  if (!codeColour)
    pointCloud.removeColors();
  assert(codeColour == pointCloud.hasColors());

  bool codeReflectance = std::any_of(
    params->encoder.sps.attributeSets.begin(),
    params->encoder.sps.attributeSets.end(),
    [](const AttributeDescription& desc) {
      return desc.attributeLabel == KnownAttributeLabel::kReflectance;
    });
  if (!codeReflectance)
    pointCloud.removeReflectances();
  assert(codeReflectance == pointCloud.hasReflectances());

  clock->start();

  if (params->colorTransform == COLOR_TRANSFORM_RGB_TO_YCBCR) {
    pointCloud.convertRGBToYUV();
  }

  if (params->reflectanceScale > 1 && pointCloud.hasReflectances()) {
    const auto pointCount = pointCloud.getPointCount();
    for (size_t i = 0; i < pointCount; ++i) {
      int val = pointCloud.getReflectance(i) / params->reflectanceScale;
      pointCloud.setReflectance(i, val);
    }
  }

  // The reconstructed point cloud
  std::unique_ptr<PCCPointSet3> reconPointCloud;
  if (!params->reconstructedDataPath.empty()) {
    reconPointCloud.reset(new PCCPointSet3);
  }

  auto bytestreamLenFrameStart = bytestreamFile.tellp();

  int ret = encoder.compress(
    pointCloud, params->encoder, this, reconPointCloud.get());
  if (ret) {
    cout << "Error: can't compress point cloud!" << endl;
    return -1;
  }

  auto bytestreamLenFrameEnd = bytestreamFile.tellp();
  int frameLen = bytestreamLenFrameEnd - bytestreamLenFrameStart;

  std::cout << "Total frame size " << frameLen << " B" << std::endl;

  clock->stop();

  if (!params->reconstructedDataPath.empty()) {
    if (params->colorTransform == COLOR_TRANSFORM_RGB_TO_YCBCR) {
      reconPointCloud->convertYUVToRGB();
    }

    if (params->reflectanceScale > 1 && reconPointCloud->hasReflectances()) {
      const auto pointCount = reconPointCloud->getPointCount();
      for (size_t i = 0; i < pointCount; ++i) {
        int val =
          reconPointCloud->getReflectance(i) * params->reflectanceScale;
        reconPointCloud->setReflectance(i, val);
      }
    }

    std::string recName{expandNum(params->reconstructedDataPath, frameNum)};
    reconPointCloud->write(recName, !params->outputBinaryPly);
  }

  return 0;
}

//----------------------------------------------------------------------------

void
SequenceEncoder::onOutputBuffer(const PayloadBuffer& buf)
{
  writeTlv(buf, bytestreamFile);
}

//----------------------------------------------------------------------------

void
SequenceEncoder::onPostRecolour(const PCCPointSet3& cloud)
{
  if (params->postRecolorPath.empty()) {
    return;
  }

  std::string plyName{expandNum(params->postRecolorPath, frameNum)};

  // todo(df): stop the clock
  if (params->colorTransform != COLOR_TRANSFORM_RGB_TO_YCBCR) {
    cloud.write(plyName);
    return;
  }

  PCCPointSet3 tmpCloud(cloud);
  tmpCloud.convertYUVToRGB();
  tmpCloud.write(plyName);
}

//============================================================================

SequenceDecoder::SequenceDecoder(const Parameters* params) : params(params)
{}

//----------------------------------------------------------------------------

int
SequenceDecoder::decompress(Stopwatch* clock)
{
  ifstream fin(params->compressedStreamPath, ios::binary);
  if (!fin.is_open()) {
    return -1;
  }

  frameNum = params->firstFrameNum;
  this->clock = clock;

  PayloadBuffer buf;

  clock->start();

  while (true) {
    PayloadBuffer* buf_ptr = &buf;
    readTlv(fin, &buf);

    // at end of file (or other error), flush decoder
    if (!fin)
      buf_ptr = nullptr;

    if (decoder.decompress(params->decoder, buf_ptr, this)) {
      cout << "Error: can't decompress point cloud!" << endl;
      return -1;
    }

    if (!buf_ptr)
      break;
  }

  fin.clear();
  fin.seekg(0, ios_base::end);
  std::cout << "Total bitstream size " << fin.tellg() << " B" << std::endl;

  clock->stop();

  return 0;
}

//----------------------------------------------------------------------------

void
SequenceDecoder::onOutputCloud(const PCCPointSet3& decodedPointCloud)
{
  // copy the point cloud in order to modify it according to the output options
  PCCPointSet3 pointCloud(decodedPointCloud);

  if (params->colorTransform == COLOR_TRANSFORM_RGB_TO_YCBCR) {
    pointCloud.convertYUVToRGB();
  }

  if (params->reflectanceScale > 1 && pointCloud.hasReflectances()) {
    const auto pointCount = pointCloud.getPointCount();
    for (size_t i = 0; i < pointCount; ++i) {
      int val = pointCloud.getReflectance(i) * params->reflectanceScale;
      pointCloud.setReflectance(i, val);
    }
  }

  // Dump the decoded colour using the pre inverse scaled geometry
  if (!params->preInvScalePath.empty()) {
    std::string filename{expandNum(params->preInvScalePath, frameNum)};
    pointCloud.write(params->preInvScalePath, !params->outputBinaryPly);
  }

  decoder.inverseQuantization(pointCloud);

  clock->stop();

  std::string decName{expandNum(params->reconstructedDataPath, frameNum)};
  if (!pointCloud.write(decName, !params->outputBinaryPly)) {
    cout << "Error: can't open output file!" << endl;
  }

  clock->start();

  // todo(df): frame number should be derived from the bitstream
  frameNum++;
}

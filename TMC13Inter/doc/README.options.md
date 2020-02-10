General options
---------------

### `--help`
Print a list of available command line (and configuration file) options
along with their default values and exit.

### `--config=FILE`, `-c`
This specifies a configuration file to be immediately loaded.

### `--mode=VALUE`
This option selects the codec's mode of operation.  A value of 0 enables
encoding functionality.  A value of 1 switches to decoding mode.


I/O parameters
--------------

### `--firstFrameNum=INT-VALUE`
The initial frame number of the input or output sequence.
The software replaces any instance of a '%d' printf format directive
with the current frame number when evaluating the following options:

- uncompressedDataPath
- reconstructedDataPath
- postRecolourPath
- preInvScalePath

NB: When decoding, this option relates only to the output file names.

In order to have the decoder produce identically numbered output ply
files as the encoder input, specify the same value of firstFrameNum for
the decoder.

### `--frameCount=INT-VALUE`
(Encoder only)
The number of frames to be encoded.

### `--uncompressedDataPath=FILE`
(Encoder only)
The input source point cloud to be compressed.  The first instance of
'%d' in FILE will be expanded with the current frame number.

### `--compressedStreamPath=FILE`
The compressed bitstream file output when encoding or input when decoding.

### `--reconstructedDataPath=FILE`
The reconstructed point cloud file.  When encoding, the output is the
locally decoded picture.  It is expected that the reconstructed output
of the encoder and decoder match exactly.

The first instance of '%d' in FILE will be expanded with the current
frame number.

### `--postRecolourPath=FILE`
(Encoder only)
As part of the encoding process, it may be necessary to re-colour the
point cloud if the point geometry is altered.  This diagnostic output
file corresponds to the re-coloured point cloud prior to attribute
coding without output geometry scaling.

The first instance of '%d' in FILE will be expanded with the current
frame number.

### `--preInvScalePath=FILE`
(Decoder only)
This diagnostic output corresponds to the decoded point cloud (geometry
and attributes) prior to output geometry scaling.

When compared to the output of `postRecolourPath`, the performance of
attribute coding may be directly measured without being confounded
by any geometry losses.

The first instance of '%d' in FILE will be expanded with the current
frame number.

### `--outputBinaryPly=0|1`
Sets the output format of PLY files (Binary=1, ASCII=0).  Reading and
writing binary PLY files is more efficient than the ASCII variant,
but are less suited to simple scripts and direct human inspection.

If outputting non-integer point co-ordinates (eg, due to the output
geometry scaling), the precision of the binary and ASCII versions are
not identical.

### `--colourTransform=0|1`
Controls the use of a colour space transformation before attribute
coding and after decoding.

  | Value | Description            |
  |:-----:| ---------------------- |
  | 0     | none                   |
  | 1     | RGB to YCbCr (Rec.709) |

### `--hack.reflectanceScale=0|1`
Some input data uses 8-bit reflectance data scaled by 255 and represented
using 16-bit attributes.  This option enables a conversion of 16-bit to
8-bit at the encoder, and the corresponding conversion from 8-bit back to
16-bit at the decoder.  If the original data has been scaled by 255, the
conversion process is lossless.


Decoder-specific options
========================
There are no decoder specific options at this time.


Encoder-specific options
========================

### `--seq_bounding_box_xyz0=INT-VALUE,INT-VALUE,INT-VALUE`
Specifies the position of the bottom-left of the sequence
bounding box.  All points are translated by this value prior
to encoding.

### `--seq_bounding_box_whd=INT-VALUE,INT-VALUE,INT-VALUE`
Specifies the width, height, and depth of the sequence
bounding box.

### `--positionQuantizationScale=REAL-FACTOR`
Prior to encoding, scale the point cloud geometry by multiplying each
co-ordinate by the real *FACTOR* and rounding to integer precision.  The
scale factor is written to the bitstream and a decoder may use it to
provide output at the original scale.

NB: when using trisoup geometry coding, use `triSoupIntToOrigScale`
instead of this option.

### `--positionQuantizationScaleAdjustsDist2=0|1`
This option simplifies the specification of the per-attribute `dist2`
parameter.

The squared distance threshold used for generating levels-of-detail in
attribute coding is dependent on the point cloud density and is therefore
affected by geometry quantization.  When this parameter is enabled,
`dist2` values are scaled by `positionQuantizationScale` squared, thereby
allowing `dist2` to be specified as an intrinsic property of the source
sequence.

### `--mergeDuplicatedPoints=0|1`
Controls the ability to code duplicate points.  When duplicate point
merging is enabled, bitstream syntax related to duplicate points is
disabled and a pre-filtering process is used to remove co-located points.

### `--disableAttributeCoding=0|1`
This option instructs the encoder to ignore all options relating to
attribute coding, as if they had never been configured.

### `--randomAccessPeriod=INT-VALUE`
Sets the distance (in coded frames) between random access points
(intra frames) when encoding a sequence.  A value of 1 implies
intra only coding.

Geometry coding
---------------

### `--geometryCodec=1|2`
Selects the geometry coding method.  Both available methods use octree
coding as a base, where points are represented by the occupancy of an
octree.  The trisoup method terminates the octree coding early and
continues by encoding triangles which are subsequently voxelised to
produce points.

  | Value | Description      |
  |:-----:| ---------------- |
  | 0     | invalid          |
  | 1     | octree only      |
  | 2     | octree + trisoup |

### `--bitwiseOccupancyCoding=0|1`
In octree geometry coding, there are both byte-wise and bit-wise tools to
encode the occupancy data.  This option selects between the two methods.

### `--neighbourContextRestriction=0|1`
Octree occupancy coding is contextualised in part by the occupancy of
neighbouring octree nodes.  The neighbour context restriction limits
the use of neighbouring nodes to direct octree siblings.

NB: This option conflicts with `neighbourAvailBoundaryLog2`.  It is
necessary to set `neighbourAvailBoundaryLog2`=0 when
`neighbourContextRestriction`=1.

### `--neighbourAvailBoundaryLog2=INT-VALUE`
Defines the volume within which octree nodes are considered available
for use in occupancy contextualisation and intra occupancy prediction.

A value of 0 indicates that no constraint is applied.

The software currently supports a maximum value of 8 or 9 when
intra occupancy prediction prediction is enabled or disabled
respectively.

### `--inferredDirectCodingMode=0|1`
Controls the use of early termination of the geometry octree for
isolated points.

### `--intra_pred_max_node_size_log2=INT-VALUE`
Intra occupancy prediction uses an octree node's neighbours to predict
its occupancy.  The prediction mode is enabled for octree nodes smaller
than or equal to the configured size.  A value of 0 disables intra
occupancy prediction.

### `--ctxOccupancyReductionFactor=INT-VALUE`
Adjusts the number of contexts used in bit-wise occupancy coding.
The total number of contexts used is 256 >> *VALUE*.

NB: the final standard is expected to define this factor as a constant.

### `--triSoupDepth=INT-VALUE`
Configures the maximum log2 node size of the geometry octree used by
trisoup geometry coding.

### `--triSoupLevel=INT-VALUE`
Configures the log2 octree node size used for triangulation (the
reconstructed surface) in trisoup geometry coding.

### `--triSoupIntToOrigScale=REAL-VALUE`
This option is the inverse of `positionQuantizationScale`, and is used to
configure input scaling when trisoup geometry coding is enabled.

NB: this option will be removed in a future software version.

### `motionParamPreset=1|2`
Preset options for controlling motion estimation and compensation
parameters.

  | Value | Description |
  |:-----:|:------------]
  |   1   | Large scale point clouds |
  |   2   | Small voxelised point clouds |

Attribute coding
----------------

The codec may be configured to represent one or more attributes.
The configuration of each attribute is independent from all others.
To configure coding of an attribute, first set the attribute options,
then save the configuration using the `attribute` option.

### `--attribute=NAME`
Saves the current attribute configuration for coding the named attribute.

  | Name        | Description |
  |:----------- |---|
  | colour      | r, g, and b properties as a tri-stimulus attribute |
  | reflectance | refc or reflectance property as a single-stimulus attribute |

This option must be specified after the options corresponding to
the attribute.

### `--bitdepth=INT-VALUE`
The bitdepth of the attribute data.  NB, this is not necessarily the
same as the bitdepth of the PLY property.  

### `--transformType=0|1|2`
Coding method to use for the current attribute:

  | Value | Description                                                |
  |:-----:| ---------------------------------------------------------- |
  | 0     | Hierarchical neighbourhood prediction                      |
  | 1     | Region Adaptive Hierarchical Transform (RAHT)              |
  | 2     | Hierarchical neighbourhood prediction as lifting transform |

### `--rahtLeafDecimationDepth=INT-VALUE`
Sets coefficients to zero in the bottom n levels of the RAHT tree.
This option provides a means to perform chroma-subsampling.  Applies
when `attribute=colour` only.

### `--rahtQuantizationStep=INT-VALUE`
Deprecated -- use `quantizationStepsLuma`.

### `--rahtDepth=INT-VALUE`
Number of bits for Morton representation of RAHT co-ordinate
components.

### `--numberOfNearestNeighboursInPrediction=INT-VALUE`
Attribute's maximum number of nearest neighbours to be used for
prediction.

### `--adaptivePredictionThreshold=INT-VALUE`
Neighbouring attribute value difference that enables choice of
single|multi predictors. Applies to transformType=2 only.
A value of -1 is replaced by 2**(bitdepth-2).

### `--attributeSearchRange=INT-VALUE`
Range for nearest neighbour search.

### `--max_num_direct_predictors=INT-VALUE`
Maximum number of nearest neighbour candidates used in direct
attribute prediction.

### `--levelOfDetailCount=INT-VALUE`
Attribute's number of levels of detail.

### `--quantizationStepLuma=INT-VALUE`
Attribute's luma quantization step size.

### `--quantizationStepChroma=INT-VALUE`
Attribute's chroma quantization step size.  Only applies when
`attribute=colour`.

### `--dist2=INT-VALUE|INT-VALUE-LIST`
Attribute's list of squared distances, or initial value for automatic
derivation.

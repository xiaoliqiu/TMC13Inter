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

#include "tables.h"

// indicates impossible values in the following table
static const int x = 0;

const uint8_t pcc::kNeighPattern64to10[64] = {
  0, 1, 1, 2, 1, 3, 3, 4, 1, 3, 3, 4, 2, 4, 4, 5, 1, 3, 3, 4, 3, 6,
  6, 7, 3, 6, 6, 7, 4, 7, 7, 8, 1, 3, 3, 4, 3, 6, 6, 7, 3, 6, 6, 7,
  4, 7, 7, 8, 2, 4, 4, 5, 4, 7, 7, 8, 4, 7, 7, 8, 5, 8, 8, 9};

const uint8_t pcc::kNeighPattern64to6[64] = {
  0, 4, 4, x, 4, 2, 2, x, 4, 2, 2, x, x, x, x, x, 5, 3, 3, x, 3, 1,
  1, x, 3, 1, 1, x, x, x, x, x, 5, 3, 3, x, 3, 1, 1, x, 3, 1, 1, x,
  x, x, x, x, x, x, x, x, x, x, x, x, x, x, x, x, x, x, x, x};

const uint8_t pcc::kNeighPattern10to7[10] = {0, 1, 2, 3, 4, 5, 3, 4, 5, 6};

const uint8_t pcc::kNeighPattern7to5[7] = {0, 1, 2, 1, 4, 4, 3};

const uint8_t pcc::kOccMapRotateXIdFromPattern[64] = {
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 2, 1, 0,
  0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 1, 1, 2, 1, 0, 0, 0, 1, 0, 0, 0,
  2, 0, 0, 0, 0, 2, 2, 1, 2, 3, 3, 0, 2, 3, 3, 0, 1, 0, 0, 0};

const uint8_t pcc::kOccMapRotateZIdFromPatternXY[16] = {
  0, 0, 2, 0, 3, 3, 2, 1, 1, 0, 1, 3, 1, 2, 0, 0};

const bool pcc::kOccMapRotateYIdFromPattern[64] = {
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0,
  0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
  1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0};

const uint8_t pcc::kOccMapRotateX270[256] = {
  0,   4,   1,   5,   8,   12,  9,   13,  2,   6,   3,   7,   10,  14,  11,
  15,  64,  68,  65,  69,  72,  76,  73,  77,  66,  70,  67,  71,  74,  78,
  75,  79,  16,  20,  17,  21,  24,  28,  25,  29,  18,  22,  19,  23,  26,
  30,  27,  31,  80,  84,  81,  85,  88,  92,  89,  93,  82,  86,  83,  87,
  90,  94,  91,  95,  128, 132, 129, 133, 136, 140, 137, 141, 130, 134, 131,
  135, 138, 142, 139, 143, 192, 196, 193, 197, 200, 204, 201, 205, 194, 198,
  195, 199, 202, 206, 203, 207, 144, 148, 145, 149, 152, 156, 153, 157, 146,
  150, 147, 151, 154, 158, 155, 159, 208, 212, 209, 213, 216, 220, 217, 221,
  210, 214, 211, 215, 218, 222, 219, 223, 32,  36,  33,  37,  40,  44,  41,
  45,  34,  38,  35,  39,  42,  46,  43,  47,  96,  100, 97,  101, 104, 108,
  105, 109, 98,  102, 99,  103, 106, 110, 107, 111, 48,  52,  49,  53,  56,
  60,  57,  61,  50,  54,  51,  55,  58,  62,  59,  63,  112, 116, 113, 117,
  120, 124, 121, 125, 114, 118, 115, 119, 122, 126, 123, 127, 160, 164, 161,
  165, 168, 172, 169, 173, 162, 166, 163, 167, 170, 174, 171, 175, 224, 228,
  225, 229, 232, 236, 233, 237, 226, 230, 227, 231, 234, 238, 235, 239, 176,
  180, 177, 181, 184, 188, 185, 189, 178, 182, 179, 183, 186, 190, 187, 191,
  240, 244, 241, 245, 248, 252, 249, 253, 242, 246, 243, 247, 250, 254, 251,
  255};

const uint8_t pcc::kOccMapRotateX090[256] = {
  0,   2,   8,   10,  1,   3,   9,   11,  4,   6,   12,  14,  5,   7,   13,
  15,  32,  34,  40,  42,  33,  35,  41,  43,  36,  38,  44,  46,  37,  39,
  45,  47,  128, 130, 136, 138, 129, 131, 137, 139, 132, 134, 140, 142, 133,
  135, 141, 143, 160, 162, 168, 170, 161, 163, 169, 171, 164, 166, 172, 174,
  165, 167, 173, 175, 16,  18,  24,  26,  17,  19,  25,  27,  20,  22,  28,
  30,  21,  23,  29,  31,  48,  50,  56,  58,  49,  51,  57,  59,  52,  54,
  60,  62,  53,  55,  61,  63,  144, 146, 152, 154, 145, 147, 153, 155, 148,
  150, 156, 158, 149, 151, 157, 159, 176, 178, 184, 186, 177, 179, 185, 187,
  180, 182, 188, 190, 181, 183, 189, 191, 64,  66,  72,  74,  65,  67,  73,
  75,  68,  70,  76,  78,  69,  71,  77,  79,  96,  98,  104, 106, 97,  99,
  105, 107, 100, 102, 108, 110, 101, 103, 109, 111, 192, 194, 200, 202, 193,
  195, 201, 203, 196, 198, 204, 206, 197, 199, 205, 207, 224, 226, 232, 234,
  225, 227, 233, 235, 228, 230, 236, 238, 229, 231, 237, 239, 80,  82,  88,
  90,  81,  83,  89,  91,  84,  86,  92,  94,  85,  87,  93,  95,  112, 114,
  120, 122, 113, 115, 121, 123, 116, 118, 124, 126, 117, 119, 125, 127, 208,
  210, 216, 218, 209, 211, 217, 219, 212, 214, 220, 222, 213, 215, 221, 223,
  240, 242, 248, 250, 241, 243, 249, 251, 244, 246, 252, 254, 245, 247, 253,
  255};

const uint8_t pcc::kOccMapRotateX270Y180[256] = {
  0,  128, 32, 160, 64, 192, 96,  224, 16, 144, 48, 176, 80, 208, 112, 240,
  8,  136, 40, 168, 72, 200, 104, 232, 24, 152, 56, 184, 88, 216, 120, 248,
  2,  130, 34, 162, 66, 194, 98,  226, 18, 146, 50, 178, 82, 210, 114, 242,
  10, 138, 42, 170, 74, 202, 106, 234, 26, 154, 58, 186, 90, 218, 122, 250,
  4,  132, 36, 164, 68, 196, 100, 228, 20, 148, 52, 180, 84, 212, 116, 244,
  12, 140, 44, 172, 76, 204, 108, 236, 28, 156, 60, 188, 92, 220, 124, 252,
  6,  134, 38, 166, 70, 198, 102, 230, 22, 150, 54, 182, 86, 214, 118, 246,
  14, 142, 46, 174, 78, 206, 110, 238, 30, 158, 62, 190, 94, 222, 126, 254,
  1,  129, 33, 161, 65, 193, 97,  225, 17, 145, 49, 177, 81, 209, 113, 241,
  9,  137, 41, 169, 73, 201, 105, 233, 25, 153, 57, 185, 89, 217, 121, 249,
  3,  131, 35, 163, 67, 195, 99,  227, 19, 147, 51, 179, 83, 211, 115, 243,
  11, 139, 43, 171, 75, 203, 107, 235, 27, 155, 59, 187, 91, 219, 123, 251,
  5,  133, 37, 165, 69, 197, 101, 229, 21, 149, 53, 181, 85, 213, 117, 245,
  13, 141, 45, 173, 77, 205, 109, 237, 29, 157, 61, 189, 93, 221, 125, 253,
  7,  135, 39, 167, 71, 199, 103, 231, 23, 151, 55, 183, 87, 215, 119, 247,
  15, 143, 47, 175, 79, 207, 111, 239, 31, 159, 63, 191, 95, 223, 127, 255};

const uint8_t pcc::kOccMapRotateX090Y180[256] = {
  0,  16, 64, 80, 32, 48, 96,  112, 128, 144, 192, 208, 160, 176, 224, 240,
  1,  17, 65, 81, 33, 49, 97,  113, 129, 145, 193, 209, 161, 177, 225, 241,
  4,  20, 68, 84, 36, 52, 100, 116, 132, 148, 196, 212, 164, 180, 228, 244,
  5,  21, 69, 85, 37, 53, 101, 117, 133, 149, 197, 213, 165, 181, 229, 245,
  2,  18, 66, 82, 34, 50, 98,  114, 130, 146, 194, 210, 162, 178, 226, 242,
  3,  19, 67, 83, 35, 51, 99,  115, 131, 147, 195, 211, 163, 179, 227, 243,
  6,  22, 70, 86, 38, 54, 102, 118, 134, 150, 198, 214, 166, 182, 230, 246,
  7,  23, 71, 87, 39, 55, 103, 119, 135, 151, 199, 215, 167, 183, 231, 247,
  8,  24, 72, 88, 40, 56, 104, 120, 136, 152, 200, 216, 168, 184, 232, 248,
  9,  25, 73, 89, 41, 57, 105, 121, 137, 153, 201, 217, 169, 185, 233, 249,
  12, 28, 76, 92, 44, 60, 108, 124, 140, 156, 204, 220, 172, 188, 236, 252,
  13, 29, 77, 93, 45, 61, 109, 125, 141, 157, 205, 221, 173, 189, 237, 253,
  10, 26, 74, 90, 42, 58, 106, 122, 138, 154, 202, 218, 170, 186, 234, 250,
  11, 27, 75, 91, 43, 59, 107, 123, 139, 155, 203, 219, 171, 187, 235, 251,
  14, 30, 78, 94, 46, 62, 110, 126, 142, 158, 206, 222, 174, 190, 238, 254,
  15, 31, 79, 95, 47, 63, 111, 127, 143, 159, 207, 223, 175, 191, 239, 255};

const uint8_t pcc::kOccMapRotateY090[256] = {
  0,   16,  1,   17,  64,  80,  65,  81,  4,   20,  5,   21,  68,  84,  69,
  85,  32,  48,  33,  49,  96,  112, 97,  113, 36,  52,  37,  53,  100, 116,
  101, 117, 2,   18,  3,   19,  66,  82,  67,  83,  6,   22,  7,   23,  70,
  86,  71,  87,  34,  50,  35,  51,  98,  114, 99,  115, 38,  54,  39,  55,
  102, 118, 103, 119, 128, 144, 129, 145, 192, 208, 193, 209, 132, 148, 133,
  149, 196, 212, 197, 213, 160, 176, 161, 177, 224, 240, 225, 241, 164, 180,
  165, 181, 228, 244, 229, 245, 130, 146, 131, 147, 194, 210, 195, 211, 134,
  150, 135, 151, 198, 214, 199, 215, 162, 178, 163, 179, 226, 242, 227, 243,
  166, 182, 167, 183, 230, 246, 231, 247, 8,   24,  9,   25,  72,  88,  73,
  89,  12,  28,  13,  29,  76,  92,  77,  93,  40,  56,  41,  57,  104, 120,
  105, 121, 44,  60,  45,  61,  108, 124, 109, 125, 10,  26,  11,  27,  74,
  90,  75,  91,  14,  30,  15,  31,  78,  94,  79,  95,  42,  58,  43,  59,
  106, 122, 107, 123, 46,  62,  47,  63,  110, 126, 111, 127, 136, 152, 137,
  153, 200, 216, 201, 217, 140, 156, 141, 157, 204, 220, 205, 221, 168, 184,
  169, 185, 232, 248, 233, 249, 172, 188, 173, 189, 236, 252, 237, 253, 138,
  154, 139, 155, 202, 218, 203, 219, 142, 158, 143, 159, 206, 222, 207, 223,
  170, 186, 171, 187, 234, 250, 235, 251, 174, 190, 175, 191, 238, 254, 239,
  255};

const uint8_t pcc::kOccMapRotateY270[256] = {
  0,  2,  32,  34,  8,  10, 40,  42,  128, 130, 160, 162, 136, 138, 168, 170,
  1,  3,  33,  35,  9,  11, 41,  43,  129, 131, 161, 163, 137, 139, 169, 171,
  16, 18, 48,  50,  24, 26, 56,  58,  144, 146, 176, 178, 152, 154, 184, 186,
  17, 19, 49,  51,  25, 27, 57,  59,  145, 147, 177, 179, 153, 155, 185, 187,
  4,  6,  36,  38,  12, 14, 44,  46,  132, 134, 164, 166, 140, 142, 172, 174,
  5,  7,  37,  39,  13, 15, 45,  47,  133, 135, 165, 167, 141, 143, 173, 175,
  20, 22, 52,  54,  28, 30, 60,  62,  148, 150, 180, 182, 156, 158, 188, 190,
  21, 23, 53,  55,  29, 31, 61,  63,  149, 151, 181, 183, 157, 159, 189, 191,
  64, 66, 96,  98,  72, 74, 104, 106, 192, 194, 224, 226, 200, 202, 232, 234,
  65, 67, 97,  99,  73, 75, 105, 107, 193, 195, 225, 227, 201, 203, 233, 235,
  80, 82, 112, 114, 88, 90, 120, 122, 208, 210, 240, 242, 216, 218, 248, 250,
  81, 83, 113, 115, 89, 91, 121, 123, 209, 211, 241, 243, 217, 219, 249, 251,
  68, 70, 100, 102, 76, 78, 108, 110, 196, 198, 228, 230, 204, 206, 236, 238,
  69, 71, 101, 103, 77, 79, 109, 111, 197, 199, 229, 231, 205, 207, 237, 239,
  84, 86, 116, 118, 92, 94, 124, 126, 212, 214, 244, 246, 220, 222, 252, 254,
  85, 87, 117, 119, 93, 95, 125, 127, 213, 215, 245, 247, 221, 223, 253, 255};

const uint8_t pcc::kOccMapMirrorXY[256] = {
  0,   2,   1,   3,   8,   10,  9,   11,  4,   6,   5,   7,   12,  14,  13,
  15,  32,  34,  33,  35,  40,  42,  41,  43,  36,  38,  37,  39,  44,  46,
  45,  47,  16,  18,  17,  19,  24,  26,  25,  27,  20,  22,  21,  23,  28,
  30,  29,  31,  48,  50,  49,  51,  56,  58,  57,  59,  52,  54,  53,  55,
  60,  62,  61,  63,  128, 130, 129, 131, 136, 138, 137, 139, 132, 134, 133,
  135, 140, 142, 141, 143, 160, 162, 161, 163, 168, 170, 169, 171, 164, 166,
  165, 167, 172, 174, 173, 175, 144, 146, 145, 147, 152, 154, 153, 155, 148,
  150, 149, 151, 156, 158, 157, 159, 176, 178, 177, 179, 184, 186, 185, 187,
  180, 182, 181, 183, 188, 190, 189, 191, 64,  66,  65,  67,  72,  74,  73,
  75,  68,  70,  69,  71,  76,  78,  77,  79,  96,  98,  97,  99,  104, 106,
  105, 107, 100, 102, 101, 103, 108, 110, 109, 111, 80,  82,  81,  83,  88,
  90,  89,  91,  84,  86,  85,  87,  92,  94,  93,  95,  112, 114, 113, 115,
  120, 122, 121, 123, 116, 118, 117, 119, 124, 126, 125, 127, 192, 194, 193,
  195, 200, 202, 201, 203, 196, 198, 197, 199, 204, 206, 205, 207, 224, 226,
  225, 227, 232, 234, 233, 235, 228, 230, 229, 231, 236, 238, 237, 239, 208,
  210, 209, 211, 216, 218, 217, 219, 212, 214, 213, 215, 220, 222, 221, 223,
  240, 242, 241, 243, 248, 250, 249, 251, 244, 246, 245, 247, 252, 254, 253,
  255};

const uint8_t pcc::kOccMapRotateZ270[256] = {
  0,   16,  32,  48,  1,   17,  33,  49,  2,   18,  34,  50,  3,   19,  35,
  51,  64,  80,  96,  112, 65,  81,  97,  113, 66,  82,  98,  114, 67,  83,
  99,  115, 128, 144, 160, 176, 129, 145, 161, 177, 130, 146, 162, 178, 131,
  147, 163, 179, 192, 208, 224, 240, 193, 209, 225, 241, 194, 210, 226, 242,
  195, 211, 227, 243, 4,   20,  36,  52,  5,   21,  37,  53,  6,   22,  38,
  54,  7,   23,  39,  55,  68,  84,  100, 116, 69,  85,  101, 117, 70,  86,
  102, 118, 71,  87,  103, 119, 132, 148, 164, 180, 133, 149, 165, 181, 134,
  150, 166, 182, 135, 151, 167, 183, 196, 212, 228, 244, 197, 213, 229, 245,
  198, 214, 230, 246, 199, 215, 231, 247, 8,   24,  40,  56,  9,   25,  41,
  57,  10,  26,  42,  58,  11,  27,  43,  59,  72,  88,  104, 120, 73,  89,
  105, 121, 74,  90,  106, 122, 75,  91,  107, 123, 136, 152, 168, 184, 137,
  153, 169, 185, 138, 154, 170, 186, 139, 155, 171, 187, 200, 216, 232, 248,
  201, 217, 233, 249, 202, 218, 234, 250, 203, 219, 235, 251, 12,  28,  44,
  60,  13,  29,  45,  61,  14,  30,  46,  62,  15,  31,  47,  63,  76,  92,
  108, 124, 77,  93,  109, 125, 78,  94,  110, 126, 79,  95,  111, 127, 140,
  156, 172, 188, 141, 157, 173, 189, 142, 158, 174, 190, 143, 159, 175, 191,
  204, 220, 236, 252, 205, 221, 237, 253, 206, 222, 238, 254, 207, 223, 239,
  255};

const uint8_t pcc::kOccMapRotateZ180[256] = {
  0,  64, 128, 192, 16, 80, 144, 208, 32, 96,  160, 224, 48, 112, 176, 240,
  4,  68, 132, 196, 20, 84, 148, 212, 36, 100, 164, 228, 52, 116, 180, 244,
  8,  72, 136, 200, 24, 88, 152, 216, 40, 104, 168, 232, 56, 120, 184, 248,
  12, 76, 140, 204, 28, 92, 156, 220, 44, 108, 172, 236, 60, 124, 188, 252,
  1,  65, 129, 193, 17, 81, 145, 209, 33, 97,  161, 225, 49, 113, 177, 241,
  5,  69, 133, 197, 21, 85, 149, 213, 37, 101, 165, 229, 53, 117, 181, 245,
  9,  73, 137, 201, 25, 89, 153, 217, 41, 105, 169, 233, 57, 121, 185, 249,
  13, 77, 141, 205, 29, 93, 157, 221, 45, 109, 173, 237, 61, 125, 189, 253,
  2,  66, 130, 194, 18, 82, 146, 210, 34, 98,  162, 226, 50, 114, 178, 242,
  6,  70, 134, 198, 22, 86, 150, 214, 38, 102, 166, 230, 54, 118, 182, 246,
  10, 74, 138, 202, 26, 90, 154, 218, 42, 106, 170, 234, 58, 122, 186, 250,
  14, 78, 142, 206, 30, 94, 158, 222, 46, 110, 174, 238, 62, 126, 190, 254,
  3,  67, 131, 195, 19, 83, 147, 211, 35, 99,  163, 227, 51, 115, 179, 243,
  7,  71, 135, 199, 23, 87, 151, 215, 39, 103, 167, 231, 55, 119, 183, 247,
  11, 75, 139, 203, 27, 91, 155, 219, 43, 107, 171, 235, 59, 123, 187, 251,
  15, 79, 143, 207, 31, 95, 159, 223, 47, 111, 175, 239, 63, 127, 191, 255};

const uint8_t pcc::kOccMapRotateZ090[256] = {
  0,  4,  8,  12, 64,  68,  72,  76,  128, 132, 136, 140, 192, 196, 200, 204,
  1,  5,  9,  13, 65,  69,  73,  77,  129, 133, 137, 141, 193, 197, 201, 205,
  2,  6,  10, 14, 66,  70,  74,  78,  130, 134, 138, 142, 194, 198, 202, 206,
  3,  7,  11, 15, 67,  71,  75,  79,  131, 135, 139, 143, 195, 199, 203, 207,
  16, 20, 24, 28, 80,  84,  88,  92,  144, 148, 152, 156, 208, 212, 216, 220,
  17, 21, 25, 29, 81,  85,  89,  93,  145, 149, 153, 157, 209, 213, 217, 221,
  18, 22, 26, 30, 82,  86,  90,  94,  146, 150, 154, 158, 210, 214, 218, 222,
  19, 23, 27, 31, 83,  87,  91,  95,  147, 151, 155, 159, 211, 215, 219, 223,
  32, 36, 40, 44, 96,  100, 104, 108, 160, 164, 168, 172, 224, 228, 232, 236,
  33, 37, 41, 45, 97,  101, 105, 109, 161, 165, 169, 173, 225, 229, 233, 237,
  34, 38, 42, 46, 98,  102, 106, 110, 162, 166, 170, 174, 226, 230, 234, 238,
  35, 39, 43, 47, 99,  103, 107, 111, 163, 167, 171, 175, 227, 231, 235, 239,
  48, 52, 56, 60, 112, 116, 120, 124, 176, 180, 184, 188, 240, 244, 248, 252,
  49, 53, 57, 61, 113, 117, 121, 125, 177, 181, 185, 189, 241, 245, 249, 253,
  50, 54, 58, 62, 114, 118, 122, 126, 178, 182, 186, 190, 242, 246, 250, 254,
  51, 55, 59, 63, 115, 119, 123, 127, 179, 183, 187, 191, 243, 247, 251, 255};

//============================================================================

// transition tables L= 10; 256 ctx
const uint8_t pcc::kCtxMapOctreeOccupancyEvolutionOn0[256] = {
  0,   1,   2,   3,   4,   5,   6,   7,   8,   9,   9,   10,  11,  12,  13,
  14,  15,  16,  17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,
  29,  29,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,
  43,  44,  45,  46,  47,  48,  48,  49,  50,  51,  52,  53,  54,  55,  56,
  57,  58,  59,  60,  61,  62,  63,  64,  65,  66,  66,  67,  68,  69,  70,
  71,  72,  73,  74,  75,  76,  77,  78,  79,  80,  81,  82,  82,  83,  84,
  85,  86,  87,  88,  89,  90,  91,  92,  93,  94,  95,  96,  96,  97,  98,
  99,  100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 109, 110, 111, 112,
  113, 114, 115, 116, 117, 118, 119, 120, 120, 121, 122, 123, 124, 125, 126,
  127, 128, 129, 130, 130, 131, 132, 133, 134, 135, 136, 137, 138, 138, 139,
  140, 141, 142, 143, 144, 145, 146, 146, 147, 148, 149, 150, 151, 152, 152,
  153, 154, 155, 156, 157, 158, 158, 159, 160, 161, 162, 163, 163, 164, 165,
  166, 167, 167, 168, 169, 170, 171, 171, 172, 173, 174, 175, 175, 176, 177,
  178, 178, 179, 180, 181, 181, 182, 183, 183, 184, 185, 185, 186, 187, 187,
  188, 189, 189, 190, 191, 191, 192, 192, 193, 194, 194, 195, 195, 196, 196,
  197, 197, 198, 198, 199, 199, 199, 200, 200, 201, 201, 201, 202, 202, 202,
  203, 203, 203, 203, 204, 204, 204, 204, 204, 204, 204, 205, 205, 205, 205,
  205};

const uint8_t pcc::kCtxMapOctreeOccupancyEvolutionOn1[256] = {
  49,  49,  49,  49,  49,  49,  49,  49,  49,  50,  50,  50,  50,  51,  51,
  51,  51,  52,  52,  52,  53,  53,  54,  54,  54,  55,  55,  56,  56,  57,
  57,  58,  58,  59,  59,  60,  61,  61,  62,  62,  63,  64,  64,  65,  66,
  66,  67,  68,  68,  69,  70,  70,  71,  72,  72,  73,  74,  75,  75,  76,
  77,  78,  78,  79,  80,  81,  82,  82,  83,  84,  85,  86,  86,  87,  88,
  89,  90,  90,  91,  92,  93,  94,  95,  95,  96,  97,  98,  99,  100, 101,
  101, 102, 103, 104, 105, 106, 107, 107, 108, 109, 110, 111, 112, 113, 114,
  115, 115, 116, 117, 118, 119, 120, 121, 122, 123, 123, 124, 125, 126, 127,
  128, 129, 130, 131, 132, 133, 133, 134, 135, 136, 137, 138, 139, 140, 141,
  142, 143, 144, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155,
  156, 157, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169,
  170, 171, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183,
  184, 185, 186, 187, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197,
  198, 199, 200, 201, 202, 203, 204, 205, 205, 206, 207, 208, 209, 210, 211,
  212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 224, 225,
  226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240,
  241, 242, 243, 244, 245, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254,
  255};

//============================================================================

const uint8_t pcc::kDualLutOccupancyCoderInit[10][32] = {
  /* [0] = */ {1,   2,   3,   4,   5,   7,   8,   10,  12,  16, 17,
               20,  32,  42,  44,  48,  49,  64,  65,  68,  76, 80,
               115, 128, 136, 160, 163, 168, 192, 207, 224, 240},
  /* [1] = */ {1,  2,  3,   4,   5,   6,   8,   10,  12,  14, 15,
               16, 17, 19,  20,  21,  32,  34,  48,  64,  68, 69,
               84, 85, 128, 136, 160, 162, 192, 194, 243, 250},
  /* [2] = */ {1,   3,   4,   8,   10,  11,  13,  16,  18,  20, 21,
               29,  40,  42,  48,  51,  64,  80,  81,  84,  85, 124,
               128, 138, 142, 144, 163, 168, 192, 197, 213, 215},
  /* [3] = */ {1,   2,   3,   4,   5,   11,  16,  17,  21,  34,  42,
               43,  64,  68,  69,  76,  80,  81,  84,  85,  128, 136,
               138, 142, 160, 168, 208, 212, 221, 223, 232, 245},
  /* [4] = */ {2,   8,   10,  17,  19,  21,  32,  34,  42,  47,  48,
               63,  64,  65,  69,  81,  84,  85,  87,  113, 136, 138,
               160, 162, 168, 170, 234, 243, 245, 248, 250, 255},
  /* [5] = */ {16,  29,  31,  64,  69,  80,  84,  85,  87,  93,  95,
               102, 113, 117, 119, 125, 127, 168, 170, 174, 213, 215,
               221, 223, 239, 242, 245, 247, 252, 253, 254, 255},
  /* [6] = */ {1,   2,   3,   4,   5,   8,   10,  16,  17,  21, 31,
               32,  34,  42,  48,  64,  68,  69,  80,  81,  84, 85,
               128, 130, 136, 138, 160, 162, 168, 170, 212, 253},
  /* [7] = */ {1,  2,   4,   5,   8,   10,  11,  16,  20,  21, 23,
               34, 42,  43,  53,  63,  64,  69,  77,  81,  84, 85,
               93, 128, 130, 138, 160, 168, 170, 171, 213, 247},
  /* [8] = */ {5,   17,  21,  23,  29,  42,  64,  65,  69,  80,  84,
               85,  87,  93,  126, 127, 130, 136, 138, 151, 162, 168,
               170, 171, 190, 205, 213, 223, 239, 252, 253, 254},
  /* [9] = */ {16,  19,  63,  68,  93,  95,  96,  118, 119, 126, 127,
               162, 174, 187, 191, 197, 216, 217, 219, 221, 223, 236,
               239, 245, 246, 247, 250, 251, 252, 253, 254, 255}};

//============================================================================

const uint32_t pcc::kMortonCode256X[256] = {
  0x00000000, 0x00000001, 0x00000008, 0x00000009, 0x00000040, 0x00000041,
  0x00000048, 0x00000049, 0x00000200, 0x00000201, 0x00000208, 0x00000209,
  0x00000240, 0x00000241, 0x00000248, 0x00000249, 0x00001000, 0x00001001,
  0x00001008, 0x00001009, 0x00001040, 0x00001041, 0x00001048, 0x00001049,
  0x00001200, 0x00001201, 0x00001208, 0x00001209, 0x00001240, 0x00001241,
  0x00001248, 0x00001249, 0x00008000, 0x00008001, 0x00008008, 0x00008009,
  0x00008040, 0x00008041, 0x00008048, 0x00008049, 0x00008200, 0x00008201,
  0x00008208, 0x00008209, 0x00008240, 0x00008241, 0x00008248, 0x00008249,
  0x00009000, 0x00009001, 0x00009008, 0x00009009, 0x00009040, 0x00009041,
  0x00009048, 0x00009049, 0x00009200, 0x00009201, 0x00009208, 0x00009209,
  0x00009240, 0x00009241, 0x00009248, 0x00009249, 0x00040000, 0x00040001,
  0x00040008, 0x00040009, 0x00040040, 0x00040041, 0x00040048, 0x00040049,
  0x00040200, 0x00040201, 0x00040208, 0x00040209, 0x00040240, 0x00040241,
  0x00040248, 0x00040249, 0x00041000, 0x00041001, 0x00041008, 0x00041009,
  0x00041040, 0x00041041, 0x00041048, 0x00041049, 0x00041200, 0x00041201,
  0x00041208, 0x00041209, 0x00041240, 0x00041241, 0x00041248, 0x00041249,
  0x00048000, 0x00048001, 0x00048008, 0x00048009, 0x00048040, 0x00048041,
  0x00048048, 0x00048049, 0x00048200, 0x00048201, 0x00048208, 0x00048209,
  0x00048240, 0x00048241, 0x00048248, 0x00048249, 0x00049000, 0x00049001,
  0x00049008, 0x00049009, 0x00049040, 0x00049041, 0x00049048, 0x00049049,
  0x00049200, 0x00049201, 0x00049208, 0x00049209, 0x00049240, 0x00049241,
  0x00049248, 0x00049249, 0x00200000, 0x00200001, 0x00200008, 0x00200009,
  0x00200040, 0x00200041, 0x00200048, 0x00200049, 0x00200200, 0x00200201,
  0x00200208, 0x00200209, 0x00200240, 0x00200241, 0x00200248, 0x00200249,
  0x00201000, 0x00201001, 0x00201008, 0x00201009, 0x00201040, 0x00201041,
  0x00201048, 0x00201049, 0x00201200, 0x00201201, 0x00201208, 0x00201209,
  0x00201240, 0x00201241, 0x00201248, 0x00201249, 0x00208000, 0x00208001,
  0x00208008, 0x00208009, 0x00208040, 0x00208041, 0x00208048, 0x00208049,
  0x00208200, 0x00208201, 0x00208208, 0x00208209, 0x00208240, 0x00208241,
  0x00208248, 0x00208249, 0x00209000, 0x00209001, 0x00209008, 0x00209009,
  0x00209040, 0x00209041, 0x00209048, 0x00209049, 0x00209200, 0x00209201,
  0x00209208, 0x00209209, 0x00209240, 0x00209241, 0x00209248, 0x00209249,
  0x00240000, 0x00240001, 0x00240008, 0x00240009, 0x00240040, 0x00240041,
  0x00240048, 0x00240049, 0x00240200, 0x00240201, 0x00240208, 0x00240209,
  0x00240240, 0x00240241, 0x00240248, 0x00240249, 0x00241000, 0x00241001,
  0x00241008, 0x00241009, 0x00241040, 0x00241041, 0x00241048, 0x00241049,
  0x00241200, 0x00241201, 0x00241208, 0x00241209, 0x00241240, 0x00241241,
  0x00241248, 0x00241249, 0x00248000, 0x00248001, 0x00248008, 0x00248009,
  0x00248040, 0x00248041, 0x00248048, 0x00248049, 0x00248200, 0x00248201,
  0x00248208, 0x00248209, 0x00248240, 0x00248241, 0x00248248, 0x00248249,
  0x00249000, 0x00249001, 0x00249008, 0x00249009, 0x00249040, 0x00249041,
  0x00249048, 0x00249049, 0x00249200, 0x00249201, 0x00249208, 0x00249209,
  0x00249240, 0x00249241, 0x00249248, 0x00249249};

const uint32_t pcc::kMortonCode256Y[256] = {
  0x00000000, 0x00000002, 0x00000010, 0x00000012, 0x00000080, 0x00000082,
  0x00000090, 0x00000092, 0x00000400, 0x00000402, 0x00000410, 0x00000412,
  0x00000480, 0x00000482, 0x00000490, 0x00000492, 0x00002000, 0x00002002,
  0x00002010, 0x00002012, 0x00002080, 0x00002082, 0x00002090, 0x00002092,
  0x00002400, 0x00002402, 0x00002410, 0x00002412, 0x00002480, 0x00002482,
  0x00002490, 0x00002492, 0x00010000, 0x00010002, 0x00010010, 0x00010012,
  0x00010080, 0x00010082, 0x00010090, 0x00010092, 0x00010400, 0x00010402,
  0x00010410, 0x00010412, 0x00010480, 0x00010482, 0x00010490, 0x00010492,
  0x00012000, 0x00012002, 0x00012010, 0x00012012, 0x00012080, 0x00012082,
  0x00012090, 0x00012092, 0x00012400, 0x00012402, 0x00012410, 0x00012412,
  0x00012480, 0x00012482, 0x00012490, 0x00012492, 0x00080000, 0x00080002,
  0x00080010, 0x00080012, 0x00080080, 0x00080082, 0x00080090, 0x00080092,
  0x00080400, 0x00080402, 0x00080410, 0x00080412, 0x00080480, 0x00080482,
  0x00080490, 0x00080492, 0x00082000, 0x00082002, 0x00082010, 0x00082012,
  0x00082080, 0x00082082, 0x00082090, 0x00082092, 0x00082400, 0x00082402,
  0x00082410, 0x00082412, 0x00082480, 0x00082482, 0x00082490, 0x00082492,
  0x00090000, 0x00090002, 0x00090010, 0x00090012, 0x00090080, 0x00090082,
  0x00090090, 0x00090092, 0x00090400, 0x00090402, 0x00090410, 0x00090412,
  0x00090480, 0x00090482, 0x00090490, 0x00090492, 0x00092000, 0x00092002,
  0x00092010, 0x00092012, 0x00092080, 0x00092082, 0x00092090, 0x00092092,
  0x00092400, 0x00092402, 0x00092410, 0x00092412, 0x00092480, 0x00092482,
  0x00092490, 0x00092492, 0x00400000, 0x00400002, 0x00400010, 0x00400012,
  0x00400080, 0x00400082, 0x00400090, 0x00400092, 0x00400400, 0x00400402,
  0x00400410, 0x00400412, 0x00400480, 0x00400482, 0x00400490, 0x00400492,
  0x00402000, 0x00402002, 0x00402010, 0x00402012, 0x00402080, 0x00402082,
  0x00402090, 0x00402092, 0x00402400, 0x00402402, 0x00402410, 0x00402412,
  0x00402480, 0x00402482, 0x00402490, 0x00402492, 0x00410000, 0x00410002,
  0x00410010, 0x00410012, 0x00410080, 0x00410082, 0x00410090, 0x00410092,
  0x00410400, 0x00410402, 0x00410410, 0x00410412, 0x00410480, 0x00410482,
  0x00410490, 0x00410492, 0x00412000, 0x00412002, 0x00412010, 0x00412012,
  0x00412080, 0x00412082, 0x00412090, 0x00412092, 0x00412400, 0x00412402,
  0x00412410, 0x00412412, 0x00412480, 0x00412482, 0x00412490, 0x00412492,
  0x00480000, 0x00480002, 0x00480010, 0x00480012, 0x00480080, 0x00480082,
  0x00480090, 0x00480092, 0x00480400, 0x00480402, 0x00480410, 0x00480412,
  0x00480480, 0x00480482, 0x00480490, 0x00480492, 0x00482000, 0x00482002,
  0x00482010, 0x00482012, 0x00482080, 0x00482082, 0x00482090, 0x00482092,
  0x00482400, 0x00482402, 0x00482410, 0x00482412, 0x00482480, 0x00482482,
  0x00482490, 0x00482492, 0x00490000, 0x00490002, 0x00490010, 0x00490012,
  0x00490080, 0x00490082, 0x00490090, 0x00490092, 0x00490400, 0x00490402,
  0x00490410, 0x00490412, 0x00490480, 0x00490482, 0x00490490, 0x00490492,
  0x00492000, 0x00492002, 0x00492010, 0x00492012, 0x00492080, 0x00492082,
  0x00492090, 0x00492092, 0x00492400, 0x00492402, 0x00492410, 0x00492412,
  0x00492480, 0x00492482, 0x00492490, 0x00492492};

const uint32_t pcc::kMortonCode256Z[256] = {
  0x00000000, 0x00000004, 0x00000020, 0x00000024, 0x00000100, 0x00000104,
  0x00000120, 0x00000124, 0x00000800, 0x00000804, 0x00000820, 0x00000824,
  0x00000900, 0x00000904, 0x00000920, 0x00000924, 0x00004000, 0x00004004,
  0x00004020, 0x00004024, 0x00004100, 0x00004104, 0x00004120, 0x00004124,
  0x00004800, 0x00004804, 0x00004820, 0x00004824, 0x00004900, 0x00004904,
  0x00004920, 0x00004924, 0x00020000, 0x00020004, 0x00020020, 0x00020024,
  0x00020100, 0x00020104, 0x00020120, 0x00020124, 0x00020800, 0x00020804,
  0x00020820, 0x00020824, 0x00020900, 0x00020904, 0x00020920, 0x00020924,
  0x00024000, 0x00024004, 0x00024020, 0x00024024, 0x00024100, 0x00024104,
  0x00024120, 0x00024124, 0x00024800, 0x00024804, 0x00024820, 0x00024824,
  0x00024900, 0x00024904, 0x00024920, 0x00024924, 0x00100000, 0x00100004,
  0x00100020, 0x00100024, 0x00100100, 0x00100104, 0x00100120, 0x00100124,
  0x00100800, 0x00100804, 0x00100820, 0x00100824, 0x00100900, 0x00100904,
  0x00100920, 0x00100924, 0x00104000, 0x00104004, 0x00104020, 0x00104024,
  0x00104100, 0x00104104, 0x00104120, 0x00104124, 0x00104800, 0x00104804,
  0x00104820, 0x00104824, 0x00104900, 0x00104904, 0x00104920, 0x00104924,
  0x00120000, 0x00120004, 0x00120020, 0x00120024, 0x00120100, 0x00120104,
  0x00120120, 0x00120124, 0x00120800, 0x00120804, 0x00120820, 0x00120824,
  0x00120900, 0x00120904, 0x00120920, 0x00120924, 0x00124000, 0x00124004,
  0x00124020, 0x00124024, 0x00124100, 0x00124104, 0x00124120, 0x00124124,
  0x00124800, 0x00124804, 0x00124820, 0x00124824, 0x00124900, 0x00124904,
  0x00124920, 0x00124924, 0x00800000, 0x00800004, 0x00800020, 0x00800024,
  0x00800100, 0x00800104, 0x00800120, 0x00800124, 0x00800800, 0x00800804,
  0x00800820, 0x00800824, 0x00800900, 0x00800904, 0x00800920, 0x00800924,
  0x00804000, 0x00804004, 0x00804020, 0x00804024, 0x00804100, 0x00804104,
  0x00804120, 0x00804124, 0x00804800, 0x00804804, 0x00804820, 0x00804824,
  0x00804900, 0x00804904, 0x00804920, 0x00804924, 0x00820000, 0x00820004,
  0x00820020, 0x00820024, 0x00820100, 0x00820104, 0x00820120, 0x00820124,
  0x00820800, 0x00820804, 0x00820820, 0x00820824, 0x00820900, 0x00820904,
  0x00820920, 0x00820924, 0x00824000, 0x00824004, 0x00824020, 0x00824024,
  0x00824100, 0x00824104, 0x00824120, 0x00824124, 0x00824800, 0x00824804,
  0x00824820, 0x00824824, 0x00824900, 0x00824904, 0x00824920, 0x00824924,
  0x00900000, 0x00900004, 0x00900020, 0x00900024, 0x00900100, 0x00900104,
  0x00900120, 0x00900124, 0x00900800, 0x00900804, 0x00900820, 0x00900824,
  0x00900900, 0x00900904, 0x00900920, 0x00900924, 0x00904000, 0x00904004,
  0x00904020, 0x00904024, 0x00904100, 0x00904104, 0x00904120, 0x00904124,
  0x00904800, 0x00904804, 0x00904820, 0x00904824, 0x00904900, 0x00904904,
  0x00904920, 0x00904924, 0x00920000, 0x00920004, 0x00920020, 0x00920024,
  0x00920100, 0x00920104, 0x00920120, 0x00920124, 0x00920800, 0x00920804,
  0x00920820, 0x00920824, 0x00920900, 0x00920904, 0x00920920, 0x00920924,
  0x00924000, 0x00924004, 0x00924020, 0x00924024, 0x00924100, 0x00924104,
  0x00924120, 0x00924124, 0x00924800, 0x00924804, 0x00924820, 0x00924824,
  0x00924900, 0x00924904, 0x00924920, 0x00924924};

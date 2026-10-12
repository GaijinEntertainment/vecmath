//
// Dagor Engine 6.5 - 1st party libs
// Copyright (C) Gaijin Games KFT.  All rights reserved.
//
#pragma once

// x86 F16C half<->float conversions, split out of dag_vecMath_pc_sse.h because _mm_cvtph_ps and
// _mm_cvtps_ph are only declared by <immintrin.h>, the umbrella intrinsic header: it costs about
// 90ms of parse per translation unit against the <smmintrin.h> an SSE4.1 target otherwise needs.
// A target built with F16C (or AVX2) sets _TARGET_HAS_FC16 and pulls this in through
// dag_vecMath.h, so v_float_to_half and friends keep using the hardware path. A plain SSE4.1 PC
// build has no F16C at compile time and converts through the v_sw_* emulation instead; code that
// dispatches to the hardware path after a cpuid check includes this header itself and pays for
// <immintrin.h> only there.

#include <vecmath/dag_vecMathDecl.h>

#if _TARGET_SIMD_SSE

#include <immintrin.h>

VECTORCALL VECMATH_FINLINE vec4f v_fc16_half_to_float_lo(vec4i v)
{
  return _mm_cvtph_ps(v);
}

VECTORCALL VECMATH_FINLINE vec4i v_fc16_float_to_half_rtne_lo(vec4f v)
{
  return _mm_cvtps_ph(v, _MM_FROUND_TO_NEAREST_INT);
}

VECTORCALL VECMATH_FINLINE vec4i v_fc16_float_to_half_down_lo(vec4f v)
{
  return _mm_cvtps_ph(v, _MM_FROUND_TO_NEG_INF);
}

VECTORCALL VECMATH_FINLINE vec4i v_fc16_float_to_half_up_lo(vec4f v)
{
  return _mm_cvtps_ph(v, _MM_FROUND_TO_POS_INF);
}

VECTORCALL VECMATH_FINLINE vec4i v_fc16_float_to_half_trunc_lo(vec4f v)
{
  return _mm_cvtps_ph(v, _MM_FROUND_TO_ZERO);
}

#endif

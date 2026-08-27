#include <assert.h>

#include "dag_vecMath.h"
#include "dag_vecMath_pc_sse.h"
#include "dag_vecMath_trig_constexpr.h"
#include "dag_vecMath_trig.h"

#pragma warning( push )
#pragma warning( disable : 4995 )		// Allow deprecated functions, cos/sin/log etc...
#pragma warning( disable : 4220 )		// Disable warning treated as an error, for this test case only.

#include <cmath>						// Using to check our compile-time math with std implementation.
void test_constexpr_trig() {
#define COMPILETIME_TEST_COMPARE_VALS(a, b, ...)                                                                                                               \
  {                                                                                                                                                            \
    constexpr auto avalue = a(__VA_ARGS__);                                                                                                                    \
    const auto bvalue = b(__VA_ARGS__);                                                                                                                        \
    assert(v_abs_c(avalue - bvalue) < 0.0000001);                                                                                                              \
  }

#define COMPILETIME_TEST_COMPARE_MVALS(a, b, ...)                                                                                                              \
  {                                                                                                                                                            \
    constexpr auto avalue = a(__VA_ARGS__);                                                                                                                    \
    const auto bvalue = b(__VA_ARGS__);                                                                                                                        \
    assert(v_abs_c<float>(avalue.m128_f32[0] - bvalue.m128_f32[0]) < 0.00001);                                                                                 \
    assert(v_abs_c<float>(avalue.m128_f32[1] - bvalue.m128_f32[1]) < 0.00001);                                                                                 \
    assert(v_abs_c<float>(avalue.m128_f32[2] - bvalue.m128_f32[2]) < 0.00001);                                                                                 \
    assert(v_abs_c<float>(avalue.m128_f32[3] - bvalue.m128_f32[3]) < 0.00001);                                                                                 \
  }

  // abs
  COMPILETIME_TEST_COMPARE_VALS(v_abs_c, std::fabs, 0.0);
  COMPILETIME_TEST_COMPARE_VALS(v_abs_c, std::fabs, -0.0);
  COMPILETIME_TEST_COMPARE_VALS(v_abs_c, std::fabs, 1.0);
  COMPILETIME_TEST_COMPARE_VALS(v_abs_c, std::fabs, -1.0);

  // cos
  constexpr vec4f point0{ 0.f, 1.f, 2.f, 3.f };
  constexpr vec4f point1{ -1.5f, -1.0f, 1.0f, 1.5f };
  constexpr vec4f point2{ 0.0f, 0.0001f, 1.0001f, 1.5f };
  constexpr vec4f point3{ 1.5f, -1.5f, 0.0001f, 0.0f };
  
  COMPILETIME_TEST_COMPARE_VALS(v_cos_c, std::cos, -1.5);
  COMPILETIME_TEST_COMPARE_VALS(v_cos_c, std::cos, 0.0);

  // cos vec4f
  auto cos_m = [](vec4f p) { return vec4f{ std::cos(p.m128_f32[0]), std::cos(p.m128_f32[1]), std::cos(p.m128_f32[2]), std::cos(p.m128_f32[3]) }; };
  COMPILETIME_TEST_COMPARE_MVALS(v_cos_m, cos_m, point1);
  COMPILETIME_TEST_COMPARE_MVALS(v_cos_m, cos_m, point0);

  // tan
  
  COMPILETIME_TEST_COMPARE_VALS(v_tan_c, std::tan, 0.0);
  COMPILETIME_TEST_COMPARE_VALS(v_tan_c, std::tan, 0.001);
  COMPILETIME_TEST_COMPARE_VALS(v_tan_c, std::tan, 1.001);
  COMPILETIME_TEST_COMPARE_VALS(v_tan_c, std::tan, 1.5);
  COMPILETIME_TEST_COMPARE_VALS(v_tan_c, std::tan, -1.5);

  // tan vec4f
  auto tan_m = [](vec4f p) { return vec4f{ std::tan(p.m128_f32[0]), std::tan(p.m128_f32[1]), std::tan(p.m128_f32[2]), std::tan(p.m128_f32[3]) }; };
  COMPILETIME_TEST_COMPARE_MVALS(v_tan_m, tan_m, point2);
  COMPILETIME_TEST_COMPARE_MVALS(v_tan_m, tan_m, point3);


  // sin
  COMPILETIME_TEST_COMPARE_VALS(v_sin_c, std::sin, -1.5f);
  COMPILETIME_TEST_COMPARE_VALS(v_sin_c, std::sin, 0.0f);
  COMPILETIME_TEST_COMPARE_VALS(v_sin_c, std::sin, 0.001f);
  COMPILETIME_TEST_COMPARE_VALS(v_sin_c, std::sin, 1.001f);
  COMPILETIME_TEST_COMPARE_VALS(v_sin_c, std::sin, 1.5f);
  COMPILETIME_TEST_COMPARE_VALS(v_sin_c, std::sin, 11.1f);

  // sin vec4f
  auto sin_m = [](vec4f p) { return vec4f{ std::sin(p.m128_f32[0]), std::sin(p.m128_f32[1]), std::sin(p.m128_f32[2]), std::sin(p.m128_f32[3]) }; };
  COMPILETIME_TEST_COMPARE_MVALS(v_sin_m, sin_m, point0);
  COMPILETIME_TEST_COMPARE_MVALS(v_sin_m, sin_m, point1);

  // last check
  // rad: 0.f        1.f         2.f         3.f 
  // val: 0.00000000 0.841471016 0.909297407 0.141120046}
  constexpr vec4f gold_sin = v_sin_m(vec4f{0.f, 1.f, 2.f, 3.f});
  assert(gold_sin.m128_f32[0] < 0.000001 
          && v_abs_c(gold_sin.m128_f32[1] - 0.841471016) < 0.000001
          && v_abs_c(gold_sin.m128_f32[2] - 0.909297407) < 0.000001
          && v_abs_c(gold_sin.m128_f32[3] - 0.141120046) < 0.000001);

  printf("All tests passed: OK");
}
#pragma warning(pop)

#pragma once

constexpr long double const_pi = 3.1415926535897932384626433832795028841972L;
constexpr long double const_half_pi = 1.5707963267948966192313216916397514420986L;

template <typename T> constexpr bool v_isnan(const T x) noexcept { return x != x; }

constexpr float v_infinity() noexcept { return __builtin_huge_val(); }

template <typename T> constexpr bool v_isneginf(const T x) noexcept { return (x == -__builtin_huge_val()); }

template <typename T> constexpr bool v_isposinf(const T x) noexcept { return (x == __builtin_huge_val()); }

template <typename T> constexpr bool v_isinf(const T x) noexcept { return (v_isneginf(x) || v_isposinf(x)); }

template <typename T> constexpr bool v_isfinite(const T x) noexcept { return (!v_isnan(x)) && (!v_isinf(x)); }

template <typename T> constexpr T v_abs_c(const T x) noexcept { return (x == T(0) ? T(0) : (x < T(0)) ? (-x) : x); }
constexpr bool							v_isodd_c(const int64_t x)	noexcept { return (x & 1U) != 0; }

namespace detail_pow_integral {
// forward declaration
template <typename T1, typename T2> constexpr T1 pow_integral_compute(const T1 base, const T2 exp_term) noexcept;

// integral-valued powers using method described in
// https://en.wikipedia.org/wiki/Exponentiation_by_squaring

template <typename T1, typename T2> constexpr T1 pow_integral_compute_recur(const T1 base, const T1 val, const T2 exp_term) noexcept {
  if (exp_term > T2(1)) {
    if (v_isodd_c(exp_term)) {
      return pow_integral_compute_recur(base * base, val * base, exp_term / 2);
    }

    return pow_integral_compute_recur(base * base, val, exp_term / 2);
  } else if (exp_term == T2(1)) {
    return (val * base);
  }

  return val;
}

template <typename T1, typename T2>
constexpr T1 pow_integral_sgn_check(const T1 base, const T2 exp_term) noexcept {
  return (pow_integral_compute_recur(base, T1(1), exp_term));
}

template <typename T1, typename T2> constexpr T1 pow_integral_compute(const T1 base, const T2 exp_term) noexcept {
  if (exp_term == T2(3)) {
    return (base * base * base);
  } else if (exp_term == T2(2)) {
    return (base * base);
  } else if (exp_term == T2(1)) {
    return base;
  } else if (exp_term == T2(0)) {
    return T1(1);
  } else if (exp_term == std::numeric_limits<T2>::min()) {
    return T1(0);
  } else if (exp_term == std::numeric_limits<T2>::max()) {
    return v_infinity();
  } else {
    return pow_integral_sgn_check(base, exp_term);
  }
}

template <typename T1, typename T2>
constexpr T1 _pow_integral(const T1 base, const T2 exp_term) noexcept {
  return pow_integral_compute(base, static_cast<int64_t>(exp_term));
}

} // namespace detail_pow_integral

namespace detail_floor {

template <typename T> constexpr int floor_resid(const T x, const T x_whole) noexcept { return ((x < T(0)) && (x < x_whole)); }

template <typename T> constexpr T floor_int(const T x, const T x_whole) noexcept { return (x_whole - static_cast<T>(floor_resid(x, x_whole))); }

template <typename T> constexpr T _floor(const T x) noexcept {
  return (v_isnan(x) ? std::numeric_limits<T>::quiet_NaN() : !v_isfinite(x) ? x : std::numeric_limits<T>::min() > v_abs_c(x) ? x : floor_int(x, T(static_cast<int64_t>(x))));
}

} // namespace detail_floor

namespace detail_tan {
// this is based on a fourth-order expansion of tan(z) using Bernoulli numbers
template <typename T> constexpr T tan_series_exp_long(const T x) noexcept {
  return (-1 / x                                                      // iter 1
          + (x / 3                                                    // iter 2
             + (detail_pow_integral::_pow_integral(x, 3) / 45         // iter 3
                + (2 * detail_pow_integral::_pow_integral(x, 5) / 945 // iter 4
                   + detail_pow_integral::_pow_integral(x, 7) / 4725) // iter 5
                )));
}

template <typename T> constexpr T tan_series_exp(const T x) noexcept {
  // the value tan(pi/2) is somewhat of a convention;
  // technically the function is not defined at EXACTLY pi/2,
  // but this is floating point pi/2
  // otherwise we use an expansion around pi/2
  return (std::numeric_limits<T>::min() > v_abs_c(x - T(const_half_pi)) ? T(1.633124e+16) : tan_series_exp_long(x - T(const_half_pi)));
}

template <typename T> constexpr T tan_cf_recur(const T xx, const int depth, const int max_depth) noexcept {
  return (depth < max_depth ? T(2 * depth - 1) - xx / tan_cf_recur(xx, depth + 1, max_depth) : T(2 * depth - 1));
}

template <typename T> constexpr T tan_cf_main(const T x) noexcept {
  return (((x > T(1.55) && x < T(1.60)))
              ? tan_series_exp(x)
              : x > T(1.4) ? x / tan_cf_recur(x * x, 1, 45) : x > T(1) ? x / tan_cf_recur(x * x, 1, 35) : x / tan_cf_recur(x * x, 1, 25));
}

template <typename T> constexpr T tan_begin(const T x, const int count = 0) noexcept {
  return (x > T(const_pi) ? count > 1 ? std::numeric_limits<T>::quiet_NaN() : (tan_begin(x - T(const_pi) * detail_floor::_floor(x / T(const_pi)), count + 1))
                          : tan_cf_main(x));
}

template <typename T> constexpr T _tan(const T x) noexcept {
  return (v_isnan(x) ? std::numeric_limits<T>::quiet_NaN() : (std::numeric_limits<T>::min() > v_abs_c(x)) ? T(0) : x < T(0) ? (-tan_begin(-x)) : tan_begin(x));
}
} // namespace detail_tan

namespace detail_cos {

template <typename T> constexpr T cos_impl(const T x) noexcept { return (T(1) - x * x) / (T(1) + x * x); }

template <typename T> constexpr T _cos(const T x) noexcept {
  return (v_isnan(x) ? std::numeric_limits<T>::quiet_NaN()
                    : (std::numeric_limits<T>::min() > v_abs_c(x))
                          ? T(1) // indistinguishable from 0
                          : (std::numeric_limits<T>::min() > v_abs_c(x - T(const_half_pi)))
                                ? T(0) // special cases: pi/2 and pi
                                : (std::numeric_limits<T>::min() > v_abs_c(x + T(const_half_pi)))
                                      ? T(0)
                                      : (std::numeric_limits<T>::min() > v_abs_c(x - T(const_pi)))
                                            ? -T(1)
                                            : (std::numeric_limits<T>::min() > v_abs_c(x + T(const_pi))) ? -T(1) : cos_impl(detail_tan::_tan(x / T(2))));
}
} // namespace detail_cos

namespace detail_sin {
template <typename T> constexpr T sin_compute(const T x) noexcept { return T(2) * x / (T(1) + x * x); }

template <typename T> constexpr T _sin(const T x) noexcept {
  return (v_isnan(x) ? std::numeric_limits<T>::quiet_NaN()
                    : std::numeric_limits<T>::min() > v_abs_c(x)
                          ? T(0)
                          : std::numeric_limits<T>::min() > v_abs_c(x - T(const_half_pi)) // special cases: pi/2 and pi
                                ? T(1)
                                : std::numeric_limits<T>::min() > v_abs_c(x + T(const_half_pi))
                                      ? -T(1)
                                      : std::numeric_limits<T>::min() > v_abs_c(x - T(const_pi))
                                            ? T(0)
                                            : std::numeric_limits<T>::min() > v_abs_c(x + T(const_pi)) ? -T(0) : sin_compute(detail_tan::_tan(x / T(2))));
}

} // namespace detail_sin

template <typename T> constexpr T v_tan_c(const T x) noexcept { return detail_tan::_tan(static_cast<T>(x)); }
                      constexpr vec4f v_tan_m(const vec4f a) noexcept { return vec4f{detail_tan::_tan(a.m128_f32[0]),
                                                                                     detail_tan::_tan(a.m128_f32[1]),
                                                                                     detail_tan::_tan(a.m128_f32[2]),
                                                                                     detail_tan::_tan(a.m128_f32[3])}; }

template <typename T> constexpr T v_cos_c(const T x) noexcept { return detail_cos::_cos(static_cast<T>(x)); }
                      constexpr vec4f v_cos_m(const vec4f a) noexcept { return vec4f{detail_cos::_cos(a.m128_f32[0]),
                                                                                     detail_cos::_cos(a.m128_f32[1]),
                                                                                     detail_cos::_cos(a.m128_f32[2]),
                                                                                     detail_cos::_cos(a.m128_f32[3])}; }

template <typename T> constexpr T v_sin_c(const T x) noexcept { return detail_sin::_sin(static_cast<T>(x)); }
                      constexpr vec4f v_sin_m(const vec4f a) noexcept { return vec4f{detail_sin::_sin(a.m128_f32[0]),
                                                                                     detail_sin::_sin(a.m128_f32[1]),
                                                                                     detail_sin::_sin(a.m128_f32[2]),
                                                                                     detail_sin::_sin(a.m128_f32[3])}; }
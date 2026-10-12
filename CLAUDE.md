# vecmath: maintainer guide

Platform-abstracted SIMD C API: SSE2/SSSE3/SSE4.1 (x86) and NEON (ARM), with a scalar per-lane
backend for no-SIMD targets. Engine code outside vecmath still assumes SSE or NEON; validate
scalar changes with `vecMathUnitTest` (`jam -sVecMathScalar=yes` in
`prog/samples/unittest/vecMathUnitTest`).

## Read by task

| Need | Read |
|---|---|
| Choose or call an API, add a helper; correctness, FP and antipattern rules | `usage.md` |
| Write, port or tune a function; lane utilization and per-instruction cost | `microarch.md` |
| Work in `dag_vecMath_double.h` | `vec4d.md` |
| Inspect emitted code, or benchmark | `compilers.md` |

Read `usage.md` before you add a helper: most of what a caller needs already exists, and a new
function must fit its conventions. Read `microarch.md` while you write the code, not after: its
layout choices are costly to retrofit.

## File map

| File | Purpose |
|---|---|
| `dag_vecMathDecl.h` | types, platform detection, alignment |
| `dag_vecMath.h` | complete platform-independent public API (~650 declarations with comments) |
| `dag_vecMath_const.h` | constants (`V_C_HALF`, `V_C_ONE`, `V_C_PI`, `V_C_UNIT_1000`, `V_CI_MASK*`, etc.) |
| `dag_vecMath_pc_sse.h` / `dag_vecMath_neon.h` / `dag_vecMath_scalar.h` | basic SSE / NEON / scalar primitives |
| `dag_vecMath_common.h` | shared bbox, frustum, quaternion and matrix code on the primitives, without hardware intrinsics |
| `dag_vecMath_trig.h` | polynomial trig approximations |
| `dag_vecMath_double.h` | `vec4d` for SSE/AVX, NEON and scalar; include it through `dag_vecMath.h` |
| `dag_vecMath_fc16.h` | x86 F16C half/float wrappers that need `<immintrin.h>`; a plain SSE4.1 build includes them only on demand |

Declare each new public function in `dag_vecMath.h` with a comment: callers grep it as the API
reference.

## Non-negotiable contracts

- **Matrix aliasing:** callers rely on a destination that is also an input
  (`v_mat44_mul(m, m, rel)`), so a new matrix function must be safe for it. Read every input
  before the first store; `v_mat44_transpose` takes its source by value, and `v_mat44_inverse43`
  copies its input first.
- **Names and precision:** the unsuffixed function is precise and identical on every platform;
  `_est`/`_unprecise` explicitly trade precision for speed, and `_safe` explicitly tolerates edge
  data (see `usage.md`). Never make an estimate the default.
- **Mirrors:** mirrored (det < 0) input has one engine-wide convention: decompose and
  quat-from-matrix code flips col2 and folds the negative into `scale.z`. The scale-removal
  primitives stay sign-neutral: do not fold a det test into them.
- **Cross-platform:** a function added to one ISA path (SSE, NEON, scalar) needs the others, and
  all paths must agree on results, not only compile.

## Implementation rules

- Most functions are `VECMATH_FINLINE` (force-inlined) and rely on dead-code elimination of
  unused results; use plain `VECMATH_INLINE` for a new large body, as ray-box and triangle tests
  do.
- On SSE, `VECTORCALL` passes `vec4f` in registers. A non-inlined call spills different vector
  state per platform: SysV (Linux, PlayStation) treats all xmm as scratch and passes up to 8
  vectors in xmm0-7; Windows x64 keeps xmm6-15 callee-saved and needs `__vectorcall` to pass
  vectors in registers; AArch64 passes in v0-v7 and keeps only the low 64 bits of v8-v15.
- Prefer branchless code (`v_sel`, `v_sel_b`, `v_btsel`, `v_bbox3_sel`). A scalar extract for a
  readable comparison is fine, but never call stdlib float functions in the SIMD backends; only
  the scalar backend maps floor/ceil/trunc/round/sqrt to `math.h`, by design.

## Compile-time budget

Header parse dominates: with nothing used, the backend is about 2% of the `dag_vecMath.h` compile
cost (clang, SSE4.1 PC flags). A large unused inline body is almost free; do not move it to a
`.cpp` for build speed. The `VECMATH_INLINE` bodies are about 6% of the header text, and removing
them saves only 1-2 ms at clang `-O2`. Emitting one costs about 1.4 ms, only in a TU that calls
it. The ones with many callers (`v_is_visible_b`, `v_screen_size_b`, `v_ray_box_intersect_dist`)
are per-object cull code in inner loops, where an out-of-line call would also spill registers.

The intrinsic headers cost more: clang parses `<immintrin.h>` in about 90 ms against 16 ms for
`<smmintrin.h>`, MSVC in about 30 against 3 ms. `dag_vecMath_pc_sse.h` includes the umbrella only
when the target compiles an intrinsic above SSE4.1, and F16C stays in `dag_vecMath_fc16.h`: keep
a new intrinsic behind the same test. AVX targets get the umbrella anyway from
`dag_vecMathDecl.h` (`vec4d` is `__m256d`), and clang rejects a direct include of the narrower
AVX/FMA/F16C headers. Do not use an umbrella-only intrinsic when SSE2/SSE4.1 gives the same code:
`v_set_all_bitsi` compares zero with itself, not `_mm_undefined_si128`, because both compilers
fold every spelling to one form (clang `pcmpeqd`, MSVC a constant load).

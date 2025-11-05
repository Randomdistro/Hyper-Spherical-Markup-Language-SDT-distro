#pragma once

#include <cmath>
#include <algorithm>
#include <numbers>
#include <optional>
#include <utility>
#include <vector>

namespace hsml {
namespace core {
namespace hps21 {

using Scalar = double;

// Constants
inline constexpr Scalar TWO_PI = std::numbers::pi_v<Scalar> * 2.0;
inline constexpr Scalar HALF_PI = std::numbers::pi_v<Scalar> * 0.5;
inline constexpr Scalar EPS = 1e-12;

// Helpers
Scalar wrap2pi(Scalar x);
Scalar clamp(Scalar v, Scalar lo, Scalar hi);
Scalar safeDiv(Scalar n, Scalar d, Scalar fallback = 0);

// 21D state
struct HPS21 {
  Scalar r{1};
  Scalar th{HALF_PI};
  Scalar ph{0};
  Scalar ga{0};
  Scalar sig{0};
  Scalar kr{0};
  Scalar kth{0};
  Scalar kph{0};
  Scalar w{0};
  Scalar phase{0};
  Scalar eps{1};
  Scalar rho_e{1};
  Scalar tau{1};
  Scalar coh{1};
  Scalar occ{0};
  Scalar lvl{0};
  Scalar sym{0};
  Scalar chi{0};
  Scalar bd{1};
  Scalar h{0};
  Scalar pi{0};
};

struct Invariants {
  std::optional<Scalar> totalThroughput{};
  std::optional<Scalar> totalOcclusion{};
  std::optional<Scalar> totalElectron{};
};

HPS21 normalize(const HPS21& s);

// S^2 geometry
Scalar gcDistance(Scalar th1, Scalar ph1, Scalar th2, Scalar ph2);
std::pair<Scalar, Scalar> gcStep(Scalar th, Scalar ph, Scalar psi, Scalar s);

// Operators
HPS21 pressurePush(const HPS21& s, Scalar biasTh, Scalar biasPh, Scalar gain = 1);
HPS21 spinAdvance(const HPS21& s, Scalar dt, Scalar k = 1);
HPS21 phaseExchange(const HPS21& s, Scalar dt, Scalar k = 1);

struct Neighbor { Scalar th; Scalar ph; Scalar solidAngle; };
HPS21 accumulateOcclusion(const HPS21& s, const std::vector<Neighbor>& nbrs);

struct BoundaryMixOut { Scalar transmitted; Scalar reflected; HPS21 s; };
BoundaryMixOut boundaryMix(const HPS21& s, Scalar incident);

// Integration
enum class OccMode { Accumulate, Zero };
struct StepOpts {
  Scalar dt{0.01};
  std::optional<Invariants> invariants{};
  Scalar gpush{1}, gspin{1}, gphase{1};
  OccMode occMode{OccMode::Accumulate};
};

HPS21 enforceInvariants(const HPS21& s, const Invariants& inv);
HPS21 step(const HPS21& s, const StepOpts& opt, Scalar biasTh, Scalar biasPh, const std::vector<Neighbor>& nbrs);

// Utils
HPS21 makeState(const HPS21& init = {});
HPS21 blend(const HPS21& a, const HPS21& b, Scalar t);

// Time-on-Target
struct TargetSpec { Scalar th; Scalar ph; Scalar alpha; Scalar phase_t{0}; };
struct ToTState { Scalar tot{0}; Scalar duty{0}; bool onTarget{false}; };
Scalar dirCos(Scalar th1, Scalar ph1, Scalar th2, Scalar ph2);
ToTState updateToT(const HPS21& s, const TargetSpec& t, const ToTState& in, Scalar dt, Scalar T_c = 1.0);

} // namespace hps21
} // namespace core
} // namespace hsml



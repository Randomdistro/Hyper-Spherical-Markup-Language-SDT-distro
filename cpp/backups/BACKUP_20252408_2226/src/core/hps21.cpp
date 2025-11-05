#include "hsml/core/hps21.h"
#include <cmath>
#include <numbers>

namespace hsml {
namespace core {
namespace hps21 {

using Scalar = double;

Scalar wrap2pi(Scalar x) {
    x = std::fmod(x, std::numbers::pi_v<Scalar> * 2.0);
    return (x < 0) ? x + std::numbers::pi_v<Scalar> * 2.0 : x;
}

Scalar clamp(Scalar v, Scalar lo, Scalar hi) {
    return v < lo ? lo : (v > hi ? hi : v);
}

Scalar safeDiv(Scalar n, Scalar d, Scalar fallback) {
    if (!std::isfinite(n) || !std::isfinite(d) || std::abs(d) < 1e-12) return fallback;
    return n / d;
}

HPS21 normalize(const HPS21& s) {
    HPS21 o = s;
    o.th = clamp(o.th, 0.0, std::numbers::pi_v<Scalar>);
    o.ph = wrap2pi(o.ph);
    o.ga = wrap2pi(o.ga);
    o.coh = clamp(o.coh, 0.0, 1.0);
    o.occ = clamp(o.occ, 0.0, 1.0);
    o.bd  = clamp(o.bd,  0.0, 1.0);
    return o;
}

Scalar gcDistance(Scalar th1, Scalar ph1, Scalar th2, Scalar ph2) {
    const Scalar dth = th2 - th1;
    const Scalar dph = ph2 - ph1;
    const Scalar a = std::pow(std::sin(dth * 0.5), 2) + std::sin(th1) * std::sin(th2) * std::pow(std::sin(dph * 0.5), 2);
    const Scalar x = std::fmin<Scalar>(1.0, std::sqrt(a));
    return 2.0 * std::asin(x);
}

std::pair<Scalar, Scalar> gcStep(Scalar th, Scalar ph, Scalar psi, Scalar s) {
    const Scalar sinth = std::sin(th);
    const Scalar kx = -std::sin(psi);
    const Scalar ky =  std::cos(psi);
    const Scalar dth = kx;
    const Scalar dph = safeDiv(ky, (std::abs(sinth) < 1e-12 ? std::copysign(1e-12, sinth) : sinth), 0.0);
    Scalar th2 = clamp(th + dth * s, 0.0, std::numbers::pi_v<Scalar>);
    Scalar ph2 = wrap2pi(ph + dph * s);
    return { th2, ph2 };
}

HPS21 pressurePush(const HPS21& s, Scalar biasTh, Scalar biasPh, Scalar gain) {
    const Scalar d = gcDistance(s.th, s.ph, biasTh, biasPh);
    const Scalar step = gain * s.occ * s.tau * safeDiv(1.0, d, 0.0);
    auto dir = gcStep(s.th, s.ph, s.ga, step);
    HPS21 o = s; o.th = dir.first; o.ph = dir.second; return normalize(o);
}

HPS21 spinAdvance(const HPS21& s, Scalar dt, Scalar k) {
    const Scalar shear = (std::abs(s.kth) + std::abs(s.kph)) * std::exp(-std::abs(s.sig));
    HPS21 o = s;
    o.ga = wrap2pi(s.ga + k * shear * dt);
    o.w  = s.w + k * shear * dt * (s.chi >= 0 ? 1.0 : -1.0);
    return normalize(o);
}

HPS21 phaseExchange(const HPS21& s, Scalar dt, Scalar k) {
    const Scalar dcoh = k * (s.eps > 0 ? 1.0 : -1.0) * std::sin(s.phase) * dt;
    HPS21 o = s;
    o.coh = clamp(s.coh + dcoh, 0.0, 1.0);
    o.eps = s.eps - dcoh;
    o.phase = wrap2pi(s.phase + dt * (1.0 + s.tau));
    return normalize(o);
}

HPS21 accumulateOcclusion(const HPS21& s, const std::vector<Neighbor>& neighbors) {
    Scalar sum = 0.0;
    for (const auto& n : neighbors) sum += n.solidAngle;
    HPS21 o = s; o.occ = clamp(s.occ + sum / (4.0 * std::numbers::pi_v<Scalar>), 0.0, 1.0);
    return normalize(o);
}

BoundaryMixOut boundaryMix(const HPS21& s, Scalar incident) {
    const Scalar t = clamp(s.bd, 0.0, 1.0);
    const Scalar transmitted = incident * t;
    const Scalar reflected   = incident * (1.0 - t);
    HPS21 o = s; o.eps = s.eps + transmitted;
    return { transmitted, reflected, normalize(o) };
}

HPS21 enforceInvariants(const HPS21& s, const Invariants& inv) {
    HPS21 o = s;
    if (inv.totalOcclusion && *inv.totalOcclusion >= 0.0) {
        const Scalar target = clamp(*inv.totalOcclusion, 0.0, 1.0);
        o.occ = clamp(0.5 * (o.occ + target), 0.0, 1.0);
    }
    if (inv.totalThroughput && *inv.totalThroughput >= 0.0) {
        const Scalar t = std::max<Scalar>(1e-9, *inv.totalThroughput);
        o.eps = std::min(o.eps, t);
    }
    if (inv.totalElectron && *inv.totalElectron >= 0.0) {
        o.rho_e = std::max<Scalar>(0.0, o.rho_e);
    }
    return normalize(o);
}

HPS21 step(const HPS21& s, const StepOpts& opt, Scalar biasTh, Scalar biasPh, const std::vector<Neighbor>& nbrs) {
    HPS21 out = s;
    const bool occZero = (opt.occMode == OccMode::Zero);
    if (occZero) { out.occ = 0; }
    out = pressurePush(out, biasTh, biasPh, opt.gpush * opt.dt);
    out = spinAdvance(out, opt.dt, opt.gspin);
    out = phaseExchange(out, opt.dt, opt.gphase);
    if (!occZero) out = accumulateOcclusion(out, nbrs); else out.occ = 0;
    if (opt.invariants) out = enforceInvariants(out, *opt.invariants);
    return normalize(out);
}

HPS21 makeState(const HPS21& init) { return normalize(init); }

HPS21 blend(const HPS21& a, const HPS21& b, Scalar t) {
    auto lerp = [&](Scalar x, Scalar y){ return x + (y - x) * t; };
    HPS21 o;
    o.r = lerp(a.r, b.r); o.th = lerp(a.th, b.th); o.ph = wrap2pi(lerp(a.ph, b.ph)); o.ga = wrap2pi(lerp(a.ga, b.ga)); o.sig = lerp(a.sig, b.sig);
    o.kr = lerp(a.kr, b.kr); o.kth = lerp(a.kth, b.kth); o.kph = lerp(a.kph, b.kph);
    o.w = lerp(a.w, b.w); o.phase = wrap2pi(lerp(a.phase, b.phase));
    o.eps = lerp(a.eps, b.eps); o.rho_e = lerp(a.rho_e, b.rho_e); o.tau = lerp(a.tau, b.tau);
    o.coh = clamp(lerp(a.coh, b.coh), 0.0, 1.0); o.occ = clamp(lerp(a.occ, b.occ), 0.0, 1.0);
    o.lvl = lerp(a.lvl, b.lvl); o.sym = lerp(a.sym, b.sym); o.chi = lerp(a.chi, b.chi);
    o.bd = clamp(lerp(a.bd, b.bd), 0.0, 1.0); o.h = lerp(a.h, b.h); o.pi = lerp(a.pi, b.pi);
    return normalize(o);
}

Scalar dirCos(Scalar th1, Scalar ph1, Scalar th2, Scalar ph2) {
    const Scalar s1 = std::sin(th1), c1 = std::cos(th1);
    const Scalar s2 = std::sin(th2), c2 = std::cos(th2);
    const Scalar dph = ph2 - ph1;
    return c1 * c2 + s1 * s2 * std::cos(dph);
}

ToTState updateToT(const HPS21& s, const TargetSpec& t, const ToTState& in, Scalar dt, Scalar T_c) {
    const Scalar cosSep = dirCos(s.th, s.ph, t.th, t.ph);
    const Scalar alpha = std::max<Scalar>(t.alpha, 1e-12);
    const bool on = (cosSep >= std::cos(alpha));
    const Scalar safeOccFactor = (1 - s.occ);
    const Scalar weight = s.eps * s.coh * safeOccFactor * (1 + s.tau);
    const Scalar incr = on ? weight * dt : 0.0;
    const Scalar dutyInst = on ? 1.0 : 0.0;
    const Scalar a = std::exp(-dt / std::max<Scalar>(1e-12, T_c));
    const Scalar duty = a * in.duty + (1.0 - a) * dutyInst;
    return { in.tot + incr, duty, on };
}

} // namespace hps21
} // namespace core
} // namespace hsml



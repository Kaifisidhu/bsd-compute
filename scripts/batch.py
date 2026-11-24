#!/usr/bin/env python3
"""Batch payload for the bsd-compute farm.

For each curve in curves.json (label, certified rank, minimal a-invariants,
conductor — drawn from the BSD campaign's certified set), compute:
  * a_p table for all primes p <= 500 (PARI ellap)
  * analytic rank + leading L-value (ellanalyticrank)
  * analytic Sha via the BSD quotient:
        Sha_an = L^(r)(1)/r!  *  #tors^2 / (Omega * Reg * prod c_p)
    with Reg = ellregulator (1 for rank 0), prod c_p = elltamagawa
    (on these MINIMAL models), Omega = real period.
All PARI work is done by GP string evaluation (campaign rule: cypari2
convenience wrappers silently cap precision; strings do not).

Usage: batch.py curves.json out.json
Prints per-curve timings and overall throughput (curves/hour).
"""
import json
import sys
import time

import cypari2

pari = cypari2.Pari()
pari.allocatemem(1 << 29)


def gp(s):
    return pari(s)


def batch(curves_path, out_path):
    curves = json.load(open(curves_path))
    results = []
    t_all = time.time()
    for c in curves:
        t0 = time.time()
        try:
            row = process_curve(c)
        except Exception as e:  # one bad curve must not kill the batch
            row = {"label": c["label"], "error": f"{type(e).__name__}: {e}"}
            print(f"{c['label']:>8}  ERROR {row['error']}", flush=True)
        row["seconds"] = round(time.time() - t0, 3)
        results.append(row)
    wall = time.time() - t_all
    summary = {
        "curves": len(results),
        "wall_seconds": round(wall, 2),
        "curves_per_hour": round(3600.0 * len(results) / wall, 1) if wall > 0 else None,
        "results": results,
    }
    with open(out_path, "w") as f:
        json.dump(summary, f, indent=1)
    print(f"TOTAL {len(results)} curves in {wall:.1f}s -> {summary['curves_per_hour']} curves/hour")
    return summary


def process_curve(c):
    ainvs = str(c["ainvs"])
    gp(f"E = ellinit({ainvs})")
    conductor = int(gp("ellglobalred(E)[1]"))
    gp("P = select(p -> p <= 500, primes(100))")
    ap = [int(x) for x in gp("vector(#P, i, ellap(E, P[i]))")]
    primes = [int(x) for x in gp("P")]
    ar = gp("ellanalyticrank(E)")
    rank_an = int(ar[0])
    llead = float(ar[1])
    tors = int(gp("elltors(E)[1]"))
    # Tamagawa product over FINITE primes only, from elllocalred's c_p.
    # Convention pairing (the E_34 lesson): on the campaign's PARI 2.17
    # build, elltamagawa(E) = (# real components) x (finite product) --
    # E_34: 32 vs 16 = c_2 * c_17 = 4 * 4; 37a1: 2 vs c_37 = 1. The BSD
    # quotient is only correct if Omega and the Tamagawa factor use
    # MATCHING conventions: below, Omega is the FULL real period
    # (component period x #components, the campaign's Omega_E = 2*w1
    # for E_34) and tam is the finite product. elltamagawa's raw value
    # is recorded for comparison.
    tam_raw = int(gp("elltamagawa(E)"))
    gp("BAD = ellglobalred(E)[4]")
    tam = int(gp("prod(i = 1, matsize(BAD)[1], elllocalred(E, BAD[i,1])[4])"))
    # Regulator: no ellregulator in this PARI build; det of the height
    # matrix over the ellrank generators instead.
    reg = float(gp("matdet(ellheightmatrix(E, ellrank(E)[4]))")) if rank_an > 0 else 1.0
    # Real period Omega = positive generator of (period lattice cap R).
    # ellperiods' basis varies: sometimes one basis vector is real,
    # sometimes the real period is a small combination (e.g. w1 + w2
    # for a conjugate basis). Search small combinations for the
    # smallest positive real lattice point. (Validated: 11a1 -> 1.)
    gp("PER = ellperiods(E)")
    gp("OMEGA = 0")
    gp("for(m = -3, 3, for(n = -3, 3, if(m || n, "
       "z = m*PER[1] + n*PER[2]; "
       "if(abs(imag(z)) < 1e-10 * max(1, abs(z)) && real(z) > 0, "
       "if(OMEGA == 0 || real(z) < OMEGA, OMEGA = real(z))))))")
    omega = float(gp("OMEGA"))
    # Full BSD real period: component period x number of real
    # components (2 iff disc > 0). Pairs with the finite Tamagawa
    # product above.
    ncomp = 2 if int(gp("sign(E.disc)")) > 0 else 1
    omega_bsd = omega * ncomp
    sha_raw = (llead * tors * tors / (omega_bsd * reg * tam)
               if rank_an == c["certified_rank"] and omega > 0 else None)
    # ellrank's returned generators can span a finite-index subgroup
    # (witnessed: 37a1, where ellrank returns 3x the true generator, so
    # the naive regulator is 9x too large and Sha_an comes out 1/9).
    # Detect an index-m inflation: smallest m with sha_raw * m^2 within
    # 1e-6 of a positive integer. Flagged, never silent; the campaign's
    # CERT pipeline remains the authority on Sha.
    sha_an, index_corr = sha_raw, 1
    if sha_raw is not None:
        for m in range(1, 8):
            v = sha_raw * m * m
            if abs(v - round(v)) < 1e-6 and round(v) >= 1:
                sha_an, index_corr = v, m
                break
    row = {
        "label": c["label"],
        "conductor": conductor,
        "conductor_ok": conductor == c["conductor"],
        "certified_rank": c["certified_rank"],
        "analytic_rank": rank_an,
        "ap_primes": primes,
        "ap": ap,
        "L_leading": llead,
        "torsion_order": tors,
        "tamagawa_product_finite": tam,
        "elltamagawa_raw": tam_raw,
        "regulator": reg,
        "real_period_component": omega,
        "real_period_bsd_full": omega_bsd,
        "sha_analytic_raw": sha_raw,
        "ellrank_index_correction": index_corr,
        "sha_analytic": sha_an,
    }
    print(f"{c['label']:>8}  N={conductor:<7} rank_an={rank_an} "
          f"Sha_an={sha_an if sha_an is None else round(sha_an, 6)}"
          f"{'' if index_corr == 1 else f' (raw {round(sha_raw, 6)}, index x{index_corr})'}",
          flush=True)
    return row


if __name__ == "__main__":
    batch(sys.argv[1], sys.argv[2])

\\ PARI-VERIFY v3 (FARM-FIXED): E34, y^2 = x^3 - 1156x, N = 18496, p = 5.
\\ Independent route: PARI overconvergent modular symbols (Pollack-Stevens).
\\ NO Sage/eclib code anywhere in this path. Eigensymbol by Hecke eigenvalues
\\ (msinit -> msnew -> msfromhecke); normalisation by the twist-sum identity
\\ as EVALUATED by mseval (verified against msfromell on 11a1 and 389a1):
\\ sum_a (D|a) x^+({a/D, oo}) = -sqrt(D) * L(E^D,1)/Omega^+  for fundamental
\\ D > 0 — the sqrt(D) is the Gauss sum tau((D|.)); the msfromell doc's
\\ displayed identity omits it, and v2's uncorrected use made its lambda
\\ irrational and D-dependent (its agreement gate could never have passed).
\\ lambda := S_D(phi_raw)/(-sqrt(D)*L(E^D,1)/Omega^+) is D-independent
\\ rational, recognised exactly by lindep. L-values via PARI lfun.
\\ RUN ONLY AFTER pari_control.gp PASSES (c2 = 221 mod 5^5, 18971 mod 5^8).
\\ CLAIM under test: c2 = 43338 mod 5^7, v5(c2) = 0; c0, c1 vanish.
\\ nprec = 10 -> c2 known mod 5^8 (one digit beyond the claim).
\\
\\ ================= REPAIR NOTES (v3) — full diagnosis in FIX-NOTES.md =====
\\ Same three repairs as pari_control.gp v3:
\\ (1) BASIS: msfromhecke(M, v, NW) returns NEW-subspace coordinates
\\     (observed PARI 2.15 farm and PARI 2.17 local). v3 detects the basis
\\     by row count, embeds via NW[1], never multiplies mismatched dims.
\\ (2) GRAMMAR: no multi-statement for()/if() bodies, no braced blocks
\\     nested in loop bodies (both mis-parse when read from a file);
\\     all logic in functions, loop bodies are single expressions.
\\ (3) FAIL-FAST: everything inside main(); any error or failed chk()
\\     aborts; "E34 DONE" prints only as main's final action; NOTHING
\\     follows the main() call. Route-validity checks (conductor, a5,
\\     symbol, eigenvector, lambda agreement, twist-2 exactness) abort.
\\     The CLAIM valuations (v5(c2 - 43338), v5(q5 - 1)) are RESULTS, not
\\     route checks: they are printed, never asserted — disagreement with
\\     the claim is data for the reader, and E34 DONE still prints.
\\     Stage markers: "STAGE eigensymbol OK", "STAGE normalisation OK",
\\     "STAGE oms OK". Gate guidance for the workflow is in FIX-NOTES.md.
default(timer, 1);
default(realprecision, 60);

\\ ---------- assertion helper: a failed check aborts the whole run ----------
chk(c, msg) = if(!c, error(msg));

\\ ---------- globals shared between stages (see pari_control.gp v3) ---------
trytwist(D) =
{
  my(Et, Lv, Ssum, ratio, w, lam);
  Et = ellinit([0,0,0, GA*D^2, GB*D^3]);
  Lv = lfun(Et, 1);
  print("twist D=", D, ": L(E^D,1) = ", Lv);
  if(abs(Lv) <= 10^(-20), return(0));
  Ssum = sum(a=0, D-1, kronecker(D,a) * mseval(G_M, G_phiraw, [a/D, oo]));
  ratio = Ssum / (-sqrt(D) * (Lv / G_Om));
  w = lindep([ratio, 1]);
  lam = -w[2]/w[1];
  print("  S_D(raw) = ", Ssum, " ; recognised lambda = ", lam);
  G_lams = concat(G_lams, [lam]);
  G_used = concat(G_used, [D]);
  return(1);
};

collecttwist(D) = if(#G_lams < 2, trytwist(D), 0);

evecstep(p) =
{
  if(G_N % p == 0, return(0));
  G_checked = G_checked + 1;
  if(mshecke(G_M, p) * G_phiraw != ellap(G_E, p) * G_phiraw, G_bad = G_bad + 1);
  return(0);
};

\\ ---------- main ------------------------------------------------------------
main() =
{
  my(E, N, M, NW, plist, v, S, dFull, dNew, r, phi_raw);
  my(Om, Ds, lambda, phi, D2, Et2, Lr2, w2, q2, S2);
  my(nprec, Mp, PHI, mu, ser, c0, c1, c2, c3, c4, claim, cnt);
  my(detH, rts, alpha, lg, Reg_gamma, q5);
  GA = -1156; GB = 0;
  E = ellinit([0,0,0, GA, GB]);
  G_E = E;
  N = ellglobalred(E)[1];
  G_N = N;
  print("E34 conductor = ", N, " (expect 18496)");
  chk(N == 18496, "unexpected conductor for E34");
  print("a5 = ", ellap(E,5), " (expect -2)");
  chk(ellap(E,5) == -2, "unexpected a5 for E34");
  print("Omega^+ = E.omega[1] = ", E.omega[1], " (frozen campaign value 0.8993583214466512...)");
  M = msinit(N, 2, 1);
  G_M = M;
  print("msinit done: dim full + space = ", msdim(M));
  NW = msnew(M);
  print("msnew done: dim new + subspace = ", matsize(NW[1])[2]);
  plist = [3,5,7,11,13,29,37,41];
  v = vector(#plist, i, [plist[i], ellap(E, plist[i])]);
  print("eigenvalue constraints (p, a_p): ", v);
  S = msfromhecke(M, v, NW);
  print("msfromhecke returned matrix size: ", matsize(S));
  dFull = msdim(M);
  dNew = matsize(NW[1])[2];
  chk(dFull != dNew, "basis ambiguity: full dim equals new-subspace dim");
  chk(matsize(S)[2] >= 1, "msfromhecke returned no columns (inconsistent eigenvalues?)");
  r = matsize(S)[1];
  if(r == dNew,
    phi_raw = NW[1] * S[,1],
    if(r == dFull,
      phi_raw = S[,1],
      error("msfromhecke returned unexpected row count")));
  chk(#phi_raw == dFull, "eigensymbol length mismatch after embedding");
  chk(msissymbol(M, phi_raw) == 1, "embedded eigensymbol is not a modular symbol");
  chk(matsize(S)[2] == 1, "eigenspace dimension is not 1");
  G_phiraw = phi_raw;
  print("eigensymbol embedded in full basis; length = ", #phi_raw);
  cnt = 0;
  forprime(p=2, 50, if(N % p != 0, cnt = cnt + 1));
  G_bad = 0; G_checked = 0;
  forprime(p=2, 50, evecstep(p));
  chk(G_checked == cnt, "eigenvector verification did not evaluate every prime");
  print("eigenvector failures at primes <= 50 not dividing N: ", G_bad);
  print("eigenvector primes genuinely checked: ", G_checked);
  chk(G_bad == 0, "eigensymbol fails a Hecke eigenvalue check");
  print("STAGE eigensymbol OK");
  \\ ---- Normalisation via two independent quadratic twists ----
  Om = E.omega[1];
  G_Om = Om;
  Ds = [5, 13, 21, 29, 37];
  G_lams = []; G_used = [];
  for(i=1, #Ds, collecttwist(Ds[i]));
  chk(#G_lams == 2, "fewer than two twists with L(E^D,1) != 0 found");
  print("twists used: ", G_used, " ; lambdas: ", G_lams);
  print("lambda agreement (exact rational equality): ", G_lams[1] == G_lams[2]);
  chk(G_lams[1] == G_lams[2], "lambda disagreement between the two twists");
  lambda = G_lams[1];
  phi = phi_raw / lambda;
  D2 = G_used[2];
  Et2 = ellinit([0,0,0, GA*D2^2, GB*D2^3]);
  Lr2 = -sqrt(D2) * (lfun(Et2, 1) / Om);
  w2 = lindep([Lr2, 1]);
  q2 = -w2[2]/w2[1];
  S2 = sum(a=0, D2-1, kronecker(D2,a) * mseval(M, phi, [a/D2, oo]));
  print("twist-2 check: normalised symbol sum = ", S2, " ; recognised -sqrt(D)*L-ratio = ", q2, " ; EXACTLY equal: ", S2 == q2);
  chk(S2 == q2, "twist-2 exact normalisation check failed");
  print("STAGE normalisation OK");
  \\ ---- Overconvergent lift and p-adic L series ----
  \\ Convention factor: PARI's OMS series from the msfromell x^+ equals
  \\ c_inf times Sage/eclib's L_ratio series, c_inf = # real components
  \\ (2 iff Delta > 0) — verified exactly on 11a1 (factor 1), 389a1 and
  \\ 64a1 (factor 2). The claim 43338 is Sage-convention, so the OMS
  \\ input symbol is phi / c_inf. See FIX-NOTES.md.
  cinf = if(E.disc > 0, 2, 1);
  print("real components c_inf = ", cinf, " (Delta = ", E.disc, ") ; OMS input symbol = phi / c_inf");
  phioms = phi / cinf;
  nprec = 10;
  Mp = mspadicinit(M, 5, nprec, 0);
  print("mspadicinit done, n = ", nprec, " (flag 0 = ordinary symbols only; a5=-2 is a 5-adic unit)");
  PHI = mstooms(Mp, phioms);
  print("mstooms done");
  mu = mspadicmoments(Mp, PHI);
  print("mspadicmoments done");
  ser = mspadicseries(mu);
  print("SERIES = ", ser);
  chk(type(ser) == "t_POL" || type(ser) == "t_SER", "mspadicseries returned unexpected type");
  c0 = polcoef(ser, 0); c1 = polcoef(ser, 1); c2 = polcoef(ser, 2);
  c3 = polcoef(ser, 3); c4 = polcoef(ser, 4);
  chk(type(c2) == "t_PADIC", "c2 is not a 5-adic number");
  print("c0 = ", c0);
  print("c1 = ", c1);
  print("c2 = ", c2);
  print("c3 = ", c3);
  print("c4 = ", c4);
  print("STAGE oms OK");
  print("v5(c2) = ", valuation(c2));
  claim = 3 + 2*5 + 3*5^2 + 5^3 + 4*5^4 + 3*5^5 + 2*5^6 + O(5^30);
  print("v5(c2 - claim43338) = ", valuation(c2 - claim), " (>=7 iff agreement with claim mod 5^7)");
  print("DIAGNOSTIC (not a gate): v5(c2 + claim43338) = ", valuation(c2 + claim), " (large ONLY if a sign/orientation flip occurred)");
  print("v5(c0) = ", valuation(c0), " ; v5(c1) = ", valuation(c1));
  \\ ---- Item 4: q5 from THIS route's c2 via eq. (3.1) ----
  \\ det H digits: kappa5-H compute_H.out (ellpadicregulator, O(5^82)); terms to 5^14.
  detH = 3*5^2 + 4*5^3 + 2*5^4 + 5^6 + 3*5^8 + 4*5^9 + 4*5^11 + 5^12 + 2*5^15 + O(5^16);
  \\ alpha: unit root of X^2 + 2X + 5 in Z_5 (alpha = 3 mod 5), own lift:
  rts = polrootspadic(x^2 + 2*x + 5, 5, 16);
  alpha = 0;
  for(i=1, #rts, if(Mod(lift(rts[i]), 5) == 3, alpha = rts[i]));
  chk(alpha != 0, "unit root alpha of X^2 + 2X + 5 not found");
  print("alpha (own polrootspadic lift) = ", alpha);
  lg = log(6 + O(5^16));
  Reg_gamma = detH / lg^2;
  q5 = c2 / ((1 - alpha^(-1))^2 * Reg_gamma);
  print("Reg_gamma = ", Reg_gamma);
  print("q5 = ", q5);
  print("v5(q5 - 1) = ", valuation(q5 - 1), " (>=7 iff q5 = 1 mod 5^7)");
  print("E34 DONE");
};

main();

\\ PARI-VERIFY CONTROL v3 (FARM-FIXED): E56, y^2 = x^3 - 56x, N = 12544, p = 5.
\\ Independent route: PARI overconvergent modular symbols (Pollack-Stevens).
\\ NO Sage/eclib code anywhere in this path. Eigensymbol built by Hecke
\\ eigenvalues (msinit -> msnew -> msfromhecke), NOT msfromell (which
\\ overflows the PARI stack at this level in this build: see log v1).
\\ Normalisation: fixed by the twist-sum identity for the msfromell symbol
\\ x^+, as EVALUATED by mseval in this build (verified empirically against
\\ msfromell itself on 11a1 and on 389a1, several fundamental D > 0):
\\   sum_a (D|a) x^+({a/D, oo}) = -sqrt(D) * L(E^D,1)/Omega^+,
\\ the sqrt(D) being the Gauss sum tau((D|.)) and the sign mseval's path
\\ orientation. (The msfromell doc's displayed identity omits the Gauss
\\ factor; v2 used it uncorrected — a third, mathematical defect: its
\\ "lambda" was -lambda*sqrt(D), irrational and D-dependent, so v2's
\\ lambda-agreement gate could NEVER have passed. See FIX-NOTES.md.)
\\ With L(E^D,1) from PARI lfun (analytic route, independent of symbols),
\\ lambda := S_D(phi_raw) / (-sqrt(D) * L(E^D,1)/Omega^+) is therefore a
\\ D-independent rational, recognised exactly by lindep.
\\ Expected (archive/Sage): c2 = 221 mod 5^5; fresh n=9 value 18971 mod 5^8.
\\
\\ ================= REPAIR NOTES (v3) — full diagnosis in FIX-NOTES.md =====
\\ (1) BASIS: msfromhecke(M, v, NW) returns coordinates in the NEW-subspace
\\     basis (row count = #columns of NW[1]) — observed on the farm's
\\     PARI 2.15 ([318,1] at msdim 1796 / newdim 318, run 45333918122) and
\\     locally on PARI 2.17 ([32,1] at msdim 33 / newdim 32, level 389).
\\     v2 instead "tested" the convention by multiplying the 318-vector by
\\     the 1796x1796 Hecke matrix; the GP error aborted its test loops
\\     BEFORE any comparison evaluated, leaving the pass-flags at their
\\     initial values, and v2 then embedded nothing and cascaded. v3
\\     detects the basis by ROW COUNT only, embeds via NW[1] when in
\\     new-coordinates, and never multiplies mismatched dimensions.
\\ (2) GRAMMAR: multi-statement bodies inside for()/if() arguments are not
\\     valid GP, and braced blocks nested inside loop bodies mis-parse when
\\     read from a file (gp binary: "unexpected closing brace"; cypari2
\\     reader: "embedded braces not yet implemented"). v2's twist loop was
\\     therefore unparseable on every build — a second, independent defect.
\\     v3 puts all logic in functions; loop bodies are single expressions
\\     (function calls); if() branches are single expressions (return()).
\\ (3) FAIL-FAST: everything runs inside main(). Any GP error or failed
\\     chk() aborts main; "CONTROL DONE" prints only as main's final
\\     action, and NOTHING follows the main() call in this file. Under the
\\     gp binary an aborted run still exits 0, so the workflow gate must
\\     grep for the DONE marker + digit lines AND for absence of "***".
\\     Stage markers: "STAGE eigensymbol OK", "STAGE normalisation OK",
\\     "STAGE oms OK". The eigenvector verification counts a prime as
\\     checked only after its comparison has been evaluated, and asserts
\\     the checked count equals the number of primes <= 50 not dividing N:
\\     a false "0 failures" line is impossible in v3.
\\ (4) VERSIONS: the code relies only on behaviour verified identical on
\\     PARI 2.15 (farm, apt pari-gp) and 2.17 (local): msinit/msnew/
\\     msfromhecke signatures, msfromhecke's new-basis return (see (1)),
\\     mseval, mspadicinit/mstooms/mspadicmoments/mspadicseries.
default(timer, 1);
default(realprecision, 60);

\\ ---------- assertion helper: a failed check aborts the whole run ----------
chk(c, msg) = if(!c, error(msg));

\\ ---------- globals shared between stages ----------
\\ GA, GB : short-model coefficients; the curve is y^2 = x^3 + GA*x + GB and
\\          its quadratic twist by D is y^2 = x^3 + GA*D^2*x + GB*D^3.
\\ G_E, G_N, G_M, G_phiraw (eigensymbol embedded in the FULL basis), G_Om,
\\ G_lams, G_used (twist results), G_bad, G_checked (eigenvector counters).

\\ ---------- one twist step: on success appends lambda and D to globals ----
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

\\ ---------- one prime step of the all-primes eigenvector verification -----
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
  my(nprec, Mp, PHI, mu, ser, c0, c1, c2, c3, c4, target5, target8, cnt);
  GA = -56; GB = 0;
  E = ellinit([0,0,0, GA, GB]);
  G_E = E;
  N = ellglobalred(E)[1];
  G_N = N;
  print("E56 conductor = ", N, " (expect 12544)");
  chk(N == 12544, "unexpected conductor for E56");
  print("a5 = ", ellap(E,5), " (expect -2)");
  chk(ellap(E,5) == -2, "unexpected a5 for E56");
  print("Omega^+ = E.omega[1] = ", E.omega[1]);
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
  \\ Basis detection by row count ONLY (never trial-multiply mismatched dims):
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
  \\ Genuine all-primes eigenvector verification: G_checked increments only
  \\ after the comparison for that prime has actually been evaluated (any
  \\ error aborts main), and the total is asserted against an independent
  \\ count of the primes <= 50 not dividing N.
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
  \\ Exact-rational verification: with phi normalised, the second twist's
  \\ symbol sum must equal -sqrt(D2)*L(E^D2,1)/Omega^+ exactly; the latter
  \\ is rational in truth (it is S_{D2}(x^+)) and is recognised by lindep.
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
  \\ Convention factor (FOURTH v2 defect, found in repair validation):
  \\ PARI's OMS series built from the msfromell x^+ equals c_inf times
  \\ Sage/eclib's L_ratio series, where c_inf = # real components of E
  \\ (2 iff Delta > 0, else 1) — verified exactly on 11a1 (Delta<0,
  \\ factor 1), 389a1 and 64a1 (Delta>0, factor 2; rank- and Manin-
  \\ constant-independent). The archive digits are Sage-convention, so
  \\ the OMS input symbol is phi / c_inf. (v5(c2) is unaffected — 2 is
  \\ a 5-adic unit — but digit-exact agreement needs the factor.)
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
  target5 = 1 + 4*5 + 3*5^2 + 5^3 + O(5^30);
  print("v5(c2 - (1+4*5+3*5^2+5^3)) = ", valuation(c2 - target5), " (>=5 iff agreement mod 5^5)");
  print("DIAGNOSTIC (not a gate): v5(c2 + (1+4*5+3*5^2+5^3)) = ", valuation(c2 + target5), " (large ONLY if a sign/orientation flip occurred)");
  chk(valuation(c2 - target5) >= 5, "CONTROL GATE FAILED: c2 disagrees with archive mod 5^5");
  target8 = 18971 + O(5^30);
  print("v5(c2 - 18971) = ", valuation(c2 - target8), " (>=8 iff agreement with fresh n=9 record mod 5^8)");
  chk(valuation(c2 - target8) >= 8, "CONTROL GATE FAILED: c2 disagrees with n=9 record mod 5^8");
  print("v5(c0) = ", valuation(c0), " ; v5(c1) = ", valuation(c1));
  print("CONTROL DONE");
};

main();

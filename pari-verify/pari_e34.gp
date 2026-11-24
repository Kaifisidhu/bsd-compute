\\ PARI-VERIFY: E34, y^2 = x^3 - 1156x, N = 18496, p = 5. (v2 pipeline)
\\ Independent route: PARI overconvergent modular symbols (Pollack-Stevens).
\\ NO Sage/eclib code anywhere in this path. Eigensymbol by Hecke eigenvalues
\\ (msinit -> msnew -> msfromhecke); normalisation by the twist-sum identity
\\ sum_a (D|a) x^+([a/D]-[oo]) = L(E^D,1)/Omega^+ (msfromell doc), L via lfun.
\\ RUN ONLY AFTER pari_control.gp PASSES (c2 = 221 mod 5^5, 18971 mod 5^8).
\\ CLAIM under test: c2 = 43338 mod 5^7, v5(c2) = 0; c0, c1 vanish.
\\ nprec = 10 -> c2 known mod 5^8 (one digit beyond the claim).
default(timer, 1);
default(realprecision, 60);
E = ellinit([0,0,0,-1156,0]);
N = ellglobalred(E)[1];
print("E34 conductor = ", N, " (expect 18496)");
print("a5 = ", ellap(E,5), " (expect -2)");
print("Omega^+ = E.omega[1] = ", E.omega[1], " (frozen campaign value 0.8993583214466512...)");
M = msinit(N, 2, 1);
print("msinit done: dim full + space = ", msdim(M));
NW = msnew(M);
print("msnew done: dim new + subspace = ", matsize(NW[1])[2]);
plist = [3,5,7,11,13,29,37,41];
v = vector(#plist, i, [plist[i], ellap(E, plist[i])]);
print("eigenvalue constraints (p, a_p): ", v);
S = msfromhecke(M, v, NW);
print("msfromhecke returned matrix size: ", matsize(S));
testvec = S[,1];
okM = msissymbol(M, testvec);
eigM = 1;
for(i=1, #plist, if(mshecke(M,plist[i])*testvec != ellap(E,plist[i])*testvec, eigM=0));
print("columns-as-M-coords: issymbol = ", okM, ", eigenvector for all p in plist: ", eigM);
testvec2 = NW[1] * S[,1];
eigM2 = 1;
for(i=1, #plist, if(mshecke(M,plist[i])*testvec2 != ellap(E,plist[i])*testvec2, eigM2=0));
print("columns-as-NW-coords (x NW[1]): eigenvector for all p in plist: ", eigM2);
if(eigM, phi_raw = testvec, phi_raw = testvec2);
print("phi_raw chosen; eigenspace dimension = ", matsize(S)[2]);
bad = 0;
forprime(p=2, 50, if(N % p != 0 && mshecke(M,p)*phi_raw != ellap(E,p)*phi_raw, bad++));
print("eigenvector failures at primes <= 50 not dividing N: ", bad);
\\ ---- Normalisation via two independent quadratic twists ----
Om = E.omega[1];
Ds = [5, 13, 21, 29, 37];
lams = [];
used = [];
for(i=1, #Ds,
  D = Ds[i];
  Et = ellinit([0,0,0,-1156*D^2,0]);
  Lv = lfun(Et, 1);
  print("twist D=", D, ": L(E^D,1) = ", Lv);
  if(abs(Lv) > 10^(-20),
    Ssum = sum(a=0, D-1, kronecker(D,a) * mseval(M, phi_raw, [a/D, oo]));
    ratio = Ssum / (Lv / Om);
    w = lindep([ratio, 1]);
    lam = -w[2]/w[1];
    print("  S_D(raw) = ", Ssum, " ; recognised lambda = ", lam);
    lams = concat(lams, [lam]); used = concat(used, [D]);
    if(#lams >= 2, break)
  )
);
print("twists used: ", used, " ; lambdas: ", lams);
print("lambda agreement (exact rational equality): ", lams[1] == lams[2]);
lambda = lams[1];
phi = phi_raw / lambda;
D2 = used[2];
Et2 = ellinit([0,0,0,-1156*D2^2,0]);
Lr2 = lfun(Et2, 1) / Om;
w2 = lindep([Lr2, 1]);
q2 = -w2[2]/w2[1];
S2 = sum(a=0, D2-1, kronecker(D2,a) * mseval(M, phi, [a/D2, oo]));
print("twist-2 check: normalised symbol sum = ", S2, " ; recognised L-ratio = ", q2, " ; EXACTLY equal: ", S2 == q2);
\\ ---- Overconvergent lift and p-adic L series ----
nprec = 10;
Mp = mspadicinit(M, 5, nprec, 0);
print("mspadicinit done, n = ", nprec, " (flag 0 = ordinary symbols only; a5=-2 is a 5-adic unit)");
PHI = mstooms(Mp, phi);
print("mstooms done");
mu = mspadicmoments(Mp, PHI);
print("mspadicmoments done");
ser = mspadicseries(mu);
print("SERIES = ", ser);
c0 = polcoef(ser, 0); c1 = polcoef(ser, 1); c2 = polcoef(ser, 2);
c3 = polcoef(ser, 3); c4 = polcoef(ser, 4);
print("c0 = ", c0);
print("c1 = ", c1);
print("c2 = ", c2);
print("c3 = ", c3);
print("c4 = ", c4);
print("v5(c2) = ", valuation(c2));
claim = 3 + 2*5 + 3*5^2 + 5^3 + 4*5^4 + 3*5^5 + 2*5^6 + O(5^30);
print("v5(c2 - claim43338) = ", valuation(c2 - claim), " (>=7 iff agreement with claim mod 5^7)");
print("v5(c0) = ", valuation(c0), " ; v5(c1) = ", valuation(c1));
\\ ---- Item 4: q5 from THIS route's c2 via eq. (3.1) ----
\\ det H digits: kappa5-H compute_H.out (ellpadicregulator, O(5^82)); terms to 5^14.
detH = 3*5^2 + 4*5^3 + 2*5^4 + 5^6 + 3*5^8 + 4*5^9 + 4*5^11 + 5^12 + 2*5^15 + O(5^16);
\\ alpha: unit root of X^2 + 2X + 5 in Z_5 (alpha = 3 mod 5), own lift:
rts = polrootspadic(x^2 + 2*x + 5, 5, 16);
alpha = 0;
for(i=1, #rts, if(Mod(lift(rts[i]), 5) == 3, alpha = rts[i]));
print("alpha (own polrootspadic lift) = ", alpha);
lg = log(6 + O(5^16));
Reg_gamma = detH / lg^2;
q5 = c2 / ((1 - alpha^(-1))^2 * Reg_gamma);
print("Reg_gamma = ", Reg_gamma);
print("q5 = ", q5);
print("v5(q5 - 1) = ", valuation(q5 - 1), " (>=7 iff q5 = 1 mod 5^7)");
print("E34 DONE");

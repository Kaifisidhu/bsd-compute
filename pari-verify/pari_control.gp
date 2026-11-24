\\ PARI-VERIFY CONTROL v2: E56, y^2 = x^3 - 56x, N = 12544, p = 5.
\\ Independent route: PARI overconvergent modular symbols (Pollack-Stevens).
\\ NO Sage/eclib code anywhere in this path. Eigensymbol built by Hecke
\\ eigenvalues (msinit -> msnew -> msfromhecke), NOT msfromell (which
\\ overflows the PARI stack at this level in this build: see log v1).
\\ Normalisation: fixed by the twist-sum identity in msfromell's own doc:
\\   sum_a (D|a) x^+([a/D]-[oo]) = L(E^D,1)/Omega^+,  Omega^+ = E.omega[1],
\\ with L(E^D,1) from PARI lfun (analytic route, independent of symbols).
\\ Expected (archive/Sage): c2 = 221 mod 5^5; fresh n=9 value 18971 mod 5^8.
default(timer, 1);
default(realprecision, 60);
E = ellinit([0,0,0,-56,0]);
N = ellglobalred(E)[1];
print("E56 conductor = ", N, " (expect 12544)");
print("a5 = ", ellap(E,5), " (expect -2)");
print("Omega^+ = E.omega[1] = ", E.omega[1]);
M = msinit(N, 2, 1);
print("msinit done: dim full + space = ", msdim(M));
NW = msnew(M);
print("msnew done: dim new + subspace = ", matsize(NW[1])[2]);
plist = [3,5,7,11,13,29,37,41];
v = vector(#plist, i, [plist[i], ellap(E, plist[i])]);
print("eigenvalue constraints (p, a_p): ", v);
S = msfromhecke(M, v, NW);
print("msfromhecke returned matrix size: ", matsize(S));
\\ Determine coordinate convention: columns in M-basis or in NW-basis?
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
\\ Verify eigenvector at ALL primes p <= 50, p not dividing N (exact, over Q):
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
  Et = ellinit([0,0,0,-56*D^2,0]);
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
\\ Exact-rational verification: with phi normalised, twist sums must equal
\\ the recognised L-ratios exactly. Recognise L-ratio for twist 2 as rational:
D2 = used[2];
Et2 = ellinit([0,0,0,-56*D2^2,0]);
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
target5 = 1 + 4*5 + 3*5^2 + 5^3 + O(5^30);
print("v5(c2 - (1+4*5+3*5^2+5^3)) = ", valuation(c2 - target5), " (>=5 iff agreement mod 5^5)");
target8 = 18971 + O(5^30);
print("v5(c2 - 18971) = ", valuation(c2 - target8), " (>=8 iff agreement with fresh n=9 record mod 5^8)");
print("v5(c0) = ", valuation(c0), " ; v5(c1) = ", valuation(c1));
print("CONTROL DONE");

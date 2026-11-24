\\ Proof of life: E_34 = [0,0,0,-1156,0], y^2 = x^3 - 1156x.
\\ Cross-check values (from the local campaign environment):
\\   conductor 18496, ellrank = [2, 2, 0, [[-16,120],[-2,48]]],
\\   a_5 = -2, elltamagawa = 32 (the known non-minimal-model artefact value;
\\   reproducing it here proves the farm matches the local PARI environment,
\\   it is NOT the true Tamagawa product, which is 16).
E = ellinit([0,0,0,-1156,0]);
print("conductor=", ellglobalred(E)[1]);
print("ellrank=", ellrank(E));
print("a5=", ellap(E,5));
print("ap_table_p_le_50:");
forprime(p=2, 50, print(p, " ", ellap(E,p)));
print("elltamagawa=", elltamagawa(E));
\\ The finite Tamagawa product, prime by prime (elllocalred's c_p):
print("elllocalred_2=", elllocalred(E, 2));
print("elllocalred_17=", elllocalred(E, 17));
print("disc_sign=", sign(E.disc));
print("torsion=", elltors(E));
quit;

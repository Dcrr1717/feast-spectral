"""Paper F, F0: validate the relativistic dipole angular coefficients.

Gates:
G0.1  Hermiticity/symmetry: A_L(k',k,m) == A_L(k,k',m) (cos theta is
      self-adjoint; the coupled matrix must be symmetric).
G0.2  Selection rule: A_L nonzero only for l' = l +- 1 (large) etc.
G0.3  NONREL LIMIT (the decisive gate).  As c->inf the two j-partners
      of a given l' merge.  The nonrel Stark couples l -> l+1 with
      Gaunt c_l = (l+1)/sqrt((2l+1)(2l+3)).  The relativistic large-
      component couplings from a fixed kappa to the two partners of
      l'=l+1 must combine (sum of squares, since the merged nonrel
      level is the incoherent pair) to the Gaunt^2 -- OR, more
      precisely, the coupling to each j-partner weighted by its
      nonrel spinor content reconstructs the Gaunt.  We test the
      clean statement: sum over j'-partners of |A_L(k'->k)|^2 (both
      partners of l'=l+-1) equals the nonrel |c_{l,l'}|^2 summed the
      same way.  Concretely for m=1/2:
        1s1/2 (kappa=-1, l=0): couples to l'=1 partners p1/2(+1),
        p3/2(-2).  sum |A_L|^2 over these must equal (1/3) = c_0^2.
G0.4  A known closed value: <p3/2 m=1/2|cos|s1/2 m=1/2> etc. against a
      hand table if available; here we assert the m=1/2 s-p values.
"""

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from angular import A_L, A_S, kappa_to_jl, kappa_ladder  # noqa: E402


def gaunt(l):
    return (l + 1) / ((2 * l + 1) * (2 * l + 3)) ** 0.5


print("F0 angular gate (m = 1/2):\n")
m2 = 1
lad = kappa_ladder(-1, 8, m2)
print("kappa-ladder:", lad,
      "  (l:", [kappa_to_jl(k)[1] for k in lad], ")\n")

# G0.1 symmetry
print("G0.1 symmetry A_L(k',k)=A_L(k,k'):")
bad = 0
for i, kp in enumerate(lad):
    for k in lad:
        if abs(A_L(kp, k, m2) - A_L(k, kp, m2)) > 1e-12:
            bad += 1
print(f"   asymmetric pairs: {bad}  ({'OK' if bad == 0 else 'FAIL'})")

# G0.2 selection rule (large): l' = l +- 1
print("\nG0.2 selection rule (nonzero A_L only for l'=l+-1):")
viol = 0
for kp in lad:
    for k in lad:
        a = A_L(kp, k, m2)
        lp, l = kappa_to_jl(kp)[1], kappa_to_jl(k)[1]
        if abs(a) > 1e-12 and abs(lp - l) != 1:
            viol += 1
print(f"   violations: {viol}  ({'OK' if viol == 0 else 'FAIL'})")

# G0.3 nonrel-limit sum rule
print("\nG0.3 nonrel-limit sum rule  sum_j' |A_L|^2 = c_l^2:")
for l in (0, 1, 2):
    # source channels of orbital l (both partners if exist for m=1/2)
    srcs = [k for k in lad if kappa_to_jl(k)[1] == l]
    for ks in srcs:
        # target partners with l'=l+1
        tgt = [k for k in lad if kappa_to_jl(k)[1] == l + 1]
        s = sum(A_L(kt, ks, m2) ** 2 for kt in tgt)
        g2 = gaunt(l) ** 2
        ok = abs(s - g2) < 1e-9
        print(f"   from kappa={ks:+d} (l={l},j={abs(ks)-0.5}): "
              f"sum|A_L|^2={s:.6f}  c_l^2={g2:.6f}  "
              f"{'OK' if ok else 'MISMATCH'}")

# G0.4 explicit s-p values
print("\nG0.4 explicit s1/2 -> p couplings (m=1/2):")
for kp in (1, -2):   # p1/2, p3/2
    jl = kappa_to_jl(kp)
    print(f"   <kappa={kp:+d} (p{'1/2' if kp==1 else '3/2'})|cos|kappa=-1 (s1/2)>: "
          f"A_L={A_L(kp,-1,m2):+.6f}  A_S={A_S(kp,-1,m2):+.6f}")

# small-component channel identity check
print("\nsmall-comp A_S selection (uses -kappa parity):")
for kp in lad[:4]:
    row = [f"{A_S(kp,k,m2):+.4f}" for k in lad[:4]]
    print(f"   kappa={kp:+d}: " + "  ".join(row))

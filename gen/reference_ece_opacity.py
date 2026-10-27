"""Independent principal-value quadrature oracle for the X1 FLR correction.

Bornatici et al., Nucl. Fusion 23, 1153 (1983), Eqs. (2.3.68), (3.1.30).
https://doi.org/10.1088/0029-5515/23/9/005
Requires mpmath; used only to reproduce tests/data/ece_x1_flr.dat.
"""
from pathlib import Path
import mpmath as mp
mp.mp.dps = 90

def fq(x, q):
    if x == 0:
        return mp.mpc(1/(q-1), 0)
    g = lambda t: t**(q-1)*mp.exp(-t)
    gx = g(x)
    # The pole subtraction has zero PV integral over [0, 2*x].
    def regular(t):
        if t == x:
            return gx*((q-1)/x-1)
        return (g(t)-gx)/(t-x)
    re = (mp.quad(regular, [0, x, 2*x]) +
          mp.quad(lambda t: g(t)/(t-x), [2*x, mp.inf]))/mp.gamma(q)
    return mp.mpc(re, -mp.pi*gx/mp.gamma(q))

root = Path(__file__).resolve().parents[1]
lines = []
for text in ['0', '.3', '1', '3', '6', '9', '12', '15']:
    x = mp.mpf(text)
    f5, f7 = fq(x, mp.mpf('2.5')), fq(x, mp.mpf('3.5'))
    for density in ['.3', '.8', '1', '1.2', '1.5', '1.9']:
        r = mp.mpf(density)
        correction = (abs(1-r*f7)**2 + mp.mpf(4)/5*r*x*
                      (f5.real-r*(f5.conjugate()*f7).real) +
                      mp.mpf(6)/35*r*r*x*x*abs(f5)**2)
        lines.append(' '.join(mp.nstr(v, 24) for v in [x, r, correction]))
(root/'tests/data/ece_x1_flr.dat').write_text('\n'.join(lines)+'\n')

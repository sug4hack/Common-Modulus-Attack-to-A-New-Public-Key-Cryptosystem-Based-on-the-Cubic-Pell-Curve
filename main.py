from sage.all import *
import pandas as pd

# =======================
# Operasi pada ℙₐ(N)
# =======================

def pell_add(P1, P2, a):
    x1, y1, z1 = P1
    x2, y2, z2 = P2
    x3 = x1 * x2 + a * (y2 * z1 + y1 * z2)
    y3 = x2 * y1 + x1 * y2 + a * z1 * z2
    z3 = y1 * y2 + x2 * z1 + x1 * z2
    return (x3, y3, z3)

def pell_scalar_mult(n, P, a):
    R = (ZN(1), ZN(0), ZN(0))  # Elemen netral
    for bit in bin(n)[2:]:
        R = pell_add(R, R, a)
        if bit == '1':
            R = pell_add(R, P, a)
    return R

# =======================
# Kurva dan Titik Cubic Pell
# =======================

class CubicPellCurve:
    def __init__(self, N, a):
        self.N = N
        self.ZN = Zmod(N)
        self.a = self.ZN(a)

    def __repr__(self):
        return f"CubicPellCurve over Zmod({self.N}) defined by x^3 + {self.a}*y^3 + {self.a^2}*z^3 - 3*{self.a}*x*y*z = 1"

    def is_on_curve(self, P):
        x, y, z = P.x, P.y, P.z
        lhs = x**3 + self.a * y**3 + self.a**2 * z**3 - 3 * self.a * x * y * z
        return lhs == 1

    def zero(self):
        return CubicPellPoint(self, self.ZN(1), self.ZN(0), self.ZN(0))

    def point(self, x, y, z):
        return CubicPellPoint(self, self.ZN(x), self.ZN(y), self.ZN(z))


class CubicPellPoint:
    def __init__(self, curve, x, y, z):
        assert isinstance(curve, CubicPellCurve)
        self.curve = curve
        self.x = x
        self.y = y
        self.z = z

    def __repr__(self):
        return f"({self.x} : {self.y} : {self.z})"

    def __eq__(self, other):
        return (self.x, self.y, self.z) == (other.x, other.y, other.z)

    def __neg__(self):
        x, y, z = self.x, self.y, self.z
        a = self.curve.a
        return CubicPellPoint(
            self.curve,
            x**2 - a * y * z,
            a * z**2 - x * y,
            y**2 - x * z
        )

    def __add__(self, Q):
        x1, y1, z1 = self.x, self.y, self.z
        x2, y2, z2 = Q.x, Q.y, Q.z
        a = self.curve.a
        x3 = x1 * x2 + a * (y2 * z1 + y1 * z2)
        y3 = x2 * y1 + x1 * y2 + a * z1 * z2
        z3 = y1 * y2 + x2 * z1 + x1 * z2
        return CubicPellPoint(self.curve, x3, y3, z3)

    def __radd__(self, other):
        return self if other == 0 else NotImplemented

    def __rmul__(self, n):
        if n < 0:
            return (-self).__rmul__(-n)
        R = self.curve.zero()
        Q = self
        for bit in bin(n)[2:]:
            R = R + R
            if bit == '1':
                R = R + Q
        return R

# =======================
# Parameter dan Kunci
# =======================

p, q, r, s = 922039, 760531, 1, 3
N = p^3 * q^3
ZN = Zmod(N)
e = ZN(190681261905711342654691)

d1 = inverse_mod(e, p**(2*(r-1)) * q**(2*(s-1)) * (p**2+p+1) * (q**2+q+1))
d2 = inverse_mod(e, p**(2*(r-1)) * q**(2*(s-1)) * (p-1)**2 * (q-1)**2)
d3 = inverse_mod(e, p**(2*(r-1)) * q**(2*(s-1)) * (p**2+p+1) * (q-1)**2)
d4 = inverse_mod(e, p**(2*(r-1)) * q**(2*(s-1)) * (p-1)**2 * (q**2+q+1))

# =======================
# Enkripsi
# =======================

xM, yM = 94727413669590175405397, 400429216716868987768230
a = mod((1 - xM**3) * inverse_mod(yM**3, N), N)
curve = CubicPellCurve(N, a)
M = curve.point(xM, yM, ZN(0))
C = e * M
xC, yC, zC = C.x, C.y, C.z

# =======================
# Lifting Hensel
# =======================

Zp, Zq, Zqs = Zmod(p**r), Zmod(q), Zmod(q**s)
R_p = PolynomialRing(Zp, 'x')
x = R_p.gen()
f = xC**3 + x * yC**3 + x**2 * zC**3 - 3 * x * xC * yC * zC - 1
sols_mod_p = [r for r in Zp if f(r) == 0]
x1, x2 = Integer(sols_mod_p[1]), Integer(sols_mod_p[0])

R = PolynomialRing(ZZ, 'y')
y = R.gen()
f_y_ZZ = xC**3 + y * yC**3 + y**2 * zC**3 - 3 * y * xC * yC * zC - 1

R_q = PolynomialRing(Zq, 'y')
f_y_modq = f_y_ZZ.change_ring(Zq)
sols_y_mod_q = [r for r in Zq if f_y_modq(r) == 0]
y1, y2 = sols_y_mod_q[0], sols_y_mod_q[1]

def hensel_lift(f, r0, p, n):
    RZ = f.parent()
    Zp = Zmod(p)
    f_modp = RZ.change_ring(Zp)(f)
    df_modp = f_modp.derivative()

    if f_modp(r0) != 0:
        raise ValueError(f"f({r0}) ≠ 0 mod {p}")
    if df_modp(r0) == 0:
        raise ValueError(f"f'({r0}) ≡ 0 mod {p}, tidak dapat diangkat.")

    r = Integer(r0)
    modulus = p
    for _ in range(2, n+1):
        modulus *= p
        Rmod = Zmod(modulus)
        f_mod = RZ.change_ring(Rmod)(f)
        df_mod = f_mod.derivative()
        num = f_mod(r)
        den = df_mod(r)
        correction = Rmod(num) * Rmod(den).inverse_of_unit()
        r = (r - correction).lift() % modulus

    return r

y1_lifted = hensel_lift(f_y_ZZ, y1, q, s)
y2_lifted = hensel_lift(f_y_ZZ, y2, q, s)

# =======================
# CRT dan Dekripsi
# =======================

ap1, ap2 = Integer(x1), Integer(x2)
aq1, aq2 = Integer(y1_lifted), Integer(y2_lifted)

a1 = crt(ap1, aq1, p**r, q**s)
a2 = crt(ap2, aq2, p**r, q**s)
a3 = crt(ap1, aq2, p**r, q**s)
a4 = crt(ap2, aq1, p**r, q**s)
a_list = [a1, a2, a3, a4]

curve_a = [CubicPellCurve(N, ai) for ai in a_list]
Ci = [curve_a[i].point(xC, yC, zC) for i in range(4)]

def is_cubic_residue(a, n):
    return power_mod(a, (n - 1) // 3, n) == 1

def pilih_Di(ai, p, q, d1, d2, d3, d4):
    is_res_p = is_cubic_residue(ai % p, p)
    is_res_q = is_cubic_residue(ai % q, q)
    if not is_res_p and not is_res_q:
        return d1
    elif is_res_p and is_res_q:
        return d2
    elif not is_res_p and is_res_q:
        return d3
    elif is_res_p and not is_res_q:
        return d4

selected_D = [pilih_Di(ai, p, q, d1, d2, d3, d4) for ai in a_list]

# =======================
# Output
# =======================

for i in range(4):
    print(f"Hasil dekripsi {i+1}:", selected_D[i] * Ci[i])

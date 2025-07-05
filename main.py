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
        self.curve = curve
        self.x = x
        self.y = y
        self.z = z
        assert isinstance(curve, CubicPellCurve)

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
        if other == 0:
            return self
        else:
            return NotImplemented

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


# Parameter
p = 922039
q = 760531
r = 1
s = 3

# Modulus
N = 405601968528411801552349
ZN = Zmod(N)

# Kunci publik
e = ZN(190681261905711342654691)

#Kunci privat
d1 = inverse_mod(e, p**(2*(r-1)) * q**(2*(s-1)) * (p**2+p+1) * (q**2+q+1))
d2 = inverse_mod(e, p**(2*(r-1)) * q**(2*(s-1)) * (p-1)**2 * (q-1)**2)
d3 = inverse_mod(e, p**(2*(r-1)) * q**(2*(s-1)) * (p**2+p+1) * (q-1)**2)
d4 = inverse_mod(e, p**(2*(r-1)) * q**(2*(s-1)) * (p-1)**2 * (q**2+q+1))

# Plaintext
xM = 94727413669590175405397
yM = 400429216716868987768230
print("Plaintext (xM, yM) =", (xM, yM))


#Hitung a = (1-xM**3)/(yM**3) Mod N
a = mod((1 - xM**3) * inverse_mod(yM**3, N), N)
print("a =", a)

#curve
curve = CubicPellCurve(N, a)

# Titik pesan M
M = curve.point(xM, yM, ZN(0))
print("Point M =", M)

#Enkripsi pesan M
C = e*M
print("Ciphertext C =", C)

xC = C.x
yC = C.y
zC = C.z

pr = p**r
Zp = Zmod(pr)

R = PolynomialRing(Zp, 'x')
x = R.gen()

f = xC**3 + x * yC**3 + x**2 * zC**3 - 3 * x * xC * yC * zC - 1

sols_mod_p = [r for r in Zp if f(r) == 0]

print("Solution x mod p:", sols_mod_p)
x1 = sols_mod_p[1]
x2 = sols_mod_p[0]  
  
qs = q**s
Zq = Zmod(q)
Zqs = Zmod(qs)

R = PolynomialRing(ZZ, 'y')
y = R.gen()
f_y_ZZ = xC**3 + y * yC**3 + y**2 * zC**3 - 3 * y * xC * yC * zC - 1

R_q = PolynomialRing(Zq, 'y')
y_q = R_q.gen()
f_y_modq = f_y_ZZ.change_ring(Zq)

sols_y_mod_q = [r for r in Zq if f_y_modq(r) == 0]
print("Solution y mod q:", sols_y_mod_q)

def hensel_lift(f, r0, p, n):
    RZ = f.parent()
    fZ = f

    Zp = Zmod(p)
    f_modp = RZ.change_ring(Zp)(fZ)
    df_modp = f_modp.derivative()

    if f_modp(r0) != 0:
        raise ValueError(f"f({r0}) ≠ 0 mod {p}")
    if df_modp(r0) == 0:
        raise ValueError(f"f'({r0}) ≡ 0 mod {p}, can't lifted.")

    r = Integer(r0)
    modulus = p
    for i in range(2, n+1):
        modulus *= p
        Rmod = Zmod(modulus)
        f_mod = RZ.change_ring(Rmod)(fZ)
        df_mod = f_mod.derivative()

        num = f_mod(r)
        den = df_mod(r)
        den_mod = Rmod(den)

        if not den_mod.is_unit():
            raise ValueError(f"f'({r}) not unit mod {modulus}")

        correction = Rmod(num) * den_mod.inverse_of_unit()
        r = (r - correction).lift() % modulus

    return r

# Lifting solusi
y1 = sols_y_mod_q[0]
y2 = sols_y_mod_q[1]

y1_lifted = hensel_lift(f_y_ZZ, y1, q, s)
y2_lifted = hensel_lift(f_y_ZZ, y2, q, s)

print("y1 lifted to q^s:", y1_lifted)
print("y2 lifted to q^s:", y2_lifted)

ap1 = Integer(x1)
ap2 = Integer(x2)
aq1 = Integer(y1_lifted)
aq2 = Integer(y2_lifted)

# CRT
a1 = crt(ap1, aq1, p**r, q**s)
a2 = crt(ap2, aq2, p**r, q**s)
a3 = crt(ap1, aq2, p**r, q**s)
a4 = crt(ap2, aq1, p**r, q**s)

print("a1 =", a1)
print("a2 =", a2)
print("a3 =", a3)
print("a4 =", a4)

curve_a1 = CubicPellCurve(N, a1)
curve_a2 = CubicPellCurve(N, a2)
curve_a3 = CubicPellCurve(N, a3)
curve_a4 = CubicPellCurve(N, a4)

C_a1 = curve_a1.point(xC, yC, zC)
C_a2 = curve_a2.point(xC, yC, zC)
C_a3 = curve_a3.point(xC, yC, zC)
C_a4 = curve_a4.point(xC, yC, zC)
Ci = [C_a1, C_a2, C_a3, C_a4]

def is_cubic_residue(a, n):
    return power_mod(a, (n - 1) // 3, n) == 1

def choose_Di(ai, p, q, d1, d2, d3, d4):
    ai_mod_p = ai % p
    ai_mod_q = ai % q

    is_res_p = is_cubic_residue(ai_mod_p, p)
    is_res_q = is_cubic_residue(ai_mod_q, q)

    if not is_res_p and not is_res_q:
        return d1
    elif is_res_p and is_res_q:
        return d2
    elif not is_res_p and is_res_q:
        return d3
    elif is_res_p and not is_res_q:
        return d4


a_list = [a1, a2, a3, a4]
selected_D = [choose_Di(ai, p, q, d1, d2, d3, d4) for ai in a_list]
print("Selected D values:", selected_D)

for i in range(4):
    print(f"decryption result {i+1}:", selected_D[i]*Ci[i])








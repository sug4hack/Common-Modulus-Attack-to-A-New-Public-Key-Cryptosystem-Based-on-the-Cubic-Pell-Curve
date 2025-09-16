from sage.all import *


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


def _ABC_from_point(C, N):
    ZN = Zmod(N)
    x, y, z = ZN(C.x), ZN(C.y), ZN(C.z)
    A = z**3
    B = y**3 - 3*x*y*z
    Cc = x**3 - 1
    return Integer(A), Integer(B), Integer(Cc)

def _solve_linear_congruence(a, b, n):
    g = gcd(a, n)
    if b % g != 0:
        return []
    a1, b1, n1 = a // g, b // g, n // g
    x0 = (Integer(b1) * inverse_mod(Integer(a1), Integer(n1))) % n1
    return [(x0 + k*n1) % n for k in range(g)]

def recover_a(C1, C2, N):
    print(f"C1: {C1}, C2: {C2}, N: {N}")
    A1, B1, C1c = _ABC_from_point(C1, N)
    A2, B2, C2c = _ABC_from_point(C2, N)
    print(f"A1: {A1}, B1: {B1}, C1: {C1c}")
    print(f"A2: {A2}, B2: {B2}, C2: {C2c}") 

    D = (A2*B1 - A1*B2) % N
    K = (A2*C1c - A1*C2c) % N
     
    print(f"D: {D}, K: {K}")

    sols = _solve_linear_congruence(int(D), int((-K) % N), int(N))
    print(f"Initial sols: {sols}")

    if not sols:
        D2 = (C2c*A1 - C1c*A2) % N
        K2 = (C2c*B1 - C1c*B2) % N
        sols = []
        if K2 % N == 0:
            sols.append(0)
        try:
            invD2 = inverse_mod(Integer(D2), Integer(N))
            sols.append(int((-K2 * invD2) % N))
        except Exception:
            pass

    return sols



# Parameters
p = 922039
q = 760531
r = 1
s = 3
N = p**r * q**s
a_actual = 402129345655132093067351
curve = CubicPellCurve(N, a_actual)

P = curve.point(
    94727413669590175405397,
    400429216716868987768230,
    0
)

e1 = 190681261905711342654691
e2 = 41920150622089930848871

C1 = e1 * P
C2 = e2 * P

# Step 1: Get 4 candidate a values
a_candidates = recover_a(C1, C2, N)
# Step 2 & 3: Check which candidate makes both C1 and C2 valid points
valid_curve = None
for a in a_candidates:
    candidate_curve = CubicPellCurve(N, a)
    if candidate_curve.is_on_curve(C1) and candidate_curve.is_on_curve(C2):
        valid_curve = candidate_curve
        break

# Step 4: Reconstruct ciphertexts on valid curve
C1_fixed = CubicPellPoint(valid_curve, C1.x, C1.y, C1.z)
C2_fixed = CubicPellPoint(valid_curve, C2.x, C2.y, C2.z)

# Step 5: Solve x*e1 + y*e2 = 1 using extended Euclidean algorithm
g, x, y = xgcd(e1, e2)
assert g == 1

# Step 6: Recover the original message point
P_recovered = x * C1_fixed + y * C2_fixed

#print all parameters and results
print(f"p: {p}")
print(f"q: {q}")
print(f"N: {N}")
print(f"a_actual: {a_actual}")
print(f"e1: {e1}")
print(f"e2: {e2}")
print(f"C1: {C1}")
print(f"C2: {C2}")
print(f"x: {x}")
print(f"y: {y}")
print(f"Recovered a candidates: {a_candidates}")
print(f"Recovered Point P: {P_recovered}")
#is P_recovered equal to P?
print(f"Is recovered P equal to original P? {P_recovered == P}")

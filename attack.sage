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


def recover_a_fast(xC, yC, zC, N):
    ZN = Zmod(N)
    x, y, z = ZN(xC), ZN(yC), ZN(zC)

    A = z**3
    B = y**3 - 3 * x * y * z
    C = x**3 - 1

    disc = B**2 - 4 * A * C
    sqrt_disc = disc.sqrt()
    inv_2A = (2 * A).inverse_of_unit()
    a1 = (-B + sqrt_disc) * inv_2A
    a2 = (-B - sqrt_disc) * inv_2A

    return a1, a2


# Parameters
N = 405601968528411801552349
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
a1, a2 = recover_a_fast(C1.x, C1.y, C1.z, N)
a3, a4 = recover_a_fast(C2.x, C2.y, C2.z, N)
a_candidates = [a1, a2, a3, a4]

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
print(P_recovered == P)
print(P_recovered)

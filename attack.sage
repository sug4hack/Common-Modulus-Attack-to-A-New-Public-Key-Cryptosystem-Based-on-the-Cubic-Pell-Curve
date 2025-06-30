###############################
# CLASS DEFINITIONS
###############################

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

###############################
# FUNGSI: Recovery Parameter a
###############################
def recover_a_fast(xC, yC, zC, N):
    ZN = Zmod(N)
    x, y, z = ZN(xC), ZN(yC), ZN(zC)

    A = z**3
    B = y**3 - 3 * x * y * z
    C = x**3 - 1

    disc = B**2 - 4 * A * C
    try:
        sqrt_disc = disc.sqrt()
    except ValueError:
        raise ValueError("There is no square root of the discriminant mod N")

    inv_2A = (2 * A).inverse_of_unit()
    a1 = (-B + sqrt_disc) * inv_2A
    a2 = (-B - sqrt_disc) * inv_2A

    return a1, a2

###############################
# FUNGSI: Identifikasi kandidat a yang benar
###############################
def identify_correct_a(C, a_candidates, N):
    for a in a_candidates:
        curve_candidate = CubicPellCurve(N, a)
        if curve_candidate.is_on_curve(C):
            return a
    return None

###############################
# COMMON MODULUS ATTACK DEMO
###############################

# Parameter sistem
N = 405601968528411801552349
a_actual = 402129345655132093067351
curve = CubicPellCurve(N, a_actual)

# Titik pesan P
P = curve.point(
    94727413669590175405397,
    400429216716868987768230,
    0
)

# Dua ciphertext dengan eksponen berbeda
e1 = 190681261905711342654691
e2 =  41920150622089930848871

C1 = e1 * P
C2 = e2 * P
print(curve.is_on_curve(C2))

# Serangan dimulai — Penyerang tidak tahu a
# Langkah 1: Recover kandidat a dari C1
a1, a2 = recover_a_fast(C1.x, C1.y, C1.z, N)
print("Kandidat a:", a1, a2)

# Langkah 2: Identifikasi a yang sesuai
a_found = identify_correct_a(C2, [a1, a2], N)
print("a yang berhasil direkonstruksi:", a_found)

# Langkah 3: Buat ulang kurva berdasarkan a yang ditemukan
if a_found is None:
    raise RuntimeError("Gagal menemukan nilai a yang valid dari ciphertext.")
else:
    curve_recovered = CubicPellCurve(N, a_found)

# Langkah 4: Gunakan Extended Euclidean Algorithm pada (e1, e2)
# x*e1 + y*e2 = 1 → x dan y diketahui (hasil luar dari egcd)
x = -1237081083026857444837
y =  5627083359451083671408


# Langkah 5: Pulihkan P dari C1 dan C2
P_recovered = x * C1 + y * C2
print("Apakah P berhasil dipulihkan?", P_recovered == P)

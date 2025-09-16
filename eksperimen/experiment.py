from sage.all import *
import time, statistics, csv, json, math, tracemalloc
from dataclasses import dataclass, asdict
from typing import List, Dict, Any, Optional

class CubicPellCurve:
    def __init__(self, N, a):
        self.N  = Integer(N)
        self.ZN = Zmod(N)
        self.a  = self.ZN(a)

    def is_on_curve(self, P):
        x, y, z = P.x, P.y, P.z
        lhs = x**3 + self.a*y**3 + (self.a**2)*z**3 - 3*self.a*x*y*z
        return lhs == 1

    def zero(self):
        return CubicPellPoint(self, self.ZN(1), self.ZN(0), self.ZN(0))

    def point(self, x, y, z):
        return CubicPellPoint(self, self.ZN(x), self.ZN(y), self.ZN(z))


class CubicPellPoint:
    def __init__(self, curve, x, y, z):
        self.curve = curve
        self.x, self.y, self.z = x, y, z

    def __repr__(self):
        return f"({self.x}:{self.y}:{self.z})"

    def __eq__(self, other):
        return (self.x, self.y, self.z) == (other.x, other.y, other.z)

    def __neg__(self):
        x, y, z = self.x, self.y, self.z
        a = self.curve.a
        return CubicPellPoint(self.curve, x**2 - a*y*z, a*z**2 - x*y, y**2 - x*z)

    def __add__(self, Q):
        x1, y1, z1 = self.x, self.y, self.z
        x2, y2, z2 = Q.x, Q.y, Q.z
        a = self.curve.a
        
        x3 = x1*x2 + a*(y2*z1 + y1*z2)
        y3 = x2*y1 + x1*y2 + a*(z1*z2)
        z3 = y1*y2 + x2*z1 + x1*z2
        return CubicPellPoint(self.curve, x3, y3, z3)

    def __radd__(self, other):
        return self if other == 0 else NotImplemented

    def __rmul__(self, n):
        if n < 0:
            return (-self).__rmul__(-n)
        R = self.curve.zero()
        Q = self
        for bit in bin(int(n))[2:]:
            R = R + R
            if bit == '1':
                R = R + Q
        return R


# ====================== Attack Utilities ======================
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
    A1, B1, C1c = _ABC_from_point(C1, N)
    A2, B2, C2c = _ABC_from_point(C2, N)

    D = (A2*B1 - A1*B2) % N
    K = (A2*C1c - A1*C2c) % N
    sols = _solve_linear_congruence(int(D), int((-K) % N), int(N))

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

    good = []
    for aa in sols:
        cc = CubicPellCurve(N, aa)
        if cc.is_on_curve(C1) and cc.is_on_curve(C2):
            good.append(aa % N)

    return list(dict.fromkeys(good))


def random_prime_bits(bits):
    lo = 1 << (bits - 1)
    hi = (1 << bits) - 1
    return Integer(random_prime(hi, lbound=lo))

def random_coprime_to(N):
    while True:
        x = randint(2, N - 2)
        if gcd(x, N) == 1:
            return x


# ====================== Benchmark Harness ======================
@dataclass
class TrialResult:
    N_bits: int
    trial_idx: int
    small_exponents: bool
    success: bool
    p_bits: int
    q_bits: int

    t_keygen_s: float
    t_encrypt_s: float
    t_recover_a_s: float
    t_cma_s: float
    t_total_s: float
    peak_mem_kb: int

    e1_bits: int
    e2_bits: int
    u_bits: int
    v_bits: int

    p: int
    q: int
    N: int
    a_true: int
    e1: int
    e2: int
    u: int
    v: int
    
    C1x: int
    C1y: int
    C1z: int
    C2x: int
    C2y: int
    C2z: int

    a_found: int


def gen_instance_attack_z0_bench(N_bits: int, small_exponents: bool = True, verbose: bool = False) -> TrialResult:
    tracemalloc.start()
    t_wall0 = time.perf_counter()

    # --- Keygen ---
    t0 = time.perf_counter()
    p_bits = N_bits // 2
    q_bits = N_bits - p_bits
    p = random_prime_bits(p_bits)
    q = random_prime_bits(q_bits)
    N = p * q
    ZN = Zmod(N)

    x = random_coprime_to(N)
    y = random_coprime_to(N)
    a_true = int(mod((1 - x**3) * inverse_mod(y**3, N), N))

    curve = CubicPellCurve(N, a_true)
    P = curve.point(x, y, 0)
    assert curve.is_on_curve(P), "P must be on-curve"
    t1 = time.perf_counter()
    t_keygen = t1 - t0

    # --- Encrypt ---
    t0 = time.perf_counter()
    if small_exponents:
        e1, e2 = 65537, 65539
    else:
        while True:
            e1 = randint(2, N - 2)
            e2 = randint(2, N - 2)
            if gcd(e1, e2) == 1:
                break
    C1 = int(e1) * P
    C2 = int(e2) * P
    t1 = time.perf_counter()
    t_encrypt = t1 - t0

    # --- Recover a ---
    t0 = time.perf_counter()
    a_candidates = recover_a(C1, C2, N)
    valid_curve = None
    a_found = -1
    for a_can in a_candidates:
        candidate_curve = CubicPellCurve(N, a_can)
        if candidate_curve.is_on_curve(C1) and candidate_curve.is_on_curve(C2):
            valid_curve = candidate_curve
            a_found = int(a_can)
            break
    t1 = time.perf_counter()
    t_recover_a = t1 - t0

    t0 = time.perf_counter()
    g, u, v = xgcd(e1, e2)
    assert g == 1
    u_bits = Integer(abs(u)).nbits()
    v_bits = Integer(abs(v)).nbits()

    success = False
    if valid_curve is not None:
        C1f = CubicPellPoint(valid_curve, C1.x, C1.y, C1.z)
        C2f = CubicPellPoint(valid_curve, C2.x, C2.y, C2.z)
        P_rec = int(u) * C1f + int(v) * C2f
        success = (P_rec == P)
    t1 = time.perf_counter()
    t_cma = t1 - t0

    t_total = time.perf_counter() - t_wall0
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    return TrialResult(
        N_bits=N_bits,
        trial_idx=0,
        small_exponents=small_exponents,
        success=bool(success),
        p_bits=p.nbits(), q_bits=q.nbits(),

        t_keygen_s=t_keygen,
        t_encrypt_s=t_encrypt,
        t_recover_a_s=t_recover_a,
        t_cma_s=t_cma,
        t_total_s=t_total,
        peak_mem_kb=int(peak // 1024),

        e1_bits=Integer(e1).nbits(),
        e2_bits=Integer(e2).nbits(),
        u_bits=int(u_bits),
        v_bits=int(v_bits),

        p=int(p),
        q=int(q),
        N=int(N),
        a_true=int(a_true),
        e1=int(e1),
        e2=int(e2),
        u=int(u),
        v=int(v),
        C1x=int(C1.x), C1y=int(C1.y), C1z=int(C1.z),
        C2x=int(C2.x), C2y=int(C2.y), C2z=int(C2.z),

        a_found=int(a_found),
    )

def run_benchmark(
    sizes: List[int],
    trials_per_size: int = 3,
    small_exponents: bool = True,
    seed: Optional[int] = 1337,
    csv_path: str = "cpc_attack_bench.csv",
    json_path: str = "cpc_attack_summary.json",
    verbose: bool = True
) -> Dict[str, Any]:
    if seed is not None:
        set_random_seed(seed)

    all_trials: List[TrialResult] = []
    rows: List[Dict[str, Any]] = []

    for nb in sizes:
        for ti in range(trials_per_size):
            res = gen_instance_attack_z0_bench(nb, small_exponents=small_exponents, verbose=False)
            res.trial_idx = ti
            all_trials.append(res)

            row = asdict(res)
            rows.append(row)

            if verbose:
                print(f"[N~{nb}][trial {ti+1}/{trials_per_size}] "
                      f"succ={res.success}  total={res.t_total_s:.3f}s  "
                      f"keygen={res.t_keygen_s:.3f}s, enc={res.t_encrypt_s:.3f}s, "
                      f"recA={res.t_recover_a_s:.6f}s, cma={res.t_cma_s:.3f}s, "
                      f"peak={res.peak_mem_kb} KB")

    # Write CSV (per-trial)
    fieldnames = list(rows[0].keys()) if rows else []
    if fieldnames:
        with open(csv_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fieldnames)
            w.writeheader()
            w.writerows(rows)

    def stats(vals):
        return {
            "min": float(min(vals)),
            "max": float(max(vals)),
            "mean": float(statistics.mean(vals)),
            "median": float(statistics.median(vals)),
            "stdev": float(statistics.pstdev(vals)) if len(vals) > 1 else 0.0
        }

    summary: Dict[str, Any] = {"params": {
        "sizes": sizes,
        "trials_per_size": trials_per_size,
        "small_exponents": small_exponents,
        "seed": seed
    }, "by_size": {}}

    for nb in sizes:
        group = [t for t in all_trials if t.N_bits == nb]
        if not group:
            continue
        summary["by_size"][str(nb)] = {
            "success_rate": float(sum(1 for t in group if t.success) / len(group)),
            "t_total_s": stats([t.t_total_s for t in group]),
            "t_keygen_s": stats([t.t_keygen_s for t in group]),
            "t_encrypt_s": stats([t.t_encrypt_s for t in group]),
            "t_recover_a_s": stats([t.t_recover_a_s for t in group]),
            "t_cma_s": stats([t.t_cma_s for t in group]),
            "peak_mem_kb": stats([t.peak_mem_kb for t in group]),
            "e_bits": {
                "e1_bits": stats([t.e1_bits for t in group]),
                "e2_bits": stats([t.e2_bits for t in group]),
                "u_bits": stats([t.u_bits for t in group]),
                "v_bits": stats([t.v_bits for t in group]),
            }
        }

    with open(json_path, "w") as f:
        json.dump(summary, f, indent=2)

    if verbose:
        print("\n=== SUMMARY (per size) ===")
        for nb in sizes:
            if str(nb) not in summary["by_size"]:
                continue
            s = summary["by_size"][str(nb)]
            print(f"N~{nb} bits | success={s['success_rate']*100:.1f}% | "
                  f"total mean={s['t_total_s']['mean']:.3f}s | "
                  f"peak mem mean={s['peak_mem_kb']['mean']:.0f} KB")

    return summary


if __name__ == "__main__":
    SIZES = [256, 512, 1024, 2048, 4096]
    TRIALS = 10                   
    SMALL_EXP = False

    _ = run_benchmark(
        sizes=SIZES,
        trials_per_size=TRIALS,
        small_exponents=SMALL_EXP,
        seed=20250913,
        csv_path="cpc_attack_bench.csv",
        json_path="cpc_attack_summary.json",
        verbose=True
    )

#!/usr/bin/env sage
# -*- coding: utf-8 -*-
"""
Statistical attack on a combining gamma generator

Attack outline (Siegenthaler's attack):
  1) Compute Walsh-Hadamard coefficients of combining function f
  2) Find ALL affine statistical analogs of f and their coincidence probabilities
  3) Compute correlation immunity cor(f)
  4) For each chosen analog g (depends on k out of n registers):
     - Compute required material T (with delta = 0.01)
     - Enumerate all 2^(mk) candidate initial states for those k registers
     - For each candidate, generate T bits via g, compare with gamma (Hamming distance)
     - The candidate with minimum Hamming distance is the true initial state
  5) Recover all 6 registers (possibly using multiple analogs + brute force)
  6) Verify by generating gamma from the full recovered state
"""

import time
import re
import math

# IMPORTANT: While sage preparsing converts ^ to ** (exponentiation) and ^^ to ^.
# To avoid any issues, we define explicit XOR and popcount functions
# using Python's int methods directly.
def bxor(a, b):
    """Bitwise XOR, safe from Sage's ^ preparser."""
    return int.__xor__(int(a), int(b))

def popcount(x):
    """Count number of 1-bits in integer x."""
    return bin(int(x)).count('1')

# ============================================================
# Step 1. Parse the combining function and build truth table
# ============================================================
# The combining function f(x0,...,x5) acts on outputs of 6 LFSRs.
# It is defined over GF(2) — all arithmetic is mod 2 (XOR for +, AND for *).
# We build a truth table: for each of the 2^6 = 64 input combinations,
# store f(x0,...,x5) in a lookup array. This enables O(1) evaluation later.

print("=" * 60)
print("Step 1. Combining function and truth table")
print("=" * 60)

n = 6   # number of LFSRs / variables in f
m = 10  # length of each LFSR (bits)

# Parse the combining function from fileS
with open('combinfuncV10.txt', 'r') as file:
    func_str = file.read().strip()

# Use Sage's BooleanPolynomialRing for symbolic representation
B = BooleanPolynomialRing(n, 'x', order='degrevlex')
B.inject_variables()
xvars = B.gens()

# Parse: split by " + ", extract x_{i} indices, build polynomial
terms = func_str.split(' + ')
f_poly = B(0)
for term in terms:
    term = term.strip()
    if term == '1':
        f_poly += B(1)
        continue
    var_indices = re.findall(r'x_\{(\d+)\}', term)
    if not var_indices:
        continue
    monomial = B(1)
    for v in var_indices:
        monomial *= xvars[int(v)]
    f_poly += monomial

print(f"f = {f_poly}")
print(f"deg(f) = {f_poly.degree()}")

# Build truth table using Sage's BooleanFunction
# f_tt[i] = f(bit_0(i), bit_1(i), ..., bit_5(i))
# where i = x0 + 2*x1 + 4*x2 + 8*x3 + 16*x4 + 32*x5
from sage.crypto.boolean_function import BooleanFunction
bf = BooleanFunction(f_poly)
f_tt = list(bf.truth_table(format='int'))  # list of 0/1 of length 64

wt_f = sum(f_tt)
print(f"Weight of f: {wt_f} / {2**n}  (balance = {float(wt_f / 2**n):.4f})")

# ============================================================
# Step 2. Walsh-Hadamard Transform
# ============================================================
# The Walsh-Hadamard coefficient f_hat(alpha) is defined as:
#   f_hat(alpha) = sum_{x in V_n} (-1)^{f(x) XOR <alpha, x>}
# where <alpha, x> = XOR of (alpha_i AND x_i) is the inner product mod 2.
#
# These coefficients reveal the correlation between f and each affine function.
# Sage's BooleanFunction computes this directly.

print("\n" + "=" * 60)
print("Step 2. Walsh-Hadamard Transform")
print("=" * 60)

wh = list(bf.walsh_hadamard_transform())

print("All Walsh-Hadamard coefficients:")
for alpha in range(2**n):
    bits = ''.join(str((alpha >> j) & 1) for j in range(n))
    vars_str = ' + '.join(f'x{j}' for j in range(n) if (alpha >> j) & 1)
    if not vars_str:
        vars_str = "const"
    status = f"f_hat = {wh[alpha]:+4d}" if wh[alpha] != 0 else "f_hat =    0"
    marker = " <-- nonzero" if wh[alpha] != 0 else ""
    print(f"  alpha = {bits}  ({vars_str:20s}): {status}{marker}")

# ============================================================
# Step 3. Affine statistical analogs and coincidence probabilities
# ============================================================
# For each alpha with f_hat(alpha) != 0:
#   - If f_hat(alpha) > 0: the analog is g = <alpha, x>       (c = 0)
#   - If f_hat(alpha) < 0: the analog is g = <alpha, x> + 1   (c = 1)
#   - Coincidence probability: Pr(f = g) = 1/2 * (1 + |f_hat(alpha)| / 2^n)
#   - The "epsilon" (bias) is: eps = |f_hat(alpha)| / 2^n
#   - Error probability: p = Pr(f != g) = 1/2 * (1 - eps)
#
# The analog g depends on k = wt(alpha) variables (registers).
# The fewer variables, then attack is "cheaper" (fewer candidates to enumerate).

print("\n" + "=" * 60)
print("Step 3. Affine statistical analogs")
print("=" * 60)

analogs = []
for alpha in range(2**n):
    fh = wh[alpha]
    if fh == 0:
        continue

    c = 0 if fh > 0 else 1
    eps = float(abs(fh)) / float(2**n)
    prob = 0.5 * (1.0 + eps)    # Pr(f = g)
    err = 0.5 * (1.0 - eps)     # Pr(f != g)

    var_indices = [j for j in range(n) if (alpha >> j) & 1]
    k = len(var_indices)

    # Build human-readable string for g
    if k == 0:
        g_str = str(c)
    else:
        parts = [f"x{j}" for j in var_indices]
        g_str = " + ".join(parts)
        if c == 1:
            g_str += " + 1"

    analogs.append({
        'alpha': alpha, 'c': c, 'fh': fh, 'eps': eps,
        'prob': prob, 'err': err,
        'vars': var_indices, 'k': k, 'g_str': g_str
    })

print(f"\n{'Analog g':30s} {'k':>3s} {'eps':>8s} {'Pr(f=g)':>8s} {'Pr(f!=g)':>9s} {'f_hat':>6s}")
print("-" * 70)
for a in sorted(analogs, key=lambda x: (x['k'], -x['eps'])):
    print(f"  g = {a['g_str']:24s} {a['k']:3d} {a['eps']:8.4f} "
          f"{a['prob']:8.4f} {a['err']:9.4f} {a['fh']:+6d}")
print(f"\nTotal nonzero coefficients: {len(analogs)}")

# ============================================================
# Step 4. Correlation immunity
# ============================================================
# cor(f) = max k such that f_hat(alpha) = 0 for every nonzero alpha
# with Hamming weight wt(alpha) <= k.
# If no such k exists, cor(f) = 0.
#
# High correlation immunity means no low-weight analogs exist,
# making the attack harder (must enumerate more registers simultaneously).

print("\n" + "=" * 60)
print("Step 4. Correlation immunity")
print("=" * 60)

cor_f = 0
for k_test in range(1, n + 1):
    all_zero = True
    for alpha in range(1, 2**n):
        if bin(alpha).count('1') <= k_test and wh[alpha] != 0:
            all_zero = False
            break
    if all_zero:
        cor_f = k_test
    else:
        break

print(f"cor(f) = {cor_f}")
if cor_f == 0:
    print("  -> f has weight-1 nonzero WHT coefficients => single-register analogs exist")
else:
    print(f"  -> minimum weight of a useful analog is {cor_f + 1}")

# ============================================================
# Step 5. Read the gamma sequence
# ============================================================
print("\n" + "=" * 60)
print("Step 5. Reading gamma sequence")
print("=" * 60)

with open('gamma2V10.txt', 'r') as file:
    gamma_str = file.read().strip()
gamma = [int(b) for b in gamma_str]
N = len(gamma)
n_zeros = gamma.count(0)
n_ones = gamma.count(1)
print(f"Gamma length: {N}")
print(f"Number of 0s: {n_zeros} ({float(100*n_zeros/N):.1f}%)")
print(f"Number of 1s: {n_ones} ({float(100*n_ones/N):.1f}%)")

# ============================================================
# Step 6. LFSR definitions
# ============================================================
# Each LFSR has length m = 10 with a given feedback polynomial.
# The state is represented as a 10-bit integer.
# At each step: output = LSB (bit 0), feedback = XOR of tap positions,
# shift right, insert feedback at MSB (bit 9).
#
# The "tap mask" encodes which coefficients c_i are nonzero in
# p(x) = x^10 + c_9*x^9 + ... + c_1*x + c_0.
# For p(x) = x^10 + x^3 + 1: c_0 = 1, c_3 = 1 => tap_mask = (1<<0)|(1<<3) = 9.

print("\n" + "=" * 60)
print("Step 6. LFSR parameters")
print("=" * 60)

lfsr_info = [
    ("R1", "x^10 + x^3 + 1",                (1 << 0) | (1 << 3)),
    ("R2", "x^10 + x^7 + 1",                (1 << 0) | (1 << 7)),
    ("R3", "x^10 + x^4 + x^3 + x + 1",     (1 << 0) | (1 << 1) | (1 << 3) | (1 << 4)),
    ("R4", "x^10 + x^8 + x^3 + x^2 + 1",   (1 << 0) | (1 << 2) | (1 << 3) | (1 << 8)),
    ("R5", "x^10 + x^8 + x^4 + x^3 + 1",   (1 << 0) | (1 << 3) | (1 << 4) | (1 << 8)),
    ("R6", "x^10 + x^8 + x^5 + x + 1",     (1 << 0) | (1 << 1) | (1 << 5) | (1 << 8)),
]

tap_masks = []
for name, poly_str, tap_mask in lfsr_info:
    tap_masks.append(tap_mask)
    print(f"  {name}: {poly_str:35s}  tap_mask = {tap_mask}")

# ============================================================
# Step 7. Compute required material T for each analog
# ============================================================
# From the theory (Proposition 1):
#   T = ceil(8 * eps^{-2} * ln(2^{mk} / delta))
# where:
#   m = 10 (register length), k = number of registers in the analog,
#   delta = 0.01 (error probability of the attack),
#   eps = |f_hat(alpha)| / 2^n (bias of the analog).

print("\n" + "=" * 60)
print("Step 7. Required material T for each analog")
print("=" * 60)

delta = 0.01

for a in analogs:
    k = a['k']
    eps = a['eps']
    if k == 0:
        a['T'] = 0
        continue
    T_val = math.ceil(8 * eps**(-2) * math.log(2**(m * k) / delta))
    a['T'] = min(T_val, N)

print(f"delta = {delta}")
print(f"\n{'Analog g':30s} {'k':>3s} {'T':>8s} {'Registers':>20s} {'Candidates':>12s}")
print("-" * 80)
for a in sorted(analogs, key=lambda x: (x['k'], x.get('T', 0))):
    if a['k'] == 0:
        continue
    reg_str = ', '.join(f"R{j+1}" for j in a['vars'])
    cand = (2**m - 1) ** a['k']
    print(f"  g = {a['g_str']:24s} {a['k']:3d} {a['T']:8d} {reg_str:>20s} {cand:>12d}")

# ============================================================
# Step 8. Precompute LFSR sequences
# ============================================================
# For efficiency, we precompute the output sequence of each LFSR for
# every possible nonzero initial state (1..1023). Each sequence is stored
# as a big Python integer where bit i = output at time i.
# This allows ultra-fast Hamming distance computation via XOR + popcount.

print("\n" + "=" * 60)
print("Step 8. Precomputing LFSR sequences")
print("=" * 60)

max_T = max((a['T'] for a in analogs if a['k'] > 0), default=N)
# Add margin for possible offset issues; cap at N
seq_len = min(max_T + 2, N + 2)

def gen_lfsr_bigint(init_state, tap_mask, length):
    """Generate LFSR sequence as a big integer (bit i = output at time i).

    At each step:
      - Output the LSB (s_0)
      - Compute feedback = XOR of bits at tap positions
      - Shift right, insert feedback at MSB position (bit m-1)

    The sequence s_0, s_1, s_2, ... corresponds to successive outputs.
    """
    state = init_state
    mask_m = (1 << m) - 1
    result = 0
    for i in range(length):
        if state & 1:
            result |= (1 << i)
        fb = bin(state & tap_mask).count('1') & 1
        state = ((state >> 1) | (fb << (m - 1))) & mask_m
    return result

print(f"Sequence length: {seq_len}")
print(f"States per register: {2**m - 1}")
t0 = time.time()

# lfsr_seqs[j][s] = big integer of sequence for register j, initial state s
lfsr_seqs = []
for j in range(n):
    seqs = {}
    tap = tap_masks[j]
    for s in range(1, 2**m):
        seqs[s] = gen_lfsr_bigint(s, tap, seq_len)
    lfsr_seqs.append(seqs)

t_precomp = time.time() - t0
print(f"Precomputation time: {t_precomp:.2f} s")
print(f"Total sequences: {n * (2**m - 1)}")

# ============================================================
# Step 9. Run statistical attacks
# ============================================================
# Strategy:
#   - Try analogs with k=1 first (cheapest: 1023 candidates each)
#   - Then k=2 if needed (1023^2 \approx 1M candidates)
#   - Brute-force any remaining registers
#
# For each analog g depending on registers i1, ..., ik:
#   g(t) = seq_{i1}(t) XOR seq_{i2}(t) XOR ... XOR c
#   As big integers: g_int = seq_i1 XOR seq_i2 XOR ... XOR c_mask
#   Hamming distance = popcount(g_int XOR gamma_int)
#
# We try TWO offset conventions (gamma might start from LFSR time 0 or 1):
#   offset=0: gamma[k] uses LFSR output at time k
#   offset=1: gamma[k] uses LFSR output at time k+1
# The one that gives a clear separation (best << T/2) is correct.

print("\n" + "=" * 60)
print("Step 9. Statistical attacks")
print("=" * 60)

recovered = {}  # register_index -> initial_state

# Build gamma as big integer (bit i = gamma[i])
gamma_int = 0
for i in range(N):
    if gamma[i]:
        gamma_int |= (1 << i)

# Sort analogs: smallest k first, then largest eps (best quality)
attack_order = sorted(
    [a for a in analogs if a['k'] > 0],
    key=lambda x: (x['k'], -x['eps'])
)

working_offset = None  # will be determined from first successful attack

for a in attack_order:
    var_indices = a['vars']

    # Skip if all registers in this analog are already recovered
    unrecovered = [j for j in var_indices if j not in recovered]
    if len(unrecovered) == 0:
        continue
    # For simplicity, only attack when ALL registers in this analog are unknown
    if len(unrecovered) != a['k']:
        continue

    k = a['k']
    T = a['T']
    c = a['c']
    eps = a['eps']

    print(f"\n--- Attacking with g = {a['g_str']} ---")
    print(f"    Registers: {[j+1 for j in var_indices]}, k={k}, T={T}, eps={eps:.4f}")

    c_mask_full = (1 << seq_len) - 1 if c == 1 else 0  # all-ones if c=1, else 0

    best_result = None  # will hold (dist, states_dict, offset_used)

    # Try both offsets (unless already determined)
    offsets_to_try = [working_offset] if working_offset is not None else [0, 1]

    for offset in offsets_to_try:
        T_mask = (1 << T) - 1
        gamma_seg = (gamma_int >> offset) & T_mask  # gamma bits [offset, offset+T)
        c_mask = (c_mask_full >> offset) & T_mask

        best_dist = T + 1
        best_states = None
        second_best = T + 1

        t0 = time.time()

        if k == 1:
            j = var_indices[0]
            for s in range(1, 2**m):
                seq = (lfsr_seqs[j][s] >> offset) & T_mask
                g_int = bxor(seq, c_mask)
                dist = popcount(bxor(g_int, gamma_seg))
                if dist < best_dist:
                    second_best = best_dist
                    best_dist = dist
                    best_states = {j: s}
                elif dist < second_best:
                    second_best = dist

        elif k == 2:
            j1, j2 = var_indices
            for s1 in range(1, 2**m):
                seq1 = (lfsr_seqs[j1][s1] >> offset) & T_mask
                for s2 in range(1, 2**m):
                    seq2 = (lfsr_seqs[j2][s2] >> offset) & T_mask
                    g_int = bxor(bxor(seq1, seq2), c_mask)
                    dist = popcount(bxor(g_int, gamma_seg))
                    if dist < best_dist:
                        second_best = best_dist
                        best_dist = dist
                        best_states = {j1: s1, j2: s2}
                    elif dist < second_best:
                        second_best = dist

        elif k == 3:
            j1, j2, j3 = var_indices
            for s1 in range(1, 2**m):
                seq1 = (lfsr_seqs[j1][s1] >> offset) & T_mask
                for s2 in range(1, 2**m):
                    seq12 = bxor(seq1, (lfsr_seqs[j2][s2] >> offset) & T_mask)
                    for s3 in range(1, 2**m):
                        seq3 = (lfsr_seqs[j3][s3] >> offset) & T_mask
                        g_int = bxor(bxor(seq12, seq3), c_mask)
                        dist = popcount(bxor(g_int, gamma_seg))
                        if dist < best_dist:
                            second_best = best_dist
                            best_dist = dist
                            best_states = {j1: s1, j2: s2, j3: s3}
                        elif dist < second_best:
                            second_best = dist

        t_attack = time.time() - t0

        expected_true = T * a['err']   # expected Hamming distance for correct state
        expected_false = T * 0.5       # expected for wrong state

        print(f"    [offset={offset}] time={t_attack:.2f}s  "
              f"best={best_dist}/{T} ({float(best_dist/T):.4f})  "
              f"2nd_best={second_best}/{T} ({float(second_best/T):.4f})  "
              f"expected_true~{expected_true:.0f}  expected_false~{expected_false:.0f}")

        # Accept if best is clearly better than second best and close to expected
        if best_dist < expected_false * 0.8:
            if best_result is None or best_dist < best_result[0]:
                best_result = (best_dist, best_states, offset)

    if best_result is not None:
        dist, states, off = best_result
        if working_offset is None:
            working_offset = off
            print(f"    => Determined working offset = {working_offset}")
        for j, s in states.items():
            recovered[j] = s
            print(f"    => R{j+1}: initial state = {s} (0b{s:010b})")
    else:
        print(f"    => Attack inconclusive (no clear winner)")

print(f"\n--- Recovered so far: {len(recovered)} / {n} registers ---")
for j in sorted(recovered.keys()):
    print(f"  R{j+1}: {recovered[j]} (0b{recovered[j]:010b})")

# ============================================================
# Step 10. Brute-force remaining registers
# ============================================================
# If some registers are still unknown, recover them by brute force:
# enumerate all candidate states, generate gamma using full f, compare.
# With recovered registers fixed, each candidate only varies the unknowns.

unrecovered_regs = [j for j in range(n) if j not in recovered]

if unrecovered_regs:
    print(f"\n{'='*60}")
    print(f"Step 10. Brute-forcing remaining registers: {[j+1 for j in unrecovered_regs]}")
    print(f"{'='*60}")

    num_unknown = len(unrecovered_regs)
    total_candidates = (2**m - 1) ** num_unknown
    print(f"Candidates to try: {total_candidates}")

    # Use first T_brute bits for verification
    T_brute = min(2000, N)
    off = working_offset if working_offset is not None else 1

    gamma_seg_brute = 0
    for i in range(T_brute):
        if gamma[i] == 1:
            gamma_seg_brute |= (1 << i)

    best_dist_bf = T_brute + 1
    best_states_bf = None
    t0 = time.time()

    # Build candidate iterator
    from itertools import product as iprod
    for combo in iprod(*(range(1, 2**m) for _ in unrecovered_regs)):
        # Set up all register states
        all_states = dict(recovered)
        for idx, j in enumerate(unrecovered_regs):
            all_states[j] = combo[idx]

        # Generate gamma using full f and truth table lookup
        gen_int = 0
        reg_states = [all_states[j] for j in range(n)]
        states_copy = list(reg_states)
        mask_m = (1 << m) - 1

        for t_idx in range(off + T_brute):
            if t_idx >= off:
                # Build input to f: bit j = LSB of register j's state
                f_input = 0
                for j in range(n):
                    if states_copy[j] & 1:
                        f_input |= (1 << j)
                if f_tt[f_input]:
                    gen_int |= (1 << (t_idx - off))

            # Clock all registers
            for j in range(n):
                fb = bin(states_copy[j] & tap_masks[j]).count('1') & 1
                states_copy[j] = ((states_copy[j] >> 1) | (fb << (m - 1))) & mask_m

        dist = popcount(bxor(gen_int, gamma_seg_brute))
        if dist < best_dist_bf:
            best_dist_bf = dist
            best_states_bf = dict(all_states)

    t_bf = time.time() - t0
    print(f"Brute-force time: {t_bf:.2f} s")
    print(f"Best Hamming distance: {best_dist_bf} / {T_brute}")

    if best_states_bf:
        for j in unrecovered_regs:
            recovered[j] = best_states_bf[j]
            print(f"  R{j+1}: {recovered[j]} (0b{recovered[j]:010b})")

print(f"\n*** All recovered initial states: ***")
for j in range(n):
    if j in recovered:
        print(f"  R{j+1}: {recovered[j]:4d} (0b{recovered[j]:010b})")
    else:
        print(f"  R{j+1}: UNKNOWN")

# ============================================================
# Step 11. Verification
# ============================================================
# Generate gamma from the recovered initial states using the full
# combining function f, and compare with the input gamma bit by bit.

print(f"\n{'='*60}")
print("Step 11. Verification")
print(f"{'='*60}")

if len(recovered) < n:
    print("Cannot verify: not all registers recovered.")
else:
    # Try both offset conventions for verification
    for off_name, off in [("offset=0", 0), ("offset=1", 1)]:
        reg_states = [recovered[j] for j in range(n)]
        states_v = list(reg_states)
        mask_m = (1 << m) - 1

        generated = []
        for t_idx in range(N + 2):
            # Evaluate f on current outputs (LSB of each register)
            f_input = 0
            for j in range(n):
                if states_v[j] & 1:
                    f_input |= (1 << j)
            generated.append(f_tt[f_input])

            # Clock all registers
            for j in range(n):
                fb = bin(states_v[j] & tap_masks[j]).count('1') & 1
                states_v[j] = ((states_v[j] >> 1) | (fb << (m - 1))) & mask_m

        # Compare with gamma (using offset)
        match_count = sum(1 for i in range(N) if gamma[i] == generated[i + off])
        pct = float(100 * match_count / N)
        print(f"  {off_name}: {match_count}/{N} ({pct:.2f}%)")

        if match_count == N:
            print(f"  *** VERIFICATION PASSED ({off_name})! ***")
            break

# ============================================================
# Summary
# ============================================================
print(f"\n{'='*60}")
print("Summary")
print(f"{'='*60}")
print(f"Combining function: f = {f_poly}")
print(f"deg(f) = {f_poly.degree()}, wt(f) = {wt_f}/{2**n}")
print(f"cor(f) = {cor_f}")
print(f"Number of affine analogs: {len(analogs)}")
print(f"Gamma length: {N}")
print(f"delta = {delta}")
print(f"LFSR length: m = {m}")
print(f"Number of LFSRs: n = {n}")
print(f"Precomputation time: {t_precomp:.2f} s")
print(f"\nRecovered initial states:")
for j in range(n):
    if j in recovered:
        print(f"  R{j+1}: {recovered[j]:4d} (0b{recovered[j]:010b})")

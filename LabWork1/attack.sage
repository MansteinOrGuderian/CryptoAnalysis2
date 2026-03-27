#!/usr/bin/env sage
# -*- coding: utf-8 -*-
"""
Algebraic attack on a filtering gamma generator
"""

import time
import re
import sys

# ============================================================
# Step 1. Define the Boolean polynomial ring B_64
# ============================================================
print("=" * 60)
print("Step 1. Initializing Boolean polynomial ring B_64")
print("=" * 60)

n = 64
B = BooleanPolynomialRing(n, 'x', order='degrevlex')
B.inject_variables()
x = B.gens()
print(f"Ring: {B}")
print(f"Number of variables: {n}")

# ============================================================
# Step 2. Feedback polynomial and companion matrix
# ============================================================
print("\n" + "=" * 60)
print("Step 2. Building the LFSR companion matrix")
print("=" * 60)

# p(x) = x^64 + x^63 + x^62 + x^60 + x^59 + x^58 + x^57 + x^56 +
#         x^53 + x^50 + x^47 + x^45 + x^44 + x^43 + x^42 + x^41 +
#         x^40 + x^39 + x^38 + x^37 + x^36 + x^34 + x^32 + x^30 +
#         x^28 + x^24 + x^18 + x^15 + x^14 + x^13 + x^11 + x^9 +
#         x^6 + x^4 + 1

R_poly = GF(2)['t']
t = R_poly.gen()
p = (t^64 + t^63 + t^62 + t^60 + t^59 + t^58 + t^57 + t^56 +
     t^53 + t^50 + t^47 + t^45 + t^44 + t^43 + t^42 + t^41 +
     t^40 + t^39 + t^38 + t^37 + t^36 + t^34 + t^32 + t^30 +
     t^28 + t^24 + t^18 + t^15 + t^14 + t^13 + t^11 + t^9 +
     t^6 + t^4 + 1)

print(f"Feedback polynomial: {p}")
print(f"Polynomial degree: {p.degree()}")

# Companion matrix over GF(2)
C_gf2 = companion_matrix(p, format='bottom')

# Feedback coefficients (c_0, c_1, ..., c_63)
# Recurrence: s_{j+64} = c_0*s_j + c_1*s_{j+1} + ... + c_63*s_{j+63}
feedback_coeffs = [int(p[i]) for i in range(64)]
feedback_indices = [i for i in range(64) if feedback_coeffs[i] == 1]
print(f"Nonzero feedback coefficient positions: {feedback_indices}")

# ============================================================
# Step 3. Read and parse the filter function
# ============================================================
print("\n" + "=" * 60)
print("Step 3. Filter function")
print("=" * 60)

with open('filterfuncV10.txt', 'r') as file:
    func_str = file.read().strip()

# Parse LaTeX notation: x_{i} -> x[i], space = multiplication, " + " = addition
terms = func_str.split(' + ')
f = B(0)
for term in terms:
    var_indices = re.findall(r'x_\{(\d+)\}', term)
    monomial = B(1)
    for v in var_indices:
        monomial *= x[int(v)]
    f += monomial

print(f"Degree of f: {f.degree()}")
print(f"Number of monomials in f: {len(f)}")

from collections import Counter
deg_distribution = Counter()
for monom in f:
    deg_distribution[monom.degree()] += 1
for d in sorted(deg_distribution.keys()):
    print(f"  Degree {d}: {deg_distribution[d]} monomials")

# ============================================================
# Step 4. Build Groebner bases for ideals <f> and <f+1>
# ============================================================
print("\n" + "=" * 60)
print("Step 4. Building Groebner bases")
print("=" * 60)

# --- Groebner basis of ideal <f> ---
print("\n--- Ideal <f> ---")
t0 = time.time()
I_f = B.ideal(f)
GB_f = I_f.groebner_basis()
t_gb_f = time.time() - t0
print(f"Build time: {t_gb_f:.4f} s")
print(f"Basis size: {len(GB_f)}")

degrees_f = [g.degree() for g in GB_f]
min_deg_f = min(degrees_f)
min_funcs_f = [g for g in GB_f if g.degree() == min_deg_f]
print(f"Minimal degree: {min_deg_f}")
print(f"Number of minimal-degree functions: {len(min_funcs_f)}")
print("Minimal-degree functions from <f>:")
for i, g in enumerate(min_funcs_f):
    print(f"  g_{i+1} = {g}")

# --- Groebner basis of ideal <f+1> ---
# For our variant, building <f+1> exceeds memory limits
# (polybori throws "Built-in matrix-size exceeded!").
# As stated in the task: if building one basis is too resource-intensive,
# we build only the other — provided it suffices for the attack.
# Since <f> yields a degree-2 function, this is more than sufficient.

print("\n--- Ideal <f + 1> ---")
print("SKIPPED: building <f+1> exceeds memory (matrix-size limit).")
print("The basis of <f> alone is sufficient for the attack,")
print("since it provides a function of degree 2.")
GB_f1 = None
t_gb_f1 = None
min_deg_f1 = float('inf')
min_funcs_f1 = []

# Algebraic immunity (upper bound)
# True AI = min(min_deg in <f>, min_deg in <f+1>)
# We only know the <f> side, so AI(f) <= min_deg_f
print(f"\n*** Algebraic immunity: AI(f) <= {min_deg_f} ***")
print(f"    (upper bound; <f+1> could not be computed)")

e = min_deg_f  # degree of equations we will use

# ============================================================
# Step 5. Read the gamma sequence
# ============================================================
print("\n" + "=" * 60)
print("Step 5. Reading gamma sequence")
print("=" * 60)

with open('gammaV10.txt', 'r') as file:
    gamma_str = file.read().strip()
gamma = [int(b) for b in gamma_str]
N = len(gamma)
n_zeros = gamma.count(0)
n_ones = gamma.count(1)
print(f"Gamma length: {N}")
print(f"Number of 0s: {n_zeros} ({float(100*n_zeros/N):.1f}%)")
print(f"Number of 1s: {n_ones} ({float(100*n_ones/N):.1f}%)")

# ============================================================
# Step 6. Build the system of lower-degree equations
# ============================================================
print("\n" + "=" * 60)
print("Step 6. Building the equation system")
print("=" * 60)

# Attack strategy:
#   We have h in <f> with deg h = 2.
#   <f> = Ann(f + 1), so h * (f + 1) = 0.
#   When gamma_i = 0  =>  f(C^i * x) = 0  =>  h(C^i * x) = 0.
#   We use ONLY the ticks where gamma_i = 0.

h = min_funcs_f[0]  # the degree-2 annihilator
print(f"Using annihilator: h = {h}")
print(f"Degree of h: {h.degree()}")

# Sanity check: h must satisfy h*(f+1) = 0 in B_64
print("\nSanity check: h*(f+1) == 0 ?")
check = h * (f + 1)
print(f"  h*(f+1) = {check}")
assert check == 0, "FATAL: annihilator check failed!"
print("  PASSED")

# Number of monomials up to degree e in n variables (linearization bound)
M = sum(binomial(n, k) for k in range(e + 1))
print(f"Number of monomials up to degree {e}: M = {M}")
print(f"Estimated equations needed: ~{M}")
print(f"Available zero-ticks in gamma: {n_zeros}")
print(f"  -> {'sufficient' if n_zeros >= M else 'INSUFFICIENT!'}")

# Generate the symbolic LFSR sequence via the recurrence relation.
# s[j] is a linear combination of x_0,...,x_{n-1} in B.
# The LFSR state at time i is: (s[i], s[i+1], ..., s[i+n-1]) = C^i * x
#
# Convention from the theory: gamma_i = f(C^i * x), i = 1,2,3,...
# In the file: gamma[0] = gamma_1, gamma[1] = gamma_2, etc.
# So gamma[k] corresponds to state (s[k+1], ..., s[k+64]).
#
# However, the exact convention is uncertain, so we will try BOTH:
#   Convention A: gamma[k] = f(C^{k+1} * x)  =>  state s[k+1 : k+65]
#   Convention B: gamma[k] = f(C^{k}   * x)  =>  state s[k   : k+64]

print("\nGenerating symbolic LFSR sequence...")
t0 = time.time()

max_ticks = N  # use ALL gamma bits

# s[0] = x_0, ..., s[63] = x_63 (initial state variables)
s = list(x)

# Extend via LFSR recurrence: s[j+64] = sum c_i * s[j+i]
for j in range(max_ticks + n + 1):
    new_bit = B(0)
    base = len(s) - n
    for idx in feedback_indices:
        new_bit += s[base + idx]
    s.append(new_bit)

t_seq = time.time() - t0
print(f"Sequence generation time: {t_seq:.4f} s")
print(f"Symbolic sequence length: {len(s)}")

# Build equations using Convention A (theory: gamma[k] <-> C^{k+1})
# If this doesn't verify, we will try Convention B.
for convention_name, offset in [("A (gamma[k] = f(C^{k+1}*x))", 1),
                                 ("B (gamma[k] = f(C^{k}*x))", 0)]:
    print(f"\n{'='*60}")
    print(f"Trying convention {convention_name}")
    print(f"{'='*60}")

    print("\nBuilding equations...")
    t0 = time.time()

    equations = []
    eq_count = 0

    for k in range(max_ticks):
        if gamma[k] != 0:
            continue  # only use ticks where gamma = 0

        i = k + offset  # state index
        state_i = s[i:i+n]
        eq = h(*state_i)

        if eq == 0:
            continue
        if eq == 1:
            print(f"  WARNING: contradiction at tick k={k}!")
            continue

        equations.append(eq)
        eq_count += 1

        if eq_count % 2000 == 0:
            print(f"  Built {eq_count} nontrivial equations (tick {k})")

    t_build = time.time() - t0
    print(f"\nEquation system build time: {t_build:.4f} s")
    print(f"Total nontrivial equations: {len(equations)}")
    print(f"Equation degree: {e}")

    # Print the first 10 equations
    print("\nFirst 10 equations:")
    for i_eq, eq in enumerate(equations[:10]):
        eq_str = str(eq)
        if len(eq_str) > 120:
            eq_str = eq_str[:120] + "..."
        print(f"  {i_eq+1}: {eq_str} = 0")

    # ============================================================
    # Step 7. Solve the equation system
    # ============================================================
    print(f"\nSolving a system of {len(equations)} equations...")
    t0 = time.time()
    sol_ideal = B.ideal(equations)
    sol_GB = sol_ideal.groebner_basis()
    t_solve = time.time() - t0
    print(f"Solve time: {t_solve:.4f} s")
    print(f"Solution basis size: {len(sol_GB)}")

    # Print the basis
    print("\nGroebner basis of the system:")
    for g in sol_GB:
        print(f"  {g}")

    # ============================================================
    # Extract solutions using variety()
    # ============================================================
    print("\nExtracting solutions via variety()...")
    t0 = time.time()
    try:
        V = sol_ideal.variety()
    except Exception as ex:
        print(f"variety() failed: {ex}")
        V = []
    t_var = time.time() - t0
    print(f"variety() time: {t_var:.4f} s")
    print(f"Number of solutions: {len(V)}")

    if len(V) == 0:
        print("No solutions found, trying next convention...")
        continue

    # Convert to bit vectors
    solutions = []
    for sol_dict in V:
        sol = [int(sol_dict[x[i]]) for i in range(n)]
        solutions.append(sol)
        sol_str = ''.join(map(str, sol))
        print(f"  Solution: {sol_str}  (Hex: {hex(int(sol_str, 2))})")

    # Filter nonzero (task requires nonzero initial state)
    nonzero_solutions = [s_sol for s_sol in solutions if any(b == 1 for b in s_sol)]
    print(f"Nonzero solutions: {len(nonzero_solutions)}")

    if len(nonzero_solutions) == 0:
        print("Only zero solution found, trying next convention...")
        continue

    # ============================================================
    # Step 8. Verification (try BOTH clock orderings)
    # ============================================================
    print(f"\n{'='*60}")
    print("Step 8. Verification")
    print(f"{'='*60}")

    verified_solution = None
    verified_mode = None

    for solution in nonzero_solutions:
        for mode_name, clock_first in [("clock-first (gamma[k]=f(C^{k+1}*x))", True),
                                        ("eval-first  (gamma[k]=f(C^{k}*x))", False)]:
            curr = vector(GF(2), solution)

            generated_gamma = []
            for i_tick in range(N):
                if clock_first:
                    curr = C_gf2 * curr
                    val = f(*list(curr))
                else:
                    val = f(*list(curr))
                    curr = C_gf2 * curr
                generated_gamma.append(int(val))

            match_count = sum(1 for i_m in range(N)
                              if gamma[i_m] == generated_gamma[i_m])
            pct = float(100 * match_count / N)
            print(f"  {mode_name}: {match_count}/{N} ({pct:.2f}%)")

            if match_count == N:
                print(f"  *** VERIFICATION PASSED with {mode_name}! ***")
                verified_solution = solution
                verified_mode = mode_name
                break

        if verified_solution is not None:
            break

    if verified_solution is not None:
        print(f"\n*** Initial state recovered correctly! ***")
        print(f"*** Verification mode: {verified_mode} ***")
        break  # exit the convention loop
    else:
        print(f"\nConvention {convention_name} did not produce a valid solution.")

# ============================================================
# Summary
# ============================================================
print("\n" + "=" * 60)
print("Summary")
print("=" * 60)
print(f"Groebner basis <f> build time:     {t_gb_f:.4f} s")
print(f"Groebner basis <f+1>:              SKIPPED (memory limit)")
print(f"Sequence generation time:          {t_seq:.4f} s")
print(f"Equation system build time:        {t_build:.4f} s")
print(f"System solve time:                 {t_solve:.4f} s")
print(f"Algebraic immunity:                AI(f) <= {e}")
print(f"Annihilator used:                  h = {h}")
print(f"Number of equations:               {len(equations)}")
if verified_solution is not None:
    vs = ''.join(map(str, verified_solution))
    print(f"Initial state:                     {vs}")
    print(f"Initial state (hex):               {hex(int(vs, 2))}")
    print(f"Verification mode:                 {verified_mode}")
else:
    print("Initial state:                     NOT FOUND")

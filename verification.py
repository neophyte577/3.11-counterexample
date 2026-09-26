"""exact finite checks for the accompanying preprint 
"""
from collections import Counter
from fractions import Fraction as Q
from itertools import combinations, product


def compose(g, f):
    return tuple(g[x] for x in f)


def fibers(f):
    return frozenset(frozenset(i for i, y in enumerate(f) if y == z)
                     for z in set(f))


def profile(f):
    return tuple(sorted(map(len, fibers(f))))


def transition(law):
    n = len(next(iter(law)))
    return tuple(tuple(sum((p for f, p in law.items() if f[i] == j), Q(0))
                       for j in range(n)) for i in range(n))


def matrix_product(a, b):
    return tuple(tuple(sum((a[i][k] * b[k][j] for k in range(len(b))), Q(0))
                       for j in range(len(b[0]))) for i in range(len(a)))


def stationary(pi, p):
    return tuple(sum((pi[i] * p[i][j] for i in range(len(pi))), Q(0))
                 for j in range(len(pi))) == pi


def reachable(p):
    n = len(p)
    for start in range(n):
        seen, todo = {start}, [start]
        while todo:
            for j, probability in enumerate(p[todo.pop()]):
                if probability and j not in seen:
                    seen.add(j)
                    todo.append(j)
        if len(seen) != n:
            return False
    return True


def check_graph_example():
    edges = {frozenset(e) for e in [(0, 1), (1, 2), (3, 4), (4, 5)]}
    oriented = [(u, v) for e in edges for u in e for v in e if u != v]
    colors = [{0, 2, 3, 5}, {0, 2, 4}]
    maps = {tuple(u if x in a else v for x in range(6))
            for a in colors for u, v in oriented}
    assert len(maps) == 16
    law = {f: Q(1, 16) for f in maps}
    pi = tuple(Q(x, 8) for x in [1, 2, 1, 1, 2, 1])
    p = transition(law)
    assert p == (pi,) * 6
    assert matrix_product(p, p) == p and stationary(pi, p)
    assert all(x > 0 for row in p for x in row)
    assert Counter(profile(f) for f in maps) == {(2, 4): 8, (3, 3): 8}
    for f, edge in product(maps, edges):
        assert frozenset(f[x] for x in edge) in edges
    for f, g in product(maps, repeat=2):
        gf = compose(g, f)
        assert gf in maps
        assert fibers(gf) == fibers(f)
        assert set(gf) == set(g)
    kernel_edges = {frozenset((i, j)) for i, j in combinations(range(6), 2)
                    if all(f[i] != f[j] for f in maps)}
    assert kernel_edges == edges
    for f in maps:
        assert all(sum((pi[x] for x in c), Q(0)) == Q(1, 2)
                   for c in fibers(f))
    print('Graph example: 16 maps, 64 edge checks, 256 compositions; all passed.')
    print('P rows:', ' '.join(map(str, pi)))
    print('Profile probabilities: (2, 4) = 1/2; (3, 3) = 1/2.')


def check_state_splitting():
    base = [(0, 1, 0, 1), (0, 1, 1, 0),
            (2, 3, 2, 3), (2, 3, 3, 2)]
    # Lifted ordering: 1a, 1b, 2, 3a, 3b, 4.
    q = (0, 0, 1, 2, 2, 3)
    sections = list(product((0, 1), (2,), (3, 4), (5,)))
    law = Counter()
    for f, s in product(base, sections):
        lift = tuple(s[f[q[y]]] for y in range(6))
        law[lift] += Q(1, len(base) * len(sections))
    assert len(law) == 8 and set(law.values()) == {Q(1, 8)}
    p = transition(law)
    base_p = transition({f: Q(1, 4) for f in base})
    multiplicity = [2, 1, 2, 1]
    expected = tuple(tuple(base_p[q[y]][q[z]] / multiplicity[q[z]]
                           for z in range(6)) for y in range(6))
    assert p == expected
    pi = tuple(Q(1, 4 * multiplicity[q[z]]) for z in range(6))
    assert stationary(pi, p) and reachable(p)
    assert any(p[i][i] > 0 for i in range(6))  # irreducibility -> period one
    assert any(x == 0 for row in p for x in row)
    profiles = Counter()
    for f, probability in law.items():
        profiles[profile(f)] += probability
        assert all(sum((pi[x] for x in c), Q(0)) == Q(1, 2)
                   for c in fibers(f))
    assert profiles == {(2, 4): Q(1, 2), (3, 3): Q(1, 2)}
    for g, f in product(law, repeat=2):
        assert set(compose(g, f)) == set(g)
        assert fibers(compose(g, f)) == fibers(f)
    print('State splitting: 8 maps; exact marginals, masses and 64 compositions passed.')


if __name__ == '__main__':
    check_graph_example()
    check_state_splitting()
    print('Fin.')

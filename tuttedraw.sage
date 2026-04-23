#!/usr/bin/env sage
# tuttedraw.sage — visualize planar and plane graphs
# (c) Manfred Scheucher 2017-2023
# Refactored with Claude based on draw.sage and draw_plane.sage

from ast import literal_eval
from copy import copy
import argparse


# ─── IPE export ──────────────────────────────────────────────────────────────

def graph_to_ipe(G, filepath, colormap):
    points = G.get_pos()

    x0 = min(x for x, y in points.values())
    y0 = min(y for x, y in points.values())
    x1 = max(1, max(x for x, y in points.values()) - x0)
    y1 = max(1, max(y for x, y in points.values()) - y0)
    M = 392 / max(x1, y1)
    points = {v: (100 + float((p[0] - x0) * M), 100 + float((p[1] - y0) * M))
              for v, p in points.items()}

    with open(filepath, 'w') as g:
        g.write('<?xml version="1.0"?>\n'
                '<!DOCTYPE ipe SYSTEM "ipe.dtd">\n'
                '<ipe version="70005" creator="Ipe 7.1.4">\n'
                '<info created="D:20150825115823" modified="D:20150825115852"/>\n')
        with open('ipestyle.txt') as f:
            for line in f:
                g.write('\t\t' + line)
        g.write('<page>\n'
                '<layer name="alpha"/>\n'
                '<layer name="beta"/>\n'
                '<view layers="alpha beta" active="alpha"/>\n')

        for i, j, color in G.edges():
            layer = 'alpha' if color is not None else 'beta'
            col = 'black' if color not in colormap else colormap[color]
            x, y = points[i]
            u, v = points[j]
            g.write(f'<path stroke="{col}" pen="heavier" layer="{layer}">\n'
                    f'{x} {y} m\n{u} {v} l\n</path>\n')

        for x, y in points.values():
            g.write(f'<use layer="alpha" name="mark/disk(sx)" '
                    f'pos="{x} {y}" size="normal" stroke="black"/>\n')

        g.write('</page>\n</ipe>')

    print(f"finished {filepath}")


# ─── Tutte embedding ──────────────────────────────────────────────────────────

def tutte_layout(G, outer_face, weights):
    V = G.vertices()
    pos = {}
    l = len(outer_face)
    a0 = pi / l + pi / 2
    for i in range(l):
        ai = a0 + pi * 2 * i / l
        pos[outer_face[i]] = (cos(ai), sin(ai))

    n = len(V)
    M = zero_matrix(RR, n, n)
    b = zero_matrix(RR, n, 2)
    for i, v in enumerate(V):
        if v in pos:
            M[i, i] = 1
            b[i, 0], b[i, 1] = pos[v]
        else:
            s = 0
            for u in G.neighbors(v):
                j = V.index(u)
                wu = weights[u, v]
                s += wu
                M[i, j] = -wu
            M[i, i] = s

    sol = M.pseudoinverse() * b
    return {V[i]: tuple(sol[i]) for i in range(n)}


def best_outer_face(G):
    """Pick the outer face with highest automorphism symmetry."""
    drawn = set()
    candidates = []
    maxsym = 0

    for face in G.faces():
        outer_face = [e[0] for e in face]
        augmented = Graph(G)
        for v in outer_face:
            augmented.add_edge((-1 - v, v))

        sym = augmented.automorphism_group().order()
        if sym < maxsym:
            continue
        if sym > maxsym:
            maxsym = sym
            candidates = []

        gstr = augmented.canonical_label().sparse6_string()
        if gstr in drawn:
            continue
        drawn.add(gstr)
        candidates.append(outer_face)

    print(f"maxsym: {maxsym}")
    outer_face = candidates[0]
    for C in candidates:
        if len(C) >= 20:
            outer_face = C
    return outer_face


def draw_tutte(args):
    for ct, line in enumerate(open(args.filename), 1):
        line = line.strip()
        print(f"graph# {ct} : {line}")

        if args.informat == 'list':
            G = Graph(literal_eval(line))
        else:
            G = Graph(line)  # s6 / g6

        G.set_pos(G.layout_planar())
        outer_face = best_outer_face(G)

        weights = {(u, v): 1 for u, v in G.edges(labels=False)}
        weights.update({(v, u): 1 for u, v in G.edges(labels=False)})
        G.set_pos(tutte_layout(G, outer_face, weights))

        colormap = {c: c for c in G.edge_labels()}
        outfile = f"{args.filename}_{ct}.tutte.{args.outformat}"

        if args.outformat in ('pdf', 'png'):
            G.plot().save(outfile)
        elif args.outformat == 'ipe':
            graph_to_ipe(G, outfile, colormap)

        print(f"wrote visualization to {outfile}\n")


# ─── Plane graph (plantri) ────────────────────────────────────────────────────

def cyclic_rotations(L):
    for i in range(len(L)):
        yield L[i:] + L[:i]


def add_edge(edges, u, v):
    if (u, v) > (v, u):
        u, v = v, u
    assert (u, v) not in edges
    edges.append((u, v))


def triangulate(f, E_plus):
    for v in f:
        if f.count(v) > 1:
            iv = f.index(v)
            f = f[iv:] + f[:iv]
            add_edge(E_plus, f[-1], f[1])
            triangulate(f[1:], E_plus)
            return

    E_induced = {e for e in E_plus if len(set(e) & set(f)) == 2}
    u = [v for v in f if len({e for e in E_induced if v in e}) == 2][0]
    iu = f.index(u)
    f = f[iu:] + f[:iu]
    for v in f[2:-1]:
        add_edge(E_plus, u, v)


def draw_plane(args):
    for ct, line in enumerate(open(args.filename), 1):
        line = line.strip()
        print(f"plane graph #{ct} : {line}")

        n_str, rot_str = line.split(' ', 1)
        n = int(n_str)
        rotations = rot_str.split(',')
        V = [chr(ord('a') + i) for i in range(n)]
        rot = {V[i]: rotations[i] for i in range(n)}

        E = [(v, w) for v in V for w in rot[v] if v < w]

        F = []
        for v in V:
            for w in rot[v]:
                face = [w]
                prv, cur = v, w
                while True:
                    nxt = rot[cur][rot[cur].index(prv) - 1]
                    prv, cur = cur, nxt
                    if prv == v and cur == w:
                        break
                    face.append(cur)
                if face == min(cyclic_rotations(face)):
                    F.append(face)

        assert len(V) - len(E) + len(F) == 2, "Euler formula check failed"

        E_plus = copy(E)
        for f in F:
            triangulate(f, E_plus)
        assert len(E_plus) == 3 * len(V) - 6

        G_plus = Graph(E_plus)
        G_plus.is_planar(set_pos=True)
        G = Graph(E)
        G.set_pos(G_plus.get_pos())

        outfile = f"{args.filename}_{ct}.{args.outformat}"
        G.plot().save(outfile)
        print(f"wrote to {outfile}\n")


# ─── CLI ──────────────────────────────────────────────────────────────────────

def parse_args():
    parser = argparse.ArgumentParser(
        description="Visualize planar and plane graphs using Tutte embeddings"
    )
    sub = parser.add_subparsers(dest='mode', required=True)

    p_tutte = sub.add_parser('tutte', help='Tutte embedding of planar graphs')
    p_tutte.add_argument('filename', help='Input file (one graph per line)')
    p_tutte.add_argument('informat', choices=['list', 's6', 'g6'],
                         help='list = Python edge list, s6 = sparse6, g6 = graph6')
    p_tutte.add_argument('outformat', choices=['png', 'pdf', 'ipe'],
                         help='Output format')

    p_plane = sub.add_parser('plane', help='Draw plane graphs (plantri output)')
    p_plane.add_argument('filename', help='Input file (plantri -a -p output)')
    p_plane.add_argument('outformat', choices=['png', 'pdf'], default='png',
                         nargs='?', help='Output format (default: png)')

    return parser.parse_args()


args = parse_args()
if args.mode == 'tutte':
    draw_tutte(args)
elif args.mode == 'plane':
    draw_plane(args)

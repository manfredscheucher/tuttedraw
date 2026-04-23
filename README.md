# tuttedraw

SageMath scripts to visualize planar and plane graphs using Tutte embeddings.

---

## Files

- `tuttedraw.sage` — combined, refactored script (refactored with Claude based on `old/draw.sage` and `old/draw_plane.sage`)
- `old/` — original scripts

## Usage

### `tutte` — planar graphs (edge list / sparse6 / graph6)

```
sage tuttedraw.sage tutte <file> <informat> <outformat>
```

Given a plain text file where each line encodes a planar graph, this computes a nice visualization for each graph. The script iteratively computes weighted Tutte embeddings and sets appropriate weights so that faces and edges become "nice". For more information see http://arxiv.org/abs/1708.06449.

| Argument | Values | Description |
|---|---|---|
| `informat` | `list`, `s6`, `g6` | `list` = Python edge list, `s6` = sparse6, `g6` = graph6 |
| `outformat` | `png`, `pdf`, `ipe` | `ipe` = XML for the [Ipe extensible drawing editor](https://ipe.otfried.org/), which can be used to further modify the visualization in a WYSIWYG style |

### `plane` — plane graphs (plantri output)

```
sage tuttedraw.sage plane <file> [outformat]
```

When enumerating plane graphs using plantri (e.g. `plantri 5 -a -c1 -p > plane5.ascii`), each line encodes a plane graph via cyclic rotations around each vertex. The script adds auxiliary edges to triangulate the plane graph so that the embedding is unique (Whitney's theorem) and then computes vertex coordinates. See: https://users.cecs.anu.edu.au/~bdm/plantri/

---

## Example

```
$ sage tuttedraw.sage tutte example_graph_edgelist.txt list png
graph# 1 : [(0, 17), (0, 18), ...]
maxsym: 4
wrote visualization to example_graph_edgelist.txt_1.tutte.png
```

## Dependencies

- [SageMath](https://www.sagemath.org/)
- `ipestyle.txt` (included, required for `ipe` output)

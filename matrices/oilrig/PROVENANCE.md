# Oil-rig matrix pair

These are the **complete original numerical matrix files**, not a truncated submatrix or a
synthetic replacement. The stiffness matrix is `bcsstk02.mtx`; the matching mass matrix is
`bcsstm02.mtx`. They are the 66 x 66 statically condensed oil-rig generalized eigenproblem
`K x = lambda M x` in the Harwell-Boeing / SuiteSparse collections. K has 2211 stored
lower-triangle entries (4356 expanded nonzeros); M has 66 diagonal entries.

Primary collection descriptions:

- https://math.nist.gov/MatrixMarket/data/Harwell-Boeing/bcsstruc1/bcsstk02.html
- https://math.nist.gov/MatrixMarket/data/Harwell-Boeing/bcsstruc1/bcsstm02.html

NIST identifies the source as John Lewis, Boeing Computer Services, and the discipline as
dynamic analyses in structural engineering. The original headers identify collection
editors I. Duff, R. Grimes and J. Lewis, and the 1982 accession.

The UTF-8 Matrix Market copies included here were obtained from this pinned public mirror:

- Repository: `x3042/ExactODEReduction.jl`
- Commit: `156926246c8621ae20460066769c1345fe6e7682`
- K path: `data/sscollection/bcsstk02/bcsstk02.mtx`
- M path: `data/sscollection/bcsstm02/bcsstm02.mtx`
- K Git blob SHA-1: `0831271feb8da04c7792b34985efb00d68a21192`
- M Git blob SHA-1: `aafd612dbe95f77c90e3624de925003400801ccf`

The included files were checked **byte for byte through their Git blob hashes**. Their
original headers are retained. `tools/check_matrices.py` repeats that check offline;
`SHA256SUMS` records additional hashes of the bundled files.

This source attribution is separate from the MIT license of the new SecularLab code.
No assertion is made that a mirror's software license changes the original collection's
terms. Retain the collection acknowledgments when using or redistributing these data.
The benchmark values are published numerical data, not a current oil-rig measurement
campaign. High precision certifies the represented rational pencil, not the physical model.

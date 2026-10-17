# Inputs and outputs

## Matrix input

The header is `%%MatrixMarket matrix coordinate real symmetric`, or one of the supported
combinations `coordinate/array`, `real/integer`, `symmetric/general`. Matrices must be square,
order 1..256. Complex, pattern, Hermitian and skew-symmetric fields are rejected explicitly.

Coordinate entries are 1-based `row column value`. Symmetric storage must use only one triangle,
plus diagonal entries. Duplicate coordinates are summed exactly; mixing upper and lower
triangles is rejected rather than double-counted. General input supplies all intended entries.
Array storage is column-major, using the lower triangle for symmetric input and **one numeric
value per data line**. Lines beginning with `%` and empty lines are ignored after the header.
Trailing unexpected data, incomplete matrices and invalid dimensions are rejected.

The numeric parser accepts integers, fractions, finite decimal and scientific notation,
including Fortran D exponents. Fraction numerators/denominators must themselves be integers.
A signed denominator is normalized. Whitespace is allowed around a token, not inside it.
Each numeric token is limited to 20000 characters; decimal exponents to magnitude 10000.

Fractions in `.mtx` files are a **documented SecularLab extension**, not part of the standard
Matrix Market real format. Some synthetic examples and the generated canonical input copies
use these exact fractions. Conventional Matrix Market readers may reject them. The original
external oil-rig files are standard real-coordinate Matrix Market files.

## `solve` output directory

- `input_K.mtx`, `input_M.mtx`: canonical exact matrices. With affine input, K is the effective
  specialized matrix; these two files alone are sufficient for replay.
- `command.txt`: invocation arguments and the exact specialized parameter, when present.
- `coefficients.tsv`: descending exact integer coefficients of det(zLM-LK), including its
  leading coefficient, with L identified in the comment.
- `result.json`: status, counts, limits, backend and phase timings.
- `eigenvalues.tsv`: written only on complete certification. Includes index, fixed-place
  decimal, exact rational lower/upper endpoints, and exact radius about that decimal.

TSV headers identify the format version and decimal-place count. Indices are zero-based,
largest eigenvalue first. A radius `p/q` is exact, not a decimal display rounded in an
uncontrolled direction. Singleton intervals are legitimate exact rational eigenvalues.
A result JSON alone is not a certificate. Match it with the matrices and TSV.

The program prints progress to stderr and accepted root values to stdout. It creates the
output directory only after parsing/validating ordinary options. Domain errors discovered
while analyzing the matrices can leave input/command files in that new directory, but
cannot produce a false successful certificate file. Preserve successful run directories
and choose a new `--out` directory when rerunning.

## Coefficient-only output

A fixed pencil produces `coefficients.tsv`. A formal affine pencil produces
`coefficients_affine.tsv`, whose columns are `z_power`, `t_power`, `integer_coefficient`.
It represents det(zLM-LK0-tLK1), not an already normalized monic polynomial. All coefficients
retain actual signs. The known determinant multiplier is L^n.

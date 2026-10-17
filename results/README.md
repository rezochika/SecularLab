# Local runs

Local run configurations write into this directory. Each run needs a new or empty subdirectory.
Recorded delivery tests are in `../verification/`; never replace those with a differently built run
without changing their provenance records. To compare platforms, send back the new `result.json`,
`eigenvalues.tsv`, and the compiler/GMP/Boost versions.

Delphy
======

Delphy is a fast, scalable, accurate and accessible tool for Bayesian phylogenetics based on Explicit Mutation-Annotated
Trees (EMATs).  EMATs are an extension of the Mutation-Annotated Trees (MATs) introduced by
[UShER](https://github.com/yatisht/usher) where nodes have explicit times, mutations are represented as explicit timed
events along branches, and missing data is explicitly represented.  EMATs open substantial simplifications and scaling
opportunities in the calculations powering Bayesian phylogenetics, at the cost of some statistical efficiency with
respect to the traditional calculations based on Felsenstein pruning.  For genomic epidemiology datasets, where the
total number of mutations on a tree is comparable to the number of samples, this is a very favorable trade-off.

These sources comprise the "core" computational engine of Delphy.  The web application that allows users to immediately
and intuitively use Delphy (currently at [https://delphy.fathom.info](https://delphy.fathom.info)) is developed in
collaboration with [Fathom Information Design](https://fathom.info), and its separately licensed sources are hosted
[here](https://github.com/fathominfo/delphy-web).

References
----------

_TO COME_: Whitepaper with overview of key ideas and accuracy+speed benchmarks (expected: late May 2024)

_TO COME_: Preprint with full details and benchmarking (expected: late June 2024)

Credits and Acknowledgements
----------------------------

Delphy draws a lot of inspiration from:

- [UShER and matOptimize](https://github.com/yatisht/usher)
- [BEAST](https://github.com/beast-dev/beast-mcmc)
- [BEAST2](https://github.com/CompEvol/beast2)
- [MAPLE](https://github.com/NicolaDM/MAPLE)
- [IQ-TREE2](github.com/iqtree/iqtree2/)

Copyright (c) 2022-2024 Broad Institute, Inc.  See [LICENSE](LICENSE) for details.

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

* [Whitepaper with overview of key ideas and accuracy+speed benchmarks](delphywp.pdf)

* _TO COME_: Preprint with full details and benchmarking (expected: late July 2024)

System Requirements
-------------------
Delphy can be compiled either as a native, standalone command-line program (`delphy` and its graphical cousin, `delphy_ui`), or as a WebAssembly bundle that lies at the core of [`delphy-web`](https://github.com/fathominfo/delphy-web).  Build instructions for both can be found in [`INSTALL.md`](INSTALL.md).  Delphy does not require special hardware, such as a GPU.

Delphy was developed and is primarily tested under Linux (Ubuntu 22.04.4 LTS, x86-64).  [Delphy-web](https://github.com/fathominfo/delphy-web), deployed at [https://delphy.fathom.info](https://delphy.fathom.info), should work in any modern web browser; see the [Delphy-web README.md](https://github.com/fathominfo/delphy-web) for specific browsers and operating systems on which it has been tested.

Google Colab Tutorials
-------------------
Ready-to-use Google Colab tutorials for downloading and formatting sequencing data from the [NCBI Virus database](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/), as well as user-provided data, can be accessed [here](https://colab.research.google.com/github/broadinstitute/delphy/blob/main/tutorials/delphy_workflow.ipynb) (blank notebook) and [here](https://colab.research.google.com/github/broadinstitute/delphy/blob/main/tutorials/delphy_workflow_sars_example.ipynb) (SARS-CoV-2 example).  These notebooks streamline the uniform formatting of sequencing data and the associated metadata for input into Delphy, and then run Delphy using a pre-compiled binary.  The output files can be visualized on the [Delphy web interface](https://delphy.fathom.info).  

Credits and Acknowledgements
----------------------------

Delphy is developed in the [Sabeti Lab](https://www.sabetilab.org/) at the [Broad
Institute](https://www.broadinstitute.org/).

Delphy draws a lot of inspiration from:

- [UShER and matOptimize](https://github.com/yatisht/usher)
- [BEAST](https://github.com/beast-dev/beast-mcmc)
- [BEAST2](https://github.com/CompEvol/beast2)
- [MAPLE](https://github.com/NicolaDM/MAPLE)
- [IQ-TREE2](github.com/iqtree/iqtree2/)

Copyright (c) 2022-2024 Broad Institute, Inc.  See [LICENSE](LICENSE) for details.

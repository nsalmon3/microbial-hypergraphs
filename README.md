# microbial-hypergraphs
A command line tool for generating statistics and visualizations of microbial hypergraphs

**WARNING** Running hypercorrelations on populations with more than 100 OTUs, group sizes larger than 4, or thresholds too low (under 0.3 for pearson hypercorrelations), can result in egregious runtimes or memory usage. The computational complexity of these computations makes it effectively impossible to investigate the full hypergraph (all group sizes) of any dataset. Please use care when running the package.
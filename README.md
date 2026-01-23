# microbial-hypergraphs
A command line tool for generating statistics and visualizations of microbial hypergraphs

**WARNING** Running hypercorrelations on populations with more than 100 OTUs, group sizes larger than 4, or thresholds too low (under 0.3 for pearson hypercorrelations), can result in egregious runtimes or memory usage. The computational complexity of these computations makes it effectively impossible to investigate the full hypergraph (all group sizes) of any dataset. Please use care when running the package.

## Using the CLI

For now, try a command like:
`python -m microbial_hypergraphs hypercorrelate --population-name=PresentInAllSamples --correlation-name=Pearson --hypercorrelation-name=Minimal --group-size=3 --threshold=0.5`

This should put a `.xlsx` file in the `exports` directory with hypercorrelation information.
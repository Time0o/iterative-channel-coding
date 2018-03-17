A quantitative comparison of different iterative decoding algorithms applied to
classical channel code descriptions.

# Prerequisites
* make
* a C++ compiler supporting C++11
* Python 3.*

# Running demos
Run `make demo` to see the decoding schemes in action.

# Running benchmarks
Run `make benchmark` to run the decoding performance benchmarks for the code
descriptions given in `resources/control_matrices.txt`. Results are placed
per code in appropriately named subfolders of the `benchmarks` directory. They
can be visualized by invoking the `tool/plot` script with one of the subfolders
as the first command line argument.

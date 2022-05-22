# SBWT

This repository is under construction.

# Compiling

```
git submodule init
git submodule update
cd build
cmake ..
make
```

**Troubleshooting**: If you run into problems involving the &lt;filesystem&gt; header, you probably need to update your compiler. The compiler `g++-8` should be sufficient. Install a new compiler and direct CMake to use it with the `-DCMAKE_CXX_COMPILER` option. For example, to set the compiler to `g++-8`, run CMake with the option `-DCMAKE_CXX_COMPILER=g++-8`. 

# Index construction

For small inputs, we have an in-memory construction algorithm that loads all k-mers in memory. To build the plain matrix variant for the example data in example_data/coli3.fna, you can use the following command:

```
./build/bin/build_plain_matrixboss --in-fasta example_data/coli3.fna -o index.matrixboss -k 30
```

For larger inputs, we provide contruction from a Themisto index. Themisto is included as a submodule in this repository. First, you need to install Themisto by going to its subdirectory `./Themisto` and following the compilation instruction in the readme of Themisto (sorry about not having automatic compilation of Themisto included in this repostory). After compiling Themisto, to build the Themisto index on our example data, run the following:

```
./Themisto/build/bin/themisto build -k 30 -i example_data/coli3.fna --temp-dir temp --no-colors -o example_data/coli3
```

This will write the index into the file example_data/coli3.tdbg. You can then build the plain matrix SBWT with:

```
./build/bin/build_plain_matrixboss --in-themisto example_data/coli3.tdbg -o example_data/coli3.matrix
```

You can use the plain matrix representation to build any of our variants. For example, to build the rrr-compressed matrix variant, run:

```
./build/bin/build_variant_from_matrix -i example_data/coli3.matrix --variant rrr-matrix --temp-dir temp -o example_data/coli3.rrrmatrix
```

Running `./build/bin/build_variant_from_matrix` without parameters gives the list of all available variants.

# Running queries

Currently the query program only accepts the plain matrix variant. To change the variant, you need to edit the template parameter in the source code and recompile. This is pretty bad -- we are working on making the k-mer search function take the variant as a parameter. Anyway, the queries on the plain matrix variant can be run as follows:

```
./build/bin/kmer-search -i example_data/coli3.matrix -q example_data/queries.fna -o out.txt -k 30 --temp-dir temp
```

This prints for each query of length n in the input a line containing n-k+1 space-separated integers, which are the ranks of the columns representing the k-mer in the index. If the k-mer is not found, -1 is printed.

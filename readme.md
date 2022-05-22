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

To build one of the SBWT variants, run the executable `sbwt_build`.

For small inputs, we have an in-memory construction algorithm that loads all k-mers in memory. This is run if the `--in-fasta` option is given, like in the example below.

```
./build/bin/sbwt_build --in-fasta example_data/coli3.fna -o index.sbwt -k 30 --variant plain-matrix
```

For larger inputs, we provide contruction from a Themisto index. We are also working on direct construction from a KMC database. Themisto is included as a submodule in this repository. First, you need to install Themisto by going to its subdirectory `./Themisto` and following the compilation instruction in the readme of Themisto. After compiling Themisto, to build the Themisto index on our example data, run the following:

```
./Themisto/build/bin/themisto build -k 30 -i example_data/coli3.fna --temp-dir temp --no-colors -o example_data/coli3
```

This will write the index into the file example_data/coli3.tdbg. You can then build the plain matrix SBWT with:

```
./build/bin/sbwt_build --in-themisto example_data/coli3.tdbg -o index.sbwt -k 30 --variant plain-matrix
```

The list of all command line options and parameters is below:

```
Construct an SBWT variant.
Usage:
  ./build/bin/sbwt_build [OPTION...]

  -o, --out-file arg     Output filename.
      --variant arg      The SBWT variant to build. Available variants: 
                         plain-matrix rrr-matrix mef-matrix plain-split 
                         rrr-split mef-split plain-concat mef-concat 
                         plain-subsetwt rrr-subsetwt
      --in-fasta arg     Build in internal memory from a FASTA file (takes 
                         a lot of memory). (default: "")
      --in-themisto arg  Build from a Themisto .tdbg file. (default: "")
  -k arg                 Value of k (must not be given if --in-themisto is 
                         given because themisto defines the k) (default: 0)
  -h, --help             Print usage
```

# Running queries

Currently the query program only accepts the plain matrix variant. To change the variant, you need to edit the template parameter in the source code and recompile. This is pretty bad -- we are working on making the k-mer search function take the variant as a parameter. Anyway, the queries on the plain matrix variant can be run as follows:

```
./build/bin/kmer-search -i example_data/coli3.matrix -q example_data/queries.fna -o out.txt -k 30 --temp-dir temp
```

This prints for each query of length n in the input a line containing n-k+1 space-separated integers, which are the ranks of the columns representing the k-mer in the index. If the k-mer is not found, -1 is printed.

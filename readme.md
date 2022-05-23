# SBWT

This is the code for the paper [Succinct k-mer Set Representations Using Subset Rank Queries on the Spectral Burrows-Wheeler Transform (SBWT)](https://www.biorxiv.org/content/10.1101/2022.05.19.492613v1). The repository includes implementations of the various SBWT variants described in the paper. Note that contrary to many other k-mer membership data structures, our code is not aware of DNA reverse complements. That is, it considers a k-mer and its reverse complement as separate k-mers.

We are currently actively working on the code. Top items on the to-do list are the following:

* Streaming queries, that is, queries with a rolling k-mer window. Currently, the code queries all k-mers separately without reusing computation from previous k-mers.
* Reverse complement aware indexing.
* Construction directly from a sorted KMC database.

# Compiling

```
git submodule init
git submodule update
cd build
cmake .. -DMAX_KMER_LENGTH=32
make
```

Change the parameter `-DMAX_KMER_LENGTH=32` to increase the maximum allowed k-mer length, up to 255.

**Troubleshooting**: If you run into problems involving the `<filesystem>` header, you probably need to update your compiler. The compiler `g++-8` should be sufficient. Install a new compiler and direct CMake to use it with the `-DCMAKE_CXX_COMPILER` option. For example, to set the compiler to `g++-8`, run CMake with the option `-DCMAKE_CXX_COMPILER=g++-8`. 

# Index construction

To build one of the SBWT variants, run `./build/bin/sbwt build`.

For small inputs, we have an in-memory construction algorithm that loads all k-mers in memory. This is run if the `--in-fasta` option is given, like in the example below.

```
./build/bin/sbwt build --in-fasta example_data/coli3.fna -o index.sbwt -k 30 --variant plain-matrix
```

For larger inputs, we provide contruction from a Themisto index. Themisto is included as a submodule in this repository. First, you need to install Themisto by going to its subdirectory `./Themisto` and following the compilation instruction in the readme of Themisto. After compiling Themisto, to build the Themisto index on our example data, run the following:

```
./Themisto/build/bin/themisto build -k 30 -i example_data/coli3.fna --temp-dir temp --no-colors -o example_data/coli3
```

This will write the index into the file example_data/coli3.tdbg. You can then build the plain matrix SBWT with:

```
./build/bin/sbwt build --in-themisto example_data/coli3.tdbg -o index.sbwt -k 30 --variant plain-matrix
```

The list of all command line options and parameters is below:

```
Construct an SBWT variant.
Usage:
  ./build/bin/sbwt build [OPTION...]

  -o, --out-file arg     Output filename.
      --variant arg      The SBWT variant to build. Available variants: 
                         plain-matrix rrr-matrix mef-matrix plain-split 
                         rrr-split mef-split plain-concat mef-concat 
                         plain-subsetwt rrr-subsetwt
      --streaming-support  Build the auxiliary bit vector for streaming 
                           query support.
      --in-fasta arg     Build in internal memory from a FASTA file (takes 
                         a lot of memory). (default: "")
      --in-themisto arg  Build from a Themisto .tdbg file. (default: "")
  -k arg                 Value of k (must not be given if --in-themisto is 
                         given because themisto defines the k) (default: 0)
  -h, --help             Print usage
```

# Running queries

To query for existence of all k-mers in an index for all sequences in a fasta-file, run the following command:

```
./build/bin/sbwt search -i temp/index.sbwt -q example_data/queries.fna -o out.txt
```

This prints for each query of length n in the input a line containing n-k+1 space-separated integers, which are the ranks of the columns representing the k-mer in the index. If the k-mer is not found, -1 is printed. The full options are:

```
Query all k-mers of all input reads.
Usage:
  ./build/bin/sbwt search [OPTION...]

  -o, --out-file arg    Output filename.
  -i, --index-file arg  Index input file.
  -q, --query-file arg  The query in FASTA format.
  -h, --help            Print usage
```

# For developers: building and running the tests 

```
git submodule init
git submodule update

# Build googletest
cd googletest
mkdir build
cd build
cmake ..
make

# Build SBWT
cd ../../build
cmake .. -DCMAKE_BUILD_TYPE=Debug -DMAX_KMER_LENGTH=32
make
```

This will build the executable `./build/bin/sbwt_tests`. Make sure to run the test executable from the root of the repository, or otherwise it will not find the example data in ./example_data.
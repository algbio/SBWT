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

The index is currently constructed from the index of the Themisto tool. Themisto is included as a submodule in this repository. First, you need to install Themisto by going to its subdirectory `./Themisto` and following the compilation instruction in the readme of Themisto (sorry about not having automatic compilation of Themisto included in this repostory). After compiling Themisto, to build the Themisto index on our example data, run the following:

```
./Themisto/build/bin/themisto build -k 30 -i example_data/coli3.fna --temp-dir temp --no-colors -o example_data/coli3
```

This will write the index into the file example_data/coli3.tdbg.

# Running queries

Example:

```
./build/bin/kmer-search -i coli3682_concat.fasta.tdbg -q queries.fna --temp-dir /tmp -o answers.txt
```


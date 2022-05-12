# Building

```
git submodule init
git submodule update
cd build
cmake ..
make
```

**Troubleshooting**: If you run into problems involving the &lt;filesystem&gt; header, you probably need to update your compiler. The compiler `g++-8` should be sufficient. Install a new compiler and direct CMake to use it with the `-DCMAKE_CXX_COMPILER` option. For example, to set the compiler to `g++-8`, run CMake with the option `-DCMAKE_CXX_COMPILER=g++-8`. 

# Running queries

Example:

```
./build/bin/kmer-search -i coli3682_concat.fasta.tdbg -q queries.fna --temp-dir /tmp -o answers.txt
```


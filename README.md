# SBWT

This is the code for the paper [Succinct k-mer Set Representations Using Subset Rank Queries on the Spectral Burrows-Wheeler Transform (SBWT)](https://www.biorxiv.org/content/10.1101/2022.05.19.492613v1). The repository includes implementations of the various SBWT variants described in the paper. The data structures answer k-mer membership queries on the input data. Note that contrary to many other k-mer membership data structures, our code is not aware of DNA reverse complements. That is, it considers a k-mer and its reverse complement as separate k-mers.

This construction algorithm is based on the lightning-fast [k-mer counter KMC](https://github.com/refresh-bio/KMC). We call the KMC binaries directly from our code. The construction is very disk-heavy, so it is recommended to run construction code off a fast SSD drive.

# Compiling

The following commands have been tested to successfully build SBWT on a clean Ubuntu 18.04 Docker image.

```
apt-get update
apt-get install -y g++ gcc cmake git python3-dev g++-8 libz-dev
git clone https://github.com/algbio/SBWT
cd SBWT/build
cmake .. -DCMAKE_CXX_COMPILER=g++-8 -DMAX_KMER_LENGTH=32
make -j8
```

Change the parameter `-DMAX_KMER_LENGTH=32` to increase the maximum allowed k-mer length, up to 255. Larger values lead to slower construction and higher disk usage during construction.

**Troubleshooting**: If you run into problems involving the `<filesystem>` header, you probably need to update your compiler. The compiler `g++-8` should be sufficient. Install a new compiler and direct CMake to use it with the `-DCMAKE_CXX_COMPILER` option. For example, to set the compiler to `g++-8`, run CMake with the option `-DCMAKE_CXX_COMPILER=g++-8`.

Note: the Elias-Fano variants make use of the `_pext_u64` instruction in the BMI2 instruction set. Older CPUs might not support this instruction. In that case, we fall back to a simple software implementation, which will ruin the performance of the Elias-Fano variants (those whose variant name starts with "mef").

# Index construction

Below is the command to build the SBWT for input data `example_data/coli3.fna` provided in this repository, with k = 30. The index is written to the file `index.sbwt`.

```
./build/bin/sbwt build -i example_data/coli3.fna -o index.sbwt -k 30
```

This builds the default variant, which is the plain matrix SBWT. Other variant can be specified with the `--variant` option.
The list of all command line options and parameters is below:

```
Construct an SBWT variant.
Usage:
  build [OPTION...]

  -i, --in-file arg             The input sequences as a FASTA or FASTQ
				file, possibly gzipped. If the file
				extension is .txt, the file is interpreted
				as a list of input files, one file on each
				line. All input files must be in the same
				format.
  -o, --out-file arg            Output file for the constructed index.
  -k, --kmer-length arg         The k-mer length.
      --variant arg             The SBWT variant to build. Available
				variants: plain-matrix rrr-matrix
				mef-matrix plain-split rrr-split mef-split
				plain-concat mef-concat plain-subsetwt
				rrr-subsetwt (default: plain-matrix)
      --add-reverse-complements
				Also add the reverse complement of every
				k-mer to the index. Warning: this creates a
				temporary reverse-complemented duplicate of
				each input file before construction. Make
				sure that the directory at --temp-dir can
				handle this amount of data. If the input is
				gzipped, the duplicate will also be
				compressed, which might take a while.
      --no-streaming-support    Save space by not building the streaming
				query support bit vector. This leads to
				slower queries.
  -t, --n-threads arg           Number of parallel threads. (default: 1)
  -a, --min-abundance arg       Discard all k-mers occurring fewer than
				this many times. By default we keep all
				k-mers. Note that we consider a k-mer
				distinct from its reverse complement.
				(default: 1)
  -b, --max-abundance arg       Discard all k-mers occurring more than this
				many times. (default: 1000000000)
  -m, --ram-gigas arg           RAM budget in gigabytes (not strictly
				enforced). Must be at least 2. (default: 2)
  -d, --temp-dir arg            Location for temporary files. (default: .)
  -h, --help                    Print usage
```

# Running queries

To query for existence of all k-mers in an index for all sequences in a fastq-file, run the following command:

```
./build/bin/sbwt search -i index.sbwt -q example_data/queries.fastq -o out.txt
```

This prints for each query of length n in the input a line containing n-k+1 space-separated integers, which are the ranks of the columns representing the k-mer in the index. If the k-mer is not found, -1 is printed. If the index was built streaming support (which is the default), the faster streaming query algorithm is automatically used. The full options are:

```
Query all k-mers of all input reads.
Usage:
  search [OPTION...]

  -o, --out-file arg    Output filename.
  -i, --index-file arg  Index input file.
  -q, --query-file arg  The query in FASTA or FASTQ format, possibly
			gzipped. Multi-line FASTQ is not supported. If the
			file extension is .txt, this is interpreted as a
			list of query files, one per line. In this case,
			--out-file is also interpreted as a list of output
			files in the same manner, one line for each input
			file.
  -z, --gzip-output     Writes output in gzipped form. This can shrink the
			output files by an order of magnitude.
  -h, --help            Print usage
```

# API

The API for the SBWT is still in the works. Do not expect a stable API at this point.

The SBWT can be constructed and queried using the [SBWT class](https://htmlpreview.github.io/?https://github.com/algbio/SBWT/blob/master/doc/html/classsbwt_1_1SBWT.html). The class is templatized by the underlying subset rank support structure. See [here](https://htmlpreview.github.io/?https://github.com/algbio/SBWT/blob/master/doc/html/variants_8hh_source.html) for types of subset rank query data structures are suitable for the template parameter.

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
cmake .. -DCMAKE_BUILD_TYPE=Debug -DBUILD_TESTS=1 -DMAX_KMER_LENGTH=32
make
```

This will build the executable `./build/bin/sbwt_tests`. Make sure to run the test executable from the root of the repository, or otherwise it will not find the example data in ./example\_data.

# Acknowledgements

The command-line parsing and the gzip support are implemented using [cxxopts](https://github.com/jarro2783/cxxopts) and [zstr](https://github.com/mateidavid/zstr) respectively.

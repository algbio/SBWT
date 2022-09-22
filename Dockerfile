FROM ubuntu:18.04

# Set some timezone or otherwise tzdata hangs the build.
ENV TZ=Asia/Dubai
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt-get update && apt-get install -y g++ gcc cmake git python3-dev g++-8 libz-dev libbz2-dev

RUN git clone https://github.com/algbio/SBWT
WORKDIR /SBWT
RUN git checkout dev

WORKDIR /SBWT/build
RUN cmake .. -DCMAKE_CXX_COMPILER=g++-8 -DMAX_KMER_LENGTH=32 -DBUILD_TESTS=1
RUN make -j8
run /SBWT/build/bin/sbwt

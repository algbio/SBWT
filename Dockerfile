FROM ubuntu:22.04

# Set some timezone or otherwise tzdata hangs the build.
ENV TZ=Asia/Dubai
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt-get update && apt-get install -y g++ gcc cmake git python3-dev g++-10 libz-dev libbz2-dev

RUN git clone https://github.com/algbio/SBWT
WORKDIR /SBWT
RUN git checkout master

WORKDIR /SBWT/build
RUN cmake .. -DCMAKE_CXX_COMPILER=g++-10 -DMAX_KMER_LENGTH=32 -DBUILD_TESTS=1
RUN make -j8
run /SBWT/build/bin/sbwt

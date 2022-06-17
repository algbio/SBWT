FROM ubuntu:18.04

# Set some timezone or otherwise tzdata hangs the build.
ENV TZ=Asia/Dubai
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt-get update && apt-get install -y g++ gcc cmake libbz2-dev git python3-dev

RUN git clone https://github.com/algbio/SBWT
WORKDIR /SBWT
RUN git checkout kmcmake

WORKDIR /SBWT/build
RUN cmake ..
RUN make -j8
run /SBWT/build/bin/sbwt

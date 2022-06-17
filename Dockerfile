FROM ubuntu 

RUN apt-get update && apt-get install -y g++ gcc cmake libbz2-dev git python3-dev

RUN git clone https://github.com/algbio/SBWT
WORKDIR /SBWT
RUN git checkout kmcmake
#RUN git submodule init
#RUN git submodule update
#WORKDIR /SBWT/KMC
#RUN make -j8

WORKDIR /SBWT/build
RUN cmake ..
RUN make -j8
run /SBWT/build/bin/sbwt

# Yeah, this include and linking process is a bit messy currently, sorry about that.

SDSL_INCLUDES=-I ../build/external/sdsl-lite/build/include
DIVSUFSORT_INCLUDES=-I ../build/external/sdsl-lite/build/external/libdivsufsort/include
SEQIO_INCLUDES=-I ../build/external/SeqIO/include
SBWT_INCLUDES=-I ../include -I ../include/sbwt 

ALL_INCLUDES=${SDSL_INCLUDES} ${DIVSUFSORT_INCLUDES} ${SEQIO_INCLUDES} ${SBWT_INCLUDES}

SBWT_LIBS=-L $(shell pwd)/../build/external/sdsl-lite/build/lib/
KMC_OBJECTS=../build/external/KMC/build/libkmc_core.a ../build/external/KMC/build/libkmc_tools.a

api_example:
	${CXX} -g -std=c++17 api_example.cpp ../build/libsbwt_static.a ${KMC_OBJECTS} ${ALL_INCLUDES} ${SBWT_LIBS} -lsdsl -lz -o api_example -Wno-deprecated-declarations -pthread




CFLAGS=-Ib
CXXFLAGS=-Ib -std=c++17 -O2

ALL=generator_B judge_B generator_A judge_A main
VPATH=b

all: ${ALL}
clean: 
	rm -f ${ALL}

CC = g++
CFLAGS = -Wall -Werror -O2 -std=c++11
LIBS = -lm -lglpk

LIBOPT = libopt.a
AR = ar
ARFLAGS = rv

SRCS = main.cpp lp_opt.cpp helper.cpp
OBJS =  $(addsuffix .o, $(basename $(SRCS)))

all: opt_fault_tol
	
%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

opt_fault_tol: $(SRCS) $(LIBOPT)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

$(LIBOPT) : $(OBJS)
	$(AR) $(ARFLAGS) $@ $+

.PHONY: clean

clean:
	rm -rf *.o opt_fault_tol $(LIBOPT)

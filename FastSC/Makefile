Dist: wscp.h main.cpp mersenne_rng.cpp
	g++ \
	-fdiagnostics-color=always \
	-I/opt/ibm/ILOG/CPLEX_Studio2211/cplex/include \
	-I/opt/ibm/ILOG/CPLEX_Studio2211/concert/include \
	-DIL_STD \
	main.cpp mersenne_rng.cpp \
	-L/opt/ibm/ILOG/CPLEX_Studio2211/cplex/lib/x86-64_linux/static_pic \
	-L/opt/ibm/ILOG/CPLEX_Studio2211/concert/lib/x86-64_linux/static_pic \
	-lconcert -lilocplex -lcplex -lm -lpthread \
	-O3 -o FastSC
	
	rm -f *~

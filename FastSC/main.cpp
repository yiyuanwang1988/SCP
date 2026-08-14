#include <fstream>
#include "wscp.h"

int main(int argc, char *argv[])
{
    ifstream infile(argv[1]);
    int  var_num, set_num;
    infile >> var_num >> set_num;
    int time_limit = atoi(argv[2]);
    WSCP wscp_solver(time_limit);
    int seed = atoi(argv[3]);
	int AC = atoi(argv[4]);
	int use_cplex = atoi(argv[5]);
    
	
    
	wscp_solver.seed = seed;
    wscp_solver.build_instance(argv[1],AC);
    wscp_solver.reduce_instance();
    wscp_solver.set_param(use_cplex);
	
    start_timing();
    wscp_solver.init();
    wscp_solver.local_search();
    wscp_solver.check_solu();
	int cnt = 0;
    for (int i = 0; i < wscp_solver.set_num; ++i)
    {
		if (wscp_solver.best_solu[i] == 1) cnt++;	
	}
    cout << "instance: " << argv[1] << " cost: " << wscp_solver.best_cost << " seed: " << argv[3] << " time: " << wscp_solver.best_time << " set number: " << cnt <<endl;
    wscp_solver.output_solu();
	wscp_solver.free_memory();
}

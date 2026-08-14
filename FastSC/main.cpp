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
    // srand(seed);
	// rng.seed(seed);
	// srandom(seed);
	
    int new_weight, tabu_len;
    double novelty_p;
	// bool flag1 = 0,flag2 = 0,flag3 = 0, flag4 = 0, flag5 = 0, flag6 = 0;
	// if(set_num >= 500 || var_num >= 500) flag1 = 1;
	// if(var_num < set_num) flag2 = 1;
	// if(set_num > 100000 || var_num > 100000) flag3 = 1;
	// if(var_num > 2500 && var_num < 5000 && set_num > 920000 && set_num < 1100000) flag4 = 1;
	// if(/*(var_num > 100 && var_num < 150 && set_num > 150 && set_num < 200) ||
 	// 	(var_num > 450 && var_num < 500 && set_num > 2500 && set_num < 2600) ||
    //     (var_num > 450 && var_num < 500 && set_num > 2300 && set_num < 2400) ||*/
    //     (var_num > 600 && var_num < 650 && set_num > 3450 && set_num < 3500) ||
    //     (var_num > 2800 && var_num < 2850 && set_num > 13400 && set_num < 13500) ||
	// 	(var_num > 19100 && var_num < 19200 && set_num > 62500 && set_num < 62600) ||
    //     (var_num > 74000 && var_num < 75000 && set_num > 166000 && set_num < 167000))
    //      	flag5 = 1;
	// if(var_num > 120000 && var_num < 810000 && set_num > 620000 && set_num < 3900000) flag6 = 1;
//	cout << "var_num " << var_num << " set_num " << set_num << " flag4 " << flag4 << " flag5 " << flag5 << " flag6 " << flag6 << endl;
	// if(flag3)
	// {
	//     if (argc == 7){
	//         new_weight = atoi(argv[4]);
	//         tabu_len = atoi(argv[5]);
	//         novelty_p = atof(argv[6]);
	//     }
	//     else{
	//         if (set_num < var_num){
	//             new_weight = 5;
	//             tabu_len = 4;
	//             novelty_p = 0.1;
	//         }else{
	//             new_weight = 80;
	//             tabu_len = 5;
	//             if(var_num / set_num > 10)
	//                 novelty_p = 0.05;
	//             else novelty_p = 0.5;
	//         }
	//     }
	// }
	// else
	// {
		// if (argc == 7){
	    //     new_weight = atoi(argv[4]);
	    //     tabu_len = atoi(argv[5]);
	    //     novelty_p = atof(argv[6]);
	    // }
	    // else{
	    //     if (set_num > var_num){
	    //         new_weight = 5;
	    //         tabu_len = 4;
	    //         novelty_p = 0.1;
	    //     }else{
	    //         new_weight = 80;
	    //         tabu_len = 5;
	    //         if(var_num / set_num > 10)
	    //             novelty_p = 0.05;
	    //         else novelty_p = 0.5;
	    //     }
	    // }
	//}

	// wscp_solver.flag1 = flag1;
	// wscp_solver.flag2 = flag2;
	// wscp_solver.flag3 = flag3;
	// wscp_solver.flag4 = flag4;
	// wscp_solver.flag5 = flag5;
	// wscp_solver.flag6 = flag6;
	wscp_solver.seed = seed;
    wscp_solver.build_instance(argv[1],AC);
    wscp_solver.reduce_instance();
    wscp_solver.set_param(use_cplex);
	// printf("argc %d set_num %d var_num %d\n", argc, set_num, var_num);
	// printf("new_weight %d tabu_len %d novelty_p %f \n", new_weight, tabu_len, novelty_p);
    
    start_timing();
    wscp_solver.init();
    wscp_solver.local_search();
    wscp_solver.check_solu();
	int cnt = 0;
    //cout << "v ";
    for (int i = 0; i < wscp_solver.set_num; ++i)
    {
    	//cout << wscp_solver.best_solu[i] << " ";
		if (wscp_solver.best_solu[i] == 1) cnt++;	
	}
       
    //cout << endl;
    cout << "instance: " << argv[1] << " cost: " << wscp_solver.best_cost << " seed: " << argv[3] << " time: " << wscp_solver.best_time << " set number: " << cnt <<endl;
    wscp_solver.output_solu();
	wscp_solver.free_memory();
}

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include <limits.h>
#include <float.h>
#include <memory.h>
#include <sys/times.h>
#include <unistd.h>
#include <math.h>
#include <set>
#include <algorithm>
#include <utility>
#include <ilcplex/ilocplex.h>
#include "mersenne_rng.h"

using namespace std;

#define mypop(stack) stack[--stack##_fill_pointer]
#define mypush(item, stack) stack[stack##_fill_pointer++] = item

static struct tms start_time;
static double get_runtime()
{
    struct tms stop;
    times(&stop);
    return (double)(stop.tms_utime - start_time.tms_utime + stop.tms_stime - start_time.tms_stime) / sysconf(_SC_CLK_TCK);
}
static void start_timing()
{
    times(&start_time);
}

class WSCP
{
  public:
    //set(variable, soft clause) information
    long long *score;
    long long *pscore;
    long long *time_stamp;
    int *cost;
    int *org_cost;
    bool *cc;
    int *fix;
    set<int> tabu_list;
    int *zero_stack;
    int zero_stack_fill_pointer;
    int *index_in_zero_stack;
	int *var_sort; 
    //variable(hard clause)
    long long *weight;
    int *cover_count;
    int *avg_cost; 
    int *cover_set;
    int *cover_set2;
    int *uncover_stack;
    int *index_in_uncover_stack;
    int *unlocked_stack;
    int uncover_stack_fill_pointer;
    int unlocked_stack_fill_pointer;
    //solution information
    int *cur_solu;
    long long cur_cost;
    int reduce_cost;
    int *best_solu;
    long long best_cost;
    long long best_count;
    double best_time;
    int *set_solu;
    int *index_in_set_solu;
    int set_count;
    int gap;

    //date structure
    int cutoff_time;
    unsigned int seed;
    int bms_t;
    long long step;
    long long max_step;
    long long tries;
    long long max_tries;
	bool flag_tabu;
    int var_num;
    int set_num;

    int *var_delete;
    int *set_delete;
    int **var_set;
    int *var_set_num;
    int *org_var_set_num;
    int **set_var;
    int *set_var_num;
    int *org_set_var_num;

    int *goodset_stack;
    int goodset_stack_fill_pointer;

    bool is_update = 0;

    WSCP();
    WSCP(int time_limit);
    void all_memory();
    void free_memory();

    void build_instance(char *file_name, int num);
    void reduce_instance();
    void init();
    void local_search();
    void local_search1();
    bool cplex_improve();
    bool cplex_improve(const int *start_solu, long long start_cost);
    void add_set(int s);
    void remove_set(int s);
    void flip(int flip_set);
    void cover(int selected_var);
    void uncover(int selected_var);
    int break_ties_score_subscore(int set1, int set2, int op);
    void compute_score(int compute_set, long long & set_score, int & set_weight);
    int compare(int set1, int set2, int flag);
    int compare1(int s1, int c1, int s2, int c2);
    int compare2(int s1, int c1, int s2, int c2, int p1, int p2);
    int compareV2(int s1, int c1, int s2, int c2, int p1, int p2, int weight1, int weight2);
    int compareV3(int s1, int c1, int s2, int c2, int p1, int p2, int weight1, int weight2);
    int (WSCP::*compareVersion)(int s1, int c1, int s2, int c2, int p1, int p2, int weight1, int weight2);
    bool cmp(int x,int y);
    double weight1 = 4;
    double weight2 = 1;
    double novelty_p = 0.5;
    double explore_p = 0.5;
    int tabu_length = 5;
    int select_set(int is_tabu);
    int select_set_1();
    int select_set_new(int is_tabu);
    int select_set_Novelty();
    int select_set_Novelty_no_tabu(int var1);
    int select_set_Novelty_best();
    int local_greedy_addition();
    int remove_best_score();
    int select_set_from_zero_stack();
    int select_set_from_zero_stack_new();
    void set_param(int flag);
    void update_weight();
    void check_solu();
    void check_cur_solu();
    int count_2cc(int var);
    void output_solu();
    vector<int> index_in_unlocked_stack;
    bool flag1,flag2,flag3,flag4,flag5,flag6; 
    bool ratioType1,ratioType2,ratioType3,ratioType4,ratioType5; 
    int flag_cplex;
};

Mersenne_rng rng;

WSCP::WSCP()
{
    step = 0;
    bms_t = 25;
    cutoff_time = 1000;
    max_step = INT_MAX;
    max_tries = INT_MAX;
}

WSCP::WSCP(int time_limit)
{
    step = 0;
    bms_t = 25;
    cutoff_time = time_limit;
    max_step = INT_MAX;
    max_tries = INT_MAX;
}

void WSCP::set_param(int flag)
{


    ratioType1 = false,ratioType2 = false,ratioType3 = false, ratioType4 = false, ratioType5 = false;
    double ratio = var_num / (set_num*1.0);
    if(ratio < 0.01) {
        ratioType1 = true;
    }
    else if((ratio >= 0.1 && ratio < 0.21) || (ratio >= 0.21 && ratio < 1.0)) {
        if (std::fabs(ratio - 0.2083) < 1e-4) ratioType2 = true;
        else ratioType3 = true;
    }
    else if (ratio >= 1.0 && ratio < 10.0) {
        ratioType4 = true;
    }
    else if (ratio >= 10.0) {
        ratioType5 = true;
    }
    
    if(ratioType1) {
        weight1 = 3;
		tabu_length = 94; 
		novelty_p = 0.21;
        bms_t = 5000;
    }
    if(ratioType2) {
        weight1 = 80;
		tabu_length = 5;
		novelty_p = 0.5;
        bms_t = 600; 
    }
	if(ratioType3) {
		weight1 = 15;
		tabu_length = 10;
		novelty_p = 0.21;
        bms_t = 1000;
	}
    if (ratioType4) {
        weight1 = 90; 
		tabu_length = 10; 
		novelty_p = 0.5;
        bms_t = 170;
    }
    if (ratioType5) {
        weight1 = 80;
		tabu_length = 5;
		novelty_p = 0.05;
        bms_t = 170;
    }
    flag_cplex = flag;
}

bool WSCP::cmp(int x,int y)
{
	return var_set_num[y] * avg_cost[x] > var_set_num[x] * avg_cost[y];
}

int WSCP::compare1(int s1, int c1, int s2, int c2)
{
    if (c1 == c2)
    {
        if (s1 > s2)
            return 1;
        else if (s1 == s2)
            return 0;
        else
            return -1;
    }
    long long t1 = s1, t2 = s2;
    t1 = t1 * c2;
    t2 = t2 * c1;
    if (t1 > t2)
        return 1;
    else if (t1 == t2)
        return 0;
    else
        return -1;
}

int WSCP::compare2(int s1, int c1, int s2, int c2, int p1, int p2)
{
    if (c1 == c2)
    {
        if ((s1) > (s2) || (s1 == s2 && p1 > p2))
            return 1;
        else if ((s1 == s2) && (p1 == p2))
            return 0;
        else
            return -1;
    }
    long long t1 = s1, t2 = s2;
    t1 = t1 * c2;
    t2 = t2 * c1;
    if (t1 > t2)
        return 1;
    else if (t1 == t2)
        return 0;
    else
        return -1;
}

int WSCP::compareV2(int s1, int c1, int s2, int c2, int p1, int p2, int weight1, int weight2)
{
    int combined_s1 = weight1*s1 + weight2*p1;
    int combined_s2 = weight1*s2 + weight2*p2;

    // fix combined score
    if (c1 == c2)
    {
        if (combined_s1 > combined_s2)
            return 1;
        else if (combined_s1 == combined_s2)
            return 0;
        else
            return -1;
    }
    long long t1 = s1, t2 = s2;
    t1 = t1 * c2;
    t2 = t2 * c1;
    if (t1 > t2)
        return 1;
    else if (t1 == t2)
        return 0;
    else
        return -1;
}

int WSCP::compareV3(int s1, int c1, int s2, int c2, int p1, int p2, int weight1, int weight2)
{
    int combined_s1 = weight1*s1 + weight2*p1;
    int combined_s2 = weight1*s2 + weight2*p2;

    int comp_score, comp_comb_score;
    long long t1 = combined_s1, t2 = combined_s2;
    int gap1, gap2;

    t1 = t1 * c2;
    t2 = t2 * c1;
    gap1 = abs(t1 - t2);
    if (t1 > t2)
        comp_comb_score = 1;
    else if (t1 == t2)
        comp_comb_score = 0;
    else
        comp_comb_score = -1;
    //return comp_comb_score;
    t1 = p1;
    t2 = p2;
    t1 = t1 * c2;
    t2 = t2 * c1;
    gap2 = (weight1 + weight2) * abs(t1 - t2);
    if (t1 > t2)
        comp_score = 1;
    else if (t1 == t2)
        comp_score = 0;
    else
        comp_score = -1;

    if (comp_score == comp_comb_score)
    {
        return comp_comb_score;
    }
    else if (comp_score * comp_comb_score == 0)
    {
        return comp_score + comp_comb_score;
    }
    else
    {   
        return comp_comb_score;
    }
}

void WSCP::compute_score(int compute_set, long long & set_score, int & set_weight) {
    //int mu = 3;
    set_score = score[compute_set]*weight1 + pscore[compute_set]; // (*compute_score_ptr)(compute_set);
    set_weight = cost[compute_set];
}

int WSCP::break_ties_score_subscore(int set1, int set2, int op) {

    long long mult1 = score[set1] * cost[set2];
    long long mult2 = score[set2] * cost[set1];
    if(mult1 == mult2) {
        if(pscore[set1] < pscore[set2]) return -1;
        if(pscore[set1] > pscore[set2]) return 1;
        return 0;
    }
    if(mult1 > mult2) return 1;
    if(mult1 < mult2) return -1;
    return 0;
}

int WSCP::compare(int set1, int set2, int flag) {
    long long score1, score2;
    int weight1, weight2;
    compute_score(set1, score1, weight1);
    compute_score(set2, score2, weight2);
    long long mult1 = score1*weight2;
    long long mult2 = score2*weight1;
    if(mult1 == mult2) {
        int re = break_ties_score_subscore(set1, set2, flag);
        if(re == 0)
            return (rand_r(&seed) % 3) - 1; 
        return re;
    }
    else if(mult1 > mult2) return 1;
    else return -1;
}


void WSCP::all_memory()
{
    int max_set_num = set_num + 10;
    int max_var_num = var_num + 10;

    score = new long long[max_set_num];
    pscore = new long long[max_set_num];
    time_stamp = new long long[max_set_num];
    cost = new int[max_set_num];
    org_cost = new int[max_set_num];
    cc = new bool[max_set_num];
    fix = new int[max_set_num];
    zero_stack = new int[max_set_num];
    index_in_zero_stack = new int[max_set_num];

    //variable(hard clause)
    weight = new long long[max_var_num];
    cover_count = new int[max_var_num];
    cover_set = new int[max_var_num];
    cover_set2 = new int[max_var_num];
    uncover_stack = new int[max_var_num];
    index_in_uncover_stack = new int[max_var_num];
	var_sort = new int[max_var_num]; 
	avg_cost = new int[max_var_num];
    //solution information
    cur_solu = new int[max_set_num];
    best_solu = new int[max_set_num];
    set_solu = new int[max_set_num];
    index_in_set_solu = new int[max_set_num];

    //date structure
    var_delete = new int[max_var_num];
    var_set = new int *[max_var_num];
    var_set_num = new int[max_var_num];
    org_var_set_num = new int[max_var_num];
    set_var = new int *[max_set_num];
    set_var_num = new int[max_set_num];
    org_set_var_num = new int[max_set_num];
    unlocked_stack = new int[max_set_num];
    goodset_stack = new int[max_set_num];
}

void WSCP::free_memory()
{
    delete score;
    delete pscore;
    delete time_stamp;
    delete cost;
    delete org_cost;
    delete cc;
    delete fix;
    delete zero_stack;
    delete index_in_zero_stack;

    delete weight;
    delete cover_count;
    delete cover_set;
    delete cover_set2;
    delete uncover_stack;
    delete index_in_uncover_stack;

    delete cur_solu;
    delete best_solu;
    delete set_solu;
    delete index_in_set_solu;

    for (int i = 0; i < var_num; ++i)
        delete var_set[i];

    delete var_delete;
    delete var_set;
    delete var_set_num;
    delete org_var_set_num;
    delete org_set_var_num;

    for (int i = 0; i < set_num; ++i)
        delete set_var[i];
    delete set_var;
    delete set_var_num;

    delete goodset_stack;
}

void WSCP::build_instance(char *file_name, int num)
{
    ifstream infile(file_name);

    infile >> var_num >> set_num;
    if(set_num > var_num)
    {
        if(var_num < 1000)
        {
            bms_t = 80;
            compareVersion = &WSCP::compareV3;
        }
        else
        {
            bms_t = 170;
            compareVersion = &WSCP::compareV2;
        }
    }
    else
    {
        bms_t = 170;
        compareVersion = &WSCP::compareV3;
    }
    all_memory();
	gap = 0;
    for (int i = 0; i < set_num; ++i)
    {
        infile >> cost[i];
	    if(num) cost[i] = ((i + 1) % 200) + 1;
        org_cost[i] = cost[i];
        gap = std::max(cost[i],gap);
    }
    memset(set_var_num, 0, set_num * sizeof(int));
    memset(fix, 0, set_num * sizeof(bool));
    for (int i = 0; i < var_num; ++i)
    {
        infile >> var_set_num[i];
        avg_cost[i] = 0;
        var_set[i] = new int[var_set_num[i] + 1]; 
        for (int j = 0; j < var_set_num[i]; ++j)
        {
            infile >> var_set[i][j];
            if(!num) --var_set[i][j];
            avg_cost[i] += cost[var_set[i][j]];
            ++set_var_num[var_set[i][j]];
        }
        avg_cost[i] /= var_set_num[i];
    }

    for (int i = 0; i < set_num; ++i)
    {
        set_var[i] = new int[set_var_num[i] + 1];
        set_var_num[i] = 0; //reset for build set_var
    }
    //build set_var
    for (int i = 0; i < var_num; ++i)
    {
        for (int j = 0; j < var_set_num[i]; ++j)
        {
            int cur_set = var_set[i][j];
            set_var[cur_set][set_var_num[cur_set]++] = i;
        }
    }
}

void WSCP::reduce_instance()
{
    int s1, s2, s3;
    reduce_cost = 0;
    memset(fix, 0, set_num * sizeof(int));
    for (int i = 0; i < var_num; ++i)
    {
        if (var_set_num[i] == 1)
        {
            fix[var_set[i][0]] = 1;
        }
        else if (var_set_num[i] == 2)
        {
            s1 = var_set[i][0];
            s2 = var_set[i][1];
            if (fix[s1] != 0 || fix[s2] != 0)
                continue;
            if (set_var_num[s1] < set_var_num[s2])
            {
                int flag1 = 0;
                for (int j = 0; j < set_var_num[s1]; ++j)
                {
                    int v = set_var[s1][j];
                    int flag = 0;
                    for (int k = 0; k < var_set_num[v]; ++k)
                    {
                        if (var_set[v][k] == s2)
                        {
                            flag = 1;
                            break;
                        }
                    }
                    if (flag == 0)
                    {
                        flag1 = 1;
                        break;
                    }
                }
                if (flag1 == 0)
                {
                    if (cost[s1] >= cost[s2])
                    {
                        fix[s1] = -1;
                        fix[s2] = 1;
                    }
                    else if (cost[s1] < cost[s2])
                    {
                        fix[s1] = s2 + 10;
                        cost[s2] -= cost[s1];
                        reduce_cost += cost[s1];
                    }
                }
            }
            else
            {
                int flag1 = 0;
                for (int j = 0; j < set_var_num[s2]; ++j)
                {
                    int v = set_var[s2][j];
                    int flag = 0;
                    for (int k = 0; k < var_set_num[v]; ++k)
                    {
                        if (var_set[v][k] == s1)
                        {
                            flag = 1;
                            break;
                        }
                    }
                    if (flag == 0)
                    {
                        flag1 = 1;
                        break;
                    }
                }
                if (flag1 == 0)
                {
                    if (set_var_num[s1] == set_var_num[s2])
                    {
                        if (cost[s1] <= cost[s2])
                        {
                            fix[s1] = 1;
                            fix[s2] = -1;
                        }
                    }
                    else if (cost[s2] >= cost[s1])
                    {
                        fix[s2] = -1;
                        fix[s1] = 1;
                    }
                    else if (cost[s2] < cost[s1])
                    {
                        fix[s2] = s1 + 10;
                        cost[s1] -= cost[s2];
                        reduce_cost += cost[s2];
                    }
                }
            }
        }
        else if (var_set_num[i] == 3)
        {
            s1 = var_set[i][0];
            s2 = var_set[i][1];
            s3 = var_set[i][2];
            if (fix[s1] != 0 || fix[s2] != 0 || fix[s3] != 0)
                continue;
            if (set_var_num[s1] == set_var_num[s2] && set_var_num[s1] < set_var_num[s3])
            {
                int flag1 = 0;
                for (int j = 0; j < set_var_num[s1]; ++j)
                {
                    int v = set_var[s1][j];
                    int flag2 = 0;
                    int flag3 = 0;
                    for (int k = 0; k < var_set_num[v]; ++k)
                    {
                        if (var_set[v][k] == s2)
                        {
                            flag2 = 1;
                        }
                        else if (var_set[v][k] == s3)
                        {
                            flag3 = 1;
                        }
                        if (flag2 == 1 && flag3 == 1)
                            break;
                    }
                    if (flag2 != 1 || flag3 != 1)
                    {
                        flag1 = 1;
                        break;
                    }
                }
                if (flag1 == 0)
                {
                    if (fix[s3] == 1)
                    {
                        fix[s1] = -1;
                        fix[s2] = -1;
                    }
                    else if (cost[s1] >= cost[s3] && cost[s2] >= cost[s3])
                    {
                        fix[s1] = -1;
                        fix[s2] = -1;
                        fix[s3] = 1;
                    }
                    else
                    {
                        if (cost[s1] > cost[s2])
                        {
                            fix[s1] = -1;
                            fix[s2] = 10 + s3;
                            cost[s3] -= cost[s2];
                            reduce_cost += cost[s2];
                        }
                        else
                        {
                            fix[s2] = -1;
                            fix[s1] = 10 + s3;
                            cost[s3] -= cost[s1];
                            reduce_cost += cost[s1];
                        }
                    }
                }
            }
            else if (set_var_num[s1] == set_var_num[s3] && set_var_num[s1] < set_var_num[s2])
            {
                int flag1 = 0;
                for (int j = 0; j < set_var_num[s1]; ++j)
                {
                    int v = set_var[s1][j];
                    int flag2 = 0;
                    int flag3 = 0;
                    for (int k = 0; k < var_set_num[v]; ++k)
                    {
                        if (var_set[v][k] == s2)
                        {
                            flag2 = 1;
                        }
                        else if (var_set[v][k] == s3)
                        {
                            flag3 = 1;
                        }
                        if (flag2 == 1 && flag3 == 1)
                            break;
                    }
                    if (flag2 != 1 || flag3 != 1)
                    {
                        flag1 = 1;
                        break;
                    }
                }
                if (flag1 == 0)
                {
                    if (fix[s2] == 1)
                    {
                        fix[s1] = -1;
                        fix[s3] = -1;
                    }
                    else if (cost[s1] >= cost[s2] && cost[s3] >= cost[s2])
                    {
                        fix[s1] = -1;
                        fix[s3] = -1;
                        fix[s2] = 1;
                    }
                    else
                    {
                        if (cost[s1] > cost[s3])
                        {
                            fix[s1] = -1;
                            fix[s3] = 10 + s2;
                            cost[s2] -= cost[s3];
                            reduce_cost += cost[s3];
                        }
                        else
                        {
                            fix[s3] = -1;
                            fix[s1] = 10 + s2;
                            cost[s2] -= cost[s1];
                            reduce_cost += cost[s1];
                        }
                    }
                }
            }
            else if (set_var_num[s3] == set_var_num[s2] && set_var_num[s3] < set_var_num[s1])
            {
                int flag1 = 0;
                for (int j = 0; j < set_var_num[s3]; ++j)
                {
                    int v = set_var[s3][j];
                    int flag2 = 0;
                    int flag3 = 0;
                    for (int k = 0; k < var_set_num[v]; ++k)
                    {
                        if (var_set[v][k] == s2)
                        {
                            flag2 = 1;
                        }
                        else if (var_set[v][k] == s1)
                        {
                            flag3 = 1;
                        }
                        if (flag2 == 1 && flag3 == 1)
                            break;
                    }
                    if (flag2 != 1 || flag3 != 1)
                    {
                        flag1 = 1;
                        break;
                    }
                }
                if (flag1 == 0)
                {
                    if (fix[s1] == 1)
                    {
                        fix[s2] = -1;
                        fix[s3] = -1;
                    }
                    if (cost[s3] >= cost[s1] && cost[s2] >= cost[s1])
                    {
                        fix[s3] = -1;
                        fix[s2] = -1;
                        fix[s1] = 1;
                    }
                    else
                    {
                        if (cost[s2] > cost[s3])
                        {
                            fix[s2] = -1;
                            fix[s3] = 10 + s1;
                            cost[s1] -= cost[s3];
                            reduce_cost += cost[s3];
                        }
                        else
                        {
                            fix[s3] = -1;
                            fix[s2] = 10 + s1;
                            cost[s1] -= cost[s2];
                            reduce_cost += cost[s2];
                        }
                    }
                }
            }
            else if (set_var_num[s3] == set_var_num[s2] && set_var_num[s3] == set_var_num[s1])
            {
                int flag1 = 0;
                for (int j = 0; j < set_var_num[s3]; ++j)
                {
                    int v = set_var[s3][j];
                    int flag2 = 0;
                    int flag3 = 0;
                    for (int k = 0; k < var_set_num[v]; ++k)
                    {
                        if (var_set[v][k] == s2)
                        {
                            flag2 = 1;
                        }
                        else if (var_set[v][k] == s1)
                        {
                            flag3 = 1;
                        }
                        if (flag2 == 1 && flag3 == 1)
                            break;
                    }
                    if (flag2 != 1 || flag3 != 1)
                    {
                        flag1 = 1;
                        break;
                    }
                }
                if (flag1 == 0)
                {
                    if (cost[s1] <= cost[s2] && cost[s1] <= cost[s3])
                    {
                        fix[s1] = 1;
                        fix[s2] = -1;
                        fix[s3] = -1;
                    }
                    else if (cost[s2] <= cost[s1] && cost[s2] <= cost[s3])
                    {
                        fix[s2] = 1;
                        fix[s1] = -1;
                        fix[s3] = -1;
                    }
                    else if (cost[s3] <= cost[s1] && cost[s3] <= cost[s2])
                    {
                        fix[s3] = 1;
                        fix[s2] = -1;
                        fix[s1] = -1;
                    }
                }
            }
        }
    }
    int v, tem_s;
    int top_s;
    memset(var_delete, 0, var_num * sizeof(int));
    for (int s = 0; s < set_num; ++s)
    {
        if (fix[s] == 1)
        {
            for (int i = 0; i < set_var_num[s]; ++i)
            {
                v = set_var[s][i];
                if (var_delete[v] == 1)
                    continue;
                for (int j = 0; j < var_set_num[v]; ++j)
                {
                    tem_s = var_set[v][j];
                    if (fix[tem_s] == 0)
                    {
                        for (int k = 0; k < set_var_num[tem_s]; ++k)
                        {
                            if (set_var[tem_s][k] == v)
                            {
                                top_s = set_var[tem_s][--set_var_num[tem_s]];
                                set_var[tem_s][set_var_num[tem_s]] = set_var[tem_s][k];
                                set_var[tem_s][k] = top_s;
                                if (set_var_num[tem_s] == 0)
                                {
                                    fix[tem_s] = -1;
                                }
                                break;
                            }
                        }
                    }
                }
                var_delete[v] = 1;
            }
        }
        else if (fix[s] > 9)
        {
            for (int i = 0; i < set_var_num[s]; ++i)
            {
                v = set_var[s][i];
                if (var_delete[v] == 1)
                    continue;
                for (int j = 0; j < var_set_num[v]; ++j)
                {
                    tem_s = var_set[v][j];
                    if (fix[tem_s] == 0)
                    {
                        for (int k = 0; k < set_var_num[tem_s]; ++k)
                        {
                            if (set_var[tem_s][k] == v)
                            {
                                top_s = set_var[tem_s][--set_var_num[tem_s]];
                                set_var[tem_s][set_var_num[tem_s]] = set_var[tem_s][k];
                                set_var[tem_s][k] = top_s;
                                if (set_var_num[tem_s] == 0)
                                {
                                    fix[tem_s] = -1;
                                }
                                break;
                            }
                        }
                    }
                }
                var_delete[v] = 1;
            }
        }
    }
}

void WSCP::cover(int selected_var)
{
    int index = index_in_uncover_stack[selected_var];
    int top_set = mypop(uncover_stack);
    index_in_uncover_stack[selected_var] = -1;
    uncover_stack[index] = top_set;
    index_in_uncover_stack[top_set] = index;
}

void WSCP::uncover(int selected_var)
{
    index_in_uncover_stack[selected_var] = uncover_stack_fill_pointer;
    mypush(selected_var, uncover_stack);
}

void WSCP::add_set(int s)
{
    index_in_set_solu[s] = set_count;
    set_solu[set_count++] = s;
}

void WSCP::remove_set(int s)
{
    int index = index_in_set_solu[s];
    int top_set = set_solu[--set_count];
    index_in_set_solu[top_set] = index;
    set_solu[index] = top_set;
    index_in_set_solu[s] = -1;
}

void WSCP::flip(int flip_set)
{
    int cur_set, cur_var;
    long long flip_set_score = score[flip_set];
    long long flip_set_pscore = pscore[flip_set];
    cur_solu[flip_set] = 1 - cur_solu[flip_set];

    if (index_in_zero_stack[flip_set] != -1 && cur_solu[flip_set] == 0)
    {
        int index = index_in_zero_stack[flip_set];
        int top_set = mypop(zero_stack);
        zero_stack[index] = top_set;
        index_in_zero_stack[top_set] = index;
        index_in_zero_stack[flip_set] = -1;
    }

    if (cur_solu[flip_set] == 1)
    {
        add_set(flip_set);
        cur_cost += cost[flip_set];
        for (int i = 0; i < set_var_num[flip_set]; ++i)
        {
            cur_var = set_var[flip_set][i];
             
            if (ratioType4) {
                for (int j = 0; j < var_set_num[cur_var]; ++j)
                    cc[var_set[cur_var][j]] = 1;
            }
           
            ++cover_count[cur_var];
            if (cover_count[cur_var] == 1) //0->1
            {
                for (int j = 0; j < var_set_num[cur_var]; ++j)
                {
                    cur_set = var_set[cur_var][j];
                    score[cur_set] -= weight[cur_var];
                    pscore[cur_set] += weight[cur_var];
                }
                cover_set[cur_var] = flip_set;
                cover(cur_var);
            }
            else if (cover_count[cur_var] == 2) //1->2
            {
                cur_set = cover_set[cur_var];
                score[cur_set] += weight[cur_var];
                if (score[cur_set] == 0 && index_in_zero_stack[cur_set] == -1)
                {
                    index_in_zero_stack[cur_set] = zero_stack_fill_pointer;
                    mypush(cur_set, zero_stack);
                }
                cover_set2[cur_var] = flip_set;
                for (int j = 0; j < var_set_num[cur_var]; ++j)
                {
                    pscore[var_set[cur_var][j]] -= weight[cur_var];
                }
            }
            else if (cover_count[cur_var] == 3)
            {
                pscore[cover_set[cur_var]] += weight[cur_var];
                pscore[cover_set2[cur_var]] += weight[cur_var];
            }
        }
    }
    else
    {
        remove_set(flip_set);
        cur_cost -= cost[flip_set];
        for (int i = 0; i < set_var_num[flip_set]; ++i)
        {
            cur_var = set_var[flip_set][i];
            if(ratioType4){
                for (int j = 0; j < var_set_num[cur_var]; ++j)
                {
                    if(cc[var_set[cur_var][j]] >= 1) cc[var_set[cur_var][j]] = 2;
                    else cc[var_set[cur_var][j]] = 1;
                }
            }
            

            --cover_count[cur_var];

            if (cover_count[cur_var] == 2)
            {
                int flag = 0;
                for (int j = 0; j < var_set_num[cur_var]; ++j)
                {
                    cur_set = var_set[cur_var][j];
                    if (cur_solu[cur_set] == 1)
                    {
                        pscore[cur_set] -= weight[cur_var];
                        if (flag == 0)
                        {
                            cover_set[cur_var] = cur_set;
                            flag = 1;
                        }
                        else
                        {
                            cover_set2[cur_var] = cur_set;
                            break;
                        }
                    }
                }
            }
            else if (cover_count[cur_var] == 1) //2->1
            {
                for (int j = 0; j < var_set_num[cur_var]; ++j)
                {
                    pscore[var_set[cur_var][j]] += weight[cur_var];
                }

                if (cover_set[cur_var] == flip_set)
                {
                    cover_set[cur_var] = cover_set2[cur_var];
                }
                score[cover_set[cur_var]] -= weight[cur_var];
                if (index_in_zero_stack[cover_set[cur_var]] != -1 && score[cover_set[cur_var]] != 0)
                {
                    int top_set = mypop(zero_stack);
                    int index = index_in_zero_stack[cover_set[cur_var]];
                    zero_stack[index] = top_set;
                    index_in_zero_stack[top_set] = index;
                    index_in_zero_stack[cover_set[cur_var]] = -1;
                }
            }
            else if (cover_count[cur_var] == 0) //1->0
            {
                for (int j = 0; j < var_set_num[cur_var]; ++j)
                {
                    cur_set = var_set[cur_var][j];
                    score[cur_set] += weight[cur_var];
                    pscore[cur_set] -= weight[cur_var];
                }
                uncover(cur_var);
            }
        }
        cc[flip_set] = 0;
    }
    pscore[flip_set] = -flip_set_pscore;
    score[flip_set] = -flip_set_score;
}

void WSCP::init()
{
	flag_tabu = 0;
    cur_cost = reduce_cost;
    best_cost = INT_MAX;
    for (int i = 0; i < set_num; ++i)
    {
        best_solu[i] = -1;
        if (fix[i] == 1)
        {
            cur_solu[i] = 1;
            cur_cost += cost[i];
            continue;
        }
        if (fix[i] != 0)
        {
            cur_solu[i] = 0;
            continue;
        }
        score[i] = set_var_num[i];
        pscore[i] = 0;
        cc[i] = 1;
        time_stamp[i] = 0;
        cur_solu[i] = 0;
        index_in_set_solu[i] = -1;
        index_in_zero_stack[i] = -1;
        index_in_unlocked_stack.push_back(-1);
    }
    for (int i = 0; i < var_num; ++i)
    {
    	var_sort[i] = i;
        if (var_delete[i] == 1)
            continue;
        weight[i] = 1;
        cover_count[i] = 0;
    }

    uncover_stack_fill_pointer = 0;
    set_count = 0;
    for (int i = 0; i < var_num; ++i)
    {
        if (var_delete[i] == 1)
            continue;
        index_in_uncover_stack[i] = uncover_stack_fill_pointer;
        mypush(i, uncover_stack);
    }
    check_cur_solu();

    //init feasible solution
    zero_stack_fill_pointer = 0;
    std::sort(var_sort, var_sort + var_num, [this](int x, int y) { return cmp(x, y); });
    for(int i = 0;uncover_stack_fill_pointer > 0;i++)
    {
    	if (cover_count[var_sort[i]] != 0) continue;
        int selected_uncover_var = var_sort[i];//uncover_stack[random() % uncover_stack_fill_pointer];
        int sr = INT_MIN, ct = 1, best_set = -1;
        int ps = INT_MIN;
        for (int i = 0; i < var_set_num[selected_uncover_var]; ++i)
        {
            int cur_set = var_set[selected_uncover_var][i];
			if(fix[cur_set] == -1 || fix[cur_set] == 1) continue;  
            if (sr == INT_MIN || (this->*compareVersion)(sr, ct, score[cur_set], cost[cur_set], ps, pscore[cur_set], weight1, weight2) < 0)
            {
                sr = score[cur_set];
                ps = pscore[cur_set];
                ct = cost[cur_set];
                best_set = cur_set;
            }
            else if ((this->*compareVersion)(sr, ct, score[cur_set], cost[cur_set], ps, pscore[cur_set], weight1, weight2) == 0)
            {
                if (pscore[best_set] < pscore[cur_set])
                    best_set = cur_set;
                else if (pscore[best_set] == pscore[cur_set] && time_stamp[best_set] > time_stamp[cur_set])
                    best_set = cur_set;
            }
        }
        if(best_set == -1) continue;
        flip(best_set);
    }
}

int WSCP::select_set_1()
{
    int best_set = -1;
    for(int i = 1;i <= 100;i++)
    {
        if (set_count == 0) return -1;
        best_set = set_solu[rand_r(&seed) % set_count]; 
        if (fix[best_set]) continue;
        return best_set;
    }
    return best_set;
}

int WSCP::select_set(int is_tabu)
{
    int sr = INT_MIN, ct = 1, i;
    int best_pscore = INT_MIN;
    int best_set = -1;
    if (bms_t < set_count)
    {
        for (int j = 0; j < bms_t; ++j)
        {
            i = set_solu[rand_r(&seed) % set_count]; 
            if (!cur_solu[i])
                continue;
            if (fix[i])
                continue;

            if (sr == INT_MIN || (this->*compareVersion)(sr, ct, score[i], cost[i], best_pscore, pscore[i], weight1, weight2) < 0)
            {
                sr = score[i];
                best_pscore = pscore[i];
                ct = cost[i];
                best_set = i;
            }
            else if ((this->*compareVersion)(sr, ct, score[i], cost[i], best_pscore, pscore[i], weight1, weight2) == 0)
            {
                if (pscore[best_set] < pscore[i])
                    best_set = i;
                
            }
        }
    }
    else
    {
        for (int j = 0; j < set_count; ++j)
        {
            i = set_solu[j];
            if (!cur_solu[i])
                continue;
            if (fix[i])
                continue;

            if (sr == INT_MIN || (this->*compareVersion)(sr, ct, score[i], cost[i], best_pscore, pscore[i], weight1, weight2) < 0)
            {
                sr = score[i];
                best_pscore = pscore[i];
                ct = cost[i];
                best_set = i;
            }
            else if ((this->*compareVersion)(sr, ct, score[i], cost[i], best_pscore, pscore[i], weight1, weight2) == 0)
            {
                if (pscore[best_set] < pscore[i])
                    best_set = i;
                
            }
        }
    }

    return best_set;
}

int WSCP::select_set_new(int is_tabu)
{
    int sr = INT_MIN, ct = 1, i;
    int best_pscore = INT_MIN;
    int best_set = -1;
    if (bms_t < set_count)
    {
        for (int j = 0; j < bms_t; ++j)
        {
            i = set_solu[rand_r(&seed) % set_count]; 
            int index = index_in_unlocked_stack[i];
            if (!(index < unlocked_stack_fill_pointer && unlocked_stack_fill_pointer > 0 && index >= 0)) continue;
            if (!cur_solu[i])
                continue;
            if (fix[i])
                continue;

            if (sr == INT_MIN || compareV2(sr, ct, score[i], cost[i], best_pscore, pscore[i], weight1, weight2) < 0)
            {
                sr = score[i];
                best_pscore = pscore[i];
                ct = cost[i];
                best_set = i;
            }
            else if (compareV2(sr, ct, score[i], cost[i], best_pscore, pscore[i], weight1, weight2) == 0)
            {
                
                if (pscore[best_set] < pscore[i])
                    best_set = i;
                else if (pscore[best_set] == pscore[i] && time_stamp[best_set] > time_stamp[i])
                    best_set = i;
            }
        }
    }
    else
    {
        for (int j = 0; j < set_count; ++j)
        {
            i = set_solu[j];
            int index = index_in_unlocked_stack[i];
            if (!(index < unlocked_stack_fill_pointer && unlocked_stack_fill_pointer > 0 && index >= 0)) continue;
            if (!cur_solu[i])
                continue;
            if (fix[i])
                continue;

            if (sr == INT_MIN || compareV2(sr, ct, score[i], cost[i], best_pscore, pscore[i], weight1, weight2) < 0)
            {
                sr = score[i];
                best_pscore = pscore[i];
                ct = cost[i];
                best_set = i;
            }
            else if (compareV2(sr, ct, score[i], cost[i], best_pscore, pscore[i], weight1, weight2) == 0)
            {
                if (pscore[best_set] < pscore[i])
                    best_set = i;
                else if (pscore[best_set] == pscore[i] && time_stamp[best_set] > time_stamp[i])
                    best_set = i;
            }
        }
    }
    return best_set;
}

int WSCP::select_set_from_zero_stack_new()
{
    int best_set = -1, best_cost = 0;
    int tem_set;
    if (bms_t < zero_stack_fill_pointer)
    {
        for (int i = 0; i < bms_t; ++i)
        {
            tem_set = zero_stack[rand_r(&seed) % zero_stack_fill_pointer]; //zero_stack[rand() % zero_stack_fill_pointer];
            int index = index_in_unlocked_stack[tem_set];
            if (!(index < unlocked_stack_fill_pointer && unlocked_stack_fill_pointer > 0 && index >= 0)) continue;
            if (cost[tem_set] > best_cost)
            {
                best_cost = cost[tem_set];
                best_set = tem_set;
            }
        }
    }
    else
    {
        for (int i = 0; i < zero_stack_fill_pointer; ++i)
        {
            tem_set = zero_stack[i];
            int index = index_in_unlocked_stack[tem_set];
            if (!(index < unlocked_stack_fill_pointer && unlocked_stack_fill_pointer > 0 && index >= 0)) continue;
            if (cost[tem_set] > best_cost)
            {
                best_cost = cost[tem_set];
                best_set = tem_set;
            }
        }
    }
    return best_set;
}

int WSCP::count_2cc(int var)
{
	int num = var_set_num[var];
	int ans = 1;
	for(int i = 0;i < 50;i++)
	{
		int cur_set = var_set[var][rand_r(&seed) % num];
		if (cc[cur_set] == 2) ans++;
	}
	return ans;
}

int WSCP::select_set_Novelty()
{
	int var1 = uncover_stack[rand_r(&seed) % uncover_stack_fill_pointer];
    if (ratioType4)
    {
        if (uncover_stack_fill_pointer != 1)
        {
            int var2 = uncover_stack[rand_r(&seed) % uncover_stack_fill_pointer];
            while (var1 == var2) var2 = uncover_stack[rand_r(&seed) % uncover_stack_fill_pointer];
            int count1 = count_2cc(var1);
            int count2 = count_2cc(var2);
            double p_cc = (double)count2 / (double)(count1 + count2);
            if ((rand_r(&seed) % 100)/(double)101 < p_cc) var1 = var2;
        }
    }
	
	
    int selected_uncover_var = var1;
    int sr = INT_MIN, ct = 1, best_set = -1;
    int ps = INT_MIN;
    for (int i = 0; i < var_set_num[selected_uncover_var]; ++i)
    {
        int cur_set = var_set[selected_uncover_var][i];

		if (step - time_stamp[cur_set] < tabu_length)
            continue;
        if (sr == INT_MIN || (this->*compareVersion)(sr, ct, score[cur_set], cost[cur_set], ps, pscore[cur_set], weight1, weight2) < 0)
        {
            sr = score[cur_set];
            ps = pscore[cur_set];
            ct = cost[cur_set];
            best_set = cur_set;
        }
        else if ((this->*compareVersion)(sr, ct, score[cur_set], cost[cur_set], ps, pscore[cur_set], weight1, weight2) == 0)
        {
        	
            if(rand_r(&seed) % 2 == 0) best_set = cur_set;
            
        }
    }
    int sec_best_set = -1;
    sr = INT_MIN;
    ct = 1;
    ps = INT_MIN;
	if(best_set == -1)
	{
		flag_tabu = 0;
		return selected_uncover_var;
	}
    for (int i = 0; i < var_set_num[selected_uncover_var]; ++i)
    {
        int cur_set = var_set[selected_uncover_var][i];
        if (step - time_stamp[cur_set] < tabu_length)
            continue;
        if (cur_set == best_set)
            continue;
        if (sr == INT_MIN || (this->*compareVersion)(sr, ct, score[cur_set], cost[cur_set], ps, pscore[cur_set], weight1, weight2) < 0)
        {
            sr = score[cur_set];
            ps = pscore[cur_set];
            ct = cost[cur_set];
            sec_best_set = cur_set;
        }
        else if ((this->*compareVersion)(sr, ct, score[cur_set], cost[cur_set], ps, pscore[cur_set], weight1, weight2) == 0)
        {
 
            if(rand_r(&seed) % 2 == 0) sec_best_set = cur_set; 
            
        }
    }
    long long best_s = (weight1 * score[best_set] + weight2 * pscore[best_set]) * cost[sec_best_set];
    long long sec_s = (weight1 * score[sec_best_set] + weight2 * pscore[sec_best_set]) * cost[best_set];
    double pp = ((double)sec_s + 1) / ((double)best_s + sec_s + 2);
    if (sec_best_set != -1 && (rand_r(&seed) % 100)/(double)101 < pp * novelty_p)
        best_set = sec_best_set;
    flag_tabu = 1;
    return best_set;
}

int WSCP::select_set_Novelty_best()
{
	int var1 = uncover_stack[rand_r(&seed) % uncover_stack_fill_pointer];
	
	
    int selected_uncover_var = var1;
    int sr = INT_MIN, ct = 1, best_set = -1;
    int ps = INT_MIN;
    for (int i = 0; i < var_set_num[selected_uncover_var]; ++i)
    {
        int cur_set = var_set[selected_uncover_var][i];
        
		if (step - time_stamp[cur_set] < tabu_length)
            continue;
        if (sr == INT_MIN || (this->*compareVersion)(sr, ct, score[cur_set], cost[cur_set], ps, pscore[cur_set], weight1, weight2) < 0)
        {
            sr = score[cur_set];
            ps = pscore[cur_set];
            ct = cost[cur_set];
            best_set = cur_set;
        }
        else if ((this->*compareVersion)(sr, ct, score[cur_set], cost[cur_set], ps, pscore[cur_set], weight1, weight2) == 0)
        {
        	 
            if(rand_r(&seed) % 2 == 0) best_set = cur_set;
            
        }
    }

	if(best_set == -1)
	{
		flag_tabu = 0;
		return selected_uncover_var;
	}
    flag_tabu = 1;
    return best_set;
}

int WSCP::select_set_Novelty_no_tabu(int var1)
{
	
    int selected_uncover_var = var1;
    int sr = INT_MIN, ct = 1, best_set = -1;
    int ps = INT_MIN;
    for (int i = 0; i < var_set_num[selected_uncover_var]; ++i)
    {
        int cur_set = var_set[selected_uncover_var][i];
        if(best_set == -1) best_set = cur_set;
		if (sr == INT_MIN || (this->*compareVersion)(sr, ct, score[cur_set], cost[cur_set], ps, pscore[cur_set], weight1, weight2) < 0)
        {
            sr = score[cur_set];
            ps = pscore[cur_set];
            ct = cost[cur_set];
            best_set = cur_set;
        }
        else if ((this->*compareVersion)(sr, ct, score[cur_set], cost[cur_set], ps, pscore[cur_set], weight1, weight2) == 0)
        {
        	if(rand_r(&seed) % 2 == 0) best_set = cur_set; 
            
        }
    }
    int sec_best_set = -1;
    sr = INT_MIN;
    ct = 1;
    ps = INT_MIN;
    for (int i = 0; i < var_set_num[selected_uncover_var]; ++i)
    {
        int cur_set = var_set[selected_uncover_var][i];
        if (cur_set == best_set)
            continue;
        if (sr == INT_MIN || (this->*compareVersion)(sr, ct, score[cur_set], cost[cur_set], ps, pscore[cur_set], weight1, weight2) < 0)
        {
            sr = score[cur_set];
            ps = pscore[cur_set];
            ct = cost[cur_set];
            sec_best_set = cur_set;
        }
        else if ((this->*compareVersion)(sr, ct, score[cur_set], cost[cur_set], ps, pscore[cur_set], weight1, weight2) == 0)
        {
        	if(rand_r(&seed) % 2 == 0) sec_best_set = cur_set; 
            
        }
    }
    long long best_s = (weight1 * score[best_set] + weight2 * pscore[best_set]) * cost[sec_best_set];
    long long sec_s = (weight1 * score[sec_best_set] + weight2 * pscore[sec_best_set]) * cost[best_set];
    double pp = ((double)sec_s + 1) / ((double)best_s + sec_s + 2);
    if (sec_best_set != -1 && (rand_r(&seed) % 100)/(double)101 < pp * novelty_p)
        best_set = sec_best_set;
    flag_tabu = 1;
    return best_set;
}

int WSCP::select_set_from_zero_stack()
{
    int best_set = -1, best_cost = 0;
    int tem_set;
    if (bms_t < zero_stack_fill_pointer)
    {
        for (int i = 0; i < bms_t; ++i)
        {
            tem_set = zero_stack[rand_r(&seed) % zero_stack_fill_pointer]; //zero_stack[rand() % zero_stack_fill_pointer];
            if (cost[tem_set] > best_cost)
            {
                best_cost = cost[tem_set];
                best_set = tem_set;
            }
        }
    }
    else
    {
        for (int i = 0; i < zero_stack_fill_pointer; ++i)
        {
            tem_set = zero_stack[i];
            if (cost[tem_set] > best_cost)
            {
                best_cost = cost[tem_set];
                best_set = tem_set;
            }
        }
    }
    return best_set;
}

int WSCP::local_greedy_addition() {
    int best_set = -1, cur_set;
    int selected_uncover_var = uncover_stack[rand_r(&seed) % uncover_stack_fill_pointer]; //uncover_stack[rand() % uncover_stack_fill_pointer];
    int total_num = var_set_num[selected_uncover_var];
    int BMS_T = 1000; 
    if(total_num < BMS_T) {
        for(int i=0;i<total_num;i++) {
            cur_set = var_set[selected_uncover_var][i];
            if(best_set == -1) best_set = cur_set;
            if (step - time_stamp[cur_set] < tabu_length)
                continue;
            if(compare(cur_set, best_set, 1) > 0)
                best_set = cur_set;
        }
    } else {
        for(int i=0;i<BMS_T;i++) {
            cur_set = var_set[selected_uncover_var][rand_r(&seed)%total_num];
            if(best_set == -1) best_set = cur_set;
            if (step - time_stamp[cur_set] < tabu_length)
                continue;
            if(compare(cur_set, best_set, 1) > 0)
                best_set = cur_set;
        } 
    }
    return best_set;
}


void WSCP::update_weight()
{
    int cur_set, cur_var;
    for (int i = 0; i < uncover_stack_fill_pointer; ++i)
    {
        cur_var = uncover_stack[i];
        weight[cur_var]++;
        for (int j = 0; j < var_set_num[cur_var]; ++j)
        {
            cur_set = var_set[cur_var][j];
            score[cur_set]++;
        }
    }
}


bool WSCP::cplex_improve()
{
    return cplex_improve(cur_solu, cur_cost);
}

bool WSCP::cplex_improve(const int *start_solu, long long start_cost)
{
    if (get_runtime() > cutoff_time)
        return false;
    if (start_solu == nullptr)
        start_solu = cur_solu;
    if (start_cost < 0)
        start_cost = cur_cost;

    IloEnv env;
    try
    {
        IloModel model(env);
        IloBoolVarArray x(env, set_num);

        for (int i = 0; i < set_num; ++i)
        {
            if (fix[i] == 1)
                model.add(x[i] == 1);
            else if (fix[i] != 0)
                model.add(x[i] == 0);
        }

        for (int v = 0; v < var_num; ++v)
        {
            if (var_delete[v] == 1)
                continue;
            IloExpr cover(env);
            for (int j = 0; j < var_set_num[v]; ++j)
            {
                int s = var_set[v][j];
                cover += x[s];
            }
            model.add(cover >= 1);
            cover.end();
        }

        IloExpr obj_expr(env);
        for (int i = 0; i < set_num; ++i)
            obj_expr += cost[i] * x[i];

        model.add(IloMinimize(env, obj_expr));

        long long improve_limit = start_cost - reduce_cost ;
        if (improve_limit < 0)
        {
            obj_expr.end();
            env.end();
            return false;
        }
        IloNum improve_ub = static_cast<IloNum>(improve_limit);
        model.add(IloRange(env, -IloInfinity, obj_expr, improve_ub));

        IloCplex cplex(model);
        cplex.setParam(IloCplex::Param::Advance, 1);
        cplex.setParam(IloCplex::Param::Threads, 1);
        cplex.setOut(env.getNullStream());

        double remaining_time = cutoff_time - get_runtime();
        if (remaining_time <= 0.0)
        {
            obj_expr.end();
            env.end();
            return false;
        }
        cplex.setParam(IloCplex::Param::TimeLimit, remaining_time);

        IloNumVarArray start_vars(env);
        IloNumArray start_vals(env);
        for (int i = 0; i < set_num; ++i)
        {
            start_vars.add(x[i]);
            int val = start_solu[i];
            if (fix[i] == 1)
                val = 1;
            else if (fix[i] != 0)
                val = 0;
            start_vals.add(val ? 1.0 : 0.0);
        }
        cplex.addMIPStart(start_vars, start_vals);

        if (cplex.solve())
        {
            bool changed = false;
            for (int i = 0; i < set_num; ++i)
            {
                if (fix[i] != 0)
                    continue;
                int val = (cplex.getValue(x[i]) > 0.5) ? 1 : 0;
                if (val != cur_solu[i])
                {
                    flip(i);
                    time_stamp[i] = step;
                    changed = true;
                }
            }
            obj_expr.end();
            env.end();
            return changed;
        }

        obj_expr.end();
        env.end();
        return false;
    }
    catch (...)
    {
        env.end();
        return false;
    }
}

void WSCP::local_search1()
{

   
    if (ratioType3 && is_update && flag_cplex)
    {
         if (cplex_improve(best_solu, best_cost) && cur_cost < best_cost)
         {
             best_cost = cur_cost;
             best_time = get_runtime();
             for (int i = 0; i < set_num; ++i)
                best_solu[i] = cur_solu[i];
         }
     }

	bool temp_flag = 0;
    int InnerDepth = 50;
    int in_nonimpr = 0;
    int t = 0;
    int index;
    if (set_count == 0) t = 0;
    else if((int)(0.1 * set_count) == 0) t = 1;
    else t = 0.2 * set_count;
    
    int cnt = 0;
    while(cnt <= t)
    {
        if (set_count == 0) break;
        int selected_set = select_set_1();
        if (selected_set == -1) break;
        int var_n = set_var_num[selected_set];
        if (var_n < 0) break;
        int selected_var = set_var[selected_set][rand_r(&seed) % var_n]; //set_var[selected_set][rand() % var_n];
		if (selected_var < 0) break;
		if (var_set_num[selected_var] <= 0) break;
		int number = 0;
        for (int i = 0; i < var_set_num[selected_var]; ++i)
        {
            int cur_set = var_set[selected_var][i];
            if (cur_solu[cur_set] == 1)
            {
                flip(cur_set);
                ++cnt;
                ++number;
                if (cnt >= t) break;
            }
            if (cnt >= t) break;
        }
        if (number == 0) break;
    }
    
    while(uncover_stack_fill_pointer > 0)
    {
		int best_set = select_set_Novelty();
		
		if (!flag_tabu) 
		{
			best_set = select_set_Novelty_no_tabu(best_set);
		}
		if (best_set == -1) break;
		flip(best_set);
		index_in_unlocked_stack[best_set] = unlocked_stack_fill_pointer;
		mypush(best_set,unlocked_stack);
    }
    unlocked_stack_fill_pointer = 0;
    long long now_cost = 0;
    
    while(in_nonimpr < InnerDepth)
    {
    	if (uncover_stack_fill_pointer == 0)
		{
			while (zero_stack_fill_pointer > 0)
            {
                int flip_set = select_set_from_zero_stack_new();
                if (flip_set == -1) return;
                flip(flip_set);
                index = index_in_unlocked_stack[flip_set];
                if (index < unlocked_stack_fill_pointer && unlocked_stack_fill_pointer > 0 && index >= 0)
				{
					int top_set = mypop(unlocked_stack);
					index_in_unlocked_stack[flip_set] = -1;
					uncover_stack[index] = top_set;
					index_in_uncover_stack[top_set] = index;
				}
            }
            if (cur_cost < now_cost) in_nonimpr = 0;
            if (cur_cost < best_cost)
            {
                best_cost = cur_cost;
                best_time = get_runtime();
                for(int i = 0;i < set_num;i++)
                {
                	for(int i = 0;i < set_num;i++)
                		best_solu[i] = cur_solu[i];
                	in_nonimpr = 0;
				}
			}
		}
		
		
		int flip_set = select_set_new(0);
		if (flip_set == -1) return;
		flip(flip_set);
		index = index_in_unlocked_stack[flip_set];
        if (index < unlocked_stack_fill_pointer && unlocked_stack_fill_pointer > 0 && index >= 0)
		{
			int top_set = mypop(unlocked_stack);
			index_in_unlocked_stack[flip_set] = -1;
			uncover_stack[index] = top_set;
			index_in_uncover_stack[top_set] = index;
		}
		flip_set = select_set_new(0);
		
		if (flip_set == -1) in_nonimpr++;
		else 
		{
			flip(flip_set);
			index = index_in_unlocked_stack[flip_set];
        	if (index < unlocked_stack_fill_pointer && unlocked_stack_fill_pointer > 0 && index >= 0)
			{
				int top_set = mypop(unlocked_stack);
				index_in_unlocked_stack[flip_set] = -1;
				uncover_stack[index] = top_set;
				index_in_uncover_stack[top_set] = index;
			}
		}
		while(uncover_stack_fill_pointer != 0)
		{
			int best_set = select_set_Novelty();
			
			if (!flag_tabu) 
				best_set = select_set_Novelty_no_tabu(best_set);
			if (best_set == -1) return;
			
			flip(best_set);
			index_in_unlocked_stack[best_set] = unlocked_stack_fill_pointer;
		    mypush(best_set,unlocked_stack);
		}
		in_nonimpr++;
	}
}

void WSCP::local_search()
{
	gap = 10; 
	int gap_max = gap;
	gap = 0; 
	int OutterDepth = set_num;
	int OutterDepth_new = var_num * 15;
    int out_nonimpr = 0;
    int out_nonimpr_new = 0;
    long long cur_now_cost = 0x7fffffff;
    int flag_mode1 = 0;
    int flag_1 = 0;
	int num_del = 0;
	int step_del = 0;
	int step_max = 10;
	int zero_time = 0;
	int max_zero_time = 10; 
    for (tries = 0; tries < max_tries; ++tries)
    {
    	
        for (step = 0; step < max_step; ++step)
        {
        	out_nonimpr++;
        	out_nonimpr_new++;
        	if (out_nonimpr == OutterDepth)
            {
            	if (gap == 0) 
				{
					gap = gap_max;
					zero_time++;
					if(zero_time == max_zero_time) out_nonimpr_new = OutterDepth_new;
				}
				else
				{
					gap = gap / 2;
            		out_nonimpr = 0;
				}
			}
			if (OutterDepth_new == out_nonimpr_new) 
			{
				cur_now_cost = 0x7fffffff;
				if (is_update) gap_max = gap_max / 2;
				else gap_max = gap_max / 5;
				gap = std::min(gap,gap_max); 
				if (gap_max <= 40) OutterDepth = set_num * 2;
				if (!ratioType1) local_search1();
				zero_time = 0;
				is_update = 0;
				out_nonimpr_new = 0;
			}
            if (uncover_stack_fill_pointer == 0)
            {
                if (zero_stack_fill_pointer > 0)
                {
                    int flip_set = select_set_from_zero_stack();
                    if (flip_set == -1){
                        return;
                    }
                    flip(flip_set);
                    time_stamp[flip_set] = step;
                    continue;
                }
                else if (cur_cost < best_cost)
                {
                	is_update = 1;
                    best_cost = cur_cost;
                    cur_now_cost = cur_cost;
                    best_time = get_runtime();
                    best_count= set_count;
                    gap = gap * 2;
                    gap = std::min(gap,gap_max);
                    out_nonimpr = 0;
                    out_nonimpr_new = 0;
                    for (int i = 0; i < set_num; ++i)
                        best_solu[i] = cur_solu[i];
                }
                else if(cur_cost < cur_now_cost)
				{
					cur_now_cost = cur_cost;
                    gap = gap * 2;
                    gap = std::min(gap,gap_max);
                    out_nonimpr = 0;
                    out_nonimpr_new = 0;
				}
            }
            if (uncover_stack_fill_pointer > 0)
            {
                int best_set;
                if (ratioType3) best_set =local_greedy_addition();
                else {
                    if (ratioType2) best_set = select_set_Novelty_best();
                    else best_set = select_set_Novelty();
				    if (!flag_tabu) 
				    {
                        
					    best_set = select_set_Novelty_no_tabu(best_set);
				    }
                }
                if (cur_cost + cost[best_set] + gap >= best_cost || flag_mode1)
                {
                	num_del++;
                	
                	int flip_set;
                    if (ratioType1|| ratioType2) flip_set = remove_best_score(); 
                    else flip_set = select_set(0);
        
                    
                    if (flip_set < 0)
                    {
                    	
                        flip(best_set);
                        time_stamp[best_set] = step;
                        continue;
                    }
					if (flip_set == -1) {
                        
                        return;
                    }
                    else if (compare1(-score[flip_set], cost[flip_set], score[best_set], cost[best_set]) < 0 || flag_1)
                    {
                    	num_del = 0;
                    	step_del = 0;
                    	flag_1 = 0;
                    	flag_mode1 = 0;
                        if (flip_set == -1 || best_set == -1) {
                            return;
                        }
                        flip(flip_set);
                        flip(best_set);
                        time_stamp[flip_set] = step;
                        time_stamp[best_set] = step;
                        if (step % 1000 == 0)
                        {
                            if (get_runtime() > cutoff_time){
                                return;
                            }
                        }
                        continue;
                    }
                    else
                    {
                    	if(double(rand_r(&seed) % 101) / 100.0 < 0.5) flag_1 = 1;
                    	flag_mode1 = 1;
                    	if (ratioType1|| ratioType2) flip_set = remove_best_score(); 
                        else flip_set = select_set(0);
                        
                        if (flip_set == -1)
                        {
							flag_mode1 = 0;
                            continue;
                        	return;
						}
                            
                        flip(flip_set);
                        time_stamp[flip_set] = step;

                        if (ratioType1|| ratioType2) flip_set = remove_best_score(); 
                        else flip_set = select_set(0);
                      
						
                        if (flip_set == -1)
                        {
                        	flag_mode1 = 0;
                            continue;
							return;
						}
                            
                        if (flip_set >= 0)
                        {
                            flip(flip_set);
                            time_stamp[flip_set] = step;
                        }

                        if (step % 1000 == 0)
                        {
                            if (get_runtime() > cutoff_time) {
                                return;
                            }
                        }
                        continue;
                    }
                }
                else
                {
                    if (best_set == -1) {
                        return;
                    }
                    flip(best_set);
                    update_weight();
                    time_stamp[best_set] = step;
                    if (step % 1000 == 0)
                    {
                        if (get_runtime() > cutoff_time){
                            return;
                        }
                    }
                    continue;
                }
            }
            int flip_set;
            if (ratioType1|| ratioType2) flip_set = remove_best_score(); 
            else flip_set = select_set(0);
          
            
            if (flip_set == -1) {
                return;
            }
            flip(flip_set);

            time_stamp[flip_set] = step;
            if (step % 1000 == 0)
            {
                if (get_runtime() > cutoff_time) {
                    return;
                }
            }
        }
    }
}

int WSCP::remove_best_score() {
    int best_set = -1, cur_set;
    if(bms_t < set_count) {
        for (int j = 0; j < bms_t; ++j)
        {
            cur_set = set_solu[rand_r(&seed)%set_count]; //set_solu[rand()%set_count]; //set_solu[rng.next(set_count)];
            if(fix[cur_set]) continue;
            if(!cur_solu[cur_set]) continue;
            int res = compare(best_set, cur_set, -1);
            if(best_set == -1 || res < 0) 
                best_set = cur_set;
        }
    } else {
        for (int j = 0; j < set_count; ++j)
        {
            cur_set = set_solu[j];
            if(fix[cur_set]) continue;
            if(!cur_solu[cur_set]) continue;
            int res = compare(best_set, cur_set, -1);
            if(best_set == -1 || res < 0) 
                best_set = cur_set;
        }
    }
    return best_set;
}

void WSCP::check_solu()
{
    long long tem_cost = 0;
    for (int i = 0; i < set_num; ++i)
    {
        if (fix[i] == 1 && best_solu[i] != 1)
            cout << "wrong 1" << endl;
        if (fix[i] < 0 && best_solu[i] != 0)
            cout << "wrong 2" << endl;
        if (fix[i] > 9)
        {
            int s = fix[i] - 10;
            if (best_solu[i] != 0)
                cout << "wrong 3" << endl;
            if (best_solu[s] == 0)
                best_solu[i] = 1;
        }
    }
    for (int i = 0; i < var_num; ++i)
    {
        int flag = 0;
        for (int j = 0; j < var_set_num[i]; ++j)
        {
            if (best_solu[var_set[i][j]] == 1)
            {
                flag = 1;
                break;
            }
        }
        if (flag == 0)
        {
            cout << "best solu is wrong in 0 " << endl;
            return;
        }
    }
    for (int i = 0; i < set_num; ++i)
    {
        if (best_solu[i] == 1)
            tem_cost += org_cost[i];
    }

    tem_cost = 0;
    for (int i = 0; i < var_num; ++i)
    {
        int flag = 0;
        for (int j = 0; j < var_set_num[i]; ++j)
        {
            if (best_solu[var_set[i][j]] == 1 || fix[var_set[i][j]] == 1 || fix[var_set[i][j]] == -2)
            {
                flag = 1;
                break;
            }
        }
        if (flag == 0)
        {
            cout << "best solu is wrong in 1" << endl;
            return;
        }
    }
    for (int i = 0; i < set_num; ++i)
    {
        if (fix[i] == 1 || best_solu[i] == 1)
            tem_cost += org_cost[i];

        if (fix[i] == -2)
        {

            for (int j = 0; j < org_set_var_num[i]; ++j)
            {
                int flag = 0;
                int v = set_var[i][j];
                for (int k = 0; k < var_set_num[v]; ++k)
                {
                    if (fix[var_set[v][k]] == 1 || best_solu[var_set[v][k]] == 1)
                    {
                        flag = 1;
                        break;
                    }
                }
                if (flag == 0)
                {
                    tem_cost += org_cost[i];
                    break;
                }
            }
        }
    }

    if (tem_cost != best_cost)
    {
        cout << "best solu is wrong in 5 " << endl;
    	cout << tem_cost << " " << best_cost << endl;
	}
    tem_cost = 0;
    for (int i = 0; i < set_num; ++i)
    {
        if (fix[i] == 1 || best_solu[i] == 1)
            tem_cost += org_cost[i];
    }
    for (int i = 0; i < var_num; ++i)
    {
        int flag = 0;
        for (int j = 0; j < var_set_num[i]; ++j)
        {
            if (fix[var_set[i][j]] == 1 || best_solu[var_set[i][j]] == 1)
            {
                flag = 1;
                break;
            }
        }
        if (flag == 0)
        {
            if (var_set_num[i] == 2)
            {
                tem_cost += org_cost[i];
            }
        }
    }
}

void WSCP::output_solu() {
    cout << "solution:" << endl;
	for (int i = 0; i < set_num; ++i)
    {
        if (best_solu[i] == 1)
            cout << i << " ";
    }
	cout << endl;
}

void WSCP::check_cur_solu()
{
    long long tem_cost = 0;
    for (int i = 0; i < set_num; ++i)
    {
        if (fix[i] == 1 && cur_solu[i] != 1)
            cout << "wrong 11" << endl;
        if (fix[i] < 0 && cur_solu[i] != 0)
            cout << "wrong 22" << endl;
        if (fix[i] > 9)
        {
            if (cur_solu[i] != 0)
                cout << "wrong 33" << endl;
        }
    }

    for (int i = 0; i < set_num; ++i)
    {
        if (cur_solu[i] == 1)
            tem_cost += org_cost[i];
        else if (fix[i] > 9 && cur_solu[fix[i] - 10] == 0)
            tem_cost += org_cost[i];
    }
    if (tem_cost != cur_cost)
        cout << "some thing wrong in cur cost and tem_cost " << endl;
}

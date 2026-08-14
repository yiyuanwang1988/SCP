# FastSC: A Fast Local Search Algorithm for the Set Covering Problem

We propose an efficient local search algorithm, FastSC, for solving the set covering problem.




## 1. CPLEX Environment Configuration
The algorithm depends on the [IBM ILOG CPLEX library](https://www.ibm.com/products/ilog-cplex-optimization-studio). Please ensure:

-   CPLEX is correctly installed at `/opt/ibm/ILOG/CPLEX_Studio2211/`  
    
-   If CPLEX is installed in a different location, you will need to modify the library paths in the Makefile.

## 2. compile

To compile the source code and generate the executable ./FastSC, use the following command:
```
make
```



## 3. Usage

To execute the code, use the following command format:
```
./FastSC <instance> <time> <seed> <is_AC> [use_cplex]
```
The above parameters are explained as follows:
-   instance: path of the input graph.
-   time: time limit (seconds).
-   seed: random seed.
-   is_AC: Whether the instance is Academic (1) or not (0). If is_AC = 1, subset costs are generated as (i+1) % 200 + 1.
-   use_cplex: Whether to use CPLEX (1) or not (0).


### **For example**

```
./FastSC ./AC/AC_01_coverC.txt 10 1 1 0
```

### **Output**

```
instance: ./AC/AC_01_coverC.txt cost: 28 seed: 1 time: 0.01 set number: 10
solution:
1 199 200 399 599 803 999 1000 1805 2403 

```

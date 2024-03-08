#include <iostream>
#include <cstdio>
#include <vector>
#include <algorithm>
#include <cmath>
#include <set>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include <utility>
#include <ctime>
#include <omp.h>
#include "Edge.h"

class Graph {
public:
    explicit Graph(FILE* file, int t);

    explicit Graph(int node_num, const std::vector<std::pair<int, int> > &edges);

    Graph(const Graph& g);

    ~Graph();

    void AC();

    void SC();

    void Shell_PDC();

    void Decom();

    void Check();

    void PDC();

    void PDC_org();


    


private:
    int n;
    double m;
    int lmax = 0;
    int kmax = 0;
    int dmax = 0;
    long long iterations = 0;

    int num_of_thread;
    std::vector<std::vector<int> > *adj;

    std::vector<int> *deg;

    int *outCdeg {nullptr};
    int *inCdeg {nullptr};

    int *outTdeg {nullptr};
    int *inTdeg {nullptr};
    
    std::vector<int> max_deg;

    std::vector<int> iH;

    std::vector<int> oH;

    std::vector<int> OC;

    std::vector<int> *iHk;

    std::vector<int> *iHk1;

    std::vector<int> *sH; //slyline core number
    std::vector<int> *tmp_results; //slyline core number
    
    omp_lock_t *lock; 

    std::vector<bool> changed; //for each vertex

    void add_edge(int u, int v);

    void sort_adj_list();

    void init_adj_deg();

    bool HIndex(int v);

    bool outHIndex(int v);

    bool oHIndex(int v);

    bool outCore(int v);

    bool Refine(int v);

    void prune(int k);
    bool Refine(int v, int k);
    void getInDegree(int k);

    bool BiSearch(std::vector<int> &NeighborSkyline, int &k, int &l);

    void update(int v);
    bool NoVerify(int v);


    /////////// for peeling ////////////
    int kCoreMax = 0;
    int Kmd;

    std::vector<int> Kdeg, kBin; //used to calculate all (k,0)-cores
    //Kdeg stores all vertices's max possible K, where there is a (K,0)-core contains it
    //kBin is used to do a bucket sort for all vertices in ascending order of their Kdeg

    std::vector<int> bin, pos;//used to calculate rows
    int *vert {nullptr};//used to calculate rows
    //give a certain K,
    //bin is a bucked sort in ascending number of vertices' max L where there is a (K,L)-core contains it
    //pos: which bin a vertex is in
    //vert: the array of all vertices in this row(namely in (K,0)-core)

    std::vector<int> subLDegree, subKDegree;           //record a vertex 's out and in degree in current subgraph
    std::vector<bool> visited, isInSub;     //isInsub used in getRow, judge whether a vetex satisfy the basic k constraint (namely in (k,0)-core)

    //we do D-core decomposition row by row
    //one row is, a certain K, all the (K, l) cores sorted by ascending order of l
    std::vector<int> rowVert;                          //store vertices in one row                   
    int rowResultPos;                    //a mark variable used to fill a row result into rowResult[]
    std::vector<int> sorted_ID;                //a strictly sorted array, sorted by out-core number
    std::vector<int> rowResult;                //a strictly sorted array, sorted by 1 core number
    std::vector<int> compulsoryVisit;
    int compulsorySubL,compulsoryNum;

    void calculate_kl(int k);
    std::vector<int> inDependentRows();
    void visit(int k, int v, int compulsory);
    std::vector<int> sort_degree();

    //functions for par_peel
    std::vector<int> Parpeel_org(int k);
    std::vector<int> Parpeel(int k);
    std::vector<int> Parpeel0();
    std::vector<int> Parpeel(int k, std::vector<int> upper);

    void pkc();
    void plc();

    //shell
    std::vector<int> kshell(std::vector<int> Din);
    

};

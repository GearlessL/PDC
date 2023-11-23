#include "Graph.h"
#include "omp.h"
#include <cassert>
#include <fstream>
#include <chrono>

// export OMP_WAIT_POLICY=passive

Graph::Graph(FILE* file, int t) {
    num_of_thread = t;
    clock_t begin = clock();
    double edge_num;
    fscanf(file, "%d%lf", &n, &edge_num);
    init_adj_deg();


    m = 0;
    for (long i = 0; i < edge_num; i++) {
        int u, v;
        fscanf(file, "%d%d", &u, &v);
        // if(u!=v)
        add_edge(u, v);
    }

    sort_adj_list();



    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    printf("graph construction: %.4f\n", elapsed_secs);
}

Graph::Graph(const int node_num, const std::vector<std::pair<int, int> > &edges) {
    n = node_num;
    init_adj_deg();

    for (auto& edge : edges) {
        add_edge(edge.first, edge.second);
    }

}

void Graph::init_adj_deg() {
    max_deg.resize(2, 0);
    iHk = new std::vector<int>[n]; //for each vertex
    iHk1 = new std::vector<int>[n]; //for each vertex
    sH = new std::vector<int>[n];
    tmp_results = new std::vector<int>[n];
    lock = new omp_lock_t[n];
    for(int i = 0; i < n; i++){
        omp_init_lock(&lock[i]);
    }

    adj = new std::vector<std::vector<int>>[2]; // 0 out edges, 1 in edges
    deg = new std::vector<int>[2]; // 0 out-degree, 1 in-degree
    m = 0;
    for (int i = 0; i < 2; i++) {
        adj[i].resize(static_cast<unsigned long>(n));
        deg[i].resize(static_cast<unsigned long>(n), 0);
    }
    changed.resize(n,true);
}

Graph::~Graph() {
    for (int i = 0; i < 2; i++) {
        adj[i].clear();
        deg[i].clear();
    }
    delete[] adj;
    delete[] deg;
    max_deg.clear();
    for(int i = 0; i < n; i++){
        omp_destroy_lock(&lock[i]);
    }
    delete[] sH;
    delete[] tmp_results;
    delete[] lock;
}

void Graph::sort_adj_list() {
    for (int i = 0; i < 2; i++) {
        for (int v = 0; v < n; v++) {
            std::sort(adj[i][v].begin(), adj[i][v].end(),
                      [&](const int& a, const int& b) -> bool
                      {
                          return deg[1 - i][a] > deg[1 - i][b];
                      });
        }
    }
}

inline void Graph::add_edge(int u, int v) {
    adj[0][u].push_back(v);
    deg[0][u] += 1;
    adj[1][v].push_back(u);
    deg[1][v] += 1;
    max_deg[0] = std::max(max_deg[0], deg[0][u]);
    max_deg[1] = std::max(max_deg[1], deg[1][v]);
    ++m;
}



bool Graph::HIndex(int v){
    changed[v] = false;
    std::vector<int> I;
    I.clear();
    for(int i=0;i<adj[1][v].size();i++){
        auto u = adj[1][v][i];
        I.push_back(iH[u]);
    }
    sort(I.begin(), I.end(), std::greater<int>());

    int h = 0;
    for(int i = 0; i < I.size(); i++){
      if(I[i] >= i + 1)
        h = i + 1;
    }

    if (iH[v] > h) {
      iH[v] = h;
      //将out-neighbor设置为true
      for(int i=0;i<adj[0][v].size();i++){
        auto u = adj[0][v][i];
        changed[u] = true;
      }
      return true;
    }
    return false;
}

bool Graph::outHIndex(int v){
    std::vector<int> I;
    I.clear();
    for(int i=0;i<adj[0][v].size();i++){
        auto u = adj[0][v][i];
        I.push_back(oH[u]);
    }
    sort(I.begin(), I.end(), std::greater<int>());

    int h = 0;
    for(int i = 0; i < I.size(); i++){
      if(I[i] >= i + 1)
        h = i + 1;
    }

    if (oH[v] > h) {
      oH[v] = h;
      return true;
    }
    return false;
}


bool Graph::oHIndex(int v){
    bool update = false;

    changed[v] = false;

    for(int i=0;i<iHk[v].size();i++){
        std::vector<int> Ik;
        Ik.clear();

        for(int j=0;j<adj[0][v].size();j++){
            auto u = adj[0][v][j];
            if((iHk[u].size()-1) < i) continue;
            else
                Ik.push_back(iHk[u][i]);
        }
        sort(Ik.begin(), Ik.end(), std::greater<int>());

        int h = 0;
        for(int j = 0; j < Ik.size(); j++) {
            if(Ik[j] >= j + 1)
                h = j + 1;
            else
                break;
        }
        if (iHk[v][i] > h){
            iHk[v][i] = h;
            update = true;
            //in-neighbor set to true
            for(int i=0;i<adj[1][v].size();i++){
                auto u = adj[1][v][i];
                changed[u] = true;
            }
        }

    }

    return update;
}


bool Graph::Refine(int v){
    bool update = false;
    changed[v] = false;

    for(int i=0;i<iHk[v].size();i++){
        if(iHk[v][i] == 0){ // avoiding iHk[v][i] == 0 but running iHk[v][i] --
            continue;
        }
        int in = 0, out = 0;
        int estimate = iHk[v][i];

        //out-neighbor
        for(int j=0;j<adj[0][v].size();j++){
            auto u = adj[0][v][j];
            if((iHk[u].size()-1)>=i && iHk[u][i]>=estimate) out++;
        }

        //in-neighbor
        for(int j=0;j<adj[1][v].size();j++){
            auto u = adj[1][v][j];
            if((iHk[u].size()-1)>=i && iHk[u][i]>=estimate) in++;
        }

        if (in < i || out < estimate){
            iHk[v][i]--;
            update = true;
        }
    }

    return update;
}

bool Graph::outCore(int v){
    std::vector<int> I;
    I.clear();
    for(int i=0;i<adj[0][v].size();i++){
        auto u = adj[0][v][i];
        I.push_back(OC[u]);
    }
    sort(I.begin(), I.end(), std::greater<int>());

    int h = 0;
    for(int i = 0; i < I.size(); i++){
      if(I[i] >= i + 1)
        h = i + 1;
    }

    if (OC[v] > h) {
      OC[v] = h;
      return true;
    }
    return false;
}

bool Graph::Refine(int v, int k){

    changed[v] = false;
    if(iH[v] < k)
        return false;


    //当前的out-core number
    int estimate = OC[v]; // OC must greater than 0
    int pre_OC = OC[v];
    bool OC_flag = false; // denote if OC flag is changed while great than zero

    

    //=======================  out neighbor of v ========================================
	int greaterequals = 0;
	std::vector<int> smallers (estimate, 0);
	bool yep_set = false;

    for(int j=0;j < adj[0][v].size();j++){
        auto u = adj[0][v][j];
        int pu = 0;
        if(iH[u] >= k)
            pu = OC[u];
        if (pu >= estimate)
			greaterequals++;
		else
			smallers[pu]++;
        
        if (greaterequals == estimate) {
			yep_set = true;
			break;
		}
    }

    // estimate must greater than 0
    if (!yep_set && estimate > 0) { // watch for degree zeros
        OC_flag = true;

		int j;
		for (j = estimate - 1; j > 0; j--) {
			greaterequals += smallers[j];
			if (greaterequals >= j)
				break;
		}

		estimate = j;
        
        OC[v] = estimate;
	}
    
    if (!yep_set && estimate == 0) { // if estimate is 0, is must in (k,0)
        changed[v] = false;
		int j;

        //out-neighbor
        for(j=0;j<adj[0][v].size();j++){
            auto u = adj[0][v][j];
            if(iH[u]>=k && pre_OC >= OC[u] && OC[u] > estimate){
                changed[u] = true;
            }
        }

        //in-neighbor
        for(j=0;j<adj[1][v].size();j++){
            auto u = adj[1][v][j];
            if(iH[u]>=k && pre_OC >= OC[u] && OC[u] > estimate){
                changed[u] = true;
            }
        }
        OC[v] = estimate;

        adj[0][v].resize(0);
        adj[1][v].resize(0);

        // changed[v] = true; // don't need since OC[v] is consisted to 0
        return true;
	}

    //=======================  in neighbor of v ========================================
    std::vector<int> neigs(estimate,0);
    int geq = 0;
    bool yes = false; // 

    for(int j=0;j<adj[1][v].size();j++){
        auto u = adj[1][v][j];
        int pu = 0;
        if(iH[u] >= k)
            pu = OC[u];
        if (pu >= estimate)
			geq++;
		else
			neigs[pu]++;

        //满足条件还是之前的值
        if (geq == k) {
			yes = true;
			break;
		}
    }

    // if (!yes && OC[v] > 0) { // watch for degree zeros
    // deg[0][v] + deg[1][v] >= OC[v] >= estimate > 0
    if (!yes && estimate > 0) { // watch for degree zeros
        OC_flag = true;

		int j;
		for (j = estimate - 1; j > 0; j--) {
			geq += neigs[j];
			if (geq >= k)
				break;
		}

		estimate = j;

        OC[v] = estimate;
	}
    
    if (!yes && estimate == 0) {
        OC[v] = estimate;

		int j;

        //out-neighbor
        for(j=0;j<adj[0][v].size();j++){
            auto u = adj[0][v][j];
            if(iH[u]>=k && pre_OC >= OC[u] && OC[u] > estimate){
                changed[u] = true;
            }
        }
        adj[0][v].resize(0);

        //in-neighbor
        for(j=0;j<adj[1][v].size();j++){
            auto u = adj[1][v][j];
            if(iH[u]>=k && pre_OC >= OC[u] && OC[u] > estimate){
                changed[u] = true;
            }
        }
        adj[1][v].resize(0);

        // changed[v] = true; // don't need since OC[v] is consisted to 0
        return true;
	}

     if (OC_flag) { // watch for degree zeros
        changed[v] = true;
        int cur_idx = 0; // use to remove edge!!!
        //out-neighbor
        for(int j=0;j<adj[0][v].size();j++){
            auto u = adj[0][v][j];
            if(iH[u]>=k && OC[u]>=estimate){
                if(pre_OC >= OC[u] ){
                    changed[u] = true;
                }
                adj[0][v][cur_idx ++] = u;
            }
        }
        adj[0][v].resize(cur_idx);

        cur_idx = 0;
        //in-neighbor
        for(int j=0;j<adj[1][v].size();j++){
            auto u = adj[1][v][j];
            if(iH[u]>=k && OC[u]>=estimate){
                if(pre_OC >= OC[u] ){
                    changed[u] = true;
                }
                adj[1][v][cur_idx ++] = u;
            }
        }
        adj[1][v].resize(cur_idx);

        return true;
	} 



	return false;

}




bool Graph::BiSearch(std::vector<int> &NeighborSkyline, int &k, int &l){
    bool exist = false; 
    int len = NeighborSkyline.size() / 2 ;
    if(len > NeighborSkyline.size() / 2){
        int a = 1;
    }
    if(NeighborSkyline[len*2-2] < k){
        return exist;
    }
        
        
    int first = 0;
    int half , middle;

    while(len > 0){
        half = len >> 1 ;
        middle = first + half;

        if(NeighborSkyline[2*middle] >= k) 
            len = half;
        else{
            first = middle +1;
            len = len - half - 1;
        }
    }
    
    if(NeighborSkyline[2*first+1] >= l) 
        exist = true;
    else 
        exist = false;
            
    return exist;
}

void Graph::update(int v){


    int lmax = sH[v][1];
    int kmax = sH[v][sH[v].size()-2];
    tmp_results[v].resize(0);


    for(int k = 0; k <= kmax; k++){
        for(int l = lmax; l >= 0; l--){

            // printf("k = %d, l = %d. ****  ",k,l);
            
            int in = 0;
            int out = 0;

            //in-neighboe
            for(int j=0;j<adj[1][v].size();j++){
                auto u = adj[1][v][j];
                if(BiSearch(sH[u], k, l)){
                    in++; 
                }
                if(in >= k) break;
            }

            //out-neighbor
            for(int j=0;j<adj[0][v].size();j++){
                auto u = adj[0][v][j];
                if(BiSearch(sH[u], k, l)){
                    out++; 
                }
                if(out >= l) break;
            }

            
            //if skyline (k,l) exist
            if(in >= k && out >= l){
                // if(k < 0 || k >kmax){
                //     int a = 1;
                // }
                int VecSize = tmp_results[v].size();
                if(VecSize == 0 || l!=tmp_results[v][VecSize-1]){
                    tmp_results[v].push_back(k);
                    tmp_results[v].push_back(l);
                    lmax = l;
                    break;
                }
                else{
                    tmp_results[v][VecSize-2] = k;
                    break;
                }
            }
        }
    }

}

bool Graph::NoVerify(int v){
    bool update = false;
    if(tmp_results[v] != sH[v]){
        
        update = true;

        sH[v] = tmp_results[v];
    }

    return update;
}



void Graph::AC(){

    omp_set_num_threads(num_of_thread);

    /////////////////phase 1: compute kmax

    iH.resize(n,0);
    #pragma omp parallel for schedule(static)
    for(int i=0;i<n;i++)
        iH[i] = deg[1][i];

    bool flag = true;

    while(flag){
        flag = false;
        #pragma omp parallel for schedule(static)
        for(int i=0;i<n;i++){
            if(changed[i]){
                if(HIndex(i))
                    flag = true;
            }     
        }
    }


    // /////////////////phase 2: compute lupp

    //initialize iHk
    for(int i=0;i<n;i++){
        int dout = deg[0][i];
        int kv = iH[i];
        for(int k=0;k<=kv;k++)
            iHk[i].push_back(dout);
    }

    flag = true;

    changed.resize(n,true);

    while(flag){
        flag = false;
        #pragma omp parallel for schedule(dynamic,1000)
        for(int i=0;i<n;i++){
            if(changed[i]){
                if(oHIndex(i))
                    flag = true;
            }
        }
    }

    

    // /////////////////phase 3: refine lupp
    flag = true;

    while(flag){
        flag = false;
        #pragma omp parallel for schedule(dynamic,1000)
        for(int i=0;i<n;i++){
            if(Refine(i))
                flag = true;
        }
    }

    printf("anchored coreness done.\n");

}

void Graph::ACP(){
    //initialize iHindex, iH[v] = din[v], oH[v] = dout[v]
    omp_set_num_threads(num_of_thread);
    iH.resize(n,0);
    oH.resize(n,0);
#pragma omp parallel for schedule(static)
    for(int i=0;i<n;i++){
        oH[i] = deg[0][i];
        iH[i] = deg[1][i];
    }

    //compute kmax, lmax
    bool flag = true;
    while(flag){
        flag = false;
#pragma omp parallel for schedule(static)
        for(int i=0;i<n;i++){
            if(HIndex(i))
                flag = true;
            if(outHIndex(i))
                flag = true;
        }
    }

    //initialize iHk
#pragma omp parallel for schedule(static)
    for(int i=0;i<n;i++){
        int dout = oH[i];
        int kv = iH[i];
        for(int k=0;k<=kv;k++)
            iHk[i].push_back(dout);
    }

    // /////////////////phase 3: refine lupp
    flag = true;

    while(flag){
        flag = false;
#pragma omp parallel for schedule(static)
        for(int i=0;i<n;i++)
            if(Refine(i))
                flag = true;
    }

    printf("anchored coreness plus done.\n");



}

//parallel compute (k, 0)-core, dout = deg[0] = 0
std::vector<int> Graph::pkc(){
    //获取线程数
    int NUM_THREADS = omp_get_max_threads();
    // omp_set_num_threads(NUM_THREADS);

    int vis_num = 0;

    std::vector<int> Din;
    Din.resize(n);
#pragma omp for
    for (int i = 0; i < n; i++) {
        Din[i] = deg[1][i];
    }


    //每个线程运行这部分代码
#pragma omp parallel
{
    int level = 0;

    // long buff_size = n/NUM_THREADS + 1;
    long buff_size = n;

    int *buff = (int *)malloc(buff_size*sizeof(int));
    assert(buff != NULL);
    
    int start = 0, end = 0;

    while(vis_num < n){

        //获取in-degree为level的顶点
        #pragma omp for schedule(dynamic,1000)
        for(int i = 0; i < n; i++){
            // printf("dddd.\n");
            if(Din[i] == level){
                buff[end] = i;
                end++;
            }
        }

        //处理buff中的顶点
        while(start < end){
            int v = buff[start];
            start++;

            for(int j = 0; j < adj[0][v].size(); j++){
                int u = adj[0][v][j];
                int din_u = Din[u];

                if(din_u > level){
                    int du = __sync_fetch_and_sub(&Din[u], 1);
                    //将下一个level的顶点放在buff末尾
                    if(du==(level+1)){
                        buff[end] = u;
                        end++;
                    }

                    if(du <= level) __sync_fetch_and_add(&Din[u], 1);
                }

            }
        }

        //累计处理过的顶点数
        __sync_fetch_and_add(&vis_num, end);

        #pragma omp barrier
        start = 0;
        end = 0;
        level = level+1;

    }
    free(buff);
}

    return Din;

}

//parallel compute (l, 0)-core, din = deg[1] = 0
std::vector<int> Graph::plc(){
    int NUM_THREADS = omp_get_max_threads();

    int vis_num = 0;

    std::vector<int> Dout;
    Dout.resize(n);

#pragma omp for
    for (int i = 0; i < n; i++) {
        Dout[i] = deg[0][i];
    }


#pragma omp parallel
{
    int level = 0;

    long buff_size = n;

    int *buff = (int *)malloc(buff_size*sizeof(int));
    assert(buff != NULL);
    
    int start = 0, end = 0;

    while(vis_num < n){

        #pragma omp for schedule(dynamic,1000)
        for(int i = 0; i < n; i++){
            if(Dout[i] == level){
                buff[end] = i;
                end++;
            }
        }

        while(start < end){
            int v = buff[start];
            start++;

            for(int j = 0; j < adj[1][v].size(); j++){
                int u = adj[1][v][j];
                int dout_u = Dout[u];

                if(dout_u > level){
                    int du = __sync_fetch_and_sub(&Dout[u], 1);
                    if(du==(level+1)){
                        buff[end] = u;
                        end++;
                    }

                    if(du <= level) __sync_fetch_and_add(&Dout[u], 1);
                }

            }
        }

        __sync_fetch_and_add(&vis_num, end);

        #pragma omp barrier
        start = 0;
        end = 0;
        level = level+1;

    }
    free(buff);
}

    return Dout;
}


void Graph::ACP_plus(){


    int NUM_THREADS = num_of_thread;
    omp_set_num_threads(NUM_THREADS);


    auto begin = std::chrono::steady_clock::now();

    OC = plc();
    iH = pkc();

    auto end1 = std::chrono::steady_clock::now();
    

    //initialize iHk
#pragma omp parallel for schedule(dynamic,1000)
    for(int i=0;i<n;i++){
        int dout = OC[i];
        int kv = iH[i];
        for(int k = 0; k <= kv; k ++)
            iHk1[i].push_back(dout);
    }

    std::vector<int> res;
    auto begin_kshell = std::chrono::steady_clock::now();
#pragma omp parallel for
    for(int i = 0; i < changed.size(); i++){
        changed[i] = false;
    }
    
    res = kshell(iH);

    auto end_kshell = std::chrono::steady_clock::now();
    double runtime_kshell = std::chrono::duration<double>(end_kshell - begin_kshell).count();
    
    printf("number of shells: %ld.\n", res.size());
    
    ///////////////////phase 3: refine lupp for each k
    int num_it = 0;
    for(int i = 1; i < res.size(); i++){
        int k = res[i];
        int k_p = res[i-1];


#pragma omp parallel for
        for(int j = kBin[k_p]; j < kBin[k]; j++){                           
            OC[vert[j]] = 0;
            changed[vert[j]] = false;
        }


#pragma omp parallel for
        for(int j = kBin[k];j < vert.size(); j++){
            changed[vert[j]] = true;
        }


        bool flag = true;

        

        auto i_start = std::chrono::steady_clock::now();

        int chunk_size = (vert.size() - kBin[k]) / NUM_THREADS / 16 + 1;
        while(flag){

            flag = false;

#pragma omp parallel for schedule(dynamic, chunk_size)
            for(int j = kBin[k]; j < vert.size(); j++){
                int u = vert[j];

                //refine out-core number of u
                // if(changed[u] && OC[u] != 0){ // add OC[u] != 0 in order to remove unnecessary steps, but it is seemed not worth
                if(changed[u]){
                    if(Refine(u,k)){
                        flag = true;
                        // printf("u: %d.\n",u);
                    }
                }
            }
            num_it += 1;
        }

        auto i_end = std::chrono::steady_clock::now();
        double i_runtime = std::chrono::duration<double>(i_end - i_start).count();

        


#pragma omp parallel for
        for(int j = kBin[k]; j < vert.size(); j++){
            int u = vert[j];
            iHk1[u][k] = OC[u];
        }
    }

    auto end2 = std::chrono::steady_clock::now();

    double runtime1 = std::chrono::duration<double>(end1 - begin).count();
    double runtime2 = std::chrono::duration<double>(end2 - end1).count();
    printf("stage 1 running time: %.4f sec;\n", runtime1);
    printf("stage kshell running time: %.4f sec;\n", runtime_kshell);
    printf("stage 2 running time: %.4f sec;\n", runtime2);


    printf("anchored coreness plus+ done.\n");

}




void Graph::SC(){
    //set number of threads
    omp_set_num_threads(num_of_thread);

    //initialize iHindex, iH[v] = din[v], oH[v] = dout[v]
    iH.resize(n,0);
    oH.resize(n,0);
    #pragma omp parallel for schedule(dynamic,1000)
    for(int i=0;i<n;i++){
        oH[i] = deg[0][i];
        iH[i] = deg[1][i];
    }

    //compute kmax, lmax
    bool flag = true;
    while(flag){
        flag = false;
        #pragma omp parallel for schedule(dynamic,1000)
        for(int i=0;i<n;i++){
            if(changed[i]){
                if(HIndex(i))
                    flag = true;
            }
        }
    }

    changed.resize(n,true);
    flag = true;
    while(flag){
        flag = false;
        #pragma omp parallel for schedule(dynamic,1000)
        for(int i=0;i<n;i++){
            if(changed[i]){
                if(outHIndex(i))
                    flag = true;
            }
        }
    }

    //initial skyline core number
    for(int i=0;i<n;i++){
        int kmax = iH[i];
        int lmax = oH[i];
        sH[i].push_back(kmax);
        sH[i].push_back(lmax);
    }

    //update sc
    flag = true;
    while(flag){
        flag = false;
#pragma omp parallel for schedule(dynamic,1000)
        for(int i=0;i<n;i++){
            update(i);
        }

#pragma omp parallel for schedule(dynamic,1000)
        for(int i=0;i<n;i++){
            if(NoVerify(i))
                flag = true;
        }
        
    }

    printf("skyline coreness done.\n");
    
}


void Graph::Check(){
    bool correctness = true;

    for(int i=0;i<n;i++){
        int kv = iH[i];
        for(int k=0;k<=kv;k++)
            if(iHk[i][k] != iHk1[i][k]){
                correctness = false;
                printf("k: %d, AC: %d, ACP: %d.\n", k, iHk[i][k], iHk1[i][k]);
                break;
            }
        
        if(!correctness){
                printf("wrong!\n");
                break;
            }
    }



    if(correctness)
        printf("right!\n");

        
}

void Graph::Decom(){
    Kmd = max_deg[1];
    visited.resize(n,false);
    isInSub.resize(n,true);
    bin.resize(n+1);
    pos.resize(n+1);
    vert.resize(n);
    kBin.resize(n);
    Kdeg.resize(n);
    subLDegree.resize(n);
    subKDegree.resize(n);
    compulsoryVisit.resize(n);

    /**
    * description: calculate all nodes's Max K-core-number when there is no out degree constraint;
     * Namely, calculate all (k,0)-core for any possible k
    */
    for (int i = 0; i < n; i++) {
        Kdeg[i] = deg[1][i];
        kBin[deg[1][i]] += 1;
        subLDegree[i] = deg[0][i];
        subKDegree[i] = deg[1][i];
    }

    int start = 0;
    for (int i = 0; i <= Kmd; i++) {
        int num = kBin[i];
        kBin[i] = start;
        start += num;
    }
    for (int i = 0; i < n; i++) {
        pos[i] = kBin[Kdeg[i]];
        vert[pos[i]] = i;
        kBin[Kdeg[i]] += 1;
    }
    for (int i = Kmd; i > 0; i--) {
        kBin[i] = kBin[i - 1];
    }
    kBin[0] = 0;

    //decompose
    for (int i = 0; i < n; i++) {
        int v = vert[i];
        for(int j=0;j<adj[0][v].size();j++){
            int u = adj[0][v][j];         //Note there we should use out link to find (v,u) and reduce u's Kdeg
            if (Kdeg[u] > Kdeg[v]) {
                int du = Kdeg[u], pu = pos[u], pw = kBin[du], w = vert[pw];
                pos[u] = pw;    vert[pu] = w;
                pos[w] = pu;    vert[pw] = u;
                kBin[du] += 1;           //pull that higher-order vertex back, and the start position consequently +1
                Kdeg[u] -= 1;
            }
        }
    }
    kCoreMax = Kdeg[vert[n-1]];  //Kdeg[v], in-degree core number of v

    std::vector<int> row = inDependentRows();



    //for each k, decompose l
    std::vector<std::vector<int>> Dindex;
    std::vector<int> sortedID;

    for (auto k : row) {
        calculate_kl(k);
        Dindex.push_back(subLDegree);   
    }

}

/**
    * description: This function get rows' No. for those rows that are independent
     *   Namely, some vertex(es) exist in this row but not any row with higher row's No.
     *   This is used to do a prune operation, avoid calculating some rows which are exactly same, mentioned in CSD paper
*/
std::vector<int> Graph::inDependentRows(){
    int record = -1;
    int* result = new int [n];
    int resultPos = 0;
    for(int i = 0; i< vert.size();i++){
        if(Kdeg[vert[i]] != record){
            result[resultPos] = Kdeg[vert[i]];
            resultPos ++;
            record  = Kdeg[vert[i]];
        }
    }
    result[0] = -1;
    int maxK = record;
    for(int i = 1;i<resultPos;i++){
        result[i] = result[i-1]+1;
    }
    result[resultPos] = maxK;

    std::vector<int> res;
    for(int i=1;i<resultPos+1;i++){
        res.push_back(result[i]);
    }
    return res;
}

/**
    * description:  Based on a list (sorted by out degree) of the nodes in (K,0) subgraph, begin to process this K-th row
*/
void Graph::calculate_kl(int k) {
    rowResultPos = 0;
    rowVert.clear();
    for(int i=kBin[k];i<vert.size();i++)
        rowVert.push_back(vert[i]);
    // sorted_ID.resize(rowVert.size());

    for(int i = 0;i<= n;i++)  {
        bin[i] = 0;       //re- initialize bin for D-core use
        pos[i] = 0;
    }
    for(int i=0;i<n;i++){
        visited[i] = false;
        isInSub[i] = true;
    }

    for(int i = 0;i<kBin[k];i++){                           //record subgraph,vert
        isInSub[vert[i]] = false;
    }
    rowVert = sort_degree();
    for (auto value : rowVert) {  //decompose
        compulsoryNum = 0; compulsorySubL = -1;
        visit(k, value, -1);
        for (int q = 0; q < compulsoryNum; q++) {
            visited[compulsoryVisit[q]] = false;
            visit(k, compulsoryVisit[q], compulsorySubL);
        }
    }
    // return sorted_ID;
}

/**
    * description: Visit and report and delete
*/
void  Graph::visit(int k, int v, int compulsory) {
    int vid = v;
    if(visited[vid]) return;
    else{
        if(compulsory != -1) subLDegree[vid] = compulsory;
        compulsorySubL = subLDegree[vid];
        visited[vid] = true;
        // sorted_ID[ rowVert.size() - 1 - rowResultPos] = vid;
        rowResultPos++;
        
        for(int i = 0;i<adj[1][vid].size();i++) {             //v's inLink is its neighbor's outlink
            int neighV = adj[1][vid][i];
            if(!isInSub[neighV]) continue;
            if (!visited[neighV] && subLDegree[neighV] > subLDegree[vid]) {
                int du = subLDegree[neighV], pu = pos[neighV], pw = bin[du], w = rowVert[pw];
                pos[neighV] = pw;   rowVert[pw]= neighV;
                pos[w] = pu;    rowVert[pu] = w;
                subLDegree[neighV] -= 1;
                bin[du] += 1;
            }
        }
        
        for(int i = 0;i<adj[0][vid].size();i++) {
            int neighV = adj[0][vid][i];
            if(!isInSub[neighV]) {
                continue;
            }
            //Here maintain a compulsory queue, all vertices in this queue should be visited and deleted before we enter net visit()method.
            if(!visited[neighV]){
                if(subLDegree[vid] < subLDegree[neighV])
                    subKDegree[neighV] -= 1;
                if(subKDegree[neighV] < k) {
                    subKDegree[neighV]++;
                    compulsoryVisit[compulsoryNum] = neighV;
                    compulsoryNum++;
                    subLDegree[neighV] = subLDegree[vid];
                    visited[neighV] = true;
                }
            }
        }
    }
}

/**
    * description: When process the K-th row of D-core Table, this fucntion is used to sort all nodes by their l degree in (K,0) core SUBGRAPH
*/
std::vector<int> Graph::sort_degree(){
    std::vector<int> tempVert;
    tempVert.resize(rowVert.size());
    int subMaxL = 0;
    for(int i = 0;i<= n;i++)  bin[i] = 0;       //re- initialize bin for D-core use
    //calculate all vertex's subK&L Degree in current subgraph
    for (auto v : rowVert) {
        int subL = deg[0][v];
        int subK = deg[1][v];
        for (int l = 0; l < adj[0][v].size(); l++) {
            if (!isInSub[adj[0][v][l]]) subL--;
        }
        for (int k = 0; k < adj[1][v].size(); k++) {
            if (!isInSub[adj[1][v][k]]) subK--;
        }
        subLDegree[v] = subL;
        subKDegree[v] = subK;
        if (subL > subMaxL) subMaxL = subL;
        bin[subL] += 1;
    }
    int start = 0;
    for(int d = 0; d <= subMaxL;d ++){
        int num = bin[d];
        bin[d] = start;
        start += num;
    }
    for (auto v : rowVert) {   //fill the array tempVert with sorted vertex
        pos[v] = bin[subLDegree[v]];
        tempVert[pos[v]] = v;
        bin[subLDegree[v]] += 1;
    }

    for(int d = subMaxL ; d > 0; d--) bin[d] = bin[d - 1];

    bin[0] = 0;
    return tempVert;
}

void Graph::PDC_org(){
    omp_set_num_threads(num_of_thread);
    int NUM_THREADS = num_of_thread;


    int vis_num = 0;

    std::vector<int> Din, Dout;
    Din.resize(n);
    Dout.resize(n);
    for (int i = 0; i < n; i++) {
        Dout[i] = deg[0][i];
        Din[i] = deg[1][i];
    }

    auto begin = std::chrono::steady_clock::now();



    
#pragma omp parallel
{
    int level = 0;

    long buff_size = n;


    int *buff = (int *)malloc(buff_size*sizeof(int));
    assert(buff != NULL);
    
    int start = 0, end = 0;

    while(vis_num < n){

        #pragma omp for schedule(static)
        for(int i = 0; i < n; i++){
            if(Din[i] == level){
                buff[end] = i;
                end++;
            }
        }

        while(start < end){
            int v = buff[start];
            start++;

            for(int j = 0; j < adj[0][v].size(); j++){
                int u = adj[0][v][j];
                int din_u = Din[u];

                if(din_u > level){
                    int du = __sync_fetch_and_sub(&Din[u], 1);
                    if(du==(level+1)){
                        buff[end] = u;
                        end++;
                    }

                    if(du <= level) __sync_fetch_and_add(&Din[u], 1);
                }

            }
        }

        __sync_fetch_and_add(&vis_num, end);

        #pragma omp barrier
        start = 0;
        end = 0;
        level = level+1;
    }

    free(buff);
}
    
    int level = *max_element(Din.begin(),Din.end());

    
    std::vector<std::vector<int>> Fres;

    Fres.push_back(Parpeel_org(0));

    auto end1 = std::chrono::steady_clock::now();

    for (int k=1; k<=level;k++) {
        Fres.push_back(Parpeel_org(k));   
    }

    auto end2 = std::chrono::steady_clock::now();

    double runtime1 = std::chrono::duration<double>(end1 - begin).count();
    double runtime2 = std::chrono::duration<double>(end2 - end1).count();
    printf("stage 1 running time: %.4f sec, stage 2 running time: %.4lf sec.\n", runtime1, runtime2);

}


std::vector<int> Graph::Parpeel_org(int k){
    int NUM_THREADS = omp_get_max_threads();

    int visited = 0;

    
    std::vector<int> Din, Dout, flag;
    Din.resize(n);
    Dout.resize(n);
    flag.resize(n);
    for (int i = 0; i < n; i++) {
        Dout[i] = deg[0][i];
        Din[i] = deg[1][i];
        flag[i] = 0;
    }

#pragma omp parallel
{
    int level = 0;

    long buff_size = n;

    int *buff = (int *)malloc(buff_size*sizeof(int));
    assert(buff != NULL);
    
    int start = 0, end = 0;

    while(visited < n){
        #pragma omp for schedule(static)
        for(int i = 0; i < n; i++){
            if(flag[i] < 1){
                if(Dout[i] == level){
                    buff[end] = i;
                    end++;
                    flag[i]++;
                }
                else if (Din[i] < k){
                    Dout[i] = level;
                    buff[end] = i;
                    end++;
                    flag[i]++;
                }
            }
        }


        while(start < end){
            int v = buff[start];
            start++;
            flag[v]++;

            for(int j = 0; j < adj[0][v].size(); j++){
                int u = adj[0][v][j];

                if(flag[u]) continue;

                int din_u = __sync_fetch_and_sub(&Din[u], 1);

                if(din_u <= k){
                    if(flag[u] == 0){
                        buff[end] = u;
                        end++; 

                        Dout[u] = level;
                        flag[u]++;
                    }                
                }

            }

            for(int j = 0; j < adj[1][v].size(); j++){
                int u = adj[1][v][j];

                if(flag[u]) continue;

                int dout_u = Dout[u];

                if(dout_u > level){
                    int du = __sync_fetch_and_sub(&Dout[u], 1);
                    if(du==(level+1)){
                        buff[end] = u;
                        end++;
                        flag[u]++;
                    }

                    if(du <= level){
                        __sync_fetch_and_add(&Dout[u], 1);
                        flag[u]++;
                    } 
                }

            }
        }

        __sync_fetch_and_add(&visited, end);

        #pragma omp barrier
        start = 0;
        end = 0;
        level = level+1;
    }

    free( buff );
}

    return Dout;  
}


void Graph::PDC(){
    int NUM_THREADS = num_of_thread;
    omp_set_num_threads(NUM_THREADS);


    int vis_num = 0;

    std::vector<int> Din, Dout;
    Din.resize(n);
    Dout.resize(n);
    for (int i = 0; i < n; i++) {
        Dout[i] = deg[0][i];
        Din[i] = deg[1][i];
    }


#pragma omp parallel
{
    int level = 0;

    long buff_size = n;


    int *buff = (int *)malloc(buff_size*sizeof(int));
    assert(buff != NULL);
    
    int start = 0, end = 0;

    while(vis_num < n){


        #pragma omp for schedule(static)
        for(int i = 0; i < n; i++){
            if(Din[i] == level){
                buff[end] = i;
                end++;
            }
        }


        while(start < end){
            int v = buff[start];
            start++;

            for(int j = 0; j < adj[0][v].size(); j++){
                int u = adj[0][v][j];
                int din_u = Din[u];

                if(din_u > level){
                    int du = __sync_fetch_and_sub(&Din[u], 1);
                    if(du==(level+1)){
                        buff[end] = u;
                        end++;
                    }

                    if(du <= level) __sync_fetch_and_add(&Din[u], 1);
                }

            }
        }

        __sync_fetch_and_add(&vis_num, end);

        #pragma omp barrier
        start = 0;
        end = 0;
        level = level+1;
    }

    
    free(buff);
}

    int level = *max_element(Din.begin(),Din.end());

    std::vector<int> res;
    res = kshell(Din);
    std::vector<std::vector<int>> Fres;

    Fres.push_back(Parpeel(0));
    for(int i=1; i<res.size(); i++){
        Fres.push_back(Parpeel(res[i],Fres.back()));
    }


    printf("decomposition done.\n");


}

//obtain all distinct k
std::vector<int> Graph::kshell(std::vector<int> Din){
    pos.resize(n+1);
    vert.resize(n);
    kBin.resize(n);

    for (int i = 0; i < n; i++) {
        kBin[Din[i]] += 1;
    }

    int level = *max_element(Din.begin(),Din.end());

    int start = 0;
    for (int i = 0; i <= level; i++) {
        int num = kBin[i];
        kBin[i] = start;
        start += num;
    }
    for (int i = 0; i < n; i++) {
        pos[i] = kBin[Din[i]];
        vert[pos[i]] = i;
        kBin[Din[i]] += 1;
    }
    for (int i = level; i > 0; i--) {
        kBin[i] = kBin[i - 1];
    }
    kBin[0] = 0;

    

    int record = -1;
    int* result = new int [n];
    int resultPos = 0;
    for(int i = 0; i< vert.size();i++){
        if(Din[vert[i]] != record){
            result[resultPos] = Din[vert[i]];
            resultPos ++;
            record  = Din[vert[i]];
        }
    }
    result[0] = -1;
    int maxK = record;
    for(int i = 1;i<resultPos;i++){
        result[i] = result[i-1]+1;
    }
    result[resultPos] = maxK;

    std::vector<int> res;
    for(int i=1;i<resultPos+1;i++){
        res.push_back(result[i]);
    }

    return res;
}

std::vector<int> Graph::Parpeel0(){
    int NUM_THREADS = num_of_thread;
    omp_set_num_threads(NUM_THREADS);

    int visited = 0;


    std::vector<int> Din, Dout, flag;
    Din.resize(n);
    Dout.resize(n);
    flag.resize(n,0);

#pragma omp for
    for (int i = 0; i < n; i++) {
        Dout[i] = deg[0][i];
        Din[i] = deg[1][i];        
    }

#pragma omp parallel
{
    int level = 0;

    long buff_size = n/NUM_THREADS + 1;

    int *buff = (int *)malloc(n*sizeof(int));
    assert(buff != NULL);
    
    int start = 0, end = 0;

    while(visited < n){
        #pragma omp for schedule(static)
        for(int u = 0; u < n; u++){
            if(flag[u] < 1){
                if(Dout[u] == level){
                    buff[end] = u;
                    end++;
                    flag[u]++;
                }
            }
        }
        while(start < end){
            int v = buff[start];
            start++;
            flag[v]++;

            for(int j = 0; j < adj[1][v].size(); j++){
                int u = adj[1][v][j];

                if(flag[u]) continue;

                int dout_u = Dout[u];

                if(dout_u > level){
                    int du = __sync_fetch_and_sub(&Dout[u], 1);
                    if(du==(level+1)){
                        buff[end] = u;
                        end++;
                        flag[u]++;
                    }

                    if(du <= level){
                        __sync_fetch_and_add(&Dout[u], 1);
                        flag[u]++;
                    } 
                }

            }
        }

        __sync_fetch_and_add(&visited, end);

        #pragma omp barrier
        start = 0;
        end = 0;
        level = level+1;
    }

    free( buff );
}

    return Dout;  
}


//parallel peeling for a given k, advanced version
std::vector<int> Graph::Parpeel(int k){
    int NUM_THREADS = num_of_thread;

    #pragma omp parallel
    {
        #pragma omp master
        NUM_THREADS = omp_get_num_threads();
    }

    int visited = 0;



    std::vector<int> Din, Dout, flag, shellVert;
    Din.resize(n);
    Dout.resize(n);
    flag.resize(n,0);

    for (int i = 0; i < n; i++) {
        Dout[i] = deg[0][i];
        Din[i] = deg[1][i];        
    }
    
    for(int i=kBin[k];i<vert.size();i++)
        shellVert.push_back(vert[i]);

    for(int i = 0;i<kBin[k];i++){                           
        flag[vert[i]] = 1;
        Dout[vert[i]] = 0;
    }

    for(int i=kBin[k];i<vert.size();i++){
        int v = vert[i];
        for (int l = 0; l < adj[0][v].size(); l++) {
            if (flag[adj[0][v][l]]) Dout[v]--;
        }
        for (int k = 0; k < adj[1][v].size(); k++) {
            if (flag[adj[1][v][k]]) Din[v]--;
        }
    }

    int n1 = shellVert.size();

    

#pragma omp parallel
{
    int level = 0;

    long buff_size = n;

    int *buff = (int *)malloc(buff_size*sizeof(int));
    assert(buff != NULL);
    
    int start = 0, end = 0;

    int test = 0;


    while(visited < n1){
        #pragma omp for schedule(static)
        for(int i = 0; i < n1; i++){
            int u = shellVert[i];
            if(flag[u] < 1){
                if(Dout[u] == level){
                    buff[end] = u;
                    end++;
                    flag[u]++;
                }
                else if (Din[u] < k){
                    Dout[u] = level;
                    buff[end] = u;
                    end++;
                    flag[u]++;
                }
            }

        }


        while(start < end){
            int v = buff[start];
            start++;
            flag[v]++;

            for(int j = 0; j < adj[0][v].size(); j++){
                int u = adj[0][v][j];

                if(flag[u]) continue;

                int din_u = __sync_fetch_and_sub(&Din[u], 1);

                if(din_u <= k){
                    if(flag[u] == 0){
                        buff[end] = u;
                        end++; 

                        Dout[u] = level;
                        flag[u]++;
                    }                
                }

            }

            for(int j = 0; j < adj[1][v].size(); j++){
                int u = adj[1][v][j];

                if(flag[u]) continue;

                int dout_u = Dout[u];

                if(dout_u > level){
                    int du = __sync_fetch_and_sub(&Dout[u], 1);
                    if(du==(level+1)){
                        buff[end] = u;
                        end++;
                        flag[u]++;
                    }

                    if(du <= level){
                        __sync_fetch_and_add(&Dout[u], 1);
                        flag[u]++;
                    } 
                }

            }

        }

        __sync_fetch_and_add(&visited, end);



        #pragma omp barrier
        start = 0;
        end = 0;
        level = level+1;

    }

    free( buff );
}

    return Dout;  
}


//parallel peeling for a given k, advanced version
std::vector<int> Graph::Parpeel(int k, std::vector<int> upper){
    int NUM_THREADS = omp_get_max_threads();

    #pragma omp parallel
    {
        #pragma omp master
        NUM_THREADS = omp_get_num_threads();
    }

    int visited = 0;



    std::vector<int> Din, Dout, flag, shellVert;
    Din.resize(n);
    Dout.resize(n);
    flag.resize(n,0);

    for (int i = 0; i < n; i++) {
        Dout[i] = deg[0][i];
        Din[i] = deg[1][i];        
    }
    
    for(int i=kBin[k];i<vert.size();i++)
        shellVert.push_back(vert[i]);

    int lmax = 0;
    for(int i = 0;i<kBin[k];i++){                           
        flag[vert[i]] = 1;
        Dout[vert[i]] = 0;
        if(upper[vert[i]]>lmax) lmax = upper[vert[i]];
    }

    for(int i=kBin[k];i<vert.size();i++){
        int v = vert[i];
        for (int l = 0; l < adj[0][v].size(); l++) {
            if (flag[adj[0][v][l]]) Dout[v]--;
        }
        for (int k = 0; k < adj[1][v].size(); k++) {
            if (flag[adj[1][v][k]]) Din[v]--;
        }
    }

    int n1 = shellVert.size();


#pragma omp parallel
{
    int level = 0;

    long buff_size = n;

    int *buff = (int *)malloc(buff_size*sizeof(int));
    assert(buff != NULL);
    
    int start = 0, end = 0;


    while(visited < n1){
        #pragma omp for schedule(static)
        for(int i = 0; i < n1; i++){
            int u = shellVert[i];
            if(flag[u] < 1){
                if(Dout[u] == level || upper[u] == level){
                    Dout[u] = level;
                    buff[end] = u;
                    end++;
                    flag[u]++;
                }
                else if (Din[u] < k){
                    Dout[u] = level;
                    buff[end] = u;
                    end++;
                    flag[u]++;
                }
            }
        }

      
        while(start < end){
            int v = buff[start];
            start++;
            flag[v]++;

            for(int j = 0; j < adj[0][v].size(); j++){
                int u = adj[0][v][j];

                if(flag[u]) continue;

                int din_u = __sync_fetch_and_sub(&Din[u], 1);

                if(din_u <= k){
                    if(flag[u] == 0){
                        buff[end] = u;
                        end++; 

                        Dout[u] = level;
                        flag[u]++;
                    }                
                }

            }

            for(int j = 0; j < adj[1][v].size(); j++){
                int u = adj[1][v][j];

                if(flag[u]) continue;

                int dout_u = Dout[u];

                if(dout_u > level){
                    int du = __sync_fetch_and_sub(&Dout[u], 1);
                    if(du==(level+1)){
                        buff[end] = u;
                        end++;
                        flag[u]++;
                    }

                    if(du <= level){
                        __sync_fetch_and_add(&Dout[u], 1);
                        flag[u]++;
                    } 
                }

            }
        }

        __sync_fetch_and_add(&visited, end);



        #pragma omp barrier
        start = 0;
        end = 0;
        level = level+1;

    }

    free( buff );
}


    return Dout;  
}


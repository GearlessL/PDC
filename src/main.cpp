#include <iostream>
#include "Graph.h"
#include <string>
#include <chrono>
#include "cmdline.h"
#include <sys/time.h>
using namespace std;


// l: out, k: in;
// 0: out, 1: in; 
int main(int argc, char *argv[]) {
    // input parameters
    cmdline::parser a;
    a.add<string>("file", 'f', "filename", true, "");
    a.add<int>("a", 'a', "algorithm", true);
    a.add<int>("t", 't', "threads", true);
  
    a.parse_check(argc, argv);

    // Read graph
    string filepath = a.get<string>("file");
    int type = a.get<int>("a");
    int t = a.get<int>("t");
    FILE* dFile = fopen(filepath.c_str(), "r");

    // // // =============================================
    // string filepath = "./materials/em.txt";
    // string ds = "em";
    // FILE* dFile = fopen(filepath.c_str(), "r");

    // int type = 6;
    // int t = 32;
    // // // =============================================

    printf("file: %s; \t algorithm: %d; \t threads: %d\n", filepath.c_str(), type, t);

    clock_t io_begin = clock();

    Graph g = Graph(dFile,t);

    clock_t io_end = clock();
    double io_secs = double(io_end - io_begin) / CLOCKS_PER_SEC;
    printf("io time: %.4f.\n", io_secs);


    struct timeval begin, end;
    gettimeofday(&begin, NULL);
    switch (type)
    {
    case 1:
        printf("done.\n");
        g.Decom();
        break;
    case 2:
        printf("done.\n");
        g.AC();
        break;
    case 3:
        printf("done.\n");
        g.SC();
        break;
    case 4:
        printf("done.\n");
        g.PDC_org();
        break;
    case 5:
        printf("done.\n");
        g.PDC();
        break;
    case 6:
        printf("done.\n");
        g.ACP_plus();
        break;
    default:
        printf("no such algorithm.\n");
        break;
    }

    gettimeofday(&end, NULL);
    double runtime = (end.tv_sec - begin.tv_sec) + (double)(end.tv_usec - begin.tv_usec)/1000000.0;
    printf("running time: %.4f sec.\n", runtime);

    return 0;
}
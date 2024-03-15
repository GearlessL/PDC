# The code of algorithms
* Peeling[1]: the state-of-the-art sequential D-core decomposition algorithm;
* AC[2]: the distributed anchored coreness-based D-core decomposition algorithm, and we parallelize it by using multi-threads;
* SC[2]: the state-of-the-art distributed skyline coreness-based D-core decomposition algorithm, and we parallelize it by using multi-threads;
* ParPeel: our proposed parallel D-core decomposition algorithm, which is depicted in Algorithm 3;
* ParPeel-Prune: our proposed parallel D-core decomposition algorithm in Algorithm 3 with pruning strategy in Lemma 3 of [1];
* Shell-PDC: our proposed parallel D-core decomposition algorithm, listed in Algorithm 5.

*[1] Yixiang Fang, Zhongran Wang, Reynold Cheng, Hongzhi Wang, Jiafeng Hu. Effective and Efficient Community Search over Large Directed Graphs. IEEE Transactions on Knowledge and Data Engineering (TKDE), 31(11): 2093-2107, 2019.*

*[2] Xuankun Liao, Qing Liu, Jiaxin Jiang, Xin Huang, Jianliang Xu, Byron Choi. Distributed D-core decomposition over large directed graphs. Proceedings of the VLDB Endowment, 15(8): 1546–1558, 2022.*

# Compiling and Running
## Compiling the program
```
g++ -std=c++17 -fopenmp ${fileDirname}/*.cpp -o ${FILE_NAME}
```

For example:
```
g++ -std=c++17 -fopenmp ./src/*.cpp -o dcore
```


## Running the program:
```
./${FILE_NAME} -t ${THREAD} -f ${GRAPH_FILE} -a ${ALGORITHM}
```

For example:
```
./dcore \
-t 64 \
-f ./materials/em.txt \
-a 6 
```

For simplicity, the program can be executed using a script named ```./run.py```.




## Input format
* THREAD:
The number of threads.

* GRAPH_FILE:
The first line contains the number of nodes and the number of directed edges in the graph, formatted as follows.
```
${NODES} ${EDGE}
``` 
Each subsequent line represents a directed edge from node $u$ to node $v$, formatted as follows.
```
${u} ${v}
```


* ALGORITHM:
The ID of the algorithm to be executed is presented as follows.

| ID | Algorithms |
| :----: | :----: |
| 1 | Peeling |
| 2 | AC |
| 3 | SC |
| 4 | ParPeel |
| 5 | ParPeel-Prune |
| 6 | Shell-PDC |

# Other
* The ten datasets used in the paper are available from:


<div style="text-align:center">
    <table>
        <tbody>
        <tr>
            <th>Name</th>
            <th>Abbr.</th>
            <th>Source</th>
        </tr>
        <tr>
            <td>Email-EuAll</td>
            <td>EM</td>
            <td rowspan="5">https://snap.stanford.edu/data/index.htm</td>
        </tr>
        <tr>
            <td>Slashdot</td>
            <td>SD</td>
        </tr>
        <tr>
            <td>Amazon</td>
            <td>AM</td>
        </tr>
        <tr>
            <td>Pokec</td>
            <td>PO</td>
        </tr>
        <tr>
            <td>Live Journal</td>
            <td>LJ</td>
        </tr>
        <tr>
            <td>Enwiki-2013</td>
            <td>EW</td>
            <td rowspan="5">https://law.di.unimi.it/index.php</td>
        </tr>
        <tr>
            <td>Hollywood</td>
            <td>HW</td>
        </tr>
        <tr>
            <td>Webbase</td>
            <td>WB</td>
        </tr>
        <tr>
            <td>IT-2004</td>
            <td>IT</td>
        </tr>
        <tr>
            <td>UK-2007</td>
            <td>UK</td>
        </tr>
        </tbody>
    </table>
</div>

* An example of the Email-EuAll dataset is provided in ```./materials/```, and the running time of all algorithms can be found in ```./result/EM/```.

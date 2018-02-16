#include <fstream>
#include <list>
#include <vector>
#include <queue>
#include <functional>
#include <pthread.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <math.h>
#include <mpi.h>

using namespace std;
#define INF 0x3f3f3f3f

typedef pair<int, int> iPair;
int E, P, V;

class Graph
{
    int V;    // No. of vertices
    list< pair<int, int> > *adj;

public:
    Graph(int V);  // Constructor
    void addEdge(int u, int v, int w);
    void shortestPath(int s);
};

Graph* graph;
int* final_graph;

Graph::Graph(int V)
{
    this->V = V;
    adj = new list<iPair> [V];
}

void Graph::addEdge(int u, int v, int w)
{
    adj[u].push_back(make_pair(v, w));
    adj[v].push_back(make_pair(u, w));
}

void Graph::shortestPath(int src)
{
    priority_queue< iPair, vector <iPair> , greater<iPair> > pq;
    vector<int> dist(V, INF);

    pq.push(make_pair(0, src));
	dist[src] = 0;

    while (!pq.empty())
    {
		int u = pq.top().second;
        pq.pop();

        list< pair<int, int> >::iterator i;
        for (i = adj[u].begin(); i != adj[u].end(); ++i)
        {
            int v = (*i).first;
            int weight = (*i).second;

            if (dist[v] > dist[u] + weight)
            {
                dist[v] = dist[u] + weight;
                pq.push(make_pair(dist[v], v));
            }
		}
	}
	for (int i = 0; i < V; i++) {
        final_graph[src*V + i] = dist[i];
    }
}

int minDistance(int dist[], bool sptSet[])
{
   int min = INT_MAX, min_index;

   for (int v = 0; v < V; v++)
     if (sptSet[v] == false && dist[v] <= min)
         min = dist[v], min_index = v;

   return min_index;
}

void dijkstra(int src)
{
	//printf("DDD %d\n", src);
     int dist[V];     // The output array.  dist[i] will hold the shortest
                      // distance from src to i
     bool sptSet[V]; // sptSet[i] will true if vertex i is included in shortest
                     // path tree or shortest distance from src to i is finalized
     for (int i = 0; i < V; i++)
        dist[i] = INT_MAX, sptSet[i] = false;

     dist[src] = 0;
     for (int count = 0; count < V-1; count++)
     {
       int u = minDistance(dist, sptSet);
       sptSet[u] = true;
		
		for (int v = 0; v < V; v++)
         if (!sptSet[v] && final_graph[u*V + v] && dist[u] != INT_MAX
                                       && dist[u]+final_graph[u*V + v] < dist[v])
            dist[v] = dist[u] + final_graph[u*V + v];
     }
     for(int i = 0; i < V; i++) {
		//if(src == 0) printf("%d ", dist[i]);
        final_graph[src*V + i] = dist[i];
    }
}

typedef struct thread_input {
    int start, end;
} thread_input;

bool isSparse;

void* sssp(void* t_inp) {
	int start = ((thread_input*)t_inp)->start;
	int end = ((thread_input*)t_inp)->end;
	for(int i = start; i < end; i++) {
		if(isSparse) graph->shortestPath(i);
		else dijkstra(i);
	}
}


int main(int argc, char** argv)
{
    assert(argc == 4);
    fstream file(argv[1], fstream::in);
    P = atoi(argv[3]);
    file >> V >> E;
	isSparse = (E*(int)(ceil(log2((double)V))) < V*V*6/10);
	//printf("isSparse: %d\n", isSparse);

    // create the graph given in above fugure
	Graph g(V);
	graph = &g;
	
    final_graph = new int [V*V];
	for(int i = 0; i < V; i++)
        for(int j = 0; j < V; j++)
            final_graph[i*V + j] = 0;
    for(int i = 0; i < E; i++) {
        int a, b, w;
        file >> a >> b >> w;
        if(isSparse) g.addEdge(a, b, w);
		final_graph[a*V + b] = final_graph[b*V + a] = w;
    }
	
    file.close();
 
	int rank;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if(P > V) P = V;
	if(rank >= P) return 0;

	//Subgraph
    int p_sub_ver[P];
    for(int i = 0; i < P; i++) {
        if(i == 0) p_sub_ver[i] = 0;
        else {
            int p_sub_batch = (i-1 < V%P)? V/P+1 : V/P;
            p_sub_ver[i] = p_sub_batch + p_sub_ver[i-1];
        }
    }
    int inner_sub_ver_cnt = (rank != P-1)? p_sub_ver[rank+1] - p_sub_ver[rank] : V - p_sub_ver[rank];
    //Subgraph


	int t_num = (inner_sub_ver_cnt > 12)? 12 : inner_sub_ver_cnt;
	//if(rank == 3) printf("sub_ver: %d\n", inner_sub_ver_cnt);
	int t_start_ver[t_num];
	pthread_t threads[t_num];
    //thread_input t_input[T];
	for(int i = 0; i < t_num; i++) {
		if(i == 0) t_start_ver[i] = 0;
		int t_batch = (i < inner_sub_ver_cnt%t_num)? inner_sub_ver_cnt/t_num+1 : inner_sub_ver_cnt/t_num;
		thread_input* t_inp = new thread_input;
		t_inp->start = p_sub_ver[rank] + t_start_ver[i];
		t_inp->end = p_sub_ver[rank] + t_start_ver[i] + t_batch;
		if(i != t_num-1) t_start_ver[i+1] = t_start_ver[i] + t_batch;
        //int* u_id = new int(p_sub_ver[rank] + i);
		//if(rank == 3) printf("t:%d,  start:%d\n", i, t_start_ver[i]);
        pthread_create(&threads[i], NULL, sssp, (void*)t_inp);
    }
    for(int i = 0; i < t_num; i++) pthread_join(threads[i], NULL);

	int recv_cnts[P];
	int displs[P];
	for(int i = 0; i < P; i++) {
		if(i == P-1) recv_cnts[i] = V * (V - p_sub_ver[i]);
		else recv_cnts[i] = V * (p_sub_ver[i+1]-p_sub_ver[i]); 
		displs[i] = p_sub_ver[i]*V;
		//printf("%d: r %d, d %d\n", rank, recv_cnts[i], displs[i]);
	}
	//for(int i = 0; i < V; i++) printf("%d ", final_graph[p_sub_ver[rank]][i]);
	//printf("HIIIIIII\n");	

	int* terminated_graph;
	if(rank == 0) terminated_graph = new int[V*V];

	MPI_Gatherv(final_graph+p_sub_ver[rank]*V, inner_sub_ver_cnt*V, MPI_INT, terminated_graph, recv_cnts, displs, MPI_INT, 0, MPI_COMM_WORLD);
	
	if(rank == 0) {
		file.open(argv[2], fstream::out);
		for(int i = 0; i < V; i++) {
			for(int j = 0; j < V; j++)
				file << terminated_graph[i*V + j] << ' ';
			file << endl;
		}
	    file.close();
	}

	

    return 0;
}


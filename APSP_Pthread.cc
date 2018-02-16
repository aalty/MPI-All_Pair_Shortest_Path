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

using namespace std;
#define INF 0x3f3f3f3f

typedef pair<int, int> iPair;
int E, T, V;

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
int** final_graph;

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
        final_graph[src][i] = dist[i];
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
         if (!sptSet[v] && final_graph[u][v] && dist[u] != INT_MAX
                                       && dist[u]+final_graph[u][v] < dist[v])
            dist[v] = dist[u] + final_graph[u][v];
     }
     for(int i = 0; i < V; i++) {
        final_graph[src][i] = dist[i];
    }
}

typedef struct thread_input {
    int start, end;
	bool isSparse;
} thread_input;

//bool isSparse;

void* sssp(void* src) {
    thread_input input = *(thread_input*)src;
    for(int i = input.start; i < input.end ; i++) {
        if(input.isSparse) {graph->shortestPath(i);}
		else {dijkstra(i);}
    }
}


int main(int argc, char** argv)
{
    assert(argc == 4);
    fstream file(argv[1], fstream::in);
    T = atoi(argv[3]);
    file >> V >> E;
	bool isSparse = (E*(int)(ceil(log2((double)V))) < V*V*6/10);
	printf("isSparse: %d, %d, %d\n", E*(int)(ceil(log2((double)V))), V*V*6/10, isSparse);

    // create the graph given in above fugure
    Graph g(V);
    graph = &g;

    final_graph = new int* [V];
    for(int i = 0; i < V; i++) final_graph[i] = new int[V];
	for(int i = 0; i < V; i++)
        for(int j = 0; j < V; j++)
            final_graph[i][j] = 0;
    for(int i = 0; i < E; i++) {
        int a, b, w;
        file >> a >> b >> w;
        g.addEdge(a, b, w);
		final_graph[a][b] = final_graph[b][a] = w;
    }
    file.close();


	pthread_t threads[T];
    thread_input t_input[T];
    int t_sub_ver[T];
    for(int t = 0; t < T; t++) {
        if(t == 0) t_sub_ver[t] = 0;
        else {
            int t_batch = (t < V%T)? V/T+1 : V/T;
            t_sub_ver[t] = t_batch + t_sub_ver[t-1];
        }
    }
	for(int i = 0; i < T; i++) {
        t_input[i].start = t_sub_ver[i];
        t_input[i].end = (i == T-1)? V : t_sub_ver[i+1];
		t_input[i].isSparse = isSparse;
        pthread_create(&threads[i], NULL, sssp, (void*)&t_input[i]);
    }
    for(int i = 0; i < T; i++) pthread_join(threads[i], NULL);

    file.open(argv[2], fstream::out);
    for(int i = 0; i < V; i++) {
        for(int j = 0; j < V; j++)
            file << final_graph[i][j] << ' ';
        file << endl;
    }
    file.close();

    return 0;
}


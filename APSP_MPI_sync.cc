#include <fstream>
#include <stdio.h>
#include <assert.h>
#include <mpi.h>

#define INF 0x3f3f3f3f

using namespace std;

int E, P, V;

int main(int argc, char** argv)
{
	assert(argc == 4);
	fstream file(argv[1], fstream::in);
	P = atoi(argv[3]);
	file >> V >> E;
	assert(P == V);
	int** graph = new int* [V];
	int* index = new int [V];
	int* edge = new int [2*E];
	for(int i = 0; i < V; i++) graph[i] = new int[V];
	for(int i = 0; i < V; i++) {
		for(int j = 0; j < V; j++) {
			if(i == j) graph[i][j] = 0;
			else graph[i][j] = INF;
		}
	}
	for(int i = 0; i < E; i++) {
		int a, b, w;
		file >> a >> b >> w;
		graph[a][b] = graph[b][a] = w;
	}
	file.close();
	for(int i = 0, ind = 0, edg = 0; i < V; i++) {
		for(int j = 0; j < V; j++) {
			if(i!=j && graph[i][j]!=INF) {
				ind++;	
				edge[edg++] = j;
			}
		}
		index[i] = ind;
	}
	
	//for(int i = 0; i < V; i++) printf("%d ", index[i]); printf("\n");
	//for(int i = 0; i < 2*E; i++) printf("%d ", edge[i]); printf("\n");
	
	int rank;
	MPI_Comm Comm_graph;
	MPI_Init(&argc, &argv);
	MPI_Graph_create(MPI_COMM_WORLD, V, index, edge, 1, &Comm_graph);
	MPI_Comm_rank(Comm_graph, &rank);
	
	/*if(rank == 49) {
		for(int i = 0; i < V; i++) {
			for(int j = 0; j < V; j++) 
				printf("%d ", graph[i][j]);
			printf("\n");
		}
		int n; 
		MPI_Graph_neighbors_count(Comm_graph, rank, &n);
		int neighbors[n];
		MPI_Graph_neighbors(Comm_graph, rank, n, neighbors);
		for(int i = 0; i < n; i++) printf("%d ", neighbors[i]); printf("\n");
	}*/

	int relax = 1, num_nb;
	MPI_Graph_neighbors_count(Comm_graph, rank, &num_nb);
	int send_buf[V], recv_buf[V*num_nb], nbs[num_nb];
	MPI_Graph_neighbors(Comm_graph, rank, num_nb, nbs);
 	for(int i = 0; i < V; i++) send_buf[i] = graph[rank][i];
	while(relax) {
		relax = 0;
		MPI_Neighbor_allgather(send_buf, V, MPI_INT, recv_buf, V, MPI_INT, Comm_graph);
		//if(rank == 2) {
		for(int i = 0; i < num_nb; i++) {
			int n = nbs[i];
			for(int j = 0; j < V; j++) {
				if(j != rank) {
					if(send_buf[j] > send_buf[n] + recv_buf[i*V + j]) {
						send_buf[j] = send_buf[n] + recv_buf[i*V + j];
						relax = 1;
					}
				}
			}
		}
		int tmp = relax;
		MPI_Allreduce(&tmp, &relax, 1, MPI_INT, MPI_BOR, Comm_graph);
				//printf(" %d ", recv_buf[i*V + j]);
			//printf("\n");
	}
	int* final_graph;
	if(rank == 0) final_graph = new int [V*V];
	MPI_Gather(send_buf, V, MPI_INT, final_graph, V, MPI_INT, 0, Comm_graph);
	if(rank == 0) {
		file.open(argv[2], fstream::out);
		for(int i = 0; i < V; i++) {
			for(int j = 0; j < V; j++) 
				file << final_graph[i*V + j] << ' ';
			file << endl;
		}
		file.close();
	}





	return 0;
}

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
	assert(P==V);
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
	
	int rank, size;
	MPI_Comm Comm_graph;
	MPI_Init(&argc, &argv);
	MPI_Graph_create(MPI_COMM_WORLD, V, index, edge, 0, &Comm_graph);
	MPI_Comm_rank(Comm_graph, &rank);
	MPI_Comm_size(Comm_graph, &size);
	
	int relax = 1, num_nb;
	MPI_Graph_neighbors_count(Comm_graph, rank, &num_nb);
	int send_buf[V], recv_buf[V*num_nb], nbs[num_nb];
	MPI_Graph_neighbors(Comm_graph, rank, num_nb, nbs);
 	for(int i = 0; i < V; i++) send_buf[i] = graph[rank][i];

	MPI_Request req;
	MPI_Status status;
	for(int i = 0; i < num_nb; i++) MPI_Isend(send_buf, V, MPI_INT, nbs[i], 0, MPI_COMM_WORLD, &req);
	int token, start = 0, end = 0;
	int vertex_color = 0, probe_flag, recv_source, recv_tag;	

	while(true) {
		MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &probe_flag, &status);
		//if(probe_flag == 1) printf("%d\n", probe_flag);
		if(probe_flag == 1) {
			recv_source = status.MPI_SOURCE;
			recv_tag = status.MPI_TAG;	//0 for data, 1 for toke
			if(recv_tag == 0) {	//send data
				MPI_Irecv(recv_buf, V, MPI_INT, recv_source, recv_tag, MPI_COMM_WORLD, &req);
				relax = 0;
				for(int i = 0; i < V; i++) {
					if(i != rank) {
						if(send_buf[i] > send_buf[recv_source] + recv_buf[i]) {
							send_buf[i] = send_buf[recv_source] + recv_buf[i];
							relax = 1;
						}
					}
				}

				if(relax == 1) {
					for(int i = 0; i < num_nb; i++) MPI_Isend(send_buf, V, MPI_INT, nbs[i], 0, MPI_COMM_WORLD, &req);
					vertex_color = 1;
				}
				else if(rank == 0) start++;
			}
			else if(recv_tag == 1) { //send token
				MPI_Irecv(&token, 1, MPI_INT, recv_source, 1, MPI_COMM_WORLD, &req);
				if(token == 2) break;

				if(rank == 0) {
					//printf("%d\n", token_tmp);
					if(token == 0 && vertex_color == 0) {
						if(end >= 10) {
							token = 2;
							for(int i = 1; i < size; i++) MPI_Isend(&token, 1, MPI_INT, i, 1, MPI_COMM_WORLD, &req);
							break;
						}
						else end++;
					}
					else if(vertex_color == 1) {
							token = 1;
							vertex_color = 0;
							//end /= 3;
					}
					else {
						token = 0;
						//end /= 3;
					}

					MPI_Isend(&token, 1, MPI_INT, size-1, 1, MPI_COMM_WORLD, &req);
				}
				else {
					if(vertex_color == 1) {
						token = 1;
						vertex_color = 0;
					}
	
					MPI_Isend(&token, 1, MPI_INT, rank-1, 1, MPI_COMM_WORLD, &req);
				}
			}
		}		

		if(start==num_nb && rank==0){
			int token_start = vertex_color;
			//printf("token_start: %d\n", token_start);
			MPI_Isend(&token_start, 1, MPI_INT, 1, 1, MPI_COMM_WORLD, &req);
			start ++;
			vertex_color = 0;
		}
	}
	//printf("fuck!\n");
	int* final_graph;
	if(rank == 0) final_graph = new int [V*V];
	MPI_Gather(send_buf, V, MPI_INT, final_graph, V, MPI_INT, 0, MPI_COMM_WORLD);
	if(rank == 0) {
		file.open(argv[2], fstream::out);
		for(int i = 0; i < V; i++) {
			for(int j = 0; j < V; j++) 
				file << final_graph[i*V + j] << ' ';
			file << endl;
		}
		file.close();
	}

	free(index);
	free(edge);
	free(graph);
	if(rank == 0) free(final_graph);
	return 0;
}

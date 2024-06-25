#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <limits.h>
#include <omp.h>
#include <mpi.h>
#include <time.h>

#define INFINITY INT_MAX

int uniqueNodes=0;
int *kMinpaths;
int k = 3; // Find K shortest paths

int** distanceMatrix(char ***dataNode, int dataSize) 
{
	int** disMatrix = malloc((dataSize * 2) * sizeof(int*));
	for (int i = 0; i < dataSize * 2; i++)
	{
        disMatrix[i] = malloc((dataSize * 2) * sizeof(int));
    }
	
	char** tempArr = malloc((dataSize * 2) * sizeof(char**));
    for (int i = 0; i < dataSize * 2; i++) 
	{
            tempArr[i] = malloc(256 * sizeof(char));
    }

    int i = 0;
	while(i<dataSize) 
	{
        bool isUnique=true;
        int j = 0;
		while(j<uniqueNodes)
		{
			if(strcmp(dataNode[i][0],tempArr[j])==0)
			{
				isUnique=false;
                break;
			}		
            j++;	
		}
        if(isUnique)
		{
            strcpy(tempArr[uniqueNodes],dataNode[i][0]);
            uniqueNodes++;
        }
        i++;
	}

    i = 0;
	while(i<dataSize)
	{
        bool isUnique=true;
        int j = 0;

		while(j<uniqueNodes) 
		{
			if(strcmp(dataNode[i][1],tempArr[j])==0) 
			{
				isUnique=false;
                break;
			}
            j++;
		}
        if(isUnique) 
		{
            strcpy(tempArr[uniqueNodes],dataNode[i][1]);
            uniqueNodes++;
        }
        i++;
	}
    i = 0;
    while( i < uniqueNodes)
	{
        int j = 0;
        while(j < uniqueNodes)
		{
            disMatrix[i][j] = INFINITY;
            if(i == j)
			{
                disMatrix[i][j] = 0;
            }
            j++;
        }
        i++;
    }

    i = 0;
    while(i<uniqueNodes)		  //adjacency matrix
    {
        int k = 0;
        while(k < dataSize)
		{
            if(strcmp(tempArr[i], dataNode[k][0])==0)
			{
                int l = 0;
                while( l < uniqueNodes)
				{
                    if(strcmp(dataNode[k][1], tempArr[l]) == 0)
					{
                        disMatrix[i][l] = atoi(dataNode[k][2]);
                        disMatrix[l][i] = atoi(dataNode[k][2]);
                    }
                    l++;
                }

            }
            k++;
        }
        i++;
    }
   return disMatrix;
}

int dataSetSize(char* filename)
{
	int counter=0;
	FILE* myfile=fopen(filename,"r");
	if(myfile==NULL)
	{
		printf("Error while opening file\n");
		return -1;
	}
	else
	{
		char line[256]="\0";
		while(fgets(line,sizeof(line),myfile)!=NULL)
		{
			counter++;
		}
		counter--;	//to remove the content's header line
	} 
	return counter;
}

char*** readFile(char* filename,int dataSize)
{
	char ***dataNode = malloc((dataSize)*sizeof(char**));
    for(int i = 0; i < dataSize; i++)
	{
        dataNode[i] = malloc((4)*sizeof(char*));
    }

    for(int i = 0; i < dataSize; i++)
	{
        for(int j = 0; j < 4; j++)
		{
            dataNode[i][j] = malloc((256)*sizeof(char));
        }
    }

	//reading data from file
	FILE* myfile=fopen(filename,"r");
	if(myfile==NULL)
	{
		printf("Error while opening file\n");
		return dataNode;
	}
	else
	{
		char line[256]="";
		//to extract and discard the contents of the header row
		fgets(line,sizeof(line),myfile);
		//reading line by line from our csv file
		int i=0;
		while(fgets(line,sizeof(line),myfile)!=NULL)
		{
			sscanf(line, "%[^,],%[^,],%[^,],%[^\n]", dataNode[i][0], dataNode[i][1], dataNode[i][2], dataNode[i][3]);
            i++;
		}	
		fclose(myfile); 
	}
    return dataNode;
}

// Declare and initialize memoization array
int *memoizationArray;


int shortestPath(int source, int destination, int **adj, int V) 
{
    // Use the memoization array to store distances calculated so far
    int* dist = memoizationArray;

    // Initialize distances to infinity for all nodes
    for (int i = 0; i < V; ++i) {
        dist[i] = INFINITY;
    }
    // Start from the source node
    dist[source] = 0; // Distance from source to source is 0

    // Declare visited array
    bool visited[uniqueNodes];
    int i = 0;
    while(i < V) {
        visited[i] = false;
        i++;
    }
    visited[source] = true;

    // Initialize the queue for BFS
    int queue[uniqueNodes];
    int front = 0;
    int rear = 0;
    queue[rear++] = source;

    // Perform BFS
    for (; front != rear; )
    {
        int u = queue[front++];

        // Check all adjacent nodes of u
        int v = 0;
        while(v < V) 
        {
            if (adj[u][v] && adj[u][v] != INFINITY) 
            {
                int newDist = dist[u] + adj[u][v];
                if (!(newDist > dist[v]))
                { 
                    // Update distance only if new distance is shorter
                    dist[v] = newDist;
                    if (visited[v] != true) 
                    {
                        visited[v] = true;
                        queue[rear++] = v; // Enqueue v only if it hasn't been visited yet
                    }
                }
            }
            v++;
        }
    }
    return dist[destination];
}

void findKthPath(int **distanceMatrix, int numNodes, int start, int end, int k) 
{   
    int mainStart = start; // Store the initial start node
    int j = 0;
    #pragma omp parallel for
    while(j < numNodes) // Loop through each node
    {
        if (distanceMatrix[start][j] != INFINITY && j != end) // Check if there's a valid path from start to the current node and it's not the end node
        {
            int localStart = j; // Set the local start node
            int i = 0;
            #pragma omp parallel for
            while(i < numNodes) // Loop through each node
            {
                int pathCost = 0; // Initialize the path cost
                if (distanceMatrix[localStart][i] && distanceMatrix[localStart][i] != INFINITY) { // Check if there's a valid path from local start to the current node and it's not the main start node
                    if(i != mainStart){
                        int tempEnd = end;
                        end = i;

                        // Compute the cost of the first segment of the path
                        #pragma omp task
                        {
                            pathCost = pathCost + shortestPath(localStart, end, distanceMatrix, numNodes); // Calculate the cost of the first segment of the path
                        }

                        end = tempEnd;
                        tempEnd = localStart;

                        // Compute the cost of the second segment of the path
                        #pragma omp task
                        {
                            if (localStart != mainStart) 
                            { // Check if the local start node is not the main start node
                                int tempEnd = end; // Store the current end node
                                int tempStart = localStart; // Store the current local start node

                                end = localStart; // Set the end node to the local start node
                                localStart = mainStart; // Set the local start node to the main start node

                                pathCost += shortestPath(localStart, end, distanceMatrix, numNodes); // Calculate the cost of the second segment of the path
                                
                                end = tempEnd; // Restore the end node
                                localStart = tempStart; // Restore the local start node
                            }
                        }


                        localStart = i; // Set the local start node to the current node

                        // Compute the cost of the third segment of the path
                        #pragma omp task
                        {
                            if (i != end) // Check if the current node is not the end node
                            {
                                pathCost += shortestPath(localStart, end, distanceMatrix, numNodes); // Calculate the cost of the third segment of the path
                            }
                        }
                        #pragma omp taskwait // Wait for all tasks to complete

                        // Update the K shortest paths array
                        #pragma omp critical
                        {
                            for (int x = 0; x < k; x++)  // Loop to maintain the K shortest paths array
                            {
                                if (pathCost < kMinpaths[x]) {
                                    int replaced = kMinpaths[x];
                                    kMinpaths[x] = pathCost;
                                    int newIndex = x + 1;

                                    for (int y = k - 1; y > newIndex; y--) {
                                        kMinpaths[y] = kMinpaths[y - 1];
                                    }
                                    kMinpaths[newIndex] = replaced;
                                    break;
                                }
                            }
                        }
                    }
                }
                i++; // Move to the next node
            }
        }
        j++; // Move to the next node
    }
}





int main(int argc, char *argv[]) 
{
	 
    MPI_Init(&argc, &argv);
    int numProcesses, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	srand(time(NULL) + rank);
    	
	int ** disMatrix = NULL;
    if (rank == 0) 		// Root process initializes the distance matrix
	{    
		char filename[]="doctorwho.csv";
		int dataSize=dataSetSize(filename);		//calling ftn to find no. of lines of our data file
		char*** dataNode = readFile(filename, dataSize);	//reading actual data from file
        disMatrix = distanceMatrix(dataNode, dataSize);		//ftn to create and initialize the adjacency matrix
    }
  
    MPI_Bcast(&uniqueNodes, 1, MPI_INT, 0, MPI_COMM_WORLD);	  // Broadcast no. of unique cities from root process to all other processes

    //initializing the memoization array with infinity
	memoizationArray = (int *)malloc(uniqueNodes * sizeof(int));
	for (int i = 0; i < uniqueNodes; ++i) 
	{
    	memoizationArray[i] = INFINITY;
	}

    if (rank != 0) 
	{
        disMatrix = (int **)malloc(uniqueNodes * sizeof(int *));
        int i = 0;
        while(i < uniqueNodes) 
		{
            disMatrix[i] = (int *)malloc(uniqueNodes * sizeof(int));
            i++;
        }
    }

    int i = 0;
	while(i<uniqueNodes)
	{
    	MPI_Bcast(disMatrix[i], uniqueNodes, MPI_INT, 0, MPI_COMM_WORLD);	// Broadcast disMatrix from root process to all other processes
        i++;
    }

	

    kMinpaths = (int *)malloc(k * sizeof(int));
    i = 0; 
    while(i < k)
	{
        kMinpaths[i] = INFINITY;
        i++;
    }

    clock_t start, stop;
    float execution_time;
    int source = (rand() % uniqueNodes) + 1; 
    int destination = source;

    while (destination == source) 
	{
        destination = (rand() % uniqueNodes) + 1;
    }

    printf("{%d ,%d}\n",source ,destination);
    // Computing K minimum paths for this source and destination pair
	start = clock();
    findKthPath(disMatrix, uniqueNodes, source, destination, k);
	stop = clock();
    printf("Shortest K(%d) Paths: ", k);

    int j = 0;
    while(j < k){
        printf("%d ", kMinpaths[j]);
        j++;
    }
    
    printf("\n");
    execution_time = ((float) (stop - start));
    execution_time /= CLOCKS_PER_SEC;
    printf("Execution time: %f seconds\n\n", execution_time);
    MPI_Finalize();
    return 0;
}
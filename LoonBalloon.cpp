//*****************************************************************************
// loonBalloon.cpp
// Author: Kevin Peters
// Date: November 20, 2017
// Class: CIS 350, Bruce Elenbogen
// Description/Purpose: Project loon is a solution to give internet access to remote areas using high-altitude balloons.
// In order for the network to function properly all the balloons must be able to communicate with all the others.
// To save batteries, each station is programmed to only transmit to its two closest stations. In the case of a tie,
// thefirst choice is to choose the westernmost station and in case that is a tie, choose the southernmost
// station.So, to test if all balloons can talk to all other balloons we need to use depth-first search.

// Input:  int loons
//		   int xCoordinate, yCoordinate
// Output: string allConnected(yes or no if all balloons are connected to each other)

//*****************************************************************************

/********************INCLUDES**************************************************/

#include <iostream>
#include <algorithm>
#include <vector>
#include <math.h>

/********************DATA TYPES************************************************/

struct Vertex {
	Vertex(int, int, int);
	int vertexNum;
	int xCoord;
	int yCoord;
	int westernMostP;
	int southernMostP;
	int sndWesternMostP;
	int sndSouthernMostP;
	int firstClosest;
	double closestDist;
	int scndClosest;
	double scndClosestDist;

};
struct DFSearch {
	DFSearch(int, int);
	int vertex;
	int touch1;
	int touch2;
};

template <class Object>
class Matrix {
public:
	Matrix(int rows = 0, int cols = 0) : Array(rows)	//constructor with two parameters
	{
		for (int i = 0; i < rows; i++)
			Array[i].resize(cols);
	}
	void resize(int rows, int cols)
	{
		Array.resize(rows);
		for (int i = 0; i < rows; i++)
			Array[i].resize(cols);
	}
	const std::vector<Object> & operator[](int row) const
	{
		return Array[row];
	}
	std::vector<Object> & operator[](int row)
	{
		return Array[row];
	}
	//getter
	int numrows() const
	{
		return Array.size();
	}
	//getter
	int numcols() const
	{
		return numrows() ? Array[0].size() : 0;	//return the size, else return 0
	}
private:
	std::vector< std::vector<Object> > Array;	//two dimensional vector called Array
};
class Graph {
public:
	Graph(int);
	void fillGraph(int loons);
	void printGraph();
	void calcDistance();
	void fillAdjMatrix();
	void printAdjMatrix();
	void depthFirstSearch();
	void destruct();
private:
	std::vector<Vertex> graphPoints;
	std::vector<DFSearch> searchTable;
	std::vector<double> distance;
	Matrix<int> adjMatrix;
};

//******************************************************************
// Vertex(int x, int y, int vertex)
// Purpose: Vertex struct contructor.
// Pre: -20 < x < 20, -20 < y < 20, -1 < vertex < 1001
// Post: Creates a node Vertex with the coordinates passed in
//******************************************************************
Vertex::Vertex(int x, int y, int vertex) {
	xCoord = x;
	yCoord = y;
	vertexNum = vertex;
	firstClosest = -1;
	scndClosest = -1;
	closestDist = -1.0;
	scndClosestDist = -1.0;
}

//******************************************************************
// DFSearch(int myVertex, int touchOne)
// Purpose: DFSearch struct contructor.
// Pre: myVertex >= 0 && myVertex <= graphPoints.size(), 0 < touchOne < graphPoints.size() * 2 - 1
// Post: Creates a struct DFSearch with the vertex and touchOne passed in
//******************************************************************
DFSearch::DFSearch(int myVertex, int touchOne) {
	vertex = myVertex;
	touch1 = touchOne;
	touch2 = -1;
}

//******************************************************************
// Graph(int loons)
// Purpose: Graph class contructor.
// Pre: 0 < loons < 1001
// Post: Creates an instance of Graph class
// and resizes the adjacency matrix to the number of loons passed in
//******************************************************************
Graph::Graph(int loons) {
	adjMatrix.resize(loons, loons);
}

//******************************************************************
// fillGraph(int loons)
// Purpose: Populate graphPoints, a vector of vertices
// Pre: 0 < loons < 1001
// Post: Fills in the vector, graphPoints with all the x,y coordinates on the graph
//******************************************************************
void Graph::fillGraph(int loons) {
	int xCoordinate, yCoordinate;
	adjMatrix.resize(loons, loons);
	for (int i = 0; i < loons; i++) {
		std::cin >> xCoordinate;
		std::cin >> yCoordinate;
		Vertex thisPoint(xCoordinate, yCoordinate, i);
		graphPoints.push_back(thisPoint);
	}
}

//******************************************************************
// printGraph()
// Purpose: Testing purposes, this will point out coordinates and other attributes of each vertex
// Pre: None
// Post: Prints out the x,y coordinates and first and second closest vertexes to each vertex in graphPoints
//******************************************************************
void Graph::printGraph() {
	for (int i = 0; i < graphPoints.size(); i++) {
		std::cout << graphPoints[i].xCoord << ", " << graphPoints[i].yCoord << ", " << graphPoints[i].firstClosest << ", " << graphPoints[i].scndClosest << std::endl;
	}
}

//******************************************************************
// calcDistance()
// Purpose: This function is used to loop through each vertex
// nested in each vertex it will loop through each loop != to itself and find the first
// and second closest points to it and save as an attribute to the vertex
// Pre:  None
// Post: Loops through all graphPoints vector and calculates each vertex's closest two vertices
// and stores them in the respective attributes in the vertex struct
//******************************************************************
void Graph::calcDistance() {
	double epsilon = .0001;
	double edgeDistance;
	for (int i = 0; i < graphPoints.size(); i++) {
		edgeDistance = 0;
		for (int j = 0; j < graphPoints.size(); j++) {
			if (i != j) {
				edgeDistance = sqrt(pow((graphPoints[i].xCoord - graphPoints[j].xCoord),2) + pow((graphPoints[i].yCoord - graphPoints[j].yCoord),2));
				//if first element is empty, set it here
				if (graphPoints[i].firstClosest == -1) {
					graphPoints[i].firstClosest = j;
					graphPoints[i].closestDist = edgeDistance;
					graphPoints[i].southernMostP = graphPoints[j].yCoord;
					graphPoints[i].westernMostP = graphPoints[j].xCoord;
				}
				//if second element is empty, set it here
				else if (graphPoints[i].scndClosest == -1) {
					if (edgeDistance < graphPoints[i].closestDist) {
						graphPoints[i].scndClosest = graphPoints[i].firstClosest;
						graphPoints[i].scndClosestDist = graphPoints[i].closestDist;
						graphPoints[i].sndSouthernMostP = graphPoints[i].southernMostP;
						graphPoints[i].sndWesternMostP = graphPoints[i].westernMostP;
						graphPoints[i].firstClosest = j;
						graphPoints[i].closestDist = edgeDistance;
						graphPoints[i].southernMostP = graphPoints[j].yCoord;
						graphPoints[i].westernMostP = graphPoints[j].xCoord;
					}
					else if (fabs(edgeDistance - graphPoints[i].closestDist) <= epsilon * fabs(edgeDistance)) {
						//do edge, first tiebreaker
						if (graphPoints[j].xCoord==graphPoints[i].westernMostP) {
							if (graphPoints[j].yCoord < graphPoints[i].southernMostP) {
								//first = edge
								graphPoints[i].firstClosest = j;
								graphPoints[i].closestDist = edgeDistance;
								graphPoints[i].southernMostP = graphPoints[j].yCoord;
								graphPoints[i].westernMostP = graphPoints[j].xCoord;
							}
						}
						else if (graphPoints[j].xCoord < graphPoints[i].westernMostP) {
							graphPoints[i].scndClosest = graphPoints[i].firstClosest;
							graphPoints[i].scndClosestDist = graphPoints[i].closestDist;
							graphPoints[i].sndSouthernMostP = graphPoints[i].southernMostP;
							graphPoints[i].sndWesternMostP = graphPoints[i].westernMostP;
							graphPoints[i].firstClosest = j;
							graphPoints[i].closestDist = edgeDistance;
							graphPoints[i].southernMostP = graphPoints[j].yCoord;
							graphPoints[i].westernMostP = graphPoints[j].xCoord;
						}
					}
					else {
						graphPoints[i].scndClosest = j;
						graphPoints[i].scndClosestDist = edgeDistance;
						graphPoints[i].sndSouthernMostP = graphPoints[j].yCoord;
						graphPoints[i].sndWesternMostP = graphPoints[j].xCoord;
					}
				}
				//first two elements are already set
				else {
					//setting new closest point
					if (edgeDistance < graphPoints[i].closestDist) {
						//moving first closest to second closest
						if (graphPoints[i].closestDist == graphPoints[i].scndClosestDist) {
							//western most tie breaker
							if (graphPoints[i].westernMostP == graphPoints[i].sndWesternMostP) {
								//southern most tie breaker
								if (graphPoints[i].southernMostP < graphPoints[i].sndSouthernMostP) {
									graphPoints[i].scndClosest = graphPoints[i].firstClosest;
									graphPoints[i].scndClosestDist = graphPoints[i].closestDist;
									graphPoints[i].sndSouthernMostP = graphPoints[i].southernMostP;
									graphPoints[i].sndWesternMostP = graphPoints[i].westernMostP;
									graphPoints[i].firstClosest = j;
									graphPoints[i].closestDist = edgeDistance;
									graphPoints[i].southernMostP = graphPoints[j].yCoord;
									graphPoints[i].westernMostP = graphPoints[j].xCoord;
								}
								else {
									//replace the first closest but not second
									graphPoints[i].firstClosest = j;
									graphPoints[i].closestDist = edgeDistance;
									graphPoints[i].southernMostP = graphPoints[j].yCoord;
									graphPoints[i].westernMostP = graphPoints[j].xCoord;
								}
							}
							else if (graphPoints[i].westernMostP < graphPoints[i].sndWesternMostP) {
								graphPoints[i].scndClosest = graphPoints[i].firstClosest;
								graphPoints[i].scndClosestDist = graphPoints[i].closestDist;
								graphPoints[i].sndSouthernMostP = graphPoints[i].southernMostP;
								graphPoints[i].sndWesternMostP = graphPoints[i].westernMostP;
								graphPoints[i].firstClosest = j;
								graphPoints[i].closestDist = edgeDistance;
								graphPoints[i].southernMostP = graphPoints[j].yCoord;
								graphPoints[i].westernMostP = graphPoints[j].xCoord;
							}
							else {
								graphPoints[i].firstClosest = j;
								graphPoints[i].closestDist = edgeDistance;
								graphPoints[i].southernMostP = graphPoints[j].yCoord;
								graphPoints[i].westernMostP = graphPoints[j].xCoord;
							}
						}
						//replace send with first and first is now j's edge
						else {
							graphPoints[i].scndClosest = graphPoints[i].firstClosest;
							graphPoints[i].scndClosestDist = graphPoints[i].closestDist;
							graphPoints[i].sndSouthernMostP = graphPoints[i].southernMostP;
							graphPoints[i].sndWesternMostP = graphPoints[i].westernMostP;
							graphPoints[i].firstClosest = j;
							graphPoints[i].closestDist = edgeDistance;
							graphPoints[i].southernMostP = graphPoints[j].yCoord;
							graphPoints[i].westernMostP = graphPoints[j].xCoord;
						}
					}
					//if edge distance = closest distance
					else if (fabs(edgeDistance - graphPoints[i].closestDist) <= epsilon * fabs(edgeDistance)) {
						if (graphPoints[j].xCoord == graphPoints[i].westernMostP) {
							if (graphPoints[j].yCoord < graphPoints[i].southernMostP) {
								//moving first closest to second closest
								if (graphPoints[i].closestDist == graphPoints[i].scndClosest) {
									//western most tie breaker
									if (graphPoints[i].westernMostP == graphPoints[i].sndWesternMostP) {
										//southern most tie breaker
										if (graphPoints[i].southernMostP < graphPoints[i].sndSouthernMostP) {
											graphPoints[i].scndClosest = graphPoints[i].firstClosest;
											graphPoints[i].scndClosestDist = graphPoints[i].closestDist;
											graphPoints[i].sndSouthernMostP = graphPoints[i].southernMostP;
											graphPoints[i].sndWesternMostP = graphPoints[i].westernMostP;
											graphPoints[i].firstClosest = j;
											graphPoints[i].closestDist = edgeDistance;
											graphPoints[i].southernMostP = graphPoints[j].yCoord;
											graphPoints[i].westernMostP = graphPoints[j].xCoord;
										}
										else {
											//replace the first closest but not second
											graphPoints[i].firstClosest = j;
											graphPoints[i].closestDist = edgeDistance;
											graphPoints[i].southernMostP = graphPoints[j].yCoord;
											graphPoints[i].westernMostP = graphPoints[j].xCoord;
										}
									}
									else if (graphPoints[i].westernMostP < graphPoints[i].sndWesternMostP) {
										graphPoints[i].scndClosest = graphPoints[i].firstClosest;
										graphPoints[i].scndClosestDist = graphPoints[i].closestDist;
										graphPoints[i].sndSouthernMostP = graphPoints[i].southernMostP;
										graphPoints[i].sndWesternMostP = graphPoints[i].westernMostP;
										graphPoints[i].firstClosest = j;
										graphPoints[i].closestDist = edgeDistance;
										graphPoints[i].southernMostP = graphPoints[j].yCoord;
										graphPoints[i].westernMostP = graphPoints[j].xCoord;
									}
									else {
										graphPoints[i].firstClosest = j;
										graphPoints[i].closestDist = edgeDistance;
										graphPoints[i].southernMostP = graphPoints[j].yCoord;
										graphPoints[i].westernMostP = graphPoints[j].xCoord;
									}
								}
								//replace send with first and first is now j's edge
								else {
									graphPoints[i].scndClosest = graphPoints[i].firstClosest;
									graphPoints[i].scndClosestDist = graphPoints[i].closestDist;
									graphPoints[i].sndSouthernMostP = graphPoints[i].southernMostP;
									graphPoints[i].sndWesternMostP = graphPoints[i].westernMostP;
									graphPoints[i].firstClosest = j;
									graphPoints[i].closestDist = edgeDistance;
									graphPoints[i].southernMostP = graphPoints[j].yCoord;
									graphPoints[i].westernMostP = graphPoints[j].xCoord;
								}
							}
						}
						else if (graphPoints[j].xCoord < graphPoints[i].westernMostP) {
							//moving first closest to second closest
							if (graphPoints[i].closestDist == graphPoints[i].scndClosest) {
								//western most tie breaker
								if (graphPoints[i].westernMostP == graphPoints[i].sndWesternMostP) {
									//southern most tie breaker
									if (graphPoints[i].southernMostP < graphPoints[i].sndSouthernMostP) {
										graphPoints[i].scndClosest = graphPoints[i].firstClosest;
										graphPoints[i].scndClosestDist = graphPoints[i].closestDist;
										graphPoints[i].sndSouthernMostP = graphPoints[i].southernMostP;
										graphPoints[i].sndWesternMostP = graphPoints[i].westernMostP;
										graphPoints[i].firstClosest = j;
										graphPoints[i].closestDist = edgeDistance;
										graphPoints[i].southernMostP = graphPoints[j].yCoord;
										graphPoints[i].westernMostP = graphPoints[j].xCoord;
									}
									else {
										//replace the first closest but not second
										graphPoints[i].firstClosest = j;
										graphPoints[i].closestDist = edgeDistance;
										graphPoints[i].southernMostP = graphPoints[j].yCoord;
										graphPoints[i].westernMostP = graphPoints[j].xCoord;
									}
								}
								else if (graphPoints[i].westernMostP < graphPoints[i].sndWesternMostP) {
									graphPoints[i].scndClosest = graphPoints[i].firstClosest;
									graphPoints[i].scndClosestDist = graphPoints[i].closestDist;
									graphPoints[i].sndSouthernMostP = graphPoints[i].southernMostP;
									graphPoints[i].sndWesternMostP = graphPoints[i].westernMostP;
									graphPoints[i].firstClosest = j;
									graphPoints[i].closestDist = edgeDistance;
									graphPoints[i].southernMostP = graphPoints[j].yCoord;
									graphPoints[i].westernMostP = graphPoints[j].xCoord;
								}
								else {
									graphPoints[i].firstClosest = j;
									graphPoints[i].closestDist = edgeDistance;
									graphPoints[i].southernMostP = graphPoints[j].yCoord;
									graphPoints[i].westernMostP = graphPoints[j].xCoord;
								}
							}
							//replace send with first and first is now j's edge
							else {
								graphPoints[i].scndClosest = graphPoints[i].firstClosest;
								graphPoints[i].scndClosestDist = graphPoints[i].closestDist;
								graphPoints[i].sndSouthernMostP = graphPoints[i].southernMostP;
								graphPoints[i].sndWesternMostP = graphPoints[i].westernMostP;
								graphPoints[i].firstClosest = j;
								graphPoints[i].closestDist = edgeDistance;
								graphPoints[i].southernMostP = graphPoints[j].yCoord;
								graphPoints[i].westernMostP = graphPoints[j].xCoord;
							}
						}
						else {
							graphPoints[i].scndClosest = j;
							graphPoints[i].scndClosestDist = edgeDistance;
							graphPoints[i].sndSouthernMostP = graphPoints[j].yCoord;
							graphPoints[i].sndWesternMostP = graphPoints[j].xCoord;
						}
					}
					else if (edgeDistance < graphPoints[i].scndClosestDist) {
						graphPoints[i].scndClosest = j;
						graphPoints[i].scndClosestDist = edgeDistance;
						graphPoints[i].sndSouthernMostP = graphPoints[j].yCoord;
						graphPoints[i].sndWesternMostP = graphPoints[j].xCoord;
					}
					//if edge = second closest using double comparision
					else if (fabs(edgeDistance - graphPoints[i].scndClosestDist) <= epsilon * fabs(edgeDistance)) {
						//tiebreaker edge and second closest point
						if (graphPoints[j].xCoord == graphPoints[i].sndWesternMostP) {
							if (graphPoints[j].yCoord < graphPoints[i].sndSouthernMostP) {
								graphPoints[i].scndClosest = j;
								graphPoints[i].scndClosestDist = edgeDistance;
								graphPoints[i].sndSouthernMostP = graphPoints[j].yCoord;
								graphPoints[i].sndWesternMostP = graphPoints[j].xCoord;
							}
						}
						else if (graphPoints[j].xCoord < graphPoints[i].sndWesternMostP) {
							graphPoints[i].scndClosest = j;
							graphPoints[i].scndClosestDist = edgeDistance;
							graphPoints[i].sndSouthernMostP = graphPoints[j].yCoord;
							graphPoints[i].sndWesternMostP = graphPoints[j].xCoord;
						}
					}
					else {

					}
				}
			}
		}
	}
}

//******************************************************************
// fillAdjMatrix()
// Purpose: This function loops through the graph points
// and fills in any connection in the adjacency matrix
// Pre:  None
// Post: Loops through all graphPoints to populate
// a '1' anywhere a vertex connects to another vertex
//******************************************************************
void Graph::fillAdjMatrix() {
	//set the IDtable to the number of vertices
	//IDtable.resize(adj.numrows(), 0);
	for (int i = 0; i < graphPoints.size(); i++) {
		for (int j = 0; j < graphPoints.size(); j++) {
			if (graphPoints[i].firstClosest == j || graphPoints[i].scndClosest == j) {
				adjMatrix[i][j] = 1;
			}
		}
	}
	return;
}

//******************************************************************
// printAdjMatrix()
// Purpose: Testing purposes, prints out the adjacency matrix
// Pre:  None
// Post: Prints out the adjacency matrix to test if it was populated correctly
//******************************************************************
void Graph::printAdjMatrix() {
	for (int i = 0; i < graphPoints.size(); i++) {
		for (int j = 0; j < graphPoints.size(); j++) {
			std::cout << adjMatrix[i][j] << " ";
		}
	std::cout << std::endl;
	}
}

//******************************************************************
// depthFirstSearch()
// Purpose: This function references the adjacency matrix to populate the search table
// and determine if there is a strong connection to all points
// Pre:  None
// Post: Populates the searchTable vector with structs "DFSearch"
// to calculate strongly connected depth-first search and will
// output a yes or no depending on if every point can communicate with every other point
//******************************************************************
void Graph::depthFirstSearch() {
	//base case is when it goes back to 0 and fills in touch point 2 or gets to the end of the line and returns one line up, recursively go through the rows
	DFSearch firstLine(0, 1);
	searchTable.push_back(firstLine);
	//touch initializes to 2 because already pushed back vertex 0 with touch 1 above
	int touch = 2;
	int maxTouch = graphPoints.size() * 2;
	bool inSearch;
	bool inLine;
	int i = 0;
	int k;
	//loop while max touch not reached our the first vertex doesn't have another touch
	while (touch <= maxTouch && searchTable[0].touch2 == -1) {
		inSearch = false;
		inLine = false;
		int j = 0;
		//loop through a line of the adj matrix until it finds a connection
		do {
			if (adjMatrix[i][j] == 1 && inLine == false) {
				//check if [i][j] is already in the depth first search table
				inSearch = false;
				k = 0;
				while (inSearch == false && k < searchTable.size()) {
					if (searchTable[k].vertex == graphPoints[j].vertexNum) {
						inSearch = true;
					}
					k++;
				}if (inSearch != true) {
					DFSearch newLine(j, touch);
					searchTable.push_back(newLine);
					touch++;
					i = j;
					inLine = true;
				}inSearch = false;
			}
			j++;
		} while (j<graphPoints.size() && inLine == false);
		if (inLine == false) {
			//we got to the end of the line without a match, increment touch2
			bool touch2Found = false;
			for (int m = 0; m < searchTable.size(); m++) {
				if (searchTable[m].touch2 != -1 && touch2Found != true) {
					touch2Found = true;
					searchTable[m - 1].touch2 = touch;
					touch++;
					if (m > 1) {
						i = searchTable[m - 2].vertex;
					}
				}
			}if (touch2Found == false) {
				searchTable[searchTable.size() - 1].touch2 = touch;
				touch++;
				i = searchTable[searchTable.size() - 2].vertex;
			}
		}
	}
	if (searchTable[0].touch2 == maxTouch) {
		std::cout << "yes" << std::endl;
	}
	else {
		std::cout << "no" << std::endl;
	}
}

//******************************************************************
// destruct()
// Purpose: This function clears the vectors and adjacency matrix used to prepare for another test
// Pre:  None
// Post: Resets graphPoints and searchTable vectors
// and the adjacency matrix to run another test
//******************************************************************
void Graph::destruct() {
	graphPoints.clear();
	searchTable.clear();
	adjMatrix.resize(0, 0);
}
int main(){
	int loons;

	//use strongly connected components algorithm to see if they're all connected
	//-20=< Xi;Yi-<20, and 1=< N =< 1000
	//weakly connected
	//break ties with WMP or SMP 
	//tree or adj matrix

	do {
		std::cin >> loons;
		Graph myGraph(loons);
		if (loons > 0) {
			if (loons <= 3) {
				myGraph.fillGraph(loons);
				std::cout << "yes" << std::endl;
			}
			else {
				myGraph.fillGraph(loons);
				//myGraph.printGraph();
				myGraph.calcDistance();
				//myGraph.printGraph();
				myGraph.fillAdjMatrix();
				//myGraph.printAdjMatrix();
				myGraph.depthFirstSearch();
				myGraph.destruct();
				//myGraph.printGraph();
				//myGraph.printAdjMatrix();
			}
		}
	//The program is terminated with loons = 0
	} while (loons != 0);
	return 0;
}

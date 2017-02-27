// Samuel Stanton
// 02/14/2017
// samuel.stanton@ucdenver.edu
// Program to determine optimal matching of refugees and host cities.

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <cmath>
#include <random>
#include <set>
#include <iomanip>
//#include <boost/heap/fibonacci_heap.hpp>

using namespace std;

// Parameters for array sizes. SOURCENUM is the number of refugee source nodes,
// HOSTNUM is the number of host cities.
const int SOURCENUM = 6;

const int HOSTNUM = 302;

const int COUNTRYNUM = 28;

// Objective weights alpha:preference, beta:balance
double alpha = .0001;

double beta = .9999;


struct edge 
{
	double cost;
	
	double cap;
	
	int head;
	
	int tail;
};


struct shortest_path {
	vector<int> path_edge_idx;
	vector<double> nodeLabels;
	double path_cost;
	int errval;
};


vector<vector<string> > readIn(vector<vector<string> > vect, string filename);

// Functions to construct graphs and adjacency lists
vector<edge> create_network(double localpref[HOSTNUM][SOURCENUM], double localbalance[HOSTNUM], double localcapacity[HOSTNUM]);

vector<edge> create_residual(vector<edge> localgraph);

vector<vector<int> > makeAdjList(vector<edge> localgraph);


// Successive Shortest Path and Dijkstra
vector<double> ssp_mincost(vector<edge> localgraph, double localgrant_totals[SOURCENUM]);

shortest_path dijkstra(int start_node, int end_node, vector<edge> localgraph, vector< vector<int> > localadj_list);


double objFuncVal(vector<edge> localgraph, vector<double> localflow);

double gradient(vector<edge> localgraph, vector<double> localflow, int index);

double meanBurden(vector<edge> localgraph, vector<double> localflow);

double prefVal(vector<edge> localgraph, vector<double> localflow);

double balanceVal(vector<edge> localgraph, vector<double> localflow);


double linesearch(vector<edge> localgraph, vector<double> localdescent, vector<double> localflow, vector<double> localdirection);

double lnSrchFncA(vector<edge> localgraph, vector<double> localflow, int index);

double lnSrchFncB(vector<edge> localgraph, vector<double> localflow, vector<double> localdescent);

double lnSrchFncC(vector<edge> localgraph, vector<double> localdescent);


int main()
{
	vector<vector<string> > cities, asylumapps, asylumgranted;


	// Specify filepaths
	cities = readIn(cities, "CityPop200K.csv");

	asylumapps = readIn(asylumapps, "AsylumAppsFirstTime.csv");

	asylumgranted = readIn(asylumgranted, "AsylumDecisionsFirstTime.csv");


	//costs
	double pref[HOSTNUM][SOURCENUM] = {};
	
	double balance[HOSTNUM] = {};
	

	// edge capacity
	double capacity[HOSTNUM] = {};
	
	double country_pref[COUNTRYNUM][SOURCENUM] = {};
	
	double app_array[COUNTRYNUM][SOURCENUM] = {};
	
	double grant_array[COUNTRYNUM][SOURCENUM] = {};
	
	double app_totals[SOURCENUM] = {};
	
	double grant_totals[SOURCENUM] = {};


	// Extracts data from vector strings, places into arrays
	for (int j = 0; j < COUNTRYNUM; j++)
	{
		for (int i = 0; i < SOURCENUM; i++)
		{
			app_array[j][i] = atof(asylumapps[j + 1][i + 1].c_str());
			grant_array[j][i] = atof(asylumgranted[j + 1][i + 1].c_str());
			app_totals[i] = app_totals[i] + app_array[j][i];
			grant_totals[i] = grant_totals[i] + grant_array[j][i];
		}

	}
	
	
	for (int j = 0; j < COUNTRYNUM; j++)
	{
		for (int i = 0; i < SOURCENUM; i++) { country_pref[j][i] = app_array[j][i] / app_totals[i]; }
	}

	
	int skip = 1;
	for (int j = 1; j < cities.size() - 1; j++)
	{
		int tempid = atoi(cities[j][2].c_str());

		if (tempid == 0) 
		{ 
			skip = skip + 1;
			continue; 
		}

		for (int i = 0; i < SOURCENUM; i++)
		{
			pref[j-skip][i] = country_pref[tempid - 1][i];
		}

		capacity[j-skip] = atof(cities[j][1].c_str()) / 100;
	}

	
	// Create network
	vector<edge> graph = create_network(pref, balance, capacity);
	
	vector<edge> staticgraph = graph;

	
	// Create adjacency list
	vector<vector<int> > adj_list = makeAdjList(graph);

	
	// Set start time for runtime calculation
	clock_t startTime = clock();
	
	double runTime;
	

	// Initialize flow
	vector<double> flow = ssp_mincost(graph, grant_totals);

	
	// Begin main loop
	double tol = 0.001; 	// Define exit tolerance;

	
	for (int rep = 0; true; rep++) 
	{
		if (rep % 100 == 99) { cout << '.'; }
		
	
		// Update graph balance edge costs with gradient
		for (int b_edge = 0; b_edge < HOSTNUM; b_edge++) 
		{
			graph[b_edge + HOSTNUM*SOURCENUM].cost = gradient(graph, flow, b_edge);
		}

		
		// Solve min cost on updated graph to find improving direction
		vector<double> direction = ssp_mincost(graph, grant_totals);
		
		vector<double> descent;
		
		for (unsigned int i = 0; i < flow.size(); i++) { descent.push_back(direction[i] - flow[i]); }

		
		// Line search for optimal step size
		double stepsize = linesearch(staticgraph, descent, flow, direction);

		
		// Update flow
		vector<double> oldflow = flow;
		
		double tempsum = 0;
		
		for (unsigned int i = 0; i < flow.size(); i++) 
		{
			flow[i] = flow[i] + stepsize*(descent[i]);
			
			tempsum = tempsum + pow((oldflow[i] - flow[i]), 2);
		}

		
		// Check stopping criteria
		if (sqrt(tempsum) < tol) 
		{ 
			runTime = (clock() - startTime) / CLOCKS_PER_SEC;
			
			cout << "\nAlgorithm has converged in " << runTime << " seconds after " << rep << " iterations. \n";
			
			break; 
		}
	
		else if (rep > pow(10, 6)) 
		{ 
			cout << "\n Algorithm failed to converge. \n";
		
			exit(-1); 
		}
	}

	
	cout << "Preference score: " << prefVal(staticgraph, flow) << endl;
	
	cout << "Variance score: " << balanceVal(staticgraph, flow) << endl;
	
	cout << "Mean Burden: " << meanBurden(staticgraph, flow) << endl;

	
	cin.get();
	cin.get();
	
}

double gradient(vector<edge> localgraph, vector<double> localflow, int index) 
{
	double gradval = 0, temp1 = 0, temp2 = 0, R = (double)HOSTNUM;
	
	int offset = HOSTNUM * SOURCENUM;

	double mean = meanBurden(localgraph, localflow);


	for (int node1 = 0; node1 < HOSTNUM; node1++) 
	{
		temp2 = mean - ( localflow[node1 + offset] / localgraph[node1 + offset].cap );

		
		if (node1 == index)
		{
			temp1 = temp1 + (1 / (R * localgraph[node1 + offset].cap) - 1 / localgraph[node1 + offset].cap) * temp2;
		}
		
		
		else
		{
			temp1 = temp1 + (1 / (R * localgraph[node1 + offset].cap)) * temp2;
		}
	}
	
	
	gradval = beta * (2 / R) * temp1;
	
	return gradval;
}

double meanBurden(vector<edge> localgraph, vector<double> localflow) 
{
	int offset = HOSTNUM * SOURCENUM;
	
	double temp2 = 0;


	for (int node2 = 0; node2 < HOSTNUM; node2++) 
	{
		temp2 = temp2 + (localflow[node2 + offset] / localgraph[node2 + offset].cap);
	}

	
	return temp2 / (double) HOSTNUM;
}

double lnSrchFncA(vector<edge> localgraph, vector<double> localflow, int index) 
{ 
	return meanBurden(localgraph, localflow) - localflow[index + HOSTNUM*SOURCENUM] / localgraph[index + HOSTNUM*SOURCENUM].cap; 
}

double lnSrchFncB(vector<edge> localgraph, vector<double> localflow, vector<double> localdescent) 
{
	double temp = 0, temp2 = 0;


	for (int i = 0; i < HOSTNUM; i++) 
	{
		temp = temp + lnSrchFncA(localgraph, localdescent, i) * lnSrchFncA(localgraph, localflow, i);
	}
	
	temp = 2 * beta / (double) HOSTNUM * temp;
	

	for (int i = 0; i < HOSTNUM * SOURCENUM; i++) 
	{
		temp2 = temp2 + localgraph[i].cost * localdescent[i];
	}


	temp2 = alpha * temp2;
	
	return temp + temp2;
}

double lnSrchFncC(vector<edge> localgraph, vector<double> localdescent) 
{
	double temp = 0;
	
	
	for (int i = 0; i < HOSTNUM; i++) 
	{
		temp = temp + pow(lnSrchFncA(localgraph, localdescent, i), 2);
	}
	
	
	temp = 2 * beta / (double) HOSTNUM *temp;
	
	return temp;
}

double linesearch(vector<edge> localgraph, vector<double> localdescent, vector<double> localflow, vector<double> localdirection) 
{
	double lambda = 1, temp = 0;

	temp = -lnSrchFncB(localgraph, localflow, localdescent) / lnSrchFncC(localgraph, localdescent);

	
	if (temp > 0 && temp < 1) { lambda = temp; }

	
	else if (objFuncVal(localgraph, localdirection) > objFuncVal(localgraph, localflow)) { lambda = 0; }

	
	return lambda;
}

double objFuncVal(vector<edge> localgraph, vector<double> localflow) 
{
	double temp1 = 0, temp2 = 0;
	
	int offset = HOSTNUM * SOURCENUM;
	
	temp1 = alpha * prefVal(localgraph, localflow);
	
	temp2 = beta * balanceVal(localgraph, localflow);
	
	return temp1 + temp2;
}

double prefVal(vector<edge> localgraph, vector<double> localflow) 
{
	double temp = 0;
	
	int offset = HOSTNUM * SOURCENUM;
	
	
	for (int i = 0; i < offset; i++) 
	{
		temp = temp + localgraph[i].cost * localflow[i] / alpha;
	}
	
	
	return temp;
}

double balanceVal(vector<edge> localgraph, vector<double> localflow) 
{
	double variance = 0, mean = meanBurden(localgraph, localflow);
	
	int offset = HOSTNUM * SOURCENUM;
	
	
	for (int i = 0; i < HOSTNUM; i++) 
	{
		variance = variance + pow((mean - localflow[i + offset] / localgraph[i + offset].cap), 2);
	}
	
	
	variance = 1 / double (HOSTNUM) * variance;
	
	return variance;
}

vector<double> ssp_mincost(vector<edge> localgraph, double localgrant_totals[SOURCENUM]) 
{
	vector<edge> residualgraph = create_residual(localgraph);
	
	vector<vector<int> > localadj_list = makeAdjList(residualgraph);
	
	vector<double> localflow(residualgraph.size(), 0);  // Initialize flow
	
	vector<double> netflow(localgraph.size(), 0);
	
	vector<double> nodePotential(HOSTNUM + SOURCENUM + 1, 0);  // Initialize node potentials
	
	vector<double> startingEdgeCosts;

	
	for (unsigned int i = 0; i < residualgraph.size(); i++)
	{
		startingEdgeCosts.push_back(residualgraph[i].cost);
	}
	

	// Define excess, demand nodes
	vector<int> excessNodes;
	
	vector<double> excessAmt;
	
	int sinkNode = HOSTNUM + SOURCENUM;
	
	double sinkAmt = 0;


	for (int i = 0; i < SOURCENUM; i++)
	{
		excessNodes.push_back(i);
		
		excessAmt.push_back(localgrant_totals[i]);
		
		sinkAmt = sinkAmt - excessAmt[i];
	}
	
	
	double MAXCOST = 0; // Nominal cost for unreachable nodes
	
	
	for (unsigned int e = 0; e < localgraph.size(); e++) 
	{
		MAXCOST = MAXCOST + abs(localgraph[e].cost);
	}


	shortest_path augmentPath, tempPath;
	
	
	while (!excessNodes.empty()) 
	{
		int sourceNode = excessNodes[0];
		

		// Update node potentials with shortest path distances from sourceNode
		augmentPath = dijkstra(sourceNode, sinkNode, residualgraph, localadj_list);

		
		for (int loopNode = 0; loopNode <= SOURCENUM + HOSTNUM; loopNode++) 
		{
			double tempCost = augmentPath.nodeLabels[loopNode];
			
			
			if (tempCost > MAXCOST) 
			{
				nodePotential[loopNode] = nodePotential[loopNode] - MAXCOST;
			}
			
			else 
			{
				nodePotential[loopNode] = nodePotential[loopNode] - tempCost;
			}
		}

		
		// Compute flow augmentation amt
		double minPathCap = DBL_MAX;
		
		for (unsigned int pathEdge = 0; pathEdge < augmentPath.path_edge_idx.size(); pathEdge++) 
		{
			double temp = residualgraph[augmentPath.path_edge_idx[pathEdge]].cap;
			
			if (temp < minPathCap) { minPathCap = temp; }
		}

		
		double delta = min(excessAmt[sourceNode], minPathCap);

		// Augment flow, update capacities
		for (unsigned int pathEdge = 0; pathEdge < augmentPath.path_edge_idx.size(); pathEdge++) 
		{
			int tempIdx = augmentPath.path_edge_idx[pathEdge];
			
			localflow[tempIdx] = localflow[tempIdx] + delta;
			
			
			if (pathEdge % 2 == 0)
			{
				residualgraph[pathEdge].cap = residualgraph[pathEdge].cap - delta;
				
				residualgraph[pathEdge + 1].cap = residualgraph[pathEdge + 1].cap + delta;
			}
			
			else
			{
				residualgraph[pathEdge].cap = residualgraph[pathEdge].cap - delta;
				
				residualgraph[pathEdge - 1].cap = residualgraph[pathEdge - 1].cap + delta;
			}
		}

		// Check if excess node exhausted
		sinkAmt = sinkAmt + delta;
		
		excessAmt[sourceNode] = excessAmt[sourceNode] - delta;
		
		
		if (excessAmt[sourceNode] <= 0) 
		{ 
			excessNodes.erase( excessNodes.begin() ); 
		}

		// Update residual edge cost to new reduced costs
		for (unsigned int i = 0; i < residualgraph.size(); i++) 
		{
			residualgraph[i].cost = startingEdgeCosts[i] - nodePotential[residualgraph[i].head] + nodePotential[residualgraph[i].tail];
		}
	}


	// Compute net flow
	for (unsigned int i = 0; i < netflow.size(); i++) 
	{
		netflow[i] = localflow[2 * i] - localflow[2 * i + 1];
	}


	return netflow;
}

shortest_path dijkstra(int start_node, int end_node, vector<edge> localgraph, vector< vector<int> > localadj_list) 
{
	shortest_path s_path;
	
	s_path.errval = 0;

	if (start_node == end_node) 
	{
		s_path.path_cost = 0;
		
		return s_path;
	}


	vector<double> min_dist(HOSTNUM + SOURCENUM + 1, DBL_MAX);
	
	vector<int> predNode(HOSTNUM + SOURCENUM + 1, -1);
	
	vector<int> predEdge(HOSTNUM + SOURCENUM + 1, -1);

	double minCostBound = 0;
	
	for (unsigned int i = 0; i < localgraph.size(); i++) 
	{
		minCostBound = minCostBound - abs(localgraph[i].cost);
	}

	
	min_dist[start_node] = 0;
	
	set< pair<double, int> > active_vertices;
	
	active_vertices.insert({ 0, start_node });

	// Future work: Implementing Fibonacci heap
	//boost::heap::fibonacci_heap < pair<double, int> >;

	
	while (!active_vertices.empty()) 
	{
		int where = active_vertices.begin()->second;
		
		active_vertices.erase( active_vertices.begin() );
		
		vector<int> edgelist = localadj_list[where];
		
		if (edgelist[0] == -1) { continue; }
	
		for (unsigned int i = 0; i < edgelist.size(); i++) 
		{
			if (localgraph[edgelist[i]].cap <= 0) { continue; }
			
			int endpt = localgraph[edgelist[i]].tail;
			
			double ecost = localgraph[edgelist[i]].cost;
			
			if (min_dist[endpt] > min_dist[where] + ecost) 
			{
				active_vertices.erase({ min_dist[endpt], endpt });
				
				min_dist[endpt] = min_dist[where] + ecost;
				
				predNode[endpt] = where;
				
				predEdge[endpt] = edgelist[i];
				
				active_vertices.insert({ min_dist[endpt], endpt });

				if (min_dist[endpt] < minCostBound) 
				{
					cout << "Negative cycle detected \n";
					cout << where << ' ' << endpt << endl;
					exit(-2);
				}
			}
		}
	}


	// Return error code -1 if end_node not reachable
	if (predNode[end_node] == -1) 
	{ 
		s_path.errval = -1;
		
		return s_path;
	}
	
	
	s_path.path_cost = min_dist[end_node];
	
	
	// Construct vector of edge indexes in shortest path
	int currentNode = end_node;
	
	while (currentNode != start_node) 
	{
		s_path.path_edge_idx.push_back(predEdge[currentNode]);
		
		currentNode = predNode[currentNode];
	}
	
	s_path.nodeLabels = min_dist;
	
	
	return s_path;
}

vector<edge> create_network(double localpref[HOSTNUM][SOURCENUM], double localbalance[HOSTNUM], double localcapacity[HOSTNUM]) 
{	
	// Initialize graph as vector of edges
	vector<edge> localgraph;
	
	int count = 0;
	
	// Add preference edges
	for (int j = 0; j < HOSTNUM; j++) 
	{
		for (int i = 0; i < SOURCENUM; i++) 
		{
			localgraph.push_back(edge());
			
			localgraph[count].cost = alpha * -localpref[j][i];
			
			localgraph[count].cap = localcapacity[j];
			
			localgraph[count].head = i;
			
			localgraph[count].tail = j + SOURCENUM;
			
			count++;
		}
	}
	
	
	// Add balance edges
	for (int j = 0; j < HOSTNUM; j++) 
	{
		localgraph.push_back(edge());
		
		localgraph[count].cost = beta * localbalance[j];
		
		localgraph[count].cap = localcapacity[j];
		
		localgraph[count].head = j + SOURCENUM;
		
		localgraph[count].tail = HOSTNUM + SOURCENUM;
		
		count++;
	}
	
	
	return localgraph;
}

vector<vector<int> > makeAdjList(vector<edge> localgraph) 
{
	vector<vector<int> > localadj_list;
	
	
	for (int n = 0; n <= HOSTNUM + SOURCENUM; n++) 
	{
		localadj_list.push_back(vector<int>());
		
		localadj_list[n].push_back(-1);
		
		for (unsigned int e = 0; e < localgraph.size(); e++) 
		{
			if (localgraph[e].head == n) 
			{
				if (localadj_list[n][0] == -1) { localadj_list[n].erase(localadj_list[n].begin()); }
				
				localadj_list[n].push_back(e);
			}
		}
	}
	
	
	return localadj_list;
}

vector<edge> create_residual(vector<edge> localgraph) 
{
	unsigned int sz = localgraph.size();
	
	
	// Create residual edges at odd indexes
	for (unsigned int i = 0; i < sz; i++) 
	{
		localgraph.insert(localgraph.begin() + 2 * i + 1, edge());
		
		localgraph[2 * i + 1].cost = -localgraph[2 * i].cost;
		
		localgraph[2 * i + 1].cap = 0;
		
		localgraph[2 * i + 1].head = localgraph[2 * i].tail;
		
		localgraph[2 * i + 1].tail = localgraph[2 * i].head;
	}
	
	
	return localgraph;
}

vector<vector<string> > readIn(vector<vector<string> > vect, string filename)
{

	fstream file;
	
	file.open(filename.c_str(), ios::in);
	
	int i = 0, k=0, colnum;
	
	string temp, temp1, temp2, temp3;


	// detect the number of elements in a row
	getline(file, temp);
	
	colnum = std::count(temp.begin(), temp.end(), ',');
	
	
	// return to beginning of file
	file.seekg(0);

	
	// make sure there are no extra lines at the end of the file nor spaces!!!!!
	while (!file.eof())
	{
		// initial push into one dimension 
		vect.push_back(vector<string>());
		
		for (k = 0; k < colnum; k++)
		{
			getline(file, temp, ',');
			
			vect[i].push_back(temp);
		}

		getline(file, temp2);
		
		vect[i].push_back(temp2);

		i++;
	}


	return vect;
}


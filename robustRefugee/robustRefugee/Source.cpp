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

const double P_HAT_SCALE = 0.1;

// Define budget of uncertainty




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
vector<edge> create_network(double localpref[HOSTNUM][SOURCENUM], double localbalance[HOSTNUM], double localcapacity[HOSTNUM], double alpha, double beta);

vector<edge> create_residual(vector<edge> localgraph);

vector<vector<int> > makeAdjList(vector<edge> localgraph);


// Successive Shortest Path and Dijkstra
vector<double> ssp_mincost(vector<edge> localgraph, double localgrant_totals[SOURCENUM]);

shortest_path dijkstra(int start_node, int end_node, vector<edge> localgraph, vector< vector<int> > localadj_list);


double objFuncVal(vector<edge> localgraph, vector<double> localflow, double alpha, double beta, double rho);

double gradient(vector<edge> localgraph, vector<double> localflow, int index, double alpha, double beta, double mean);

double gradient2(vector<edge> localgraph, vector<double> localflow, int index, double alpha, double affineP, double rho, double normV);

double meanBurden(vector<edge> localgraph, vector<double> localflow);

double prefValRobust(vector<edge> localgraph, vector<double> localflow, double alpha, double beta, double rho);

double balanceVal(vector<edge> localgraph, vector<double> localflow);


double linesearch(vector<edge> localgraph, vector<double> localdescent, vector<double> localflow, vector<double> localdirection, double alpha, double beta, double rho);

double lnSrchFncA(vector<edge> localgraph, vector<double> localflow, int index);

double lnSrchFncB(vector<edge> localgraph, vector<double> localflow, vector<double> localdescent, double alpha, double beta);

double lnSrchFncC(vector<edge> localgraph, vector<double> localdescent, double alpha, double beta);


double binarySearch(vector<double> fixedTerm, double alpha, double rho);

vector<double> preComputefPrime(vector<edge> localgraph, vector<double> currentFlow, vector<double> descentDir);

double fPrime(double lambda, vector<double> fixedTerm, double alpha, double rho);



void solver(double alpha, double beta, double rho);


int main()
{
	double rho;
	
	// Main routine, writes data to output file

	//for (int i = 1; i < 31; i++)
	//{
	//	// Objective weights: alpha for preference, beta for balance

	//	rho = (double) i / 10 + 5;
	//	
	//	solver(0.5, 0.5, rho);
	//}

	rho = 1;
	
	solver(0.5, 0.5, rho);


	cin.get();
	cin.get();

}

void solver(double alpha, double beta, double rho)
{

	double affinePerturb = 0.1; // Perturbation factor

	double tol = 1;		// Define exit tolerance;
	
	vector<vector<string> > cities, asylumapps, asylumgranted;

	fstream objVals;
	fstream	appAssignments;

	objVals.open("RA-RobustobjValsE-00alpha50R1.txt", ios::app);

	appAssignments.open("RA-RobustassignmentsE-00alpha50R1.txt", ios::app);

	// Specify filepaths
	cities = readIn(cities, "CityPop200K.csv");

	asylumapps = readIn(asylumapps, "AsylumAppsFirstTime.csv");

	asylumgranted = readIn(asylumgranted, "AsylumDecisionsFirstTime.csv");


	int countryID = 0;

	vector<double> countryflow(COUNTRYNUM, 0);

	//costs
	double pref[HOSTNUM][SOURCENUM] = { 0 };

	double balance[HOSTNUM] = { 0 };


	// edge capacity
	double capacity[HOSTNUM] = { 0 };

	double capTotal = 0;

	double country_pref[COUNTRYNUM][SOURCENUM] = { 0 };

	double app_array[COUNTRYNUM][SOURCENUM] = { 0 };

	double grant_array[COUNTRYNUM][SOURCENUM] = { 0 };

	double app_totals[SOURCENUM] = { 0 };

	double grant_totals[SOURCENUM] = { 0 };

	double scaleTerm = 5 * pow(10, 6);




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


	vector<int> changeOfCountry(COUNTRYNUM - 1, 0);

	vector<string> countryName(COUNTRYNUM);

	vector<double> countryUrbPop(COUNTRYNUM, 0);

	int j = 0;
	int country = 0;

	while (j < HOSTNUM)
	{
		int tempid = atoi(cities[j + country + 1][2].c_str());

		if (tempid == 0)
		{
			countryName[country] = cities[j + country + 1][0];

			if (j != 0) { changeOfCountry[country - 1] = j; }

			country++;

			continue;
		}

		for (int i = 0; i < SOURCENUM; i++)
		{
			pref[j][i] = country_pref[country - 1][i] / scaleTerm;
		}

		double urbPop = atof(cities[j + country + 1][1].c_str());

		capacity[j] = urbPop / 100;

		capTotal = capTotal + capacity[j];

		countryUrbPop[country - 1] = countryUrbPop[country - 1] + urbPop;

		j++;
	}


	// Create network
	vector<edge> graph = create_network(pref, balance, capacity, alpha, beta);

	vector<edge> staticgraph = graph;


	// Create adjacency list
	vector<vector<int> > adj_list = makeAdjList(graph);


	// Set start time for runtime calculation
	clock_t startTime = clock();

	double runTime;


	// Initialize flow with perfectly balanced flow to improve solving time

	vector<double> flow(graph.size(), 0);

	int offset = HOSTNUM * SOURCENUM;

	for (int j = 0; j < HOSTNUM; j++)
	{
		for (int i = 1; i < SOURCENUM; i++)
		{
			flow[j + offset] = flow[j + offset] + app_totals[i] * capacity[j] / capTotal;
		}

	}

	for (int j = 0; j < HOSTNUM; j++)
	{
		for (int i = 0; i < offset; i++)
		{
			if (graph[i].tail == j)
			{
				flow[i] = app_totals[graph[i].head] * capacity[j] / capTotal;
			}
		}
	}


	// Begin main loop


	for (int rep = 0; true; rep++)
	{
		if (rep % 100 == 99) {std::cout << '.'; }


		// Update graph balance edge costs with gradient
		
		double normV = 0;
		for (int i = 0; i < HOSTNUM*SOURCENUM; i++)
		{
			normV = normV + pow(affinePerturb * staticgraph[i].cost * flow[i], 2);
		}

		normV = sqrt(normV);

		double mean = meanBurden(staticgraph, flow);
		

		for (int i = 0; i < HOSTNUM*SOURCENUM; i++)
		{
			graph[i].cost = gradient2(staticgraph, flow, i, alpha, affinePerturb, rho, normV);
		}

		
		for (int b_edge = 0; b_edge < HOSTNUM; b_edge++)
		{
			graph[b_edge + HOSTNUM*SOURCENUM].cost = gradient(staticgraph, flow, b_edge, alpha, beta, mean);
		}

		// Solve min cost on updated graph to find improving direction
		vector<double> direction = ssp_mincost(graph, grant_totals);

		vector<double> descent;

		for (unsigned int i = 0; i < flow.size(); i++) { descent.push_back(direction[i] - flow[i]); }


		// Binary Search for optimal stepsize

		vector<double> fixedTerm = preComputefPrime(staticgraph, flow, descent);
		
		double stepsize = binarySearch(fixedTerm, alpha, rho);
		
		// In the deterministic case we used this command for stepsize

		// double stepsize = linesearch(staticgraph, descent, flow, direction, alpha, beta);


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

			cout << "Algorithm has converged in " << runTime << " seconds after " << rep << " iterations. \n";

			break;
		}

		else if (rep > pow(10, 6))
		{
			cout << "\n Algorithm failed to converge. \n";

			exit(-1);
		}
	}

	setprecision(12);

	//Aggregate city flows back to country flows



	for (int i = 0; i < HOSTNUM; i++)
	{
		if (i == changeOfCountry[countryID])
		{
			countryID++;
		}

		countryflow[countryID] = countryflow[countryID] + flow[i + HOSTNUM*SOURCENUM];
	}

	//cout << "Country		Urban Pop.		Asylum Apps Assigned\n";
	for (int i = 0; i < COUNTRYNUM; i++)
	{
		appAssignments << countryName[i] << ',' << countryUrbPop[i] << ',' << round(countryflow[i]) << endl;
	}

	appAssignments << endl;

	appAssignments.close();

	cout << "rho " << rho << endl;

	cout << "Preference score: " << prefValRobust(staticgraph, flow, alpha, beta, rho) << endl;

	cout << "Variance score: " << balanceVal(staticgraph, flow) << endl;

	cout << "Objective function value " << objFuncVal(staticgraph, flow, alpha, beta, rho) << endl;

	objVals << rho << ',' << prefValRobust(staticgraph, flow, alpha, beta, rho) << ',' << balanceVal(staticgraph, flow) << ',' << objFuncVal(staticgraph, flow, alpha, beta, rho) << ',' << runTime << endl;

	objVals.close();

	return;
}


double gradient(vector<edge> localgraph, vector<double> localflow, int index, double alpha, double beta, double mean)
{
	double gradval = 0, temp1 = 0, temp2 = 0, R = (double)HOSTNUM;

	int offset = HOSTNUM * SOURCENUM;

	for (int node1 = 0; node1 < HOSTNUM; node1++)
	{
		temp2 = mean - (localflow[node1 + offset] / localgraph[node1 + offset].cap);


		if (node1 == index)
		{
			temp1 = temp1 + (1 / (R * localgraph[node1 + offset].cap) - 1 / localgraph[node1 + offset].cap) * temp2;
		}


		else
		{
			temp1 = temp1 + (1 / (R * localgraph[node1 + offset].cap)) * temp2;
		}
	}


	gradval = (1-alpha)*(2 / R) * temp1;

	return gradval;
}

double gradient2(vector<edge> localgraph, vector<double> localflow, int index, double alpha, double affineP, double rho, double normV)
{
	double gradval = alpha*(localgraph[index].cost + rho * pow(affineP * localgraph[index].cost, 2) * localflow[index] / normV);

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


	return temp2 / (double)HOSTNUM;
}


double lnSrchFncA(vector<edge> localgraph, vector<double> localflow, int index)
{
	return meanBurden(localgraph, localflow) - localflow[index + HOSTNUM*SOURCENUM] / localgraph[index + HOSTNUM*SOURCENUM].cap;
}

double lnSrchFncB(vector<edge> localgraph, vector<double> localflow, vector<double> localdescent, double alpha, double beta)
{
	double temp = 0, temp2 = 0;


	for (int i = 0; i < HOSTNUM; i++)
	{
		temp = temp + lnSrchFncA(localgraph, localdescent, i) * lnSrchFncA(localgraph, localflow, i);
	}

	temp = 2 * beta / (double)HOSTNUM * temp;


	for (int i = 0; i < HOSTNUM * SOURCENUM; i++)
	{
		temp2 = temp2 + localgraph[i].cost * localdescent[i];
	}


	temp2 = alpha * temp2;

	return temp + temp2;
}

double lnSrchFncC(vector<edge> localgraph, vector<double> localdescent, double alpha, double beta)
{
	double temp = 0;


	for (int i = 0; i < HOSTNUM; i++)
	{
		temp = temp + pow(lnSrchFncA(localgraph, localdescent, i), 2);
	}


	temp = 2 * beta / (double)HOSTNUM *temp;

	return temp;
}

double linesearch(vector<edge> localgraph, vector<double> localdescent, vector<double> localflow, vector<double> localdirection, double alpha, double beta, double rho)
{
	double lambda = 1, temp = 0;

	temp = -lnSrchFncB(localgraph, localflow, localdescent, alpha, beta) / lnSrchFncC(localgraph, localdescent, alpha, beta);


	if (temp > 0 && temp < 1) { lambda = temp; }


	else if (objFuncVal(localgraph, localdirection, alpha, beta, rho) > objFuncVal(localgraph, localflow, alpha, beta, rho)) { lambda = 0; }


	return lambda;
}


double binarySearch(vector<double> fixedTerm, double alpha, double rho)
{
	double upperBd = 1, lowerBd = 0, lambda = 0.5, exitTol = pow(10, -2);
	int count = 0;

		if (fPrime(upperBd, fixedTerm, alpha, rho) < 0)
		{
			lambda = upperBd;
		}
		
		else if (fPrime(lowerBd, fixedTerm, alpha, rho) > 0)
		{
			lambda = lowerBd;
		}

		while (abs(fPrime(lambda, fixedTerm, alpha, rho)) >= exitTol || count > 99)
		{
		if (fPrime(lambda, fixedTerm, alpha, rho) > 0)
		{
			upperBd = lambda;
		}
		
		else if (fPrime(lambda, fixedTerm, alpha, rho) < 0)
		{
			lowerBd = lambda;
		}

		lambda = (lowerBd + upperBd) / 2;

		count++;
	}

	return lambda;
}

vector<double> preComputefPrime(vector<edge> localgraph, vector<double> currentFlow, vector<double> descentDir)
{
	vector<double> fixedTerm(6, 0);
	
	vector<double> pHat(HOSTNUM*SOURCENUM, 0);

	double meanDescentB = meanBurden(localgraph, descentDir);

	double meanFlowB = meanBurden(localgraph, currentFlow);


	double temp = 0;

	for (int i = 0; i < HOSTNUM*SOURCENUM; i++)
	{
		pHat[i] = localgraph[i].cost * P_HAT_SCALE;
		
		// Terms in robust p(x) without lambda dependence
		
		fixedTerm[0] = fixedTerm[0] + localgraph[i].cost * descentDir[i];

		fixedTerm[1] = fixedTerm[1] + pow(pHat[i], 2) * descentDir[i] * currentFlow[i];

		fixedTerm[2] = fixedTerm[2] + pow(pHat[i], 2) * pow(descentDir[i], 2);

		fixedTerm[3] = fixedTerm[3] + pow(pHat[i], 2) * pow(currentFlow[i], 2);


	}
	
	int offset = HOSTNUM * SOURCENUM;

	for (int i = offset; i < offset + HOSTNUM; i++)
	{
		// THETA FROM BALANCE DERIVATIVE

		fixedTerm[4] = fixedTerm[4] + (meanDescentB - descentDir[i] / localgraph[i].cap)*(meanFlowB - currentFlow[i] / localgraph[i].cap);

		// PHI FROM BALANCE DERIVATIVE
		
		fixedTerm[5] = fixedTerm[5] + pow(meanDescentB - descentDir[i] / localgraph[i].cap, 2);
	}

	fixedTerm[4] = 2 / (double) HOSTNUM * fixedTerm[4];

	fixedTerm[5] = 2 / (double)HOSTNUM * fixedTerm[5];

	return fixedTerm;
}

double fPrime(double lambda, vector<double> fixedTerm, double alpha, double rho)
{
	double primeVal = 0;

	
	// Calculate p'(x)

	primeVal = fixedTerm[0] + rho * (fixedTerm[1] + lambda * fixedTerm[2]) / sqrt(fixedTerm[3] + 2 * fixedTerm[1] * lambda + fixedTerm[2] * pow(lambda, 2));

	primeVal = alpha * primeVal;

	
	// f'(x) = a * p'(x) + (1-a) b'(x)
	
	primeVal = primeVal + (1 - alpha) * (fixedTerm[4] + lambda * fixedTerm[5]);
	
	
	return primeVal;
}


double objFuncVal(vector<edge> localgraph, vector<double> localflow, double alpha, double beta, double rho)
{
	double temp1 = 0, temp2 = 0;

	int offset = HOSTNUM * SOURCENUM;

	temp1 = alpha * prefValRobust(localgraph, localflow, alpha, beta, rho);

	temp2 = beta * balanceVal(localgraph, localflow);

	return temp1 + temp2;
}

double prefValRobust(vector<edge> localgraph, vector<double> localflow, double alpha, double beta, double rho)
{
	double temp = 0;

	int offset = HOSTNUM * SOURCENUM;

	double normVal = 0;

	for (int i = 0; i < offset; i++)
	{
		double cost = localgraph[i].cost, flowVal = localflow[i];

		temp = temp + cost * flowVal;

		normVal = normVal + pow(cost*P_HAT_SCALE * flowVal, 2);
	}

	normVal = rho * sqrt(normVal);

	return temp + normVal;
}

double balanceVal(vector<edge> localgraph, vector<double> localflow)
{
	double variance = 0, mean = meanBurden(localgraph, localflow);

	int offset = HOSTNUM * SOURCENUM;


	for (int i = 0; i < HOSTNUM; i++)
	{
		variance = variance + pow((mean - localflow[i + offset] / localgraph[i + offset].cap), 2);
	}


	variance = 1 / double(HOSTNUM) * variance;

	return variance;
}


vector<double> ssp_mincost(vector<edge> localgraph, double localapptotals[SOURCENUM])
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

		excessAmt.push_back(localapptotals[i]);

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
			excessNodes.erase(excessNodes.begin());
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

		active_vertices.erase(active_vertices.begin());

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


vector<edge> create_network(double localpref[HOSTNUM][SOURCENUM], double localbalance[HOSTNUM], double localcapacity[HOSTNUM], double alpha, double beta)
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

			localgraph[count].cost = -localpref[j][i];

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

		localgraph[count].cost = localbalance[j];

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

	if (!file.good())
	{
		cout << filename << " not found.\n";

		exit(-1);
	}

	int i = 0, k = 0, colnum;

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


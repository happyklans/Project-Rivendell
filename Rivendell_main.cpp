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
#include "stdafx.h"


using namespace std;

// Parameters for array sizes. SOURCENUM is the number of refugee source nodes,
// HOSTNUM is the number of host cities.
const int SOURCENUM = 6;
const int HOSTNUM = 302;

struct edge {
	double cost;
	double cap;
	int head;
	int tail;
};

class edge_data //Class for extracting data from the elements of graph
{
public:
	edge_data(int, int, double, double);
	void update_data(int, int, double, double);
	int head_val() { return head; };
	vector<int> return_position_values();
	vector<double> return_contextual_values();
private:
	double cost, cap;
	int head, tail;
};

edge_data::edge_data(int i1, int i2, double d1, double d2) //default constructor
{
	cost = d1;
	cap = d2;
	head = i1;
	tail = i2;
}
void edge_data::update_data(int newi1, int newi2, double newd1, double newd2) // updates info for iteration
{
	cost = newd1;
	cap = newd2;
	head = newi1;
	tail = newi2;
}
vector<int> edge_data::return_position_values()
{
	std::vector<int> position_values(2);
	position_values[0] = head;
	position_values[1] = tail;
	return position_values;
}
vector<double> edge_data::return_contextual_values()
{
	std::vector<double> contextual_values(2);

	contextual_values[0] = cost;
	contextual_values[1] = cap;
	return contextual_values;
}

vector<vector<int> > populate_adj_lst(vector<edge> graph, int node_num) // returns a vector in which all the elements are int vectors which represent 
{																	//the adjacent edges atached to each node. 
	edge_data iter_info(0, 0, 0, 0); //calls the default constructor

	int node_ID;
	int edge_ID;
	int edges;
	std::vector<int> adjlist;
	std::vector<std::vector<int>> adjlistholder;

	edges = graph.size();
	node_ID = 0;
	edge_ID = 0;
	adjlistholder.begin();
	adjlist.begin();

	iter_info.update_data(graph[edge_ID].head, graph[edge_ID].tail, graph[edge_ID].cost, graph[edge_ID].cap);

	for (; node_ID <= node_num; node_ID++);
	{
		for (; edge_ID < edges; edge_ID++)
		{
			if (iter_info.head_val() == node_ID)
			{
				adjlist[node_ID] += edge_ID;

			}
			iter_info.update_data(graph[edge_ID].head, graph[edge_ID].tail, graph[edge_ID].cost, graph[edge_ID].cap); // grab new info

			adjlistholder[node_ID] = adjlist; //add list of adjacent edges to adjacent edge superstructure
			adjlist.clear(); // empty adjlist for next iteration
		}
	}

	return(adjlistholder);
}

vector<vector<string> > readIn(vector<vector<string> > vect, string filename)
{

	fstream file;
	file.open(filename, ios::in);
	int i = 0, k = 0;
	string temp;
	string temp1;
	string temp2;
	string temp3;
	int colnum;

	// detect the number of elements in a row
	getline(file, temp);
	colnum = std::count(temp.begin(), temp.end(), ',');
	// return to beginning of file
	file.seekg(0);

	// make sure there are no extra lines at the end of the file nor spaces!!!!!
	while (!file.eof())
	{

		// initial push into one dimention 
		vect.push_back(vector<string>());
		for (k = 0; k < colnum; k++)
		{
			getline(file, temp, ',');
			vect[i].push_back(temp);
		}

		getline(file, temp2);
		vect[i].push_back(temp2);

		//if (temp[0] == " ")
		//{
		//	// swap each character and pop one off the end
		//}



		//stringstream ss;
		//ss << i;
		//temp3 = ss.str();

		//vect[i].push_back(temp3);



		i++;

	}








	return vect;

}

//
//[city name][][][][][][][][][][][][][][][][][][][][][][][][][][][]
//[population][][][][][][][][][][][][][][][][][][][][][][][][][][][]
//[city number][][][][][][][][][][][][][][][][][][][][][][][][][][][]
//[index number][][][][][][][][][][][][][][][][][][][][][][][][][][][]...........n
//
//



int main()

{

	vector<vector<string> > cities, asylumapps, asylumgranted;

	cities = readIn(cities, "CityPop200K.csv");
	asylumapps = readIn(asylumapps, "AsylumAppsFirstTime.csv");
	asylumgranted = readIn(asylumgranted, "AsylumDecisionsFirstTime.csv");

	// displaying the fields
	// note that if you swap the j parameter for 0,1,2,3 
	// you will be able to search on that specific field/ column etc. 
	// please check formatting of file as the delimiting character is a comma!!!!!!!!!!!!

	for (int j = 0; j < cities[1].size(); j++)
	{
		cout << cities[1][j] << ' ';
	}
	cout << "\n";

	for (int j = 0; j < asylumapps[1].size(); j++)
	{
		cout << asylumapps[1][j] << ' ';
	}
	cout << "\n";

	for (int j = 0; j < asylumgranted[1].size(); j++)
	{
		cout << asylumgranted[1][j] << ' ';
	}
	cout << "\n";

	//costs
	double pref[HOSTNUM][SOURCENUM] = {};
	double balance[HOSTNUM] = {};

	// edge capacity
	double capacity[HOSTNUM] = {};

	// edge flow decision variables
	double xflow[HOSTNUM][SOURCENUM] = {};
	double yflow[HOSTNUM] = {};

	double country_pref[28][SOURCENUM] = {};
	double app_array[28][SOURCENUM] = {};
	double grant_array[28][SOURCENUM] = {};
	double app_totals[SOURCENUM] = {};

	for (int j = 0; j < 28; j++)
	{
		for (int i = 0; i < SOURCENUM; i++)
		{
			app_array[j][i] = atof(asylumapps[j + 1][i + 1].c_str());
			grant_array[j][i] = atof(asylumgranted[j + 1][i + 1].c_str());
			app_totals[i] = app_totals[i] + app_array[j][i];
		}
	}

	for (int j = 0; j < 28; j++)
	{
		for (int i = 0; i < SOURCENUM; i++)
		{
			country_pref[j][i] = app_array[j][i] / app_totals[i];
		}
	}

	int skip = 1;
	for (int j = 1; j < cities.size(); j++)
	{
		int tempid = atoi(cities[j][2].c_str());

		if (tempid == 0)
		{
			skip = skip + 1;
			continue;
		}

		for (int i = 0; i < SOURCENUM; i++)
		{
			pref[j - skip][i] = country_pref[tempid - 1][i];
		}

		capacity[j - skip] = atof(cities[j][1].c_str());
	}

	// Initialize graph as vector of edges
	vector<edge> graph;
	int count = 0;
	for (int j = 0; j < HOSTNUM; j++)
	{
		for (int i = 0; i < SOURCENUM; i++)
		{
			graph.push_back(edge());
			graph[count].cost = pref[j][i];
			graph[count].cap = capacity[j];
			graph[count].head = i;
			graph[count].tail = j + SOURCENUM;
			count++;
		}
	}

	for (int j = 0; j < HOSTNUM; j++)
	{
		graph.push_back(edge());
		graph[count].cost = balance[j];
		graph[count].cap = capacity[j];
		graph[count].head = j + SOURCENUM;
		graph[count].tail = HOSTNUM + SOURCENUM;
		count++;
	}

	// Create residual edges at odd indexes
	unsigned int sz = graph.size();
	for (unsigned int i = 0; i < sz; i++)
	{
		graph.insert(graph.begin() + 2 * i + 1, edge());
		graph[2 * i + 1].cost = -graph[2 * i].cost;
		graph[2 * i + 1].cap = 0;
		graph[2 * i + 1].head = graph[2 * i].tail;
		graph[2 * i + 1].tail = graph[2 * i].head;
	}
	// for (int i = 0; i < 10; i++) {
	// cout << graph[i].cost << ' ' << graph[i].cap << ' ' << graph[i].head << ' ' << graph[i].tail << endl;
	// }

	cin.get();
	cin.get();

	vector<vector<int> > populate_adj_lst(graph, HOSTNUM + SOURCENUM);



}

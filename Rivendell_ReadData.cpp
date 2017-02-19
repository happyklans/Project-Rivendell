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


vector<vector<string> > readIn(vector<vector<string> > vect, string filename);


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
			pref[j-skip][i] = country_pref[tempid - 1][i];
		}

		capacity[j-skip] = atof(cities[j][1].c_str());
	}

	// Initialize graph as vector of edges
	vector<edge> graph;
	int count = 0;
	for (int j = 0; j < HOSTNUM; j++) {
		for (int i = 0; i < SOURCENUM; i++) {
			graph.push_back(edge());
			graph[count].cost = pref[j][i];
			graph[count].cap = capacity[j];
			graph[count].head = i;
			graph[count].tail = j + SOURCENUM;
			count++;
		}
	}
	for (int j = 0; j < HOSTNUM; j++) {
		graph.push_back(edge());
		graph[count].cost = balance[j];
		graph[count].cap = capacity[j];
		graph[count].head = j + SOURCENUM;
		graph[count].tail = HOSTNUM + SOURCENUM;
		count++;
	}

	// Create residual edges at odd indexes
	unsigned int sz = graph.size();
	for (unsigned int i = 0; i < sz; i++) {
		graph.insert(graph.begin() + 2*i + 1, edge());
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
	
}

vector<vector<string> > readIn(vector<vector<string> > vect, string filename)
{

	fstream file;
	file.open(filename, ios::in);
	int i = 0, k=0;
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


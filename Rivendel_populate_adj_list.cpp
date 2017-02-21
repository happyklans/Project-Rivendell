#include <vector>
#include <iostream>
#include <string>
using namespace std;

class edge_data //Class for extracting data from the elements of graph
{
public:
	edge_data(int, int, double, double, );
	void update_data(int, int, double, double);
	int head_val{ return(head); }
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

edge_data::update_data(int newi1, int newi2, double newd1, double newd2) // updates info for iteration
{
	cost = newd1;
	cap = newd2;
	head = newi1;
	tail = newi2;
}

vector<vector<int> > pop_adj_lst(vector<edge> graph, int node_num) // returns a vector in which all the elements are int vectors which represent 
{																	//the adjacent edges atached to each node. 
	new edge_data iter_info (0,0,0,0); //calls the default constructor
	int nodes;
	int node_ID;
	int edge_ID;
	int edges;
	struct temp_passer;
	std::vector<int> adjlist;

	edges = graph.size();
	node_ID = 0;
	adjlist.begin();

	

	for (node_ID <= node_num; node_ID++;)
	{
		
		
		for (edge_ID = 0;  edge_ID < edges; edge_ID++;)
		{
			temp_passer = graph[edge_ID]; //holder which allows the extracter class to get data
			iter_info.update_data(temp_passer.head, temp_passer.tail, temp_passer.cost, temp_passer.cap); // updates the data being iterated over

			if (iter_info.head_val == node_ID)
			{
				adjlist.insert(node_ID, edge_ID);
				cout << adjlist[edge_ID] << "\n";
				
			}
			iter_info.update_data()

		}
	}
	
	return(adjlist);
}
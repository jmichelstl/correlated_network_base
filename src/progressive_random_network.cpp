/*
 * Read in a periodic lattice file, then randomly shuffle the edges, and
 * keep a fraction of edges within a specified range, such that, when the
 * number of retained edges is incremented, the edges already present are kept.
 */

#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <tuple>
#include <fstream>
#include <ctime>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <set>
#include <map>
#include <sstream>

using namespace std;

#define SUCCESS 0
#define FAILURE 1

#define vec2d tuple<double, double>
#define edge_data tuple<int, int, short>

//Read in information about the number of points and edges, and the horizontal
//offset for wrapping bonds from the right to the left edge
int read_periodic_lattice_file(ifstream &infile, vector<vec2d> &points, vector<edge_data> &edges, double &offset){

	string nextline;
	int num_read, num_points, num_edges, idx1, idx2;
	double x, y;
	short factor;

	//If a header can be read, obtain the number of points and edges.
	//Otherwise, return and indicate failure.
	if(!infile.eof() && infile.is_open()){
		getline(infile, nextline);
		num_read = sscanf(nextline.c_str(), "%d %d %lf", &num_points, &num_edges, &offset);
		if(num_read == 3){
			if(num_points <= 0 || num_edges <= 0 || offset < 0){
				cerr << "There was a bad header.\n";
				infile.close();
				return FAILURE;
			}
		}
		else return FAILURE;
	}
	else{
		if(infile.is_open()) infile.close();
		return FAILURE;
	}

	//If the file is still open, read points until all expected point have
	//been read, or the end of the file is reached.
	while(! infile.eof() && points.size() < num_points){
		getline(infile, nextline);
		num_read = sscanf(nextline.c_str(), "%lf %lf", &x, &y);
		if(num_read == 2){
			points.push_back(make_tuple(x, y));
		}
	}
	
	//If there are too few points, indicate failure and abort.
	if(points.size() < num_points){
		infile.close();
		return FAILURE;
	}

	//If the file is still open, read edges until all expected edges have
	//been read, or until the end of the file is reached.
	while(!infile.eof() && edges.size() < num_edges){
		getline(infile, nextline);
		num_read = sscanf(nextline.c_str(), "%d %d %hd", &idx1, &idx2, &factor);
		if(num_read == 3){
			edges.push_back(make_tuple(idx1, idx2, factor));
		}
	}
	infile.close();

	//If the number of edges is equal to the expected number, indicate
	//success. Otherwise, indicate failure.
	if(edges.size() == num_edges) return SUCCESS;
	else return FAILURE;
}

//Perform a random shuffling of the bonds
void shuffle_edges(vector<edge_data> &bonds){

	boost::random::mt19937 gen(time(0));
	int pos, swap, length;
	edge_data temp;

	length = bonds.size();
	boost::random::uniform_int_distribution<> dist(0, length - 1);

	//At each bond index, pick a random bond with which to swap that bond.
	for(pos = 0; pos < length; pos ++){
		swap = dist(gen);
		temp = bonds[swap];
		bonds[swap] = bonds[pos];
		bonds[pos] = temp;
	}
}

//When a subset of bonds is retained, produce a condensed list of just those
//points that participate in a retained bond. Reassign point indices within
//retained bonds.
void reassign_indices(vector<vec2d> &points, vector<edge_data> &edges){

	set<int> kept_points;
	map<int, int> assignments;
	int pindex;
	vector<vec2d> new_points;
	vector<edge_data> new_edges;

	//Make a record of the retained points
	for(edge_data next_edge : edges){
		kept_points.insert(get<0>(next_edge));
		kept_points.insert(get<1>(next_edge));
	}

	//Map indices of retained points to the sequence in which they occur
	//in the set of all retained points.
	pindex = 0;
	for(auto iter = kept_points.begin(); iter != kept_points.end(); iter++){
		assignments.insert(make_pair(*iter, pindex));
		pindex ++;
	}

	for(int next_index : kept_points){
		new_points.push_back(points[next_index]);
	}

	//Create a new set of edges with new point indices.
	for(edge_data next_edge : edges){
		new_edges.push_back(make_tuple(assignments[get<0>(next_edge)], assignments[get<1>(next_edge)], get<2>(next_edge)));
	}

	points = new_points;
	edges = new_edges;
}

//Given a set of points and edges, and an offset, keep a subset of edges,
//cull all points not included in a retained edge, and reassign point indices.
//Finally, write all data to a new network file.
void write_plf(vector<vec2d> points, vector<edge_data> edges, double offset, double bond_p, string file_name){

	ofstream output;

	//Produce a subset of all edges in the network containing the specified
	//portion of the original edges.
	vector<edge_data> subset(edges.begin(),edges.begin()+(int)(edges.size()*bond_p));

	//Cull unused points from the description of the network
	reassign_indices(points, subset);
	
	//Write the header, followed by a list of points, and then a list of
	//edges.
	output.open(file_name, ofstream::out);
	output << points.size() << "\t" << subset.size() << "\t" <<offset<<"\n";

	for(vec2d next_point : points){
		output << get<0>(next_point) <<"\t"<< get<1>(next_point) <<"\n";
	}

	for(edge_data next_edge : subset){
		output << get<0>(next_edge) << "\t" << get<1>(next_edge) << "\t" << get<2>(next_edge) << "\n";
	}

	output.close();
}

//Read an argument containing a file name, a minimum and maximum portion of
//bonds to keep, an increment for the bond portion, and a base name. Shuffle the
//edges, and generate periodic lattice files in which each desired portion of 
//edges is retained.
int main(int argc, char **argv){

	vector<vec2d> points;
	vector<edge_data> edges, shuffled_edges;
	double min_p, max_p, inc, curr_p, offset, num_read;
	ifstream input;
	ostringstream name_stream(ostringstream::ate);

	//If there are too few arguments, abort
	if(argc < 6){
		cerr << "To few arguments were provided.\n";
		return FAILURE;
	}
	
	//Validate input
	string file_name(argv[1]);
	
	num_read = sscanf(argv[2], "%lf", &min_p);
	if(num_read < 1 || min_p < 0){
		cerr << "Invalid input.\n";
		return FAILURE;
	}

	num_read = sscanf(argv[3], "%lf", &max_p);
	if(num_read < 1 || max_p < 0 || max_p < min_p || max_p > 1){
		cerr << "Invalid input.\n";
		return FAILURE;
	}

	num_read = sscanf(argv[4], "%lf", &inc);
	if(num_read < 1 || inc <= 0){
		cerr << "Invalid input.\n";
		return FAILURE;
	}

	string base_name(argv[5]);

	//Read in network data, and generate networks with the desired bond
	//portions. Generate names by appending a bond fraction to the base
	//name.
	
	input.open(file_name, ifstream::in);
	if(! input.is_open()){
		cerr << "The file could not be read.\n";
		return FAILURE;
	}

	if(read_periodic_lattice_file(input, points, edges, offset) != SUCCESS){
		cerr << "Not all information could be read.\n";
		return FAILURE;
	}
	
	shuffle_edges(edges);
	
	name_stream << setprecision(3);

	for(curr_p = min_p; curr_p <= max_p; curr_p += inc){
		name_stream.str(base_name);
		name_stream.flush();
		name_stream << "_p" << (int) (curr_p*100) << ".dat";
		name_stream.flush();
		write_plf(points, edges, offset, curr_p, name_stream.str());
	}

	return SUCCESS;
}

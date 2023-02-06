#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <tuple>
#include <string>
#include <stdio.h>
#include <unordered_set>
#include <fftw3.h>
#include "file_lists.hpp"
#include <unistd.h>
#include <stdlib.h>
#include <cstring>
#include <cmath>

using namespace std;

#define TOL 1e-8
#define min(a, b) a < b ? a : b
#define max(a, b) a > b ? a : b

//Knuthian recipe for komputing a kwality integer hash function
#define knuth_hash(a) (2654435769*a) >> 16

//Data type to store two-dimensional vectors
typedef tuple<double, double> vec2d;

//Data type to represent a point
struct Point{

	Point(){
	}

	Point(double myX, double myY) : x(myX), y(myY) {
	}

	bool operator == (const Point& point) const{
		return(abs(point.x - x) < TOL && abs(point.y - y) < TOL);
	}

	double x, y;

	//Define an order relation for sorting
	friend bool operator < (const Point &p1, const Point &p2){
		if(p2.y > p1.y + TOL) return true;
		else if(p1.y > p2.y + TOL) return false;
		else if(p2.x > p1.x + TOL) return true;
		else return false;
    	}
};

//Data type to represent an edge, specified by the integer indices of its
//endpoints in a master list.
struct IPair{

	IPair(){}

	IPair(int arg1, int arg2) : idx1(arg1), idx2(arg2) {}

	int idx1, idx2;

	bool operator == (const IPair &ip) const{
		return ip.idx1 == idx1 && ip.idx2 == idx2;
	}
};

//Hashing function for storing index pairs
namespace std {

    template<> struct hash<IPair>{
        typedef size_t result_type;
        typedef IPair argument_type;

        size_t operator() (const IPair& ip) const;
    };

    size_t hash<IPair>::operator() (const IPair& ip) const {
        return (size_t) knuth_hash(ip.idx1) ^ knuth_hash(ip.idx2);
    }
}

//Split a line into pieces, using whitespace characters as delimiters
vector<string> split_line(string line){
    vector<string> result;
    char *full_line, *token;

    full_line = (char *) malloc((1 + line.size()) * sizeof(char));
    strcpy(full_line, line.c_str());

    token = strtok(full_line, " \t");
    while(token != NULL){
        result.push_back(string(token));
        token = strtok(NULL, " \t");
    }

    free(full_line);

    return result;
}

//Read a header, a point list and an edge list from a network description file.
//Produce a map from points to integer indices, and a list of edges specified
//as tuples of integer indices.
bool read_periodic_network_file(string file_name, map<Point, int> &pmap, unordered_set<IPair> &edges, double &period){

	ifstream input;
	int num_read, npoints, nedges, pindex, idx1, idx2;
	string nextline;
	map<int, int> translation;
	double x, y;
	Point p;

	pindex = 0;

	//Make sure the file can be read
	input.open(file_name, ifstream::in);
	if(! input.is_open()){
		cerr << "The file " << file_name << " could not be read.\n";
		return false;
	}

	else if(input.eof()){
		cerr << "The end of the file was reached prematurely.\n";
		input.close();
		return false;
	}

	//Empty all data structures
	pmap.clear();
	edges.clear();

	//Read the header
	getline(input, nextline); 
	num_read = sscanf(nextline.c_str(), "%d %d %lf", &npoints, &nedges, &period);
	if(num_read < 3 || npoints <= 0 || nedges <= 0 || period < 0){
		cerr << "The header was invalid.\n";
		input.close();
		return false;
	}

	//Obtain points
	while(! input.eof() && pmap.size() < npoints){
		getline(input, nextline);
		num_read = sscanf(nextline.c_str(), "%lf %lf", &x, &y);
		if(num_read == 2){
			p.x = x;
			p.y = y;

			pmap.insert(make_pair(p, pindex++));
		}
	}

	if(pindex < npoints){
		cerr << "Too few valid points were read.\n";
		input.close();
		return false;
	}

	//Obtain edges
	while(! input.eof() && edges.size() < nedges){
		getline(input, nextline);
		num_read = sscanf(nextline.c_str(), "%d %d", &idx1, &idx2);
		if(num_read == 2){
			if(idx1 >= 0 && idx2 >= 0 && idx1 < npoints && idx2 < npoints){
				edges.emplace(IPair(min(idx1,idx2),max(idx1,idx2)));
			}
			else{
				cerr << "Out of bounds indices were found.\n";
				input.close();
				return false;
			}
		}
	}

	input.close();

	if(edges.size() < nedges){
		cerr << "Too few edges were read.\n";
		return false;
	}


	return true;
}

//Given a starting point, primitive lattice vectors, a number of
//steps along each lattice vector direction, and a list of Wigner-Seitz cell
//edges, find the edges in the Wigner-Seitz cell for each Bravais lattice
//point.
void tally_edges(map<Point, int> pmap, unordered_set<IPair> edge_set, vector<vector<vec2d>> ws_edges, vector<double> &tallies, double maxx, double period, Point begin, vec2d v1, vec2d v2, int n1, int n2){

	int idx1, idx2, eIdx1, eIdx2, tally;
	Point neighbor;

	//Coordinates of a Bravais lattice site and a neighbor site
	double bx, by, nx, ny, mean = 0;

	//Iterate over each Bravais site under consideration, and
	//counter the number of bonds in its Wigner-Seitz cell
	for(idx1 = 0; idx1 < 2*n1; idx1 ++){
		for(idx2 = 0; idx2 < 2*n2; idx2 ++){
			tally = 0;
			bx = begin.x + get<0>(v1)*idx1 + get<0>(v2)*idx2;
			by = begin.y + get<1>(v1)*idx1 + get<1>(v2)*idx2;

			for(vector<vec2d> next_pair : ws_edges){
				neighbor.x = bx + get<0>(next_pair[0]);
				if(neighbor.x > maxx){
					neighbor.x -= floor((neighbor.x-maxx)/period)*period;
				}
				neighbor.y = by + get<1>(next_pair[0]);

				if(pmap.find(neighbor) != pmap.end()){
					eIdx1 = pmap[neighbor];
				}
				else continue;

				neighbor.x = bx + get<0>(next_pair[1]);
				if(neighbor.x > maxx){
					neighbor.x -= floor((neighbor.x-maxx)/period)*period;
				}
				neighbor.y = by + get<1>(next_pair[1]);

				if(pmap.find(neighbor) != pmap.end()){
					eIdx2 = pmap[neighbor];
				}
				else continue;

				if(edge_set.find(IPair(min(eIdx1,eIdx2),max(eIdx1,eIdx2))) != edge_set.end()){
					tally ++;
				}
			}
			tallies[idx1*n2 + idx2] = tally;
			mean += tally;
		}
	}

	mean /= tallies.size();
	for(idx1 = 0; idx1 < tallies.size(); idx1++){
		tallies[idx1] -= mean;
	}
}

//Given a list of numbers of edges in a set of Wigner-Seitz cells, find the
//two-point correlation function for edge count as a function of the separation
//of two Bravais sites, expressed in terms of the coefficients by which
//the primtive lattice vectors must be multiplied to yield the displacement
//vector between the two sites.
void get_corr_funcs(vector<double> counts, double *corr_table, int n1, int n2, bool normalize){

	int idx1, idx2, offset1, offset2, currIdx, neighIdx, corrIdx, count = 0;
	double normFact;

	//Initialize all correlations to 0
	for(idx1 = 0; idx1 < n1; idx1 ++){
		for(idx2 = 0; idx2 < n2; idx2++){
			corr_table[idx1*n2 + idx2] = 0;
		}

	}

	//Normalize totals by dividing by the number of Bravais sites
	//considered, and then by dividing each correlation function by
	//the value at zero displacement.
	for(idx1 = 0; idx1 < n1; idx1 ++){
		for(idx2 = 0; idx2 < n2; idx2++){
			currIdx = idx1*n2 + idx2;

			for(offset1 = 0; offset1 < n1; offset1++){
				for(offset2 = 0; offset2 < n2; offset2++){
					count ++;
					corrIdx = offset1*n2 + offset2;
					neighIdx = (offset1+idx1)*n2 + offset2 + idx2;

					corr_table[corrIdx] += counts[currIdx]*counts[neighIdx];
				}
			}
		}
	}

	//Normalize totals by dividing by the number of Bravais sites
	//considered, and then by dividing each correlation function by
	//the value at zero displacement.
	normFact = normalize ? 1 / corr_table[0] : 1. / (n1*n2);
	for(idx1 = 0; idx1 < n1; idx1++){
		for(idx2 = 0; idx2 < n2; idx2++){
			corr_table[idx1*n2 + idx2] *= normFact; 
		}
	}
}

//Attempt to read a line from a file, and print an error message if the
//attempt is unsuccessful.
//
//Return value
//	True, if the attempt succeeded
//	False, otherwise
bool try_to_read(ifstream &input, string &line){

	while(! input.eof()){
		getline(input, line);
		//Ignore blank lines, and lines beginning with a # sign
		if(! line.compare("") == 0 && line[0] != '#') return true;
	}
	return false;
}

//Parse an input file to harvest information about the network description
//file, the primitive lattice vectors used to generate the network, the
//total number of steps to take along each primitive vector's direction,
//the coordinates of the first Bravais lattice site to consider, and the
//bonds within a Wigner-Seitz cell.
bool parse_job_file(string name, string &file_desc, string &tag, vec2d &v1, vec2d &v2, int &n1, int &n2, Point &offset, vector<vector<vec2d>> &ws_edges){

	ifstream input;
	string nextline;
	double dx1, dy1, dx2, dy2;
	int num_read, idx;
	vector<string> tokens;

	input.open(name, ifstream::in);
	if(! input.is_open()){
		cerr << "The input file could not be read.\n";
		return false;
	}

	//Find the network description file naming pattern
	if(try_to_read(input, nextline)){

	        tokens = split_line(nextline);

        	//If at least one string of non-whitespace characters is found,
	        //initialize the base name, and, optionally, a tag to append to
		//the base name.
		if(tokens.size() >= 1){
			file_desc = tokens[0];
			if(tokens.size() > 1) tag = tokens[1];
        	}
		else{
			cerr << "No valid file name pattern was read.\n";
			input.close();
			return false;
		}
	}
	else return false;

	//Read the primitive lattice vectors
	if(try_to_read(input, nextline)){
		num_read = sscanf(nextline.c_str(), "%lf %lf", &dx1, &dy1);
		if(num_read == 2){
			v1 = make_tuple(dx1, dy1);
		}
		else{
			cerr << "No valid lattice vector was read.\n";
			input.close();
			return false;
		}
	}
	else return false;
	
	if(try_to_read(input, nextline)){
		num_read = sscanf(nextline.c_str(), "%lf %lf", &dx1, &dy1);
		if(num_read == 2){
			v2 = make_tuple(dx1, dy1);
		}
		else{
			cerr << "No valid lattice vector was read.\n";
			input.close();
			return false;
		}
	}
	else return false;
	
	//Read the number of steps along each primitive lattice vector
	if(try_to_read(input, nextline)){
		num_read = sscanf(nextline.c_str(), "%d %d", &n1, &n2);
		if(num_read < 2 || n1 <= 0 < n2 <= 0){
			cerr << "Valid step counts could not be read.\n";
			input.close();
			return false;
		}
	}
	else return false;
	
	//Read the offset from the origin at which the first Bravais
	//lattice site should be found
	if(try_to_read(input, nextline)){
		num_read = sscanf(nextline.c_str(), "%lf %lf", &dx1, &dy1);
		if(num_read == 2){
			offset.x = dx1;
			offset.y = dy1;
		}
		else{
			cerr << "No valid vector was read.\n";
			input.close();
			return false;
		}
	}
	else return false;
	
	//Read lists of edges in the W-S cell, specified as displacements
	//of end points from a Bravais lattice site
	while(! input.eof()){
		if(! try_to_read(input, nextline)){
			break;
		}

		num_read = sscanf(nextline.c_str(), "%lf %lf %lf %lf", &dx1, &dy1, &dx2, &dy2);
		if(num_read == 4){
			ws_edges.push_back(vector<vec2d>());
			idx = ws_edges.size() - 1;
			ws_edges[idx].push_back(make_tuple(dx1, dy1));
			ws_edges[idx].push_back(make_tuple(dx2, dy2));
		}
		else cerr << "Invalid input found where edge end points were expected.\n";
	}

	input.close();

	if(ws_edges.size() > 0) return true;

	else{
		cerr << "No edges within a Wigner-Seitz cell were specified.\n";
		return false;
	}
}

//Read an input file describing necessary parameters. Read each network
//description file satisfying a given pattern and find the two-point
//correlation function for the bond density.
int main(int argc, char **argv){

	vector<string> in_files, descs;
	string file_desc, tag;
       	vec2d v1, v2, b1, b2;
	Point offset;
	int n1, n2, iter, idx1, idx2, tstride;
	vector<vector<vec2d>> ws_edges;
	map<Point, int> pmap;
	unordered_set<IPair> edges;
	double period, maxx, norm_fact, kx, ky, pspec, real, cplx;
	char c;
	char *batch_file = NULL, *ext = NULL;
	vector<double> edge_counts;
	fftw_plan plan;
	vector<vector<string>> out_names;
	bool normalize = true;
	FILE *report;

	//Arrays to hold the real-space representation of the two-point
	//bond density correlation function, and its discrete power spectrum
	double *corr_table;
	fftw_complex *transform;

	while((c = getopt(argc, argv, "b:e:m")) != -1){
		switch(c) {
			case 'b':
				batch_file = (char *) malloc((1 + strlen(optarg))*sizeof(char));
				strcpy(batch_file, optarg);
				break;
			case 'e':
				ext = (char *) malloc((1 + strlen(optarg))*sizeof(char));
				strcpy(ext, optarg);
				break;
			case 'm':
				//Find means of correlation products, rather
				//than normalizing the maximum value to 1
				normalize = false;
				break;
			case '?' :
				if(optopt == 'b'){
					cerr << "Option \"b\" requires a file name.\n";
					return -1;
				}
				if(optopt == 'e'){
					cerr << "Option \"e\" requires an extension name.\n";
					return -1;
				}
				else if(isprint(optopt)){
					fprintf(stderr, "Unrecognized option: %c\n", optopt);
				}
				else{
					cerr << stderr, "Unknown option character.\n";
                		}
				break;
			default:
				break;
		}
	}


	//Ensure a job file was provided, and ensure all information can
	//be read from it.	
	if(! batch_file){
		cerr << "No job file was specified.\n";
		return -1;
	}

	if(! parse_job_file(string(batch_file), file_desc, tag, v1, v2, n1, n2, offset, ws_edges)){
		cerr << "Parsing of the job script failed.\n";
		free(batch_file);
		return -2;
	}
	else free(batch_file);

	//Enumerate all network description files matching the pattern given
	//in the job file, and generate accompanying names for output files
	//reporting correlation functions and their power spectra.
	descs.push_back("gor");
	descs.push_back("pwr");
	if(ext){
		get_file_lists(file_desc, tag, string(ext), in_files, descs, out_names);
		free(ext);
	}
	else{
		get_file_lists(file_desc, tag, "dat", in_files, descs, out_names);
	}

	if(in_files.size() == 0){
		cerr << "No files matching the given pattern were found.\n";
		exit(-3);
	}

	//Initialize a list to hold the edge count for each Wigner-Seitz cell
	edge_counts = vector<double>(4*n1*n2, 0.);

	//Find the reciprocal lattice vectors associated with the direct-space
	//primtive lattice vectors
	norm_fact = 2*M_PI/(get<1>(v1)*get<0>(v2) - get<0>(v1)*get<1>(v2));
	b1 = make_tuple(norm_fact/n1*(-get<1>(v2)), norm_fact/n1*get<0>(v2));

	norm_fact = 2*M_PI/(get<1>(v2)*get<0>(v1) - get<0>(v2)*get<1>(v1));
	b2 = make_tuple(norm_fact/n2*(-get<1>(v1)), norm_fact/n2*get<0>(v1));

	//Set up the workspace for computing the density correlation function
	//and its power spectrum
	tstride = n2/2 + 1;
	corr_table = (double *) malloc(sizeof(double) * n1*n2);
	transform =(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*n1*tstride);
	plan = fftw_plan_dft_r2c_2d(n1, n2, corr_table, transform, FFTW_ESTIMATE);

	//Iterate over each network description file. Find the WS cell edges
	//for each point to be considered, compute an averaged set of
	//bond density correlation functions, and write a report of findings
	//to a file.
	for(iter = 0; iter < in_files.size(); iter ++){

		if(! read_periodic_network_file(in_files[iter],pmap, edges, period)){
			cerr << "The network file " << in_files[iter] << " was invalid.\n";
			continue;
		}

		//Find edge counts, and compute the two-point edge count
		//correlator and its discrete Fourier transform
		tally_edges(pmap, edges, ws_edges, edge_counts, maxx, period, offset, v1, v2, n1, n2);
		get_corr_funcs(edge_counts, corr_table, n1, n2, normalize);
		//We've planned the work; now, let's work the plan.
		fftw_execute(plan);

		//Report the two-point correlator as a function of displacement,
		//and the power spectrum as a function of wave vector.
		report = fopen(out_names[0][iter].c_str(), "w");
		for(idx1 = 0; idx1 < n1; idx1 ++){
			for(idx2 = 0; idx2 < n2; idx2 ++){
				fprintf(report, "%d\t%d\t%1.12lf\n",idx1,idx2,corr_table[idx1*n2 + idx2]);
			}
		}
		fclose(report);

		report = fopen(out_names[1][iter].c_str(), "w");
		for(idx1 = 0; idx1 < n1/2 + 1; idx1 ++){
			for(idx2 = 0; idx2 < n2/2 + 1; idx2 ++){
				kx = idx1*get<0>(b1) + idx2*get<0>(b2);
				ky = idx1*get<1>(b1) + idx2*get<1>(b2);
				real = transform[idx1*tstride + idx2][0];
				cplx = transform[idx1*tstride + idx2][1];
				pspec = real*real + cplx*cplx;
				fprintf(report, "%1.8lf\t%1.8lf\t%1.12lf\n", kx, ky, pspec);
			}
		}
		fclose(report);
	}

	fftw_destroy_plan(plan);
	fftw_free(transform);
	free(corr_table);

	return 0;
}

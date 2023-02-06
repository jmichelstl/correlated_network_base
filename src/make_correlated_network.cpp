/*
Author: Jonathan Michel

This program first forms a lattice, consisting of ordered pairs of connected
Cartesian points in two dimensions. The program then optionally removes bonds
such that a certain fraction is retained. Bonds are randomly selected for
retention, and candidates are kept with a probability that increases with
the number of adjacent bonds that have already been retained. Each row is
specified as a repeating pattern of points, and the overall lattice is specified
as a repeating set of rows. The first piece of information in a row 
specification is a horizontal offset from the left edge of the space to
be filled with a lattice to the left-most point in the row. The next piece of
information is the vertical spacing from the current row to the next row in
the pattern. Finally, each point in the pattern is then indicated by a number 
specifying the horizontal spacing from this point to the next point in the 
pattern.
*/

#include "network_utils.h"
#include <iostream>
#include <vector>
#include <map>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <cmath>
#include <fstream>
#include <float.h>
#include <set>
#include <tuple>
#include <unistd.h>
#include "make_network_mesh.hpp"

using namespace std;

typedef tuple<int, int, char> edge_datum;
typedef tuple<int, int, int> itriple;

//The default and maximum values for the minimum angle constraint which may
//be optionally specified when the user chooses to create a network-FEM hybrid
#define DEF_MIN_ANG 30
#define MAX_MIN_ANG 60

//Dilute a network in a correlated manner, with a correlation coefficient c,
//such that a bond selected as a candidate for retention is kept with 
//probability p = (1 - c)^(zmax - nn), where zmax is the maximum number of
//bonds to which a given bond can be adjacent, and nn is the current number
//neighbor bonds.
void correlated_dilution(vector<edge_datum> &edges, vector<Point> &points, double portion, double correlation, int maxz, bool mesh){

     //Edges to retain upon dilution
     vector<edge_datum> keepers;
     vector<Point> new_points;
     edge_datum edat;
     set<int> kept_points;
     int num_points = points.size();
     map<int, int> pmap;

     //Structures to track which edges have been retained, and how many
     //retained edges contain each point
     vector<int> candidates, valence(num_points, 0); 

     //Quantities to determine whether to keep a candidate bond
     double probability, randval;
     int neighbors, rand_pos, edge_index, to_keep, iter;

     //Structures to generate pseudorandom numbers
     const gsl_rng_type *T = gsl_rng_mt19937;
     gsl_rng *r = gsl_rng_alloc(T);

     to_keep = (int) (portion * edges.size());

     //Populate the list of candidates
     for(iter = 0; iter < edges.size(); iter ++) candidates.push_back(iter);

     while(keepers.size() < to_keep){

         //Select a point at random, then find the number of retained neighbor 
         //bonds so far, compute the probability of acceptance, and compare
         //this probability to a random number
         rand_pos = gsl_rng_uniform_int(r, candidates.size());

         edge_index = candidates[rand_pos];
	 edat = edges[edge_index];
         neighbors = valence[get<0>(edat)] + valence[get<1>(edat)];
         probability = pow(1 - correlation, maxz - neighbors);

         //If a bond is selected for retention, take it out of the pool of
         //possibilities, and indicate that it has been kept.
         if(probability > gsl_rng_uniform(r)){
             keepers.push_back(edges[edge_index]);
             candidates.erase(candidates.begin() + rand_pos);
             valence[get<0>(edat)] += 1;
             valence[get<1>(edat)] += 1;
	     if(! mesh){
                 kept_points.insert(get<0>(edat));
                 kept_points.insert(get<1>(edat));
             }
         }
     }

     free(r);

    if(! mesh){
        iter = 0;
        for(int next_index : kept_points){
            pmap.insert(make_pair(next_index, iter));
	    new_points.push_back(points[next_index]);
	    iter ++;
        }

        points = new_points;

        edges.clear();
        for(edge_datum next_edge : keepers){
            edges.push_back(make_tuple(pmap[get<0>(next_edge)], pmap[get<1>(next_edge)], get<2>(next_edge)));
        }
    }

    else edges = keepers;

}

//Multiply each component of an ordered list of doubles by a constant factor
void scale_vector(vector<double>& in, double scale){
    int index;
    for(index = 0; index < in.size(); index++){
        in.at(index) = in.at(index)*scale;
    }
}

/*
Import the rules for creating a lattice from a file. The file should specify
a 2D lattice a set of lines of points. Each line of points should specify
an inital x offset from the left-most x coordinate, a y offset for the
subsequent row of points, and then one x offset for each point in a repeating
pattern defining the row. Each line of points should be defined on a single line
of text, with a new line signifying a new line is now being defined. A blank
line signifies the end of line definitions. After all lines of  points are 
defined, a series of displacements should be defined from points to neighbors 
with which those points share a bond. Displacements for each point with in
each line should be specified as pairs of numbers on the same line, separated
by a space. A blank line should be used to signify that all neighbors for a
point have been specified.

Arguments:
	-rules: A vector of vectors of double-precision floating point numbers
	specifying the means for laying out lattice sites
	-nns: A map from an integer index to a vector of nearest neighbor
	displacements
	-scale: A multiplier for uniformly contracting or dilating the
	lattice, to shrink or expand the separation between lattice sites
	relative to what is specified in the lattice description file

Return value:
	A boolean indicating whether a valid set of lattice rules has been
	read
*/
bool import_lattice(vector<vector<double>>& rules, map<int,vector<vector<double>>>& nns, double scale){
    ifstream latfile;
    string name, nextline;
    bool again, blank, fileopen = false;
    vector<double> rule, nn;
    int lcount = 0, pointiter, nncount = 0, ruleiter = 0;

    //Prompt for file name
    do{
        printf("Enter the lattice file name: ");
        getline(cin, nextline);
        name = split(nextline, ' ')[0];

        if(!name.empty()) latfile.open(name);

        if(!latfile.is_open()){
            again = yesno("No file was read. Try again? ");
            if(!again) return false;
        }
    }while(!latfile.is_open());

    //Process rules until a blank line is reached
    blank = false;
    while(!latfile.eof()){
        lcount ++;
        getline(latfile, nextline);
        if(nextline.empty()) break;

        rule = parse_doubles(split(nextline, ' '));
        if(rule.size() < 3){
            fprintf(stderr, "Too few numbers in rule on line line %d\n",lcount);
            latfile.close();
            return false;
        }
        scale_vector(rule, scale);
        rules.push_back(rule);
    }
    
    //Read nearest neighbor rules
    for(vector<double> nextrule : rules){
       ruleiter ++;
       for(pointiter = 1; pointiter <= nextrule.size() - 2; pointiter ++){
           vector<vector<double>> nextset;
           while(!latfile.eof()){
               lcount ++;
               getline(latfile, nextline);
               if(nextline.empty()) break;
               nn = parse_doubles(split(nextline, ' '));
               if(nn.size() != 2) fprintf(stderr, "Insufficient information for nearest neighbor rule on line %d.\n", lcount);
               else{
                   scale_vector(nn, scale);
                   nextset.push_back(nn);
               }
           }
           nns.insert(pair<int, vector<vector<double>>>(nncount++, nextset));
       }
    }


    latfile.close();
    return true;
}

//Add a single row to a network
void add_row(vector<double> rule, map<int, vector<vector<double>>> nns, int nx, int ny, double cycle_len, double max_x, vector<Point> &points, vector<edge_datum> &edges, double y, int base, set<Point> &pset, double ymax){

    double x, x2, y2;
    char x_offset;
    int row_iter, rule_iter;

    Point p1, p2;
    pair<set<Point>::iterator, bool> ret;


    x = rule[0];
    //Iterate over the repeating unit for each row nx times to create
    //the number of desired repeating units in each row
    for(row_iter = 0; row_iter < nx; row_iter++){
        for(rule_iter = 0; rule_iter < rule.size() - 2; rule_iter++){

            //Given the current x and y positions, create a point,
	    //and if it is not already in the list of lattice points,
	    //update the list of points
	    p1 = Point(x, y, points.size());
	    ret = pset.insert(p1);
	    p1 = *(ret.first);
	    if(ret.second) points.push_back(p1);

	    //For each site in the repeated pattern for a given row,
	    //iterate over the set of nearest neighbors
	    for(vector<double> offset : nns[base + rule_iter]){
                x2 = x + offset[0];
        	y2 = y + offset[1];

                if(y2 > ymax + FLOAT_TOL) continue;

		//Add an offset as needed in the x direction
		//if the neighbor point is out of bounds, such that
		//the lattice creation process wraps around from the
		//right to the left edge.
		if(x2 > max_x + FLOAT_TOL){
		    x2 -= cycle_len;
		    x_offset = 1;
		}
		else if(x2 < -FLOAT_TOL){
		    x2 += cycle_len;
		    x_offset = -1;
		}
		else x_offset = 0;

		p2 = Point(x2, y2, points.size());
		ret = pset.insert(p2);
		p2 = *(ret.first);
		if(ret.second) points.push_back(p2);

	        //Update the list of edges with an edge connecting
	        //p1 to p2
	        edges.push_back(make_tuple(p1.index, p2.index, x_offset));
	    }
	    x += rule[2 + rule_iter];
		
	}
    }
}

//Given a set of rules for lattice sites and neighbors, produce a set of
//edges, consisting of pairs of Cartesian points
void make_edges(vector<vector<double>> rules, map<int,vector<vector<double>>> nns, int nx, int ny, double cycle_len, double max_x, vector<Point> &points, vector<edge_datum> &edges){

    double y = 0, ymax = 0;
    int col_iter, num_rules, base, nn_count;

    set<Point> pset;

    base = 0;
    nn_count = nns.size();

    for(vector<double> next_rule : rules){
        ymax += next_rule[1];
    }
    ymax *= ny;

    //Iterate over the set of row rules ny times to create the desired number
    //of repeating units in the vertical direction
    for(col_iter = 0; col_iter < ny; col_iter++){

	for(vector<double> rule : rules){

            add_row(rule,nns,nx,ny,cycle_len,max_x,points,edges,y,base,pset,ymax);

            y += rule[1];
	    base = (base + rule.size() - 2) % nn_count;
	}
    }

    //Cap off the network with one more row following the first row pattern
    add_row(rules[0],nns,nx,ny,cycle_len,max_x,points,edges,y,base, pset,y);
}


//Given a list of integers, generate all shufflings such that, when partitioned
//into consecutive pairs, the first member of each pair occurs before the
//second element of the pair in the original sequence.
vector<vector<int>> get_shufflings(vector<int> list){

        vector<vector<int>> shufflings;
        int iter, inner;

        if(list.size() <= 2){
                shufflings.push_back(list);
        }

        else{
                for(iter = 1; iter < list.size(); iter++){

                        vector<int> rest;
                        for(inner = 1; inner < list.size(); inner ++){
                                if(inner != iter){
                                        rest.push_back(list[inner]);
                                }
                        }

                        for(vector<int> next_list : get_shufflings(rest)){
                                next_list.insert(next_list.begin(), list[0]);
                                next_list.insert(next_list.begin()+1, list[iter]);
                                shufflings.push_back(next_list);
                        }
                }
        }

        return shufflings;
}

//Find the geometric mean of the angles between pairs of bonds implied by a
//set of bending triples.
double ang_geom_mean(vector<Point> pts, int center, vector<int> neigh_list){

        int iter;
        Point cpt, p1, p2;
        double dx1, dx2, dy1, dy2, l1, l2, product;

        product = 1;
        cpt = pts[center];
        for(iter = 0; iter < neigh_list.size() - 1; iter += 2){
                p1 = pts[neigh_list[iter]];
                p2 = pts[neigh_list[iter + 1]];

                dx1 = p1.x - cpt.x;
                dy1 = p1.y - cpt.y;
                dx2 = p2.x - cpt.x;
                dy2 = p2.y - cpt.y;

                product *= acos((dx1*dx2+dy1*dy2)/sqrt(dx1*dx1+dy1*dy1)/sqrt(dx2*dx2+dy2*dy2));
        }

        return pow(product, 1 / floor(neigh_list.size() / 2));
}

vector<itriple> generate_triples(vector<Point> pts, vector<edge_datum> edges){

        map<int, vector<int>> neigh_map;
        int idx1, idx2, niter;
        double gmean, highmean;
        vector<int> optimum;
        vector<itriple> btriples;

        for(idx1 = 0; idx1 < pts.size(); idx1++){
                neigh_map.insert(make_pair(idx1, vector<int>()));
        }

        //Make a mapping from each vertex index to the a list of the indices of
        //its neighbors.
        for(edge_datum dat : edges){
                idx1 = get<0>(dat);
                idx2 = get<1>(dat);
                neigh_map[idx1].push_back(idx2);
                neigh_map[idx2].push_back(idx1);
        }

        for(auto iter = neigh_map.begin(); iter != neigh_map.end(); iter ++){

                //If the neighbor list has just one element, there is nothing
                //to do.
                if(iter->second.size() < 2) continue;

                //Generate a shuffling of the current neighor list for each
                //unique way to partition the neighbors into pairs.
                vector<vector<int>> shufflings = get_shufflings(iter->second);
                highmean = FLT_MIN;

                for(vector<int> shuffling : shufflings){
                        gmean = ang_geom_mean(pts, iter->first, shuffling);
                        if(gmean > highmean){
                                highmean = gmean;
                                optimum = shuffling;
                        }
                }

                for(niter = 0; niter < optimum.size() - 1; niter += 2){
                        btriples.push_back(make_tuple(iter->first,optimum[niter], optimum[niter+1]));
                }
        }

        return btriples;
}

/*
 *Given a set of points and a set of edges, split each edge into two
 *edges. Augment the list of points to include the midpoints of the previous
 *set of edges.
 */
void add_hinges(vector<Point> &points, vector<edge_datum> &edges, double width){

	vector<edge_datum> new_edges;
	int idx1, idx2, idx3;
	double x1, y1, x2, y2, x3, y3;
	char offset;

	for(edge_datum edat : edges){
		idx1 = get<0>(edat);
		idx2 = get<1>(edat);
		offset = get<2>(edat);

		x1 = points[idx1].x;
		y1 = points[idx1].y;
		x2 = points[idx2].x;
		y2 = points[idx2].y;

		x3 = (x1 + x2 + offset*width) / 2;
		y3 = (y1 + y2) / 2;
		idx3 = points.size();

		points.push_back(Point(x3, y3));
		new_edges.push_back(make_tuple(idx1, idx3, 0));
		new_edges.push_back(make_tuple(idx3, idx2, offset));
	}

	edges = new_edges;
}

//Obtain a Delaunay finite element mesh that includes the points in the network.
//Reassign the point list and the indices defining edges, as well as creating
//triangular facets and computing triangle areas.
void network_to_mesh(vector<Point> &in_pts, double offset, double min_ang, double max_area, vector<itriple> &facets, vector<double> &areas, vector<int> &lookup){

    vector<Point> out_pts;
    int idx, v1, v2, v3;
    map<Point, int> pmap;
    double dx1, dy1, dx2, dy2;

    //Create the new network-FEM hybrid
    make_network_mesh(in_pts, out_pts, facets, offset, min_ang, max_area);

    //Make a mapping from old to new point indices
    idx = 0;
    for(Point p : out_pts){
        pmap.insert(make_pair(p, idx++));
    }

    for(idx = 0; idx < in_pts.size(); idx++){
        lookup.push_back(pmap[in_pts[idx]]);
    }

    //Replace the point set
    in_pts.assign(out_pts.begin(), out_pts.end());

    //Compute triangle areas
    for(itriple next_facet : facets){
        v1 = get<0>(next_facet);
	v2 = get<1>(next_facet);
	v3 = get<2>(next_facet);

	dx1 = in_pts[v2].x - in_pts[v1].x;
	dy1 = in_pts[v2].y - in_pts[v1].y;
	dx2 = in_pts[v3].x - in_pts[v1].x;
	dy2 = in_pts[v3].y - in_pts[v1].y;

	areas.push_back(.5 * abs(dx1*dy2 - dy1*dx2));
    }
}

//If a network-mesh hybrid has been created, point indices in edge 
//specifications should be updated.
void replace_edges(vector<edge_datum> &in_edges, vector<int> lookup){

    vector<edge_datum> new_edges;

    for(edge_datum edat : in_edges){
        new_edges.push_back(make_tuple(lookup[get<0>(edat)], lookup[get<1>(edat)], get<2>(edat)));
    }

    in_edges.assign(new_edges.begin(), new_edges.end());
}

//Report a periodic network to a file
void report_plf(FILE *report_file, vector<Point> points, vector<edge_datum> edges, double offset, vector<itriple> triples){

    //First report a header containing the number of points, the number of
    //edges, and the dimensions in the x and y directions
    fprintf(report_file, "%ld\t%ld\t%2.12lf", points.size(), edges.size(), offset);
    
    if(triples.size() > 0){
	    fprintf(report_file, "\t%ld", triples.size());
    }
    fprintf(report_file, "\n");

    //If the portion of bonds retained is less than 1, only report those points
    //that are still included in the network. Otherwise, report all points
    for(Point p : points){
        fprintf(report_file, "%2.12lf\t%2.12lf\n", p.x, p.y);
    }

    //Finally report each edge as a pair of points, followed by a pair of
    //offsets to apply to the second point, in units of the x and y dimensions
    //of the space in which the lattice is created.

    for(edge_datum next_edge : edges){
        fprintf(report_file, "%d\t%d\t%hhd\n", get<0>(next_edge), get<1>(next_edge), get<2>(next_edge));
    }

    if(triples.size() > 0){
	    for(itriple t : triples){
		    fprintf(report_file, "%d\t%d\t%d\n", get<0>(t),get<1>(t),get<2>(t));
	    }
    }
}

//Report a network-mesh hybrid to a file
void report_mesh(FILE *fh, vector<Point> points, vector<edge_datum> edges, vector<itriple> triples, vector<itriple> triangles, vector<double> areas, double offset){

    //Print header explaining the number of each attribute present
    fprintf(fh, "%ld\t%ld\t%ld", points.size(), edges.size(), triangles.size());

    if(triples.size() > 0){
        fprintf(fh, "\t%ld", triples.size());
    }

    fprintf(fh, "\t%12.8lf\n", offset);

    //Report vertices, edges, mesh cells, and mesh cell areas
    for(Point p : points){
        fprintf(fh, "%12.8lf %12.8lf\n", p.x, p.y);
    }

    for(edge_datum next_edge : edges){
        fprintf(fh, "%d %d %d\n", get<0>(next_edge), get<1>(next_edge), get<2>(next_edge));
    }

    for(itriple next_cell : triangles){
        fprintf(fh, "%d %d %d\n", get<0>(next_cell), get<1>(next_cell), get<2>(next_cell));
    }

    for(double next_area : areas){
        fprintf(fh, "%12.8lf\n", next_area);
    }
}

int main(int argc, char **argv){

    int nx, ny, num_read, maxz, rule_index;
    double scale, portion, correlation, cycle_len, far_right, max_x;
    size_t size;
    bool make_a_mesh = false;
    double min_ang = DEF_MIN_ANG, max_area = 0;

    //Data structures to hold rules for lattice generation
    vector<vector<double>> rules;
    map<int, vector<vector<double>>> nns;

    //Lists of edges and points
    vector<edge_datum> edges;
    vector<Point> points;
    vector<itriple> triples;

    //Data structures for optional FEM mesh
    vector<int> lookup;
    vector<double> areas;
    vector<itriple> facets;

    //Data values associated with creating a seed for pseudorandom number
    //generation
    FILE *ranfile;
    unsigned seed;
    char num_string[11];

    //Map from original point indices to the indices of points after dilution,
    //which may exclude some points from the network
    //map<int, int> pmap;

    string response;

    FILE *report_file;

    char c;
    bool get_triples = false, hinge = false;

    //Check to see if a triples should be determined
    while((c = getopt(argc, argv, "a:hmq:t")) != -1){
        switch(c) {
            case 'a':
                if(sscanf(optarg, "%lf", &max_area) == 1){
		    if(max_area <= 0){
		        cerr << "The maximum area must be positive.\n";
			max_area = 0;
		    }
		}
		else cerr << "The option \"-a\" requires a positive number.\n";
                break;
            case 'h':
                hinge = true;
                break;
            case 'm':
		make_a_mesh = true;
		break;
            case 'q':
                if(sscanf(optarg, "%lf", &min_ang) == 1){
		    if(min_ang <= 0 || min_ang > MAX_MIN_ANG){
		        cerr << "The minimum angle must be positive and less then or equal to " << MAX_MIN_ANG << ".\n";
			min_ang = DEF_MIN_ANG;
		    }
		}
		else cerr << "The option \"-q\" requires a number between 0 and " << MAX_MIN_ANG << ".\n";
		break;
            case 't':
                get_triples = true;
                break;
            case '?' :
		if(optopt == 'a'){
		    cerr << "The option \"-a\" requires a positive number.\n";
		}
		else if(optopt == 'q'){
		    cerr << "The option \"-q\" requires a number between 0 and " << MAX_MIN_ANG << ".\n";
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


    //Read in the number of repeating units per row and per column
    while(true){
	cout << "Enter the number of horizontal and vertical repeating units: ";
	getline(cin, response);
	num_read = sscanf(response.c_str(), "%d %d", &nx, &ny);
	if(num_read > 0 && nx > 0 && ny > 0) break;
	else cerr << "Enter two positive integers.\n";
    }

    //Prompt for scale factor
    scale = get_double("Enter the scale factor: ", FLT_MIN, FLT_MAX);

    //Import information from a network description file
    import_lattice(rules, nns, scale);

    //Now that the lattice has been specified and the numbers of repeating
    //horizontal and vertical units have been specified, define the bounds
    //of the space in which the lattice is to be formed
    cycle_len = 0;

    for(rule_index = 2; rule_index < rules[0].size(); rule_index ++){
        cycle_len += rules[0][rule_index];
    }

    //The maximum x coordinate is the number of horizontal repeating units
    //times the width of one such repeating unit, plus the horizontal offset
    //between the left bound of the nextwork and the left-most site in the
    //row
    max_x = FLT_MIN;
    for(vector<double> next_rule : rules){
	far_right = next_rule[0] + cycle_len*nx - *next_rule.rbegin();
	if(far_right > max_x) max_x = far_right;
    }

    //Prompt for bond portion
    portion = get_double("Enter the bond portion: ", 0, 1);

    //Generate the lattice
    make_edges(rules, nns, nx, ny, nx*cycle_len, max_x, points, edges);

    //If the portion of bonds to be retained is less than 1, carry out
    //dilution
    if(portion < 1){
        //Prepare for random number generation    
        ranfile = fopen("/dev/urandom", "r");
        size = fread(&seed, sizeof(unsigned), 1, ranfile);
        fclose(ranfile);

		if(size < 1){
			cerr << "Error: too few bits were read.\n";
		}

        sprintf(num_string, "%u", seed);
        setenv("GSL_RNG_SEED", num_string, 1);
        gsl_rng_env_setup();


        //Prompt for correlation strength
        correlation = get_double("Enter the correlation strength: ", 0, 1);

        //Prompt for the coordination number
        do{
            cout << "Enter the nominal coordination number: ";
            getline(cin, response);
            num_read = sscanf(response.c_str(), "%d", &maxz);
            if(num_read < 1) cerr << "Enter an integer.\n";
            else if(num_read < 1) cerr << "Enter a positive integer\n";
        }while(maxz < 1);

        //Use the protocol discussed in the preamble to this code to dilute
        //bonds in a spatially correlated manner
        correlated_dilution(edges, points, portion, correlation, maxz, make_a_mesh);
    }

    //Add hinges, if desired
    if(hinge) add_hinges(points, edges, nx*cycle_len);

    if(make_a_mesh){
	if(max_area == 0) max_area = .5 * scale * scale;
        network_to_mesh(points, nx*cycle_len, min_ang, max_area, facets, areas, lookup);
        replace_edges(edges, lookup);
    }

    if(get_triples){
        triples = generate_triples(points, edges);
    }

    //Prompt for file to report the network
    do{
        cout << "Enter the file name for the network: ";
        getline(cin, response);
    }while(response.compare("") == 0);

    report_file = fopen(response.c_str(), "w");

    //Report the created structures to an appropriate file
    if(! make_a_mesh){
        report_plf(report_file, points, edges, nx*cycle_len, triples);
    }
    else{
        report_mesh(report_file, points, edges, triples, facets, areas, nx*cycle_len);
    }

    fclose(report_file);

    return 0;
}

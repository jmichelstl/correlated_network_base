/*
 * This program produces a "kagomization" of a random point set. First, a
 * random Poisson packing is made, with a specified minimum separation, and a
 * specified set of bounds on x and y coordinates. Next, eight copies of this 
 * distribution are produced, and the nine instances of the packing are arranged
 * in a 3x3 grid. Next, a Voronoi diagram of the point set is produced. Finally,
 * lines connecting the midpoints of adjacent Voronoi edges that are contained 
 * within the prescribed x and y bounds are kept, with edges protruding beyond
 * the left and right boundaries wrapping around to the opposite side.
 */

#define REAL double

#include <tuple>
#include <set>
#include <algorithm>
#include <boost/polygon/voronoi.hpp>
#include <boost/polygon/voronoi_builder.hpp>
#include <boost/polygon/voronoi_diagram.hpp>
#include <boost/random.hpp>
#include <ctime>
#include <cstdint>
#include <map>
#include <stdio.h>
#include <string>
#include <iostream>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <set>

using namespace std;
using boost::polygon::voronoi_builder;
using boost::polygon::voronoi_diagram;

#define MAX_ATTEMPTS 30
#define MAX_INVALID_ENTRIES 10
#define READ_SUCCESS 0
#define READ_FAILURE -1
#define INT_FLOOR 0xF8000000
#define INT_CEIL 0x07FFFFFF
#define DIFLOOR -134217728
#define DIRANGE 268435455

#define min(a, b) a < b ? a : b
#define max(a, b) a > b ? a : b

#define ipair tuple<int, int>
#define itriple tuple<int, int, int>
#define edge_triple tuple<int, int, short>
#define TOL 1e-8

struct Point {
	double x, y;
	Point(double my_x, double my_y) : x(my_x), y(my_y) {}

	Point(){
		x = 0;
		y = 0;
	}

	bool operator == (const Point& point) const{
		return(abs(point.x - x) < TOL && abs(point.y - y) < TOL);
	}

	friend bool operator < (const Point &p1, const Point &p2){
		if(p2.y > p1.y + TOL) return true;
		else if(p1.y > p2.y + TOL) return false;
		else if(p2.x > p1.x + TOL) return true;
		else return false;
	}
};

/*
 * Define structures for Voronoi diagram construction
 */
struct IPoint {
	int x;
	int y;
	IPoint(int my_x, int my_y) : x(my_x), y(my_y) {}
};

struct Segment {
	IPoint p0;
	IPoint p1;

	Segment(int x1, int y1, int x2, int y2) : p0(x1, y1), p1(x2, y2) {}
};

namespace boost {
namespace polygon {

template <>
struct geometry_concept<IPoint> {
	typedef point_concept type;
};

template <>
struct point_traits<IPoint> {
	typedef int coordinate_type;

	static inline coordinate_type get(const IPoint& point, orientation_2d orient) {
		return (orient == HORIZONTAL) ? point.x : point.y;
	}
};

template <>
struct geometry_concept<Segment> {
	typedef segment_concept type;
};

template <>
struct segment_traits<Segment> {
	typedef int coordinate_type;
	typedef IPoint point_type;

	static inline point_type get(const Segment& segment, direction_1d dir) {
		return dir.to_int() ? segment.p1 : segment.p0;
	}
};
}
}

//Find the distance between two points
double distance(Point p1, Point p2){
	return sqrt((p2.x - p1.x)*(p2.x - p1.x) + (p2.y - p1.y)*(p2.y - p1.y));
}

//Given x and y bounds, and a minimum separation, create a Poisson disk packing.
//Create a grid whose cells have diagonal lengths equal to the minimum spacing.
//Seed the grid with an initial point, and add this point to a list of
//"active points". Pick an active point at random, and create a candidate
//neighbor separated from this point by an amount between the minimum distance
//and twice the minimum distance, either until the candidate is separated from
//all other points by the requisite distance, or the maximum number of attempts
//have been made. If no neighbor is successfully added, eject the current
//active point from the list.
vector<Point> poisson_packing(double minx, double maxx, double miny, double maxy, double min_spacing){

	int attempts, p_count, x_stride, y_stride, grindex;
	int x_idx, y_idx, randex, x_itr, y_itr, xoff_iter, yoff_iter;
	int minx_idx, miny_idx, maxx_idx, maxy_idx;
	double grid_len, x_range, y_range, rand_x, rand_y, rlen, rangle;
	double w_minx, w_maxx, w_miny, w_maxy, w_xrange, w_yrange;
	bool valid;

	//A list of all points created, and a list of points which can acquire
	//neighbors.
	vector<Point> poisson_points, active_points;

	Point generator, candidate, neighbor;

	//Structures generating random numbers
	time_t now = time(0);
	boost::random::lagged_fibonacci9689 dgen{static_cast<uint32_t>(now)};
	now = time(0);
	boost::random::mt19937 igen{static_cast<uint32_t>(now)};
	boost::random::uniform_int_distribution<> udist(0, 0);

	//Set up wide boundaries for the extended grid

	//Initialize the occupancy table. At first, set the entry for each cell
	//to -1, to indicate the cell is vacant.
	grid_len = min_spacing / sqrt(2);
	x_range = maxx - minx;
	y_range = maxy - miny;

	w_minx = minx - x_range;
	w_maxx = maxx + x_range;
	w_miny = miny - y_range;
	w_maxy = maxy + y_range;
	w_xrange = 3 * x_range;
	w_yrange = 3 * y_range;

	x_stride = (int) (w_xrange / grid_len) + 1;
	y_stride = (int) (w_yrange / grid_len) + 1;

	vector<int> occupancy_table(x_stride * y_stride, -1);

	//Seed the grid with an initial point
	rand_x = dgen() * x_range + minx;
	rand_y = dgen() * y_range + miny;
	candidate = Point(rand_x, rand_y);
	active_points.push_back(candidate);

	for(xoff_iter = -1; xoff_iter <= 1; xoff_iter ++){
		for(yoff_iter = -1; yoff_iter <=1; yoff_iter ++){
			x_idx = (int) ((rand_x+x_range*xoff_iter-w_minx)/grid_len);
			y_idx = (int) ((rand_y+y_range*yoff_iter-w_miny)/grid_len);
			occupancy_table[y_idx * x_stride + x_idx] = poisson_points.size();
			poisson_points.push_back(Point(rand_x+x_range*xoff_iter,rand_y+y_range*yoff_iter));
		}
	}

	//Add additional points until the list of active points is empty. An
	//empty list of active points indicates no more free grid cells are
	//available
	while(active_points.size() > 0){

		attempts = 0;
			
		//Reset the uniform distribution to set the range
		//according to the size of the list of active points
		boost::random::uniform_int_distribution<>::param_type new_params(0, active_points.size() - 1);
		udist.param(new_params);

		//Create a candidate point
		randex = udist(igen);
		generator = active_points[randex];

		while(attempts < MAX_ATTEMPTS){

			valid = true;

			do{
				rlen = grid_len * (dgen() + 1);
				rangle = 2 * M_PI * dgen();
				rand_x = generator.x + rlen * cos(rangle);
				rand_y = generator.y + rlen * sin(rangle);
			}while(!(rand_x>=minx && rand_x<=maxx && rand_y>=miny && rand_y<=maxy));
			x_idx = (int) ((rand_x - w_minx) / grid_len); 
			y_idx = (int) ((rand_y - w_miny) / grid_len);
			grindex = y_idx * x_stride + x_idx;
			if(occupancy_table[grindex] >= 0){
				attempts ++;
				continue;
			}

			candidate = Point(rand_x, rand_y);

			//Find the grid cell in which the candidate lies, and
			//check that grid cell and the neighboring grid cells
			//for points.
			minx_idx = max(x_idx - 2, 0);
			miny_idx = max(y_idx - 2, 0);
			maxx_idx = min(x_idx + 2, x_stride - 1);	
			maxy_idx = min(y_idx + 2, y_stride - 1);

			for(x_itr = minx_idx; x_itr <= maxx_idx; x_itr ++){
				for(y_itr = miny_idx;y_itr<= maxy_idx;y_itr ++){
					grindex = x_itr + y_itr * x_stride;
					if(occupancy_table[grindex] >= 0){
						neighbor = poisson_points[occupancy_table[grindex]];
						if(distance(candidate, neighbor) < min_spacing){
							valid = false;
							break;
						}
					}
				}
				if(! valid) break;
			}
			
			if(valid){
				for(xoff_iter = -1; xoff_iter <= 1; xoff_iter ++){
					for(yoff_iter = -1; yoff_iter <=1; yoff_iter ++){
						x_idx = (int) ((rand_x+x_range*xoff_iter-w_minx)/grid_len);
						y_idx = (int) ((rand_y+y_range*yoff_iter-w_miny)/grid_len);
						occupancy_table[y_idx * x_stride + x_idx] = poisson_points.size();
						poisson_points.push_back(Point(rand_x+x_range*xoff_iter,rand_y+y_range*yoff_iter));
					}
				}

				active_points.push_back(candidate);
				break;
			}

			else{
				attempts ++;
			}

		}

		if(attempts == MAX_ATTEMPTS){
			active_points.erase(active_points.begin()+randex);
		}
	}

	return poisson_points;
}

//Find the Voronoi diagram of a point set
void get_voronoi_edges(vector<IPoint> points, vector<Point> &vPoints, vector<ipair> &edge_data, double minx, double xrange, double miny, double yrange){

	Point p1, p2;
	vector<Point> v_points;
	map<Point, int> pmap;
	int p_index = 0, idx1, idx2;
	double x1, y1, x2, y2;
	vector<Segment> segments;

	//Create a Voronoi diagram of the point packing
	voronoi_diagram<double> vd;
        construct_voronoi(points.begin(), points.end(), segments.begin(), segments.end(), &vd);

	//Iterate over edges in the Voronoi diagram, and add a record of each
	//to data structures that track points and connectivity.
        for(voronoi_diagram<double>::const_edge_iterator it = vd.edges().begin(); it != vd.edges().end(); it++){
                if(it->is_primary() && it->is_finite()){

			//Obtain end points, and reverse the mapping of the
			//points used to generate the Voronoi diagram to the
			//integer grid.
			x1 = (it->vertex0()->x()-DIFLOOR)*xrange/DIRANGE + minx;
			y1 = (it->vertex0()->y()-DIFLOOR)*yrange/DIRANGE + miny;
			x2 = (it->vertex1()->x()-DIFLOOR)*xrange/DIRANGE + minx;
			y2 = (it->vertex1()->y()-DIFLOOR)*yrange/DIRANGE + miny;

			p1 = Point(x1, y1);
                        p2 = Point(x2, y2);

			//If either point has not been encountered, add it
			//to the map from points to indices, and assign it a
			//unique index.
			if(pmap.find(p1) == pmap.end()){
				pmap.insert(make_pair(p1, p_index++));
				vPoints.push_back(p1);

			}
			if(pmap.find(p2) == pmap.end()){
				pmap.insert(make_pair(p2, p_index++));
				vPoints.push_back(p2);
			}

			idx1 = pmap[p1];
			idx2 = pmap[p2];

			//Add a record of the edge
			if(idx1 < idx2) edge_data.push_back(make_tuple(idx1,idx2));
                }
        }

}

//Convert points with double precision floating point coordinates to points
//with integer coordinates. Rescale the range of the x and y coordinates to
//the range representable by single precision integers.
vector<IPoint> snap_to_grid(vector<Point> in, double xmin, double xrange, double ymin, double yrange){

	vector<IPoint> out_points;
	int int_x, int_y;

	for(Point in_point : in){
		int_x = (int) ((in_point.x - xmin) / xrange * (INT_CEIL-INT_FLOOR)) + INT_FLOOR;
		int_y = (int) ((in_point.y - ymin) / yrange * (INT_CEIL-INT_FLOOR)) + INT_FLOOR;
		out_points.push_back(IPoint(int_x, int_y));
	}

	return out_points;
}

bool in_bounds(double minx, double maxx, double miny, double maxy, double x, double y){
	return x >= minx && x <= maxx && y >= miny && y <= maxy;
};

//Find the edge along the rectangular bounding box which is intersected by a
//given segment
Point intersection(Point p1, Point p2, double minx, double maxx, double miny, double maxy){

	double x, slope;

	slope = abs(p1.x - p2.x) > TOL ? (p2.y - p1.y)/(p2.x - p1.x) : FLT_MAX;

	//If the y coordinate of the second point is below the lower y bound,
	//look for an intersection with the bottom edge. Otherwise, look for
	//an intersection with the top edge.
	if(p2.y < miny){
		x = (miny - p1.y) / slope + p1.x;
		if(x >= minx && x <= maxx) return Point(x, miny);
	}
	else{
		x = (maxy - p1.y) / slope + p1.x;
		if(x >= minx && x <= maxx) return Point(x, maxy);
	}
	
	//If there is not an intersection with the top or bottom edge, then,
	//if the x coordinate is to the left of the lower x bound, look for
	//an intersection with the left edge. Otherwise, look for an 
	//intersection with the right edge.
	if(x < minx){
		return Point(minx, p1.y + slope * (minx - p1.x));
	}
	else{
		return Point(maxx, p1.y + slope * (maxx - p1.x));
	}
}

//Sort neighbors of a central point according to the angles between the
//positive x axis and the edges from the central point to each of the neighbors.
void sort_neighbors(vector<Point> points, int center, vector<int> &indices){

	Point cPoint =  points[center];
	int next_index, sort_pos, dest, replace_pos;
	Point p1, p2;
	double ang1, ang2;

	for(sort_pos = 1; sort_pos < indices.size(); sort_pos ++){
		next_index = indices[sort_pos];
		p1 = points[next_index];
		ang1 = atan2(p1.x - cPoint.x, p1.y - cPoint.y);

		dest = 0;
		while(dest < sort_pos){
			p2 = points[indices[dest]];
			ang2 = atan2(p2.x - cPoint.x, p2.y - cPoint.y);
			if(ang1 < ang2) break;
			dest ++;
		}

		for(replace_pos = sort_pos; replace_pos > dest; replace_pos --){
			indices[replace_pos] = indices[replace_pos - 1];
		}

		indices[dest] = next_index;
	}
		
}

//Connect midpoints of adjacent Voronoi edges. Keep edges connecting midpoints
//of adjacent Voronoi edges if they are at least partly within the bounding
//box. If an edge extends below the bottom or top of the bounding box, truncate
//it at the point where it intersects the bounding box. If an edge extends
//beyond the left or right edge, use one of two polices:
//	-If both end points are within the y bounds, add a horizontal offset
//	to the end point that lies out of bounds, so that the point lies in
//	bounds.
//
//	-If one end point has a y coordinate that is out of bounds, truncate
//	the edge at the edge of the bounding box it intersects.
void kagomize_voronoi_diagram(vector<Point> v_points, vector<ipair> v_edges, vector<Point> &k_points, vector<edge_triple> &edge_data, double minx, double maxx, double miny, double maxy, double spacing){

	map<int, vector<int>> neighbor_lists;
	map<Point, int> pmap;
	int piter, niter, idx1, idx2, pIndex = 0;
	Point mid1, mid2;
	vector<int> neighbors;
	double x1, x2, y1, y2;
	size_t len;
	bool p1In, p2In;
	Point p1, p2;
	short offset;

	for(piter = 0; piter < v_points.size(); piter ++){
		neighbor_lists.insert(make_pair(piter, vector<int>()));
	};

	//Create data structure mapping points' indices to a list of indices
	//of neighboring points.
	for(ipair next_pair : v_edges){
		idx1 = get<0>(next_pair);
		idx2 = get<1>(next_pair);

		neighbor_lists[idx1].push_back(idx2);
		neighbor_lists[idx2].push_back(idx1);
	}

	//For each point, follow one of two policies:
	//	-If the point has three neighbors, connect each pair of 
	//	midpoints of edges connected to that point
	//
	//	-If the point has more than three neighbors, sort the neighbors
	//	according to the angle between the positive x axis and a line
	//	from the central point to the neighor. Then connect midpoint
	//	of each edge to the midpoint of the edge that comes next in
	//	cyclical order.
	for(piter = 0; piter < v_points.size(); piter ++){

		//Include any edge if both end points are in bounds. If just one
		//point is in bounds, apply the aforementioned protocol to
		//truncate the edge.
		
		if(neighbor_lists[piter].size() > 3){
			sort_neighbors(v_points, piter, neighbor_lists[piter]);
		}

		vector<int> curr_list = neighbor_lists[piter];
		len = neighbor_lists[piter].size();
		for(niter = 0; niter < len; niter++){
			x1 = (v_points[piter].x + v_points[curr_list[niter]].x) / 2;
			y1 = (v_points[piter].y + v_points[curr_list[niter]].y) / 2;
			x2 = (v_points[piter].x + v_points[curr_list[(niter+1)%len]].x)/2;
			y2 = (v_points[piter].y + v_points[curr_list[(niter+1)%len]].y)/2;

			p1In = in_bounds(minx, maxx, miny, maxy, x1, y1);
			p2In = in_bounds(minx, maxx, miny, maxy, x2, y2);

			if(! p1In && ! p2In) continue;

			if(p1In && p2In){
				p1 = Point(x1, y1);
				p2 = Point(x2, y2);
				offset = 0;
			}

			else{
				p1 = p1In ? Point(x1, y1) : Point(x2, y2);
				p2 = p1In ? Point(x2, y2) : Point(x1, y1);

				if(p2.y >= miny && p2.y <= maxy){
					//Only include a wraparound edge if
					//it wraps from the right edge back
					//to the left edge. Otherwise, there
					//will be duplicates.
					if(p2.x > maxx){
						p2.x = p2.x + minx - maxx;
						offset = 1;
					}

					else continue;
				}

				else{
					p2 = intersection(p1,p2,minx,maxx,miny,maxy);
					offset = 0;
				}
			}

			if(pmap.find(p1) == pmap.end()){
				pmap.insert(make_pair(p1,pIndex++));
				k_points.push_back(p1);
			}
			if(pmap.find(p2) == pmap.end()){
				pmap.insert(make_pair(p2,pIndex++));
				k_points.push_back(p2);
			}

			idx1 = pmap[p1];
			idx2 = pmap[p2];
			edge_data.push_back(make_tuple(idx1,idx2,offset));
		}
	}

}

int read_bounds(string prompt, double &min, double &max){

	int bad_entries = 0;
	string response;

        do{
                cout << prompt;
                getline(cin, response);
                if(sscanf(response.c_str(), "%lf %lf" ,&min, &max) == 2){
                        if(min < max) return READ_SUCCESS;

                        else{
                                cerr << "The min must be less than the max.\n";
                                bad_entries ++;
                        }
                }

                else{
                        cerr << "Enter two numbers.\n";
                        bad_entries ++;
                }

        }while(bad_entries < MAX_INVALID_ENTRIES);

	cerr << "The maximum number of invalid entries has been reached.\n";

	return READ_FAILURE;
}

/*
 * Remove a portion of the bonds in a network completely at random
 */
void simple_dilute(vector<Point> &out_pts, vector<edge_triple> &out_trpls, vector<itriple> &out_bts, double p){

	int to_keep, randex;
	map<int, int> change_map;
	int idx1, idx2, idx3;
	set<ipair> kept_edges;

	vector<Point> in_pts(out_pts.begin(), out_pts.end());
	vector<edge_triple> in_trpls(out_trpls.begin(), out_trpls.end());
	vector<itriple> in_bts(out_bts.begin(), out_bts.end());

	to_keep = p * in_trpls.size();

	//Structures generating random numbers
	time_t now = time(0);
	boost::random::mt19937 igen{static_cast<uint32_t>(now)};
	boost::random::uniform_int_distribution<> udist(0, in_trpls.size() - 1);

	while(out_trpls.size() < to_keep){

		//Choose a random edge from the initial list of edges for
		//inclusion.
		randex = udist(igen);
		idx1 = get<0>(in_trpls[randex]);
		idx2 = get<1>(in_trpls[randex]);

		if(change_map.find(idx1) == change_map.end()){
			change_map.insert(make_pair(idx1, out_pts.size()));
			out_pts.push_back(in_pts[idx1]);
		}

		if(change_map.find(idx2) == change_map.end()){
			change_map.insert(make_pair(idx2, out_pts.size()));
			out_pts.push_back(in_pts[idx2]);
		}

		out_trpls.push_back(make_tuple(change_map[idx1], change_map[idx2], get<2>(in_trpls[randex])));
		kept_edges.insert(make_tuple(min(change_map[idx1], change_map[idx2]), max(change_map[idx1], change_map[idx2])));

		in_trpls.erase(in_trpls.begin() + randex);
		
		//Reset the uniform distribution to set the range
		//according to the size of the list of unchosen edges
		boost::random::uniform_int_distribution<>::param_type new_params(0, in_trpls.size() - 1);
		udist.param(new_params);

	}

	//If triples are provided, cull those triples in which not all three
	//points have not been retained.
	out_bts.clear();
	for(itriple bt : in_bts){
		idx1 = get<0>(bt);
		idx2 = get<1>(bt);
		idx3 = get<2>(bt);

		if(change_map.find(idx1)==change_map.end()||change_map.find(idx2)==change_map.end()||change_map.find(idx3)==change_map.end()) continue;

		idx1 = change_map[idx1];
		idx2 = change_map[idx2];
		idx3 = change_map[idx3];

		if(kept_edges.find(make_pair(min(idx1,idx2), max(idx1,idx2)))!=kept_edges.end()&&kept_edges.find(make_pair(min(idx1,idx3), max(idx1,idx3)))!=kept_edges.end()){
			out_bts.push_back(make_tuple(idx1, idx2, idx3));
		}
	}
}

//Dilute the network with optional structural correlation and polarization
void cor_pol_dilute(vector<Point> &out_pts, vector<edge_triple> &out_trpls, vector<itriple> &out_bts, double p, double corr, double nx, double ny, double pstren){

	map<edge_triple, vector<edge_triple>> neigh_lists;
	map<edge_triple, int> neigh_cnts;
	map<int, vector<edge_triple>> point_to_edge;
	edge_triple e1, e2;
	double ret_prob, cosine, xdiff, ydiff;
	map<int, int> change_map;
	int to_keep, randex, idx1, idx2, idx3, iter, inner;
	Point p1, p2;
	set<ipair> kept_edges;

	vector<Point> in_pts(out_pts.begin(), out_pts.end());
	vector<edge_triple> in_trpls(out_trpls.begin(), out_trpls.end());
	vector<itriple> in_bts(out_bts.begin(), out_bts.end());

	//Structures generating random numbers
	time_t now = time(0);
	boost::random::mt19937 igen{static_cast<uint32_t>(now)};
	boost::random::uniform_int_distribution<> udist(0, in_trpls.size() - 1);
	now = time(0);
	boost::random::lagged_fibonacci9689 dgen{static_cast<uint32_t>(now)};

	to_keep = p * in_trpls.size();

	//Create neighbor lists from edges to the edges with which they share
	//vertices, and a map from each edge to the number of neighbors that
	//have been retained.
	if(corr > 0){

		//Initialize an auxilliary map from point indices to edges
		//which have those points as end points
		for(iter = 0; iter < in_pts.size(); iter++){
			point_to_edge.insert(make_pair(iter, vector<edge_triple>()));
		}

		for(edge_triple t : in_trpls){
			neigh_lists.insert(make_pair(t, vector<edge_triple>()));
			idx1 = get<0>(t);
			idx2 = get<1>(t);
			point_to_edge[idx1].push_back(t);
			point_to_edge[idx2].push_back(t);
			neigh_cnts.insert(make_pair(t, 0));
		}

		for(auto miter=point_to_edge.begin(); miter!=point_to_edge.end(); miter++){
			for(iter = 0; iter < miter->second.size(); iter++){
				e1 = miter->second[iter];
				for(inner = iter+1; inner < miter->second.size();inner ++){
					e2 = miter->second[inner];

					neigh_lists[e1].push_back(e2);
					neigh_lists[e2].push_back(e1);
				}
			}
		}

	}

	//Add edges with a probability (1 - c)^(6 - nn) * (|n . p|)^(2*e),
	//for correlation strength c, retained nearest neighbor bonds nn,
	//bond normal vector n, polarization normal vector p, and polarization
	//strength e.
	out_pts.clear();
	out_trpls.clear();
	do{
		do{
			randex = udist(igen);
			ret_prob = 1;
			e1 = in_trpls[randex];
			idx1 = get<0>(e1);
			idx2 = get<1>(e1);
			p1 = in_pts[idx1];
			p2 = in_pts[idx2];

			//If correlated dilution is desired, find the number
			//of nearest neighbors of the candidate edge and modify
			//the retention probability accordingly.
			if(corr > 0){
				ret_prob *= pow(1 - corr, max(0,6-neigh_cnts[e1]));
			}

			//If polarization is desired, calculate the cosine
			//between the candidate bond and the polarization
			//direction, and modify the retention probality.
			if(pstren > 0){
				xdiff = p2.x - p1.x;
				ydiff = p2.y - p1.y;
				cosine = abs((xdiff*nx+ydiff*ny)/sqrt(xdiff*xdiff+ydiff*ydiff));
				ret_prob *= pow(cosine, 2*pstren);
			}

		}while(dgen() >= ret_prob);

		//If a random number on [0,1] is less than the retention
		//probability, keep the candidate bond. Update the
		//nearest neighbor counts for the neighbors of the kept
		//bond, add any end point to the list of output points
		//if it has not been added, and remove the candidate
		//bond from the list of candidates, while adding it
		//to the list of output bonds.
		if(change_map.find(idx1) == change_map.end()){
			change_map.insert(make_pair(idx1, out_pts.size()));
			out_pts.push_back(p1);
		}

		if(change_map.find(idx2) == change_map.end()){
			change_map.insert(make_pair(idx2, out_pts.size()));
			out_pts.push_back(p2);
		}

		for(edge_triple t : neigh_lists[e1]){
			neigh_cnts[t] = neigh_cnts[t] + 1;
		}

		out_trpls.push_back(make_tuple(change_map[idx1],change_map[idx2],get<2>(e1)));
		kept_edges.insert(make_tuple(min(change_map[idx1], change_map[idx2]), max(change_map[idx1], change_map[idx2])));
		in_trpls.erase(in_trpls.begin() + randex);

		boost::random::uniform_int_distribution<>::param_type new_params(0, in_trpls.size() - 1);
		udist.param(new_params);
	}while(out_trpls.size() < to_keep);

	//If triples are provided, cull those triples in which not all three
	//points have not been retained.
	out_bts.clear();
	for(itriple bt : in_bts){
		idx1 = get<0>(bt);
		idx2 = get<1>(bt);
		idx3 = get<2>(bt);

		if(change_map.find(idx1)==change_map.end()||change_map.find(idx2)==change_map.end()||change_map.find(idx3)==change_map.end()) continue;

		idx1 = change_map[idx1];
		idx2 = change_map[idx2];
		idx3 = change_map[idx3];

		if(kept_edges.find(make_pair(min(idx1,idx2), max(idx1,idx2)))!=kept_edges.end()&&kept_edges.find(make_pair(min(idx1,idx3), max(idx1,idx3)))!=kept_edges.end()){
			out_bts.push_back(make_tuple(idx1, idx2, idx3));
		}
	}
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

//Find bond triples for calculating bending terms. Where edges converge on
//a vertex, pair off edges such that the geometric mean of the angles between
//paired angles is maximized.
vector<itriple> generate_triples(vector<Point> pts, vector<edge_triple> edges){

	map<int, vector<int>> neigh_map;
	int idx1, idx2, niter;
	double gmean, highmean;
	vector<int> optimum;
	vector<itriple> btriples;

	for(idx1 = 0; idx1 < pts.size(); idx1++){
		neigh_map.insert(make_pair(idx1, vector<int>()));
	}

	//Make a mapping from each vertex index to a list of the indices of
	//its neighbors.
	for(edge_triple t : edges){
		idx1 = get<0>(t);
		idx2 = get<1>(t);
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

//Prompt a user for x and y bounds, and a minimum spacing. With these 
//parameters, create a random point set, find the Voronoi diagram for this 
//random point set, and from there produce a "kagomized" network, to be reported
//as a periodic lattice file (.plf).
int main(int argc, char **argv){

	//Parameters for creating a random point packing
	double minx, miny, maxx, maxy, xrange, yrange, spacing, bondp;
	double corr, nx, ny, pstren, len;
	//Initial random point set to be kagomized
	vector<Point> poisson_points, voronoi_points, kag_points;
	//Points snapped to an integer grid
	vector<IPoint> grid_points;
	//Integer index pairs to represent Voronoi edges
	vector<ipair> vor_edges;
	//List of data structures representing eges for a periodic lattice file,
	//with end points and, if needed, a horizontal offset
	vector<edge_triple> periodic_edges;
	//File for reporting the periodic lattice file
	FILE *output;
	string response, file_name;
	int bad_entries;
	bool get_triples = false;
	char c;
	vector<itriple> bond_triples;

	//Process an optional command line argument indicating triples should
	//be written.
	while((c = getopt(argc, argv, "t")) != -1){
		switch(c) {
			case 't':
				get_triples = true;
				break;
			case '?' :
	        	        if(isprint(optopt)){
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

	//Prompt for x and y bounds
	if(read_bounds("Enter minumum and maximum x bounds: ",minx,maxx) == READ_FAILURE){
		cerr << "The program will now exit.\n";
		return READ_FAILURE;
	}

	if(read_bounds("Enter minumum and maximum y bounds: ",miny,maxy) == READ_FAILURE){
		cerr << "The program will now exit.\n";
		return READ_FAILURE;
	}

	//Prompt for minimum separation
	bad_entries = 0;
	do{
		cout << "Enter the minimum separation: ";
		getline(cin, response);
		if(sscanf(response.c_str(), "%lf", &spacing) == 1){
			if(spacing > 0) break;

			else{
				cerr << "Enter a positive number.\n";
				bad_entries ++;
			}
		}

		else{
			cerr << "Enter a positive number.\n";
			bad_entries ++;
		}

		cout << "Running.\n";
	}while(bad_entries < MAX_INVALID_ENTRIES);
	if(bad_entries == MAX_INVALID_ENTRIES){
		cerr << "No separation could be read. The program will exit.\n";
		return -1;
	}

	//Prompt for the retained portion of bonds. If the portion is less than
	//1, offer the option for correlation. If correlation is chosen, prompt
	//for the correlation strength.
	
	//Offer the option to add polarization. If polarizationis desired,
	//prompt for the preferred direction and strength of polarization.
	
	//Prompt for the name of the file to which the network should be written
	bad_entries = 0;
	do{
		cout << "Enter a file name: ";
		getline(cin, file_name);
		if(! file_name.compare("") == 0) break;
		else{
			cerr << "Enter a name.\n";
			bad_entries ++;
		}
	}while(bad_entries < MAX_INVALID_ENTRIES);
	if(bad_entries == MAX_INVALID_ENTRIES){
		cerr << "No file name could be read. The program will exit.\n";
		return -1;
	}

	//Offer the opportunity for random dilution. If dilution is desired,
	//prompt for bond portion, dilution correlation strength, polarization
	//direction, and polarization strength.
	bad_entries = 0;
	do{
		cout << "Enter the portion of bonds to keep: ";
		getline(cin, response);
		if(sscanf(response.c_str(), "%lf", &bondp) == 1){
			if(bondp >= 0 && bondp <= 1){
				break;
			}
		}

		cerr << "Enter a number between 0 and 1.\n";
		bad_entries ++;

	}while(bad_entries < MAX_INVALID_ENTRIES);
	if(bad_entries == MAX_INVALID_ENTRIES){
		cerr << "No valid entries could be read.\n";
		exit(-1);
	}

	if(bondp < 1){
		//Get correlation strength
		bad_entries = 0;
		do{
			cerr << "Enter the correlation strength: ";
			getline(cin, response);
			if(sscanf(response.c_str(), "%lf", &corr) == 1){
				if(corr >= 0 && corr <= 1){
					break;
				}
			}

			bad_entries ++;
			cerr << "Enter a number between 0 and 1.\n";
		}while(bad_entries < MAX_INVALID_ENTRIES);
		if(bad_entries == MAX_INVALID_ENTRIES){
			cerr << "No valid respone could be read.\n";
			exit(-1);
		}

		//Get polarization strength
		bad_entries = 0;
		do{
			cerr << "Enter the polarization strength: ";
			getline(cin, response);
			if(sscanf(response.c_str(), "%lf", &pstren) == 1){
				break;
			}

			bad_entries ++;
			cerr << "Enter a number between 0 and 1.\n";
		}while(bad_entries < MAX_INVALID_ENTRIES);
		
		//Get the polarization direction
		bad_entries = 0;
		do{
			cout << "Enter the polarization vector components: ";
			getline(cin, response);
			if(sscanf(response.c_str(), "%lf %lf", &nx, &ny) == 2){
				if(!(nx == 0 && ny == 0)){
					len = sqrt(nx*nx + ny*ny);
					nx /= len;
					ny /= len;
					break;
				}

				else{
					cerr << "At least one component must be non-zero\n";
				}
			}

			cerr << "Enter two real numbers, not both zero.\n";

			bad_entries ++;
		}while(bad_entries < MAX_INVALID_ENTRIES);
	}

	xrange = maxx - minx;
	yrange = maxy - miny;

	//Construct the random point set, and produce a 3x3 grid with 9 copies
	poisson_points = poisson_packing(minx, maxx, miny, maxy, spacing);
	//poisson_points = duplicate_grid(poisson_points, xrange, yrange);

	/*output = fopen("ppacking.txt", "w");
	for(Point p : poisson_points){
		fprintf(output, "%2.12lf\t%2.12lf\n", p.x, p.y);
	}
	fclose(output);*/

	//Obtain the Voronoi diagram
	grid_points = snap_to_grid(poisson_points, minx - xrange, 3*xrange, miny - yrange, 3*yrange);

	/*output = fopen("ipoints.txt", "w");
	for(IPoint p : grid_points){
		fprintf(output, "%d\t%d\n", p.x, p.y);
	}
	fclose(output);*/

	get_voronoi_edges(grid_points, voronoi_points, vor_edges, minx - xrange, 3*xrange, miny - yrange, 3*yrange);

	/*output = fopen("vpoints.txt", "w");
	for(Point p : voronoi_points){
		fprintf(output, "%2.12lf\t%2.12lf\n", p.x, p.y);
	}
	fclose(output);

	output = fopen("vedges.txt", "w");
	for(tuple<int,int> t : vor_edges){
		fprintf(output, "%d\t%d\n", get<0>(t), get<1>(t));
	}
	fclose(output);*/

	//Kagomize the point set
	kagomize_voronoi_diagram(voronoi_points, vor_edges, kag_points, periodic_edges, minx, maxx, miny, maxy, spacing);

	//If desired, generate bending triples
	if(get_triples){
		bond_triples = generate_triples(kag_points, periodic_edges);
	}

	//If the portion of bonds to be retained is less than 1, dilute the
	//network. Pass as arguments copies of the initial kagomized network's 
	//points and edges, and overwrite these arrays with the points and 
	//edges of the diluted network.
	if(bondp < 1){
		if(corr > 0 || pstren > 0){
			cor_pol_dilute(kag_points, periodic_edges, bond_triples, bondp, corr, nx, ny, pstren);
		}

		else{
			simple_dilute(kag_points, periodic_edges, bond_triples, bondp);
		}
	}

	//Write the newly created network to a periodic lattice file.
	output = fopen(file_name.c_str(), "w");

	if(output != NULL){
		
		fprintf(output, "%ld\t%ld\t%2.12lf", kag_points.size(), periodic_edges.size(), maxx - minx);
	
		if(get_triples) fprintf(output, "\t%ld", bond_triples.size());

		fprintf(output, "\n");

		for(Point p : kag_points){
			fprintf(output, "%2.12lf\t%2.12lf\n", p.x, p.y);
		}

		for(edge_triple et : periodic_edges){
			fprintf(output, "%d\t%d\t%hd\n", get<0>(et), get<1>(et), get<2>(et));
		}

		if(get_triples){
			for(itriple it : bond_triples){
				fprintf(output, "%d\t%d\t%d\n", get<0>(it), get<1>(it), get<2>(it));
			}
		}

		fclose(output);

		return 0;
	}

	else{
		cerr << "The report file could not be opened.\n";
		return -1;
	}
}

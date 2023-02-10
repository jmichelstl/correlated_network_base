/*
 * This header file handles the creation of triangular finite element
 * meshes with an embedded network. The meshing is accomplished with
 * Triangle, copyrighted by Jonathan Richard Shewchuck. Any use of
 * this software should retain this attribution, and any commercial
 * deployment requires express permission from Jonathan Richard Shewchuck.
 */

#define VOID int
#define REAL double
#define COMP_TOL 1e-10

#include "triangle.hpp"
#include "network_utils.h"
#include <vector>
#include <tuple>
#include <string.h>
#include <stdlib.h>
//#include <unordered_set>
#include <set>
#include <math.h>

//Add boundary points to the point set such that the mesh will be suitable for
//Lee-Edwards boundary conditions, the convex hull will be a rectangle, and
//boundary points are suitable closely spaced to ensure that a quality
//triangular mesh can be formed without adding Steiner points on the convex 
//hull.
void add_boundary_points(std::vector<Point> &points, double offset, double min_ang, double max_area){

    std::set<Point> p_set, left, right, bottom, top;
    double xmin=FLT_MAX, xmax=FLT_MIN, ymin = FLT_MAX, ymax = FLT_MIN, spacing;
    double low, high, inc, coord;
    Point new_point;

    //Make a record of points, and find those points that lie on the left,
    //bottom, and top edges
    for(Point p : points){

	p_set.insert(p);

	//Update the minimal x coordinate and the left edge
        if(p.x < xmin - COMP_TOL){
            xmin = p.x;
	    left.clear();
	    left.insert(p);
	}
	else if(p.x < xmin + COMP_TOL){
            left.insert(p);
	}

	//Update the maximal x coordinate and the right edge
        if(p.x > xmax + COMP_TOL){
            xmax = p.x;
	    right.clear();
	    right.insert(p);
	}
	else if(p.x > xmax - COMP_TOL){
            right.insert(p);
	}

	//Update the minimal y coordinate and the bottom edge
        if(p.y < ymin - COMP_TOL){
            ymin = p.y;
	    bottom.clear();
	    bottom.insert(p);
	}
	else if(p.y < ymin + COMP_TOL){
            bottom.insert(p);
	}

	//Update the maximal y coordinate and the top edge
        if(p.y > ymax + COMP_TOL){
            ymax = p.y;
	    top.clear();
	    top.insert(p);
	}
	else if(p.y > ymax - COMP_TOL){
            top.insert(p);
	}
    }

    //Add points to the boundaries so that boundary points along the right edge
    //are simply boundary points along the left edge, translated to the right
    //by a specified horizonal offset, and such that boundary points are
    //suitably closely spaced to ensure a quality triangulation can be
    //produced without adding additional boundary points.

    //Look for the corners, and add any corner that is missing
    if(p_set.find(Point(xmin, ymin)) == p_set.end()){
        new_point = Point(xmin, ymin);
        points.push_back(new_point);
	p_set.insert(new_point);
	left.insert(new_point);
	bottom.insert(new_point);
    }

    if(p_set.find(Point(xmin + offset, ymin)) == p_set.end()){
        new_point = Point(xmin + offset, ymin);
        points.push_back(new_point);
	p_set.insert(new_point);
	bottom.insert(new_point);
    }

    if(p_set.find(Point(xmin, ymax)) == p_set.end()){
        new_point = Point(xmin, ymax);
        points.push_back(new_point);
	p_set.insert(new_point);
	left.insert(new_point);
	top.insert(new_point);
    }

    if(p_set.find(Point(xmin + offset, ymax)) == p_set.end()){
        new_point = Point(xmin + offset, ymax);
        points.push_back(new_point);
	p_set.insert(new_point);
	top.insert(new_point);
    }

    //Find the spacing between successive supplementary points
    spacing = 2*sqrt(max_area * tan(min_ang * M_PI / 360));

    //Insert left-hand points, and their translated versions on the right
    //edge
    for(auto iter = left.begin(); iter != prev(left.end()); iter++){
        low = (*iter).y;
	high = (*(next(iter))).y;
	inc = (high - low) / floor((high - low) / spacing);
	for(coord = low + inc; coord < high; coord += inc){

            new_point = Point(xmin, coord);
	    if(p_set.find(new_point) == p_set.end()){
	        p_set.insert(new_point);
                points.push_back(new_point);
            }
	    
	    new_point = Point(xmin + offset, coord);
	    if(p_set.find(new_point) == p_set.end()){
	        p_set.insert(new_point);
	        points.push_back(new_point);
	    }
	}

	new_point = Point(xmin + offset, high);
	if(p_set.find(new_point) == p_set.end()){
	    p_set.insert(new_point);
	    points.push_back(new_point);
	}
    }

    //If a point on the right edge has no counterpart along the left edge,
    //add such a counterpart.
    if(xmax > xmin + offset - COMP_TOL){
        for(auto iter = right.begin(); iter != right.end(); iter ++){
            new_point = Point(xmin, (*iter).y);
	    if(p_set.find(new_point) == p_set.end()){
                p_set.insert(new_point);
		points.push_back(new_point);
	    }
        }
    }

    //Insert bottom points
    for(auto iter = bottom.begin(); iter != prev(bottom.end()); iter++){
        low = (*iter).x;
	high = (*(next(iter))).x;
	inc = (high - low) / floor((high - low) / spacing);
	for(coord = low + inc; coord < high; coord += inc){
            new_point = Point(coord, ymin);
            if(p_set.find(new_point) == p_set.end()){
                points.push_back(new_point);
            }
	}
    }

    //Insert top points
    for(auto iter = top.begin(); iter != prev(top.end()); iter++){
        low = (*iter).x;
	high = (*(next(iter))).x;
	inc = (high - low) / floor((high - low) / spacing);
	for(coord = low + inc; coord < high; coord += inc){
            new_point = Point(coord, ymax);
            if(p_set.find(new_point) == p_set.end()){
                points.push_back(new_point);
            }
	}
    }
}

//Given a set of points and mesh quality constraints, create a conforming
//Delaunay triangulation of a point set.
void make_network_mesh(std::vector<Point> in_pts, std::vector<Point> &out_pts, std::vector<tuple<int, int, int>> &facets, double offset, double min_ang, double max_area){

    char flags[31];
    struct triangulateio in, out;
    int idx;

    //Prepare flags for triangulation
    //sprintf(flags, "pBzq%10.8lfa%10.8lfQ", min_ang, max_area);
    sprintf(flags, "pcBzq%10.8lfa%10.8lfQY", min_ang, max_area);

    out_pts.clear();
    facets.clear();

    //Prepare input triangulation struct to describe the point set to be
    //triangulated.
    add_boundary_points(in_pts, offset, min_ang, max_area);
    in.numberofpoints = in_pts.size();
    in.numberofpointattributes = 0;
    in.pointlist = (double *) malloc(in.numberofpoints * 2 * sizeof(double));

    for(idx = 0; idx < in_pts.size(); idx ++){
        in.pointlist[2 * idx] = in_pts[idx].x;
	in.pointlist[2 * idx + 1] = in_pts[idx].y;
    }

    in.pointmarkerlist = (int *) NULL;
    in.numberofsegments = 0;
    in.numberofholes = 0;
    in.numberofregions = 0;

    //Prepare the triangulation data structure for the result of
    //triangulation.
    out.pointlist = (double *) NULL;
    out.pointmarkerlist = (int *) NULL;
    out.trianglelist = (int *) NULL;
    out.numberofpointattributes = 0;
    out.segmentlist = (int *) NULL;
    out.segmentmarkerlist = (int *) NULL;

    //Prepare the triangulation
    triangulate(flags, &in, &out, NULL);

    //Fill the lists of output points and triangles
    for(idx = 0; idx < out.numberofpoints; idx ++){
        out_pts.push_back(Point(out.pointlist[idx*2], out.pointlist[idx*2+1]));
    }

    for(idx = 0; idx < out.numberoftriangles; idx++){
        facets.push_back(make_tuple(out.trianglelist[idx*3],out.trianglelist[idx*3+1],out.trianglelist[idx*3+2]));
    }
   
    free(in.pointlist);
    free(out.pointlist);
    free(out.trianglelist);
    free(out.segmentlist);
}

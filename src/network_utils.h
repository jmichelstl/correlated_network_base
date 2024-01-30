#ifndef NET_UTILS_H
#define NET_UTILS_H

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <cmath>
#include <set>
#include <cfloat>
#include <unordered_map>

using namespace std;


#define PRIME1 (int) 1021297
#define PRIME2 (int) 1021651 //Primes for computing hash functions
#define FLOAT_TOL 1E-8 //Tolerance for deeming floating point numbers different
#define BIG_SLOPE 1E10
#define INF 1E12

struct Point{
    Point() {
        this->tol = 0.0000001;
    }

    Point(double xval, double yval) : x(xval), y(yval){
        this->tol = 0.0000001;
    }

    Point(double xval, double yval, int i) : x(xval), y(yval), index(i){
        this->tol = 0.0000001;
    }

    double x,y;
    double tol;
    int index;

    bool operator == (const Point& point) const{
        return(abs(point.x - x) < tol && abs(point.y - y) < tol);
    }

    friend bool operator < (const Point &p1, const Point &p2){
        if(p2.y > p1.y + p1.tol) return true;
        else if(p1.y > p2.y + p2.tol) return false;
        else if(p2.x > p1.x + p1.tol) return true;
        else return false;
    }
};

struct Edge{
    Edge() {}

    Edge(Point start, Point end) : p1(start), p2(end){
        top = p1.y > p2.y ? p1.y : p2.y;
        bottom = p1.y < p2.y ? p1.y : p2.y;
        left = p1.x < p2.x ? p1.x : p2.x;
        right = p1.x > p2.x ? p1.x : p2.x;
    }

    Point p1, p2;
    double top, bottom, left, right;

    void reset(Point newp1, Point newp2){
        p1 = newp1;
        p2 = newp2;

        top = p1.y > p2.y ? p1.y : p2.y;
        bottom = p1.y < p2.y ? p1.y : p2.y;
        left = p1.x < p2.x ? p1.x : p2.x;
        right = p1.x > p2.x ? p1.x : p1.x;
    }

    bool operator == (const Edge& edge) const {
        return((edge.p1==p1 && edge.p2==p2)||(edge.p2==p1 && edge.p1==p2));
    }

    bool operator < (const Edge& edge) const {
        if(top > edge.top) return false;
        else if(top == edge.top && right > edge.right) return false;
        else return true;
    }

    double length(){
        return sqrt((p2.x - p1.x)*(p2.x - p1.x) + (p2.y - p1.y)*(p2.y - p1.y));
    }

    Point midpoint(){
        return Point((p1.x + p2.x) / 2, (p1.y + p2.y) / 2);
    }
};

bool intersection(Edge e1, Edge e2){
    double slope1, slope2, x;
    Point p1 = e1.p1, p2 = e1.p2, p3 = e2.p1, p4 = e2.p2;
    bool inxrange1, inxrange2;

    slope1 = p1.x == p2.x ? BIG_SLOPE : (p2.y - p1.y)/(p2.x - p1.x);
    slope2 = p3.x == p4.x ? BIG_SLOPE : (p4.y - p3.y)/(p4.x - p3.x);

    x = (p1.y - p3.y + slope2*p3.x - slope1*p1.x)/(slope2 - slope1);

    inxrange1 = x >= min(p1.x,p2.x) && x <= max(p1.x, p2.x);
    inxrange2 = x >= min(p3.x,p4.x) && x <= max(p3.x, p4.x);

    return inxrange1 && inxrange2;
}

bool intersect_exclude_end(Edge e1, Edge e2){
    double slope1, slope2, x;
    Point p1 = e1.p1, p2 = e1.p2, p3 = e2.p1, p4 = e2.p2;
    bool inxrange1, inxrange2;

    slope1 = p1.x == p2.x ? BIG_SLOPE : (p2.y - p1.y)/(p2.x - p1.x);
    slope2 = p3.x == p4.x ? BIG_SLOPE : (p4.y - p3.y)/(p4.x - p3.x);

    x = (p1.y - p3.y + slope2*p3.x - slope1*p1.x)/(slope2 - slope1);

    inxrange1 = x > min(p1.x,p2.x) && x < max(p1.x, p2.x);
    inxrange2 = x > min(p3.x,p4.x) && x < max(p3.x, p4.x);

    return inxrange1 && inxrange2;
}

struct PolyData{

    PolyData(vector<Edge> edges) : polyEdges(edges) {
        left = INF;
        right = -INF;
        low = INF;
        high = -INF;

        for(Edge e : polyEdges){
            //See either end point is furthest to the right or left
            left = e.p1.x < left ? e.p1.x : left;
            right = e.p1.x > right ? e.p1.x : right;
            left = e.p2.x < left ? e.p2.x : left;
            right = e.p2.x > right ? e.p2.x : right;

            //See if either end point is lowest or highest
            low = e.p1.y < low ? e.p1.y : low;
            high = e.p1.y > high ? e.p1.y : high;
            low = e.p2.y < low ? e.p2.y : low;
            high = e.p2.y > high ? e.p2.y : high;
        }
    }

    double left, right, low, high;
    vector<Edge> polyEdges;

    bool operator < (const PolyData& pdat) const {
        if(high < pdat.low || high < pdat.high || low < pdat.low) return true;
        else return false;
    }

    bool isLeft(Point p){
        return p.x < left;
    }

    bool isRight(Point p){
        return p.x > right;
    }

    bool isBelow(Point p){
        return p.y < low;
    }

    bool isAbove(Point p){
        return p.y > high;
    }

    bool contains(Point p){
        bool inpoly = false;
        Point p1, p2;
        double x = p.x, y = p.y;

        if(isLeft(p) || isRight(p) || isBelow(p) || isAbove(p)) return false;

        for(Edge next : polyEdges){
            p1 = next.p1;
            p2 = next.p2;
            if(p1.y < y && p2.y >= y || p2.y < y && p1.y >= y){
                if(p1.x + (y - p1.y)/(p2.y-p1.y)*(p2.x-p1.x) < x){
                    inpoly = !inpoly;
                }
            }
        }

        return inpoly;
    }

    //First see if a point lies on any line segment of the polygon
    bool contains_permissive(Point p){
        double det, minx, maxx;

        for(Edge next : polyEdges){
            minx = min(next.p1.x, next.p2.x);
            maxx = max(next.p1.x, next.p2.x);
            det = (next.p1.x - p.x)*(next.p2.y - p.y);
            det -= (next.p1.y - p.y)*(next.p2.x  - p.x);
            if(abs(det)<FLOAT_TOL && p.x< maxx+FLOAT_TOL && p.x>minx-FLOAT_TOL){
                return true;
            }
        }

        return this->contains(p);
    }

};

namespace std {

    template<> struct hash<Point>{
        typedef size_t result_type;
        typedef Point argument_type;

        size_t operator() (const Point& p) const;
    };

    size_t hash<Point>::operator() (const Point& p) const {
        return (size_t) ((p.x+p.y)*(p.x+p.y)*(p.x*p.y+p.x*p.x+p.y*p.y*PRIME1))^PRIME2;
    }

    template<> struct hash<Edge>{
        typedef size_t result_type;
        typedef Edge argument_type;

        size_t operator() (const Edge& e) const;
    };

    size_t hash<Edge>::operator() (const Edge& e) const{
        size_t h1, h2;
        Point p1 = e.p1, p2 = e.p2;

        h1 = (size_t) ((p1.x+p1.y)*(p1.x+p1.y)*(p1.x*p1.y+p1.x*p1.x+p1.y*p1.y*PRIME1))^PRIME2;
        h2 = (size_t) ((p2.x+p2.y)*(p2.x+p2.y)*(p2.x*p2.y+p2.x*p2.x+p2.y*p2.y*PRIME1))^PRIME2;

        return h1*h1 + h2*h2;
    }
/*
    size_t hash<PolyData>::operator() (const PolyData& pd) const{
        size_t hashval = 1;

        for(Point p : pd.getPoints()){
            hashval *= hash<Point>(p);
        }
        return hashval;
    }
*/
}


bool yesno(string message){
    string response;

    cout << message << "(y/n): ";

    while(true){
        getline(cin, response);
        if(!response.compare("y")) return true;
        else if(!response.compare("n")) return false;
        else printf("Respond with \"y\" or \"n\": ");
    }
}

vector<string> split(string input, char delim){
    vector<string> result;
    size_t start, iter, len;
    bool reading_token = false;

    for(iter = 0; iter < input.size(); iter++){
        if(delim != input[iter]){
            if(!reading_token){
                reading_token = true;
                start = iter;
            }
        }

        else{
            if(reading_token){
                reading_token = false;
                result.push_back(input.substr(start, iter - start));
            }
        }
    }

    if(reading_token) result.push_back(input.substr(start, iter - start));


    return result;
}

/*
vector<string> split(string input, char delim){
    vector<string> result;
    size_t start, curr, len;

    curr = 0;
    do{
        len = 0;
        start = curr;
        while(input[curr] != delim && curr < input.size()){
            curr++;
            len ++;
        }
        result.push_back(input.substr(start, len));
        curr ++;
    }while(curr < input.size());

    return result;
}
*/
vector<double> parse_doubles(vector<string> numstring){
    double nextnum;
    vector<double> result;

    for(string s : numstring){
        if(sscanf(s.c_str(), "%lf", &nextnum)) result.push_back(nextnum);
    }

    return result;
}

vector<double> getdoubles(string prompt){
    double nextnum;
    vector<double> result;
    string response;

    printf("%s", prompt.c_str());
    getline(cin, response);

    return parse_doubles(split(response, ' '));
}

double get_double(string prompt, double min, double max){

    double result;
    string response;
    int num_read;

    do{
        cout << prompt;
        getline(cin, response);
        num_read = sscanf(response.c_str(), "%lf", &result);
        if(num_read == 0){
            cerr << "Enter a number between " << min << " and " << max << "\n";
        }
        else if(result < min || result > max){
            cerr << "The number must be between "<<min<< " and " << max << "\n";
        }
        else return result;

    }while(true);

}

void makeunion(vector<int>& setvec, int root1, int root2){
    if(setvec[root2] < setvec[root1]) setvec[root1] = root2;
    else{
        if(setvec[root1] == setvec[root2]) setvec[root1]--;

        setvec[root2] = root1;
    }
}

int find_root(vector<int> setvec, int elem){
    int pos = elem;

    if(elem >= setvec.size()) return -1;

    while(setvec[pos] >= 0) pos = setvec[pos];

    return pos;
}

char get_choice(string message, map<char, string> c_map, vector<char> choices){

    bool valid = false;
    string response;
    char choice;

    do{
        cout << message;
        for(auto it = c_map.begin(); it != c_map.end(); it++){
            cout << it->first << ": " << it->second << "\n";
        }
        getline(cin, response);
        if(response.compare("") == 0){
            cerr << "Please make a choice.\n";
        }

        else{
            choice = response[0];
            for(char next_choice : choices){
                if(choice == next_choice) return choice;
            }
            cerr << "Choice \"" << choice << "\"" << "not recognized.\n";
        }
    }while(!valid);

    return ' ';
}

//Functions to sort edges prioritizing alternatively the highest or lowest
//point of an edge

bool sort_by_bottom(Edge e1, Edge e2){
    return e1.bottom < e2.bottom;
}

bool sort_by_top(Edge e1, Edge e2){
    return e1.top < e2.top;
}

string enter_decline(string message){
    string response;
    cout << message << ", or return to decline: ";
    getline(cin, response);
    return response;
}

void open_dat_file(string prompt, ifstream& file_stream){
    string response;

    do{
        cout << prompt;
        getline(cin, response);
        file_stream.open(response);
        if(! file_stream.is_open()){
            if(! yesno("The file could not be read. Try again?")) break;
        }
    }while(! file_stream.is_open());
}

void open_output_file(string prompt, ofstream& file_stream){
    string response;

    do{
        cout << prompt;
        getline(cin, response);
        file_stream.open(response, ofstream::out);
        if(! file_stream.is_open()){
            if(! yesno("The file could not be read. Try again?")) break;
        }
    }while(! file_stream.is_open());
}

//Given a set of edges, get the underlying set of points
set<Point> get_points(vector<Edge> edges){
    set<Point> points;

    for(Edge next : edges){
        points.insert(next.p1);
        points.insert(next.p2);
    }

    return points;
}

//Perform a uniform displacement of all edges in a set
void displace(vector<Edge> &edges, double xdisp, double ydisp){
    int iter;
    set<Point> edge_points;
    unordered_map<Point, Point> replace_map;
    Point p1, p2;

    edge_points = get_points(edges);
    for(Point p : edge_points){
       replace_map.insert(make_pair(p, Point(p.x + xdisp, p.y + ydisp)));
    }
    //delete &edge_points;

    for(iter = 0; iter < edges.size(); iter++){
        p1 = edges[iter].p1;
        p2 = edges[iter].p2;
        edges[iter] = Edge(replace_map[p1], replace_map[p2]);
    }
}

Point rotate_point(Point in, double angle, Point pivot){

    double sine = sin(angle), cosine = cos(angle);
    double newx, newy;

    newx = cosine*in.x - sine*in.y + pivot.x*(1 - cosine) + pivot.y*sine;
    newy = cosine*in.y + sine*in.x + pivot.y*(1 - cosine) - pivot.x*sine;
    return Point(newx, newy);
}

void rotate_edges(vector<Edge> &edges, double angle, Point pivot){
    double sine = sin(angle), cosine = cos(angle);
    double newx, newy;
    set<Point> edge_points;
    unordered_map<Point, Point> replace_map;
    int iter;
    Point p1, p2;

    //Get the original points, then rotate them by the specified angle about
    //the specified pivot point
    edge_points = get_points(edges);
    for(Point p : edge_points){
        if(p == pivot) replace_map.insert(make_pair(p, p));
        else{
            newx = cosine*p.x - sine*p.y + pivot.x*(1 - cosine) + pivot.y*sine;
            newy = cosine*p.y + sine*p.x + pivot.y*(1 - cosine) - pivot.x*sine;
            replace_map.insert(make_pair(p, Point(newx, newy)));
        }
    }
    //delete &edge_points;

    //Look up the images of an edge's end points after rotation and change
    //the edge's end points to these new points
    for(iter = 0; iter < edges.size(); iter++){
        p1 = edges[iter].p1;
        p2 = edges[iter].p2;
        edges[iter].reset(replace_map[p1], replace_map[p2]);
    }
}

void get_extremes(vector<Edge> edges, double &minx, double &miny, double &maxx, double &maxy){
    minx = FLT_MAX;
    maxx = FLT_MIN;
    miny = FLT_MAX;
    maxy = FLT_MIN;
    Point p1, p2;

    for(Edge e : edges){
        p1 = e.p1;
        p2 = e.p2;

        minx = p1.x < minx ? p1.x : minx;
        minx = p2.x < minx ? p2.x : minx;
        maxx = p1.x > maxx ? p1.x : maxx;
        maxx = p2.x > maxx ? p2.x : maxx;

        miny = p1.y < miny ? p1.y : miny;
        miny = p2.y < miny ? p2.y : miny;
        maxy = p1.y > maxy ? p1.y : maxy;
        maxy = p2.y > maxy ? p2.y : maxy;
    }
}

double x_intersect(Edge e, double yval){
    double slope;
    Point p1 = e.p1, p2 = e.p2;

    slope = p1.x == p2.x ? BIG_SLOPE : (p2.y - p1.y) / (p2.x - p1.x);
    return p1.x + (yval - p1.y) / slope;
}

double y_intersect(Edge e, double xval){
    double slope;
    Point p1 = e.p1, p2 = e.p2;

    slope = p1.x == p2.x ? BIG_SLOPE : (p2.y - p1.y) / (p2.x - p1.x);
    return p1.y + slope * (xval - p1.x);
}

#endif

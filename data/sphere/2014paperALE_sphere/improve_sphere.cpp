/*! \file improve_sphere.cpp

    \brief <Program Name>

     Revision history:

          Revised by Christian Power dd.mm.yyyy
          Originally written by Christian Power
               (power22c@gmail.com) 15.06.2015

     Read a dgf file and improve the nodes. This codes only works for the unit sphere.
	 With a minor modification we can count vertices.

     Created by Christian Power on 15.06.2015
     Copyright (c) 2015 Christian Power. All rights reserved.
 */

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <algorithm>
#include <cmath>

using namespace std;

vector<string> load_file(const string ifname){
	ifstream ifs(ifname);
	if (!ifs) throw runtime_error("load_file()");
	vector<string> lines;
	for(string line; getline(ifs,line); )
		lines.push_back(line);
	ifs.close();
	return lines;
}

typedef vector<string>::iterator Line_pt;

void write_file(const string ofname, const vector<string>& content){
	ofstream ofs(ofname);
	for(const auto& line : content)
		ofs << line << endl;
	ofs.close();
}

void improve_sphere(Line_pt first, Line_pt last){
	while (first != last){
		istringstream line {*first};
		double d1,d2,d3;
		if( !(line >> d1 >> d2 >> d3) || !(line >> ws).eof() )
			throw runtime_error("improve_sphere()");

		// improving nodes
		const double node_norm = sqrt(d1*d1 + d2*d2 + d3*d3);
		d1 /= node_norm;
		d2 /= node_norm;
		d3 /= node_norm;

		ostringstream improved_line;
		improved_line << scientific
					  << d1 << ' ' << d2 << ' ' << d3;
		*first = improved_line.str();
		++first;
	}
}

int main() try{
	// improving sphere mesh code
	// const string filename = "Sphere_6.dgf";
	// vector<string> lines = load_file(filename);
	// 
	// Line_pt first_vertex = find(lines.begin(), lines.end(), "VERTEX");
	// if (first_vertex == lines.end()) throw runtime_error("main()");
	// ++first_vertex;	// prev. *first_vertex == VERTEX, now it is the first vertex 
	// Line_pt last_vertex = find(first_vertex, lines.end(), "#");
	// if (last_vertex == lines.end()) throw runtime_error("main()");
	//
	// improve_sphere(first_vertex, last_vertex);
	//
	// write_file(filename, lines);

	// counting vertices code
	ofstream notepad {"nr_of_vertices.txt"};

	const string basename = "Sphere_";
	const string suffix = ".dgf";

	for(int i = 0; i<6; ++i){
		const string dgfname = basename + to_string(i+1) + suffix;
		vector<string> lines = load_file(dgfname);

		Line_pt first_vertex = find(lines.begin(), lines.end(), "VERTEX");
		if (first_vertex == lines.end()) throw runtime_error("main()");
		++first_vertex;	// prev. *first_vertex == VERTEX, now it is the first vertex 
		Line_pt last_vertex = find(first_vertex, lines.end(), "#");
		if (last_vertex == lines.end()) throw runtime_error("main()");		
		
		notepad << dgfname << ": " << distance(first_vertex,last_vertex) << endl;
	}
	notepad.close();
}
catch(exception& e){
	cerr << e.what() << '\n';
	return 1;
}
catch(...){
	cerr << "An error occured!\n";
	return 2;
}

/*! Log:
 */

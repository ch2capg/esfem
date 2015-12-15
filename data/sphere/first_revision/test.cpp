#include <string>
#include <sstream>
#include <fstream>
#include <iostream> 
#include <exception>
using namespace std;

bool readLine_until_string(istream& is, string& s, const string& condition){
	getline(is,s);
	if(s == condition)
		return false;
	return true;
}

int main(int argc, char *argv[])
try{
	ifstream ifs {argv[1]};
	if(!ifs) throw runtime_error("File not found!");
	ofstream ofs {argv[2]};

	for(string s; readLine_until_string(ifs, s, "SIMPLEX");){
		ofs << s << endl;
	}
	ofs << "SIMPLEX" << endl;

	for(string s; readLine_until_string(ifs, s, "#");){
		double d1,d2,d3;
		stringstream ss {s};
		ss >> d1 >> d2 >> d3;
		int i1 = d1-1, i2 = d2-1, i3 = d3-1;
		ofs << i1 << ' ' << i2 << ' ' << i3 << endl;
	}
	ifs.close();
	ofs << "#\n#\n";
	ofs.close();
    return 0;
}
catch(exception& e){
	cerr << e.what() << '\n';
	return 1;
}

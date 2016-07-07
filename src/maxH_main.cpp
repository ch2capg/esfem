/*! \file maxH_main.cpp
    \brief Calculate the maximum h of a finite element grid

     Revision history
     --------------------------------------------------

          Revised by Christian Power July 2016
          Originally written by Christian Power
               (power22c@gmail.com) July 2016

    \author Christian Power 
    \date 07. July 2016
    \copyright Copyright (c) 2016 Christian Power.  All rights reserved.
 */

#include <iostream>
#include <valarray>
#include <vector>
#include <fstream>
#include <sstream>
#include <dassert.h>

using namespace std;

//! File as a vector of lines
struct file : vector<string>{
  //! Report an error
  struct bad : runtime_error{
    using runtime_error::runtime_error;
  };

  //! For error report
  const string filename;

  //! Read file
  explicit file(const string& fname)
    :filename {fname}

  {
    ifstream ifs {filename};
    if(!ifs) throw bad {Assert::compose(__FILE__, __LINE__, fname)};
    for(string s; getline(ifs, s); ) push_back(s);
    ifs.close();
  }
};

//! Nodes as a vector of valarrays of size 3
template<class T>
struct nodes : vector<valarray<T> >{
  //! Report an error
  struct bad : runtime_error{
    using runtime_error::runtime_error;
  };
  //! Read nodes from file
  /*! \warning Potentially narrowing double to int! */
  explicit nodes(const file& f){
    istringstream iss;
    iss.exceptions(ios_base::badbit | ios_base::failbit);
    size_t lineno {1};
    try{
      for(const auto& line : f){
	iss.str(line);
	double d1, d2, d3;
	iss >> d1 >> d2 >> d3;
	T t1(d1), t2(d2), t3(d3);  // this may narrow
	this->push_back({t1, t2, t3});
	iss.clear();
	++lineno;
      }
    }
    catch(...){
      ostringstream oss;
      oss << f.filename << ", line " << lineno;
      throw_with_nested(bad {Assert::compose(__FILE__, __LINE__, oss.str())});
    }
  }
};
//! Balazs node list
using nodelist = nodes<double>;
//! Balazs element list
using elemlist = nodes<int>;

//! Nodelist with corresponding elemlist
struct fempair{
  //! Report an error
  struct bad : runtime_error{
    using runtime_error::runtime_error;
  };

  //! The node list
  nodelist n;
  //! Zero based element list
  elemlist e;
  
  explicit fempair(const int& no) try
    :n {init_nodelist(no)}, e {init_elemlist(no)} {}
  catch(...){
    throw_with_nested(bad {Assert::compose(__FILE__, __LINE__ , "Constructor")});
  }
private:
  //! Constructor helper
  nodelist init_nodelist(const int& no){
    ostringstream oss;
    oss << SPHEREPATH << "/raw_data/Sphere_nodes" << no << ".txt";
    file f {oss.str()};
    return nodelist {f};
  }
  //! Constructor helper
  elemlist init_elemlist(const int& no){
    ostringstream oss;
    oss << SPHEREPATH << "/raw_data/Sphere_elements" << no << ".txt";
    file f {oss.str()};
    elemlist e {f};
    // make e zero based
    for(auto& vi : e) vi -= 1;
    return e;
  }
};

template<class F>
void evolve_nodes(fempair& fm, F f){
  for(auto& vd : fm.n) f(vd);
}

/*! \pre valarray has size of 3 */
vector<pair<int, int> > all_combinations(const valarray<int>& v){
  Assert::dynamic(v.size() == 3, Assert::compose(__FILE__, __LINE__, "!= 3"));
  return {{v[0], v[1]}, {v[0], v[2]}, {v[1], v[2]}}; // has size 3
}

/*! \pre valarray has size of 3*/
double distance(const valarray<double>& vd1, const valarray<double>& vd2){
  Assert::dynamic(vd1.size() == 3 && vd2.size() == 3,
		  Assert::compose(__FILE__, __LINE__, "!=3"));
  double rv {};
  for(int i = 0; i < 3; ++i) rv += pow(vd1[i] - vd2[i],2);
  return sqrt(rv);
}

double max_h(const fempair& fm){
  constexpr int no_edges {3};
  valarray<double> all_h(fm.e.size() * no_edges);
  auto ptr = begin(all_h);
  for(const auto& elem : fm.e){
    const auto comb = all_combinations(elem);
    for(const auto& p : comb){
      *ptr = distance(fm.n.at(p.first), fm.n.at(p.second));
      ++ptr;
    }
  }
  return all_h.max();
}

class evolution{
  double t {2.};
  double rE {2.};
  double rA {1.};
  double k {.5};
  double rt {(rE * rA)/( rE * exp(-k*t) + rA *( 1 - exp(-k*t)))};
public:
  template<class C>
  void operator()(C& c){ c *= rt; }
};

//! For testing nodes<T>
template<class T>
void print(const valarray<T>& v){
  cout << "(";
  for(size_t i = 0; i < v.size() - 1; ++i) cout << v[i] << ", ";
  cout << v[v.size()-1] << ")" << endl;
}

void print_what(const exception& e){
  cerr << e.what() << '\n';
  try{
    rethrow_if_nested(e);
  }
  catch(const exception& e){
    print_what(e);
  }
  catch(...){
    cerr << "Unknown nested error\n";
  }
}

int main() try{
  for(int i = 1; i <= 5; ++i){
    fempair fp {i};
    evolve_nodes(fp, evolution {});
    const auto h = max_h(fp);
    cout << "h_max " << i << ": " << h << endl;
  }
 }
 catch(const exception& e){
   print_what(e);
   return 1;
 }
 catch(...){
   cerr << "Unknown error\n";
   return 2;
 }

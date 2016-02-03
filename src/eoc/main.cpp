#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <valarray>

class Cmd_args;
class Data_reader;

class Cmd_args{
public:
  const std::vector<std::string>& filelist() const;
private:
  std::vector<std::string> all_errFiles {};
};
class Data_reader{
public:
  Data_reader(const std::string& filename);
  
private:
  std::valarray<double> l2_vector {};
  std::valarray<double> h1_vector {};
  double dT {0.};

};

void print_errMsg(const std::exception&);
Cmd_args check_usage(const int argc, const char* const* argv);

int main(int argc, char** argv) try{
  const auto cmd_args = check_usage(argc, argv);
  const auto filelist = cmd_args.filelist();
  for(auto& file : filelist){
    
  }
}
catch(const std::exception& e){
  print_errMsg(e);
  return 1;
}
catch(...){
  std::cerr << "Unkown error in main() detected\n";
  return 2;
}

// ----------------------------------------------------------------------
// Internal Implementation

void print_errMsg(const std::exception& e){
  std::cerr << e.what() << '\n';
  try{
    std::rethrow_if_nested(e);
  }
  catch(const std::exception& nested){
    print_errMsg(nested);
  }
}

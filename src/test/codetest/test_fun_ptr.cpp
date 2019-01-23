// short example to test function callback objects
// compile with:
// $ g++ -Wall --std=c++11 test_fun_ptr.cpp -o fun

#include <string>
#include <functional>
#include <vector>
#include <iostream>
#include <sstream>

using namespace std::placeholders; // for _1, _2 etc.

// define our container structure
struct data_def {
  std::string name;
  float * array_ptr;
  int array_size;
  std::function<std::string(int)> fun;
};


// define a class that knows how to print stuff
class knowsaboutdata
{
public:
  knowsaboutdata()
  {

  };
  // define some functions

  // how to handle 1D index
  std::string index_1d(int idx)
  {
    char buff[100];
    snprintf(buff, sizeof(buff), "%d", idx);
    std::string buff_string = buff;
    
    return buff_string;
  }

  // how to handle 2d index
  std::string index_2d(int c)
  {
    std::ostringstream string_stream;
    
    string_stream << c/array2d_size_x << "\t" << c%array2d_size_y;
  
    std::string copy_of_str = string_stream.str();
    
    return copy_of_str;
  }

  // some data arrays
  float array1d[13] = {1,2,3,4,5,6,7,8,9,10,11,12,13};
  const int array1d_size = 13;
  float array2d[9] = {1,2,3,4,5,6,7,8,9};
  const int array2d_size_x = 3;
  const int array2d_size_y = 3;
};


void print_def(data_def& d, int i)
{
  std::string coord = "";
  if (d.fun != nullptr)
    coord = d.fun(i);
  
  std::cout << d.name << "\t" <<i << "\t"<<coord<< std::endl;
  
}


int main(int argc, char * argv[])
{
  knowsaboutdata kad;
  // define array of definitions
  std::vector<data_def> data_definitions  = {
					     { std::string("array1d"), kad.array1d, kad.array1d_size,std::bind(&knowsaboutdata::index_1d, &kad, _1)},
					     { std::string("array2d"), kad.array2d, kad.array2d_size_x* kad.array2d_size_y, std::bind(&knowsaboutdata::index_2d, &kad, _1)},
					     // the empty function case
					     { "array1d_nofun", kad.array1d, kad.array1d_size, std::function<std::string(int)>()}
  };

  // loop on definitions and print out index 4
  for (auto & def : data_definitions)
  {
    print_def(def,7);
  }
  
  return 0;
}




 

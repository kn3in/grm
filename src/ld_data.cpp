#include <fstream>
#include <stdexcept>
#include <vector>
#include "ld_data.hpp"

ld_data::ld_data(std::string path_to_ld_file) {
   std::ifstream ld_file(path_to_ld_file.c_str());
   if(!ld_file) throw std::runtime_error("Runtime error: can not open ld file (.bed as per ucsc bed format) [" + path_to_ld_file + "]");

      int my_int;
      
      while(ld_file) {
      
      ld_file >> my_int; // attempt to read
      if(ld_file.eof()) break; // are we past EOF?
      chr.push_back(my_int);

      ld_file >> my_int;
      start.push_back(my_int);

      ld_file >> my_int;
      stop.push_back(my_int);

   }
   
   ld_file.close();

}

ld_data::~ld_data() {

}

#ifndef INPUT_H
#define INPUT_H

#include "Point.h"
#include <string>
#include <vector>
class ReadInput

{
public:
	ReadInput(const std::string & filename, std::vector<Point> & coord, std::vector<std::vector<int>> & inc);

	

private:
	std::vector<std::string> split(std::string line, std::string delimiter = ",");

};



#endif // !INPUT_H

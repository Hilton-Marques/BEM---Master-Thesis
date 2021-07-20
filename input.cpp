#include "input.h"
#include <fstream>

ReadInput::ReadInput(const std::string& filename, std::vector<Point>& coord, std::vector<std::vector<int>> & inc)
{
  std::ifstream  input(filename);
  std::string firstLine;
  std::getline(input, firstLine);
  std::vector<std::string> fistLineSplited = split(firstLine);
  int nPts = std::stoi(fistLineSplited[0]);
  int nEls = std::stoi(fistLineSplited[1]);

  for (int i = 0; i < nPts; i++)
  {
    std::string line;
    std::getline(input, line);
    std::vector<std::string> lineSplited = split(line);
    coord.push_back( Point(std::stod(lineSplited[0]), std::stod(lineSplited[1]), std::stod(lineSplited[2])));
  }
  for (int i = 0; i < nEls; i++)
  {
    std::string line;
    std::getline(input, line);
    std::vector<std::string> lineSplited = split(line);
    std::vector<int> incI = { std::stoi(lineSplited[0]) , std::stoi(lineSplited[1]), std::stoi(lineSplited[2]) };
    inc.push_back(incI);
  }
}
std::vector<std::string> ReadInput::split(std::string line, std::string delimeter)
{
  std::vector<std::string> res;

  int start = 0;
  int end = line.find(",");
  while (end != -1) {
    res.push_back(line.substr(start, end - start));
    start = end + delimeter.size();
    end = line.find(delimeter, start);
  }
  res.push_back(line.substr(start,line.size()));

  return res;
}

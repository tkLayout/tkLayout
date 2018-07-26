#include "ITCabling/HvLine.hh"


HvLine::HvLine(const std::string name) :
  name_(name)
{};


HvLine::~HvLine() {
  delete powerChain_;    // TO DO: switch to smart pointers and remove this!
  powerChain_ = nullptr;
}

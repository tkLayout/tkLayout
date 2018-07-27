#include "OuterCabling/DTC.hh"


DTC::DTC(const std::string name, const double phiSectorWidth, const int phiSectorRef, const Category& type, const int slot, const bool isPositiveCablingSide) :
  name_(name),
  phiSectorWidth_(phiSectorWidth),
  phiSectorRef_(phiSectorRef),
  type_(type),
  slot_(slot),
  isPositiveCablingSide_(isPositiveCablingSide)
{
  plotColor_ = computePlotColor(phiSectorRef, type, slot);
};


const int DTC::computePlotColor(const int phiSectorRef, const Category& type, const int slot) const {
  int plotColor = 0;
  if (type == Category::PS10G || type == Category::PS5G) plotColor = slot;
  else if (type == Category::SS) plotColor = 6 + slot;
  plotColor += 12 * phiSectorRef;
  return plotColor;
}

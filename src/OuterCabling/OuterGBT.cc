#include "OuterCabling/OuterGBT.hh"


OuterGBT::OuterGBT(const int GBTId) :
  myGBTCMSSWId_(0) // Want to have consecutive integers, ordered by detID, so done after the full map is created
{
  myGBTId_ = GBTId;
};


/*
 *  Connect a Module to the GBT.
 */
void OuterGBT::addModule(Module* m) { 
  modules_.push_back(m);
}

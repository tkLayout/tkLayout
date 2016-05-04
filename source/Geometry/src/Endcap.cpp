#include "Endcap.h"
#include "MessageLogger.h"
#include "SupportStructure.h"

using material::SupportStructure;

void Endcap::cutAtEta(double eta) { 
  for (auto& d : disks_) d.cutAtEta(eta); 
  disks_.erase_if([](const Disk& d) { return d.numRings() == 0; });
  numDisks(disks_.size());
}

vector<double> Endcap::findMaxDsDistances() { // drill down into the property trees to find the maximum dsDistance
  vector<double> maxDsDistances; 
  double endcapDsDistance = propertyTree().get("dsDistance", 0.);

  PropertyNode<int> ringNode("");
  for (auto& tel : pair2range(propertyTree().equal_range("Ring"))) // scan Ring subtrees outside Disks
    ringNode.fromPtree(tel.second);
  for (auto& rnel : ringNode) {
    Property<double, NoDefault> ringDsDistance; //endcapDsDistance);
    for (auto& tel : pair2range(rnel.second.equal_range("dsDistance"))) ringDsDistance.fromPtree(tel.second);
    if (ringDsDistance.state()) {
      if (maxDsDistances.size() < rnel.first) maxDsDistances.resize(rnel.first);//, endcapDsDistance);
      maxDsDistances[rnel.first-1] = MAX(maxDsDistances[rnel.first-1], ringDsDistance()); 
    }
  }

  for (auto& dnel : diskNode) {
    double diskDsDistance = dnel.second.get("dsDistance", endcapDsDistance);
    ringNode.clear();
    for (auto& tel : pair2range(dnel.second.equal_range("Ring"))) ringNode.fromPtree(tel.second); // scan Ring subtrees inside Disks
    for (auto& rnel : ringNode) {
      Property<double, NoDefault> ringDsDistance; // (diskDsDistance);
      for (auto& tel : pair2range(rnel.second.equal_range("dsDistance"))) ringDsDistance.fromPtree(tel.second);
      if (maxDsDistances.size() < rnel.first) maxDsDistances.resize(rnel.first); //, diskDsDistance);
      if (ringDsDistance.state()) {
        maxDsDistances[rnel.first-1] = MAX(maxDsDistances[rnel.first-1], ringDsDistance()); 
      } else {
        maxDsDistances[rnel.first-1] = MAX(maxDsDistances[rnel.first-1], diskDsDistance);
      }
    }
  }
  maxDsDistances.push_back(endcapDsDistance); // adds a default element to be used in case disks need more rings than the vector specifies dsDistances for
  return maxDsDistances;
}

void Endcap::build() {
  try {
    logINFO(Form("Building %s", fullid(*this).c_str()));
    check();

    if (!innerZ.state()) innerZ(barrelMaxZ() + barrelGap());
    else if(barrelGap.state()) logWARNING("'innerZ' was set, ignoring 'barrelGap'");

    vector<double> maxDsDistances = findMaxDsDistances();
    vector<Disk*> tdisks;

    double alpha = pow(outerZ()/innerZ(), 1/double(numDisks()-1)); // geometric progression factor

    for (int i = 1; i <= numDisks(); i++) {
      Disk* diskp = GeometryFactory::make<Disk>();
      diskp->myid(i);

      // Standard is to build & calculate parameters for the central disc (in the middle)
      diskp->buildZ((innerZ() + outerZ())/2);

      // Apply correct offset for each disc versus middle position
      double offset = pow(alpha, i-1) * innerZ();
      diskp->placeZ(offset);

      // Store parameters in a tree
      diskp->store(propertyTree());
      if (diskNode.count(i) > 0) diskp->store(diskNode.at(i));

      // To test the extreme cases -> one needs to test either first or last layer (based on parity)
      diskp->zHalfLength((outerZ()-innerZ())/2.);

      // Build
      diskp->build(maxDsDistances);

      // Mirror discs
      Disk* diskn = GeometryFactory::clone(*diskp);
      diskn->mirrorZ();

      tdisks.push_back(diskp);
      tdisks.push_back(diskn);
    }
    std::stable_sort(tdisks.begin(), tdisks.end(), [](Disk* d1, Disk* d2) { return d1->minZ() < d2->maxZ(); });
    for (Disk* d : tdisks) disks_.push_back(d);
    
  } catch (PathfulException& pe) { pe.pushPath(fullid(*this)); throw; }

  // Supports defined within a Barrel
  for (auto& mapel : supportNode) {
    SupportStructure* s = new SupportStructure();
    s->store(propertyTree());
    s->store(mapel.second);
    s->buildInEndcap(*this);
    supportStructures_.push_back(s);
  }

  cleanup();
  builtok(true);
}

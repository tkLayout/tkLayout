#include "Endcap.h"

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
    std::cout << ">>> Building " << fullid(*this) << " <<<" << std::endl;
    check();

    if (!minZ.state()) minZ(barrelMaxZ() + barrelGap());

    vector<double> maxDsDistances = findMaxDsDistances();
    vector<Disk*> tdisks;

    double alpha = pow(maxZ()/minZ(), 1/double(numDisks()-1)); // geometric progression factor

    for (int i = 1; i <= numDisks(); i++) {
      Disk* diskp = new Disk();
      diskp->setup();
      diskp->myid(i);

      diskp->buildZ((minZ() + maxZ())/2);
      double offset = pow(alpha, i-1) * minZ();
      diskp->placeZ(offset);

      diskp->store(propertyTree());
      if (diskNode.count(i) > 0) diskp->store(diskNode.at(i));

      diskp->zError(diskp->zError() + (maxZ()-minZ())/2);

      diskp->build(maxDsDistances);

      Disk* diskn = new Disk(*diskp);
      diskn->setup();
      diskn->mirrorZ();

      tdisks.push_back(diskp);
      tdisks.push_back(diskn);
    }
    std::stable_sort(tdisks.begin(), tdisks.end(), [](Disk* d1, Disk* d2) { return d1->minZ() < d2->maxZ(); });
    for (Disk* d : tdisks) disks_.push_back(d);
    
  } catch (PathfulException& pe) { pe.pushPath(fullid(*this)); throw; }

  cleanup();
  builtok(true);
}

define_enum_strings(Endcap::MinZType) = { "absolute", "barrelgap" };

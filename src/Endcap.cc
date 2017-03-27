#include "Endcap.hh"
#include "MessageLogger.hh"
#include "SupportStructure.hh"

using material::SupportStructure;

void Endcap::cutAtEta(double eta) { 
  for (auto& d : disks_) d.cutAtEta(eta); 
  disks_.erase_if([](const Disk& d) { return d.numRings() == 0; });
  numDisks(disks_.size());
}


/** Scans the Disk property tree.
    This info is needed to build any disk in the Endcap.
*/
const ScanDiskInfo Endcap::scanDiskPropertyTree(int diskNumber) const {

  Disk* diskTemplate = GeometryFactory::make<Disk>();
  diskTemplate->myid(diskNumber);
  diskTemplate->store(propertyTree());
  if (diskNode.count(diskNumber) > 0) diskTemplate->store(diskNode.at(diskNumber));
  const ScanDiskInfo& diskInfo = diskTemplate->scanPropertyTree();
  return diskInfo;
}

/** Scans the Endcap property tree and gathers info from the innermost and outermost disks in the Endcap.
    This info is needed to build any disk in the Endcap.
*/
const ScanEndcapInfo Endcap::scanPropertyTree() const {

  const ScanDiskInfo& innermostDiskInfo = scanDiskPropertyTree(1);
  const ScanDiskInfo& outermostDiskInfo = scanDiskPropertyTree(numDisks());
  const ScanEndcapInfo& extremaDisksInfo = std::make_pair(innermostDiskInfo, outermostDiskInfo);
  return extremaDisksInfo;
}


void Endcap::build() {
  try {
    logINFO(Form("Building %s", fullid(*this).c_str()));
    check();

    if (!innerZ.state()) innerZ(barrelMaxZ() + barrelGap());
    else if(barrelGap.state()) logWARNING("'innerZ' was set, ignoring 'barrelGap'");

    // Before any disk is build, needs to take info from PropertyTree from innermost and outermost disks in the Endcap.
    // (In +Z side, 'innermost' corresponds to 'lower Z', and 'outermost' corresponds to 'bigger Z').
    ScanEndcapInfo extremaDisksInfo = scanPropertyTree();

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
      diskp->build(extremaDisksInfo);

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



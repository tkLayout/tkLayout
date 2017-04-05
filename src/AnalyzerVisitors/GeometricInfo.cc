#include <AnalyzerVisitors/GeometricInfo.hh>


   //**************************************//
    //*             Visitor                *//
    //*  Layers and disks : global info    *//
    //*                                    *//
    //**************************************//

void LayerDiskSummaryVisitor::preVisit() {
  layerTable->setContent(0, 0, "Barrel :");
  layerTable->setContent(1, 0, "Layer");
  layerTable->setContent(2, 0, "r");
  layerTable->setContent(3, 0, "z_max");
  layerTable->setContent(4, 0, "# rods");
  layerTable->setContent(5, 0, "# mods");
  diskTable->setContent(0, 0, "Endcap :");
  diskTable->setContent(1, 0, "Disk");
  diskTable->setContent(2, 0, "z");
  diskTable->setContent(3, 0, "# rings");
  diskTable->setContent(4, 0, "# mods");
  nMB = SimParms::getInstance().numMinBiasEvents();
}

void LayerDiskSummaryVisitor::visit(const Barrel& b) {
  layerTable->setContent(0, 1 + nBarrelLayers, b.myid());
}

void LayerDiskSummaryVisitor::visit(const Layer& l) {
  if (l.maxZ() < 0.) return;
  ++nBarrelLayers;
  totalBarrelModules += l.totalModules();
  layerTable->setContent(1, nBarrelLayers, l.myid());
  layerTable->setContent(2, nBarrelLayers, l.placeRadius(), coordPrecision);
  layerTable->setContent(3, nBarrelLayers, l.maxZ(), coordPrecision);
  layerTable->setContent(4, nBarrelLayers, l.numRods());
  layerTable->setContent(5, nBarrelLayers, l.totalModules());
}

void LayerDiskSummaryVisitor::visit(const Endcap& e) {
  nEndcaps++;
  diskTable->setContent(0, 1 + nDisks, e.myid());

  RootWTable* diskName = new RootWTable();
  diskName->setContent(0, 0, e.myid() + ",  Disc 1 :");
  diskNames.push_back(diskName);

  RootWTable* ringTable = new RootWTable();
  ringTable->setContent(0, 0, "Ring :");
  ringTable->setContent(1, 0, "r"+subStart+"min"+subEnd);
  ringTable->setContent(2, 0, "r"+subStart+"low"+subEnd);
  ringTable->setContent(3, 0, "r"+subStart+"centre"+subEnd);
  ringTable->setContent(4, 0, "r"+subStart+"high"+subEnd);
  ringTable->setContent(5, 0, "r"+subStart+"max"+subEnd);
  ringTable->setContent(6, 0, "# mods");
  ringTables.push_back(ringTable);
}

void LayerDiskSummaryVisitor::visit(const Disk& d) {
  nRings = 0;
  if (d.averageZ() < 0.) return;
  ++nDisks;
  totalEndcapModules += d.totalModules();
  diskTable->setContent(1, nDisks, d.myid());
  diskTable->setContent(2, nDisks, d.averageZ(), coordPrecision);
  diskTable->setContent(4, nDisks, d.totalModules());
}

void LayerDiskSummaryVisitor::visit(const Ring& r) {
  if (r.averageZ() < 0. || r.numModules() == 0) return;
  ++nRings;
  diskTable->setContent(3, nDisks, nRings);
  ringTables.at(nEndcaps-1)->setContent(0, nRings, r.myid());
  ringTables.at(nEndcaps-1)->setContent(6, nRings, r.numModules());
}

void LayerDiskSummaryVisitor::visit(const Module& m) {
  TagMaker tmak(m);

  std::string aSensorTag = tmak.sensorGeoTag;
  tagMapPositions[aSensorTag].insert(tmak.posTag);
  tagMapCount[aSensorTag]++;
  tagMapCountChan[aSensorTag] += m.totalChannels();
  tagMapMaxStripOccupancy[aSensorTag] = MAX(m.stripOccupancyPerEvent()*nMB, tagMapMaxStripOccupancy[aSensorTag]);
  tagMapMaxHitOccupancy[aSensorTag] = MAX(m.hitOccupancyPerEvent()*nMB, tagMapMaxHitOccupancy[aSensorTag]);
  tagMapAveStripOccupancy[aSensorTag] += m.stripOccupancyPerEvent()*nMB;
  tagMapAveHitOccupancy[aSensorTag] += m.hitOccupancyPerEvent()*nMB;
  // modules' spatial resolution along the local X axis is not parametrized
  if (!m.hasAnyResolutionLocalXParam()) {
    tagMapIsParametrizedXResolution[aSensorTag] = false;
    tagMapAveRphiResolution[aSensorTag] += m.nominalResolutionLocalX();
    tagMapAveRphiResolutionRmse[aSensorTag] += 0.;
  }
  // modules' spatial resolution along the local X axis is parametrized
  else {
    tagMapIsParametrizedXResolution[aSensorTag] = true;
    if (boost::accumulators::count(m.rollingParametrizedResolutionLocalX) > 0) { // calculation only on hit modules
      tagMapSumXResolution[aSensorTag] += sum(m.rollingParametrizedResolutionLocalX);
      tagMapSumSquaresXResolution[aSensorTag] += moment<2>(m.rollingParametrizedResolutionLocalX) * boost::accumulators::count(m.rollingParametrizedResolutionLocalX);
      tagMapCountXResolution[aSensorTag] += double(boost::accumulators::count(m.rollingParametrizedResolutionLocalX));
    }
  }
  // modules' spatial resolution along the local Y axis is not parametrized
  if (!m.hasAnyResolutionLocalYParam()) {
    tagMapIsParametrizedYResolution[aSensorTag] = false;
    tagMapAveYResolution[aSensorTag] += m.nominalResolutionLocalY();
    tagMapAveYResolutionRmse[aSensorTag] += 0.;
  }
  // modules' spatial resolution along the local Y axis is parametrized
  else {
    tagMapIsParametrizedYResolution[aSensorTag] = true;
    if (boost::accumulators::count(m.rollingParametrizedResolutionLocalY) > 0) { // calculation only on hit modules
      tagMapSumYResolution[aSensorTag] += sum(m.rollingParametrizedResolutionLocalY);
      tagMapSumSquaresYResolution[aSensorTag] += moment<2>(m.rollingParametrizedResolutionLocalY) * boost::accumulators::count(m.rollingParametrizedResolutionLocalY);
      tagMapCountYResolution[aSensorTag] += double(boost::accumulators::count(m.rollingParametrizedResolutionLocalY));
    }
  }
  //tagMapAveRphiResolutionTrigger[aSensorTag] += m.resolutionRPhiTrigger();
  //tagMapAveYResolutionTrigger[aSensorTag] += m.resolutionYTrigger();
  tagMapSensorPowerAvg[aSensorTag] += m.sensorsIrradiationPowerMean();
  if (tagMapSensorPowerMax[aSensorTag] < m.sensorsIrradiationPowerMean()) tagMapSensorPowerMax[aSensorTag] = m.sensorsIrradiationPowerMean();
  totCountMod++;
  totCountSens += m.numSensors();
  totChannel += m.totalChannels();
  totArea += m.area()*m.numSensors();
  totalSensorPower += m.sensorsIrradiationPowerMean();
  if (tagMap.find(aSensorTag)==tagMap.end()){
    // We have a new sensor geometry
    tagMap[aSensorTag] = &m;
  }
}

void LayerDiskSummaryVisitor::visit(const EndcapModule& m) {
  if (m.side() != 1 || m.disk() != 1) return;

  ringTables.at(nEndcaps-1)->setContent(1, nRings, m.minR(), coordPrecision);
  ringTables.at(nEndcaps-1)->setContent(2, nRings, sqrt(pow(m.minR(),2)+pow(m.minWidth()/2.,2)), coordPrecision); // Ugly, this should be accessible as a method
  ringTables.at(nEndcaps-1)->setContent(3, nRings, m.center().Rho(), coordPrecision);
  ringTables.at(nEndcaps-1)->setContent(4, nRings, m.minR()+m.length(), coordPrecision);
  ringTables.at(nEndcaps-1)->setContent(5, nRings, m.maxR(), coordPrecision);
}

void LayerDiskSummaryVisitor::postVisit() {
  layerTable->setContent(0, nBarrelLayers+1, "Total");
  layerTable->setContent(5, nBarrelLayers+1, totalBarrelModules);
  diskTable->setContent(0, nDisks+1, "Total");
  diskTable->setContent(4, nDisks+1, totalEndcapModules*2);
}



    //***************************************//
    //*                Visitor              *//
    //* Automatic-placement tilted layers : *//
    //*            Additional info          *//
    //*                                     *//
    //***************************************//

/**
   Visits tilted layers and gathers info on its flat and tilted parts.
*/
void TiltedLayersVisitor::visit(const Layer& l) {
  // Only for tilted layers with automatic placement.
  if (l.isTilted() && l.isTiltedAuto()) {
    numTiltedLayers++;

    // Initializes layer name
    RootWTable* tiltedLayerName = new RootWTable();
    tiltedLayerName->setContent(0, 0, "Layer " + std::to_string(l.myid()) + " :");
    tiltedLayerNames.push_back(tiltedLayerName);

    // Initializes flat part name
    RootWTable* flatPartName = new RootWTable();
    flatPartName->setContent(0, 0, "Flat part :");
    flatPartNames.push_back(flatPartName);

    // Initializes tilted part name
    RootWTable* tiltedPartName = new RootWTable();
    tiltedPartName->setContent(0, 0, "Tilted part :");
    tiltedPartNames.push_back(tiltedPartName);

    // FILLS TILTED PART TABLE
    RootWTable* tiltedPartTable = new RootWTable();
    for (int i=0; i < l.tiltedRingsGeometry().size(); i++) {
      int ringNumber = l.buildNumModulesFlat() + 1 + i;
      tiltedPartTable->setContent(0, 0, "Ring");
      tiltedPartTable->setContent(0, i+1, ringNumber);
      tiltedPartTable->setContent(1, 0, "tiltAngle (°)");
      tiltedPartTable->setContent(1, i+1, l.tiltedRingsGeometry()[ringNumber]->tiltAngle(), anglePrecision);
      tiltedPartTable->setContent(2, 0, "tiltAngleIdeal" + subStart + "Inner" + subEnd + " (°)");
      tiltedPartTable->setContent(2, i+1, l.tiltedRingsGeometry()[ringNumber]->tiltAngleIdealInner(), anglePrecision);
      tiltedPartTable->setContent(3, 0, "deltaTiltIdeal" + subStart + "Inner" + subEnd + " (°)");
      tiltedPartTable->setContent(3, i+1, l.tiltedRingsGeometry()[ringNumber]->deltaTiltIdealInner(), anglePrecision);
      tiltedPartTable->setContent(4, 0, "tiltAngleIdeal" + subStart + "Outer" + subEnd + " (°)");
      tiltedPartTable->setContent(4, i+1, l.tiltedRingsGeometry()[ringNumber]->tiltAngleIdealOuter(), anglePrecision);
      tiltedPartTable->setContent(5, 0, "deltaTiltIdeal" + subStart + "Outer" + subEnd + " (°)");
      tiltedPartTable->setContent(5, i+1, l.tiltedRingsGeometry()[ringNumber]->deltaTiltIdealOuter(), anglePrecision);
      tiltedPartTable->setContent(6, 0, "theta_g (°)");
      tiltedPartTable->setContent(6, i+1, l.tiltedRingsGeometry()[ringNumber]->theta_g(), anglePrecision);
      tiltedPartTable->setContent(7, 0, "r" + subStart + "Inner" + subEnd);
      tiltedPartTable->setContent(7, i+1, l.tiltedRingsGeometry()[ringNumber]->innerRadius(), coordPrecision);
      tiltedPartTable->setContent(8, 0, "r" + subStart + "Outer" + subEnd);
      tiltedPartTable->setContent(8, i+1, l.tiltedRingsGeometry()[ringNumber]->outerRadius(), coordPrecision);
      tiltedPartTable->setContent(9, 0, "averageR (on Ring)");
      tiltedPartTable->setContent(9, i+1, l.tiltedRingsGeometry()[ringNumber]->averageR(), coordPrecision);
      tiltedPartTable->setContent(10, 0, "gapR");
      tiltedPartTable->setContent(10, i+1, l.tiltedRingsGeometry()[ringNumber]->gapR(), coordPrecision);
      tiltedPartTable->setContent(11, 0, "z" + subStart + "Inner" + subEnd);
      tiltedPartTable->setContent(11, i+1, l.tiltedRingsGeometry()[ringNumber]->zInner(), coordPrecision);
      tiltedPartTable->setContent(12, 0, "z" + subStart + "Outer" + subEnd);
      tiltedPartTable->setContent(12, i+1, l.tiltedRingsGeometry()[ringNumber]->zOuter(), coordPrecision);
      tiltedPartTable->setContent(13, 0, "averageZ (on Ring)");
      tiltedPartTable->setContent(13, i+1, l.tiltedRingsGeometry()[ringNumber]->averageZ(), coordPrecision);
      tiltedPartTable->setContent(14, 0, "deltaZ" + subStart + "Inner" + subEnd + " (Ring i & i-1)");
      tiltedPartTable->setContent(14, i+1, l.tiltedRingsGeometryInfo().deltaZInner()[ringNumber], coordPrecision);
      tiltedPartTable->setContent(15, 0, "deltaZ" + subStart + "Outer" + subEnd + " (Ring i & i-1)");
      tiltedPartTable->setContent(15, i+1, l.tiltedRingsGeometryInfo().deltaZOuter()[ringNumber], coordPrecision);
      tiltedPartTable->setContent(16, 0, "phiOverlap");
      tiltedPartTable->setContent(16, i+1, l.tiltedRingsGeometry()[ringNumber]->phiOverlap(), coordPrecision);
      tiltedPartTable->setContent(17, 0, "zOverlap" + subStart + "Outer" + subEnd);
      tiltedPartTable->setContent(17, i+1, l.tiltedRingsGeometry()[ringNumber]->ringZOverlap(), zOverlapPrecision);
      double zErrorInner = l.tiltedRingsGeometryInfo().zErrorInner()[ringNumber];
      tiltedPartTable->setContent(18, 0, "zError" + subStart + "Inner" + subEnd + " (Ring i & i-1)");
      if (!std::isnan(zErrorInner)) tiltedPartTable->setContent(18, i+1, zErrorInner, coordPrecision);
      else tiltedPartTable->setContent(18, i+1, "n/a");
      double zErrorOuter = l.tiltedRingsGeometryInfo().zErrorOuter()[ringNumber];
      tiltedPartTable->setContent(19, 0, "zError" + subStart + "Outer" + subEnd + " (Ring i & i-1)");
      if (!std::isnan(zErrorOuter)) tiltedPartTable->setContent(19, i+1, zErrorOuter, coordPrecision);
      else tiltedPartTable->setContent(19, i+1, "n/a");
    }
    tiltedPartTables.push_back(tiltedPartTable);

    // FILLS FLAT PART TABLE
    RootWTable* flatPartTable = new RootWTable();

    StraightRodPair* minusBigDeltaRod = (l.bigParity() > 0 ? l.flatPartRods().at(1) : l.flatPartRods().front());
    const auto& minusBigDeltaModules = minusBigDeltaRod->modules().first;
    StraightRodPair* plusBigDeltaRod = (l.bigParity() > 0 ? l.flatPartRods().front() : l.flatPartRods().at(1));
    const auto& plusBigDeltaModules = plusBigDeltaRod->modules().first;

    int i = 0;
    for (const auto& m : minusBigDeltaModules) {
      int ringNumber = i + 1;
      flatPartTable->setContent(0, 0, "Ring");
      flatPartTable->setContent(0, i+1, ringNumber);
      flatPartTable->setContent(1, 0, "r" + subStart + "Inner" + subEnd);
      flatPartTable->setContent(1, i+1, m.center().Rho(), coordPrecision);
      flatPartTable->setContent(3, 0, "averageR (on Flat part)");
      flatPartTable->setContent(3, i+1, l.flatPartAverageR(), coordPrecision);
      flatPartTable->setContent(4, 0, "bigDelta");
      flatPartTable->setContent(4, i+1, l.bigDelta(), coordPrecision);
      flatPartTable->setContent(5, 0, "smallDelta");
      flatPartTable->setContent(5, i+1, l.smallDelta(), coordPrecision);
      flatPartTable->setContent(6, 0, "z");
      flatPartTable->setContent(6, i+1, m.center().Z(), coordPrecision);
      flatPartTable->setContent(7, 0, "phiOverlap");
      flatPartTable->setContent(7, i+1, (((minusBigDeltaRod->zPlusParity() * pow(-1, (i%2))) > 0) ? l.flatPartPhiOverlapSmallDeltaPlus() : l.flatPartPhiOverlapSmallDeltaMinus()), coordPrecision);
      // In case beamSpotCover == false, zOverlap is the only parameter used as a Z-coverage constraint in the geometry construction process.
      // (There is then no zError taken into account).
      // As a result, it is interesting to display zOverlap !
      int extraLine = 0;
      if (!l.flatPartRods().front()->beamSpotCover()) {
	extraLine = 1;
	flatPartTable->setContent(8, 0, "zOverlap");
	flatPartTable->setContent(8, i+1, l.flatPartRods().front()->zOverlap(), coordPrecision);
	// WARNING : zOverlap, in the geometry construction process, is one value common for all flat part (or straight rod)
      }
      // calculates zError for ring (i) with respect to ring (i-1)
      if (i > 0) {
	// Inner modules in the tilted rings
	double zErrorInner = l.flatRingsGeometryInfo().zErrorInner()[i];
	flatPartTable->setContent(8 + extraLine, 0, "zError" + subStart + "Inner" + subEnd + " (Ring i & i-1)");
	if (!std::isnan(zErrorInner)) flatPartTable->setContent(8 + extraLine, i+1, zErrorInner, coordPrecision);
	else flatPartTable->setContent(8 + extraLine, i+1, "n/a");
	// Outer modules in the tilted rings
	double zErrorOuter = l.flatRingsGeometryInfo().zErrorOuter()[i];
	flatPartTable->setContent(9 + extraLine, 0, "zError" + subStart + "Outer" + subEnd + " (Ring i & i-1)");
	if (!std::isnan(zErrorOuter)) flatPartTable->setContent(9 + extraLine, i+1, zErrorOuter, coordPrecision);
	else flatPartTable->setContent(9 + extraLine, i+1, "n/a");
      }
      i++;
    }
    i = 0;
    for (const auto& m : plusBigDeltaModules) {
      flatPartTable->setContent(2, 0, "r" + subStart + "Outer" + subEnd);
      flatPartTable->setContent(2, i+1, m.center().Rho(), coordPrecision);
      i++;
    }
    flatPartTables.push_back(flatPartTable);
  } // end of 'fills flat part table'
}


    //************************************//
    //*               Visitor             //
    //*            AllModulesCsv          //
    //*                                   //
    //************************************//
void TrackerVisitor::preVisit() {
  //output_ << "Section/C:Layer/I:Ring/I:r_mm/D:z_mm/D:tiltAngle_deg/D:phi_deg/D:meanWidth_mm/D:length_mm/D:sensorSpacing_mm/D:sensorThickness_mm/D, DetId/I" << std::endl;
  output_ << "DetId/U, BinaryDetId/B, Section/C, Layer/I, Ring/I, r_mm/D, z_mm/D, tiltAngle_deg/D, phi_deg/D, meanWidth_mm/D, length_mm/D, sensorSpacing_mm/D, sensorThickness_mm/D" << std::endl;
}

void TrackerVisitor::visit(const Barrel& b) {
  sectionName_ = b.myid();
}

void TrackerVisitor::visit(const Endcap& e) {
  sectionName_ = e.myid();
}

void TrackerVisitor::visit(const Layer& l) {
  layerId_ = l.myid();
}

void TrackerVisitor::visit(const Disk& d) {
  layerId_ = d.myid();
}

void TrackerVisitor::visit(const Module& m) {
  output_ << m.myDetId() << ","
	  << m.myBinaryDetId() << ","
	  << sectionName_ << ", "
	  << layerId_ << ", "
	  << m.moduleRing() << ", "
	  << std::fixed << std::setprecision(6)
	  << m.center().Rho() << ", "
	  << m.center().Z() << ", "
	  << m.tiltAngle() * 180. / M_PI << ", "
	  << m.center().Phi() * 180. / M_PI << ", "
	  << m.meanWidth() << ", "
	  << m.length() << ", "
	  << m.dsDistance() << ", "
	  << m.sensorThickness()
	  << std::endl;
}


    //************************************//
    //*               Visitor             //
    //*            BarrelModulesCsv       //
    //*                                   //
    //************************************//
void BarrelVisitor::preVisit() {
  output_ << "DetId, BinaryDetId, Barrel-Layer name, r(mm), z(mm), tiltAngle(deg), num mods, meanWidth(mm) (orthoradial), length(mm) (along Z), sensorSpacing(mm), sensorThickness(mm)" << std::endl;
}
void BarrelVisitor::visit(const Barrel& b) {
  barName_ = b.myid();
}
void BarrelVisitor::visit(const Layer& l) {
  layId_ = l.myid();
  numRods_ = l.numRods();
}
void BarrelVisitor::visit(const BarrelModule& m) {
  if (m.posRef().phi > 2) return;
  output_ << m.myDetId() << ", "
	  << m.myBinaryDetId() << ","
	  << barName_ << "-L" << layId_ << ", "
	  << std::fixed << std::setprecision(6)
	  << m.center().Rho() << ", "
	  << m.center().Z() << ", "
	  << m.tiltAngle() * 180. / M_PI << ", "
	  << numRods_/2. << ", "
	  << m.meanWidth() << ", "
	  << m.length() << ", "
	  << m.dsDistance() << ", "
	  << m.sensorThickness()
	  << std::endl;
}

std::string BarrelVisitor::output() const { return output_.str(); }
 

    //************************************//
    //*               Visitor             //
    //*            EndcapModulesCsv       //
    //*                                   //
    //************************************//
void EndcapVisitor::preVisit() {
  output_ << "DetId, BinaryDetId, Endcap-Disc name, Ring, r(mm), z(mm), tiltAngle(deg), phi(deg),  meanWidth(mm) (orthoradial), length(mm) (radial), sensorSpacing(mm), sensorThickness(mm)" << std::endl;
}

void EndcapVisitor::visit(const Endcap& e) {
  endcapName_ = e.myid();
}

void EndcapVisitor::visit(const Disk& d)  {
  diskId_ = d.myid();
}

void EndcapVisitor::visit(const EndcapModule& m) {
  if (m.minZ() < 0.) return;

  output_	<< m.myDetId() << ", "
		<< m.myBinaryDetId() << ","
		<< endcapName_ << "-D" << diskId_ << ", "
		<< m.ring() << ", "
		<< std::fixed << std::setprecision(6)
		<< m.center().Rho() << ", "
		<< m.center().Z() << ", "
		<< m.tiltAngle() * 180. / M_PI << ", "
		<< m.center().Phi() * 180. / M_PI << ", "
		<< m.meanWidth() << ", "
		<< m.length() << ", "
		<< m.dsDistance() << ", "
		<< m.sensorThickness()
		<< std::endl;
}

std::string EndcapVisitor::output() const { return output_.str(); }
    

    //************************************//
    //*               Visitor             //
    //*            Sensors DetIds         //
    //*                                   //
    //************************************//
void TrackerSensorVisitor::visit(Barrel& b) {
  sectionName_ = b.myid();
}

void TrackerSensorVisitor::visit(Endcap& e) {
  sectionName_ = e.myid();
}

void TrackerSensorVisitor::visit(Layer& l)  {
  layerId_ = l.myid();
}

void TrackerSensorVisitor::visit(Disk& d)  {
  layerId_ = d.myid();
}

void TrackerSensorVisitor::visit(Module& m)  {
  moduleRing_ = m.moduleRing();
}

void TrackerSensorVisitor::visit(Sensor& s) {
  output_ << s.myDetId() << ","
	  << s.myBinaryDetId() << ","
	  << sectionName_ << ", "
	  << layerId_ << ", "
	  << moduleRing_ << ", "
	  << std::fixed << std::setprecision(6)
	  << s.hitPoly().getCenter().Rho() << ", "
	  << s.hitPoly().getCenter().Z() << ", "
	  << s.hitPoly().getCenter().Phi() * 180. / M_PI
	  << std::endl;
}

std::string TrackerSensorVisitor::output() const { return output_.str(); }
   

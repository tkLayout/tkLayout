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
  endcapId = e.myid();
  diskTable->setContent(0, 1 + nDisks, endcapId);

  RootWTable* endcapName = new RootWTable();
  endcapName->setContent(0, 0, endcapId + ",  Disc 1 :");
  endcapNames.push_back(endcapName);

  RootWTable* endcapTable = new RootWTable();
  endcapTable->setContent(0, 0, "Ring :");
  endcapTable->setContent(1, 0, "r"+subStart+"min"+subEnd);
  endcapTable->setContent(2, 0, "r"+subStart+"low"+subEnd);
  endcapTable->setContent(3, 0, "r"+subStart+"centre"+subEnd);
  endcapTable->setContent(4, 0, "r"+subStart+"high"+subEnd);
  endcapTable->setContent(5, 0, "r"+subStart+"max"+subEnd);
  endcapTable->setContent(6, 0, "phiOverlap");
  endcapTable->setContent(7, 0, "# mods");
  endcapTables.push_back(endcapTable);
}

void LayerDiskSummaryVisitor::visit(const Disk& d) {
  nRings = 0;
  nRingsTotal = d.numRings() - d.numEmptyRings();
  if (d.centerZ() < 0.) return;
  ++nDisks;
  totalEndcapModules += d.totalModules();
  diskTable->setContent(1, nDisks, d.myid());
  diskTable->setContent(2, nDisks, d.centerZ(), coordPrecision);
  diskTable->setContent(4, nDisks, d.totalModules());

  RootWTable* diskName = new RootWTable();
  diskId = endcapId + ",  Disc " + std::to_string(d.myid()) + " :";
  diskName->setContent(0, 0, diskId);
  diskNames.push_back(diskName);

  RootWTable* zErrorTable = new RootWTable();
  zErrorTable->setContent(0, 0, "Ring :");
  zErrorTable->setContent(1, 0, "zError (Ring i & i+1)");
  zErrorTables.push_back(zErrorTable);
}

void LayerDiskSummaryVisitor::visit(const Ring& r) {
  if (r.averageZ() < 0. || r.numModules() == 0) return;
  ++nRings;
  diskTable->setContent(3, nDisks, nRings);
  endcapTables.at(nEndcaps-1)->setContent(0, nRings, r.myid());
  endcapTables.at(nEndcaps-1)->setContent(6, nRings, r.actualPhiOverlap(), coordPrecision);
  endcapTables.at(nEndcaps-1)->setContent(7, nRings, r.numModules());
  zErrorTables.at(nDisks-1)->setContent(0, nRings, r.myid());
  if (nRings != nRingsTotal) zErrorTables.at(nDisks-1)->setContent(1, nRings, r.actualZError(), coordPrecision);
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
  totalSensorPower += m.sensorsIrradiationPowerMean();
  if (tagMap.find(aSensorTag)==tagMap.end()){
    // We have a new sensor geometry
    tagMap[aSensorTag] = &m;
  }
  // Summing sensor areas calculated from their respective polygons
  for (int i = 0; i < m.numSensors(); i++) {
    totArea += m.sensors().at(i).hitPoly().getDoubleArea()/2.;
  }
}

void LayerDiskSummaryVisitor::visit(const EndcapModule& m) {
  if (m.side() != 1 || m.disk() != 1) return;

  endcapTables.at(nEndcaps-1)->setContent(1, nRings, m.minR(), coordPrecision);
  endcapTables.at(nEndcaps-1)->setContent(2, nRings, sqrt(pow(m.minR(),2)+pow(m.minWidth()/2.,2)), coordPrecision); // Ugly, this should be accessible as a method
  endcapTables.at(nEndcaps-1)->setContent(3, nRings, m.center().Rho(), coordPrecision);
  endcapTables.at(nEndcaps-1)->setContent(4, nRings, m.minR()+m.length(), coordPrecision);
  endcapTables.at(nEndcaps-1)->setContent(5, nRings, m.maxR(), coordPrecision);
}

void LayerDiskSummaryVisitor::postVisit() {
  layerTable->setContent(0, nBarrelLayers+1, "Total");
  layerTable->setContent(5, nBarrelLayers+1, totalBarrelModules);
  diskTable->setContent(0, nDisks+1, "Total");
  diskTable->setContent(4, nDisks+1, totalEndcapModules*2);
}


    //***************************************//
    //*                Visitor              *//
    //*             Skewed layers:          *//
    //*            Additional info          *//
    //*                                     *//
    //***************************************//
/**
   Visits skewed layers.
*/
void SkewedLayersVisitor::visit(const Layer& l) {
  // Only for skewed layers.
  if (l.isSkewedForInstallation()) {
    numSkewedLayers++;

    // FILLS LAYER TABLE
    RootWTable* layerTable = new RootWTable();
    layerTable->setContent(0, 0, "Layer " + std::to_string(l.myid()) + " :");
    layerTable->setContent(1, 0, "layer Rho [mm]");
    layerTable->setContent(1, 1, l.placeRadiusHint(), coordPrecision);
    layerTable->setContent(2, 0, "bigDelta [mm]");
    layerTable->setContent(2, 1, l.bigDelta(), coordPrecision);
    layerTable->setContent(3, 0, "smallDelta [mm]");
    layerTable->setContent(3, 1, l.smallDelta(), coordPrecision);
    layerTable->setContent(4, 0, "inner module Rho" + subStart + "center" + subEnd + " [mm]");
    layerTable->setContent(4, 1, l.placeRadiusHint() - l.bigDelta(), coordPrecision);
    layerTable->setContent(5, 0, "outer module Rho" + subStart + "center" + subEnd + " [mm]");
    layerTable->setContent(5, 1, l.placeRadiusHint() + l.bigDelta(), coordPrecision);
    layerTable->setContent(6, 0, "skewed module Rho" + subStart + "center" + subEnd + " [mm]");
    layerTable->setContent(6, 1, l.skewedModuleCenterRho(), coordPrecision);
    layerTable->setContent(7, 0, "skewed module Rho" + subStart + "min" + subEnd + " [mm]");
    layerTable->setContent(7, 1, l.skewedModuleMinRho(), coordPrecision);
    layerTable->setContent(8, 0, "skewed module Rho" + subStart + "max" + subEnd + " [mm]");
    layerTable->setContent(8, 1, l.skewedModuleMaxRho(), coordPrecision);
    layerTable->setContent(9, 0, "skewed module edge Shift [mm]");
    layerTable->setContent(9, 1, l.skewedModuleEdgeShift(), coordPrecision);
    layerTable->setContent(10, 0, "skew angle [°]");
    layerTable->setContent(10, 1, l.skewAngle() * 180. / M_PI, coordPrecision);
    layerTable->setContent(11, 0, "phiOverlap [mm]");
    layerTable->setContent(11, 1, l.unitPhiOverlapLength(), coordPrecision);
    layerTable->setContent(12, 0, "horizontal overlap with (X=0) plane (each IT half) [mm]");
    layerTable->setContent(12, 1, l.installationHorizontalOverlapLength(), coordPrecision);
    layerTable->setContent(13, 0, "phiOverlap angular Ratio");
    layerTable->setContent(13, 1, l.installationOverlapRatio(), coordPrecision);
   
    tables.push_back(layerTable);
  } // end of 'fills layer table'
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
    const int numTiltedRings = l.tiltedRingsGeometry().size();
    for (int i=0; i < numTiltedRings; i++) {
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
  output_ << "DetId/U, BinaryDetId/B, Section/C, Layer/I, Ring/I, SensorCenter rho(mm), SensorCenter z(mm), tiltAngle_deg/D, skewAngle_deg/D, yawAngle_deg/D, phi_deg/D, vtxOneX_mm/D, vtxOneY_mm/D,vtxTwoX_mm/D,vtxTwoY_mm/D,vtxThreeX_mm/D,vtxThreeY_mm/D,vtxFourX_mm/D,vtxFourY_mm/D, meanWidth_mm/D, length_mm/D, sensorSpacing_mm/D, sensorThickness_mm/D" << std::endl;
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
	  << m.skewAngle() * 180. / M_PI << ", "
	  << m.yawAngle() * 180. / M_PI << ", "
	  << m.center().Phi() * 180. / M_PI << ", "
          << m.getVertex(0).X() << ", "
          << m.getVertex(0).Y() << ", "
          << m.getVertex(1).X() << ", "
          << m.getVertex(1).Y() << ", "
          << m.getVertex(2).X() << ", "
          << m.getVertex(2).Y() << ", "
          << m.getVertex(3).X() << ", "
          << m.getVertex(3).Y() << ", "
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
  output_ << "DetId, BinaryDetId, Barrel-Layer name, SensorCenter rho(mm), SensorCenter z(mm), tiltAngle(deg), num mods, meanWidth(mm) (orthoradial), length(mm) (along Z), sensorSpacing(mm), sensorThickness(mm)" << std::endl;
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
  output_ << "DetId, BinaryDetId, Endcap-Disc name, Ring, SensorCenter rho(mm), SensorCenter z(mm), tiltAngle(deg), yawAngle(deg), phi(deg), vtxOneX_mm/D, vtxOneY_mm/D,vtxTwoX_mm/D,vtxTwoY_mm/D,vtxThreeX_mm/D,vtxThreeY_mm/D,vtxFourX_mm/D,vtxFourY_mm/D, meanWidth(mm) (orthoradial), length(mm) (radial), sensorSpacing(mm), sensorThickness(mm)" << std::endl;
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
		<< m.yawAngle() * 180. / M_PI << ", "
		<< m.center().Phi() * 180. / M_PI << ", "
                << m.getVertex(0).X() << ", "
                << m.getVertex(0).Y() << ", "
                << m.getVertex(1).X() << ", "
                << m.getVertex(1).Y() << ", "
                << m.getVertex(2).X() << ", "
                << m.getVertex(2).Y() << ", "
                << m.getVertex(3).X() << ", "
                << m.getVertex(3).Y() << ", "
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
   

    //************************************//
    //*               Visitor             //
    //*            ModulesToDTCsCsv       //
    //*                                   //
    //************************************//
ModulesToDTCsVisitor::ModulesToDTCsVisitor(bool isPositiveCablingSide) {
  isPositiveCablingSide_ = isPositiveCablingSide;
}

void ModulesToDTCsVisitor::preVisit() {
  output_ << "Module_DetId/i, Module_Section/C, Module_Layer/I, Module_Ring/I, Module_phi_deg/D, MFB/I, OPT_Services_Channel/I, PWR_Services_Channel/I, MFC/I, MFC_type/C, DTC_name/C, DTC_CMSSW_Id/i, DTC_Phi_Sector_Ref/I, type_/C, DTC_Slot/I, DTC_Phi_Sector_Width_deg/D" << std::endl;
}

void ModulesToDTCsVisitor::visit(const Barrel& b) {
  sectionName_ = b.myid();
}

void ModulesToDTCsVisitor::visit(const Endcap& e) {
  sectionName_ = e.myid();
}

void ModulesToDTCsVisitor::visit(const Layer& l) {
  layerId_ = l.myid();
}

void ModulesToDTCsVisitor::visit(const Disk& d) {
  layerId_ = d.myid();
}

void ModulesToDTCsVisitor::visit(const Module& m) {
  const OuterBundle* myBundle = m.getBundle();
  if (myBundle != nullptr) {
    if (myBundle->isPositiveCablingSide() == isPositiveCablingSide_) {
      std::stringstream moduleInfo;
      moduleInfo << m.myDetId() << ", "
		 << sectionName_ << ", "
		 << layerId_ << ", "
		 << m.moduleRing() << ", "
		 << std::fixed << std::setprecision(6)
		 << m.center().Phi() * 180. / M_PI << ", ";

      std::stringstream bundleInfo;
      bundleInfo << myBundle->myid() << ", ";

      const OuterCable* myCable = myBundle->getCable();
      if (myCable != nullptr) {
	std::stringstream cableInfo;
	cableInfo << myCable->myid() << ", "
		  << any2str(myCable->type()) << ", ";
	bundleInfo << myCable->opticalChannelSection()->channelNumber() << " " 
		   << any2str(myCable->opticalChannelSection()->channelSlot()) << ", "
		   << myBundle->powerChannelSection()->channelNumber() << " " 
		   << any2str(myBundle->powerChannelSection()->channelSlot()) << ", ";
	
	const OuterDTC* myDTC = myCable->getDTC();
	if (myDTC != nullptr) {
	  std::stringstream DTCInfo;
	  DTCInfo << myDTC->name() << ", "
		  << myDTC->getCMSSWId() << ", "
		  << myDTC->phiSectorRef() << ", "
		  << any2str(myDTC->type()) << ", "
		  << myDTC->slot() << ", "
		  << std::fixed << std::setprecision(6)
		  << myDTC->phiSectorWidth() * 180. / M_PI;
	  output_ << moduleInfo.str() << bundleInfo.str() << cableInfo.str() << DTCInfo.str() << std::endl;
	}
	else output_ << moduleInfo.str() << bundleInfo.str() << cableInfo.str() << std::endl;
      }
      else output_ << moduleInfo.str() << bundleInfo.str() << std::endl;
    }
  }
}


    //************************************//
    //*               Visitor             //
    //*     CMSSWOuterTrackerCablingMap   //
    //*                                   //
    //************************************//
void CMSSWOuterTrackerCablingMapVisitor::preVisit() {
  output_ << "Module_DetId/U, GBT_CMSSW_IdPerDTC/U, DTC_CMSSW_Id/U" << std::endl;
}

void CMSSWOuterTrackerCablingMapVisitor::visit(const Module& m) {
  std::stringstream moduleInfo;
  moduleInfo << m.myDetId() << ", ";
  const OuterGBT* myGBT = m.getOuterGBT();
  if (myGBT) {
    std::stringstream GBTInfo;
    GBTInfo << myGBT->getCMSSWId() << ", ";
    const OuterDTC* myDTC = m.getDTC();
    if (myDTC) {
      std::stringstream DTCInfo;
      DTCInfo << myDTC->getCMSSWId();	 
      output_ << moduleInfo.str() << GBTInfo.str() << DTCInfo.str() << std::endl;
    }
  }
  else output_ << moduleInfo.str() << std::endl;
}


    //************************************//
    //*               Visitor             //
    //*    InnerTrackerModulesToDTCsCsv   //
    //*                                   //
    //************************************//

void InnerTrackerModulesToDTCsVisitor::preVisit() {
  output_ << "Module_DetId/i, Module_Section/C, Module_Layer/I, Module_Ring/I, Module_phi_deg/D, N_Chips_Per_Module/I, N_Channels_Per_Module/I, Is_LongBarrel/O, Power_Chain/I, Power_Chain_Type/C, N_ELinks_Per_Module/I, LpGBT_Id/C, LpGBT_CMSSW_IdPerDTC/U, MFB/I, DTC_Id/I, DTC_CMSSW_Id/U, IsPlusZEnd/O, IsPlusXSide/O" << std::endl;
}

void InnerTrackerModulesToDTCsVisitor::visit(const Barrel& b) {
  sectionName_ = b.myid();
}

void InnerTrackerModulesToDTCsVisitor::visit(const Endcap& e) {
  sectionName_ = e.myid();
}

void InnerTrackerModulesToDTCsVisitor::visit(const Layer& l) {
  layerId_ = l.myid();
}

void InnerTrackerModulesToDTCsVisitor::visit(const Disk& d) {
  layerId_ = d.myid();
}

void InnerTrackerModulesToDTCsVisitor::visit(const Module& m) {
  const PowerChain* myPowerChain = m.getPowerChain();
  if (myPowerChain != nullptr) {
    std::stringstream moduleInfo;
    moduleInfo << m.myDetId() << ","
	       << sectionName_ << ", "
	       << layerId_ << ", "
	       << m.moduleRing() << ", "
	       << std::fixed << std::setprecision(6)
	       << m.center().Phi() * 180. / M_PI << ", "
	       << m.outerSensor().totalROCs() << ", "
	       << m.totalChannels() << ", ";

    std::stringstream powerChainInfo;
    powerChainInfo << any2str(myPowerChain->isLongBarrel()) << ","
		   << myPowerChain->myid() << ","
		   << any2str(myPowerChain->powerChainType()) << ",";

    const GBT* myGBT = m.getGBT();
    if (myGBT != nullptr) {
      std::stringstream GBTInfo;
      GBTInfo << myGBT->numELinksPerModule() << ","
	      << any2str(myGBT->GBTId()) << ","
	      << myGBT->getCMSSWId() << ",";

      const InnerBundle* myBundle = myGBT->getBundle();
      if (myBundle != nullptr) {
	std::stringstream bundleInfo;
	bundleInfo << myBundle->myid() << ",";
	
	const InnerDTC* myDTC = myBundle->getDTC();
	if (myDTC != nullptr) {
	  std::stringstream DTCInfo;
	  DTCInfo << myDTC->myid() << ","
		  << myDTC->getCMSSWId() << ","
		  << myDTC->isPositiveZEnd() << ","
		  << myDTC->isPositiveXSide();
	  output_ << moduleInfo.str() << powerChainInfo.str() << GBTInfo.str() << bundleInfo.str() << DTCInfo.str() << std::endl;
	}
	else output_ << moduleInfo.str() << powerChainInfo.str() << GBTInfo.str() << bundleInfo.str() << std::endl;
      }
      else output_ << moduleInfo.str() << powerChainInfo.str() << GBTInfo.str() << std::endl;
    }
    else output_ << moduleInfo.str() << powerChainInfo.str() << std::endl;
  }
}


    //************************************//
    //*               Visitor             //
    //*     CMSSWInnerTrackerCablingMap   //
    //*                                   //
    //************************************//
void CMSSWInnerTrackerCablingMapVisitor::preVisit() {
  output_ << "Module_DetId/U, GBT_CMSSW_IdPerDTC/U, DTC_CMSSW_Id/U" << std::endl;
}

void CMSSWInnerTrackerCablingMapVisitor::visit(const Module& m) {
  std::stringstream moduleInfo;
  moduleInfo << m.myDetId() << ", ";

  const GBT* myGBT = m.getGBT();
  if (myGBT) {
    std::stringstream GBTInfo;
    GBTInfo << myGBT->getCMSSWId() << ",";

    const InnerDTC* myDTC = m.getInnerDTC();
    if (myDTC) {
      std::stringstream DTCInfo;
      DTCInfo << myDTC->getCMSSWId();	 
      output_ << moduleInfo.str() << GBTInfo.str() << DTCInfo.str() << std::endl;
    }
    else output_ << moduleInfo.str() << GBTInfo.str() << std::endl;
  }
  else output_ << moduleInfo.str() << std::endl;
}

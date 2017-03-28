
#include "DetectorModule.hh"
#include "ModuleCap.hh"

void DetectorModule::build() {
  check();

  if (!decorated().builtok()) {
    decorated().store(propertyTree());
    decorated().build();
  }
  if (numSensors() > 0) {
    for (int i = 0; i < numSensors(); i++) {
      Sensor* s = GeometryFactory::make<Sensor>();
      s->parent(this);
      s->myid(i+1);
      s->store(propertyTree());
      if (sensorNode.count(i+1) > 0) s->store(sensorNode.at(i+1));
      if (numSensors() == 1) s->innerOuter(SensorPosition::NO);
      else {
	if (i == 0) s->innerOuter(SensorPosition::LOWER);
	else if (i == 1) s->innerOuter(SensorPosition::UPPER);
	else s->innerOuter(SensorPosition::NO);
      }
      s->build();
      sensors_.push_back(s);
      materialObject_.sensorChannels[i+1]=s->numChannels();
    }
  } else {
    Sensor* s = GeometryFactory::make<Sensor>();  // fake sensor to avoid defensive programming when iterating over the sensors and the module is empty
    s->parent(this);
    s->myid(1);
    s->innerOuter(SensorPosition::NO);
    s->build();
    sensors_.push_back(s);
  }

  materialObject_.store(propertyTree());
  materialObject_.build();
}


void DetectorModule::setup() {
  nominalResolutionLocalX.setup([this]() {
      // only set up this if no model parameter specified
      //std::cout <<  "hasAnyResolutionLocalXParam() = " <<  hasAnyResolutionLocalXParam() << std::endl;

      if (!hasAnyResolutionLocalXParam()) {
	//std::cout << "nominalResolutionLocalX and resolutionLocalXBarrel parameters are all unset. Use of default formulae." << std::endl;
	double res = 0;
	for (const Sensor& s : sensors()) res += pow(meanWidth() / s.numStripsAcrossEstimate() / sqrt(12), 2);
	return sqrt(res)/numSensors();
      }
      // if model parameters specified, return -1
      else return -1.0;
    });
  nominalResolutionLocalY.setup([this]() {
      // only set up this if no model parameters not specified
      if (!hasAnyResolutionLocalYParam()) {
	//std::cout << "resolutionLocalY and resolutionLocalYBarrel parameters are all unset. Use of default formulae." << std::endl;
	if (stereoRotation() != 0.) return nominalResolutionLocalX() / sin(stereoRotation());
	else {
	  return length() / maxSegments() / sqrt(12); // NOTE: not combining measurements from both sensors. The two sensors are closer than the length of the longer sensing element, making the 2 measurements correlated. considering only the best measurement is then a reasonable approximation (since in case of a PS module the strip measurement increases the precision by only 0.2% and in case of a 2S the sensors are so close that they basically always measure the same thing)
	}
      }
      // if model parameters specified, return -1
      else return -1.0;
    });
  for (Sensor& s : sensors_) s.parent(this); // set the parent for the sensors once again (in case the module's been cloned)
};






std::map<std::string, double> DetectorModule::extremaWithHybrids() const {

    std::map<std::string, double> extrema;

    //                                                   OUTER TRACKER MODULE
    //
    //  Top View
    //  ------------------------------
    //  |            L(5)            |  
    //  |----------------------------|     y
    //  |     |                |     |     ^
    //  |B(4) |     Between    | F(3)|     |
    //  |     |       (7)      |     |     +----> x
    //  |----------------------------|
    //  |            R(6)            |     
    //  ------------------------------     
    //                                            z
    //  Side View                                 ^
    //         ---------------- OuterSensor(2)    |
    //  ====== ================ ====== Hybrids    +----> x
    //         ---------------- InnerSensor(1)
    //  ============================== 
    //          SupportPlate(8)                      
    //
    //  R(6) and L(5) are Front-End Hybrids.
    //  B(4) and F(3) are Service Hybdrids.
    //
    //  SupportPlate(8) thickness is of course null for 2S modules.


    //                                                      PIXEL MODULE
    //
    //  Top View 
    //        ------------------           y
    //        |                |           ^
    //        |     Hybrid     |           |
    //        |       (1)      |           +----> x
    //        ------------------    
    //                                             z
    //                                             ^
    //  Side View                                  |
    //         ================ Hybrid  (1)        +----> x
    //         ---------------- Sensor  (2)
    //         ================ Chip    (3)
    //
    // Chip(3) volume can contain Bumps and any other material for simplification.



    // =========================================================================================================
    // Finding Xmin/Xmax/Ymin/Ymax/Zmin/Zmax/Rmin/Rmax/RminatZmin/RmaxatZmax, taking hybrid volumes into account
    // =========================================================================================================
    //
    // Module polygon
    //   top view
    //   v1                v2
    //    *---------------*
    //    |       ^ my    |
    //    |       |   mx  |
    //    |       *------>|
    //    |     center    |
    //    |               |
    //    *---------------*
    //   v0                v3
    //  (v4)
    //
    //   side view
    //    ----------------- top
    //    ----------------- bottom


    double         rmin;
    double         rmax;
    double         xmin;
    double         xmax;
    double         ymin;
    double         ymax;
    double         zmin;
    double         zmax;
    double         rminatzmin;
    double         rmaxatzmax;
    std::vector<XYZVector> vertex; 



    double width = area() / length();
    double expandedModWidth = width + 2 * serviceHybridWidth();
    double expandedModLength = length() + 2 * frontEndHybridWidth();
    double expandedModThickness;
    if (!isPixelModule()) { expandedModThickness = dsDistance() + 2.0 * (supportPlateThickness() + sensorThickness()); }
    //double expandedModThickness = dsDistance() + supportPlateThickness()+sensorThickness(); SHOULD BE THIS !!!!
    else { expandedModThickness = sensorThickness() + 2.0 * MAX(chipThickness(), hybridThickness()); }

    
    vector<double> xv; // x list (in global frame of reference) from which we will find min/max.
    vector<double> yv; // y list (in global frame of reference) from which we will find min/max.
    vector<double> zv; // z (in global frame of reference) list from which we will find min/max.
    vector<double> rv; // radius list (in global frame of reference) from which we will find min/max.
    vector<double> ratzminv; // radius list (in global frame of reference) at zmin from which we will find min/max.
    vector<double> ratzmaxv; // radius list (in global frame of reference) at zmax from which we will find min/max.

    // mx: (v2+v3)/2 - center(), my: (v1+v2)/2 - center()
    XYZVector mx = 0.5*( basePoly().getVertex(2) + basePoly().getVertex(3) ) - center() ;
    XYZVector my = 0.5*( basePoly().getVertex(1) + basePoly().getVertex(2) ) - center() ;

    // new vertexes after expansion due to hybrid volumes
    const int npoints = 5; // v0,v1,v2,v3,v4(=v0)
    XYZVector v[npoints-1];
    v[0] = center() - (expandedModWidth / width ) * mx - (expandedModLength / length()) * my;
    v[1] = center() - (expandedModWidth / width ) * mx + (expandedModLength / length()) * my;
    v[2] = center() + (expandedModWidth / width ) * mx + (expandedModLength / length()) * my;
    v[3] = center() + (expandedModWidth / width ) * mx - (expandedModLength / length()) * my;

    // Calculate all vertex candidates (8 points)
    XYZVector v_top[npoints];    // module's top surface
    XYZVector v_bottom[npoints]; // module's bottom surface

    for (int ip = 0; ip < npoints-1; ip++) {
      v_top[ip]    = v[ip] + 0.5*expandedModThickness * normal();
      v_bottom[ip] = v[ip] - 0.5*expandedModThickness * normal();

      // for debuging
      vertex.push_back(v_top[ip]);
      vertex.push_back(v_bottom[ip]);

      // Calculate xmin, xmax, ymin, ymax, zmin, zmax
      xv.push_back(v_top[ip].X());
      xv.push_back(v_bottom[ip].X());
      yv.push_back(v_top[ip].Y());
      yv.push_back(v_bottom[ip].Y());
      zv.push_back(v_top[ip].Z());
      zv.push_back(v_bottom[ip].Z());
    }
    // Find min and max
    xmin = *std::min_element(xv.begin(), xv.end());
    xmax = *std::max_element(xv.begin(), xv.end());
    ymin = *std::min_element(yv.begin(), yv.end());
    ymax = *std::max_element(yv.begin(), yv.end());
    zmin = *std::min_element(zv.begin(), zv.end());
    zmax = *std::max_element(zv.begin(), zv.end());


    // Calculate module's mid-points (8 points)
    XYZVector v_mid_top[npoints]; // module's top surface mid-points
    XYZVector v_mid_bottom[npoints]; // module's bottom surface mid-points

    v_top[npoints-1] = v_top[0]; // copy v0 as v4 for convenience
    v_bottom[npoints-1] = v_bottom[0]; // copy v0 as v4 for convenience

    for (int ip = 0; ip < npoints-1; ip++) {
      v_mid_top[ip] = (v_top[ip] + v_top[ip+1]) / 2.0;
      v_mid_bottom[ip] = (v_bottom[ip] + v_bottom[ip+1]) / 2.0;
    }

    // Calculate rmin, rmax, rminatzmin, rmaxatzmax...
    for (int ip = 0; ip < npoints-1; ip++) {

      // module's bottom surface
      if (fabs(v_bottom[ip].Z() - zmin) < 0.001) {
	v_bottom[ip].SetZ(0.); // projection to xy plan.
	ratzminv.push_back(v_bottom[ip].R());
      }
      if (fabs(v_bottom[ip].Z() - zmax) < 0.001) {
	v_bottom[ip].SetZ(0.); // projection to xy plan.
	ratzmaxv.push_back(v_bottom[ip].R());
      }
      v_bottom[ip].SetZ(0.); // projection to xy plan.
      rv.push_back(v_bottom[ip].R());

      // module's top surface
      if (fabs(v_top[ip].Z() - zmin) < 0.001) {
	v_top[ip].SetZ(0.); // projection to xy plan.
	ratzminv.push_back(v_top[ip].R());
      }
      if (fabs(v_top[ip].Z() - zmax) < 0.001) {
	v_top[ip].SetZ(0.); // projection to xy plan.
	ratzmaxv.push_back(v_top[ip].R());
      }
      v_top[ip].SetZ(0.); // projection to xy plan.
      rv.push_back(v_top[ip].R());
    
      // module's bottom surface mid-points
      if (fabs(v_mid_bottom[ip].Z() - zmin) < 0.001) {
	v_mid_bottom[ip].SetZ(0.); // projection to xy plan.
	ratzminv.push_back(v_mid_bottom[ip].R());
      }
      v_mid_bottom[ip].SetZ(0.); // projection to xy plan.
      rv.push_back(v_mid_bottom[ip].R());
    
      // module's top surface mid-points
      if (fabs(v_mid_top[ip].Z() - zmin) < 0.001) {
	v_mid_top[ip].SetZ(0.); // projection to xy plan.
	ratzminv.push_back(v_mid_top[ip].R());
      }
      v_mid_top[ip].SetZ(0.); // projection to xy plan.
      rv.push_back(v_mid_top[ip].R());
    }
    // Find min and max
    rmin = *std::min_element(rv.begin(), rv.end());
    rmax = *std::max_element(rv.begin(), rv.end());
    rminatzmin = *std::min_element(ratzminv.begin(), ratzminv.end());
    rmaxatzmax = *std::max_element(ratzmaxv.begin(), ratzmaxv.end());


    extrema["minZ"] = zmin;
    extrema["maxZ"] = zmax;
    extrema["minR"] = rmin;
    extrema["maxR"] = rmax;

    return extrema;

  }







std::pair<double, double> DetectorModule::minMaxEtaWithError(double zError) const {
  if (cachedZError_ != zError) {
    cachedZError_ = zError;
    double eta1 = (XYZVector(0., maxR(), maxZ() + zError)).Eta();
    double eta2 = (XYZVector(0., minR(), minZ() - zError)).Eta();
    double eta3 = (XYZVector(0., minR(), maxZ() + zError)).Eta();
    double eta4 = (XYZVector(0., maxR(), minZ() - zError)).Eta();
    cachedMinMaxEtaWithError_ = std::minmax({eta1, eta2, eta3, eta4});
    //cachedMinMaxEtaWithError_ = std::make_pair(MIN(eta1, eta2), MAX(eta1, eta2));
  }
  return cachedMinMaxEtaWithError_;
}

bool DetectorModule::couldHit(const XYZVector& direction, double zError) const {

  double eta       = direction.Eta();
  double phi       = direction.Phi();
  double shiftPhi  = phi + 2*M_PI;
  bool   withinEta = false;
  bool   withinPhi = false;

  // Eta region covered by module
  if (eta > minEtaWithError(zError) && eta < maxEtaWithError(zError)) withinEta = true;

  // Phi region is from <-pi;+3*pi> due to crossline at +pi -> need to check phi & phi+2*pi
  if ( (phi     >=minPhi() && phi     <=maxPhi()) ||
       (shiftPhi>=minPhi() && shiftPhi<=maxPhi()) ) withinPhi = true;

  // Checking that hit within a module region works for barrel-type modules only!!!
  if (this->shape()==ModuleShape::RECTANGULAR) return (withinEta && withinPhi);
  // ATTENTION: For wedge shaped modules, min, max procedure will not work correctly -> return true to avoid errors --> will be implemented in the future
  else return true;
}

//bool DetectorModule::couldHit(const XYZVector& direction, double zError) const {
//  double eta = direction.Eta(), phi = direction.Phi();
//  bool withinEta = eta > minEtaWithError(zError) && eta < maxEtaWithError(zError);
//  bool withinPhi;
//  if (minPhi() < 0. && maxPhi() > 0. && maxPhi()-minPhi() > M_PI) // across PI
//    withinPhi = phi < minPhi() || phi > maxPhi();
//  else 
//    withinPhi = phi > minPhi() && phi < maxPhi();
//  //bool withinPhiSub = phi-2*M_PI > minPhi() && phi-2*M_PI < maxPhi();
//  //bool withinPhiAdd = phi+2*M_PI > minPhi() && phi+2*M_PI < maxPhi();
//  return withinEta && (withinPhi /*|| withinPhiSub || withinPhiAdd*/);
//}


double DetectorModule::resolutionEquivalentRPhi(double hitRho, double trackR, double resolutionLocalX, double resolutionLocalY) const {
  double A = hitRho/(2*trackR); 
  double B = A/sqrt(1-A*A);
  return sqrt(pow((B*sin(skewAngle())*cos(tiltAngle()) + cos(skewAngle())) * resolutionLocalX,2) + pow(B*sin(tiltAngle()) * resolutionLocalY,2));
}

double DetectorModule::resolutionEquivalentZ(double hitRho, double trackR, double trackCotgTheta, double resolutionLocalX, double resolutionLocalY) const {
  double A = hitRho/(2*trackR); 
  double D = trackCotgTheta/sqrt(1-A*A);
  return sqrt(pow(((D*cos(tiltAngle()) + sin(tiltAngle()))*sin(skewAngle())) * resolutionLocalX,2) + pow((D*sin(tiltAngle()) + cos(tiltAngle())) * resolutionLocalY,2));
}


double DetectorModule::totalPower() const {
  double result = powerPerModule();
  for (const auto& sen : sensors()) {
    result += sen.powerPerChannel() * sen.numChannels();
  }
  return result;
}


double DetectorModule::stripOccupancyPerEventBarrel() const {
  double rho = center().Rho()/10.;
  double theta = center().Theta();
  double myOccupancyBarrel=(1.63e-4)+(2.56e-4)*rho-(1.92e-6)*rho*rho;
  double factor = fabs(sin(theta))*2; // 2 is a magic adjustment factor
  double dphideta = phiAperture() * etaAperture();
  double minNSegments = minSegments();
  int numStripsAcrossEstimate = sensors().begin()->numStripsAcrossEstimate();
  double modWidth = (maxWidth() + minWidth())/2.;

  double occupancy = myOccupancyBarrel / factor / (90/1e3) * (dphideta / minNSegments) * (modWidth / numStripsAcrossEstimate);

  return occupancy;
}

double DetectorModule::stripOccupancyPerEventEndcap() const {
  double rho = center().Rho()/10.;
  double theta = center().Theta();
  double z = center().Z()/10.;
  double myOccupancyEndcap = (-6.20e-5)+(1.75e-4)*rho-(1.08e-6)*rho*rho+(1.50e-5)*(z);
  double factor=fabs(cos(theta))*2; // 2 is a magic adjustment factor
  double dphideta = phiAperture() * etaAperture();
  double minNSegments = minSegments();
  int numStripsAcrossEstimate = sensors().begin()->numStripsAcrossEstimate();
  double modWidth = (maxWidth() + minWidth())/2.;

  double occupancy = myOccupancyEndcap / factor / (90/1e3) * (dphideta / minNSegments) * (modWidth / numStripsAcrossEstimate);

  return occupancy;
}

double DetectorModule::stripOccupancyPerEvent() const {
  if (fabs(tiltAngle()) < 1e-3) return stripOccupancyPerEventBarrel();
  else if (fabs(tiltAngle()) - M_PI/2. < 1e-3) return stripOccupancyPerEventEndcap();
  else return stripOccupancyPerEventBarrel()*pow(cos(tiltAngle()),2) + stripOccupancyPerEventEndcap()*pow(sin(tiltAngle()),2);
};

double DetectorModule::geometricEfficiency() const {
  double inefficiency = fabs(dsDistance() / (zCorrelation()==SAMESEGMENT ? outerSensor().stripLength() : length()) / tan(center().Theta()+tiltAngle())); // fabs prevents the inefficiency from becoming negative when theta+tilt > 90 deg (meaning the geometrically inefficient area is on the other end of the module)
  return 1-inefficiency;
}

double DetectorModule::effectiveDsDistance() const {
  if (fabs(tiltAngle()) < 1e-3) return dsDistance();
  else return dsDistance()*sin(center().Theta())/sin(center().Theta()+tiltAngle());
}

std::pair<XYZVector, HitType> DetectorModule::checkTrackHits(const XYZVector& trackOrig, const XYZVector& trackDir) {
  HitType ht = HitType::NONE;
  XYZVector gc; // global coordinates of the hit
  if (numSensors() == 1) {
    auto segm = innerSensor().checkHitSegment(trackOrig, trackDir);
    // <SMe>The following line used to return HitType::BOTH. Changing to INNER in order to avoid double hit counting</SMe>
    if (segm.second > -1) { gc = segm.first; ht = HitType::INNER; } 
  } else {
    auto inSegm = innerSensor().checkHitSegment(trackOrig, trackDir);
    auto outSegm = outerSensor().checkHitSegment(trackOrig, trackDir);
    if (inSegm.second > -1 && outSegm.second > -1) { 
      gc = inSegm.first; // in case of both sensors are hit, the inner sensor hit coordinate is returned
      ht = ((zCorrelation() == SAMESEGMENT && (inSegm.second / (maxSegments()/minSegments()) == outSegm.second)) || zCorrelation() == MULTISEGMENT) ? HitType::STUB : HitType::BOTH;
    } else if (inSegm.second > -1) { gc = inSegm.first; ht = HitType::INNER; }
    else if (outSegm.second > -1) { gc = outSegm.first; ht = HitType::OUTER; }
  }
  //basePoly().isLineIntersecting(trackOrig, trackDir, gc); // this was just for debug
  if (ht != HitType::NONE) numHits_++;
  return std::make_pair(gc, ht);
};

std::string DetectorModule::summaryType() const  {
  std::string result;
  result+=moduleType();
  if (dsDistance()!=0) result+=" "+any2str(dsDistance(), 1)+" mm";
  return result;
};

std::string DetectorModule::summaryFullType() const  {
  std::string result;
  if (dsDistance()!=0) result=any2str(dsDistance(), 1)+" mm ";
  result+=moduleType();
  if (subdet()==BARREL) {
    if (isTilted()) result += " Tilted";
    else result += " Flat";
    result += " Barrel";
  } else if (subdet()==ENDCAP) {
    result += " Endcap";
  } else {
    std::cerr << "HELP! I am a bit lost here..." << std::endl;
  }
  return result;
};

//BarrelModule::BarrelModule(Decorated* decorated) : DetectorModule(decorated) {
//setup();
                                    // this was already commented
                                    //myModuleCap_ = new ModuleCap(this);
                                    //myModuleCap_->setCategory(MaterialProperties::b_mod);
//}

void BarrelModule::check() {
  PropertyObject::check();

  //std::cout <<  "hasAnyResolutionLocalYParam() = " <<  hasAnyResolutionLocalYParam() << std::endl;

  if (nominalResolutionLocalX.state() && hasAnyResolutionLocalXParam()) throw PathfulException("Only one between resolutionLocalX and resolutionLocalXBarrelParameters can be specified.");

  if (nominalResolutionLocalY.state() && hasAnyResolutionLocalYParam()) throw PathfulException("Only one between resolutionLocalY and resolutionLocalYBarrelParameters can be specified.");
}


void BarrelModule::build() {
  try {
    DetectorModule::build();
    //myModuleCap_->setCategory(MaterialProperties::b_mod);
    decorated().rotateY(M_PI/2);
    rAxis_ = normal();
    tiltAngle_ = 0.;
    skewAngle_ = 0.;
    for (auto& s : sensors_) { s.subdet(ModuleSubdetector::BARREL); }
  }
  catch (PathfulException& pe) { pe.pushPath(*this, myid()); throw; }
  cleanup();
  builtok(true);
}

//EndcapModule::EndcapModule(Decorated* decorated) : DetectorModule(decorated) { 
//setup();
                                         // this was already commented
                                         //myModuleCap_ = new ModuleCap(this);
                                         //myModuleCap_->setCategory(MaterialProperties::e_mod);
//} 


void EndcapModule::check() {
  PropertyObject::check();

 if (nominalResolutionLocalX.state() && hasAnyResolutionLocalXParam()) throw PathfulException("Only one between resolutionLocalX and resolutionLocalXEndcapParameters can be specified.");

 if (nominalResolutionLocalY.state() && hasAnyResolutionLocalYParam()) throw PathfulException("Only one between resolutionLocalY and resolutionLocalYEndcapParameters can be specified.");
}


void EndcapModule::build() {
  try {
    DetectorModule::build();
    //myModuleCap_->setCategory(MaterialProperties::e_mod);
    rAxis_ = (basePoly().getVertex(0) + basePoly().getVertex(3)).Unit();
    tiltAngle_ = M_PI/2.;
    skewAngle_ = 0.;
    for (auto& s : sensors_) { s.subdet(ModuleSubdetector::ENDCAP); }
  }
  catch (PathfulException& pe) { pe.pushPath(*this, myid()); throw; }
  cleanup();
  builtok(true);
}




define_enum_strings(SensorLayout) = { "nosensors", "mono", "pt", "stereo" };
define_enum_strings(ZCorrelation) = { "samesegment", "multisegment" };
define_enum_strings(ReadoutType) = { "strip", "pixel", "pt" };
define_enum_strings(ReadoutMode) = { "binary", "cluster" };


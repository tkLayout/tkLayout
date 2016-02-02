
#include "DetectorModule.h"
#include "ModuleCap.h"

/*
DetectorModule* DetectorModule::assignType(const string& type, DetectorModule* m) {
  using namespace boost::property_tree;
  if (type.empty()) return new TypedModule(m); // CUIDADO default type should be NullTypeModule or something similar, not the parent class!
  ptree pt;
  TypedModule* tmod;
  try {
    info_parser::read_info(join<string>({Paths::stdconfig, Paths::moduleTypes, type}, Paths::sep), pt);
    string layout = pt.get<string>("sensorLayout", "");
    if (layout == "mono") tmod = new SingleSensorModule(m);
    else if (layout == "pt") tmod = new PtModule(m);
    else if (layout == "stereo") tmod = new StereoModule(m);
    else throw InvalidPropertyValue("sensorLayout", layout);
  } 
  catch(info_parser::info_parser_error& e) { throw InvalidPropertyValue("moduleType", type); }
  catch(ptree_bad_path& e) { throw CheckedPropertyMissing("sensorLayout"); }

  tmod->store(pt);

  return tmod;
}
*/


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
      s->build();
      sensors_.push_back(s);
      materialObject_.sensorChannels[i+1]=s->numChannels();
    }
  } else {
    Sensor* s = GeometryFactory::make<Sensor>();  // fake sensor to avoid defensive programming when iterating over the sensors and the module is empty
    s->parent(this);
    s->myid(1);
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


//double DetectorModule::calculateAlphaAndBeta(double hitRho, double trackR) const {
  //}

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
  }
  catch (PathfulException& pe) { pe.pushPath(*this, myid()); throw; }
  cleanup();
  builtok(true);
}




define_enum_strings(SensorLayout) = { "nosensors", "mono", "pt", "stereo" };
define_enum_strings(ZCorrelation) = { "samesegment", "multisegment" };
define_enum_strings(ReadoutType) = { "strip", "pixel", "pt" };
define_enum_strings(ReadoutMode) = { "binary", "cluster" };


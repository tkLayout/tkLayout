#include <algorithm>
#include <cmath>
#include <vector>
#include <array>
#include <map>

#include <Math/Vector3D.h>

#include "DetectorModule.hh"
#include "ModuleCap.hh"
#include "OuterCabling/OuterBundle.hh"
#include "OuterCabling/OuterCable.hh"
#include "OuterCabling/OuterDTC.hh"

#include "InnerCabling/PowerChain.hh"
#include "InnerCabling/HvLine.hh"
#include "InnerCabling/GBT.hh"
#include "InnerCabling/InnerBundle.hh"
#include "InnerCabling/InnerDTC.hh"


void DetectorModule::setup() {
  nominalResolutionLocalX.setup([this]() {
      // only set up this if no model parameter specified
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


void DetectorModule::check() {
  PropertyObject::check();

  // LOCAL X RESOLUTION PARAMETERS
  if (nominalResolutionLocalX.state() && hasAnyResolutionLocalXParam()) throw PathfulException("Only one between resolutionLocalX and resolutionLocalXParameters can be specified.");

  if (hasAnyResolutionLocalXParam()) {
    if ( !resolutionLocalXParam0.state() 
	 || !resolutionLocalXParam1.state() 
	 || !resolutionLocalXParam2.state() 
	 || !resolutionLocalXParam3.state() 
	 || !resolutionLocalXParam4.state()
	 || !resolutionLocalXParam5.state() 
	 || !resolutionLocalXParam6.state() 
	 || !resolutionLocalXParam7.state()
	 || !resolutionLocalXParam8.state() 
	 || !resolutionLocalXParam9.state()
	 ) throw PathfulException("Local Y spatial resolution, tilted module. Resolution cfg file not properly chosen. You did not assign a resolution cfg file specific to a tilted module.");
  }

  // LOCAL Y RESOLUTION PARAMETERS
  if (nominalResolutionLocalY.state() && hasAnyResolutionLocalYParam()) throw PathfulException("Only one between resolutionLocalY and resolutionLocalYParameters can be specified.");

  if (hasAnyResolutionLocalYParam()) {
    if ( !resolutionLocalYParam0.state() 
	 || !resolutionLocalYParam1.state() 
	 || !resolutionLocalYParam2.state() 
	 || !resolutionLocalYParam3.state() 
	 || !resolutionLocalYParam4.state()
	 || !resolutionLocalYParam5.state() 
	 || !resolutionLocalYParam6.state() 
	 || !resolutionLocalYParam7.state()
	 || !resolutionLocalYParam8.state() 
	 || !resolutionLocalYParam9.state()
	 ) throw PathfulException("Local Y spatial resolution, tilted module. Resolution cfg file not properly chosen. You did not assign a resolution cfg file specific to a tilted module.");
  }

}


void DetectorModule::build() {
  check();

  if (!decorated().builtok()) {
    decorated().store(propertyTree());
    decorated().build();
  }

  if (numSensors() > 0) {
    int nSensors = numSensors();
    if (nSensors > 2) {
      logERROR("More than 2 sensors found. Only the first 2 will be used.");
      nSensors = 2;
      numSensors.scaleByUnit(2.0001 / static_cast<double>(numSensors())); // cheating the read-only property
    }

    for (int i = 0; i < nSensors; i++) {
      // Build the sensor
      Sensor* s = GeometryFactory::make<Sensor>();
      s->parent(this);
      s->myid(i+1);
      s->store(propertyTree());
      if (sensorNode.count(i+1) > 0)
        s->store(sensorNode.at(i+1));
      if (nSensors == 1)
        s->innerOuter(SensorPosition::NO);
      else
        s->innerOuter(i == 0 ? SensorPosition::LOWER : SensorPosition::UPPER);
      s->build();

      // Add it to the module's list of sensors
      sensors_.push_back(s);
      materialObject_.sensorChannels[i+1] = s->numChannels();
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


std::map<std::string, double> DetectorModule::extremaWithHybrids() const {
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
    //   v1                v0
    //    *---------------*
    //    |       ^ my    |
    //    |       |   mx  |
    //    |       *------>|
    //    |     center    |
    //    |               |
    //    *---------------*
    //   v2                v3
    //
    //   side view
    //    ----------------- top
    //    ----------------- bottom

    double width = area() / length();
    double expandedModWidth = width + 2 * serviceHybridWidth();
    double expandedModLength = length() + outerSensorExtraLength() + 2 * frontEndHybridWidth();
    double expandedModThickness;
    if (!isPixelModule())
      expandedModThickness = dsDistance() + sensorThickness() + supportPlateThickness();
    else
      expandedModThickness = sensorThickness() + chipThickness() + hybridThickness();

    const ROOT::Math::XYZVector& modCenter = center();
    const ROOT::Math::XYZVector modX(getLocalX() * expandedModWidth * 0.5);
    const ROOT::Math::XYZVector modY(getLocalY() * expandedModLength * 0.5);
    const ROOT::Math::XYZVector modZ(normal() * expandedModThickness * 0.5);

    // Expanded module corners
    std::array<ROOT::Math::XYZVector, 4> xyCorners;
    for (std::size_t i = 0; i < 4; ++i) {
      double signX = (i == 0 || i == 3) ? 1. : -1;
      double signY = (i == 0 || i == 1) ? 1. : -1;
      xyCorners[i] = modCenter + signX * modX + signY * modY;
    }

    // Calculate all vertex candidates (corners + mid points)
    std::array<ROOT::Math::XYZVector, 16> allPoints;
    for (std::size_t i = 0; i < 4; ++i) {
      std::size_t nextIdx = (i + 1) % 4;

      // Top and bottom corners
      auto botVtx = xyCorners[i] - modZ;
      auto topVtx = xyCorners[i] + modZ;
      allPoints[i*4]     = topVtx;
      allPoints[i*4 + 1] = botVtx;

      // Top and bottom midpoints between this corner and the next
      allPoints[i*4 + 2] = (topVtx + (xyCorners[nextIdx] + modZ)) * 0.5;
      allPoints[i*4 + 3] = (botVtx + (xyCorners[nextIdx] - modZ)) * 0.5;
    }

    // Find min and max Z and R
    const auto [zMin, zMax] = std::minmax_element(allPoints.begin(), allPoints.end(),
      [](const ROOT::Math::XYZVector& a, const ROOT::Math::XYZVector& b) { return a.Z() < b.Z(); });
    const auto [rMin, rMax] = std::minmax_element(allPoints.begin(), allPoints.end(),
      [](const ROOT::Math::XYZVector& a, const ROOT::Math::XYZVector& b) { return a.Rho() < b.Rho(); });

    std::map<std::string, double> extrema {
      {"minZ", zMin->Z()},
      {"maxZ", zMax->Z()},
      {"minR", rMin->Rho()},
      {"maxR", rMax->Rho()}
    };

    return extrema;
  }

bool DetectorModule::couldHit(double trackPhi, double trackSlope, double zError) const {
  // Checking that hit within a module region works for barrel-type modules only!!!
  // ATTENTION: For wedge shaped modules, min, max procedure will not work correctly 
  // -> return true to avoid errors --> will be implemented in the future
  if (shape() != ModuleShape::RECTANGULAR) return true;

  // Phi region is from <-pi;+3*pi> due to crossline at +pi -> need to check phi & phi+2*pi
  const double modMinPhi = minPhi();
  const double modMaxPhi = maxPhi();

  bool validPhi = (trackPhi >= modMinPhi && trackPhi <= modMaxPhi);
  if (!validPhi) {
    const double shiftPhi = trackPhi + 2.0 * M_PI;
    validPhi = (shiftPhi >= modMinPhi && shiftPhi <= modMaxPhi);
  }
  if (!validPhi) return false;

  // Eta covered by module
  const double zMinErr = minZ() - zError;
  const double zMaxErr = maxZ() + zError;
  const double rMin = minR();
  const double rMax = maxR();

  // Calculate the 4 corner slopes and the track's
  const auto [minSlope, maxSlope] = std::minmax({ zMinErr / rMin, zMinErr / rMax, 
                                                  zMaxErr / rMin, zMaxErr / rMax });

  return (trackSlope >= minSlope || trackSlope <= maxSlope);
}


/*
 * RETURN GLOBAL SPATIAL RESOLUTION ON X AXIS.
 */
double DetectorModule::resolutionEquivalentRPhi(double hitRho, double trackR, double resolutionLocalX, double resolutionLocalY) const {
  double A = hitRho/(2*trackR); 
  double B = A/sqrt(1-A*A);
  return sqrt(pow((B*sin(skewAngle())*cos(tiltAngle()) + cos(skewAngle())) * resolutionLocalX,2) + pow(B*sin(tiltAngle()) * resolutionLocalY,2));
}


/*
 * RETURN GLOBAL SPATIAL RESOLUTION ON Y AXIS.
 */
double DetectorModule::resolutionEquivalentZ(double hitRho, double trackR, double trackCotgTheta, double resolutionLocalX, double resolutionLocalY) const {
  double A = hitRho/(2*trackR); 
  double D = trackCotgTheta/sqrt(1-A*A);
  return sqrt(pow(((D*cos(tiltAngle()) + sin(tiltAngle()))*sin(skewAngle())) * resolutionLocalX,2) + pow((D*sin(tiltAngle()) + cos(tiltAngle())) * resolutionLocalY,2));
}


/*
 * RETURN LOCAL SPATIAL RESOLUTION ON X AXIS.
 * This resolution is nominal by default, and parametrized if the paramnetrization model is called.
 */
const double DetectorModule::resolutionLocalX(const TVector3& trackDirection) const {
  if (!hasAnyResolutionLocalXParam()) { return nominalResolutionLocalX(); }
  else { return calculateParameterizedResolutionLocalX(trackDirection); }
}


/*
 * RETURN LOCAL SPATIAL RESOLUTION ON Y AXIS.
 * This resolution is nominal by default, and parametrized if the paramnetrization model is called.
 */
const double DetectorModule::resolutionLocalY(const TVector3& trackDirection) const {
  if (!hasAnyResolutionLocalYParam()) { return nominalResolutionLocalY(); }
  else { return calculateParameterizedResolutionLocalY(trackDirection); }
}


/*
 * Compute parametrized local spatial resolution on X axis.
 */
const double DetectorModule::calculateParameterizedResolutionLocalX(const TVector3& trackDirection) const {
  const double tanLorentzAngle = (is3DPixelModule() ? 0. :
				  0.053 * SimParms::getInstance().magField() * cos(tiltAngle())  // dependancy on tilt angle is here!!! :)
				  );
  const double cotanAlpha = 1./tan(alpha(trackDirection));         // Riccardo's theta = alpha - Pi/2    => than(theta) = -cotan(alpha)
  const double fabsTanDeepAngle = fabs(-cotanAlpha - tanLorentzAngle);  

  const bool isLocalXAxis = true;
  const double resolutionLocalX = calculateParameterizedResolutionLocalAxis(fabsTanDeepAngle, isLocalXAxis);

  return resolutionLocalX;
}


/*
 * Compute parametrized local spatial resolution on Y axis.
 */
const double DetectorModule::calculateParameterizedResolutionLocalY(const TVector3& trackDirection) const {
  const double cotanBeta = 1./tan(beta(trackDirection));           // Riccardo's theta = beta - Pi/2    => than(theta) = -cotan(beta)
  const double fabsTanDeepAngle = fabs(-cotanBeta); 
                 
  const bool isLocalXAxis = false;
  const double resolutionLocalY = calculateParameterizedResolutionLocalAxis(fabsTanDeepAngle, isLocalXAxis);

  return resolutionLocalY;
}


/*
 * Compute spatial resolution on a local axis as a function of |tan(deepAngle)|.
 * This is the parametrization model from Riccardo del Burgo, ETZ.
 * The model itself is the same for both X axis and Y axis.
 */
const double DetectorModule::calculateParameterizedResolutionLocalAxis(const double fabsTanDeepAngle, const bool isLocalXAxis) const {
  // Parameter values
  const double param0 = (isLocalXAxis ? resolutionLocalXParam0() : resolutionLocalYParam0());
  const double param1 = (isLocalXAxis ? resolutionLocalXParam1() : resolutionLocalYParam1());
  const double param2 = (isLocalXAxis ? resolutionLocalXParam2() : resolutionLocalYParam2());
  const double param3 = (isLocalXAxis ? resolutionLocalXParam3() : resolutionLocalYParam3());
  const double param4 = (isLocalXAxis ? resolutionLocalXParam4() : resolutionLocalYParam4());
  const double param5 = (isLocalXAxis ? resolutionLocalXParam5() : resolutionLocalYParam5());
  const double param6 = (isLocalXAxis ? resolutionLocalXParam6() : resolutionLocalYParam6());
  const double param7 = (isLocalXAxis ? resolutionLocalXParam7() : resolutionLocalYParam7());
  const double param8 = (isLocalXAxis ? resolutionLocalXParam8() : resolutionLocalYParam8());
  const double param9 = (isLocalXAxis ? resolutionLocalXParam9() : resolutionLocalYParam9());

  // |tan(deepAngle)|
  const double x = fabsTanDeepAngle;

  // Use of parametrization model
  const double resolutionLocalAxis = param0 + param1 * x 
    + param2 * exp(-param9 * x) * cos(param3 * x + param4)
    + param5 * exp(-0.5 * pow(((x - param6) / param7), 2))
    + param8 * pow(x, 0.5);

  return resolutionLocalAxis;
}


/*
 * Chekc whether the value of at least one parameter to compute a parametrized spatial resolution in X, is specified.
 */
const bool DetectorModule::hasAnyResolutionLocalXParam() const { 
  return ( resolutionLocalXParam0.state() 
	   || resolutionLocalXParam1.state() 
	   || resolutionLocalXParam2.state() 
	   || resolutionLocalXParam3.state() 
	   || resolutionLocalXParam4.state()
	   || resolutionLocalXParam5.state() 
	   || resolutionLocalXParam6.state() 
	   || resolutionLocalXParam7.state()
	   || resolutionLocalXParam8.state() 
	   || resolutionLocalXParam9.state()
	   );
}


/*
 * Chekc whether the value of at least one parameter to compute a parametrized spatial resolution in Y, is specified.
 */
const bool DetectorModule::hasAnyResolutionLocalYParam() const { 
  return ( resolutionLocalYParam0.state() 
	   || resolutionLocalYParam1.state() 
	   || resolutionLocalYParam2.state() 
	   || resolutionLocalYParam3.state() 
	   || resolutionLocalYParam4.state()
	   || resolutionLocalYParam5.state() 
	   || resolutionLocalYParam6.state() 
	   || resolutionLocalYParam7.state()
	   || resolutionLocalYParam8.state() 
	   || resolutionLocalYParam9.state()
	   );
}


/*
 * Compute the alpha incident angle.
 * See README for definition of alpha angle. 
 * alpha = (X, projectedTrack) (oriented angle between the 2 vectors).
 * X is the vector of the Lorentz drift (= localX).
 * projectedTrack is the projection of the track in the plane of normal localY.
 */
const double DetectorModule::alpha(const TVector3& trackDirection) const {

  // Sensor local X vector (= Lorentz drift in TBPX).
  const TVector3& localX = getLocalX();

  // Sensor local Y vector.
  const TVector3& localY = getLocalY();

  // Project the track vector into the plane of normal localY.
  const TVector3& projectedTrack = CoordinateOperations::projectv1OnPlaneOfNormalUnitv2(trackDirection, localY);

  // alpha = (localX, projectedTrack) (oriented angle between the 2 vectors).
  const double alpha = femod(localX.Angle(projectedTrack), 2. * M_PI);

  if (alpha >= M_PI) { 
    logERROR("alpha angle should be in ]0 180[, but found alpha = " + any2str(alpha * 180. / M_PI)); 
  }

  /* Keep old formula.
     This old formula is only true in (CMS_X, CMS_Y) plane.
     It is is a projection of the new formula onto the (CMS_X, CMS_Y) plane (perfect match).

     const double trackPhi = trackDirection.Phi();

     // Oriented angle: between the normal to the sensor (oriented towards outer radius) and the track vector.
     const double sensorNormalToTrackDeltaPhi = trackPhi - (center().Phi() + skewAngle());

     // Oriented angle: between local X vector (= Lorentz drift) and the track vector.
     const double sensorXToTrackDeltaPhi = (!flipped() ? M_PI / 2. + sensorNormalToTrackDeltaPhi : M_PI / 2. - sensorNormalToTrackDeltaPhi);

     const double alpha = femod(sensorXToTrackDeltaPhi, 2.*M_PI);
  */

  return alpha;
}


/*
 * Compute the beta incident angle.
 * See README for definition of beta angle.
 * alpha = (Y, projectedTrack) (oriented angle between the 2 vectors).
 * Y is the sensor localY vector.
 * projectedTrack is the projection of the track in the plane of normal localX.
 */
const double DetectorModule::beta(const TVector3& trackDirection) const {

  // Sensor local X vector
  const TVector3& localX = getLocalX();

  // Sensor local Y vector.
  const TVector3& localY = getLocalY();

  // Project the track vector into the plane of normal localX.
  const TVector3& projectedTrack = CoordinateOperations::projectv1OnPlaneOfNormalUnitv2(trackDirection, localX);

  // beta = (localY, projectedTrack) (oriented angle between the 2 vectors).
  const double beta = femod(localY.Angle(projectedTrack), 2. * M_PI);

  if (beta >= M_PI) { 
    logERROR("beta angle should be in ]0 180[, but found beta = " + any2str(beta * 180. / M_PI)); 
  }

  /* Keep old formula.
     This old formula is not really correct.
     It does not take into account a dependency in track vector's phi (true story!!), which occurs for both the barrel and the forward.
     This is because we are interested in the PROJECTION of the track vector into the plane of normal localX.
     When we project the track vector into the plane of normal localX:
     In the event the hit is not in the (sensor center, CMS_Z) plane, or if the module is skewed: 
     * the component in (CMS_X, CMS_Y) plane of the track vector gets a smaller norm. 
     * the component along CMS_Z of the track vector stays the same.
     Hence, the projected track vector is more 'parallel to CMS_Z' than the track vector.
     This implies that beta angle is further from M_PI/2 than eta is the barrel, and closer to M_PI/2 than eta in the forward.
     The effect is quite negligible though, except for skewed ladders. 
     It is way cleaner to have a vectorial formula amyway.

     const double beta = theta + tiltAngle();

     // Oriented angle: between local Y vector and the track vector.
     const double orientedBeta = (fabs(tiltAngle() - M_PI/2.) < insur::geom_zero ? beta : (M_PI - beta));
  */

  return beta;
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
  if (numSensors() == 1 || (numSensors() > 1 && sensorLayout() == SensorLayout::MONO)) {
    // Assuming 1 or more non-overlapping sensors on the same plane
    for (auto& s : sensors_) {
      auto segm = s.checkHitSegment(trackOrig, trackDir);
      // <SMe>The following line used to return HitType::BOTH. Changing to INNER in order to avoid double hit counting</SMe>
      if (segm.second > -1) { gc = segm.first; ht = HitType::INNER; break; }
    }
  } else {
    auto inSegm = innerSensor().checkHitSegment(trackOrig, trackDir);
    auto outSegm = outerSensor().checkHitSegment(trackOrig, trackDir);
    if (inSegm.second > -1 && outSegm.second > -1) { 
      gc = inSegm.first; // in case of both sensors are hit, the inner sensor hit coordinate is returned
      ht = ((zCorrelation() == SAMESEGMENT && (inSegm.second / (maxSegments()/minSegments()) == outSegm.second)) || zCorrelation() == MULTISEGMENT) ? HitType::STUB : HitType::BOTH;
    } else if (inSegm.second > -1) { gc = inSegm.first; ht = HitType::INNER; }
    else if (outSegm.second > -1) { gc = outSegm.first; ht = HitType::OUTER; }
  }
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


// OT CABLING
const int DetectorModule::isPositiveCablingSide() const {
  int isPositiveCablingSide = 0;
  const OuterBundle* myBundle = getBundle();
  if (myBundle) {
    isPositiveCablingSide = (myBundle->isPositiveCablingSide() ? 1 : -1);
  }
  return isPositiveCablingSide;
}


const int DetectorModule::bundlePlotColor() const {
  int bundlePlotColor = 0;
  const OuterBundle* myBundle = getBundle();
  if (myBundle) {
    bundlePlotColor = myBundle->plotColor();
  }
  return bundlePlotColor;
}


const int DetectorModule::opticalChannelSectionPlotColor() const {
  int opticalChannelPlotColor = 0;
  const OuterBundle* myBundle = getBundle();
  if (myBundle) {
    const ChannelSection* mySection = myBundle->opticalChannelSection();
    if (mySection) {
      opticalChannelPlotColor = mySection->plotColor();
    }
  }
  return opticalChannelPlotColor;
}


const int DetectorModule::powerChannelSectionPlotColor() const {
  int powerChannelPlotColor = 0;
  const OuterBundle* myBundle = getBundle();
  if (myBundle) {
    const ChannelSection* mySection = myBundle->powerChannelSection();
    if (mySection) {
      powerChannelPlotColor = mySection->plotColor();
    }
  }
  return powerChannelPlotColor;
}


const OuterDTC* DetectorModule::getDTC() const {
  const OuterDTC* myDTC = nullptr;
  const OuterBundle* myBundle = getBundle();
  if (myBundle) {
    const OuterCable* myCable = myBundle->getCable();
    if (myCable) {
      myDTC = myCable->getDTC();
    }
  }
  return myDTC;
}


const int DetectorModule::dtcPlotColor() const {
  int dtcPlotColor = 0;
  const OuterDTC* myDTC = getDTC();
  if (myDTC) {
    dtcPlotColor = myDTC->plotColor();
  }
  return dtcPlotColor;
}


const int DetectorModule::dtcPhiSectorRef() const {
  int dtcPhiSectorRef = 0;
  const OuterDTC* myDTC = getDTC();
  if (myDTC) {
    dtcPhiSectorRef = myDTC->phiSectorRef();
  }
  return dtcPhiSectorRef;
}


// IT CABLING
const int DetectorModule::isPositiveZEnd() const {
  int isPositiveZEnd = 0;
  const PowerChain* myPowerChain = getPowerChain();
  if (myPowerChain) {
    isPositiveZEnd = (myPowerChain->isPositiveZEnd() ? 1 : -1);
  }
  return isPositiveZEnd;
}


const bool DetectorModule::isPositiveXSide() const {
  bool isPositiveXSide;
  const PowerChain* myPowerChain = getPowerChain();
  if (myPowerChain) {
    isPositiveXSide = myPowerChain->isPositiveXSide();
  } else {
    logERROR("Module found not connected to any power chain. I do not know whether it is positive x side or not");
    return true;
  }
  return isPositiveXSide;
}


const int DetectorModule::powerChainPlotColor() const {
  int powerChainPlotColor = 0;
  const PowerChain* myPowerChain = getPowerChain();
  if (myPowerChain) {
    powerChainPlotColor = myPowerChain->plotColor();
  }
  return powerChainPlotColor;
}


const int DetectorModule::gbtPlotColor() const {
  int gbtPlotColor = 0;
  const GBT* myGBT = getGBT();
  if (myGBT) {
    gbtPlotColor = myGBT->plotPowerChainColor();
  }
  return gbtPlotColor;
}


const InnerBundle* DetectorModule::getInnerBundle() const {
  const InnerBundle* myBundle = nullptr;
  const GBT* myGBT = getGBT();
  if (myGBT) {
    myBundle = myGBT->getBundle();
  }
  return myBundle;
}


const int DetectorModule::innerBundlePlotColor() const {
  int bundlePlotColor = 0;
  const InnerBundle* myBundle = getInnerBundle();
  if (myBundle) {
    bundlePlotColor = myBundle->plotColor();
  }
  return bundlePlotColor;
}


const InnerDTC* DetectorModule::getInnerDTC() const {
  const InnerDTC* myDTC = nullptr;
  const InnerBundle* myBundle = getInnerBundle();
  if (myBundle) {
    myDTC = myBundle->getDTC();
  }
  return myDTC;
}


const int DetectorModule::innerDTCPlotColor() const {
  int dtcPlotColor = 0;
  const InnerDTC* myDTC = getInnerDTC();
  if (myDTC) {
    dtcPlotColor = myDTC->plotColor();
  }
  return dtcPlotColor;
}



void BarrelModule::build() {
  try {
    DetectorModule::build();

    rAxis_ = XYZVector(1., 0., 0.);

    // tilt
    tiltAngle_ = 0.;

    for (auto& s : sensors_) { s.subdet(ModuleSubdetector::BARREL); }
  }
  catch (PathfulException& pe) { pe.pushPath(*this, myid()); throw; }
  cleanup();
  builtok(true);
}


void EndcapModule::build() {
  try {
    DetectorModule::build();

    rAxis_ = XYZVector(0., -1., 0.);

    // tilt
    tiltAngle_ = M_PI_2;

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


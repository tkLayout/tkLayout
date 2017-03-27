#include "GeometricModule.hh"

double ModuleHelpers::polygonAperture(const Polygon3d<4>& poly) { 
  auto minmax = std::minmax_element(poly.begin(), poly.end(), [](const XYZVector& v1, const XYZVector& v2) { return v1.Phi() < v2.Phi(); }); 
  return minmax.second->Phi() - minmax.first->Phi(); 
}


// Distance between origin and module if hits
// Otherwise -1
double GeometricModule::trackCross(const XYZVector& PL, // Base line point
                                   const XYZVector& PU) // Line direction
{
  double distance;

  distance = triangleCross(basePoly().getVertex(0),
                           basePoly().getVertex(1),
                           basePoly().getVertex(2),
                           PL, PU);

  if (distance<0) {
    distance = triangleCross(basePoly().getVertex(0),
                             basePoly().getVertex(2),
                             basePoly().getVertex(3),
                             PL, PU);
  }

  if (distance>=0) {
    numHits_++;
    distance /= PU.r();
  } else {
    return -1;
  }

  return distance;

}



double GeometricModule::triangleCross(const XYZVector& P1, // Triangle points
                                      const XYZVector& P2,
                                      const XYZVector& P3,
                                      const XYZVector& PL, // Base line point
                                      const XYZVector& PU) // Line direction
{
  bool moduleHit = false;

  // Triangle coordinates
  // t - P1 = alpha * (P2-P1) + beta * (P3-P1)

  // Line coordinates:
  // r = PL + gamma * PU


  // How to solve the generic problem:

  // d = PL - P1
  // A = (P2-P1, P3-P1, -PU)
  // v = intersection point
  TVectorD d = TVectorD(3);
  TMatrixD A = TMatrixD(3, 3);
  TVectorD v = TVectorD(3);

  d(0)=PL.x()-P1.x();
  d(1)=PL.y()-P1.y();
  d(2)=PL.z()-P1.z();

  A(0, 0)=P2.x()-P1.x();
  A(0, 1)=P3.x()-P1.x();
  A(0, 2)=-1*PU.x();
  A(1, 0)=P2.y()-P1.y();
  A(1, 1)=P3.y()-P1.y();
  A(1, 2)=-1*PU.y();
  A(2, 0)=P2.z()-P1.z();
  A(2, 1)=P3.z()-P1.z();
  A(2, 2)=-1*PU.z();

  Double_t determ;
  A.InvertFast(&determ);

  // The matrix is invertible
  if (determ!=0) {
    // v = A^{-1} * d
    v = A * d;
    // v(0) = alpha : triangle local coordinate
    // v(1) = beta  : triangle local coordinate
    // v(2) = gamma : line coordinate

    //     std::cout << "Alpha: " << v(0) << std::endl;
    //     std::cout << "Beta:  " << v(1) << std::endl;
    //     std::cout << "Gamma: " << v(2) << std::endl;

    if (
      (v(0)>=0)&&
      (v(1)>=0)&&
      ((v(0)+v(1))<=1)
      ) {
      moduleHit = true;
    } else {
      moduleHit = false;
    }

  } else {
    // It does not cross the triangle
    moduleHit=false;
    std::cout << "Matrix is not invertible" << std::endl; // debug
  }


  if (moduleHit) {
    return v(2);
  } else {
    return -1.;
  }

}






void RectangularModule::check() {
  GeometricModule::check();

  if (length_.state() && width.state()) {
    aspectRatio(length_()/width());
  } else if (!length_.state() && width.state() && aspectRatio.state()) { 
    length_(width() * aspectRatio());
  } else if (length_.state() && !width.state() && aspectRatio.state()) {
    width(length_() / aspectRatio());
  } else if (!length_.state() && !width.state() && aspectRatio.state()) {
    length_(waferDiameter() * sin(atan(aspectRatio())));
    width(waferDiameter() * cos(atan(aspectRatio())));
  } else {
    throw PathfulException("Module geometry is inconsistently specified");
  }
}


void RectangularModule::build() {
  try { 
    check();
    
    int iContourPoint=0;
    ContourPoint p;
    while (contourPointNode.count(iContourPoint) > 0) {
      p.store(contourPointNode.at(iContourPoint));
      ROOT::Math::XYZVector contourPoint;
      contourPoint.SetXYZ(p.pointX(), p.pointY(), 0);
      contour_.push_back(contourPoint);
      iContourPoint++;
    }

    float l = length(), w = width();
    basePoly_ << XYZVector( l/2, w/2, 0)
              << XYZVector(-l/2, w/2, 0)
              << XYZVector(-l/2,-w/2, 0)
              << XYZVector( l/2,-w/2, 0);
    cleanup();
    builtok(true);
  }
  catch (PathfulException& pe) { pe.pushPath(*this, myid()); throw; }
}




void WedgeModule::build() {
  try {
    check();

    //////// BEGIN COPY-PASTE WITH MINIMAL ADJUSTMENTS ////////
    cropped_ = false;

    double r = waferDiameter()/2.;
    double phi = buildAperture()/2.;// We need the half angle covered by the module
    double d = buildDistance();

    //double gamma1; // alternate
    double gamma2;
    double h1, h2, b1, b2, dfar, l;
    h1 = d * tan(phi);// The short (half)base
    // Distance of short base from the wafer center
    b1 = sqrt(pow(r, 2)-pow(h1, 2)); // main + alternate

    // y coordinate of the wafer center
    l = b1 + d; // main

    // Distance of the far angle form the z axis
    gamma2 = l*cos(phi) + sqrt(pow(r, 2)-pow(l*sin(phi), 2)); // main

    h2 = gamma2 * sin(phi);// The long (half)base
    dfar = gamma2 * cos(phi);// Distance of long base from the z axis

    // The distance of the long base from the wafer center
    //b2 =  sqrt(pow(r,2)-pow(h2,2)); // old
    b2 = dfar - d - b1;

    // NOTE: in principle we don't need to compute b2 to get the
    // module's corner coordinates. Still we use this way of computing
    // b2 to ease the computation of the module's area

    // Add a check: if the module overcomes the max rho
    // it must be cut.
    if (buildCropDistance.state() && dfar > buildCropDistance()) {
      amountCropped_ = dfar - buildCropDistance();
      b1 = 0;
      b2 = buildCropDistance() - d;
      h2 = h1/d * buildCropDistance();
      cropped_ = true;
    }

    // Some member variable computing:
    area_     = fabs((b1+b2) * (h2+h1));
    length_   = (b1 + b2);
    //phiWidth_ = 2*phi;
    minWidth_  = 2 * h1;
    maxWidth_  = 2 * h2;
  //dist_     = d;
  //aspectRatio_ = length_/(h1+h2);

    // Right-handed drawing, (not left-handed as previously)
    basePoly_ << (XYZVector( length_/2., maxWidth_/2., 0))
              << (XYZVector(-length_/2., minWidth_/2., 0))
              << (XYZVector(-length_/2.,-minWidth_/2., 0))
              << (XYZVector( length_/2.,-maxWidth_/2., 0));

    cleanup();
    builtok(true);
  }
  catch (PathfulException& pe) { pe.pushPath(*this, myid()); throw; }
}

define_enum_strings(ModuleShape) = { "rectangular", "wedge" };

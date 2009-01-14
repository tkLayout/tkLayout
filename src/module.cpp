#include <iostream>
#include <cmath>
#include "module.hh"

#include "Math/RotationZ.h"
#include "Math/Vector3D.h"

#include "TMatrixD.h"
#include "TVectorD.h"

#include "TGeoManager.h"
#include "TGeoArb8.h"
#include "TGeoMatrix.h"

#include <boost/lexical_cast.hpp>

using namespace ROOT::Math;

/******************/
/*                */
/* Generic module */
/*                */
/******************/

Module::~Module() {

}

Module::Module() {
  waferDiameter_= defaultWaferDiameter_;
  setDefaultParameters();
}

Module::Module(double waferDiameter) {
  waferDiameter_ = waferDiameter;
  setDefaultParameters();
}

void Module::setDefaultParameters() {
  // See the header file for the values of these
  height_             = defaultHeight_;
  nHits_              = defaultNHits_;
  thickness_          = defaultThickness_;
  area_               = defaultArea_;
  id_                 = "NoId";
  tag_                = "";
  type_               = "NoType";
  thisVolume_         = NULL;
  color_              = defaultColor_;
  inSection_          = defaultInSection_;
  nChannelsPerFace_   = defaultChannelsPerFace_;
  nSegments_          = defaultSegments_;
  nStripAcross_       = defaultStripAcross_;
  nFaces_             = defaultFaces_;
  readoutType_        = Strip;
  computedBoundaries_ = false;
}

void Module::print() {
  std::cout <<
    "(" <<
    corner_[3].X() << ", " << 
    corner_[3].Y() << ", " << 
    corner_[3].Z() << ")\t\t" << 
    "(" <<
    corner_[2].X() << ", " << 
    corner_[2].Y() << ", " << 
    corner_[2].Z() << ")" << std::endl <<
    "(" <<
    corner_[0].X() << ", " << 
    corner_[0].Y() << ", " << 
    corner_[0].Z() << ")\t\t" << 
    "(" <<
    corner_[1].X() << ", " << 
    corner_[1].Y() << ", " << 
    corner_[1].Z() << ")" << std::endl ;
}

void Module::translate(XYZVector Delta) {
  for (int i=0; i<4; i++) {
//     corner_[i].x+=Delta.x;
//     corner_[i].y+=Delta.y;
//     corner_[i].z+=Delta.z;
    corner_[i]+=Delta;
  }
}

void Module::rotatePhi(double phi) {
  RotationZ phiRot(phi);

  for (int i=0; i<4; i++) {
    corner_[i]=phiRot(corner_[i]);
    //    newx=corner_[i].x*cos(phi)-corner_[i].y*sin(phi);
    //    newy=corner_[i].x*sin(phi)+corner_[i].y*cos(phi);
    //    corner_[i].x=newx;
    //    corner_[i].y=newy;
  }
}

void Module::rotateY_PI() {
  for (int i=0; i<4; i++) {
    corner_[i].SetX(corner_[i].X()*-1.);
    corner_[i].SetZ(corner_[i].Z()*-1.);
  }
}

void Module::reflectZ() {
  for (int i=0; i<4; i++) {
    corner_[i].SetZ(corner_[i].Z()*-1.);
  }
}


XYZVector* Module::marginBorderSide(double margin, int side){
  
  XYZVector* result;

  if ((side<0)||(side>3)) {
    std::cerr << "Module::marginBorderSide was invoked with side out of range" << std::endl;
    return NULL;
  };

  int otherSideIndex = (side+2) % 4;
  int nextPoint=(side+1) % 4;
  int otherNextPoint=(otherSideIndex+1) % 4;

  XYZVector thisSide=(corner_[side]+corner_[nextPoint])/2.;
  XYZVector otherSide=(corner_[otherSideIndex]+corner_[otherNextPoint])/2.;

  XYZVector delta=otherSide-thisSide;
  delta *= margin/delta.R();
  
  result = new XYZVector(thisSide + delta);

  return result;
}


XYZVector* Module::marginBorder(double widthMargin, double lengthMargin, int corner ) {

  XYZVector* result;

  if ((corner<0)||(corner>3)) {
    std::cerr << "Module::marginBorder was invoked with corner out of range" << std::endl;
    return NULL;
  };

  
  XYZVector *thisCorner, *neighbourWidthCorner, *neighbourLengthCorner;
  XYZVector deltaWidth, deltaLength;

  int parity=corner%2;
  int neighbourWidthIndex, neighbourLengthIndex;

  if (parity==0) {
    neighbourWidthIndex=corner+1;
    neighbourLengthIndex=corner-1;
  } else {
    neighbourWidthIndex=corner-1;
    neighbourLengthIndex=corner+1;
  }

  if (neighbourWidthIndex<0)neighbourWidthIndex+=4;
  if (neighbourLengthIndex<0)neighbourLengthIndex+=4;
  if (neighbourWidthIndex>3)neighbourWidthIndex-=4;
  if (neighbourLengthIndex>3)neighbourLengthIndex-=4;

  neighbourWidthCorner=&corner_[neighbourWidthIndex];
  neighbourLengthCorner=&corner_[neighbourLengthIndex];
  thisCorner=&corner_[corner];

  deltaWidth=(*neighbourWidthCorner)-(*thisCorner);
  deltaLength=(*neighbourLengthCorner)-(*thisCorner);

  deltaWidth*=(widthMargin/deltaWidth.R());
  deltaLength*=(lengthMargin/deltaLength.R());

  result=new XYZVector((*thisCorner)+deltaWidth+deltaLength);

  return result;
}


int Module::projectSideZ(int side, double destZ, double displaceZ /*=0*/ ) {
  if ((side<0)||(side>3)) {
    std::cerr << "Module::projectSideZ was invoked with a nonsense side index" << std::endl;
    return -1;
  }

  XYZVector* displacement = NULL;

  if (displaceZ!=0) {
    displacement = new XYZVector(0,0,displaceZ);
    for (int i=0; i<4; i++) {
      corner_[i]-=(*displacement);
    }
    destZ-=displacement->Z();
  }


  int firstPoint=side;
  int secondPoint=(side+1) % 4;

//   std::cerr << "Refpoints:" << firstPoint << " " << secondPoint << std::endl;

  XYZVector refPoint = (corner_[firstPoint]+corner_[secondPoint])/2.;

//   std::cerr << "refPoint"
// 	    <<  refPoint.X() << " " 
// 	    <<  refPoint.Y() << " " 
// 	    <<  refPoint.Z()
// 	    << std::endl;
  
  double ratio = destZ/refPoint.Z();

//   std::cerr << "ratio:" << ratio << std::endl;

  XYZVector Delta = refPoint*(ratio-1);

//   std::cerr << "Delta:"
// 	    <<  Delta.X() << " " 
// 	    <<  Delta.Y() << " " 
// 	    <<  Delta.Z()
// 	    << std::endl;
  
  for (int i=0; i<4; i++) {
    corner_[i]+=Delta;
  }

  if (displacement!=NULL) {
    for (int i=0; i<4; i++) {
      corner_[i]+=(*displacement);
    }
    delete displacement;
  }
  
  return 0;
}

int Module::projectSideRho(int side, double destRho, double displaceZ /*=0*/ ) {
  if ((side<0)||(side>3)) {
    std::cerr << "Module::placeSideRho was invoked with a nonsense side index" << std::endl;
    return -1;
  }

  XYZVector* displacement = NULL;

  if (displaceZ!=0) {
    displacement = new XYZVector(0,0,displaceZ);
    for (int i=0; i<4; i++) {
      corner_[i]-=(*displacement);
    }
    destRho-=displacement->Rho();
  }

  int firstPoint=side;
  int secondPoint=(side+1) % 4;

//   std::cerr << "Refpoints: " << firstPoint << "," << secondPoint << std::endl;
//   std::cerr << firstPoint << "= "
// 	    << corner_[firstPoint].X() << " " 
// 	    << corner_[firstPoint].Y() << " " 
// 	    << corner_[firstPoint].Z()
// 	    << std::endl;  
//   std::cerr << secondPoint << "= "
// 	    << corner_[secondPoint].X() << " " 
// 	    << corner_[secondPoint].Y() << " " 
// 	    << corner_[secondPoint].Z()
// 	    << std::endl;
  

  XYZVector refPoint = (corner_[firstPoint]+corner_[secondPoint])/2.;

//   std::cerr << "refPoint= "
// 	    <<  refPoint.X() << " " 
// 	    <<  refPoint.Y() << " " 
// 	    <<  refPoint.Z()
// 	    << std::endl;
  
  
  double ratio = destRho/refPoint.Rho();

//   std::cerr << "ratio:" << ratio << std::endl;

  XYZVector Delta = refPoint*(ratio-1);

//   std::cerr << "Delta:"
// 	    <<  Delta.X() << " " 
// 	    <<  Delta.Y() << " " 
// 	    <<  Delta.Z()
// 	    << std::endl;
  
  for (int i=0; i<4; i++) {
    corner_[i]+=Delta;
  }

  if (displacement!=NULL) {
    for (int i=0; i<4; i++) {
      corner_[i]+=(*displacement);
    }
    delete displacement;
  }
  
  return 0;
}


// Distance between origin and module if hits
// Otherwise -1
double Module::trackCross(const XYZVector& PL, // Base line point
			  const XYZVector& PU) // Line direction
{
  double distance;

  distance = triangleCross(corner_[0],
			   corner_[1],
			   corner_[2], 
			   PL, PU);

  if (distance<0) {
    distance = triangleCross(corner_[0],
			     corner_[2],
			     corner_[3], 
			     PL, PU);
  }
  
  if (distance>=0) {
    nHits_++;
    distance /= PU.r();
  } else {
    return -1;
  }

  return distance;
  
}



double Module::triangleCross(const XYZVector& P1, // Triangle points
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
  TMatrixD A = TMatrixD(3,3);
  TVectorD v = TVectorD(3);

  d(0)=PL.x()-P1.x();
  d(1)=PL.y()-P1.y();
  d(2)=PL.z()-P1.z();

  A(0,0)=P2.x()-P1.x();
  A(0,1)=P3.x()-P1.x();
  A(0,2)=-1*PU.x();
  A(1,0)=P2.y()-P1.y();
  A(1,1)=P3.y()-P1.y();
  A(1,2)=-1*PU.y();
  A(2,0)=P2.z()-P1.z();
  A(2,1)=P3.z()-P1.z();
  A(2,2)=-1*PU.z();

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

// Creates a new object in the volume container
// of geometry geom
void Module::shapeVolume(TGeoVolume* container,
			 TGeoMedium* medium,
			 TGeoManager* geom) {


  // Fist of all: find the local coordinates of the module
  XYZVector ex, ey, ez;
  XYZVector b, c, d, p;
  
  b = corner_[1]-corner_[0];
  c = corner_[2]-corner_[0];
  d = corner_[3]-corner_[0];
  
  ex = b/b.R();
  p  = (d.Dot(ex) * ex);
  ey = d - p;
  ey /= ey.R();
  ez = ex.Cross(ey);

  // Local coordinates of module's corners
  // x'           y'
  // 0            0
  // b.R()        0
  // c.Dot(ex)    c.Dot(ey)
  // d.Dot(ex)    d.Dot(ey)

  // TODO: we should check that {b, c, d}.Dot(ez) < 10E-5 or so
//   std::cerr << b.Dot(ez) << "==0 ?" << std::endl;
//   std::cerr << c.Dot(ez) << "==0 ?" << std::endl;
//   std::cerr << d.Dot(ez) << "==0 ?" << std::endl;
    
  // This is the rotation (columns of matrix myRot are
  // the module axes expressed in the global frame)
  Double_t matrixElements[9];
  matrixElements[0]=ex.X();
  matrixElements[1]=ey.X();
  matrixElements[2]=ez.X();
  matrixElements[3]=ex.Y();
  matrixElements[4]=ey.Y();
  matrixElements[5]=ez.Y();
  matrixElements[6]=ex.Z();
  matrixElements[7]=ey.Z();
  matrixElements[8]=ez.Z();


  TGeoRotation *myRot    = new TGeoRotation();
  myRot->SetMatrix(matrixElements); // Ecco la matrice

  // We add the translation (point 0 of module)
  TGeoCombiTrans *myCombi = new TGeoCombiTrans(corner_[0].X(),
					       corner_[0].Y(),
					       corner_[0].Z(),
					       myRot);

  // Delete the volume if already existing
  if (thisVolume_) delete thisVolume_;

  // Create and get a grip on th object
  // TODO: place a global counter here...
  TGeoVolume *thisVolume_ = geom->MakeArb8(id_.c_str(), medium, thickness_);
  thisVolume_->SetLineColor(color_);
  TGeoArb8 *arb = (TGeoArb8*)thisVolume_->GetShape();

  // Add the points of the module in its local coords
  // with respect to point 0
  for (int i=0; i<5; i+=4) {
    arb->SetVertex(0+i, 0        , 0);
    arb->SetVertex(1+i, b.R()    , 0);
    arb->SetVertex(2+i, c.Dot(ex), c.Dot(ey));
    arb->SetVertex(3+i, d.Dot(ex), d.Dot(ey));
  }
  
  // And finally place it into the empty space!
  container->AddNode(thisVolume_, 1, myCombi);

}


TPolyLine3D* Module::getContour() {
  TPolyLine3D *contour = new TPolyLine3D(5);
  for (Int_t i=0;i<4;i++) {
    contour->SetPoint(i,
		      corner_[i].X(),
		      corner_[i].Y(),
		      corner_[i].Z());
  }
  contour->SetPoint(4,
		    corner_[0].X(),
		    corner_[0].Y(),
		    corner_[0].Z());
  contour->SetLineColor(color_);
  return contour;
}





edge Module::getEdgeRhoSide(int direction) {
  Module fakeModule;
  XYZVector* marginBorder;
  
  int nextIndex;
  for (int i=0; i<4; i++) {
    nextIndex=(i+1) % 4;
    marginBorder = marginBorderSide(0, i);
    fakeModule.corner_[i]=(*marginBorder);
    delete marginBorder;
  }

  return fakeModule.getEdgeRho(direction);
}


edge Module::getEdgeRho(int direction) {
  edge result;
  
  if (direction==0) {
    std::cerr << "Module::getEdgeRho direction was chosen as zero." << std::endl;
    result.first=-1;
    result.second=-1;
    return result;
  }

  XYZVector* referencePoint[4];
  
  for (int i=0; i<4; i++) {
    referencePoint[i]=&corner_[i];      
  }


  if (direction<0) {
    // Get the point with lowest Rho
    result.first=referencePoint[0]->Rho();
    result.second=0;
    for (int i=1; i<4; i++) {
      if (referencePoint[i]->Rho()<result.first) {
	result.first=referencePoint[i]->Rho();
	result.second=i;
      }
    }
  } else {
    // Get the point with highest Rho
    result.first=referencePoint[0]->Rho();
    result.second=0;
    for (int i=1; i<4; i++) {
      if (referencePoint[i]->Rho()>result.first) {
	result.first=referencePoint[i]->Rho();
	result.second=i;
      }
    }
  }


  return result;
  
}

int Module::setEdgeRho(double newRho, int direction) {
  if (direction==0) {
    std::cerr << "Module::setEdgeRho direction was chosen as zero." << std::endl;
    return -1;
  }

  edge edgePoint = getEdgeRho((-1)*direction);

  XYZVector start = corner_[edgePoint.second];
  XYZVector dest  = start/start.Rho()*newRho;
  XYZVector delta = dest-start;  
  
  for (int i=0; i<4; i++) {
    corner_[i]+=delta;
  }

  return 1;
}

int Module::setEdgeRhoSide(double newRho, int direction) {
  if (direction==0) {
    std::cerr << "Module::setEdgeRho direction was chosen as zero." << std::endl;
    return -1;
  }

  edge edgePoint = getEdgeRhoSide((-1)*direction);
  XYZVector* anEdge = marginBorderSide(0,edgePoint.second);

  XYZVector start = *(anEdge);
  XYZVector dest  = start/start.Rho()*newRho;
  XYZVector delta = dest-start;  
  
  for (int i=0; i<4; i++) {
    corner_[i]+=delta;
  }

  delete anEdge;

  return 1;
}

double Module::getMinTheta() {
  double result, aTheta;
  result = corner_[0].Theta();
  for (int i=1; i<4; i++) {
    aTheta = corner_[i].Theta();
    if (aTheta<result) result=aTheta;
  }
  return result;
}

double Module::getMaxTheta() {
  double result, aTheta;
  result = corner_[0].Theta();
  for (int i=1; i<4; i++) {
    aTheta = corner_[i].Theta();
    if (aTheta>result) result=aTheta;
  }
  return result;
}

double Module::getMeanTheta() {

  XYZVector meanPoint(0,0,0);

  for (int i=0; i<4; i++) {
    meanPoint += corner_[i];
  }
  meanPoint /= 4;

  return meanPoint.Theta();
}

double Module::getMaxRho() {
  double maxRho=-1;
  for (uint i = 0; i < 4 ; i++) {
    if (corner_[i].Rho()>maxRho) maxRho=corner_[i].Rho();
  }
  return maxRho;
}

double Module::getMinRho() {
  double minRho=corner_[0].Rho();
  for (uint i = 1; i < 4 ; i++) {
    if (corner_[i].Rho()<minRho) minRho=corner_[i].Rho();
  }
  return minRho;
}

double Module::getMaxZ() {
  double maxZ=corner_[0].Z();
  for (uint i = 1; i < 4 ; i++) {
    if (corner_[i].Z()>maxZ) maxZ=corner_[i].Z();
  }
  return maxZ;
}

double Module::getMinZ()  {
  double minZ=corner_[0].Z();
  for (uint i = 1; i < 4 ; i++) {
    if (corner_[i].Z()<minZ) minZ=corner_[i].Z();
  }
  return minZ;
}


XYZVector Module::getMeanPoint() {

  XYZVector meanPoint(0,0,0);

  for (int i=0; i<4; i++) {
    meanPoint += corner_[i];
  }
  meanPoint /= 4;

  return meanPoint;
}

void Module::computeStripArea() {
  XYZVector meanPoint = getMeanPoint();
  // double meanEta = meanPoint.eta(); // TODO: remove this line
  double meanPhi = meanPoint.phi();

  double maxTheta, minTheta, maxPhi, minPhi, aPhi;
  maxTheta = getMaxTheta();
  minTheta = getMinTheta();
  //  std::cerr << "theta:\t" << minTheta << "\t" << maxTheta << std::endl;
  
  Module* anotherModule = new Module(*this);
  anotherModule->rotatePhi(-1*meanPhi);
  maxPhi=0;
  minPhi=0;
  for (int i=0; i<4 ; i++) {
    aPhi=anotherModule->corner_[i].phi();
    if (aPhi>maxPhi) maxPhi=aPhi;
    if (aPhi<minPhi) minPhi=aPhi;
  }

  // std::cerr << "phi:\t" << minPhi << "\t" << maxPhi << std::endl;
  
  double phiWidth, etaWidth;
  phiWidth = maxPhi-minPhi;
  etaWidth = log(tan(maxTheta/2.))-1*log(tan(minTheta/2.));
  //  std::cerr << "widths: theta:\t" << etaWidth << " phi:\t" << phiWidth << std::endl;

  stripArea_ = phiWidth * etaWidth / double(nChannelsPerFace_);
  //  std::cerr << "area:\t" << stripArea_ << std::endl;
}

double Module::getOccupancyPerEvent() {
  XYZVector meanPoint = getMeanPoint();
  double meanEta = fabs(meanPoint.eta());

  computeStripArea();
  double spOcc;

  if (meanEta<1) {
    spOcc = 3.7 * (1+meanEta);
  } else {
    spOcc = 3.7 * 2;
  }

  // Per ottenere l'occupanza vera basta moltiplicare spOcc
  // per l'area della strip espressa in unità di (phi, eta).
  // Il risultato deve essere in seguito moltiplicato per il numero di eventi di minimum bias per evento.
  // (5, 24 o 400, per bassa, alta e super luminosità). 

  //  std::cerr << "occupancy: " << spOcc*stripArea_ << std::endl;

  return spOcc*stripArea_;
}


double Module::getLowPitch() {
  XYZVector acrossV = corner_[0] - corner_[3];
  return (acrossV.R()/nStripAcross_);
}

double Module::getHighPitch() {
  XYZVector acrossV = corner_[2] - corner_[1];
  return (acrossV.R()/nStripAcross_);
}

void Module::computeBoundaries(double zError) {
  double minEta, maxEta;
  double minPhi, maxPhi;
  double thisEta, thisPhi;
  uint i,j;
  double z;

  // Compute eta boundaries:
  minEta=corner_[0].Eta();
  maxEta=corner_[0].Eta();
  for (z=-zError*5; z<=zError*5; z+=zError*5) {
    for (i=0; i<4; i++) {
      thisEta = (corner_[i]+XYZVector(0,0,z)).Eta();
      if (thisEta>maxEta) maxEta=thisEta;
      if (thisEta<minEta) minEta=thisEta;
    }
    for (i=0; i<4; i++) {
      j=i+1; if (j==4) j=0;
      XYZVector sideMiddle((corner_[i]+corner_[j])/2.+XYZVector(0,0,z));
      thisEta = sideMiddle.Eta();
      if (thisEta>maxEta) maxEta=thisEta;
      if (thisEta<minEta) minEta=thisEta;    
    }
  }

  // Compute Phi boundaries
  minPhi=corner_[0].Phi();
  maxPhi=corner_[0].Phi();
  for (i=1; i<4; i++) {
    thisPhi = corner_[i].Phi();
    if (thisPhi>maxPhi) maxPhi=thisPhi;
    if (thisPhi<minPhi) minPhi=thisPhi;
  }
  for (i=0; i<4; i++) {
    j=i+1; if (j==4) j=0;
    XYZVector sideMiddle((corner_[i]+corner_[j])/2.);
    thisPhi = sideMiddle.Phi();
    if (thisPhi>maxPhi) maxPhi=thisPhi;
    if (thisPhi<minPhi) minPhi=thisPhi;
  }


  // Assign boundaries to the module's
  // member variables
  boundaryMinEta_=minEta;
  boundaryMaxEta_=maxEta;
  boundaryMinPhi_=minPhi;
  boundaryMaxPhi_=maxPhi;

  // Assign logical member variables
  if ((maxPhi-minPhi)>M_PI) isAcrossPi_=true; else isAcrossPi_=false;
  computedBoundaries_ = true;

}


bool Module::couldHit(double eta, double phi) {

  if (!computedBoundaries_) {
    std::cerr << "ERROR: missing boundaries after computeBoundaries()" << std::endl;
    return true;
  }
  
  // If track's eta is outside the boundaries, then the result is certainly false
  if ((eta<boundaryMinEta_)||(eta>boundaryMaxEta_)) return false;

  // withinPhi represents whether the track's phi
  // is comprised within minimum and maximum phi of the module
  bool withinPhi = ((phi>boundaryMinPhi_)&&(phi<boundaryMaxPhi_));

  // Explanation:
  // if the modules lies across the M_PI, then in order for the track to cross the module
  // we need the track's phi to be greater than maximum phi (M_PI-f)
  // AND less then minimum phi (-(M_PI-f))

  // Table of truth:
  // withinPhi   isAcrossPi_    withinPhi XOR isAcrossPi_
  // 0           0              0
  // 0           1              1
  // 1           0              1
  // 1           1              0

  // Thus XOR (^) is our operator...
  return ((withinPhi^isAcrossPi_));
  
}


/******************/
/*                */
/* Barrel module  */
/*                */
/******************/

BarrelModule::~BarrelModule() {
}

BarrelModule::BarrelModule(double waferDiameter, double heightOverWidth) : Module(waferDiameter) {
  setSensorGeometry(heightOverWidth);
}

BarrelModule::BarrelModule(double heightOverWidth /*=1*/ ) : Module() {
  setSensorGeometry(heightOverWidth);
}

void BarrelModule::setSensorGeometry(double heightOverWidth) {
  double cornerAngle = atan(heightOverWidth);
  height_ = waferDiameter_ * sin (cornerAngle);
  width_  = waferDiameter_ * cos (cornerAngle);
  area_ = height_*width_;

  // Upon creation the module is placed horizontally
  // with its center laying on 0,0,0
  for (int i=0; i<4; i++) corner_[i].SetY(0.);
  corner_[0].SetX((-1)*width_/2);
  corner_[1].SetX(width_/2);
  corner_[2].SetX(width_/2);
  corner_[3].SetX((-1)*width_/2);

  corner_[0].SetZ(height_/2);
  corner_[1].SetZ(height_/2);
  corner_[2].SetZ((-1)*height_/2);
  corner_[3].SetZ((-1)*height_/2);  
}


// This would not measure correctly the margin! TODO
edge BarrelModule::getEdgeZSide(int direction, double margin /*= 0*/) {
  BarrelModule fakeModule;



  XYZVector* marginBorder;

  int nextIndex;
  for (int i=0; i<4; i++) {
    nextIndex=(i+1) % 4;
    marginBorder = marginBorderSide(margin, i);
    fakeModule.corner_[i]=(*marginBorder);
    delete marginBorder;
  }

  return fakeModule.getEdgeZ(direction, 0);
}

edge BarrelModule::getEdgeZ(int direction, double margin /*= 0*/ ) {
  edge result;
  
  if (direction==0) {
    std::cerr << "BarrelModule::getEdgeZ direction was chosen as zero." << std::endl;
    result.first=-1;
    result.second=-1;
    return result;
  }

  XYZVector* referencePoint[4];
  
  if (margin==0) {
    for (int i=0; i<4; i++) {
      referencePoint[i]=&corner_[i];      
    }
  } else {
    for (int i=0; i<4; i++) {
      referencePoint[i]= marginBorder(margin,margin,i);
    }    
  }



  if (direction<0) {
    // Get the point with lowest Z
    result.first=referencePoint[0]->Z();
    result.second=0;
    for (int i=1; i<4; i++) {
      if (referencePoint[i]->Z()<result.first) {
	result.first=referencePoint[i]->Z();
	result.second=i;
      }
    }
  } else {
    // Get the point with highest Z
    result.first=referencePoint[0]->Z();
    result.second=0;
    for (int i=1; i<4; i++) {
      if (referencePoint[i]->Z()>result.first) {
	result.first=referencePoint[i]->Z();
	result.second=i;
      }
    }
  }


  if (margin!=0) {
    for (int i=0; i<4; i++)
      delete referencePoint[i];
  }


  return result;
  
}

int BarrelModule::setEdgeZ(double newZ, int direction) {
  if (direction==0) {
    std::cerr << "BarrelModule::setEdgeZ direction was chosen as zero." << std::endl;
    return -1;
  }

  edge edgePoint = getEdgeZ((-1)*direction);
  XYZVector delta(0,
		  0,
		  newZ-edgePoint.first);
  
  
  for (int i=0; i<4; i++) {
    corner_[i]+=delta;
  }

  return 1;
}

// This would not measure correctly the margin! TODO
edge BarrelModule::getEdgePhiSide(int direction, double margin /*= 0*/) {


  BarrelModule fakeModule;
  XYZVector* marginBorder;

  int nextIndex;
  for (int i=0; i<4; i++) {
    nextIndex=i+1;
    if (nextIndex==4) nextIndex=0;
    marginBorder = marginBorderSide(margin, i);
    fakeModule.corner_[i]=(*marginBorder);
    delete marginBorder;
  }

  return fakeModule.getEdgePhi(direction, 0);
}

edge BarrelModule::getEdgePhi(int direction, double margin /*= 0*/) {
  bool thereIsThirdQuad=false;
  bool thereIsFourthQuad=false;
  edge result;

  if (direction==0) {
    std::cerr << "In BarrelModule::getEdgePhi direction was chosen as zero." << std::endl;
    result.first=-1;
    result.second=-1;
    return result;
  }


  XYZVector* referencePoint[4];

  if (margin==0) {
    for (int i=0; i<4; i++) {
      referencePoint[i]=&corner_[i];      
    }
  } else {
    for (int i=0; i<4; i++) {
      referencePoint[i]=marginBorder(margin,margin,i);
    }    
  }

  double phi[4];


  for (int i=0; i<4; i++) {
    phi[i]=referencePoint[i]->Phi();
    if ((phi[i]>=M_PI/2.)&&(phi[i]<=M_PI)) thereIsThirdQuad=true;
    if ((phi[i]<=(-1)*M_PI/2.)&&(phi[i]>=(-1)*M_PI)) thereIsFourthQuad=true;
  }

  // Adjust negative values to be in the range (0, 2Pi)
  // in case we have phis in the third AND fourth quadrant
  if ((thereIsThirdQuad) && (thereIsFourthQuad)) {
    for (int i=0; i<4; i++) {
      if (phi[i]<0) phi[i]+=2*M_PI;
    }
  }

  // If direction is greater than zero we are looking for the highest phi
  if (direction>0) {
    result.first=phi[0];
    result.second=0;
    for (int i=1; i<4; i++) {
      if (phi[i]>result.first) {
	result.first=phi[i];
	result.second=i;
      }
    }
  }

  // If direction is lower than zero we are looking for the lowest phi
  if (direction<0) {
    result.first=phi[0];
    result.second=0;
    for (int i=1; i<4; i++) {
      if (phi[i]<result.first) {
	result.first=phi[i];
	result.second=i;
      }
    }
  }


  if (margin!=0) {
    for (int i=0; i<4; i++)
      delete referencePoint[i];
  }


  return result;
}

int BarrelModule::setEdgePhi(double newPhi, int direction) {

  if (direction==0) {
    std::cerr << "In BarrelModule::getEdgePhi direction was chosen as zero." << std::endl;
    return -1;
  }

  double startEdgePhi = getEdgePhi(-1*direction).first;
  rotatePhi(newPhi-startEdgePhi);

  return 1;
}




/******************/
/*                */
/* Endcap module  */
/*                */
/******************/

EndcapModule::~EndcapModule() {
}

// Just to make an endcapmodule out of a module, without touching the geometry
EndcapModule::EndcapModule(const Module& aModule) : Module(aModule) {
  cut_ = false;

  phiWidth_ = 0;
  widthLo_  = 0;
  widthHi_  = 0;
  dist_     = 0;
}

EndcapModule::EndcapModule(const Module& sampleModule, double alpha, double d, double maxRho /* = -1 */) : Module(sampleModule) {
  setSensorGeometry(alpha, d, maxRho);
}

EndcapModule::EndcapModule(const EndcapModule& sampleModule, double alpha, double d, double maxRho /* = -1 */) : Module(sampleModule){
  
  phiWidth_ = sampleModule.phiWidth_;
  widthLo_ = sampleModule.widthLo_;
  widthHi_ = sampleModule.widthHi_;
  dist_ = sampleModule.dist_;
  cut_ = sampleModule.cut_;
  lost_ = sampleModule.lost_;
  ring_ = sampleModule.ring_;  
  disk_ = sampleModule.disk_;
  
  setSensorGeometry(alpha, d, maxRho);
}

EndcapModule::EndcapModule(double alpha, double d, double maxRho /* = -1 */) : Module() {
  setSensorGeometry(alpha, d, maxRho);
}

void EndcapModule::setSensorGeometry(double alpha, double d, double maxRho /* = -1 */) {
  
  cut_=false;

  double r = waferDiameter_/2.;

  // We need the half angle covered by the module
  double phi;
  phi = alpha/2.;
  
  double gamma2;
  double h1, h2, b1, b2, dfar, l;

  // The short (half)base
  h1 = d * tan(phi);
  // Distance of short base from the wafer center
  b1 = sqrt(pow(r,2)-pow(h1,2));

  // y coordinate of the wafer center
  l = b1 + d;
  // Distance of the far angle form the z axis
  gamma2 = l*cos(phi) + sqrt(pow(r,2)-pow(l*sin(phi),2));

  // The long (half)base
  h2 = gamma2 * sin(phi);
  // Distance of long base from the z axis
  dfar = gamma2 * cos(phi);

  // The distance of the long base from the wafer center
  b2 =  sqrt(pow(r,2)-pow(h2,2));
  // NOTE: in principle we don't need to compute b2 to get the
  // module's corner coordinates. Still we use this way of computing
  // b2 to ensure maximum precision to the corner being placed on
  // the wafer's edge

  // Add a check: if the module overcomes the max rho
  // it must be cut.
  if (maxRho>0) {
    if ((d+b1+b2)>maxRho) {
      lost_=d+b1+b2-maxRho;
      b1=0;
      b2=maxRho-d;
      h2=h1/d*maxRho;
      cut_=true;
    }
  }

  // Some member variable computing:
  area_     = (b1+b2)*(h2+h1);
  height_   = b1+b2;
  phiWidth_ = 2*phi;
  widthLo_  = 2*h1;
  widthHi_  = 2*h2;
  dist_     = d;

  // Some checks
  // std::cerr << "r1-r: " << sqrt(pow(h1,2)+pow(b1,2)) - r << std::endl;
  // std::cerr << "r2-r: " << sqrt(pow(h2,2)+pow(b2,2)) - r << std::endl;
  // std::cerr << "Area: " << area_ << std::endl;

  // The four corners:
  // std::cerr << "h1: " << h1   << std::endl;
  // std::cerr << "d1: " << d    << std::endl;
  // std::cerr << "h2: " << h2   << std::endl;
  // std::cerr << "d2: " << d+b1+b2 << std::endl;
  
  corner_[0] = XYZVector(h1,d,0);
  corner_[1] = XYZVector(h2,d+b1+b2,0);
  corner_[2] = XYZVector(-1*h2,d+b1+b2,0);
  corner_[3] = XYZVector(-1*h1,d,0);

}



// ############################
// #                          #
// #  Sensor TAG should match #
// #                          #
// ############################

std::string BarrelModule::getSensorTag() {
  std::string result ;
  
  result = "Barrel";
  if (tag_!="") {
    result += "/Tag=" + tag_;
  }
  result += "/Width=" +
    boost::lexical_cast<std::string>(width_) +
    "/Height=" +
    boost::lexical_cast<std::string>(height_) +
    "/Thickness=" +
    boost::lexical_cast<std::string>(thickness_) +
    "/Type=" + type_;
  


  return result;
}


std::string EndcapModule::getSensorTag() {
  std::string result ;

  result = "Endcap";
  if (tag_!="") {
    result += "/Tag=" + tag_;
  }
  result += "/WidthLo=" +
    boost::lexical_cast<std::string>(widthLo_) +
    "/WidthHi=" +
    boost::lexical_cast<std::string>(widthHi_) +
    "/Height=" +
    boost::lexical_cast<std::string>(height_) +
    "/Thickness=" +
    boost::lexical_cast<std::string>(thickness_) +
    "/Type=" + type_;
  
  return result;
}

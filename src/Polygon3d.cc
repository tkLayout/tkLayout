#include "Polygon3d.hh"

template<> XYZVector crossProduct<XYZVector>(const XYZVector& v1, const XYZVector& v2) { return v1.Cross(v2); }
template<> XYZVector unitVector<XYZVector>(const XYZVector& vector) { return vector.Unit(); }






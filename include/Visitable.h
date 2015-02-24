/**
 * @file Visitable.h
 *
 * @date 14/mag/2014
 * @author Stefano Martina
 */

#ifndef VISITABLE_H_
#define VISITABLE_H_

class Visitable {
public:
  Visitable(){};
  virtual ~Visitable(){};

  virtual void accept(GeometryVisitor& v) = 0;
  virtual void accept(ConstGeometryVisitor& v) const = 0;
};

#endif /* VISITABLE_H_ */

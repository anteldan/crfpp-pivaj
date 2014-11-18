//
//  CRF++PIVAJ -- A CRF toolkit derived from CRF++ for 
//  the PIVAJ project (see http://plair.univ-rouen.fr).
//
//
//  Copyright(C) 2014 Julien Lerouge <julien@lerouge.me>
//  Copyright(C) 2005-2007 Taku Kudo <taku@chasen.org>
//
#include <cmath>
#include "path.h"
#include "common.h"

namespace CRFPP {

void Path::calcExpectation(double *expected, double Z, size_t size) const {
  const double c = std::exp(lnode->alpha + cost + rnode->beta - Z);
  for (const int *f = fvector; *f != -1; ++f) {
    expected[*f + lnode->y * size + rnode->y] += c;
  }
}

void Path::add(Node *_lnode, Node *_rnode) {
  lnode = _lnode;
  rnode = _rnode;
  lnode->rpath.push_back(this);
  rnode->lpath.push_back(this);
}
}

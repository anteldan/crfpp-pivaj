//
//  CRF++PIVAJ -- A CRF toolkit derived from CRF++ for 
//  the PIVAJ project (see http://plair.univ-rouen.fr).
//
//
//  Copyright(C) 2014 Julien Lerouge <julien@lerouge.me>
//  Copyright(C) 2005-2007 Taku Kudo <taku@chasen.org>
//
#ifndef CRFPP_FEATURE_CACHE_H_
#define CRFPP_FEATURE_CACHE_H_

#include <vector>
#include <map>
#include "freelist.h"

namespace CRFPP {

class FeatureCache: public std::vector <int *> {
 public:
  void clear() {
    std::vector<int *>::clear();
    feature_freelist_.free();
  }

  void add(const std::vector<int> &);
  void shrink(std::map<int, int> *);

  explicit FeatureCache(): feature_freelist_(8192 * 16) {}
  virtual ~FeatureCache() {}

 private:
  FreeList<int> feature_freelist_;
};
}
#endif

//
//  CRF++PIVAJ -- A CRF toolkit derived from CRF++ for 
//  the PIVAJ project (see http://plair.univ-rouen.fr).
//
//
//  Copyright(C) 2014 Julien Lerouge <julien@lerouge.me>
//  Copyright(C) 2005-2007 Taku Kudo <taku@chasen.org>
//
#ifndef CRFPP_ENCODER_H_
#define CRFPP_ENCODER_H_

#include "crfpp.h"
#include "common.h"
#include "param.h"

namespace CRFPP {
class EncoderImpl : public Encoder {
 public:
     
   bool configure(const Param &param);

   bool configure(int argc,  char** argv);

   bool configure(const char* arg);
  
   bool learn(const char *trainfile,
              const char *templfile,
              const char *modelfile);
//  bool learn(const char *, const char *,
//             const char *,
//             bool, size_t, size_t,
//             double, double,
//             unsigned short,
//             unsigned short, int);

   bool convert(const char *text_file,
                const char* binary_file);
   
 private:
    bool textmodelfile_;
    size_t maxitr_;
    size_t freq_;
    double eta_;
    double C_;
    unsigned short thread_num_;
    unsigned short shrinking_size_;
    int algorithm_;
};
}
#endif

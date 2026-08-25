#include "classconfig.h"

// node_class

int node_cnt = 0;

namespace cc {

node_class::node_class(int number_,double x_,double y_):number(number_),x(x_),y(y_){node_cnt++;}

}
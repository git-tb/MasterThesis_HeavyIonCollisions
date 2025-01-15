#ifndef UNINITFUNC_H
#define UNINITFUNC_H

#include <functional>
#include "debugmsg.h"

std::function<double(double)> UNINITIALIZED_FUNCTION = [](double x){
    DEBUGMSG("THIS FUNCTION HAS NOT BEEN INITIALIZED");
    return -1e10;
};

#endif
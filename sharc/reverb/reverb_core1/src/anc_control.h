#pragma once

#include "panel_ir.h"
#include "fxlms.hpp"

#define FILTER_ORDER 1024
#define NLMS true
#define LAG 102

extern volatile bool anc_off;

extern FxLMS<FILTER_ORDER, IR_LEN, NLMS> anc;

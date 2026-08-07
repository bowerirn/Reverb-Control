#pragma once

#include "panel_ir.h"
#include "fxlms.hpp"

#define FILTER_ORDER 512
#define NLMS true

extern volatile bool anc_off;

extern FxLMS<FILTER_ORDER, IR_LEN, NLMS> anc;

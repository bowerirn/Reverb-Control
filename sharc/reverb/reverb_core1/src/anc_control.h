#pragma once

#include "panel_ir.h"
#include "fxlms.hpp"

#define FILTER_ORDER 512
#define NLMS true

extern FxLMS<FILTER_ORDER, IR_LEN, NLMS> anc;

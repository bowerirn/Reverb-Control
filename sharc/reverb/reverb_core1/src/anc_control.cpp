#include "anc_control.h"
#include "panel_ir.h"

FxLMS<FILTER_ORDER, IR_LEN, NLMS> anc(panel_ir);

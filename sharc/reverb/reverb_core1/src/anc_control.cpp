#include "anc_control.h"
#include "panel_ir.h"

volatile bool anc_off = true;

FxLMS<FILTER_ORDER, IR_LEN, NLMS> anc(error_ir, ref_ir);

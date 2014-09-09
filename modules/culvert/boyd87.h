/* Implementation of culvert model based on
 *
 * "Generalised Head-Discharge Equations for Culverts", M. J. Boyd,
 * Fourth national local government engineering conference, Perth, 17-20
 * August 1987.  
 *
 */

/**
 * Flow rate for box culvert, inlet control.
 * HW: headwater depth
 * B: box width.
 * D: box height.
 * type: entrance type.
 * g: acceleration of gravity.
 */
double Q_inlet_box (double HW, double B, double D, int type, double g);

/**
 * Flow rate for circular pipe culvert, inlet control.
 * HW: headwater depth
 * D: diameter.
 * type: entrance type.
 * g: acceleration of gravity.
 */
double Q_inlet_pipe (double HW, double D, int type, double g);

/**
 * Flow rate for box culvert, outlet control.
 * HW: headwater depth
 * TW: tailwater depth.
 * B: box width.
 * D: box height.
 * S0: slope.
 * L: length.
 * n: Manning friction coefficient.
 * ke: entrance loss coefficient.
 * g: acceleration of gravity.
 */
double Q_outlet_box (double HW, double TW,
		     double B, double D,
		     double S0, double L,
		     double n, double ke,
		     double g);

/**
 * Flow rate for circular pipe culvert, outlet control.
 * HW: headwater depth
 * TW: tailwater depth.
 * D: diameter.
 * S0: slope.
 * L: length.
 * n: Manning friction coefficient.
 * ke: entrance loss coefficient.
 * g: acceleration of gravity.
 */
double Q_outlet_pipe (double HW, double TW,
		      double D,
		      double S0, double L,
		      double n, double ke,
		      double g);

/**
 * Flow rate for circular pipe culvert.
 * HW: headwater depth
 * TW: tailwater depth.
 * D: diameter.
 * type: entrance type.
 * S0: slope.
 * L: length.
 * n: Manning friction coefficient.
 * ke: entrance loss coefficient.
 * g: acceleration of gravity.
 */
double Q_pipe (double HW, double TW,
	       double D,
	       int type,
	       double S0, double L,
	       double n, double ke,
	       double g);

/**
 * Flow rate for circular pipe culvert, outlet control.
 * HW: headwater depth
 * TW: tailwater depth.
 * B: box width.
 * D: box height.
 * type: entrance type.
 * S0: slope.
 * L: length.
 * n: Manning friction coefficient.
 * ke: entrance loss coefficient.
 * g: acceleration of gravity.
 */
double Q_box (double HW, double TW,
	      double B, double D,
	      int type,
	      double S0, double L,
	      double n, double ke,
	      double g);

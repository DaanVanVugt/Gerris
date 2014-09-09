/* Gerris - The GNU Flow Solver
 * Copyright (C) 2009 National Institute of Water and Atmospheric Research
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 * 02111-1307, USA.  
 */

#include "cfortrantypes.h"

/**
 * Initialisation of wavewatch
 */
extern void GFSW3INIT (void);

/**
 * imported from w3srcemd.f90:W3SRCE
 * @SPEC: Spectrum (action) in 1-D form (NK*NTH elements): INPUT/OUTPUT
 * @SP2D: Spectrum (action) in 2-D form (NK*NTH elements): INPUT/OUTPUT
 * @ALPHA: Nondimensional 1-D spectrum (NK)                OUTPUT
 * @WN1: Discrete wavenumbers (NK)                         INPUT
 * @CG1: Discrete group velocities (NK)                    INPUT
 * @DEPTH: Depth                                           INPUT
 * @U10ABS: Wind speed at reference height                 INPUT
 * @U10DIR: Wind direction at reference height             INPUT
 * @USTAR: Friction velocity                               INPUT/OUTPUT
 * @EMEAN: Mean energy                                     OUTPUT (maybe model dependent?)
 * @FMEAN: Mean frequency                                  OUTPUT (maybe model dependent?)
 * @WMEAN: Mean wavenumber                                 OUTPUT (maybe model dependent?)
 * @AMAX: Maximum energy                                   OUTPUT
 * @FPI: Peak-input frequency                              INPUT/OUTPUT
 * @CD: Drag coefficient                                   OUTPUT (maybe model dependent?)
 * @Z0: Roughness length                                   OUTPUT (maybe model dependent?)
 * @DTDYN: Average dynamic time step                       OUTPUT
 * @FCUT: Cut-off frequency for tail                       OUTPUT
 * @DTG: Global time step                                  INPUT
 * @DTMIN: Minimum (partial) time step                     INPUT
 * @DTMAX: Maximum (partial) time step                     INPUT
 */
extern void W3SRCE (REAL * SPEC,
		    REAL * SP2D,
		    REAL * ALPHA,
		    REAL * WN1, 
		    REAL * CG1, 
		    REAL * DEPTH, 
		    REAL * U10ABS, 
		    REAL * U10DIR, 
		    REAL * USTAR,
		    REAL * EMEAN,
		    REAL * FMEAN,
		    REAL * WMEAN,
		    REAL * AMAX,
		    REAL * FPI,
		    REAL * CD,
		    REAL * Z0,
		    REAL * DTDYN,
		    REAL * FCUT,
		    REAL * DTG,
		    REAL * DTMIN,
		    REAL * DTMAX);

/*
    Copyright (c) 2023 Sam Pluta
    Built upon the improved FFT and IFFT UGens for SuperCollider 3
    Copyright (c) 2007-2008 Dan Stowell, incorporating code from
    SuperCollider 3 Copyright (c) 2002 James McCartney.
    All rights reserved.

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA
*/

#include "FFT_UGens.h"
#include "SC_PlugIn.h"
#include "SC_PlugIn.hpp"
#include "SC_fftlib.h"
#include <cstdio>

#if defined(__APPLE__) && !defined(SC_IPHONE)
    #include <Accelerate/Accelerate.h>
#endif


InterfaceTable *ft;



//////////////////////////////////////////////////////////////////////////////////////////////////
//constants and functions

constexpr double PI = 3.14159265358979323846;
constexpr double TWOPI = 2 * PI;

float phaseminus(float phase1, float phase2) {
  float dif = phase1 - phase2;
  if(dif < -pi)
    return dif + twopi;
  if(dif > pi)
    return dif - twopi;
  return dif;
}

float interpolateUpDown(float in, float fb, float atkCoeff, float dcyCoeff) {
	return (fb - in) * ((in > fb) ? atkCoeff : dcyCoeff) + in;
}



//////////////////////////////////////////////////////////////////////////////////////////////////

struct PV_Freezish : public PV_Unit {
	int m_numbins;
	float *m_mags, m_dc, m_nyq;
	float *m_prevPhases, *m_difPhases;
	SndBuf *m_buf;
	int m_stage;
};

extern "C" {
	#include "SC_fftlib.h"

	void PV_Freezish_Ctor(PV_Freezish *unit);
	void PV_Freezish_next(PV_Freezish* unit, int inNumSamples);
	void PV_Freezish_Dtor(PV_Freezish *unit);
}

//////////////////////////////////////////////////////////////////////////////////////////////////

//// PV_Freezish:
//// Like Josh Parmenter's PV_Freeze,
//// except 'freeze' is a continuous coefficient:
//// 0.0 = all energy from input bins
//// 1.0 = all energy from the last output bin
//// or blend these (IIR style)


void PV_Freezish_next(PV_Freezish *unit, int inNumSamples)
{
	// printf("ZIN0(0) = %f\n", ZIN0(0));
	PV_GET_BUF

	// printf("numbins = %d, buf = %d, m_stage = %d\n", numbins, ibufnum, unit->m_stage);

	float *mags = unit->m_mags;
	float *prevPhases = unit->m_prevPhases;
	SCPolarBuf *p = ToPolarApx(buf);
	float *difPhases = unit->m_difPhases;
	float atkCoeff = ZIN0(1);
	float dcyCoeff = ZIN0(2);

	switch(unit->m_stage) {
	case 0:
		// no FFT frames are ready yet. Only allocate storage
		unit->m_mags = (float*)RTAlloc(unit->mWorld, numbins * sizeof(float));
		unit->m_difPhases = (float*)RTAlloc(unit->mWorld, numbins * sizeof(float));
		unit->m_prevPhases = (float*)RTAlloc(unit->mWorld, numbins * sizeof(float));
		unit->m_numbins = numbins;
		unit->m_stage = 1;
		break;
	case 1:
		// exactly one FFT frame is ready. Save its mags and phases. output = input
		for (int i=0; i<numbins; ++i) {
			mags[i] = p->bin[i].mag;
			prevPhases[i] = p->bin[i].phase;
		}
		unit->m_dc = p->dc;
		unit->m_nyq = p->nyq;
		unit->m_stage = 2;
		break;
	case 2:
		// second FFT frame: now it's meaningful to populate difPhases
		// difPhases = difference in phase between the current and last frame
		// if freezing: keep old mags and new phases
		// below, we will do outphase = prevphase + difphase
		// but now, difphase will be inphase - prevphase so outphase = inphase
		for (int i=0; i<numbins; ++i) {
		  float in = p->bin[i].mag;
		  float fb = mags[i];
		  float out = interpolateUpDown(in, fb, atkCoeff, dcyCoeff);
		  p->bin[i].mag = out;
		  mags[i] = out;

		  difPhases[i] = phaseminus(p->bin[i].phase, prevPhases[i]);
		  prevPhases[i] = p->bin[i].phase;
		}	      
		p->dc = unit->m_dc;
		p->nyq = unit->m_nyq;
		unit->m_stage = 3;
		break;
	case 3:
		// normal case
		// freeze previous magnitudes
		for (int i=0; i<numbins; ++i) {
		  float in = p->bin[i].mag;
		  float fb = mags[i];
		  float coeff = (in > fb) ? atkCoeff : dcyCoeff;
		  float out = interpolateUpDown(in, fb, atkCoeff, dcyCoeff);
		  p->bin[i].mag = out;
		  mags[i] = out;

		  // with this, however, pretty much anything with a 'freeze' < 1.0
		  // converges onto a color that seems related to the fft period
		  // I guess perhaps that difPhases decays toward 0
		  float difPhaseIn = phaseminus(p->bin[i].phase, prevPhases[i]);
		  float difPhaseFb = difPhases[i];
		  // interpolateUpDown() can't be used b/c 'mag' comparison is relevant, not phase

		  difPhases[i] = (difPhaseFb - difPhaseIn) * coeff + difPhaseIn;
		  float phase = prevPhases[i] + difPhases[i];
		  // if(i == 5) printf("difPhases[5] = %f\n", difPhases[i]);
		  
		  if (phase > pi) /* wrap the phase */
		      phase -= twopi;
		  if (phase < -pi)
		      phase += twopi;
		  p->bin[i].phase = phase;
		  prevPhases[i] = phase;
		}
		unit->m_dc = p->dc = interpolateUpDown(p->dc, unit->m_dc, atkCoeff, dcyCoeff);
		unit->m_nyq = p->nyq = interpolateUpDown(p->nyq, unit->m_nyq, atkCoeff, dcyCoeff);
		break;
	}
}


void PV_Freezish_Ctor(PV_Freezish* unit)
{
	SETCALC(PV_Freezish_next);
	ZOUT0(0) = ZIN0(0);
	unit->m_stage = 0;
}

void PV_Freezish_Dtor(PV_Freezish* unit)
{
	RTFree(unit->mWorld, unit->m_mags);
	RTFree(unit->mWorld, unit->m_difPhases);
	RTFree(unit->mWorld, unit->m_prevPhases);
}



void init_SCComplex(InterfaceTable *inTable);

#define DefinePVUnit(name) (*ft->fDefineUnit)(#name, sizeof(PV_Unit), (UnitCtorFunc)&name##_Ctor, 0, 0);


PluginLoad(PV_Freezy)
{
    // InterfaceTable *inTable implicitly given as argument to the load function
    ft = inTable; // store pointer to InterfaceTable
    init_SCComplex(inTable);

    DefineDtorUnit(PV_Freezish);
}


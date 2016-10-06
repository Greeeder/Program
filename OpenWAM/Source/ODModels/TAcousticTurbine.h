// ---------------------------------------------------------------------------

#ifndef TAcousticTurbineH
#define TAcousticTurbineH
// ---------------------------------------------------------------------------

#include "TTubo.h"

/**
 * @file TAcousticTurbine.h
 * @author Francisco Jose Arnau Martinez <farnau@mot.upv.es>
 * @version 0.1
 *
 * @section LICENSE
 *
 * This file is part of OpenWAM.
 *
 * OpenWAM is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * OpenWAM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
 *
 * @section DESCRIPTION
 * The TurboBearings class represent the bearing system in a turbocharger.
 *
 * It has methods to compute the power losses in the bearing system. This file
 * has the interface of TurboBearings.
 *
 */

class TAcousticTurbine {
  private:

	std::vector<TTubo*> FInletPipe;
	std::vector<TTubo*> FVolute;
	TTubo *FOutletPipe;

	iVector FInletPipeID;
	iVector FVoluteID;
	int FOutletPipeID;

  public:
	TAcousticTurbine(iVector InletPipeID, iVector VoluteID, int OutletPipeID);

	TAcousticTurbine();

	~TAcousticTurbine();

	virtual double T3(int i) const;

	virtual double T3() const;

	virtual double T30(int i) const;

	virtual double P3(int i) const;

	virtual double P3() const;

	virtual double P30(int i) const;

	virtual double ExpRatio(int i) const;

	virtual double ExpRatio() const;

	virtual double DiabEfficiency(int i) const;

	virtual double P4() const;

	virtual double T4() const;

	virtual double R(int i) const {
		return FInletPipe[i]->GetRMezcla(0);
	}
	;

	virtual double SIn(int i) const {
		return FInletPipe[i]->GetArea(0);
	}
	;

	virtual double DIn(int i) const {
		return FInletPipe[i]->GetDiametro(0);
	}
	;

	virtual double DInTot() const;

	virtual double DOut() const {
		return FOutletPipe->GetDiametro(FOutletPipe->getNin() - 1);
	}
	;

	virtual double SOut() const {
		return FOutletPipe->GetArea(FOutletPipe->getNin() - 1);
	}
	;

	virtual double MassIn(int i) const;

	virtual double MassIn() const;

	virtual double MassOut() const;

	void AsignInPipe(TTubo **Pipe, int i) {
		FInletPipe[i] = Pipe[FInletPipeID[i] - 1];
	}
	;

	void AsignOutPipe(TTubo **Pipe) {
		FOutletPipe = Pipe[FOutletPipeID - 1];
	}
	;

	/*!
	 * \brief Volute outlet pressure.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * 
	 * \return Volute outlet pressure. [Pa]
	 */
	virtual RowVector VoluteOutletp() const;

};

#endif


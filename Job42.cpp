// Mi headers
#include "/pbs/home/m/mpetro/PROGRAMS/MiModule/include/MiEvent.h" 
#include "/pbs/home/m/mpetro/PROGRAMS/MiModule/include/MiSDVisuHit.h" 
#include "/pbs/home/m/mpetro/PROGRAMS/MiModule/include/MiFilters.h"

#include "TLatex.h"
#include "TRandom.h"
#include "TVector3.h"
#include "TMath.h"


#include <string>
#include <iostream>
#include <sstream>
#include <algorithm>

R__LOAD_LIBRARY(/pbs/home/m/mpetro/PROGRAMS/MiModule/lib/libMiModule.so);

////////////// Function used in script
/////////////////////////////////////////////////////////////

const double ELECTRON_MASS_MEV = 0.5109989461; // in [MeV]
const double LIGHT_SPEED = 299792458 * 1e-9 * 1000; // in [mm/ns]


TVector3* get_vertex_vector(MiEvent*  _eve, string _position, int _trackID) // returns the step position of the hit_step when it first enters tracking gas
{
	TVector3* vertexVector;
	if ( _position == "calo" )
	{
		for(int j = 0;j < _eve->getPTD()->getpart(_trackID)->getvertexv()->size();j++)
		{
			if(
				_eve->getPTD()->getpart(_trackID)->getvertex(j)->getpos() == "xcalo" || 
				_eve->getPTD()->getpart(_trackID)->getvertex(j)->getpos() == "calo"  ||
				_eve->getPTD()->getpart(_trackID)->getvertex(j)->getpos() == "gveto" 
			)
			{
				vertexVector = new TVector3(
					_eve->getPTD()->getpart(_trackID)->getvertex(j)->getr()->getX(),
					_eve->getPTD()->getpart(_trackID)->getvertex(j)->getr()->getY(),
					_eve->getPTD()->getpart(_trackID)->getvertex(j)->getr()->getZ()
				);
			}
		}
	}
	else
	{
		for(int j = 0;j < _eve->getPTD()->getpart(_trackID)->getvertexv()->size();j++)
		{
			if(_eve->getPTD()->getpart(_trackID)->getvertex(j)->getpos() == _position)
			{
				vertexVector = new TVector3(
					_eve->getPTD()->getpart(_trackID)->getvertex(j)->getr()->getX(),
					_eve->getPTD()->getpart(_trackID)->getvertex(j)->getr()->getY(),
					_eve->getPTD()->getpart(_trackID)->getvertex(j)->getr()->getZ()
				);
			}
		}
	}

	return vertexVector; //return -1 if the particle never left source foil (happens when we get "fakeItTillYouMakeIt events")
}


float get_distance(TVector3* v1, TVector3* v2)
{
	float xDifSquared = pow(v1->X() - v2->X(), 2);  // (x1 - x2)^2
	float yDifSquared = pow(v1->Y() - v2->Y(), 2);  // (y1 - y2)^2
	float zDifSquared = pow(v1->Z() - v2->Z(), 2);  // (z1 - z2)^2

	float distance = sqrt( xDifSquared + yDifSquared + zDifSquared );

	return distance; // sqrt( x^2 + y^2 + z^2 )
}

bool is_same_calo_gid(MiGID* cdGID,  MiGID* sdGID )
{
	if(
		cdGID->gettype() ==  sdGID->gettype() &&
		cdGID->getmodule() ==  sdGID->getmodule() &&
		cdGID->getside() ==  sdGID->getside() &&
		cdGID->getwall() ==  sdGID->getwall() &&
		cdGID->getcolumn() ==  sdGID->getcolumn() &&
		cdGID->getrow() ==  sdGID->getrow()
	)
	{
		return true;
	}
	return false;
}

float_t get_SD_energy(MiEvent* _eve, int _trackID)
{
	float_t	E = 0.0;
	for ( auto & SDCaloHit : *_eve->getSD()->getcalohitv() )
	{
		if( is_same_calo_gid(_eve->getCD()->getcalohit(_trackID)->getGID(), SDCaloHit.getGID() ) )
		{
			E	+= SDCaloHit.getE(); // there are sometimes multiple SD hits where one will be large and then a few very small (fraction of MeV)
		}
	}
	return E;
}


double max_energy(double e1, double e2)
{
	if(e1 >= e2)
	{
		return e1;
	}
	else
	{
		return e2;
	}
}

double min_energy(double e1, double e2)
{
	if(e1 <= e2)
	{
		return e1;
	}
	else
	{
		return e2;
	}
}


double get_beta(double _E)
{
    return TMath::Sqrt(_E * (_E + 2 * ELECTRON_MASS_MEV)) / (_E + ELECTRON_MASS_MEV);
}

double get_tTOF(double _l, double _beta)
{
    return _l/(_beta * LIGHT_SPEED);
}


////////////// MAIN BLOCK OF CODE
/////////////////////////////////////////////////////////////
void Job42()
{
	
////////////// Initialize File names/paths
/////////////////////////////////////////////////////////////
	const char* inFileName                 = "Default.root";													 //FOR TESTING PURPOSES USING ONLY FOLDER 0/
	const char* outFileName  = "output.root";
	TFile* 	    outTFile  = new TFile(outFileName, "RECREATE");

	int nFiles = 1;

////////////// Initialize variables to be saved
/////////////////////////////////////////////////////////////
	float   phi, sameSide;
	float   x1Reconstructed, y1Reconstructed, z1Reconstructed, x2Reconstructed, y2Reconstructed, z2Reconstructed;
	float   trackLength1, trackLength2;
	float   dz, dy, r;
	float   caloTime1, caloTime2, deltaCaloTime;
	float   beta1, beta2;
	float   tInt, tExt;



	TVector3 dir1;
	TVector3 dir2;

	float   reconstructedEnergy1, reconstructedEnergy2;
	float   simulatedEnergy1, simulatedEnergy2;
	float   sumE, maxE, minE, avgE, singleE;
	float   sumEsimu;

	float   Pint, Pext;
	float   lPint, lPext;

////////////// Saving Data
/////////////////////////////////////////////////////////////
	TTree* tree 			= new TTree("tree","tree");

	tree->Branch("phi", &phi, "phi/F");
	tree->Branch("sameSide", &sameSide, "sameSide/F");
	tree->Branch("dz", &dz, "dz/F");
	tree->Branch("dy", &dy, "dy/F");
	tree->Branch("r", &r, "r/F");
	tree->Branch("reconstructedEnergy1", &reconstructedEnergy1, "reconstructedEnergy1/F");
	tree->Branch("reconstructedEnergy2", &reconstructedEnergy2, "reconstructedEnergy2/F");
	tree->Branch("simulatedEnergy1", &simulatedEnergy1, "simulatedEnergy1/F");
	tree->Branch("simulatedEnergy2", &simulatedEnergy2, "simulatedEnergy2/F");
	tree->Branch("sumEsimu", &sumEsimu, "sumEsimu/F");

	tree->Branch("sumE", &sumE, "sumE/F");
	tree->Branch("maxE", &maxE, "maxE/F");
	tree->Branch("minE", &minE, "minE/F");
	tree->Branch("avgE", &avgE, "avgE/F");
	tree->Branch("singleE", &singleE, "singleE/F");
	tree->Branch("trackLength1", &trackLength1, "trackLength1/F");
	tree->Branch("trackLength2", &trackLength2, "trackLength2/F");

	tree->Branch("Pint", &Pint, "Pint/F");
	tree->Branch("Pext", &Pext, "Pext/F");
	tree->Branch("lPint", &lPint, "lPint/F");
	tree->Branch("lPext", &lPext, "lPext/F");

	tree->Branch("x1Reconstructed", &x1Reconstructed, "x1Reconstructed/F");
	tree->Branch("y1Reconstructed", &y1Reconstructed, "y1Reconstructed/F");
	tree->Branch("z1Reconstructed", &z1Reconstructed, "z1Reconstructed/F");
	tree->Branch("x2Reconstructed", &x2Reconstructed, "x2Reconstructed/F");
	tree->Branch("y2Reconstructed", &y2Reconstructed, "y2Reconstructed/F");
	tree->Branch("z2Reconstructed", &z2Reconstructed, "z2Reconstructed/F");

	tree->Branch("caloTime1", &caloTime1, "caloTime1/F");
	tree->Branch("caloTime2", &caloTime2, "caloTime2/F");
	tree->Branch("beta1", &beta1, "beta1/F");
	tree->Branch("beta2", &beta2, "beta2/F");
	tree->Branch("tInt", &tInt, "tInt/F");
	tree->Branch("tExt", &tExt, "tExt/F");
	tree->Branch("deltaCaloTime", &deltaCaloTime, "deltaCaloTime/F");

////////////// Initialize counters
/////////////////////////////////////////////////////////////
	int stepBeforeGas 				= -1;   // Represents the step just before exitting to the tracker gas volume
	int stepBeforeOM 				= -1;

	int nPassed = 0;
	double nPassedFraction = 0.0;

	TRandom rand;
////////////// Loop over files
/////////////////////////////////////////////////////////////
	for( int file = 0; file < nFiles; file ++)
	{
		stringstream ssInPath;
		ssInPath << inFileName;

		if(gSystem->AccessPathName(ssInPath.str().c_str())) // check whether Default.root exists
		{
	    	cout << "Default.root DOESNT EXIST - PATH: " << ssInPath.str().c_str() << endl;
		} 
		else 
		{
			cout << ssInPath.str().c_str() << endl;

	
	////////////// Initialize reading data
	/////////////////////////////////////////////////////////////
			TFile* 	  		inFile 						= new TFile(ssInPath.str().c_str());
			TTree* 			s 		  					= (TTree*) inFile->Get("Event");
			MiEvent*  		eve 						= new MiEvent();

			s->SetBranchAddress("Eventdata", &eve);
			int nEntries = s->GetEntries();

	/////////////////////////////////////////////////////////////
			for (int e = 0; e < nEntries; ++e) 
			{ 	
				s->GetEntry(e);

					
				cout << " event number = " << e << endl;

				nPassed += 1;
				nPassedFraction += 1.0/double(nEntries);

				reconstructedEnergy1	= eve->getPTD()->getpart(0)->getcalohit(0)->getE();
				reconstructedEnergy2	= eve->getPTD()->getpart(1)->getcalohit(0)->getE();

				simulatedEnergy1 = get_SD_energy(eve, 0);
				simulatedEnergy2 = get_SD_energy(eve, 1);

				sumEsimu = simulatedEnergy1 + simulatedEnergy2;


				sumE = reconstructedEnergy1 + reconstructedEnergy2;        
				maxE = max_energy(reconstructedEnergy1, reconstructedEnergy2);
				minE = min_energy(reconstructedEnergy1, reconstructedEnergy2);
				avgE = (reconstructedEnergy1 + reconstructedEnergy2) / 2.0;

				if (rand.Rndm() < 0.5)  // there's a weird way the falaise orders energies for background sources, so we need to randomize
					singleE = reconstructedEnergy2;
				else
					singleE = reconstructedEnergy1;

				TVector3* r1Reconstructed = get_vertex_vector(eve, "source foil", 0);  	// position vector of the foil vertex
				TVector3* r2Reconstructed = get_vertex_vector(eve, "source foil", 1);
				TVector3* r1AtOM = get_vertex_vector(eve, "calo", 0);     			// position vector where electron hits OM. Tracklength is calculated as sqrt(r1Reconstructed^2 + r1AtOM^2)
				TVector3* r2AtOM = get_vertex_vector(eve, "calo", 1);  

				x1Reconstructed = r1Reconstructed->X();
				y1Reconstructed = r1Reconstructed->Y();
				z1Reconstructed = r1Reconstructed->Z();
				x2Reconstructed = r2Reconstructed->X();
				y2Reconstructed = r2Reconstructed->Y();
				z2Reconstructed = r2Reconstructed->Z();

				caloTime1 = eve->getPTD()->getpart(0)->getcalohit(0)->gett();
				caloTime2 = eve->getPTD()->getpart(1)->getcalohit(0)->gett();
				deltaCaloTime = TMath::Abs(caloTime2 - caloTime1);

				dir1 = eve->getPTD()->getpart(0)->getdirectionfromfoil();
				dir2 = eve->getPTD()->getpart(1)->getdirectionfromfoil();

				phi   = dir1.Angle(dir2)*180/TMath::Pi();

				if( dir1.X() * dir2.X() > 0 )
				{
					sameSide = 1.0;
				}
				else
				{
					sameSide = 0.0;
				}

				trackLength1 = eve->getPTD()->getpart(0)->getTrackLength();
				trackLength2 = eve->getPTD()->getpart(1)->getTrackLength();

				dz = TMath::Abs(z1Reconstructed - z2Reconstructed);
				dy = TMath::Abs(y1Reconstructed - y2Reconstructed);
				r = get_distance(r1Reconstructed, r2Reconstructed);

				Pint = eve->getPint();
				Pext = eve->getPext();

				if( Pint == 0.0 )
					lPint = 110.0;
				else
					lPint = TMath::Abs(TMath::Log10(Pint)); //taking the absolute value of the log10 to avoid negative values and stretch the space
				
				if( Pext == 0.0 )
					lPext = 110.0;
				else
					lPext = TMath::Abs(TMath::Log10(Pext));

				beta1 = get_beta(simulatedEnergy1 / 1000); // convert energy from keV to MeV for beta calculation
				beta2 = get_beta(simulatedEnergy2 / 1000);

				tInt = pow(TMath::Abs(caloTime2 - caloTime1)  - (trackLength2 / (beta2 * LIGHT_SPEED)  + trackLength1 / (beta1 * LIGHT_SPEED)), 2);
				tExt = pow((caloTime2 - trackLength2 / (beta2 * LIGHT_SPEED)) - (caloTime1 - trackLength1 / (beta1 * LIGHT_SPEED)), 2);

				tree->Fill();

			} 	
		}
	}
	// cout << "fakeItTillYouMakeItCounter = " << fakeItTillYouMakeItCounter << endl; // these events are skipped in the simulation 
	cout <<" nPassed = " << nPassed << endl;
	cout <<" nPassedFraction = " << nPassedFraction*100 << "%" << endl;

	outTFile->Write();
	outTFile->Close();
}


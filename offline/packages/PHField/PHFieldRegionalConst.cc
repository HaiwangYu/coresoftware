
#include "PHFieldRegionalConst.h"

//root framework
#include <TFile.h>
#include <TNtuple.h>

#include <set>
#include <iostream>

using namespace std;
using namespace CLHEP;// units

PHFieldRegionalConst::PHFieldRegionalConst( const string &filename, const int verb, const float magfield_rescale) :
    PHField(verb),
		maxy_(4*cm), miny_(-4*cm),
		maxr_(15*cm), minr_(0*cm),
		field_val_(5*tesla)
{
}

void PHFieldRegionalConst::GetFieldValue(const double point[4], double *Bfield ) const
{
	double r = sqrt(point[0]*point[0] + point[2]*point[2]);
	double y = point[1];

	Bfield[0] = 0;
	Bfield[1] = 0;
	Bfield[2] = 0;

	if( r > minr_ and r < maxr_ and
			y > miny_ and y < maxy_) {
		Bfield[1] = field_val_;
	}

  return;
}

// debug function to print key/value pairs in map
void PHFieldRegionalConst::print_map() const
{
	std::cout
	<< "PHFieldRegionalConst::print_map: "
	<< "field_val_: " << field_val_/tesla << " tesla"
	<<std::endl;
}


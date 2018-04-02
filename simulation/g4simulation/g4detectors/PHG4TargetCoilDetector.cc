#include "PHG4TargetCoilDetector.h"

#include <phparameter/PHParameters.h>

#include <g4main/PHG4Utils.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/getClass.h>

#include <Geant4/G4Colour.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4PhysicalConstants.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4Polycone.hh>
#include <Geant4/G4Tubs.hh>
#include <Geant4/G4SubtractionSolid.hh>
#include <Geant4/G4UserLimits.hh>
#include <Geant4/G4VisAttributes.hh>

#include <cmath>
#include <sstream>
#include <algorithm>

using namespace std;

//_______________________________________________________________
PHG4TargetCoilDetector::PHG4TargetCoilDetector(PHCompositeNode *Node, PHParameters *parameters, const std::string &dnam, const int lyr)
  : PHG4Detector(Node, dnam)
  , params(parameters)
  , cylinder_physi(nullptr)
  , layer(lyr)
{
}

//_______________________________________________________________
bool PHG4TargetCoilDetector::IsInCylinder(const G4VPhysicalVolume *volume) const
{
  if (volume == cylinder_physi || volume->GetMotherLogical() == cylinder_physi->GetLogicalVolume())
  {
    return true;
  }
  return false;
}

bool PlaceHollowTube(
		G4LogicalVolume *mother,
		double place,
		G4UserLimits *g4userlimits,
		bool overlapcheck,
		const std::string& name,
		G4Material *Mouter,
		double in,
		double out,
		double half,
		double t,
		double phi0 = 0,
		double phi = twopi
		)
{
  G4VSolid *all_solid = new G4Tubs((name+"_all").c_str(),
  		in,
			out,
			half,
			phi0,
			phi);

  G4VSolid *inner_solid = new G4Tubs((name+"_inner").c_str(),
  		in+t,
			out-t,
			half-t,
			phi0,
			phi);

  G4VSolid * outer_solid = new G4SubtractionSolid((name+"_outer").c_str(),
  		all_solid,
			inner_solid,
			0,
			G4ThreeVector(0,0,0)
  		);


  G4VisAttributes *vis_outer = new G4VisAttributes();
	PHG4Utils::SetColour(vis_outer, Mouter->GetName());
	vis_outer->SetVisibility(true);
	vis_outer->SetForceSolid(true);

  G4LogicalVolume *outer_logic = new G4LogicalVolume(outer_solid,
  		Mouter,
			(name+"_outer").c_str(),
			nullptr, nullptr, g4userlimits);

  outer_logic->SetVisAttributes(vis_outer);

  new G4PVPlacement(
  		0,
			G4ThreeVector(0, 0, place),
			outer_logic,
			(name+"_outer").c_str(),
			mother, 0, false, overlapcheck);

  return true;
}

bool PlaceLayeredTube(
		G4LogicalVolume *mother,
		double place,
		G4UserLimits *g4userlimits,
		bool overlapcheck,
		const std::string& name,
		G4Material *Minner,
		G4Material *Mouter,
		double in,
		double out,
		double half,
		double t,
		double phi0 = 0,
		double phi = twopi
		)
{
  G4VSolid *all_solid = new G4Tubs((name+"_all").c_str(),
  		in,
			out,
			half,
			phi0,
			phi);

  G4VSolid *inner_solid = new G4Tubs((name+"_inner").c_str(),
  		in+t,
			out-t,
			half-t,
			phi0,
			phi);

  G4VSolid * outer_solid = new G4SubtractionSolid((name+"_outer").c_str(),
  		all_solid,
			inner_solid,
			0,
			G4ThreeVector(0,0,0)
  		);


  G4VisAttributes *vis_inner = new G4VisAttributes();
	PHG4Utils::SetColour(vis_inner, Minner->GetName());
	vis_inner->SetVisibility(true);
	vis_inner->SetForceSolid(true);

  G4LogicalVolume *inner_logic = new G4LogicalVolume(inner_solid,
  		Minner,
			(name+"_inner").c_str(),
			nullptr, nullptr, g4userlimits);

  inner_logic->SetVisAttributes(vis_inner);

  new G4PVPlacement(
  		0,
			G4ThreeVector(0, 0, place),
			inner_logic,
			(name+"_inner").c_str(),
			mother, 0, false, overlapcheck);


  G4VisAttributes *vis_outer = new G4VisAttributes();
	PHG4Utils::SetColour(vis_outer, Mouter->GetName());
	vis_outer->SetVisibility(true);
	vis_outer->SetForceSolid(true);

  G4LogicalVolume *outer_logic = new G4LogicalVolume(outer_solid,
  		Mouter,
			(name+"_outer").c_str(),
			nullptr, nullptr, g4userlimits);

  outer_logic->SetVisAttributes(vis_outer);

  new G4PVPlacement(
  		0,
			G4ThreeVector(0, 0, place),
			outer_logic,
			(name+"_outer").c_str(),
			mother, 0, false, overlapcheck);

  return true;
}

//_______________________________________________________________
void PHG4TargetCoilDetector::Construct(G4LogicalVolume *logicWorld)
{
	G4double z;
	G4double a;
	G4String symbol;
	G4String name;
	G4double density;
	G4int ncomponents;
	G4int natoms;

	G4Element *elHe = new G4Element(name="Helium", symbol="He" , z=2.,  a = 4.003*g/mole);
	G4Material* lHe = new G4Material(name = "G4_lHe", density = 0.145 * g/cm3, ncomponents = 1);
	lHe->AddElement(elHe, natoms = 1);

//	G4Element *elFe = new G4Element(name="Iron", symbol="Fe" , z=26., a = 55.845*g/mole);
//	G4Material* sFe = new G4Material(name = "G4_sFe", density = 7.87  * g/cm3, ncomponents = 1);
//	sFe->AddElement(elFe, natoms = 1);
//
//	G4Element *elNb = new G4Element(name="Niobium", symbol="Nb" , z=41., a = 92.906*g/mole);
//	G4Material* sNb = new G4Material(name = "G4_sNb", density = 8.57  * g/cm3, ncomponents = 1);
//	sNb->AddElement(elNb, natoms = 1);
//
//	G4Element *elTi = new G4Element(name="Titanium", symbol="Ti" , z=22., a = 47.867*g/mole);
//	G4Material* sTi = new G4Material(name = "G4_sTi", density = 4.506  * g/cm3, ncomponents = 1);
//	sTi->AddElement(elTi, natoms = 1);
//
//	G4Element *elCu = new G4Element(name="Copper", symbol="Cu" , z=29., a = 63.546*g/mole);
//	G4Material* sCu = new G4Material(name = "G4_sCu", density = 8.96  * g/cm3, ncomponents = 1);
//	sCu->AddElement(elCu, natoms = 1);
//
//	G4Element *elCr = new G4Element(name="Chromium", symbol="Cr" , z=24., a = 51.996*g/mole);
//	G4Material* sCr = new G4Material(name = "G4_sCr", density = 7.19*g/cm3, ncomponents = 1);
//	sCr->AddElement(elCr, natoms = 1);
//
//	G4Element *elNi = new G4Element(name="Nickel", symbol="Ni" , z=28., a = 28.6934*g/mole);
//	G4Material* sNi = new G4Material(name = "G4_sNi", density = 8.908*g/cm3, ncomponents = 1);
//	sNi->AddElement(elNi, natoms = 1);
//
//	G4Element *elMo = new G4Element(name="Molybdenum", symbol="Mo" , z=42., a = 95.95*g/mole);
//	G4Material* sMo = new G4Material(name = "G4_sMo", density = 10.28*g/cm3, ncomponents = 1);
//	sMo->AddElement(elMo, natoms = 1);


	// 1/(0.45/8.57+0.45/4.506+0.1/8.96) = 6.11 g/cm3
	G4Material* Coil = new G4Material(name = "Coil", density = 6.11*g/cm3, ncomponents = 3);
	Coil->AddMaterial(G4Material::GetMaterial("G4_Nb"), 45 * perCent);
	Coil->AddMaterial(G4Material::GetMaterial("G4_Ti"), 45 * perCent);
	Coil->AddMaterial(G4Material::GetMaterial("G4_Cu"), 10 * perCent);
	std::cout<< "DEBUG: " << Coil << std::endl;

	// 1/(0.6/7.87+0.2/7.18+0.15/8.902+0.05/10.22) = 7.95 g/cm3
	G4Material* SS316L = new G4Material(name = "SS316L", density = 7.95*g/cm3, ncomponents = 4);
	SS316L->AddMaterial(G4Material::GetMaterial("G4_Fe"), 60 * perCent);
	SS316L->AddMaterial(G4Material::GetMaterial("G4_Cr"), 20 * perCent);
	SS316L->AddMaterial(G4Material::GetMaterial("G4_Ni"), 15 * perCent);
	SS316L->AddMaterial(G4Material::GetMaterial("G4_Mo"),  5 * perCent);
	std::cout<< "DEBUG: " << SS316L << std::endl;


  G4VisAttributes *cylinder_vis = new G4VisAttributes();
	PHG4Utils::SetColour(cylinder_vis, "G4_He");
	cylinder_vis->SetVisibility(true);
	cylinder_vis->SetForceSolid(true);

  // determine length of cylinder using PHENIX's rapidity coverage if flag is true
  
  double l = 22.7 * cm; // length
  double ri = 6.0 * cm;
  double ro = 22.225 * cm;
  double ts = 0.3 * cm; // shell thickness

  //double cl = 1e-2 * cm;

  G4VSolid *cylinder_solid = new G4Tubs(G4String(GetName().c_str()),
  		ri,
			ro,
			l/2,
			0*deg,
			twopi);

  double steplimits = params->get_double_param("steplimits") * cm;
  G4UserLimits *g4userlimits = nullptr;
  if (isfinite(steplimits))
  {
    g4userlimits = new G4UserLimits(steplimits);
  }

  G4LogicalVolume *cylinder_logic = new G4LogicalVolume(cylinder_solid,
  		lHe,
  		G4String(GetName().c_str()),
			nullptr, nullptr, g4userlimits);

  cylinder_logic->SetVisAttributes(cylinder_vis);

  G4RotationMatrix *rotm  = new G4RotationMatrix();
  rotm->rotateX(params->get_double_param("rot_x")*deg);
  rotm->rotateY(params->get_double_param("rot_y")*deg);
  rotm->rotateZ(params->get_double_param("rot_z")*deg);
  params->Print();
  rotm->print(std::cout);
  cylinder_physi = new G4PVPlacement(
  		rotm,
  		G4ThreeVector(
  				params->get_double_param("place_x") * cm,
  				params->get_double_param("place_y") * cm,
					params->get_double_param("place_z") * cm),
					cylinder_logic,
					G4String(GetName().c_str()),
					logicWorld, 0, false, overlapcheck);

  PlaceHollowTube(
  		cylinder_logic,
			0,
			g4userlimits,
			overlapcheck,
			"Shell",
			SS316L,
			ri,
			ro,
			l/2,
			ts,
			0,
			twopi
  );

  double c1_l = 4.5 * cm;
  double c1_ri =  12.5 *cm;
  double c1_ro = 17.2 *cm;
  double c1_t = 0.5 *cm;

  PlaceLayeredTube(
  		cylinder_logic,
			-7.5 *cm,
			g4userlimits,
			overlapcheck,
			"C1",
			Coil,
			SS316L,
			c1_ri-c1_t,
			c1_ro+c1_t,
			c1_l/2+c1_t,
			c1_t,
			0,
			twopi
  );


  double c2_l = 5.7 * cm;
  double c2_ri =  7.6 *cm;
  double c2_ro = 9.4 *cm;
  double c2_t = 0.5 *cm;

  PlaceLayeredTube(
  		cylinder_logic,
			-2.8 *cm,
			g4userlimits,
			overlapcheck,
			"C2",
			Coil,
			SS316L,
			c2_ri-c2_t,
			c2_ro+c2_t,
			c2_l/2+c2_t,
			c2_t,
			0,
			twopi
  );


  double c3_l = 1 * cm;
  double c3_ri =  12.7 *cm;
  double c3_ro = 13.7 *cm;
  double c3_t = 0.5 *cm;

  PlaceLayeredTube(
  		cylinder_logic,
			0.9 *cm,
			g4userlimits,
			overlapcheck,
			"C3",
			Coil,
			SS316L,
			c3_ri-c3_t,
			c3_ro+c3_t,
			c3_l/2+c3_t,
			c3_t,
			0,
			twopi
  );
}







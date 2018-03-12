#include "PHG4CylinderDetector.h"

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
#include <Geant4/G4Tubs.hh>
#include <Geant4/G4UserLimits.hh>
#include <Geant4/G4VisAttributes.hh>

#include <cmath>
#include <sstream>

using namespace std;

//_______________________________________________________________
PHG4CylinderDetector::PHG4CylinderDetector(PHCompositeNode *Node, PHParameters *parameters, const std::string &dnam, const int lyr)
  : PHG4Detector(Node, dnam)
  , params(parameters)
  , cylinder_physi(nullptr)
  , layer(lyr)
{
}

//_______________________________________________________________
bool PHG4CylinderDetector::IsInCylinder(const G4VPhysicalVolume *volume) const
{
  if (volume == cylinder_physi)
  {
    return true;
  }
  return false;
}

//_______________________________________________________________
void PHG4CylinderDetector::Construct(G4LogicalVolume *logicWorld)
{
	G4Material *TrackerMaterial = nullptr;
	if (params->get_string_param("material").find("Target")
			!= std::string::npos) {
		G4double z;
		G4double a;
		G4String symbol;
		G4String name;
		G4double density;
		G4int ncomponents;
		G4int natoms;

		G4Element *elH = new G4Element(name="Hydrogen", symbol="H" , z=1., a = 1.01*g/mole);
		G4Element *elN = new G4Element(name="Nitrogen", symbol="N" , z=7., a = 14.0*g/mole);

		G4Material* N2 = new G4Material(name = "G4_N2", density = 0.25 * g/cm3, ncomponents = 1);
		N2->AddElement(elN, natoms = 2);

		G4Material* NH3 = new G4Material(name = "G4_NH3", density = 0.867 * g/cm3, ncomponents = 2);
		NH3->AddElement(elN, natoms = 1);
		NH3->AddElement(elH, natoms = 3);

		G4Material* Target = new G4Material(name = "Target", density = 0.59 * g/cm3, ncomponents = 2);
		Target->AddMaterial(NH3, 60 * perCent);
		Target->AddMaterial(N2,  40 * perCent);

		TrackerMaterial = Target;

		std::cout<< "DEBUG: " << TrackerMaterial << std::endl;
	} else {
		TrackerMaterial = G4Material::GetMaterial(
				params->get_string_param("material"));
	}

  if (!TrackerMaterial)
  {
    std::cout << "Error: Can not set material" << std::endl;
    exit(-1);
  }

  G4VisAttributes *siliconVis = new G4VisAttributes();
  if (params->get_int_param("blackhole"))
  {
    PHG4Utils::SetColour(siliconVis, "BlackHole");
    siliconVis->SetVisibility(false);
    siliconVis->SetForceSolid(false);
  }
  else
  {
    PHG4Utils::SetColour(siliconVis, params->get_string_param("material"));
    siliconVis->SetVisibility(true);
    siliconVis->SetForceSolid(true);
  }

  // determine length of cylinder using PHENIX's rapidity coverage if flag is true
  double radius = params->get_double_param("radius") * cm;
  double thickness = params->get_double_param("thickness") * cm;
  G4VSolid *cylinder_solid = new G4Tubs(G4String(GetName().c_str()),
                                        radius,
                                        radius + thickness,
                                        params->get_double_param("length") * cm / 2., 0, twopi);
  double steplimits = params->get_double_param("steplimits") * cm;
  G4UserLimits *g4userlimits = nullptr;
  if (isfinite(steplimits))
  {
    g4userlimits = new G4UserLimits(steplimits);
  }

  G4LogicalVolume *cylinder_logic = new G4LogicalVolume(cylinder_solid,
                                                        TrackerMaterial,
                                                        G4String(GetName().c_str()),
                                                        nullptr, nullptr, g4userlimits);
  cylinder_logic->SetVisAttributes(siliconVis);

  G4RotationMatrix *rotm  = new G4RotationMatrix();
  rotm->rotateX(params->get_double_param("rot_x")*deg);
  rotm->rotateY(params->get_double_param("rot_y")*deg);
  rotm->rotateZ(params->get_double_param("rot_z")*deg);
  params->Print();
  rotm->print(std::cout);
  cylinder_physi = new G4PVPlacement(rotm, G4ThreeVector(params->get_double_param("place_x") * cm,
                                                      params->get_double_param("place_y") * cm,
                                                      params->get_double_param("place_z") * cm),
                                     cylinder_logic,
                                     G4String(GetName().c_str()),
                                     logicWorld, 0, false, overlapcheck);
}

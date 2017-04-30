#ifndef PHG4SiliconTrackerDetector_h
#define PHG4SiliconTrackerDetector_h

#include <g4main/PHG4Detector.h>

#include <Geant4/G4RotationMatrix.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4Types.hh>
#include <Geant4/globals.hh>

#include <set>
#include <utility>
#include <vector>

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4VSolid;
class PHG4ParametersContainer;

class PHG4SiliconTrackerDetector : public PHG4Detector
{
 public:
  typedef std::vector<std::pair<int, int>> vpair;

  //! constructor
  PHG4SiliconTrackerDetector(PHCompositeNode *Node, PHG4ParametersContainer *parameters, const std::string &dnam = "SILICON_TRACKER", const vpair &layerconfig = vpair(0));

  //! destructor
  virtual ~PHG4SiliconTrackerDetector();

  //! construct
  virtual void Construct(G4LogicalVolume *world);

  //!@name volume accessors
  //@{
  int IsInSiliconTracker(G4VPhysicalVolume *) const;
  //@}

  void SuperDetector(const std::string &name)
  {
    superdetector = name;
  }
  const std::string SuperDetector() const
  {
    return superdetector;
  }
  void Detector(const std::string &name)
  {
    detector_type = name;
  }
  const std::string Detector() const
  {
    return detector_type;
  }

  const int arr_nladders_layer[4] = {20, 26, 32, 38};
  const G4double arr_radius[4] = {60.0 * mm, 80.0 * mm, 100.0 * mm, 120.0 * mm};
  const G4double arr_offsetphi[4] = {0.0 / 180. * CLHEP::pi, 0.0 / 180. * CLHEP::pi, 0.0 / 180. * CLHEP::pi, 0.0 / 180. * CLHEP::pi};
  const G4double arr_offsetrot[4] = {14.0 / 180. * CLHEP::pi, 14.0 / 180. * CLHEP::pi, 12.0 / 180. * CLHEP::pi, 11.5 / 180. * CLHEP::pi};
  const int arr_nstrips_z_sensor[2][2] = {/*Layer0*/ {5, 5}, /*Layer1-3*/ {8, 5}};
  const int arr_nstrips_phi_cell[4] = {128, 128, 128, 128};
  const G4double arr_strip_x[4] = {0.200 * mm * 0.5, 0.200 * mm * 0.5, 0.200 * mm * 0.5, 0.200 * mm * 0.5};  // 200 micron
  //const G4double arr_strip_x[4] = {0.240*mm*0.5, 0.240*mm*0.5, 0.240*mm*0.5, 0.240*mm*0.5}; // 240 micron
  //const G4double arr_strip_x[4] = {0.320*mm*0.5, 0.320*mm*0.5, 0.320*mm*0.5, 0.320*mm*0.5}; // 320 micron
  const G4double arr_strip_y[4] = {0.078 * mm * 0.5, 0.086 * mm * 0.5, 0.086 * mm * 0.5, 0.086 * mm * 0.5};
  const G4double arr_strip_z[2][2] = {/*Layer0*/ {18.0 * mm * 0.5, 18.0 * mm * 0.5}, /*Layer1-3*/ {16.0 * mm * 0.5, 20.0 * mm * 0.5}};

 private:
  void AddGeometryNode();
  int ConstructSiliconTracker(G4LogicalVolume *sandwich);
  int DisplayVolume(G4VSolid *volume, G4LogicalVolume *logvol, G4RotationMatrix *rotm = NULL);

  PHG4ParametersContainer *paramscontainer;
  int active;
  int absorberactive;
  int blackhole;

  vpair layerconfig_;
  unsigned int nlayer_;
  int layermin_;
  int layermax_;

  G4double overlap_fraction;

  std::string detector_type;
  std::string superdetector;

  const G4double sensor_edge_phi = 1.305 * mm;
  const G4double sensor_edge_z = 0.98 * mm;
  const G4double hdi_edge_z = 0.10 * mm;

  G4double eff_radius[4];
  G4double posz[4][2];
  G4double strip_x_offset[4];

  //const double hdi_x  = 0.473*mm * 0.5;
  const double hdi_x = 0.38626 * mm * 0.5;  // updated design
  const double fphx_x = 0.32 * mm * 0.5;
  const double fphx_y = 2.7 * mm * 0.5;
  const double fphx_z = 9.0 * mm * 0.5;

  //const double pgs_x = 0.35*mm * 0.5; // 70micron * 5layers
  const double pgs_x = 0.21 * mm * 0.5;    // 70micron * 3layers
  const double stave_x = 0.23 * mm * 0.5;  // too thick? TODO

  const G4double arr_hdi_y[4] = {38. * mm * 0.5, 43. * mm * 0.5, 43. * mm * 0.5, 43. * mm * 0.5};
  const G4double arr_halfladder_z[4] = {220. * mm * 0.5, 268. * mm * 0.5, 268. * mm * 0.5, 268. * mm * 0.5};
  std::set<G4LogicalVolume *> absorberlogvols;
  std::set<G4LogicalVolume *> activelogvols;
};

#endif

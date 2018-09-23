#include "PHG4TPCElectronDrift.h"
#include "PHG4TPCPadPlaneReadout.h"
//#include "PHG4TPCPadPlane.h"

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>

#include <g4detectors/PHG4Cellv1.h>
#include <g4detectors/PHG4CellContainer.h>
#include <g4detectors/PHG4CylinderCellGeom.h>
#include <g4detectors/PHG4CylinderCellGeomContainer.h>

#include <phparameter/PHParametersContainer.h>

#include <pdbcalbase/PdbParameterMapContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>

#include <phool/PHNodeIterator.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>

#include <TSystem.h>
#include <TH1.h>
#include <TFile.h>
#include <TNtuple.h>

#include <Geant4/G4SystemOfUnits.hh>

#include <gsl/gsl_randist.h>
#include <iostream>

using namespace std;

PHG4TPCElectronDrift::PHG4TPCElectronDrift(const std::string& name):
  SubsysReco(name),
  PHParameterInterface(name),
  g4cells(nullptr),
  dlong(nullptr),
  dtrans(nullptr),
  diffusion_trans(NAN),
  diffusion_long(NAN),
  drift_velocity(NAN),
  electrons_per_gev(NAN),
  min_active_radius(NAN),
  max_active_radius(NAN),
  min_time(NAN),
  max_time(NAN)
{  
  cout << "Constructor of PHG4TPCElectronDrift" << endl;
  InitializeParameters();
  RandomGenerator = gsl_rng_alloc(gsl_rng_mt19937);
  set_seed(PHRandomSeed()); // fixed seed is handled in this funtcion
  return;
}

PHG4TPCElectronDrift::~PHG4TPCElectronDrift()
{
  gsl_rng_free(RandomGenerator);
}

int PHG4TPCElectronDrift::Init(PHCompositeNode *topNode)
{
  padplane->Init(topNode);

  return Fun4AllReturnCodes::EVENT_OK;
}


int PHG4TPCElectronDrift::InitRun(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
    {
      std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
      exit(1);
    }
  PHCompositeNode *runNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "RUN" ));
  PHCompositeNode *parNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "PAR" ));
  string paramnodename = "G4CELLPARAM_" + detector;
  string geonodename = "G4CELLPAR_" + detector;
  string tpcgeonodename = "G4GEO_" + detector;
  hitnodename = "G4HIT_" + detector;
  PHG4HitContainer *g4hit = findNode::getClass<PHG4HitContainer>(topNode, hitnodename.c_str());
  if (!g4hit)
    {
      cout << Name() << " Could not locate G4HIT node " << hitnodename << endl;
      topNode->print();
      gSystem->Exit(1);
      exit(1);
    }
  cellnodename = "G4CELL_SVTX";  // + detector;
  g4cells = findNode::getClass<PHG4CellContainer>(topNode,cellnodename);
  if (! g4cells)
  {
   PHNodeIterator dstiter(dstNode);
    PHCompositeNode *DetNode = dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", detector));

    if (!DetNode)
    {
      DetNode = new PHCompositeNode(detector);
      dstNode->addNode(DetNode);
    }
    g4cells = new PHG4CellContainer();
    PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(g4cells, cellnodename.c_str(), "PHObject");
    DetNode->addNode(newNode);
  }

  seggeonodename = "CYLINDERCELLGEOM_SVTX"; // + detector;
  PHG4CylinderCellGeomContainer *seggeo = findNode::getClass<PHG4CylinderCellGeomContainer>(topNode, seggeonodename.c_str());
  if (!seggeo)
  {
    seggeo = new PHG4CylinderCellGeomContainer();
    PHCompositeNode *runNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "RUN"));
    PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(seggeo, seggeonodename.c_str(), "PHObject");
    runNode->addNode(newNode);
  }

   UpdateParametersWithMacro();
  PHNodeIterator runIter(runNode);
  PHCompositeNode *RunDetNode =  dynamic_cast<PHCompositeNode*>(runIter.findFirst("PHCompositeNode",detector));
  if (! RunDetNode)
    {
      RunDetNode = new PHCompositeNode(detector);
      runNode->addNode(RunDetNode);
    }
  SaveToNodeTree(RunDetNode,paramnodename);

  // save this to the parNode for use
  PHNodeIterator parIter(parNode);
  PHCompositeNode *ParDetNode =  dynamic_cast<PHCompositeNode*>(parIter.findFirst("PHCompositeNode",detector));
  if (! ParDetNode)
    {
      ParDetNode = new PHCompositeNode(detector);
      parNode->addNode(ParDetNode);
    }
  PutOnParNode(ParDetNode,geonodename);

  diffusion_long = get_double_param("diffusion_long");
  diffusion_trans = get_double_param("diffusion_trans");
  drift_velocity = get_double_param("drift_velocity");
  electrons_per_gev = get_double_param("electrons_per_gev");
  min_active_radius = get_double_param("min_active_radius");
  max_active_radius = get_double_param("max_active_radius");

// find TPC Geo
  PHNodeIterator tpcpariter(ParDetNode);
  PHParametersContainer *tpcparams = findNode::getClass<PHParametersContainer>(ParDetNode,tpcgeonodename);
  if (!tpcparams)
  {
    string runparamname = "G4GEOPARAM_" + detector;
    PdbParameterMapContainer *tpcpdbparams = findNode::getClass<PdbParameterMapContainer>(RunDetNode,runparamname);
    if (tpcpdbparams)
    {
      tpcparams = new PHParametersContainer(detector);
      tpcpdbparams->print();
      tpcparams->CreateAndFillFrom(tpcpdbparams,detector);
      ParDetNode->addNode(new PHDataNode<PHParametersContainer>(tpcparams,tpcgeonodename));
    }
  }
  tpcparams->Print();

  Fun4AllServer *se = Fun4AllServer::instance();
  dlong = new TH1F("difflong","longitudinal diffusion",100,diffusion_long-diffusion_long/2.,diffusion_long+diffusion_long/2.);
  se->registerHisto(dlong);
  dtrans = new TH1F("difftrans","transversal diffusion",100,diffusion_trans-diffusion_trans/2.,diffusion_trans+diffusion_trans/2.);
  se->registerHisto(dtrans);
  nt = new TNtuple("nt","stuff","hit:ts:tb:tsig:rad:zstart:zfinal");
  nthit = new TNtuple("nthit","hit stuff","hit:layer:phi:phicenter:z_gem:zcenter:weight");
  ntpad = new TNtuple("ntpad","padplane stuff","layer:phi:phiclus:z:zclus");
  se->registerHisto(nt);
  se->registerHisto(nthit);
  se->registerHisto(ntpad);
  padplane->InitRun(topNode);
  padplane->CreateReadoutGeometry(topNode,seggeo);
 
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4TPCElectronDrift::process_event(PHCompositeNode *topNode)
{
  //int verbosity = 500;

  PHG4HitContainer *g4hit = findNode::getClass<PHG4HitContainer>(topNode, hitnodename.c_str());
  if (!g4hit)
    {
      cout << "Could not locate g4 hit node " << hitnodename << endl;
      gSystem->Exit(1);
    }

  PHG4HitContainer::ConstIterator hiter;
  PHG4HitContainer::ConstRange hit_begin_end = g4hit->getHits();
  double tpc_length = 211.;
  double ihit = 0;
  for (hiter = hit_begin_end.first; hiter != hit_begin_end.second; ++hiter)
  {
    double t0 = fmax(hiter->second->get_t(0),hiter->second->get_t(1));
    if (t0 > max_time)
    {
      continue;
    }
    double eion = hiter->second->get_eion();
    unsigned int n_electrons = gsl_ran_poisson(RandomGenerator,eion*electrons_per_gev);
    if(verbosity > 100) 
      cout << "  new hit with t0, " <<  t0 << " g4hitid " << hiter->first 
	   << " eion " << eion << " n_electrons " << n_electrons 
	   << " entry z " << hiter->second->get_z(0) << " exit z " << hiter->second->get_z(1) << " avg z" << (hiter->second->get_z(0) + hiter->second->get_z(1))/2.0 
	   << endl;

    //nthit->Fill(ihit,n_electrons,eion,hiter->second->get_edep(),hiter->second->get_t(0),hiter->second->get_x(0),hiter->second->get_y(0),hiter->second->get_z(0));
    if (n_electrons <= 0)
    {
      if (n_electrons < 0)
      {
	cout << "really bad number of electrons: " << n_electrons
	     << ", eion: " << eion
	     << endl;
      }
      continue;
    }
    double dx = (hiter->second->get_x(1) - hiter->second->get_x(0))/n_electrons;
    double dy = (hiter->second->get_y(1) - hiter->second->get_y(0))/n_electrons;
    double dz = (hiter->second->get_z(1) - hiter->second->get_z(0))/n_electrons;
    double dt = (hiter->second->get_t(1) - hiter->second->get_t(0))/n_electrons;
    double x_start = hiter->second->get_x(0) + dx/2.;
    double y_start = hiter->second->get_y(0) + dy/2.;
    double z_start = hiter->second->get_z(0) + dz/2.;
    double t_start = hiter->second->get_t(0) + dt/2.;
    // cout << "g4hit created electrons: " << n_electrons 
    // 	   << " from " << eion*1000000 << " keV" << endl;
    
    for (unsigned int i=0; i<n_electrons; i++)
      {
	double radstart = sqrt(x_start*x_start + y_start*y_start);
	double r_sigma =  diffusion_trans*sqrt(tpc_length/2. - fabs(z_start));
	double rantrans = gsl_ran_gaussian(RandomGenerator,r_sigma);

	double t_path = (tpc_length/2. - fabs(z_start))/drift_velocity;
	// now the drift
	double t_sigma =  diffusion_long*sqrt(tpc_length/2. - fabs(z_start))/drift_velocity;	
	double rantime = gsl_ran_gaussian(RandomGenerator,t_sigma);

	// ADF: note that adding t_start to the path results in an error of t_start * drift_velocity in the final position 
	// The intention here seems to be to account for the propagation time of the particle to this point in the TPC
	// But the absolute time at this point in the TPC is not the relevant quantity. It is the time offset relative to the clock. 
	// We presumably calibrate the TPC to get the average value of z correct, so we should be using 1/2 of the time spread of hits in the TPC
	// which is about +/- 25% of the absolute time. So ignore t_start for now to make it easier to compare with truth information.
	// double t_final = t_start + t_path + rantime;
	double t_final = t_path + rantime;

	double z_final;
	if(z_start < 0)
	  z_final = -tpc_length/2. + t_final * drift_velocity;
	else
	  z_final = tpc_length/2. - t_final * drift_velocity;

	if (t_final < min_time || t_final > max_time)
	  {
	    cout << "skip this, t_final out of range" << endl;
	    continue;
	  }
	double ranphi = gsl_ran_flat(RandomGenerator,-M_PI,M_PI);
	double x_final = x_start + rantrans*cos(ranphi);
	double y_final = y_start + rantrans*sin(ranphi);
	double rad_final = sqrt(x_final*x_final + y_final*y_final);
	// remove electrons outside of our acceptance
	if (rad_final<min_active_radius || rad_final >max_active_radius)
	  {
	    //cout << " skip - outside active area of TPC" << endl;
	    continue;
	  }

	//if(radstart > 70.0 && radstart < 71.0)
	if(verbosity > 100)
	  {
	    if(verbosity > 1000)
	      {
		cout << "electron " << i << " g4hitid " << hiter->first << endl; 
		cout << "radstart " << radstart  << " x_start: " << x_start
		     << ", y_start: " << y_start
		     << ",z_start: " << z_start 
		     << " t_start " << t_start
		     << " t_path " << t_path
		     << " t_sigma " << t_sigma
		     << " rantime " << rantime
		     << endl;
	      }
	    //cout << "rad_final " << rad_final << " x_final " << x_final << " y_final " << y_final 
	    //	 << " z_final " << z_final << " t_final " << t_final << " zdiff " << z_final - z_start << endl; 

	    cout << "rad_final " << rad_final << " z_start " << z_start 
		 << " z_final " << z_final << " zdiff " << z_final - z_start << endl; 
	  }

	nt->Fill(ihit,t_start,t_final,t_sigma,rad_final,z_start,z_final);

	// this fills the cells and updates them on the node tree for this drifted electron hitting the GEM stack
	MapToPadPlane(x_final,y_final,z_final, hiter,ntpad,nthit);
	x_start += dx;
	y_start += dy;
	z_start += dz;
	t_start += dt;
      }
    ihit++;
  }
  if (Verbosity()>1)
    {
      cout << " Cells associated with this hit:" << endl; 
      {
	PHG4CellContainer::ConstRange cells = g4cells->getCells();
	PHG4CellContainer::ConstIterator celliter;
	for (celliter=cells.first;celliter != cells.second; ++celliter)
	  {
	    //celliter->second->identify();
	    if(celliter->second->get_layer() == 50)
	      {
		int phibin = PHG4CellDefs::SizeBinning::get_phibin(celliter->second->get_cellid());//cell->get_binphi();
		int zbin = PHG4CellDefs::SizeBinning::get_zbin(celliter->second->get_cellid());//cell->get_binz();
		cout << " cellid " << celliter->second->get_cellid() 
		     << " layer " << celliter->second->get_layer()
		     << " zbin " << zbin
		     << " phibin " << phibin
		     << " edep " << celliter->second->get_edep()
		     << endl;
	      }
	  }
      }
    }

  return Fun4AllReturnCodes::EVENT_OK;
}

void PHG4TPCElectronDrift::MapToPadPlane(const double x_gem, const double y_gem, const double t_gem, PHG4HitContainer::ConstIterator hiter, TNtuple *ntpad, TNtuple *nthit)
{

  padplane->MapToPadPlane(g4cells,x_gem,y_gem,t_gem, hiter,ntpad,nthit);
 
  return;
}

int PHG4TPCElectronDrift::End(PHCompositeNode *topNode)
{
  TFile *outf=new TFile("nt_out.root","recreate");
  outf->WriteTObject(nt);
  outf->WriteTObject(ntpad);
  outf->WriteTObject(nthit);
  outf->Close();

  return Fun4AllReturnCodes::EVENT_OK;
}

void PHG4TPCElectronDrift::set_seed(const unsigned int iseed)
{
  seed = iseed;
  gsl_rng_set(RandomGenerator,seed);
}

void PHG4TPCElectronDrift::SetDefaultParameters()
{
// Data on gasses @20 C and 760 Torr from the following source:
// http://www.slac.stanford.edu/pubs/icfa/summer98/paper3/paper3.pdf
// diffusion and drift velocity for 400kV for NeCF4 90/10 from calculations:
// https://www.phenix.bnl.gov/WWW/p/draft/prakhar/tpc/HTML_Gas_Linear/Ne_CF4_90_10.html
double Ne_dEdx = 1.56;    // keV/cm
double CF4_dEdx = 7.00;   // keV/cm
// double Ne_NPrimary = 12;    // Number/cm
// double CF4_NPrimary = 51;   // Number/cm
double Ne_NTotal = 43;     // Number/cm
double CF4_NTotal = 100;   // Number/cm
double TPC_NTot = 0.9*Ne_NTotal + 0.1*CF4_NTotal;
double TPC_dEdx = 0.90 * Ne_dEdx + 0.10 * CF4_dEdx;
double  TPC_ElectronsPerKeV = TPC_NTot / TPC_dEdx;
  set_default_double_param("diffusion_long",0.015); // cm/SQRT(cm)
  set_default_double_param("diffusion_trans",0.006); // cm/SQRT(cm)
  set_default_double_param("drift_velocity",8.0 / 1000.0); // cm/ns
  set_default_double_param("electrons_per_gev",TPC_ElectronsPerKeV*1000000.);
  set_default_double_param("min_active_radius",30.); // cm
  set_default_double_param("max_active_radius",78.); // cm
  set_default_double_param("min_time",0.); // ns
  set_default_double_param("max_time",14000.); // ns

  return;
}

void PHG4TPCElectronDrift::registerPadPlane(PHG4TPCPadPlane *inpadplane)
{
  cout << "Registering padplane " << endl;
  padplane = inpadplane;
  padplane->Detector(Detector());
  padplane->UpdateInternalParameters();
  cout << "padplane registered and parameters updated" << endl;

return;
}

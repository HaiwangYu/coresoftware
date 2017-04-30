#include "PHG4DetectorGroupSubsystem.h"
#include "PHG4Parameters.h"
#include "PHG4ParametersContainer.h"

#include <pdbcalbase/PdbParameterMap.h>
#include <pdbcalbase/PdbParameterMapContainer.h>

#include <phool/getClass.h>
#include <phool/phool.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>

#include <boost/format.hpp>

#include <iostream>
#include <sstream>

using namespace std;

PHG4DetectorGroupSubsystem::PHG4DetectorGroupSubsystem(const std::string &name, const int lyr): 
  PHG4Subsystem(name),
  params(new PHG4Parameters(Name())),
  paramscontainer(new PHG4ParametersContainer(Name())),
  savetopNode(NULL),
  overlapcheck(false),
  layer(lyr),
  usedb(0),
  beginrunexecuted(0),
  filetype(PHG4DetectorGroupSubsystem::none),
  superdetector("NONE"),
  calibfiledir("./")
{
  // put the layer into the name so we get unique names
  // for multiple layers
  ostringstream nam;
  nam << name << "_" << lyr;
  Name(nam.str().c_str());
  paramsmap[lyr] = new PHG4Parameters(Name());
}

int 
PHG4DetectorGroupSubsystem::Init(PHCompositeNode* topNode)
{

  savetopNode = topNode;
  cout << "setting param name to " << Name() << endl;
  params->set_name(Name());
  int iret = InitSubsystem(topNode);
  return iret;
}

int 
PHG4DetectorGroupSubsystem::InitRun( PHCompositeNode* topNode )
{
  PHNodeIterator iter( topNode );
  PHCompositeNode *parNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "PAR" ));
  PHCompositeNode *runNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "RUN" ));

  string g4geonodename = "G4GEO_";
  string paramnodename = "G4GEOPARAM_";
  string calibdetname;
  int isSuperDetector = 0;
  if (superdetector != "NONE")
    {
      g4geonodename += SuperDetector();
      paramscontainer = findNode::getClass<PHG4ParametersContainer>(parNode,g4geonodename);
      if (! paramscontainer)
	{
	  PHNodeIterator parIter(parNode);
          PHCompositeNode *DetNode = dynamic_cast<PHCompositeNode*>(parIter.findFirst("PHCompositeNode",SuperDetector()));
	  if (! DetNode)
	    {
	      DetNode = new PHCompositeNode(SuperDetector());
	      parNode->addNode(DetNode);
	    }
	  paramscontainer = new PHG4ParametersContainer(superdetector);
	  DetNode->addNode(new PHDataNode<PHG4ParametersContainer>(paramscontainer,g4geonodename));
	}
      paramscontainer->AddPHG4Parameters(layer,params);
      paramnodename += superdetector;
      calibdetname = superdetector;
      isSuperDetector = 1;
    }
  else
    {
      g4geonodename += params->Name();
      parNode->addNode(new PHDataNode<PHG4Parameters>(params,g4geonodename));
      paramnodename += params->Name();
      calibdetname = params->Name();
    }



  // ASSUMPTION: if we read from DB and/or file we don't want the stuff from
  // the node tree
  // We leave the defaults intact in case there is no entry for
  // those in the object read from the DB or file
  // Order: read first DB, then calib file if both are enabled
  if (ReadDB() || get_filetype() != PHG4DetectorGroupSubsystem::none)
    {
      if (ReadDB())
	{
	   ReadParamsFromDB(calibdetname,isSuperDetector);
	}
      if (get_filetype() != PHG4DetectorGroupSubsystem::none)
	{
	  ReadParamsFromFile(calibdetname, get_filetype(),isSuperDetector );
	}
    }
  else
    {
      PdbParameterMapContainer *nodeparams = findNode::getClass<PdbParameterMapContainer>(topNode,paramnodename);
      if (nodeparams)
	{
	  params->FillFrom(nodeparams, layer);
	}
    }
  // parameters set in the macro always override whatever is read from
  // the node tree, DB or file
  UpdateParametersWithMacro();
  // save updated persistant copy on node tree
  PHCompositeNode *RunDetNode = runNode;
  if (superdetector != "NONE")
    {
      PHNodeIterator runIter(runNode);
      RunDetNode = dynamic_cast<PHCompositeNode*>(runIter.findFirst("PHCompositeNode",SuperDetector()));
      if (! RunDetNode)
	{
	  RunDetNode = new PHCompositeNode(SuperDetector());
	  runNode->addNode(RunDetNode);
	}
    }
  params->SaveToNodeTree(RunDetNode,paramnodename,layer);
  int iret = InitRunSubsystem(topNode);
  if (Verbosity() > 0)
    {
      PdbParameterMapContainer *nodeparams = findNode::getClass<PdbParameterMapContainer>(topNode,paramnodename);
      cout << Name() << endl;
      nodeparams->print();
    }
  beginrunexecuted = 1;
  return iret;
}

void
PHG4DetectorGroupSubsystem::SuperDetector(const std::string &name)
{
  superdetector = name;
  return;
}

void
PHG4DetectorGroupSubsystem::set_double_param(const int detid, const std::string &name, const double dval)
{
  map<int, map<const std::string, double>>::const_iterator iter = default_double.find(detid);
  if (iter == default_double.end())
  {
    cout << "detid " << detid << " not implemented" << endl;
    cout << "implemented detector ids: " << endl;
    for (map<int, map<const std::string, double>>::const_iterator iter2 = default_double.begin(); iter2 != default_double.end(); ++iter2)
    {
      cout << "detid: " << iter2->first << endl;
    }
    return;
  }
  if (iter->second.find(name) == iter->second.end())
  {
      cout << "double parameter " << name << " not implemented" << endl;
      cout << "implemented double parameters are:" << endl;
      for (map<const string, double>::const_iterator iter3 = iter->second.begin(); iter3 != iter->second.end(); ++iter3)
      {
	cout << iter3->first << endl;
      }
      return;
  }
// here we know we have entries for the detector id and the variable name exists
// in the defaults, so now lets set it
//  iter->second[name] = dval;
  return;
}
/*
  if (default_double.find(name) == default_double.end())
    {
      cout << "double parameter " << name << " not implemented" << endl;
      cout << "implemented double parameters are:" << endl;
      for (map<const string, double>::const_iterator iter = default_double.begin(); iter != default_double.end(); ++iter)
	{
	  cout << iter->first << endl;
	}
      return;
    }
  dparams[name] = dval;
*/


double
PHG4DetectorGroupSubsystem::get_double_param(const int detid, const std::string &name) const
{
  return params->get_double_param(name);
}

void
PHG4DetectorGroupSubsystem::set_int_param(const std::string &name, const int ival)
{
/*
  if (default_int.find(name) == default_int.end())
    {
      cout << "integer parameter " << name << " not implemented" << endl;
      cout << "implemented integer parameters are:" << endl;
      for (map<const string, int>::const_iterator iter = default_int.begin(); iter != default_int.end(); ++iter)
	{
	  cout << iter->first << endl;
	}
      return;
    }
  iparams[name] = ival;
*/
}

int
PHG4DetectorGroupSubsystem::get_int_param(const std::string &name) const
{
  return params->get_int_param(name);
}

void
PHG4DetectorGroupSubsystem::set_string_param(const std::string &name, const string &sval)
{
/*
  if (default_string.find(name) == default_string.end())
    {
      cout << "string parameter " << name << " not implemented" << endl;
      cout << "implemented string parameters are:" << endl;
      for (map<const string, string>::const_iterator iter = default_string.begin(); iter != default_string.end(); ++iter)
	{
	  cout << iter->first << endl;
	}
      return;
    }
  cparams[name] = sval;
*/
}

string
PHG4DetectorGroupSubsystem::get_string_param(const std::string &name) const
{
  return params->get_string_param(name);
}

void
PHG4DetectorGroupSubsystem::UpdateParametersWithMacro()
{
  for (map<const string,double>::const_iterator iter = dparams.begin(); iter != dparams.end(); ++iter)
    {
      params->set_double_param(iter->first,iter->second);
    }
  for (map<const string,int>::const_iterator iter = iparams.begin(); iter != iparams.end(); ++iter)
    {
      params->set_int_param(iter->first,iter->second);
    }
  for (map<const string,string>::const_iterator iter = cparams.begin(); iter != cparams.end(); ++iter)
    {
      params->set_string_param(iter->first,iter->second);
    }
  return;
}

void
PHG4DetectorGroupSubsystem::set_default_double_param(const int detid, const std::string &name, const double dval)
{
  map<int, map<const std::string, double>>::iterator dmapiter = default_double.find(detid);
  if (dmapiter == default_double.end())
  {
    map<const std::string, double> newdmap;
    newdmap[name] = dval;
    default_double[detid] = newdmap;
  }
  else
  {
    if (dmapiter->second.find(name) != dmapiter->second.end())
    {
      cout << "trying to overwrite default double " << name << " " 
	   << dmapiter->second.find(name)->second << " with " << dval << endl;
      exit(1);
    }
    else
    {
      dmapiter->second[name] = dval;
    }
  }
  return;
}

void
PHG4DetectorGroupSubsystem::set_default_int_param(const int detid, const std::string &name, const int ival)
{
  map<int, map<const std::string, int>>::iterator imapiter = default_int.find(detid);
  if (imapiter == default_int.end())
  {
    map<const std::string, int> newimap;
    newimap[name] = ival;
    default_int[detid] = newimap;
  }
  else
  {
    if (imapiter->second.find(name) != imapiter->second.end())
    {
      cout << "trying to overwrite default int " << name << " " 
	   << imapiter->second.find(name)->second << " with " << ival << endl;
      exit(1);
    }
    else
    {
      imapiter->second[name] = ival;
    }
  }
  return;
}

void
PHG4DetectorGroupSubsystem::set_default_string_param(const int detid, const std::string &name, const string &sval)
{
  map<int, map<const std::string, string>>::iterator smapiter = default_string.find(detid);
  if (smapiter == default_string.end())
  {
    map<const std::string, string> newsmap;
    newsmap[name] = sval;
    default_string[detid] = newsmap;
  }
  else
  {
    if (smapiter->second.find(name) != smapiter->second.end())
    {
      cout << "trying to overwrite default string " << name << " " 
	   << smapiter->second.find(name)->second << " with " << sval << endl;
      exit(1);
    }
    else
    {
      smapiter->second[name] = sval;
    }
  }
  return;
}

void
PHG4DetectorGroupSubsystem::InitializeParameters()
{
  for (set<int>::const_iterator iter=layers.begin(); iter != layers.end(); ++iter)
  {
    set_default_int_param(*iter, "absorberactive", 0);
  set_default_int_param(*iter, "absorbertruth", 0);
  set_default_int_param(*iter, "active", 0);
  set_default_int_param(*iter, "blackhole", 0);
  }
  SetDefaultParameters(); // call method from specific subsystem
  // now load those parameters to our params class
  map<int, map<const string, double>>::const_iterator diter;
  for (diter = default_double.begin(); diter != default_double.end(); ++diter)
  {
    PHG4Parameters *detidparams = paramscontainer->GetParametersToModify(diter->first);
    if (!detidparams)
    {
      detidparams = new PHG4Parameters(boost::str(boost::format("%s_%d") %Name() %diter->first));
      paramscontainer->AddPHG4Parameters(diter->first,detidparams);
    }
    map<const string, double>::const_iterator diter2;
    for (diter2 = diter->second.begin(); diter2 != diter->second.end(); ++diter2)
    {
      detidparams->set_double_param(diter2->first, diter2->second);
    }
  }

  map<int, map<const string, int>>::const_iterator iiter;
  for (iiter = default_int.begin(); iiter != default_int.end(); ++iiter)
  {
    PHG4Parameters *detidparams = paramscontainer->GetParametersToModify(iiter->first);
    if (!detidparams)
    {
      detidparams = new PHG4Parameters(boost::str(boost::format("%s_%d") %Name() %iiter->first));
      paramscontainer->AddPHG4Parameters(iiter->first,detidparams);
    }
    map<const string, int>::const_iterator iiter2;
    for (iiter2 = iiter->second.begin(); iiter2 != iiter->second.end(); ++iiter2)
    {
      detidparams->set_int_param(iiter2->first, iiter2->second);
    }
  }

  map<int, map<const string, string>>::const_iterator siter;
  for (siter = default_string.begin(); siter != default_string.end(); ++siter)
  {
    PHG4Parameters *detidparams = paramscontainer->GetParametersToModify(siter->first);
    if (!detidparams)
    {
      detidparams = new PHG4Parameters(boost::str(boost::format("%s_%d") %Name() %siter->first));
      paramscontainer->AddPHG4Parameters(siter->first,detidparams);
    }
    map<const string, string>::const_iterator siter2;
    for (siter2 = siter->second.begin(); siter2 != siter->second.end(); ++siter2)
    {
      detidparams->set_string_param(siter2->first, siter2->second);
    }
  }
}

int
PHG4DetectorGroupSubsystem::SaveParamsToDB()
{
  int iret = 0;
  if (paramscontainer)
    {
      iret = paramscontainer->WriteToDB();
    }
  else
    {
      iret = params->WriteToDB();
    }
  if (iret)
    {
      cout << "problem committing to DB" << endl;
    }
  return iret;
}

int
PHG4DetectorGroupSubsystem::ReadParamsFromDB(const string &name, const int issuper)
{
  int iret = 0;
  if (issuper)
    {
      iret = params->ReadFromDB(name,layer);
    }
  else
    {
      iret = params->ReadFromDB();
    }
  if (iret)
    {
      cout << "problem reading from DB" << endl;
    }
  return iret;
}

int
PHG4DetectorGroupSubsystem::SaveParamsToFile(const PHG4DetectorGroupSubsystem::FILE_TYPE ftyp)
{
  string extension;
  switch(ftyp)
    {
    case xml:
      extension = "xml";
      break;
    case root:
      extension = "root";
      break;
    default:
      cout << PHWHERE << "filetype " << ftyp << " not implemented" << endl;
      exit(1);
    }
  int iret = 0;
  if (paramscontainer)
    {
      iret = paramscontainer->WriteToFile(extension,calibfiledir);
    }
  else
    {
      iret = params->WriteToFile(extension,calibfiledir);
    }
  if (iret)
    {
      cout << "problem saving to " << extension << " file " << endl;
    }
  return iret;
}

int
PHG4DetectorGroupSubsystem::ReadParamsFromFile(const string &name, const PHG4DetectorGroupSubsystem::FILE_TYPE ftyp, const int issuper)
{
  string extension;
  switch(ftyp)
    {
    case xml:
      extension = "xml";
      break;
    case root:
      extension = "root";
      break;
    default:
      cout << PHWHERE << "filetype " << ftyp << " not implemented" << endl;
      exit(1);
    }
  int iret = params->ReadFromFile(name, extension, layer, issuper, calibfiledir);
  if (iret)
    {
      cout << "problem reading from " << extension << " file " << endl;
    }
  return iret;
}

void
PHG4DetectorGroupSubsystem::SetActive(const int detid, const int i)
{
  iparams["active"] = i;
}

void
PHG4DetectorGroupSubsystem::SetAbsorberActive(const int detid, const int i)
{
  iparams["absorberactive"] = i;
}

void
PHG4DetectorGroupSubsystem::BlackHole(const int detid, const int i)
{
  iparams["blackhole"] = i;
}

void
PHG4DetectorGroupSubsystem::SetAbsorberTruth(const int detid, const int i)
{
  iparams["absorbertruth"] = i;
}

void PHG4DetectorGroupSubsystem::PrintDefaultParams() const
{
  cout << "int values: " << endl;
  map<int, map<const std::string, int>>::const_iterator iiter;
  for (iiter = default_int.begin(); iiter != default_int.end(); ++iiter)
  {
    cout << "Detector id: " << iiter->first << endl;
    map<const string, int>::const_iterator iiter2;
    for (iiter2 = iiter->second.begin(); iiter2 != iiter->second.end(); ++iiter2)
    {
      cout << iiter2->first << ": " << iiter2->second << endl;
    }
  }
  cout << "double values: " << endl;
  map<int, map<const std::string, double>>::const_iterator diter;
  for (diter = default_double.begin(); diter != default_double.end(); ++diter)
  {
    cout << "Detector id: " << diter->first << endl;
    map<const string, double>::const_iterator diter2;
    for (diter2 = diter->second.begin(); diter2 != diter->second.end(); ++diter2)
    {
      cout << diter2->first << ": " << diter2->second << endl;
    }
  }
  cout << "string values: " << endl;
  map<int, map<const std::string, string>>::const_iterator siter;
  for (siter = default_string.begin(); siter != default_string.end(); ++siter)
  {
    cout << "Detector id: " << siter->first << endl;
    map<const string, string>::const_iterator siter2;
    for (siter2 = siter->second.begin(); siter2 != siter->second.end(); ++siter2)
    {
      cout << siter2->first << ": " << siter2->second << endl;
    }
  }
  return;
}

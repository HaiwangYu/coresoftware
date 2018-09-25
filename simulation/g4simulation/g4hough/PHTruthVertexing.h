/*!
 *  \file		PHTruthVertexing.h
 *  \brief		Vertexing using truth info
 *  \author		Haiwang Yu <yuhw@nmsu.edu>
 */

#ifndef __H_PHTruthVertexing_H__
#define __H_PHTruthVertexing_H__

// PHENIX includes
#include <fun4all/SubsysReco.h>

// STL includes
#include <string>
#include <vector>

// forward declarations
class PHCompositeNode;
class SvtxClusterMap;
class SvtxVertexMap;

class PHG4TruthInfoContainer;


/// \class PHTruthVertexing
///
/// \brief Vertexing using truth info
///
class PHTruthVertexing : public SubsysReco {

public:

	PHTruthVertexing(const std::string &name = "PHTruthVertexing");
  virtual ~PHTruthVertexing() {}

	void set_vertex_error(const float & x_err, const float & y_err, const float & z_err) {
		_vertex_error.resize(3);
		_vertex_error[0] = x_err;
		_vertex_error[1] = y_err;
		_vertex_error[2] = z_err;
	}

	const std::vector<float>& get_vertex_error() const {
		return _vertex_error;
	}

protected:

	int Setup(PHCompositeNode *topNode);

	int Process();

private:
	/// create new node output pointers
	int CreateNodes(PHCompositeNode *topNode);

	/// fetch node pointers
	int GetNodes(PHCompositeNode *topNode);

  SvtxClusterMap *_cluster_map;
  SvtxVertexMap *_vertex_map;

	PHG4TruthInfoContainer* _g4truth_container;

	/// manually assigned vertex error (standard dev), cm
	std::vector<float> _vertex_error;

};

#endif //__H_PHTruthVertexing_H__

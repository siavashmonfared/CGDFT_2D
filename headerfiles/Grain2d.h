#ifndef GRAIN2D_H_
#define GRAIN2D_H_

#include "definitions.h"
#include "linterp.h"

class Grain2d {
public:
	// constructors
	Grain2d() {
		_radius = 0; _mass = 0; _momentInertia = 0; _theta = 0; _id = 0; _density = 1;
		_morphologyID = 0; _kn = 0; _ks = 0; _mu = 0; _ncontacts = 0; _gridStep = 1.;
	}

	Grain2d(const double & mass, const Vector2d & position,
			  const double & momentInertia, const double & theta,
			  const Vector2d & cmLset, const Vector2d & cmDims, const vector<double> & lsetvec, const vector<Vector2d> & pointList,
			  const double & radius, const int & id,
			  const int & morphologyID, const double & kn, const double & ks,
			  const double & mu, const double & gridStep):
			  _mass(mass), _position(position),  _momentInertia(momentInertia),
			  _theta(theta), _cmLset(cmLset), _cmDims(cmDims), _radius(radius), _id(id),
			  _morphologyID(morphologyID), _kn(kn), _ks(ks), _mu(mu), _gridStep(gridStep) {
		
		_lsetvec = lsetvec;		  
		_pointList = pointList;
		_ncontacts = 0;
		_density = 1;
		_nodeShears.resize(_pointList.size());
		_nodeContact.resize(_pointList.size());

		// Compute rotation matrix required for updating the pointlist from ref. to actual config. below
		Matrix2d rotMatrix;
		rotMatrix << cos(_theta), -sin(_theta), sin(_theta), cos(_theta);

		for (size_t i = 0; i < _pointList.size(); i++) {
			_nodeShears[i] = 0;
			_nodeContact[i] = 0;
			_pointList[i] = rotMatrix*_pointList[i] + _position;
			// _nodeNormals[i] << 0,0;
		}
			
		double gsLs = 1.;
    	double Lx = _cmDims(0)-1;
    	double Ly = _cmDims(1)-1;
    	int nx = Lx / gsLs + 1;
    	int ny = Ly / gsLs + 1;

    	vector<double> fvals_interp( (Lx*(1./_gridStep)+1) * (Ly*(1./_gridStep)+1) );
    	for (int y = 0 ; y < (Ly*(1./_gridStep)+1) ; y++){  
    	for (int x = 0 ; x < (Lx*(1./_gridStep)+1) ; x++){
    
    	Vector2d point(x*_gridStep,y*_gridStep);
    	double interpVal = bilin_interp( point, _lsetvec, Lx, Ly, _gridStep, nx);
    	fvals_interp[x+(Lx*(1./_gridStep)+1)*y] = interpVal;
    	}
    	}
				  

  		double m, mx, my;
  		m = mx = my = 0.; 
    	for (int y = 0 ; y < (Ly*(1./_gridStep)+1) ; y++){  
    	for (int x = 0 ; x < (Lx*(1./_gridStep)+1) ; x++){
  		if( fvals_interp[x+(Lx*(1./_gridStep)+1)*y] < 0 ){
  		m += 1.; 
  		mx += x;
  		my += y;
  		}
  		}
		}
  		mx = mx / m;
  		my = my / m;

		_cmLset(0) = mx;
		_cmLset(1) = my;
		_lsetvec = fvals_interp;
		_cmDims(0) = Lx*(1./_gridStep)+1;
		_cmDims(1) = Ly*(1./_gridStep)+1;
				  
				  
	}
	
	double getGridValue(vector<double>& fvals, size_t & x, size_t & y, size_t nx)  {
		return fvals[y*nx + x];
	}

	double bilin_interp( Vector2d & point, vector<double> & fvals, double Lx, double Ly, double gridStep, size_t nxLS )  {
    
		double x = point(0);
		double y = point(1);

		// if inside the level set, do bilinear interpolation
		size_t xf 	= floor(x);
		size_t yf 	= floor(y);
		size_t xc 	= ceil(x);
		size_t yc 	= ceil(y);
		double dx 	= x - xf;
		double dy 	= y - yf;
		double b1 	= getGridValue(fvals, xf, yf,nxLS);
		double b2 	= getGridValue(fvals, xc, yf,nxLS) - b1;
		double b3 	= getGridValue(fvals, xf, yc,nxLS) - b1;
		double b4 	= -b2 - getGridValue(fvals, xf, yc,nxLS) + getGridValue(fvals, xc, yc,nxLS);

		double value = b1 + b2*dx + b3*dy + b4*dx*dy;

    	return value;
	}

	// Checks if bounding circles intersect between this and other (1st level contact check)
	bool bcircleGrainIntersection(const Grain2d & other) const {
		if ((other.getPosition()(0)-_position(0))*(other.getPosition()(0)-_position(0)) +
			 (other.getPosition()(1)-_position(1))*(other.getPosition()(1)-_position(1)) <
			 (other.getRadius()+_radius)*(other.getRadius()+ _radius) ) {
			return true;
		}
		return false;
	}
	
	

	
	
	void moveGrain(const Vector2d & amount) {
		_position = _position + amount;
		for (size_t ptid = 0; ptid < _pointList.size(); ptid++) {
			_pointList[ptid] += amount;
		}
	}
	
	// Change methods
	void changeMu(const double & newmu) {
		_mu = newmu;
	}

	void changePos(const Vector2d & pos) {
		Vector2d disp = pos - _position;
		for (size_t ptid = 0; ptid < _pointList.size(); ptid++) {
			_pointList[ptid] += disp;
		}
		_position = pos;	
	}

	void changeRot(const double & rot) {
		double dtheta = rot - _theta;
		_theta = rot;
		double cosd = cos(dtheta);
		double sind = sin(dtheta);
		for (size_t ptid = 0; ptid < _pointList.size(); ptid++) {
			_pointList[ptid] << (_pointList[ptid](0)-_position(0))*cosd - (_pointList[ptid](1)-_position(1))*sind, 
									  (_pointList[ptid](0)-_position(0))*sind + (_pointList[ptid](1)-_position(1))*cosd;
			_pointList[ptid] += _position;
		}
	}

	void changeKn(const double & kn) {
		_kn = kn;
	}

	void changeKs(const double & ks) {
		_ks = ks;
	}

	void changeDensity(const double & density) {
		_mass *= density/_density;
		_momentInertia *= density/_density;
		_density = density;
	}
    
	void changeId(const size_t & id) {
		_id = id;
	}


	// Helper methods
	const double & getMass() const {
		return _mass;
	}
	const Vector2d & getPosition() const {
		return _position;
	}

	const double & getTheta() const {
		return _theta;
	}

	const Vector2d & getCmLset() const {
		return _cmLset;
	}
	
	const Vector2d & getCmDims() const {
		return _cmDims;
	}
	
	const double & getRadius() const {
		return _radius;
	}

	const double & getKn() const {
		return _kn;
	}
	const double & getKs() const {
		return _ks;
	}
	const double & getMu() const {
		return _mu;
	}
	const size_t & getId() const {
		return _id;
	}
    const size_t & getmorphologyID()const {
        return _morphologyID;
    }
    const vector<Vector2d> getPointList() const {
    	return _pointList;
    }
	const vector<double> getLsetVec() const {
    	return _lsetvec;
    }
	const vector<size_t> & getNodeContact() const {
		return _nodeContact;
	}
	const vector<double> & getNodeShears() const {
		return _nodeShears;
	}

	 
	// non const
	vector<double> & getNodeShearsNonConst() {
		return _nodeShears;
	}
	vector<size_t> & getNodeContactNonConst() {
		return _nodeContact;
	}

private:

	double 		_mass;				// particle mass
	Vector2d 	_position; 			// location of center of mass in real space
	double 		_momentInertia;		// moment of inertia in principal frame (purely diagonal terms)
	double 		_theta;				// particle orientation
	Vector2d	_cmLset; 			// center of mass wrt the level set reference configuration (cm: at (0,0), I: diagonal)
	Vector2d 	_cmDims;
	vector<Vector2d>  _pointList; 	// list of points comprising the grain in real space (translated and rotated)
	double 		_radius;			// radius of bounding sphere
	double		_kn;				// normal stiffness
	double		_ks;				// shear stiffness
	double		_mu;				// interparticle friction
	double 		_gridStep;
	double		_density;			// particle density (default 1?)
	size_t      _morphologyID;		// ID representing the morphology type of the particle
	size_t 		_ncontacts;			// number of contacts of the grain (wals + other grains)
	size_t 		_id;				// ID (numbering) of the grain
	vector<double>	_nodeShears;	// shears at each node
	vector<double> _lsetvec;
	vector<size_t>	_nodeContact;  	// index of grain the node is contacting
};

#endif /* Grain2D_H_ */

#ifndef READINPUTFILE_H_
#define READINPUTFILE_H_

#include "definitions.h"
#include "Grain2d.h"


// creates a vector of grain objects from input files
vector<Grain2d> generateGrainsFromFiles(string morphologyMaterialIDs,
										string morphologyDir,
	                                   	string Positions, double gridStep
									    ) {
									

	string line;
	string line_property;
	string line_velocity;

	string line_position;
	string partial;
	string propertyfile;
	istringstream iss;
	ifstream file_information(morphologyMaterialIDs.c_str());		// construct an ifstream and open the grain morphology file
	ifstream file_position(Positions.c_str());						// construct another ifstream and open the grain position file

	
	getline(file_information, line);								// read a line (first) of the morphology file (number of particles)
    iss.str(line);													// copy line string to into iss (we basically bind iss to the line we just read)
    getline(iss, partial, ' ');										// extracts characters until delimiter ' ' is found; the latter is extracted and discarded 
	size_t numberOfGrains = atoi(line.c_str());						// converts the string to an integer
	iss.clear();													// clear the error state of the stream (e.g. end-of-file -> no error)
    char tempfilename[100];
    // Initialize the vector of grain objects
	vector<Grain2d> grainList(numberOfGrains);
	// temp stuff
	Vector2d point;
	Vector2d position;
	double theta;


	// Go through each grain 
	for (size_t grainidx = 0; grainidx < numberOfGrains; grainidx++) {

        // Read morphology index for each particle - each index has each own property .dat file
        getline(file_information, line);
        iss.str(line);
        getline(iss, partial, ' ');
		int morphologyID = atoi(partial.c_str());
        iss.clear();
        sprintf(tempfilename, "grainproperty%d.dat", morphologyID);
        
        propertyfile = morphologyDir + tempfilename;
        ifstream file_gp(propertyfile.c_str());

        // mass
        getline(file_gp, line_property);
        double mass = atof(line_property.c_str());
	
        // moment of inertia
        getline(file_gp, line_property);
		double momentOfInertia = atof(line_property.c_str());

		// cmLset (center of mass)
		getline(file_gp, line_property);
		Vector2d cmLset;
		iss.str(line_property);
		getline(iss, partial, ' ');
		cmLset(0) = atof(partial.c_str());
		getline(iss, partial, ' ');
		cmLset(1) = atof(partial.c_str());
		iss.clear();

		// number of points on the grain surface (INTEGER)
		getline(file_gp, line_property);
		int npoints = atoi(line_property.c_str());

		// the point coordinates
		getline(file_gp, line_property);
		vector<Vector2d> pointList(npoints);
		iss.str(line_property);
		for (int ptidx = 0; ptidx < npoints; ptidx++) {
			getline(iss, partial, ' ');
			point(0) = atof(partial.c_str());
			getline(iss, partial, ' ');
			point(1) = atof(partial.c_str());
			pointList[ptidx] = point;
		}
		iss.clear();

		// bounding box radius
		getline(file_gp, line_property);
		double bboxRadius = atof(line_property.c_str());

		// level set dimensions (INTEGERS)
		getline(file_gp, line_property);
		Vector2d cmDims;
		iss.str(line_property);
		getline(iss, partial, ' ');
		cmDims(0) = atof(partial.c_str());		
		// int xdim = atoi(partial.c_str());
		getline(iss, partial, ' ');
		cmDims(1) = atof(partial.c_str());				
		// int ydim = atoi(partial.c_str());
		iss.clear();

		// level set
		getline(file_gp, line_property);
		vector<double> lsetvec(cmDims(0)*cmDims(1));
		iss.str(line_property);
		for (int i = 0; i < ( cmDims(0)*cmDims(1) ); i++) {
			getline(iss, partial, ' ');
			lsetvec[i] = atof(partial.c_str());
		}
		iss.clear();

		// kn
		getline(file_gp, line_property);
		double kn = atof(line_property.c_str());

        // ks
		getline(file_gp, line_property);
		double ks = atof(line_property.c_str());

        // mu
		getline(file_gp, line_property);
		double mu = atof(line_property.c_str());

        // Clear string for property file ready for next grain - do we need this?
        propertyfile.clear();
        line_property.clear();

	    // Read position and theta from position file
	    getline(file_position, line_position);
		iss.str(line_position);
		getline(iss, partial, ' ');
		position(1) = atof(partial.c_str());
		getline(iss, partial, ' ');
		position(0) = atof(partial.c_str());
		getline(iss, partial, ' ');
		theta = atof(partial.c_str());
		iss.clear();

		// Read velocity and omega from velocity file
		// getline(file_velocity, line_velocity);
		// iss.str(line_velocity);
		// getline(iss, partial, ' ');
		// velocity(0) = 0;//atof(partial.c_str());
		// getline(iss, partial, ' ');
		// velocity(1) = 0; // atof(partial.c_str());
		// getline(iss, partial, ' ');
		// omega = 0;//atof(partial.c_str());
		// iss.clear();

		// Levelset2d lset(lsetvec, xdim, ydim);

		// Update grain object in the vector that was created at the beginning of this function
		grainList[grainidx] = Grain2d(mass, position, momentOfInertia, theta, cmLset, cmDims, lsetvec, pointList, bboxRadius, grainidx, morphologyID, kn, ks, mu, gridStep);
	}

	return grainList;
} // end generateGrainsFromFiles



#endif // READINPUTFILE_H_ 

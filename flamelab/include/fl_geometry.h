#ifndef FL_GEOMETRY_H
#define FL_GEOMETRY_H
#include "fl_params.h"
namespace FlameLab{

	class Geometry{

		int nCell; // number of max finite volume cells
		int axialPosition;
		real rLength; // length of the reactor
		real aspectRatio;
		
	public:
		Geometry():nCell(10){// default number of max u,ber of cells			
			dz.resize(nCell);
		}

		Geometry(int n);		
		~Geometry();
		void setnCells(int n);
		void setLength(real len);
		void setAspectRatio(real ar);
		int getnCells() const;
		real getLength() const;
		void descretize();
		vector<real> getGeometry() const;
		real getAspectRatio() const;
		void setAxialPosition(int n);
		int getAxialPosition()const;

	protected:
		vector<real> dz; // vector holding descretization info


	};
};

#endif
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <stdlib.h>
#include <set>
#include <iterator>
#include <numeric> 
#include <chrono>
#include <random>
#include <complex>
#include <algorithm>
#include <mpi.h>

#include "../headerfiles/definitions.h"
#include "../headerfiles/Grain2d.h"
#include "../headerfiles/readInputFile.h"
#include "../headerfiles/linterp.h"

using namespace std;

#define QLG 4
#define HLG 2

struct node { 

	int I[QLG]; 
	size_t ilg[QLG]; 

	int nblg;
	int nbolg;
	
	node():nblg(0),nbolg(0) { }
};



/*for CG*/
int   ncom;
double *pcom, *xicom, (*nrfunc)(double []);
/*for computation*/
vector<node> nodes;  
/*for MPI*/
MPI_Comm comm2d;
int NeighBor[4];
int S=0, E=1, N=2, W=3;
int ndims = 2;
int ndof = 2;

MPI_Datatype RtScalFieldBotDOUB, RtScalFieldTopDOUB, RtScalFieldBot, RtScalFieldTop, RtVecFieldBot, RtVecFieldTop, CtVecField, CtScalField, CtScalFieldDOUB;

MPI_Status status;
int flag = 1; 
int* lengthX;
int* lengthY;
double* p;
double* force;
int* Id;
int nx_loc, ny_loc, xsum, ysum, nx_fill, ny_fill;
/*DFT*/
double* Eta;
double* rho;
double wff,wmf,mu,Temp, Boltz,yvar;
double Lx_lt_shift, Lx_rt_shift, Ly_lt_shift, Ly_rt_shift;
double* sigF;
/*readInput*/
size_t type, nit;
double Tstar, Tcw;
int nx, ny;	   
int nxp, nyp;
double gridStep, K1,K0; 	   		   
double strain, ftol, Rpore; 	
double XC, YC, obj_tol;
double rho_thresh_ub, rho_thresh_lb;
size_t nu_part;
int Lnx, Lny;

struct tensor
{
	double xx, xy;
	double yx, yy;

	void reset() {
		xx = xy = yx = yy = 0.0;
	}
};






int ipx(int x, int y, int xsize)
{
	return 2 * (x + xsize * y);
}

void decompose1d (int length[], int ndr ,int npr){

int rem = ndr % npr;

if(rem == 0){
int intL = ndr / npr;
for(int k = 0 ; k < npr; k++) length[k] = intL;
}
	
if (rem == 1){
int intL = (ndr - rem) / npr;
for(int k = 0 ; k < (npr-1); k++) length[k] = intL;
length[npr-1] = rem + intL;
}
	
if (rem > 1){
for(int k = 0 ; k < (npr-1); k++) length[k] = (ndr - rem) / npr;
length[npr-1] = (ndr - rem) / npr + rem ;
}
}


void init(double GridStep, int ndof, size_t nx_loc, size_t ny_loc, int xsum, int ysum){

	bool xulF, yulF, xllF, yllF;
	bool xulT, yulT, xllT, yllT;
	bool xug, xlg, yug, ylg;
	
	nodes.resize(nx_loc * ny_loc);

	for (size_t y = 0 ; y < ny_loc ; y++) {
			for (size_t x = 0 ; x < nx_loc ; x++) {
				
				int x_global = xsum + x;
				int y_global = ysum + y;

				xug = (x_global == nx - 1) ? false : true;
				xlg = (x_global == 0) ? false : true;
				yug = (y_global == ny - 1) ? false : true;
				ylg = (y_global == 0) ? false : true;

				xulF = (x == nx_loc - 1) ? false : true;//if false - continue
				yulF = (y == ny_loc - 1) ? false : true;//if false - continue 
				xllF = (x == 0) ? false : true; //if false - continue 
				yllF = (y == 0) ? false : true; //if false - continue

				xulT = (x == nx_loc - 1) ? true : false;//if true - continue
				yulT = (y == ny_loc - 1) ? true : false;//if true - continue 
				xllT = (x == 0) ? true : false; //if true - continue 
				yllT = (y == 0) ? true : false; //if true - continue
				
				int ivlg = 0;
				if (xug && xulF){ 
				nodes[ipx(x,y,nx_loc)/2].ilg[ivlg++] 	= 	ipx(x+1,y,nx_loc);
				}
				if (xulT){ 
				nodes[ipx(x,y,nx_loc)/2].ilg[ivlg++] 	= (nx_loc * ny_loc + y) * ndof;			
				}
				
				
				if (yug && yulF){
				nodes[ipx(x,y,nx_loc)/2].ilg[ivlg++] 	= ipx(x,y+1,nx_loc);
				}
				if (yulT){
				nodes[ipx(x,y,nx_loc)/2].ilg[ivlg++] 	= (nx_loc * ny_loc + ny_loc * 2 + (nx_loc + 2) + (x+1)) * ndof;		
				}
													
				int ivolg = 0;
				
				if (xlg && xllF) nodes[ipx(x,y,nx_loc)/2].ilg[HLG+ivolg++] 	= ipx(x-1,y,nx_loc);
				if (xllT){
				nodes[ipx(x,y,nx_loc)/2].ilg[HLG+ivolg++] = (nx_loc * ny_loc + ny_loc + (y)) * ndof;				
				}

				if (ylg && yllF) nodes[ipx(x,y,nx_loc)/2].ilg[HLG+ivolg++] 	= ipx(x,y-1,nx_loc);
				if (yllT){ 
				nodes[ipx(x,y,nx_loc)/2].ilg[HLG+ivolg++] = (nx_loc * ny_loc + ny_loc * 2 + (x+1)) * ndof;
				}
				
										
				nodes[ipx(x,y,nx_loc)/2].nblg = ivlg;
				nodes[ipx(x,y,nx_loc)/2].nbolg = ivolg;
		}
	}
}



void readSimParam(){
	
    const char * name = "../input/input_12152020.dat";
    ifstream file(name);
    file >> Lnx >> Lny;
	file >> nxp >> nyp;
	file >> Tstar >> gridStep >> Tcw >> yvar;
	file >> obj_tol >> nu_part;
	file >> rho_thresh_ub >> rho_thresh_lb;
	file >> Lx_lt_shift >> Lx_rt_shift;
	file >> Ly_lt_shift >> Ly_rt_shift;
}








void  Import_Structure_arb_particles_fromPts( vector<Grain2d> grains, int idincl, double ideta, size_t nx_loc, size_t ny_loc, int xsum, int ysum, int xshift, int yshift, int rrkk ){
	
	for (size_t k = 0 ; k < grains.size() ; k++){
	vector<Vector2d> _pointList = grains[k].getPointList();
		for (size_t j = 0 ; j < _pointList.size() ; j++){
			
			int xc_global = floor( (_pointList[j](0) + xshift) / gridStep + 1.);
			int yc_global = floor( (_pointList[j](1) + yshift) / gridStep + 1.);
			
			size_t x = xc_global - xsum + 1;
			size_t y = yc_global - ysum + 1;
			
			if ( rrkk == 0 && xc_global < lengthX[rrkk] && yc_global < lengthY[rrkk] ){
   			Id[ipx(x,y,nx_loc)/2] = idincl;
   			Eta[ipx(x,y,nx_loc)/2] = ideta;
			}
			if ( rrkk > 0 && xc_global < lengthX[rrkk-1] && xc_global > lengthX[rrkk] && yc_global < lengthY[rrkk-1] && yc_global > lengthY[rrkk]){
   			Id[ipx(x,y,nx_loc)/2] = idincl;
   			Eta[ipx(x,y,nx_loc)/2] = ideta;
			}
		}
	}
	
	MPI_Barrier(comm2d);
	/** scalar fields **/
	/** to left from right **/
	MPI_Sendrecv (&(Id[0]),1,CtScalField,NeighBor[W],flag,&(Id[nx_loc*ny_loc*1]),1 * ny_loc,MPI_INT,NeighBor[E],flag,comm2d,&status);
	/** to right from left **/
	MPI_Sendrecv (&(Id[1 * nx_loc - 1]),1,CtScalField,NeighBor[E],flag,&(Id[1 * nx_loc * ny_loc + 1 * ny_loc]), 1 * ny_loc, MPI_INT,NeighBor[W],flag,comm2d,&status);	
	/** to top from bottom **/
	MPI_Sendrecv (&(Id[0]),1,RtScalFieldTop,NeighBor[N],flag,&(Id[1*nx_loc*ny_loc + 2 * 1 * ny_loc + 1 * nx_loc + 2 * 1]),(1 * nx_loc + 2 * 1),MPI_INT,NeighBor[S],flag,comm2d,&status);
	/** to bottom from top **/
	MPI_Sendrecv (&(Id[0]),1,RtScalFieldBot,NeighBor[S],flag,&(Id[1 * nx_loc*ny_loc + 1 * 2 * ny_loc]), 1 * nx_fill,MPI_INT,NeighBor[N],flag,comm2d,&status);
	
	/** scalar fields **/
	/** to left from right **/
	MPI_Sendrecv (&(Eta[0]),1,CtScalFieldDOUB,NeighBor[W],flag,&(Eta[nx_loc*ny_loc*1]),1 * ny_loc,MPI_DOUBLE,NeighBor[E],flag,comm2d,&status);
	/** to right from left **/
	MPI_Sendrecv (&(Eta[1 * nx_loc - 1]),1,CtScalFieldDOUB,NeighBor[E],flag,&(Eta[1 * nx_loc * ny_loc + 1 * ny_loc]), 1 * ny_loc, MPI_DOUBLE,NeighBor[W],flag,comm2d,&status);	
	/** to top from bottom **/
	MPI_Sendrecv (&(Eta[0]),1,RtScalFieldTopDOUB,NeighBor[N],flag,&(Eta[1*nx_loc*ny_loc + 2 * 1 * ny_loc + 1 * nx_loc + 2 * 1]),(1 * nx_loc + 2 * 1),MPI_DOUBLE,NeighBor[S],flag,comm2d,&status);
	/** to bottom from top **/
	MPI_Sendrecv (&(Eta[0]),1,RtScalFieldBotDOUB,NeighBor[S],flag,&(Eta[1 * nx_loc*ny_loc + 1 * 2 * ny_loc]), 1 * nx_fill,MPI_DOUBLE,NeighBor[N],flag,comm2d,&status);
		
}

	
void  Import_Structure_arb_particles_fromLS( vector<Grain2d> grains, int idincl, double ideta, size_t nx_loc, size_t ny_loc, int xsum, int ysum, int xshift, int yshift, int xhi, int yhi, int rrkk ){

	for (size_t k = 0 ; k < grains.size() ; k++){
		Vector2d cmLS = grains[k].getCmLset();
		Vector2d cmDims = grains[k].getCmDims();
		Vector2d ppos = grains[k].getPosition();
		vector<double> lsvec = grains[k].getLsetVec();
		Matrix2d rotMatrix;
		double theta = grains[k].getTheta();
		rotMatrix << cos(theta), -sin(theta), sin(theta), cos(theta);		


		size_t count = 0;
		for(size_t j = 0 ; j < cmDims(1) ; j++){
			for(size_t i = 0 ; i < cmDims(0) ; i++){
				double lsval = lsvec[count];
				count++;
				if(lsval <= 0){
					
					Vector2d gridPoint(i-cmLS(0),j-cmLS(1));
					gridPoint = rotMatrix*gridPoint;
					double iRot = gridPoint(0);
					double jRot = gridPoint(1);
					
					size_t xc_global = floor ( ( iRot+cmLS(0) ) + xshift + (ppos(1)/gridStep+1)  );
					size_t yc_global = floor ( ( jRot+cmLS(1) ) + yshift + (ppos(0)/gridStep+1)  );


					size_t x = xc_global - xsum ;
					size_t y = yc_global - ysum ;
					
					if (  xc_global >= xsum && xc_global < (xsum+xhi) && yc_global >= ysum && yc_global < (ysum+yhi) ){
						
   						Id[ipx(x,y,nx_loc)/2] = idincl;
   						Eta[ipx(x,y,nx_loc)/2] = ideta;
					}

				}
			}
		}
	}
	
	
	MPI_Barrier(comm2d);
	/** scalar fields **/
	/** to left from right **/
	MPI_Sendrecv (&(Id[0]),1,CtScalField,NeighBor[W],flag,&(Id[nx_loc*ny_loc*1]),1 * ny_loc,MPI_INT,NeighBor[E],flag,comm2d,&status);
	/** to right from left **/
	MPI_Sendrecv (&(Id[1 * nx_loc - 1]),1,CtScalField,NeighBor[E],flag,&(Id[1 * nx_loc * ny_loc + 1 * ny_loc]), 1 * ny_loc, MPI_INT,NeighBor[W],flag,comm2d,&status);	
	/** to top from bottom **/
	MPI_Sendrecv (&(Id[0]),1,RtScalFieldTop,NeighBor[N],flag,&(Id[1*nx_loc*ny_loc + 2 * 1 * ny_loc + 1 * nx_loc + 2 * 1]),(1 * nx_loc + 2 * 1),MPI_INT,NeighBor[S],flag,comm2d,&status);
	/** to bottom from top **/
	MPI_Sendrecv (&(Id[0]),1,RtScalFieldBot,NeighBor[S],flag,&(Id[1 * nx_loc*ny_loc + 1 * 2 * ny_loc]), 1 * nx_fill,MPI_INT,NeighBor[N],flag,comm2d,&status);
	
	/** scalar fields **/
	/** to left from right **/
	MPI_Sendrecv (&(Eta[0]),1,CtScalFieldDOUB,NeighBor[W],flag,&(Eta[nx_loc*ny_loc*1]),1 * ny_loc,MPI_DOUBLE,NeighBor[E],flag,comm2d,&status);
	/** to right from left **/
	MPI_Sendrecv (&(Eta[1 * nx_loc - 1]),1,CtScalFieldDOUB,NeighBor[E],flag,&(Eta[1 * nx_loc * ny_loc + 1 * ny_loc]), 1 * ny_loc, MPI_DOUBLE,NeighBor[W],flag,comm2d,&status);	
	/** to top from bottom **/
	MPI_Sendrecv (&(Eta[0]),1,RtScalFieldTopDOUB,NeighBor[N],flag,&(Eta[1*nx_loc*ny_loc + 2 * 1 * ny_loc + 1 * nx_loc + 2 * 1]),(1 * nx_loc + 2 * 1),MPI_DOUBLE,NeighBor[S],flag,comm2d,&status);
	/** to bottom from top **/
	MPI_Sendrecv (&(Eta[0]),1,RtScalFieldBotDOUB,NeighBor[S],flag,&(Eta[1 * nx_loc*ny_loc + 1 * 2 * ny_loc]), 1 * nx_fill,MPI_DOUBLE,NeighBor[N],flag,comm2d,&status);
	
}


void  Import_Structure_spherical_particles( const string& _name, int idincl, double ideta, size_t npart, size_t nx_loc, size_t ny_loc, int xsum, int ysum, int xshift, int yshift ){
	
	double xc,yc,crit,rad;

	fstream file(_name);

	for (size_t k = 0 ; k < npart ; k++){
	file >> xc >> yc >> rad;
		
	xc = xc + xshift;
	yc = yc + yshift;
		
	for (size_t y = 0 ; y < ny_loc ; y++) {
			for (size_t x = 0 ; x < nx_loc ; x++) {

				int x_global = xsum + x;
				int y_global = ysum + y;

				x_global = x_global - 1;
				y_global = y_global - 1;

				crit=(1.*x_global*gridStep-1.*xc)*(1.*x_global*gridStep-1.*xc)+(1.*y_global*gridStep-1.*yc)*(1.*y_global*gridStep-1.*yc);

				if(crit<=1.*rad*rad){
   				Id[ipx(x,y,nx_loc)/2] 	= idincl;
				Eta[ipx(x,y,nx_loc)/2] 	= ideta;
    				}
    			}
    		}
	}
			
	MPI_Barrier(comm2d);
	/** scalar fields **/
	/** to left from right **/
	MPI_Sendrecv (&(Id[0]),1,CtScalField,NeighBor[W],flag,&(Id[nx_loc*ny_loc*1]),1 * ny_loc,MPI_INT,NeighBor[E],flag,comm2d,&status);
	/** to right from left **/
	MPI_Sendrecv (&(Id[1 * nx_loc - 1]),1,CtScalField,NeighBor[E],flag,&(Id[1 * nx_loc * ny_loc + 1 * ny_loc]), 1 * ny_loc, MPI_INT,NeighBor[W],flag,comm2d,&status);	
	/** to top from bottom **/
	MPI_Sendrecv (&(Id[0]),1,RtScalFieldTop,NeighBor[N],flag,&(Id[1*nx_loc*ny_loc + 2 * 1 * ny_loc + 1 * nx_loc + 2 * 1]),(1 * nx_loc + 2 * 1),MPI_INT,NeighBor[S],flag,comm2d,&status);
	/** to bottom from top **/
	MPI_Sendrecv (&(Id[0]),1,RtScalFieldBot,NeighBor[S],flag,&(Id[1 * nx_loc*ny_loc + 1 * 2 * ny_loc]), 1 * nx_fill,MPI_INT,NeighBor[N],flag,comm2d,&status);
	
	/** scalar fields **/
	/** to left from right **/
	MPI_Sendrecv (&(Eta[0]),1,CtScalFieldDOUB,NeighBor[W],flag,&(Eta[nx_loc*ny_loc*1]),1 * ny_loc,MPI_DOUBLE,NeighBor[E],flag,comm2d,&status);
	/** to right from left **/
	MPI_Sendrecv (&(Eta[1 * nx_loc - 1]),1,CtScalFieldDOUB,NeighBor[E],flag,&(Eta[1 * nx_loc * ny_loc + 1 * ny_loc]), 1 * ny_loc, MPI_DOUBLE,NeighBor[W],flag,comm2d,&status);	
	/** to top from bottom **/
	MPI_Sendrecv (&(Eta[0]),1,RtScalFieldTopDOUB,NeighBor[N],flag,&(Eta[1*nx_loc*ny_loc + 2 * 1 * ny_loc + 1 * nx_loc + 2 * 1]),(1 * nx_loc + 2 * 1),MPI_DOUBLE,NeighBor[S],flag,comm2d,&status);
	/** to bottom from top **/
	MPI_Sendrecv (&(Eta[0]),1,RtScalFieldBotDOUB,NeighBor[S],flag,&(Eta[1 * nx_loc*ny_loc + 1 * 2 * ny_loc]), 1 * nx_fill,MPI_DOUBLE,NeighBor[N],flag,comm2d,&status);
		
}




double vol_frac( double tmp[], int idp, int rrkk, int xls, int xrs, int yls, int yrs ){

	double count = 0.0;
	double count_A = 0.0;
	double count_B = 0.0;
	bool res_cond, box_cond;

	for (size_t y = 0 ; y < ny_loc ; y++) {
	for (size_t x = 0 ; x < nx_loc ; x++) {

	int x_global = (xsum + x);
	int y_global = (ysum + y);

    res_cond = ( x_global >= xls ) && ( x_global <  (nx) - xrs ) && ( y_global >= yls ) && ( y_global < (ny) - yrs ) ? true: false;//if true -> continue 
    box_cond = ( x_global >= xls ) && ( x_global <  (nx) - xrs ) && ( y_global >= yls ) && ( y_global < (ny) - yrs ) ? false: true;//if false -> continue 

	if(Id[ipx(x,y,nx_loc)/2] == idp){
	count++;
	}

	if(Id[ipx(x,y,nx_loc)/2] == idp && res_cond){
	count_A++;
	}

	if(Id[ipx(x,y,nx_loc)/2] == idp && box_cond){
	count_B++;
	}

	}
	}

	double fs,fp, count_tot, count_tot_A, count_tot_B;
	MPI_Barrier(comm2d);
	MPI_Allreduce(&count, &count_tot, 1, MPI_DOUBLE, MPI_SUM, comm2d);
	MPI_Allreduce(&count_A, &count_tot_A, 1, MPI_DOUBLE, MPI_SUM, comm2d);
	MPI_Allreduce(&count_B, &count_tot_B, 1, MPI_DOUBLE, MPI_SUM, comm2d);
	double tot = nx * ny;
	double den = pow((cbrt(nx * ny)-1),3);
	fp = (count_tot)/den;
	fs = 1.0 - fp;

	tmp[0] = count_tot;
	tmp[1] = count_tot_A;
	tmp[2] = count_tot_B;
	return fp;
}



void save_MPI_info( int rrkk,int nx_loc,int ny_loc, int xcoor,int ycoor, int xs, int xsp, int ys, int ysp ){
    FILE * sortie;
    char nomfic[256];
    sprintf(nomfic, "../output/mpi_info.%04i", rrkk);
    sortie = fopen(nomfic, "w+");
    fprintf(sortie,"%i %i %i %i %i %i %i %i %i\n",rrkk,nx_loc,ny_loc,xcoor,ycoor,xs,xs+xsp,ys,ys+ysp);
    fclose(sortie);	
}

/*
void save_MPI_info(int rrkk,int nx_loc,int ny_loc, int nz_loc, int xcoor,int ycoor,int zcoor,int xs,int xsp,int ys,int ysp,int zs,int zsp){
    FILE * sortie;
    char nomfic[256];
    sprintf(nomfic, "mpi_info.%04i", rrkk);
    sortie = fopen(nomfic, "w+");
    fprintf(sortie,"%i %i %i %i %i %i %i %i %i %i %i %i %i\n",rrkk,nx_loc,ny_loc,nz_loc,xcoor,ycoor,zcoor,xs,xs+xsp,ys,ys+ysp,zs,zs+zsp);
    fclose(sortie);	
}
*/


void initialize_rho( double ri, int idp, int nx_loc, int ny_loc ){

	for (size_t y = 0 ; y < ny_loc ; y++) {
	for (size_t x = 0 ; x < nx_loc ; x++) {

	if(Id[ipx(x,y,nx_loc)/2] == idp) rho[ipx(x,y,nx_loc)/2] = ri;
	if(Id[ipx(x,y,nx_loc)/2] != idp) rho[ipx(x,y,nx_loc)/2] = 0.0;

	}	
	}

	/** scalar fields **/
	/** to left from right **/
	MPI_Sendrecv (&(rho[0]),1,CtScalFieldDOUB,NeighBor[W],flag,&(rho[nx_loc*ny_loc*1]),1 * ny_loc,MPI_DOUBLE,NeighBor[E],flag,comm2d,&status);
	/** to right from left **/
	MPI_Sendrecv (&(rho[1 * nx_loc - 1]),1,CtScalFieldDOUB,NeighBor[E],flag,&(rho[1 * nx_loc * ny_loc + 1 * ny_loc]), 1 * ny_loc, MPI_DOUBLE,NeighBor[W],flag,comm2d,&status);	
	/** to top from bottom **/
	MPI_Sendrecv (&(rho[0]),1,RtScalFieldTopDOUB,NeighBor[N],flag,&(rho[1*nx_loc*ny_loc + 2 * 1 * ny_loc + 1 * nx_loc + 2 * 1]),(1 * nx_loc + 2 * 1),MPI_DOUBLE,NeighBor[S],flag,comm2d,&status);
	/** to bottom from top **/
	MPI_Sendrecv (&(rho[0]),1,RtScalFieldBotDOUB,NeighBor[S],flag,&(rho[1 * nx_loc*ny_loc + 1 * 2 * ny_loc]), 1 * nx_fill,MPI_DOUBLE,NeighBor[N],flag,comm2d,&status);

}



void DFT_ObjFunc( double tol, double Temp, double mu, int nx_loc, int ny_loc,int nx_fill, int ny_fill, int rrkk, int xls, int xrs, int yls, int yrs ){
	
	double GP, err_tot, rho_tot;
	double err = 1000.0;
	double errsum_loc = 0.0;
	double den = nx * ny;
  	 bool res_cond;

    double* rhom;
	rhom = new double[nx_fill * ny_fill];

	while (err > tol){

	for (size_t i = 0 ; i < (nx_fill * ny_fill) ; i++) rhom[i] = rho[i];

	double rho_loc = 0.0;
	for (size_t i = 0 ; i < (nx_loc * ny_loc) ; i++) rho_loc += rhom[i];

	for (size_t y = 0 ; y < ny_loc ; y++){
	for (size_t x = 0 ; x < nx_loc ; x++){

	int x_global = (xsum + x);
	int y_global = (ysum + y);

    res_cond = ( x_global >= xls ) && ( x_global <  (nx) - xrs ) && ( y_global >= yls ) && ( y_global < (ny) - yrs ) ? false: true; 

	GP = 0.0;
	
	for (int l = 0 ; l < nodes[ipx(x,y,nx_loc)/2].nblg ; ++l) GP += rhom[nodes[ipx(x,y,nx_loc)/2].ilg[l]/2] + yvar * (1.0 - Eta[nodes[ipx(x,y,nx_loc)/2].ilg[l]/2]);
	for (int l = HLG ; l < HLG + nodes[ipx(x,y,nx_loc)/2].nbolg ; ++l) GP += rhom[nodes[ipx(x,y,nx_loc)/2].ilg[l]/2] + yvar * (1.0 - Eta[nodes[ipx(x,y,nx_loc)/2].ilg[l]/2]);
	rho[ipx(x,y,nx_loc)/2] = Eta[ipx(x,y,nx_loc)/2] / (1.0 + exp( (-1.0/Boltz/Temp) * (mu + wff * GP) ) );
	//if (res_cond) rho[ipx(x,y,z,nx_loc,ny_loc)/6] = 0.999;
	}
	}

	/** scalar fields **/
	/** to left from right **/
	MPI_Sendrecv (&(rho[0]),1,CtScalFieldDOUB,NeighBor[W],flag,&(rho[nx_loc*ny_loc*1]),1 * ny_loc,MPI_DOUBLE,NeighBor[E],flag,comm2d,&status);
	/** to right from left **/
	MPI_Sendrecv (&(rho[1 * nx_loc - 1]),1,CtScalFieldDOUB,NeighBor[E],flag,&(rho[1 * nx_loc * ny_loc + 1 * ny_loc]), 1 * ny_loc, MPI_DOUBLE,NeighBor[W],flag,comm2d,&status);	
	/** to top from bottom **/
	MPI_Sendrecv (&(rho[0]),1,RtScalFieldTopDOUB,NeighBor[N],flag,&(rho[1*nx_loc*ny_loc + 2 * 1 * ny_loc + 1 * nx_loc + 2 * 1]),(1 * nx_loc + 2 * 1),MPI_DOUBLE,NeighBor[S],flag,comm2d,&status);
	/** to bottom from top **/
	MPI_Sendrecv (&(rho[0]),1,RtScalFieldBotDOUB,NeighBor[S],flag,&(rho[1 * nx_loc*ny_loc + 1 * 2 * ny_loc]), 1 * nx_fill,MPI_DOUBLE,NeighBor[N],flag,comm2d,&status);

	for (size_t i = 0 ; i < (nx_loc * ny_loc) ; i++) errsum_loc += pow( (rhom[i] - rho[i]) , 2);
	
	MPI_Barrier(comm2d);
	MPI_Allreduce(&errsum_loc, &err_tot, 1, MPI_DOUBLE, MPI_SUM, comm2d);
	MPI_Allreduce(&rho_loc, &rho_tot, 1, MPI_DOUBLE, MPI_SUM, comm2d);
	errsum_loc = 0.0;
	err = err_tot/den;
	}
}



void density_gradient(int IDP, int IDS, double trho[], double drx[], double dry[], double gs){

	int p0,p1,p2,p3,p4,p5,p6;

	double Dx1, Dy1, Dx0, Dy0, Dx, Dy;
	bool xulF, yulF, xllF, yllF;
	bool xulT, yulT, xllT, yllT;
	bool xug, xlg, yug, ylg;

	for (size_t y = 0 ; y < ny_loc ; y++) {
	for (size_t x = 0 ; x < nx_loc ; x++) {
		
	int x_global = (xsum + x);
	int y_global = (ysum + y);

	if (Id[ipx(x,y,nx_loc)/2] == IDP){

	xug = (x_global == nx - 1) ? false : true;
	xlg = (x_global == 0) ? false : true;
	yug = (y_global == ny - 1) ? false : true;
	ylg = (y_global == 0) ? false : true;

	xulF = (x == nx_loc - 1) ? false : true;//if false - continue
	yulF = (y == ny_loc - 1) ? false : true;//if false - continue 
	xllF = (x == 0) ? false : true; //if false - continue 
	yllF = (y == 0) ? false : true; //if false - continue

	xulT = (x == nx_loc - 1) ? true : false;//if true - continue
	yulT = (y == ny_loc - 1) ? true : false;//if true - continue 
	xllT = (x == 0) ? true : false; //if true - continue 
	yllT = (y == 0) ? true : false; //if true - continue

	if (xulF)	p1 = Id[ipx(x+1,y,nx_loc)/2];
	if (xulT)	p1 = Id[ (nx_loc * ny_loc + y) ];

	if (yulF)	p2 = Id[ipx(x,y+1,nx_loc)/2];
	if (yulT)	p2 = Id[ (nx_loc * ny_loc + ny_loc * 2 + (nx_loc + 2) + (x+1)) ];
		
	if (xllF)	p4 = Id[ipx(x-1,y,nx_loc)/2];
	if (xllT)	p4 = Id[(nx_loc * ny_loc + ny_loc + (y))];
				
	if (yllF)	p5 = Id[ipx(x,y-1,nx_loc)/2];
	if (yllT)	p5 = Id[(nx_loc * ny_loc + ny_loc * 2 + (x+1))];
				
	if (xulF)	Dx1 = trho[ipx(x+1,y,nx_loc)/2];
	if (xulT)	Dx1 = trho[ (nx_loc * ny_loc + y) ];

	if (yulF)	Dy1 = trho[ipx(x,y+1,nx_loc)/2];
	if (yulT)	Dy1 = trho[ (nx_loc * ny_loc + ny_loc * 2 + (nx_loc + 2) + (x+1)) ];
		
	if (xllF)	Dx0 = trho[ipx(x-1,y,nx_loc)/2];
	if (xllT)	Dx0 = trho[(nx_loc * ny_loc + ny_loc + (y))];
				
	if (yllF)	Dy0 = trho[ipx(x,y-1,nx_loc)/2];
	if (yllT)	Dy0 = trho[(nx_loc * ny_loc + ny_loc * 2 + (x+1))];
				
	if (p1 == IDP && p4 == IDP) Dx = (Dx1 - Dx0) / (2.0*gs);
	if (p1 == IDP && p4 == IDS) Dx = (Dx1 - trho[ipx(x,y,nx_loc)/2]) / (gs);
	if (p1 == IDS && p4 == IDP) Dx = (trho[ipx(x,y,nx_loc)/2] - Dx0) / (gs);

	if (p2 == IDP && p5 == IDP) Dy = (Dy1 - Dy0) / (2.0*gs);
	if (p2 == IDP && p5 == IDS) Dy = (Dy1 - trho[ipx(x,y,nx_loc)/2]) / (gs);
	if (p2 == IDS && p5 == IDP) Dy = (trho[ipx(x,y,nx_loc)/2] - Dy0) / (gs);

	drx[ipx(x,y,nx_loc)/2] = Dx;
	dry[ipx(x,y,nx_loc)/2] = Dy;

	}

	if(Id[ipx(x,y,nx_loc)/2] == IDS){
	drx[ipx(x,y,nx_loc)/2] = 0.0;
	dry[ipx(x,y,nx_loc)/2] = 0.0;
	}
		
	}
	}
	
}

void smooth_particles(int IDP, int IDS){

	int p1,p2,p4,p5;
	p1 = p2 = p4 = p5 = 0;
	bool xulF, yulF, xllF, yllF;
	bool xulT, yulT, xllT, yllT;
	bool xug, xlg, yug, ylg;

	for (size_t y = 0 ; y < ny_loc ; y++) {
	for (size_t x = 0 ; x < nx_loc ; x++) {
		
	int x_global = (xsum + x);
	int y_global = (ysum + y);

	if (Id[ipx(x,y,nx_loc)/2] == IDP){

	xug = (x_global == nx - 1) ? false : true;
	xlg = (x_global == 0) ? false : true;
	yug = (y_global == ny - 1) ? false : true;
	ylg = (y_global == 0) ? false : true;

	xulF = (x == nx_loc - 1) ? false : true;//if false - continue
	yulF = (y == ny_loc - 1) ? false : true;//if false - continue 
	xllF = (x == 0) ? false : true; //if false - continue 
	yllF = (y == 0) ? false : true; //if false - continue

	xulT = (x == nx_loc - 1) ? true : false;//if true - continue
	yulT = (y == ny_loc - 1) ? true : false;//if true - continue 
	xllT = (x == 0) ? true : false; //if true - continue 
	yllT = (y == 0) ? true : false; //if true - continue

	if (xulF)	p1 = Id[ipx(x+1,y,nx_loc)/2];
	if (xulT)	p1 = Id[ (nx_loc * ny_loc + y) ];

	if (yulF)	p2 = Id[ipx(x,y+1,nx_loc)/2];
	if (yulT)	p2 = Id[ (nx_loc * ny_loc + ny_loc * 2 + (nx_loc + 2) + (x+1)) ];
		
	if (xllF)	p4 = Id[ipx(x-1,y,nx_loc)/2];
	if (xllT)	p4 = Id[(nx_loc * ny_loc + ny_loc + (y))];
				
	if (yllF)	p5 = Id[ipx(x,y-1,nx_loc)/2];
	if (yllT)	p5 = Id[(nx_loc * ny_loc + ny_loc * 2 + (x+1))];
				
	if (p1 + p2 + p4 + p5 >= 3){
	Id[ipx(x,y,nx_loc)/2] = IDS;
	Eta[ipx(x,y,nx_loc)/2] = 0.0;
	}
	}
	}
	}
	
	MPI_Barrier(comm2d);
	/** scalar fields **/
	/** to left from right **/
	MPI_Sendrecv (&(Id[0]),1,CtScalField,NeighBor[W],flag,&(Id[nx_loc*ny_loc*1]),1 * ny_loc,MPI_INT,NeighBor[E],flag,comm2d,&status);
	/** to right from left **/
	MPI_Sendrecv (&(Id[1 * nx_loc - 1]),1,CtScalField,NeighBor[E],flag,&(Id[1 * nx_loc * ny_loc + 1 * ny_loc]), 1 * ny_loc, MPI_INT,NeighBor[W],flag,comm2d,&status);	
	/** to top from bottom **/
	MPI_Sendrecv (&(Id[0]),1,RtScalFieldTop,NeighBor[N],flag,&(Id[1*nx_loc*ny_loc + 2 * 1 * ny_loc + 1 * nx_loc + 2 * 1]),(1 * nx_loc + 2 * 1),MPI_INT,NeighBor[S],flag,comm2d,&status);
	/** to bottom from top **/
	MPI_Sendrecv (&(Id[0]),1,RtScalFieldBot,NeighBor[S],flag,&(Id[1 * nx_loc*ny_loc + 1 * 2 * ny_loc]), 1 * nx_fill,MPI_INT,NeighBor[N],flag,comm2d,&status);
	
	/** scalar fields **/
	/** to left from right **/
	MPI_Sendrecv (&(Eta[0]),1,CtScalFieldDOUB,NeighBor[W],flag,&(Eta[nx_loc*ny_loc*1]),1 * ny_loc,MPI_DOUBLE,NeighBor[E],flag,comm2d,&status);
	/** to right from left **/
	MPI_Sendrecv (&(Eta[1 * nx_loc - 1]),1,CtScalFieldDOUB,NeighBor[E],flag,&(Eta[1 * nx_loc * ny_loc + 1 * ny_loc]), 1 * ny_loc, MPI_DOUBLE,NeighBor[W],flag,comm2d,&status);	
	/** to top from bottom **/
	MPI_Sendrecv (&(Eta[0]),1,RtScalFieldTopDOUB,NeighBor[N],flag,&(Eta[1*nx_loc*ny_loc + 2 * 1 * ny_loc + 1 * nx_loc + 2 * 1]),(1 * nx_loc + 2 * 1),MPI_DOUBLE,NeighBor[S],flag,comm2d,&status);
	/** to bottom from top **/
	MPI_Sendrecv (&(Eta[0]),1,RtScalFieldBotDOUB,NeighBor[S],flag,&(Eta[1 * nx_loc*ny_loc + 1 * 2 * ny_loc]), 1 * nx_fill,MPI_DOUBLE,NeighBor[N],flag,comm2d,&status);		
	
	
}



void write_output_id(int num, int nx_loc, int ny_loc){
	
	int globe_ind;
 	FILE * sortie;
 	char nomfic[256];
 	sprintf(nomfic, "../output/id_%i", num);
 	sortie = fopen(nomfic, "w+");

	for (size_t y = 0 ; y < ny_loc ; y++) {
	for (size_t x = 0 ; x < nx_loc ; x++) {

	int x_global = (xsum + x);
	int y_global = (ysum + y);
				
	globe_ind = (x_global + nx * y_global);
	fprintf(sortie,"%i %i\n",globe_ind,Id[ipx(x,y,nx_loc)/2]);
	
	}	
	}
	
fclose(sortie);

}



void write_output_results(const string& _name, double a1, double a2, double a3, double a4, double a5, int a6, double a7, double a8){

	const char * c = _name.c_str();
 	FILE * sortie;
 	sortie = fopen(c, "a");

	fprintf(sortie,"%g %g %g %g %g %i %g %g\n",a1,a2,a3,a4,a5,a6,a7,a8);
	fclose(sortie);

}


void write_output_sim_param( double a1, double a2, double a3, double a4, double a5, double a6, double a7, double a8, double a9, double a10, double a11, double a12, double a13, double a14, double a15, double a16 ){

FILE * sortie;
sortie = fopen("../output/sim_param", "w+");
fprintf(sortie,"%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16);
fclose(sortie);

}



double post_processing( int idp, double tmp[], double mu_sat , int xls, int xrs, int yls, int yrs ){

	double GP_loc, rhof_loc;
	GP_loc = rhof_loc = 0.0;

	double GP_loc_A, rhof_loc_A, rhof_loc_B;
	GP_loc_A = rhof_loc_A = rhof_loc_B = 0.0;

	double count = 0.0;
    bool res_cond, box_cond;

	for (size_t y = 0 ; y < ny_loc ; y++) {
	for (size_t x = 0 ; x < nx_loc ; x++) {

	int x_global = (xsum + x);
	int y_global = (ysum + y);

    res_cond = ( x_global >= xls ) && ( x_global <  (nx) - xrs ) && ( y_global >= yls ) && ( y_global < (ny) - yrs ) ? true: false;//if true -> continue 
    box_cond = ( x_global >= xls ) && ( x_global <  (nx) - xrs ) && ( y_global >= yls ) && ( y_global < (ny) - yrs ) ? false: true;//if false -> continue 

	if(Id[ipx(x,y,nx_loc)/2] == idp){

	GP_loc += Boltz * Temp * (rho[ipx(x,y,nx_loc)/2] * log(rho[ipx(x,y,nx_loc)/2]) + (Eta[ipx(x,y,nx_loc)/2] - rho[ipx(x,y,nx_loc)/2]) * log(Eta[ipx(x,y,nx_loc)/2] - 
	rho[ipx( x,y,nx_loc )/2]));
		
	for (int l = 0 ; l < nodes[ipx(x, y,nx_loc)/2].nblg ; ++l){
	GP_loc += -wff * rho[nodes[ipx(x, y,nx_loc)/2].ilg[l]/2] * rho[ipx(x,y,nx_loc)/2];
	}
		
	for (int l = HLG ; l < HLG + nodes[ipx(x,y,nx_loc)/2].nbolg ; ++l){
	GP_loc += -wff * rho[nodes[ipx(x,y,nx_loc)/2].ilg[l]/2] * rho[ipx(x,y,nx_loc)/2];
	}
		
	GP_loc -= (mu - mu_sat) * rho[ipx(x,y,nx_loc)/2];
	rhof_loc += rho[ipx(x,y,nx_loc)/2];
	
	}

	if(Id[ipx(x,y,nx_loc)/2] == idp && res_cond){

	GP_loc_A += Boltz * Temp * (rho[ipx(x,y,nx_loc)/2] * log(rho[ipx(x,y,nx_loc)/2]) + (Eta[ipx(x,y,nx_loc)/2] - rho[ipx(x,y,nx_loc)/2]) * log(Eta[ipx(x,y,nx_loc)/2] - rho[ipx(x,y,nx_loc)/2]));
		
	for (int l = 0 ; l < nodes[ipx(x,y,nx_loc)/2].nblg ; ++l){
	GP_loc_A += -wff * rho[nodes[ipx(x,y,nx_loc)/2].ilg[l]/2] * rho[ipx(x,y,nx_loc)/2];
	}
		
	for (int l = HLG ; l < HLG + nodes[ipx(x,y,nx_loc)/2].nbolg ; ++l){
	GP_loc_A += -wff * rho[nodes[ipx(x,y,nx_loc)/2].ilg[l]/2] * rho[ipx(x,y,nx_loc)/2];
	}
		
	GP_loc_A -= (mu - mu_sat) * rho[ipx(x,y,nx_loc)/2];
	rhof_loc_A += rho[ipx(x,y,nx_loc)/2];
	count = count + 1.0;
	}

	if(Id[ipx(x,y,nx_loc)/2] == idp && box_cond) rhof_loc_B += rho[ipx(x,y,nx_loc)/2];

	}
	}
	

	double GP_tot, rhof_tot;
	double GP_tot_A, rhof_tot_A, rhof_tot_B;
	double count_g, countB_g;
	
	MPI_Barrier(comm2d);
	MPI_Allreduce(&GP_loc, &GP_tot, 1, MPI_DOUBLE, MPI_SUM, comm2d);
	MPI_Allreduce(&rhof_loc, &rhof_tot, 1, MPI_DOUBLE, MPI_SUM, comm2d);
	MPI_Allreduce(&GP_loc_A, &GP_tot_A, 1, MPI_DOUBLE, MPI_SUM, comm2d);
	MPI_Allreduce(&rhof_loc_A, &rhof_tot_A, 1, MPI_DOUBLE, MPI_SUM, comm2d);
	MPI_Allreduce(&rhof_loc_B, &rhof_tot_B, 1, MPI_DOUBLE, MPI_SUM, comm2d);
	MPI_Allreduce(&count, &count_g, 1, MPI_DOUBLE, MPI_SUM, comm2d);

	tmp[0] = rhof_tot;//rhof   box+reservoir
	tmp[1] = rhof_tot_A;//rhof   box
	tmp[2] = rhof_tot_B;

	tmp[3] = GP_tot;//local GP box+reservoir
	tmp[4] = GP_tot_A / count_g;//local GP box
  	return GP_tot;
}




void stress_calculation(size_t k, double nn, double MU, double gs, double Temp, tensor & stress, double trho[], double drx[], double dry[], double mu_sat ){

	double p0 = 0.0;
	double p1 = 0.0;
	double cn = 6.0;
	
	p0 += -Boltz * Temp * (	trho[k] * log(trho[k]) + (1.0 - trho[k]) * log(1.0 - trho[k])	) + (mu) * trho[k];
	p0 += (cn * wff / 2.0) * pow(	trho[k],2	);

	stress.xx =  p0 - (cn * pow(gs,2) * wff / 4.0	) * (	drx[k] * drx[k] + dry[k] * dry[k] 	) + (cn * pow(gs,2) / 2.0 * wff * (	drx[k] * drx[k]	)	);
	stress.xy =  cn * pow(gs,2) / 2.0 * wff * drx[k] * dry[k];
	stress.yy =  p0 - (cn * pow(gs,2) * wff / 4.0	) * (	drx[k] * drx[k] + dry[k] * dry[k] 	) + (cn * pow(gs,2) / 2.0 * wff * (	dry[k] * dry[k]	)	);

	stress.xx = stress.xx/pow(gs,2);
	stress.xy = stress.xy/pow(gs,2);
	stress.yy = stress.yy/pow(gs,2);
	stress.yx = stress.xy;
}



void write_output_fluid_stress_data(const string& _name, int num, int rrkk, double NearNeigh, double ChemPot, double LatSpac, double Temp, int idp, int ids, double trho[], double drx[], double dry[], double mu_sat, int xls, int xrs, int yls, int yrs){
	
	tensor sigma, sFres;
	bool res_cond; 
	double count, sigxx,sigyy,sigxy;
	sigxx = sigyy = sigxy = count = 0.0;

	string result = _name + to_string(num);
	result.append("_");
	result = result + to_string(rrkk);
	const char * c = result.c_str();
 	FILE * sortie;
 	sortie = fopen(c, "w+");

	for (size_t y = 0 ; y < ny_loc ; y++) {
	for (size_t x = 0 ; x < nx_loc ; x++) {

	int x_global = (xsum + x);
	int y_global = (ysum + y);

    res_cond = ( x_global >= xls ) && ( x_global <  (nx) - xrs ) && ( y_global >= yls ) && ( y_global < (ny) - yrs ) ? false: true;//if false -> continue 
				
	int globe_ind = (x_global + nx * y_global);
	stress_calculation(ipx(x,y,nx_loc)/2, NearNeigh, ChemPot, LatSpac, Temp, sigma, trho, drx, dry, mu_sat);

	if(res_cond){
	count++;
	sigxx += sigma.xx;
	sigyy += sigma.yy;

	sigxy += sigma.xy;

	}

	if(Id[ipx(x,y,nx_loc)/2] == idp) fprintf(sortie,"%i %.6g\n",globe_ind,(sigma.xx+sigma.yy)/2.0);
	if(Id[ipx(x,y,nx_loc)/2] == ids) fprintf(sortie,"%i %.6g\n",globe_ind,0.0);
			
	}	
	}
	

	fclose(sortie);

	double count_tot,sigxx_tot,sigyy_tot,sigxy_tot;
	
	MPI_Barrier(comm2d);
	MPI_Allreduce(&count, &count_tot, 1, MPI_DOUBLE, MPI_SUM, comm2d);
	MPI_Allreduce(&sigxx, &sigxx_tot, 1, MPI_DOUBLE, MPI_SUM, comm2d);
	MPI_Allreduce(&sigyy, &sigyy_tot, 1, MPI_DOUBLE, MPI_SUM, comm2d);
	MPI_Allreduce(&sigxy, &sigxy_tot, 1, MPI_DOUBLE, MPI_SUM, comm2d);

	sFres.xx = sigxx_tot / count_tot;	
	sFres.yy = sigyy_tot / count_tot;
	sFres.xy = sigxy_tot / count_tot;

	if (rrkk == 0){
	FILE * sortieA;
	sortieA = fopen("../output/sigF_res","a");
	fprintf(sortieA,"%g %g %g\n",sFres.xx,sFres.yy,sFres.xy);
	fclose(sortieA);
	}

}



void write_output_data(const string& _name, int num, int rrkk, double data[] ){
	
	string result = _name + to_string(num);
	result.append("_");
	result = result + to_string(rrkk);
	const char * c = result.c_str();
 	FILE * sortie;
 	sortie = fopen(c, "w+");

	for (size_t y = 0 ; y < ny_loc ; y++) {
	for (size_t x = 0 ; x < nx_loc ; x++) {

	int x_global = (xsum + x);
	int y_global = (ysum + y);

	int globe_ind = ( x_global + nx * y_global );

	fprintf( sortie,"%i %g\n",globe_ind, data[ipx(x,y,nx_loc)/2]);
			
	}	
	}

fclose(sortie);

}





void  reservoir_creation( int idincl, double ideta, int xls, int xrs, int yls, int yrs ){

   bool res_cond;

   for (size_t y = 0 ; y < ny_loc ; y++) {
   for (size_t x = 0 ; x < nx_loc ; x++) {

	int x_global = (xsum + x);
	int y_global = (ysum + y);
		
   res_cond = ( x_global >= xls ) && ( x_global <  (nx) - xrs ) && ( y_global >= yls ) && ( y_global < (ny) - yrs ) ? false: true;//if true -> continue 

   if (res_cond){
   Id[ipx(x,y,nx_loc)/2] = idincl;
   Eta[ipx(x,y,nx_loc)/2] = ideta;
   }

   }
   }

	MPI_Barrier(comm2d);
	/** scalar fields **/
	/** to left from right **/
	MPI_Sendrecv (&(Id[0]),1,CtScalField,NeighBor[W],flag,&(Id[nx_loc*ny_loc*1]),1 * ny_loc,MPI_INT,NeighBor[E],flag,comm2d,&status);
	/** to right from left **/
	MPI_Sendrecv (&(Id[1 * nx_loc - 1]),1,CtScalField,NeighBor[E],flag,&(Id[1 * nx_loc * ny_loc + 1 * ny_loc]), 1 * ny_loc, MPI_INT,NeighBor[W],flag,comm2d,&status);	
	/** to top from bottom **/
	MPI_Sendrecv (&(Id[0]),1,RtScalFieldTop,NeighBor[N],flag,&(Id[1*nx_loc*ny_loc + 2 * 1 * ny_loc + 1 * nx_loc + 2 * 1]),(1 * nx_loc + 2 * 1),MPI_INT,NeighBor[S],flag,comm2d,&status);
	/** to bottom from top **/
	MPI_Sendrecv (&(Id[0]),1,RtScalFieldBot,NeighBor[S],flag,&(Id[1 * nx_loc*ny_loc + 1 * 2 * ny_loc]), 1 * nx_fill,MPI_INT,NeighBor[N],flag,comm2d,&status);
	
	/** scalar fields **/
	/** to left from right **/
	MPI_Sendrecv (&(Eta[0]),1,CtScalFieldDOUB,NeighBor[W],flag,&(Eta[nx_loc*ny_loc*1]),1 * ny_loc,MPI_DOUBLE,NeighBor[E],flag,comm2d,&status);
	/** to right from left **/
	MPI_Sendrecv (&(Eta[1 * nx_loc - 1]),1,CtScalFieldDOUB,NeighBor[E],flag,&(Eta[1 * nx_loc * ny_loc + 1 * ny_loc]), 1 * ny_loc, MPI_DOUBLE,NeighBor[W],flag,comm2d,&status);	
	/** to top from bottom **/
	MPI_Sendrecv (&(Eta[0]),1,RtScalFieldTopDOUB,NeighBor[N],flag,&(Eta[1*nx_loc*ny_loc + 2 * 1 * ny_loc + 1 * nx_loc + 2 * 1]),(1 * nx_loc + 2 * 1),MPI_DOUBLE,NeighBor[S],flag,comm2d,&status);
	/** to bottom from top **/
	MPI_Sendrecv (&(Eta[0]),1,RtScalFieldBotDOUB,NeighBor[S],flag,&(Eta[1 * nx_loc*ny_loc + 1 * 2 * ny_loc]), 1 * nx_fill,MPI_DOUBLE,NeighBor[N],flag,comm2d,&status);

}

void correct_rho(double tnrho[], double lb, double ub){

	for (size_t y = 0 ; y < ny_loc ; y++) {
	for (size_t x = 0 ; x < nx_loc ; x++) {
	double lrho = rho[ipx(x,y,nx_loc)/2];
	tnrho[ipx(x,y,nx_loc)/2] = lrho;
	if(lrho < lb) tnrho[ipx(x,y,nx_loc)/2] = 1e-20;
	if(lrho > ub) tnrho[ipx(x,y,nx_loc)/2] = 0.9999999;

	}	
	}
	
	/** scalar fields **/
	/** to left from right **/
	MPI_Sendrecv (&(tnrho[0]),1,CtScalFieldDOUB,NeighBor[W],flag,&(tnrho[nx_loc*ny_loc*1]),1 * ny_loc,MPI_DOUBLE,NeighBor[E],flag,comm2d,&status);
	/** to right from left **/
	MPI_Sendrecv (&(tnrho[1 * nx_loc - 1]),1,CtScalFieldDOUB,NeighBor[E],flag,&(tnrho[1 * nx_loc * ny_loc + 1 * ny_loc]), 1 * ny_loc, MPI_DOUBLE,NeighBor[W],flag,comm2d,&status);	
	/** to top from bottom **/
	MPI_Sendrecv (&(tnrho[0]),1,RtScalFieldTopDOUB,NeighBor[N],flag,&(tnrho[1*nx_loc*ny_loc + 2 * 1 * ny_loc + 1 * nx_loc + 2 * 1]),(1 * nx_loc + 2 * 1),MPI_DOUBLE,NeighBor[S],flag,comm2d,&status);
	/** to bottom from top **/
	MPI_Sendrecv (&(tnrho[0]),1,RtScalFieldBotDOUB,NeighBor[S],flag,&(tnrho[1 * nx_loc*ny_loc + 1 * 2 * ny_loc]), 1 * nx_fill,MPI_DOUBLE,NeighBor[N],flag,comm2d,&status);	

}


int main(int argc, char *argv[]){

	readSimParam();

	/** adding reservoirs to the sides **/
	int nx_lt_shift = Lx_lt_shift / gridStep + 1;
	int nx_rt_shift = Lx_rt_shift / gridStep + 1;
	int ny_lt_shift = Ly_lt_shift / gridStep + 1;
	int ny_rt_shift = Ly_rt_shift / gridStep + 1;

	if (Lx_lt_shift == 0) nx_lt_shift = 0;
	if (Lx_rt_shift == 0) nx_rt_shift = 0;
	if (Ly_lt_shift == 0) ny_lt_shift = 0;
	if (Ly_rt_shift == 0) ny_rt_shift = 0;

	/** total number of nodes **/
	nx = Lnx / gridStep + 1 + nx_lt_shift + nx_rt_shift;
	ny = Lny / gridStep + 1 + ny_lt_shift + ny_rt_shift;
	
	/** MPI stuff **/
	MPI_Comm comm;
	int iter;
	double fret;
	int nproc, rank;     
    int reorder = 0;
	
	/** dimensions, degree of freedom, periodicity **/
	int dims[ndims];
	int periods[ndims];
	int coord[ndims];

	MPI_Init(&argc, &argv);
    comm = MPI_COMM_WORLD;
    MPI_Comm_size(comm,&nproc);
    MPI_Comm_rank(comm,&rank);

	if(rank == 0) cout<<"nx: "<<nx<<" ny: "<<ny<<endl;
	
	lengthX = new int[nxp];
	lengthY = new int[nyp];
		
	/** decomposing domain according to number of processors **/
	decompose1d(lengthX,nx,nxp);
	decompose1d(lengthY,ny,nyp);
	
	if(rank == 0){
	for (size_t i = 0 ; i < nxp ; i++) cout<<"lengthX[rank]: "<<lengthX[i]<<endl;
	for (size_t i = 0 ; i < nyp ; i++) cout<<"lengthY[rank]: "<<lengthY[i]<<endl;
	}
	

	periods[0] = 1;
    periods[1] = 1;
	
    dims[0] = nxp;
    dims[1] = nyp;
	
    MPI_Cart_create(comm, ndims, dims, periods, reorder, &comm2d);
	MPI_Cart_get(comm2d , ndims , dims , periods , coord);

	NeighBor[0] = MPI_PROC_NULL;
    NeighBor[1] = MPI_PROC_NULL;
    NeighBor[2] = MPI_PROC_NULL;
    NeighBor[3] = MPI_PROC_NULL;

    MPI_Cart_shift(comm2d,0,1,&NeighBor[W],&NeighBor[E]);
  	MPI_Cart_shift(comm2d,1,1,&NeighBor[S],&NeighBor[N]);
	
	nx_loc = lengthX[coord[0]];
	ny_loc = lengthY[coord[1]];
		
	nx_fill = (nx_loc + 2);
	ny_fill = (ny_loc + 2);
		
	xsum = 0;
	ysum = 0;
	int xx,yy;
	
	for(xx = 0 ; xx < coord[0] ; xx++) xsum += lengthX[xx];
	for(yy = 0 ; yy < coord[1] ; yy++) ysum += lengthY[yy];
	

		
	p = new double[nx_fill * ny_fill * ndof];
	sigF = new double[nx_fill * ny_fill * ndof];
		
	for(size_t i = 0 ; i < (nx_fill * ny_fill * ndof) ; i++){
	p[i] = 0.0;
	sigF[i] = 0.0;
	}
	
	
	int count = 0;
	for (int j = 0 ; j < ny_loc ; j++) {
	for (int i = 0 ; i < nx_loc ; i++) {
	int x = xsum + i;
	int y = ysum + j;
	p[ipx(i,j,nx_loc)] = x * gridStep;
	p[ipx(i,j,nx_loc) + 1] = y * gridStep;
	// Id[ipx(i,j,nx_loc)/2] = ipx(x,y,nx)/2;
	count++;
	}
	}
	
	Id = new int[nx_fill * ny_fill];
	Eta = new double[nx_fill * ny_fill];	
	
		
	MPI_Type_vector(ny_loc, ndims, ndims * nx_loc, MPI_DOUBLE, &CtVecField);/** column type vector field **/  
	MPI_Type_commit(&CtVecField);
	
	MPI_Type_vector(ny_loc, 1, 1 * nx_loc, MPI_INT, &CtScalField);/** column type vector field **/  
	MPI_Type_commit(&CtScalField);
	
	MPI_Type_vector(ny_loc, 1, 1 * nx_loc, MPI_DOUBLE, &CtScalFieldDOUB);/** column type vector field **/  
	MPI_Type_commit(&CtScalFieldDOUB);	
	
	/** this is for double type vector field **/
	int BlockLengthTopVecField[3] = {ndims,ndims * nx_loc,ndims};/** block length top field **/
	int DispTopVecField[3] = {ndims*nx_loc*ny_loc + ndims*ny_loc, 0 , ndims*nx_loc*ny_loc};/** disp top field **/
	// MPI_Datatype RtVecFieldTop;/** row type vector field top **/
	MPI_Type_indexed(3,BlockLengthTopVecField,DispTopVecField,MPI_DOUBLE,&RtVecFieldTop);
	MPI_Type_commit(&RtVecFieldTop);
	
	int BlockLengthBotVecField[3] = {2,ndims*nx_loc,2};/** block length bottom field **/
	int DispBotVecField[3] = {ndims*nx_loc*ny_loc + ndims* 2 * ny_loc - ndims, ndims * nx_loc*(ny_loc-1) , ndims* nx_loc*ny_loc + ndims * ny_loc - ndims};/** disp bottom field **/
	// MPI_Datatype RtVecFieldBot;/** row type vector field bottom **/
	MPI_Type_indexed(3,BlockLengthBotVecField,DispBotVecField,MPI_DOUBLE,&RtVecFieldBot);
	MPI_Type_commit(&RtVecFieldBot);
	
	/** this is for int type scalar field **/
	int BlockLengthTopScalField[3] = {1,1 * nx_loc,1};/** block length top field **/
	int DispTopScalField[3] = {1*nx_loc*ny_loc + 1*ny_loc, 0 , 1*nx_loc*ny_loc};/** disp top field **/
	// MPI_Datatype RtScalFieldTop;/** row type vector field top **/
	MPI_Type_indexed(3,BlockLengthTopScalField,DispTopScalField,MPI_INT,&RtScalFieldTop);
	MPI_Type_commit(&RtScalFieldTop);
	
	int BlockLengthBotScalField[3] = {1,1*nx_loc,1};/** block length bottom field **/
	int DispBotScalField[3] = {1*nx_loc*ny_loc + 1 * 2 * ny_loc - 1, 1 * nx_loc*(ny_loc-1) , 1* nx_loc*ny_loc + 1 * ny_loc - 1};/** disp bottom field **/
	// MPI_Datatype RtScalFieldBot;/** row type vector field bottom **/
	MPI_Type_indexed(3,BlockLengthBotScalField,DispBotScalField,MPI_INT,&RtScalFieldBot);
	MPI_Type_commit(&RtScalFieldBot);	
	
	/** this is for double type scalar field **/
	// MPI_Datatype RtScalFieldTopDOUB;/** row type vector field top **/
	MPI_Type_indexed(3,BlockLengthTopScalField,DispTopScalField,MPI_DOUBLE,&RtScalFieldTopDOUB);
	MPI_Type_commit(&RtScalFieldTopDOUB);	

	// MPI_Datatype RtScalFieldBotDOUB;/** row type vector field bottom **/
	MPI_Type_indexed(3,BlockLengthBotScalField,DispBotScalField,MPI_DOUBLE,&RtScalFieldBotDOUB);
	MPI_Type_commit(&RtScalFieldBotDOUB);		
	
	
	
	int IDP = 0;
	int IDS = 1;
	
	save_MPI_info( rank,nx_loc,ny_loc,coord[0],coord[1],xsum,lengthX[coord[0]],ysum,lengthY[coord[1]] );	
	init(gridStep, ndof, nx_loc, ny_loc, xsum, ysum);
	
	for(size_t i = 0 ; i < (nx_fill * ny_fill) ; i++){
	Id[i] = IDP;	
	Eta[i] = 1.0;
	}
	
	rho = new double[nx_fill * ny_fill];
	double* drhoX;
	double* drhoY;
	double* nrho;

	drhoX = new double[nx_loc * ny_loc];
	drhoY = new double[nx_loc * ny_loc];
	nrho = new double[nx_fill * ny_fill];
	
	for (size_t i = 0 ; i < (nx_loc * ny_loc) ; i++){
	drhoX[i] = 0.0;
	drhoY[i] = 0.0;
	nrho [i] = 0.0;
	}
	
	
    // Get morphology, init. position and init. velocity input files
    char tempfname[100];    
    sprintf(tempfname, "../input/trial_0/morphIDs.dat");
    string file_morph = tempfname;
    sprintf(tempfname, "../input/trial_0/positionsEquil.dat");
    string file_pos = tempfname;
    sprintf(tempfname, "../input/shapes/");
    string morph_dir = tempfname;
	
    // Generate grains
    vector<Grain2d> grains = generateGrainsFromFiles(file_morph, morph_dir, file_pos, gridStep);
    // size_t ngrains = grains.size();	
	
	Import_Structure_arb_particles_fromLS(grains, IDS, 0.0,nx_loc, ny_loc, xsum, ysum, Lx_lt_shift, Ly_lt_shift,lengthX[coord[0]],lengthY[coord[1]], rank);
	smooth_particles(IDP, IDS);
    	// Import_Structure_spherical_particles("input/FinalCoor.txt", IDS, 0.0, nu_part,nx_loc, ny_loc, xsum, ysum, Lx_lt_shift, Ly_lt_shift);

	reservoir_creation(IDP,1.0, nx_lt_shift,nx_rt_shift,ny_lt_shift,ny_rt_shift);
	write_output_id(rank,nx_loc,ny_loc);
	// write_output_data("eta_",0,rank,Eta);	
	
	
	double NN = 4.0; // number of nearest neighbors 2D simple cubic lattice 
	double surf_ten = 0.072;// N*m^-1
	Boltz = 1.38064852 * pow(10,-23);//m^2*kg*s^-2*K^-1
    // Boltz = 8.3144621;//J*mol^-1*K^-1
	wff = Boltz * Tcw * 4.0 / NN ;//J*mol^-1//m^2*kg*s^-2=N.m=J
	wmf = yvar * wff;//J*mol^-1//m^2*kg*s^-2=N.m=J
  	Temp = Tstar * wff / Boltz;
	
	double mu_sat = -wff * NN / 2.0;
	double vp0 = exp(mu_sat/Boltz/Temp);
	double a_lg = sqrt(wff/2.0/surf_ten);//m	
	
	if (rank == 0) write_output_sim_param(Boltz,wff,wmf,yvar,Tcw,Tstar,Temp,mu_sat,vp0,gridStep,Lnx,Lny, Lx_lt_shift, Lx_rt_shift, Ly_lt_shift, Ly_rt_shift);
	
	// non-dimensional values 
	Temp = Boltz*Temp/wff;
	wmf = wmf/wff;
	Boltz = Boltz/Boltz;
	mu_sat = mu_sat / wff;
	wff = wff / wff;
	
	double h,vp;
	double* res;
	double* NoPore;
	res = new double[5];
	NoPore = new double[3];		
	
	h = 0.001;
	mu = mu_sat + Boltz * Temp * log(h);
	double rhoInit = exp((1/Boltz/Temp) * mu);
	initialize_rho(0.001,IDP,nx_loc,ny_loc);
    // initialize_reservoir(0.999,nx_lt_shift,nx_rt_shift,ny_lt_shift,ny_rt_shift,nz_dn_shift,nz_up_shift);	

	count = 0;
	for (int i = 1 ; i < 101; ++i){
	double step = static_cast<double>(i) / 100.0;
	h = step;
	mu = mu_sat + Boltz * Temp * log(h);
	vp = exp(mu/Boltz/Temp);
	mu = mu / wff;	
		
	DFT_ObjFunc(obj_tol,Temp,mu,nx_loc,ny_loc,nx_fill,ny_fill,rank,nx_lt_shift,nx_rt_shift,ny_lt_shift,ny_rt_shift);	
	double GrandPot = post_processing(IDP,res, 0, nx_lt_shift,nx_rt_shift,ny_lt_shift,ny_rt_shift);
	double phi = vol_frac(NoPore,IDP,rank,nx_lt_shift,nx_rt_shift,ny_lt_shift,ny_rt_shift);
	
	correct_rho(nrho,rho_thresh_lb, rho_thresh_ub);
	density_gradient(IDP,IDS,nrho,drhoX,drhoY,gridStep);

	if (rank == 0) write_output_results("../output/sim_data_tot",res[0]/NoPore[0],res[1]/NoPore[1],res[2]/NoPore[2],mu,h,vp/vp0,res[3]/NoPore[0],res[4]);
	if (rank == 0) cout<<"step = " <<count<<" mu: "<<mu<<endl;
		
	write_output_fluid_stress_data("../output/f_",count, rank, NN, mu, gridStep,  Temp, IDP, IDS, nrho, drhoX, drhoY,mu_sat, nx_lt_shift,nx_rt_shift,ny_lt_shift,ny_rt_shift);
	write_output_data("../output/sat_",count,rank,rho);
		
	count++;
	}	
	
	
delete [] rho;
// delete [] nrho;
// delete [] res;
// delete [] NoPore;
delete [] drhoX;
delete [] drhoY;
delete [] nrho;
delete [] Id;
delete [] Eta;
delete [] p;
delete [] sigF;

	
MPI_Comm_free(&comm2d);
MPI_Finalize();

return 0;
	
}

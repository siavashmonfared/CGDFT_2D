#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <set>

using namespace std;

const double pi = 3.1415926535897; 

int Lnx,Lny;
int nxp, nyp;
int npart;

double gridStep, obj_tol; 
double Tstar, Tcw, yvar;
double rho_thresh_ub, rho_thresh_lb, rho_thresh;
double Lx_lt_shift, Lx_rt_shift, Ly_lt_shift, Ly_rt_shift;

int nx, ny;

/*
double fp, fs;
vector<int>VOL;
size_t nx, ny;	   
size_t type, nit;
double strain, ftol, Rpore; 	
double XC, YC, obj_tol;
double rho_thresh;
double wff,wmf,mu,Temp, Boltz,yvar;
double NP,NS; 
*/








int ipx(int x, int y, int xsize)
{
	return (x + xsize * y);
}



void  Import_data( const string& _name, int num, double tmp[], size_t xl, size_t xu, size_t yl, size_t yu ){
	
	double local_var;
	string result = _name + to_string(num);
	fstream file(result);
	
	size_t count = 0;
				
		for (size_t y = 0 ; y < ny ; y++) {
			for (size_t x = 0 ; x < nx ; x++) { 
			file >> local_var;
			if( x >= xl && x < xu && y >= yl && y < yu ){
			tmp[count] = local_var;
			count++;
			}
			}
		}
}


void  Import_press(const string& _name, int num, double tmp[], size_t xl, size_t xu, size_t yl, size_t yu, double npress){
	
	double local_var;
	string result = _name + to_string(num);
	fstream file(result);
	
	size_t count = 0;
				
		for (size_t y = 0 ; y < ny ; y++) {
			for (size_t x = 0 ; x < nx ; x++) { 
			file >> local_var;
			if( x >= xl && x < xu && y >= yl && y < yu ){
			tmp[count] = local_var - npress;
			count++;
			}
			}
		}
	
}


int  Import_id( const string& _name, int tmpID[], int IDP, size_t xl, size_t xu, size_t yl, size_t yu ){

	int tag;
	fstream file(_name);

	size_t count = 0;
	int NoPo = 0;
					
		for (size_t y = 0 ; y < ny ; y++) {
			for (size_t x = 0 ; x < nx ; x++) { 
			file >> tag;
			if( x >= xl && x < xu && y >= yl && y < yu ){
			tmpID[count] = tag;
			count++;
			if ( tag == IDP ) NoPo++;
			}
			}
		}
	
return NoPo;
}

/*
size_t volfrac_calculation(int tid[]) 
{
	cout<<"computing volume fractions"<<endl;
	size_t N = nx * ny;
	double tot = N;
	NP = 0.0;
	NS = 0.0;

	int f;
	for (size_t z = 0 ; z < N ; ++z) { 
	f = tid[z];
	if(f==0){NP++;};
	if(f>0){NS++;};
	}

	fp = NP/tot;
	fs = NS/tot;

	cout<<"fp: "<<fp<<" fs: "<<fs<<endl;
	size_t nsi = NS;
	return nsi;
}
*/

double find_min(const string& _name, double tmp[], size_t tot){

	const char * c = _name.c_str();
	
	double small = tmp[0];
	for (size_t i = 0 ; i < tot ; i++){
	if(tmp[i] < small) small = tmp[i];
	}

	FILE * sortie;
	sortie = fopen(c,"a");
	fprintf(sortie,"%g\n",small);
	fclose(sortie);

	return small;

}

double find_max(const string& _name, double tmp[], size_t tot){
	
   	const char * c = _name.c_str();

	double large = tmp[0];
	for (size_t i = 0 ; i < tot ; i++){
	if(tmp[i] > large) large = tmp[i];
	}

	FILE * sortie;
	sortie = fopen(c,"a");
	fprintf(sortie,"%g\n",large);
	fclose(sortie);

	return large;
}


void histogram(const string& _name, int num, double tmp[], size_t tot, double val_min, double val_max, double binwidth){
	
   string result = _name + to_string(num);
   string resultOUT = result;		
   const char * c = resultOUT.c_str();
	
	if(val_min >= 0.0 && val_max >= 0.0){
	
	int nbins;
	int* hist;
	double* intv;
	nbins = (int)ceil ( (val_max) / binwidth);
	hist = new int[nbins];
	intv = new double[nbins];
	double intv0 = val_min + (binwidth)/2.0;
	for(int i = 0 ; i < nbins ; i++) hist[i] = 0;
	for(int i = 0 ; i < nbins ; i++){
	intv0 += binwidth;
	intv[i] = intv0;
	}

	for (size_t i = 0 ; i < tot ; i++){
    int bucket = (int)floor(tmp[i] / binwidth);
    hist[bucket - (int)floor(val_min/binwidth) ] += 1;
	}	

	FILE * sortie;
	sortie = fopen(c,"w+");
	for(int i = 0 ; i < nbins ; i++) fprintf(sortie,"%g %i\n",intv[i],hist[i]);
	fclose(sortie);
	delete [] hist;
	delete [] intv;
	
	}
	
	if(val_min < 0.0 && val_max < 0.0){
		
	int nbins;
	int* hist;
	double* intv;
	nbins = (int)ceil ( (-val_min) / binwidth);
	hist = new int[nbins];
	intv = new double[nbins];
	double intv0 = val_min + (binwidth)/2.0;
	for(int i = 0 ; i < nbins ; i++) hist[i] = 0;	
	for(int i = 0 ; i < nbins ; i++){
	intv0 += binwidth;
	intv[i] = intv0;
	}

	for (size_t i = 0 ; i < tot ; i++){
    int bucket = (int)floor(-tmp[i] / binwidth);
    hist[bucket - (int)floor(-val_max/binwidth)] += 1;
	}
	FILE * sortie;
	sortie = fopen(c,"w+");
	for(int i = 0 ; i < nbins ; i++) fprintf(sortie,"%g %i\n",intv[i],hist[nbins - 1 - i]);
	fclose(sortie);
	delete [] hist;
	delete [] intv;
	}
	
	if(val_min < 0.0 && val_max >= 0.0){

	int nbins_pos, nbins_neg;
	int* hist_pos; 
	int* hist_neg;
	double* intv;
   	nbins_pos = (int)ceil( (val_max)/binwidth);
   	nbins_neg = (int)ceil( (-val_min)/binwidth);

		
   	hist_pos = new int[nbins_pos];	
   	hist_neg = new int[nbins_neg];
	intv = new double[nbins_pos + nbins_neg];
   	for(int i = 0 ; i < nbins_pos ; i++) hist_pos[i] = 0;
   	for(int i = 0 ; i < nbins_neg ; i++) hist_neg[i] = 0;
	double intv0 = val_min + (binwidth)/2.0;
	for(int i = 0 ; i < (nbins_pos + nbins_neg) ; i++){
	intv0 += binwidth;
	intv[i] = intv0;

	}

	for (size_t i = 0 ; i < tot ; i++){
	if(tmp[i] >= 0){
    	int bucket = (int)floor(tmp[i] / binwidth);

    	hist_pos[bucket] += 1;
	}
	if(tmp[i] < 0){
    	int bucket = (int)floor(-tmp[i] / binwidth);

    	hist_neg[bucket] += 1;
	}
	}
	FILE * sortie;
	sortie = fopen(c,"w+");
	for(int i = 0 ; i < nbins_neg ; i++) fprintf(sortie,"%g %i\n",intv[i],hist_neg[nbins_neg - 1 - i]);
	for(int i = 0 ; i < nbins_pos ; i++) fprintf(sortie,"%g %i\n",intv[i + nbins_neg],hist_pos[i]);
	fclose(sortie);

	delete [] hist_pos;
	delete [] hist_neg;
	delete [] intv;	
	}
		
}



void write_output(const string& _name, int npart, int nbins, double tbins[], double tcounts[], double rho){

   const char * c = _name.c_str();
   FILE * sortie;
   sortie = fopen(c, "w+");	

   for (size_t i = 0 ; i < nbins ; i++){
   double interval = 0.5 * (tbins[i] + tbins[i+1]);
   //double freq = tcounts[i]/npart/(4.0*pi*pow(interval,2)*(tbins[i+1]-tbins[i]))/rho;
   double freq = tcounts[i]/(4.0*pi*pow(interval,2)*(tbins[i+1]-tbins[i]))/rho;
   fprintf(sortie,"%.7e %.7e %.7e %.7e %.7e\n",tbins[i],tbins[i+1], interval, freq,rho);
   }
   fclose(sortie);
}

void compute_cumulants(const string& _name,int lamb, double tmp[], int tid[] , int tot, int IDP) 
{
	
	string result = _name;
	string resultOUT = result;
	result.append("_");
	string resultIN = result + to_string(lamb);
	
	const char * cA = resultIN.c_str();

	double m1, m2, count, ExprA, ExprB, ExprC;
	m1 = m2 = count = ExprA = ExprB = ExprC = 0.0;

	vector<double> tmp_storage;
	
	for (int i = 0 ; i < tot ; i++){
	if (tid[i] == IDP){
	count++;
	m1 += tmp[i];
	m2 += pow(tmp[i],2);
	tmp_storage.push_back(tmp[i]);
	}
	}
	
	m1 = m1/(count);
	m2 = m2/(count);
	
	for (int i = 0 ; i < tmp_storage.size() ; i++){
	ExprA += pow((tmp_storage[i] - m1),2);
	ExprB += pow((tmp_storage[i] - m1),3);
	ExprC += pow((tmp_storage[i] - m1),4);
	}

	double Mean = m1;//first cumulant
	double Std = sqrt(ExprA/(count - 1));
	//double Var = m2 - pow(m1,2);//second cumulant
	double Var = pow(Std,2);
	double Skew = (ExprB/count) / pow(sqrt(ExprA/count),3);
	double Kurt = (ExprC/count) / pow(ExprA/count,2);
	
	FILE * sortie;
	sortie = fopen(cA,"a");
	fprintf(sortie,"%g %g %g %g\n",Mean, Var, Skew, Kurt);
	fclose(sortie);	
}

void compute_cumulants_f(const string& _name,int lamb, double tmp[], int tid[] , int tot, int IDP, int pindex [] , int nopo) 
{
	double* tp;
	tp = new double[nopo];

	for (int i = 0 ; i < nopo ; i++) tp[i] = tmp[pindex[i]];

	string result = _name;
	string resultOUT = result;
	result.append("_");
	string resultIN = result + to_string(lamb);
	
	const char * cA = resultIN.c_str();

	double m1, m2, count, ExprA, ExprB, ExprC;
	m1 = m2 = count = ExprA = ExprB = ExprC = 0.0;
	
	for (int i = 0 ; i < nopo ; i++){
	count++;
	m1 += tp[i];
	m2 += pow(tp[i],2);
	}
	
	m1 = m1/(count);
	m2 = m2/(count);
	
	for (int i = 0 ; i < nopo ; i++){
	ExprA += pow((tp[i] - m1),2);
	ExprB += pow((tp[i] - m1),3);
	ExprC += pow((tp[i] - m1),4);
	}

	double Mean = m1;//first cumulant
	double Std = sqrt(ExprA/(count - 1));
	//double Var = m2 - pow(m1,2);//second cumulant
	double Var = pow(Std,2);
	double Skew = (ExprB/count) / pow(sqrt(ExprA/count),3);
	double Kurt = (ExprC/count) / pow(ExprA/count,2);
	
	FILE * sortie;
	sortie = fopen(cA,"a");
	fprintf(sortie,"%g %g %g %g\n",Mean, Var, Skew, Kurt);
	fclose(sortie);	

	delete [] tp;
}




int Mask_id(int tid[] , int tid_mask[] , int tot , int IDP , int IDS) 
{
	
	bool xn,yn,zn,xp,yp,zp;
	int* n;
	n = new int[4];
	for (int i = 0 ; i < 4 ; i++) n[i] = 0;
	int count = 0;
					
	for (int y = 0 ; y < ny ; y++) {
	for (int x = 0 ; x < nx ; x++) {

	xn = (x == nx - 1) ? false : true;
	yn = (y == ny - 1) ? false : true;

	xp = (x == 0) ? false : true;
	yp = (y == 0) ? false : true;

	int p1 = tid[ipx(x,y,nx)];

	if (xn) n[0] = tid[ipx(x+1,y,nx)];
	if (yn) n[1] = tid[ipx(x,y+1,nx)];

	if (xp) n[2] = tid[ipx(x-1,y,nx)];
	if (yp) n[3] = tid[ipx(x,y-1,nx)];



	int sn = 0;
	for (int i = 0 ; i < 4 ; i++) sn += n[i];
	if (p1 == IDP && sn > 0){
	tid_mask[ipx(x,y,nx)] = IDS;
	count++;
	}
		
	}
	}
	
		
	delete [] n;
	return count;
}


void gr(const string& _name, int lamb, double tdata[], int tid[], double xdata[], double ydata[], double tbins[], double tc_ll[], double tc_gg[], double tc_gl[], double tc_lg[], int tot, int IDP, double thresh, int nbins, double rcut, int volcount, int pvol, int nll, int ngg){
	
	double r;

	for (int i = 0 ; i < tot-1 ; i++){
	for (int j = i ; j < tot ; j++){

	double drx = xdata[i]-xdata[j];
	double dry = ydata[i]-ydata[j];

	r = sqrt( pow(drx,2)+pow(dry,2) );

	if (tid[i] == IDP && tid[j] == IDP && tdata[i] >= thresh && tdata[j] >= thresh){//liquid-liquid
	for (int k = 0 ; k < nbins ; k++){
	if (r <= rcut && r  > 0.0 && r >= tbins[k] && r < tbins[k+1]){
	tc_ll[k] += 1.0;
	}
	}
	}

	if (tid[i] == IDP && tid[j] == IDP && tdata[i] < thresh && tdata[j] < thresh){//gas-gas
	for (int k = 0 ; k < nbins ; k++){
	if (r <= rcut && r  > 0.0 && r >= tbins[k] && r < tbins[k+1]){
	tc_gg[k] += 1.0;
	}
	}
	}

	if (tid[i] == IDP && tid[j] == IDP && tdata[i] >= thresh && tdata[j] < thresh){//liquid-gas
	for (int k = 0 ; k < nbins ; k++){
	if (r <= rcut && r  > 0.0 && r >= tbins[k] && r < tbins[k+1]){
	tc_lg[k] += 1.0;
	}
	}
	}

	if (tid[i] == IDP && tid[j] == IDP && tdata[i] < thresh && tdata[j] >= thresh){//gas-liquid
	for (int k = 0 ; k < nbins ; k++){
	if (r <= rcut && r  > 0.0 && r >= tbins[k] && r < tbins[k+1]){
	tc_gl[k] += 1.0;
	}
	}
	}
	}
	}


	string result = _name;
	string resultOUT = result;
	result.append("_");

	string resultLL = result + to_string(lamb);
	resultLL.append("_ll");

	string resultGG = result + to_string(lamb);
	resultGG.append("_gg");

	string resultGL = result + to_string(lamb);
	resultGL.append("_gl");

	string resultLG = result + to_string(lamb);
	resultLG.append("_lg");

	write_output(resultLL,nll,nbins,tbins,tc_ll,1.*nll/1.*pvol);
	write_output(resultGG,ngg,nbins,tbins,tc_gg,1.*ngg/1.*pvol);
	write_output(resultGL,ngg,nbins,tbins,tc_gl,1.*ngg/1.*pvol);
	write_output(resultLG,nll,nbins,tbins,tc_lg,1.*nll/1.*pvol);
	
}

void  Export(const string& _name, int num, double tmp[], size_t tot){
	
	
   string result = _name + to_string(num);
   string resultOUT = result;		
   const char * c = resultOUT.c_str();
   FILE * sortie;
   sortie = fopen(c, "w+");	
   for (size_t k = 0 ; k < (tot) ; k++) fprintf(sortie,"%g\n",tmp[k]);
   fclose(sortie);
}

void  Export_coordinates(const string& _name, int npart, double xc[], double yc[], double rad){
	
	
   const char * c = _name.c_str();
   FILE * sortie;
   sortie = fopen(c, "w+");	
   for (int k = 0 ; k < (npart) ; k++) fprintf(sortie,"%g %g %g\n",xc[k],yc[k],rad);
   fclose(sortie);
}



void coarse_graining_func(	const string& _name, int num, double lambda, double L0, double data[], int id[], int IDP)	{

	int ilambda = lambda;
	int iL0 = L0;
	int tmp = (ilambda-1) * gridStep;
		
	FILE * sortie;
	if (iL0 % tmp == 0.0){
	int nvol = L0 / ((lambda-1)*gridStep) + 1;
	int* LengthX;
	int* LengthY;

	int* intX;
	int* intY;

	int tmp = lambda;
	LengthX = new int[nvol];
	LengthY = new int[nvol];

	intX = new int[tmp];
	intY = new int[tmp];

	int lbx, ubx, lby, uby;
		
	double* val_data;
	val_data = new double[ilambda * ilambda * ilambda];

	int* id_data;
	id_data = new int[ilambda * ilambda * ilambda];		
		
	for (size_t i = 0 ; i < nvol ; i++){
	//for overlapping walls
	LengthX [i] = i * (lambda - 1);
	LengthY [i] = i * (lambda - 1);

	/*
	//for non-overlapping walls
	LengthX [i] = i * (lambda );
	LengthY [i] = i * (lambda );
	LengthZ [i] = i * (lambda );*/
	}
		
	int vol_count = 0;

	for (size_t j = 0 ; j < nvol-1 ; j++){
	lby =  LengthY[j];
	uby =  LengthY[j + 1];
	for (size_t k = 0 ; k < nvol-1 ; k++){
	lbx =  LengthX[k];
	ubx =  LengthX[k + 1];

	for (size_t ii = 0 ; ii < tmp ; ii++){
	intX[ii] = lbx + ii;	
	intY[ii] = lby + ii;	
	}

	size_t count = 0;
	for (size_t yy = 0 ; yy < tmp ; yy++){
	for (size_t xx = 0 ; xx < tmp ; xx++){	
	size_t index = ipx( intX[xx],intY[yy],tmp );
	val_data[count] = data[index];
	id_data[count] = id[index];
	count++;
	}
	}
	
		
	string result = _name + to_string(num);	
	compute_cumulants(result, (ilambda-1)*gridStep, val_data,id_data, pow(ilambda,3), IDP);

	vol_count++;
	}
	}
	
	
	delete [] LengthX;
	delete [] LengthY;

	delete [] intX;
	delete [] intY;

	delete [] val_data;
	delete [] id_data;
	}					
					
							
	if (iL0 % tmp > 0.0){
	int L1 = iL0 - (iL0 % tmp);
	int nvol = L1/ ((ilambda-1)*gridStep) + 1;

	int* LengthX;
	int* LengthY;

	int* intX;
	int* intY;

	int tmp = lambda;
	LengthX = new int[nvol];
	LengthY = new int[nvol];

	intX = new int[tmp];
	intY = new int[tmp];

	int lbx, ubx, lby, uby;	
		
	double* val_data;
	val_data = new double[ilambda * ilambda];

	int* id_data;
	id_data = new int[ilambda * ilambda];		
		
	for (size_t i = 0 ; i < nvol ; i++){
	//for overlapping walls
	LengthX [i] = i * (lambda - 1);
	LengthY [i] = i * (lambda - 1);
	
	//for non-overlapping walls
	//LengthX [i] = i * (lambda );
	//LengthY [i] = i * (lambda );

	}
			
	int vol_count = 0;
	for (size_t j = 0 ; j < nvol-1 ; j++){
	lby =  LengthY[j];
	uby =  LengthY[j + 1];
	for (size_t k = 0 ; k < nvol-1 ; k++){
	lbx =  LengthX[k];
	ubx =  LengthX[k + 1];

	for (size_t ii = 0 ; ii < tmp ; ii++){
	intX[ii] = lbx + ii;	
	intY[ii] = lby + ii;	
	}
	
	size_t count = 0;
	for (size_t yy = 0 ; yy < tmp ; yy++){
	for (size_t xx = 0 ; xx < tmp ; xx++){	
	size_t index = ipx(intX[xx],intY[yy],tmp);
	val_data[count] = data[index];
	id_data[count] = id[index];
	count++;
	}
	}
	

	string result = _name + to_string(num);	
	compute_cumulants(result, (ilambda-1)*gridStep, val_data, id_data, pow(ilambda,3), IDP);
	vol_count++;
	}
	}
	

	for (size_t ii = 0 ; ii <= 	(iL0 % (ilambda - 1)) ; ii++){
	intX[ii] = intX[tmp-1] + ii;
	intY[ii] = intY[tmp-1] + ii;	
	}
		
	for (size_t ii = 	(iL0 % (ilambda - 1))+1 ; ii <= (ilambda - 	(iL0 % (ilambda - 1))+1) ; ii++){
	intX[ii] = ii - 	(iL0 % (ilambda - 1))-1;
	intY[ii] = ii - 	(iL0 % (ilambda - 1))-1;
	}

	size_t count = 0;
	for (size_t yy = 0 ; yy < tmp ; yy++){
	for (size_t xx = 0 ; xx < tmp ; xx++){	
	size_t index = ipx(intX[xx],intY[yy],tmp );
	val_data[count] = data[index];
	id_data[count] = id[index];

	count++;
	}
	}
	

	string result = _name + to_string(num);	
	compute_cumulants(result, (ilambda-1)*gridStep, val_data, id_data, pow(ilambda,3), IDP);
	vol_count++;

	delete [] LengthX;
	delete [] LengthY;

	delete [] intX;
	delete [] intY;

	delete [] val_data;
	delete [] id_data;
		
	}
}

void compute_cumulants_th(const string& _name,int lamb, double tmp[], double trho[], int tid[] , int tot, int IDP, double rth) 
{
	
	string result = _name;
	string resultOUT = result;
	result.append("_");
	string resultIN = result + to_string(lamb);
	
	const char * cA = resultIN.c_str();

	double m1, m2, count, ExprA, ExprB, ExprC;
	m1 = m2 = count = ExprA = ExprB = ExprC = 0.0;

	vector<double> tmp_storage;
	
	for (int i = 0 ; i < tot ; i++){
	if (tid[i] == IDP && trho[i] > rth){
	count++;
	m1 += tmp[i];
	m2 += pow(tmp[i],2);
	tmp_storage.push_back(tmp[i]);
	}
	}
	
	m1 = m1/(count);
	m2 = m2/(count);
	
	for (int i = 0 ; i < tmp_storage.size() ; i++){
	ExprA += pow((tmp_storage[i] - m1),2);
	ExprB += pow((tmp_storage[i] - m1),3);
	ExprC += pow((tmp_storage[i] - m1),4);
	}

	double Mean = m1;//first cumulant
	double Std = sqrt(ExprA/(count - 1));
	//double Var = m2 - pow(m1,2);//second cumulant
	double Var = pow(Std,2);
	double Skew = (ExprB/count) / pow(sqrt(ExprA/count),3);
	double Kurt = (ExprC/count) / pow(ExprA/count,2);
	
	FILE * sortie;
	sortie = fopen(cA,"a");
	fprintf(sortie,"%g %g %g %g\n",Mean, Var, Skew, Kurt);
	fclose(sortie);	
}

/*
void compute_gr(const string& _name, int num, double lambda, double L0, double data[], int id[], int IDP, double thresh, double rcut, int nbins)	{

	
	int ilambda = lambda;
	int iL0 = L0;
	int tmp = (ilambda-1) * gridStep;
	int vol_count = 0;

	if (iL0 % tmp == 0.0){
	int nvol = L0 / ((lambda-1)*gridStep) + 1;
	int* LengthX;
	int* LengthY;
	int* LengthZ;
	int* intX;
	int* intY;
	int* intZ;
	int tmp = lambda;
	LengthX = new int[nvol];
	LengthY = new int[nvol];
	LengthZ = new int[nvol];
	intX = new int[tmp];
	intY = new int[tmp];
	intZ = new int[tmp];
	int lbx, ubx, lby, uby, lbz, ubz;

	double* bins;
	double* counts_ll;
	double* counts_gg;
	double* counts_gl;
	double* counts_lg;
	bins = new double[nbins+1];
	counts_ll = new double[nbins+1];
	counts_gg = new double[nbins+1];
	counts_gl = new double[nbins+1];
	counts_lg = new double[nbins+1];

	for (int i = 0 ; i < (nbins + 1); i++) {	
    	bins[i] = 1.*i * (rcut/(1.*nbins));
	counts_ll[i] = 0.0;
	counts_gg[i] = 0.0;
	counts_gl[i] = 0.0;
	counts_lg[i] = 0.0;
	}
		
	double* val_data;
	double* x_data;
	double* y_data;
	double* z_data;
	val_data = new double[ilambda * ilambda * ilambda];
	x_data = new double[ilambda * ilambda * ilambda];
	y_data = new double[ilambda * ilambda * ilambda];
	z_data = new double[ilambda * ilambda * ilambda];

	int* id_data;
	id_data = new int[ilambda * ilambda * ilambda];		
		
	for (size_t i = 0 ; i < nvol ; i++){
	LengthX [i] = i * (lambda - 1);
	LengthY [i] = i * (lambda - 1);
	LengthZ [i] = i * (lambda - 1);
	}
		
	for (size_t i = 0 ; i < nvol-1 ; i++){
	lbz =  LengthZ[i];
	ubz =  LengthZ[i + 1];
	for (size_t j = 0 ; j < nvol-1 ; j++){
	lby =  LengthY[j];
	uby =  LengthY[j + 1];
	for (size_t k = 0 ; k < nvol-1 ; k++){
	lbx =  LengthX[k];
	ubx =  LengthX[k + 1];

	for (size_t ii = 0 ; ii < tmp ; ii++){
	intX[ii] = lbx + ii;	
	intY[ii] = lby + ii;	
	intZ[ii] = lbz + ii;
	}

	size_t count = 0;
	int ll = 0;
	int gg = 0;
	int pp = 0;
	for (size_t zz = 0 ; zz < tmp ; zz++){
	for (size_t yy = 0 ; yy < tmp ; yy++){
	for (size_t xx = 0 ; xx < tmp ; xx++){	
	size_t index = ipx(intX[xx],intY[yy],intZ[zz]);
	val_data[count] = data[index];
	id_data[count] = id[index];
	x_data[count] = intX[xx] * gridStep;
	y_data[count] = intY[yy] * gridStep;
	z_data[count] = intZ[zz] * gridStep;
	if (data[index] >= thresh && id[index] == IDP) ll++;
	if (data[index] < thresh && id[index] == IDP) gg++;
	if (id[index] == IDP) pp++;
	count++;
	}
	}
	}

	//double pvol = (pp * pow(gridStep,3)) / pow( (ilambda-1)*gridStep , 3);
	//double pvol = (pp * pow(gridStep,3));
		
	string result = _name + to_string(num);	
	result.append("_");
	result = result + to_string(vol_count);
	gr(result, (ilambda-1)*gridStep, val_data, id_data, x_data, y_data, z_data, bins, counts_ll, counts_gg, counts_gl,counts_lg, pow(ilambda,3), IDP, thresh,nbins,rcut,vol_count,pp,ll,gg);
	cout<<"lamb: "<<(ilambda-1)*gridStep<<" vol: "<<vol_count<<endl;
	vol_count++;

	}
	}
	}	
	
	delete [] LengthX;
	delete [] LengthY;
	delete [] LengthZ;
	delete [] intX;
	delete [] intY;
	delete [] intZ;
	delete [] val_data;
	delete [] id_data;
	delete [] x_data;
	delete [] y_data;
	delete [] z_data;
	delete [] bins;
	delete [] counts_ll;
	delete [] counts_gg;
	delete [] counts_gl;
	delete [] counts_lg;
	}					
					
							
	if (iL0 % tmp > 0.0){
	int L1 = iL0 - (iL0 % tmp);
	int nvol = L1/ ((ilambda-1)*gridStep) + 1;

	int* LengthX;
	int* LengthY;
	int* LengthZ;
	int* intX;
	int* intY;
	int* intZ;
	int tmp = lambda;
	LengthX = new int[nvol];
	LengthY = new int[nvol];
	LengthZ = new int[nvol];
	intX = new int[tmp];
	intY = new int[tmp];
	intZ = new int[tmp];
	int lbx, ubx, lby, uby, lbz, ubz;	

	double* bins;
	double* counts_ll;
	double* counts_gg;
	double* counts_gl;
	double* counts_lg;
	bins = new double[nbins+1];
	counts_ll = new double[nbins+1];
	counts_gg = new double[nbins+1];
	counts_gl = new double[nbins+1];
	counts_lg = new double[nbins+1];

	for (int i = 0 ; i < (nbins + 1); i++) {	
    	bins[i] = 1.*i * (rcut/(1.*nbins));
	counts_ll[i] = 0.0;
	counts_gg[i] = 0.0;
	counts_gl[i] = 0.0;
	counts_lg[i] = 0.0;
	}
		
	double* val_data;
	double* x_data;
	double* y_data;
	double* z_data;
	val_data = new double[ilambda * ilambda * ilambda];
	x_data = new double[ilambda * ilambda * ilambda];
	y_data = new double[ilambda * ilambda * ilambda];
	z_data = new double[ilambda * ilambda * ilambda];

	int* id_data;
	id_data = new int[ilambda * ilambda * ilambda];		
		
	for (size_t i = 0 ; i < nvol ; i++){
	//for overlapping walls
	LengthX [i] = i * (lambda - 1);
	LengthY [i] = i * (lambda - 1);
	LengthZ [i] = i * (lambda - 1);
	}
			
	for (size_t i = 0 ; i < nvol-1 ; i++){
	lbz =  LengthZ[i];
	ubz =  LengthZ[i + 1];
	for (size_t j = 0 ; j < nvol-1 ; j++){
	lby =  LengthY[j];
	uby =  LengthY[j + 1];
	for (size_t k = 0 ; k < nvol-1 ; k++){
	lbx =  LengthX[k];
	ubx =  LengthX[k + 1];

	for (size_t ii = 0 ; ii < tmp ; ii++){
	intX[ii] = lbx + ii;	
	intY[ii] = lby + ii;	
	intZ[ii] = lbz + ii;
	}
	
	size_t count = 0;
	int ll = 0;
	int gg = 0;
	int pp = 0;
	for (size_t zz = 0 ; zz < tmp ; zz++){
	for (size_t yy = 0 ; yy < tmp ; yy++){
	for (size_t xx = 0 ; xx < tmp ; xx++){	
	size_t index = ipx(intX[xx],intY[yy],intZ[zz]);
	val_data[count] = data[index];
	id_data[count] = id[index];
	x_data[count] = intX[xx] * gridStep;
	y_data[count] = intY[yy] * gridStep;
	z_data[count] = intZ[zz] * gridStep;
	if (data[index] >= thresh && id[index] == IDP) ll++;
	if (data[index] < thresh && id[index] == IDP) gg++;
	if (id[index] == IDP) pp++;
	count++;
	}
	}
	}

	//double pvol = (pp * pow(gridStep,3)) / pow( (ilambda-1)*gridStep , 3);
	//double pvol = (pp * pow(gridStep,3));
		
	string result = _name + to_string(num);	
	result.append("_");
	result = result + to_string(vol_count);
	gr(result, (ilambda-1)*gridStep, val_data, id_data, x_data, y_data, z_data, bins, counts_ll, counts_gg, counts_gl,counts_lg, pow(ilambda,3), IDP, thresh,nbins,rcut,vol_count,pp,ll,gg);
	cout<<"lamb: "<<(ilambda-1)*gridStep<<" vol: "<<vol_count<<endl;
	vol_count++;

	}
	}
	}

	for (size_t ii = 0 ; ii <= 	(iL0 % (ilambda - 1)) ; ii++){
	intX[ii] = intX[tmp-1] + ii;
	intY[ii] = intY[tmp-1] + ii;	
	intZ[ii] = intZ[tmp-1] + ii;
	}
		
	for (size_t ii = 	(iL0 % (ilambda - 1))+1 ; ii <= (ilambda - 	(iL0 % (ilambda - 1))+1) ; ii++){
	intX[ii] = ii - 	(iL0 % (ilambda - 1))-1;
	intY[ii] = ii - 	(iL0 % (ilambda - 1))-1;
	intZ[ii] = ii - 	(iL0 % (ilambda - 1))-1;
	}

	size_t count = 0;
	int ll = 0;
	int gg = 0;
	int pp = 0;
	for (size_t zz = 0 ; zz < tmp ; zz++){
	for (size_t yy = 0 ; yy < tmp ; yy++){
	for (size_t xx = 0 ; xx < tmp ; xx++){	
	size_t index = ipx(intX[xx],intY[yy],intZ[zz]);
	val_data[count] = data[index];
	id_data[count] = id[index];
	x_data[count] = intX[xx] * gridStep;
	y_data[count] = intY[yy] * gridStep;
	z_data[count] = intZ[zz] * gridStep;
	if (data[index] >= thresh && id[index] == IDP) ll++;
	if (data[index] < thresh && id[index] == IDP) gg++;
	if (id[index] == IDP) pp++;
	count++;
	}
	}
	}

	//double pvol = (pp * pow(gridStep,3)) / pow( (ilambda-1)*gridStep , 3);
	//double pvol = (pp * pow(gridStep,3));
		
	string result = _name + to_string(num);	
	result.append("_");
	result = result + to_string(vol_count);
	gr(result, (ilambda-1)*gridStep, val_data, id_data, x_data, y_data, z_data, bins, counts_ll, counts_gg, counts_gl,counts_lg, pow(ilambda,3), IDP, thresh,nbins,rcut,vol_count,pp,ll,gg);
	cout<<"lamb: "<<(ilambda-1)*gridStep<<" vol: "<<vol_count<<endl;
	vol_count++;

	delete [] LengthX;
	delete [] LengthY;
	delete [] LengthZ;
	delete [] intX;
	delete [] intY;
	delete [] intZ;
	delete [] val_data;
	delete [] id_data;
	delete [] x_data;
	delete [] y_data;
	delete [] z_data;
	delete [] bins;
	delete [] counts_ll;
	delete [] counts_gg;
	delete [] counts_gl;
	delete [] counts_lg;
		
	}
}
*/

int readSimData(const string& _name, double rh, int TOT){
	
	int a1;
	double a2,a3,a4,a5,a6;
	int index_out;
	ifstream file(_name);

	for (int k = 0 ; k < TOT ; ++k){
	file >> a2 >> a3 >> a4 >> a5 >> a6 >> a1;
	if(a5 == rh) index_out = k;
	}

	return index_out;
}

void  voronoi_volume_frac(const string& _name, int npart, double rad, double fi[] ){
	int ipart;
	double xc, yc, cvol;
	fstream file(_name);
	for (int i = 0 ; i < npart ; i++){
	file >> ipart >> xc >> yc >> cvol;	
	fi[i] = ( 4.0 * M_PI * pow(rad,3)  / 3.0 ) / cvol;
	}
}


	
void  Import_data_cell( const string& _name, int cpart, int idpart[], double gs, double dist[], double xpar[], double ypar[] ){
	int ipart, index;
	double d, xc, yc;
	fstream file(_name);
	for (int i = 0 ; i < cpart ; i++){
	file >> ipart >> index >> d >> xc >> yc ;	
	dist[i] = d;
	idpart[i] = ipart;
	xpar[ipart - 1] = xc;
	ypar[ipart - 1] = yc;

	}
}

void  isolate_data(double tdata[], int tdata_size , double isodata[], int tid[], int id_flag){

	int j = 0;
	for (int i = 0 ; i < tdata_size ; i++){
	if(tid[i] == id_flag){
	isodata[j] = tdata[i];
	j++;
	}	
	}
}


void compute_statistics(const string& _name, vector<double> &vect){
	
	const char * cA = _name.c_str();
	double m1, m2, count, ExprA, ExprB, ExprC;

	m1 = m2 = count = ExprA = ExprB = ExprC = 0.0;
	
	for (int i = 0 ; i < vect.size() ; i++){
	count++;
	m1 += vect[i];
	m2 += pow(vect[i],2);
	}
		
	m1 = m1/(count);
	m2 = m2/(count);
	
	for (int i = 0 ; i < vect.size() ; i++){
	ExprA += pow((vect[i] - m1),2);
	ExprB += pow((vect[i] - m1),3);
	ExprC += pow((vect[i] - m1),4);
	}

	double Mean = m1;//first cumulant
	double Std = sqrt(ExprA/(count - 1));
	//double Var = m2 - pow(m1,2);//second cumulant
	double Var = pow(Std,2);
	double Skew = (ExprB/count) / pow(sqrt(ExprA/count),3);
	double Kurt = (ExprC/count) / pow(ExprA/count,2);
	
	FILE * sortie;
	sortie = fopen(cA,"a");
	if (vect.size() > 0.0)  fprintf(sortie,"%g %g %g %g\n",Mean, Var, Skew, Kurt);
	if (vect.size() == 0.0) fprintf(sortie,"0 0 0 0\n");
	fclose(sortie);
}



void create_xsection(const string& _name, int num, double data[], int tid[], int idp, int ids, int Ny, int Nx, double flag){
	
   string result = _name + to_string(num);
   string resultOUT = result;		
   const char * c = resultOUT.c_str();
   FILE * sortie;
   sortie = fopen(c, "w+");	

   int count = 0;

	for (int y = 0 ; y < Ny ; y++){ 
	 for (int x = 0 ; x < Nx ; x++){ 
		 if ( tid[count] == idp ) fprintf(sortie,"%g\n",data[count]);
		 if ( tid[count] == ids ) fprintf(sortie,"%g\n",flag);
	 count++;
	 }
	}

 fclose(sortie);
	
}

void  findGLsites(const string& _name, int tid[], double tdata[], int N, int IDP, double thresh){
	
   int ll = 0;
   int gg = 0;
   int pp = 0;
   for ( int i = 0 ; i < N ; i++){
   if (tid[i] == IDP && tdata[i] >= thresh) ll++;
   if (tid[i] == IDP && tdata[i] < thresh) gg++;
   if (tid[i] == IDP) pp++;
   }
	
   const char * c = _name.c_str();
   FILE * sortie;
   sortie = fopen(c, "a");	
   fprintf(sortie,"%i %i %i\n",pp,ll,gg);
   fclose(sortie);
}

void readSimParam(){
	
    const char * name = "../input/input_12152020.dat";
    ifstream file(name);
    file >> Lnx >> Lny;
	file >> nxp >> nyp;
	file >> Tstar >> gridStep >> Tcw >> yvar;
	file >> obj_tol >> npart;
	file >> rho_thresh_ub >> rho_thresh_lb;
	file >> Lx_lt_shift >> Lx_rt_shift;
	file >> Ly_lt_shift >> Ly_rt_shift;
}


void  Import_sres(const string& _name, double tdata[], int N){
	
	double msxx,msyy,msxy;

	fstream file(_name);
	
	for (int i = 0 ; i < N ; i++){					
	file >> msxx >> msyy >> msxy;
	tdata[i] = (1.0/2.0)*(msxx+msyy);
 	}
}





int main ()
{
	
	readSimParam();
	
	int IDP = 0;
	int IDS = 1;
	int nsim = 100;

	nx = Lnx / gridStep + 1;
	ny = Lny / gridStep + 1;

	double L0 = (nx - 1) * gridStep;
	
	int* id;
	int* id_mask;
	double* press;
	double* sat;

	double* part_x;
	double* part_y;
	int* pindex;

   	id 	= new int[ nx * ny ];
   	id_mask 	= new int[ nx * ny ];
   	press 	= new double[ nx * ny ];
   	sat 	= new double[ nx * ny ];

   	int NoPo_r = Import_id("../output/id_ordered", id, IDP, 0, nx, 0, ny );
	for (int i = 0 ; i < ( nx*ny ) ; i++) id_mask[i] = id[i];
	int NoPo_sub = Mask_id(id , id_mask , ( nx*ny ) , IDP , IDS);
	int NoPo_p = NoPo_r - NoPo_sub;
	pindex = new int[NoPo_p];

	double* press_pore;
	double* sat_pore;
	press_pore = new double[NoPo_p];
	sat_pore = new double[NoPo_r];

	for (int i = 0 ; i < nsim ; i++){
		
	cout<<"importing data :"<<i<<endl;
	Import_press("../output/p_", i, press,0,nx,0,ny,0 );
	compute_cumulants("../output/p_tot_stat",99, press, id , (nx*ny), IDP);
	compute_cumulants("../output/p_pore_tot_stat",99, press, id_mask , (nx*ny), IDP);
	
	compute_cumulants_th("../output/p_pore_tot_stat_th",1, press,sat, id_mask , (nx*ny), IDP, 0.50);
	compute_cumulants_th("../output/p_pore_tot_stat_th",2, press,sat, id_mask , (nx*ny), IDP, 0.55);
	compute_cumulants_th("../output/p_pore_tot_stat_th",3, press,sat, id_mask , (nx*ny), IDP, 0.60);
	compute_cumulants_th("../output/p_pore_tot_stat_th",4, press,sat, id_mask , (nx*ny), IDP, 0.65);
	compute_cumulants_th("../output/p_pore_tot_stat_th",5, press,sat, id_mask , (nx*ny), IDP, 0.70);
	compute_cumulants_th("../output/p_pore_tot_stat_th",6, press,sat, id_mask , (nx*ny), IDP, 0.75);
	compute_cumulants_th("../output/p_pore_tot_stat_th",7, press,sat, id_mask , (nx*ny), IDP, 0.80);
	compute_cumulants_th("../output/p_pore_tot_stat_th",8, press,sat, id_mask , (nx*ny), IDP, 0.85);
	compute_cumulants_th("../output/p_pore_tot_stat_th",9, press,sat, id_mask , (nx*ny), IDP, 0.90);	
		
	isolate_data(press,(nx*ny),press_pore,id_mask,IDP);
	double min_data = find_min("../output/p_min", press_pore, NoPo_p);
	double max_data = find_max("../output/p_max", press_pore, NoPo_p);

	Import_data("../output/rho_",i, sat,0,nx,0,ny);
	compute_cumulants("../output/r_tot_stat",99, sat, id , (nx*ny), IDP);
	
	isolate_data(sat,(nx*ny),sat_pore,id,IDP);
	min_data = find_min("../output/r_min", sat_pore, NoPo_r);
	max_data = find_max("../output/r_max", sat_pore, NoPo_r);

	}

	delete [] press_pore;
	delete [] sat_pore;	
	delete [] id;
	delete [] id_mask;
	delete [] sat;
	delete [] press;
	delete [] part_x;
	delete [] part_y;
	delete [] pindex;

	cout<<"computation done."<<endl;

	return 0;
}


	

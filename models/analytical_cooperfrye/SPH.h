#ifndef _SPH_H_
#define _SPH_H_

#include <string>
#include <iostream>
#include <vector>
#include "vector.h"
#include "vnvar.h"
#include "spectra.h"
#include "bessel.h"
#include "tables.h"

	

template <int D,int DD>
class SPH {
private:
	
	
	double facc;
	static const int Nmax=10;
	
	double I1c,I2c,I1sc,I2sc;
	int h1tot;
	
	
	double Eperp(double p, double m) {return sqrt(p*p+m*m);}
	double pperp(double p, double phi, Vector<double,DD> u) {return p*(u.x[1]*cos(phi)+u.x[2]*sin(phi));}
	//double nuprod( int nsph) {return par[nsph].n&par[nsph].u;   }
	
	
	 static const string bulk;
	 static const string ideal ;	
	 static const string shear ;
	 static const string shearbulk;	
	 static const string bulkshear;
     static const string idealBSQ;
	 static const string v2;
	
	
	
public:
	static const double scale=0.1973;
	static const double sc3=0.1973*0.1973*0.1973;
	double T,muB,muS,muQ; //I added muBSQ here
	double m3;
	int N,NHAD;
	int spectra;
	int vnmax;
	int ngl;
	double *w,*z;
	int ntab[6];
	double ** vout, ** voutc;
	 
	
	struct PAR
	{
	Vector<double,DD> u,n,r;
	double vol,bulk,bulkpi,nu,tau,s,T,energy;
	double pi00,pi11,pi22,pi33,pi12;
    double muB,muS,muQ;// added the muBSQ info for each sph particle as well
	double enthalpy,cs2;
	char eosname;
	double spec;
	}; // structure that contains the basic info for each SPH particle e.g. normal vecotrs, four-velocity etc
	
	struct HAD
	{
	double mass,deg,theta,spar;
	double form3,halm2;
	double E0,B0,D0;
	double vfac,svfac;
	int id,baryon,strange,electric;//here these refer to the BSQ-charge qm#s of each hadron
	int var8,var9,var10,var11,var12;
	double C00,mT;
	  int anti,null;
	  vector<int> sids;
	}; // structure that is set up for each observed hadron. 

	struct HLIST
	{
	double E0,B0,D0;
	double C00,mT;
	int T;
	}; // structure that is set up for each observed hadron. 
	
	HLIST *hlist;

	int evn,evncor;  // number of sph particles, file number of event
	PAR *par;
	VNVAR *ide,*cor;
	Vector<double,DD> *qv;
	
	
	double outc;
	int typ;
	vector<HAD> had,nmes;
	string before,after,rnum;
	int start, end;
	string folder;
	string flist2;
	double *Ia,*Iac;
	
	string neg;
	
	double pt_s,pt_e,pt_step;
	double *v;
	double *vc;
	vector<double> pt,phi;
	
	
	SPH<D,DD>();
	~SPH<D,DD>();
	void readin(_inputIC ics);
	void readin2(int cev);
	double dNdpdphi(double p, double phi, HAD cur);
	void Iout(double &I1, double &I2, double p, double phi, HAD cur,int nsph);
	string convertInt(int number);
	void flist();
	void calcF2(HAD cur, int nsph,double pd,double &F0,double &F1, double &F2);
	void checknu( );
	void setcoef(int i);
	void calcsPI(int h);
	double totmul(int ht);
	void delavg(double avg, int &n, PAR * par2);
	void restart(int pTmax);
	void restartc(int pTmax);
	void setvns(int ptmax );
	void destroyvout(int pTmax) ;
	void destroyvoutc(int pTmax);
	void printrun(string ofolder);
	
	void readweights(string grid,string grid2);
	string grid,grid2;
};

        template <int D,int DD> const string SPH<D,DD>::bulk="bulk";
	template <int D,int DD> const string SPH<D,DD>::ideal="ideal";	
	template <int D,int DD> const string SPH<D,DD>::shear="shear";
	template <int D,int DD> const string SPH<D,DD>::shearbulk="shear+bulk";	
	template <int D,int DD> const string SPH<D,DD>::bulkshear="bulk+shear";
        template<int D,int DD> const string SPH<D,DD>::idealBSQ="idealBSQ";
		template<int D,int DD> const string SPH<D,DD>::v2="v2";

template <int D,int DD>
SPH<D,DD>::SPH()
{
	facc=1/(2*PI*PI);
	
}

template <int D,int DD>
SPH<D,DD>::~SPH()
{
	
}


// reads in the basic information from "input.dat" such at the number of events, events folder, particles to observe etc.
template <int D,int DD>
void SPH<D,DD>::readin(_inputIC ics)
{

	string infile="input/"+ics.man;
	string resolist,hadclist;
        FILE * myfile = fopen (infile.c_str(),"r");
        if(myfile== NULL)
        {
  		cout << "Error: "<< infile << " does not exist. \n";
		exit(1);	
  	}	
		
		
	   char charin[100];
           fscanf(myfile,"%*s %s \n",charin); // type of equations (ideal, bulk etc)
           string type=charin;
           fscanf(myfile,"%*s  %s \n",charin); // folder that contains the events
           folder=charin; 
           fscanf(myfile,"%*s  %s ",charin); // reads in pt points
           grid=charin;
           fscanf(myfile,"  %s ",charin);  //reads in phi points
           grid2=charin;
           fscanf(myfile,"%*s  %*s \n");
           fscanf(myfile,"%*s  %i %i \n",&start,&end); // range of events e.g. start=0 end=199
           if (ics.on==1){
  	   	start=ics.start; 
  	   	end=ics.end;
  	   }
           fscanf(myfile,"%*s %lf \n",&T); // here's our problem
           
                   
  	   fscanf(myfile,"%*s %s \n",charin); // file that contains df corrections
           flist2=charin;
           fscanf(myfile,"%*s  %*s \n");  // file that contains phi/pt grid
          
           fscanf(myfile,"%*s  %s",charin);  // file that contains list of hadrons and decays
           resolist=charin;
            fscanf(myfile,"%*s  %s",charin);  // file that contains list of hadrons to actually compute
           hadclist=charin;
           fclose(myfile);
           
           
          FILE * myfile3 = fopen (hadclist.c_str(),"r");
        if(myfile3== NULL)
        {
  		cout << "Error: input/numbers.dat does not exist. \n";
		exit(1);	
  	}	
  		
	 
           vector <double> hcheck;
           int subc;
           while(fscanf(myfile3,"%i", &subc)==1)
           {
           	  
           	hcheck.push_back(subc);
           }
           fclose(myfile3); 
           
           
        // reads in the resonances and their properties
	FILE * myfile2 = fopen (resolist.c_str(),"r");
        if(myfile2== NULL)
        {
  		cout << "Error: decay/input/"<< resolist<< " does not exist. \n";
		exit(1);	
  	}	
	 
           int decays;
           HAD had2;
           vector <HAD> had3;//here SQ of each hadron not just the baryon number? 6th=S; 11th=Q int
           while(fscanf(myfile2,"%i %*s %lf %*f %lf %d %d %*i%*i %*f %d %i ", &had2.id, &had2.mass,&had2.deg, &had2.baryon,&had2.strange,&had2.electric,&decays)==7)
           {
			   
           	if (had2.baryon==0) had2.theta=-1; // if it's a meson, -1 for fermi dirac distribution
           	else had2.theta=1; //if it's a baryon, +1 for fermi diract distribution
			
			//cout <<had2.id << " " << had2.mass << " " << had2.baryon << " "  << had2.strange << " " << had2.electric << " " << endl;
			//getchar();
			
           	int hc=0,hsize=hcheck.size();
           	while (hc<hsize){
				if (had2.id==hcheck[hc])
				{
	/*            		if (had2.baryon<0) had3.push_back(had2);
					else if (had2.id<0&&had2.baryon==0) nmes.push_back(had2);
				else {
				  had2.anti=0;
				  had2.null=0;
				  had.push_back(had2);//we need to come back to this */
				  had.push_back(had2);
				  break;
				}
				else hc++;
				}

           	for (int l=0;l<decays;l++)
           	{
           	fscanf(myfile2,"%*i%*i%*f%*i%*i%*i%*i%*i");
           	}
           	
           }
           cout << "size " <<  had.size() << endl;
           fclose(myfile2);
	cout << "Resonances read in!\n"; 
	
           NHAD=had.size();
        
	if (type==ideal){
		typ=0;
		before="freezeout_ev";}
	else if (type==bulk){
		typ=1;
		before="bvfreezeout_ev";}
	else if (type==shear){
		typ=2;
		before="svfreezeout_ev";}
	else if ((type==shearbulk)||(type==bulkshear)){
		typ=3;
		before="sbvfreezeout_ev";}
        else if (type==idealBSQ){
        typ=4;
        before="bsqsbvfreezeout_ev";} //added the freezeout files for BSQ: name needs to match ccake
    else if(type == v2){
		typ=5;
		before="freeze_out";

	}	
	after=".dat";
	N=end-start+1; // determines the total number of events       
        
        readweights(grid,grid2);
  	
	
}

// reads in the basic information from "input.dat" such at the number of events, events folder, particles to observe etc.
template <int D,int DD>
void SPH<D,DD>::readin2(int cev)
{
	string event="input/"+folder+"/"+before;
	event+=convertInt(cev)+after;
	FILE * myfile = fopen (event.c_str(),"r");
	if (myfile==NULL) 
	{
	cout << "Error: Can't open Event " << cev << endl;
	cout << "resonance number=" << event<< endl;
	exit(1);
	}
	
	PAR * par2=new PAR [200000];//why 200000 here?. Should not this depend on the event or is it a maximum?
	int i=0;
	int j=0;
	double avgnu=0;
	
	if (typ==0) //ideal
	{
//	fscanf(myfile,"%lf\n",&s); 
	while (!feof(myfile))   // runs over all the SPH particles in the event
	{
	fscanf(myfile,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",&par2[i].n.x[0],&par2[i].n.x[1],&par2[i].n.x[2]
	,&par2[i].u.x[0],&par2[i].u.x[1],&par2[i].u.x[2],&par2[i].vol,&par2[i].tau,&par2[i].r.x[0],&par2[i].r.x[1],&par2[i].s
	,&par2[i].T); 

	par2[i].vol/=sc3;		
	
	
	par2[i].nu=par2[i].n&par2[i].u;
	
	if (par2[i].nu<0) j++;
	else
	{
		avgnu+=par2[i].nu;
		i++;
	}

	}
	
	}
	else if (typ==1) // bulk
	{	

	while (!feof(myfile))   // runs over all the SPH particles in the event
	{
	fscanf(myfile,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf  %lf  %lf\n",&par2[i].n.x[0],&par2[i].n.x[1],&par2[i].n.x[2]
	,&par2[i].u.x[0],&par2[i].u.x[1],&par2[i].u.x[2],&par2[i].vol,&par2[i].bulkpi,&par2[i].tau,&par2[i].r.x[0],&par2[i].r.x[1]
	,&par2[i].s,&par2[i].T); 
	par2[i].vol/=sc3;
	
	par2[i].nu=par2[i].n&par2[i].u;
	if (par2[i].nu<0)
	{
		
		j++;
	}
	else
	{
		
		avgnu+=par2[i].nu;
		i++;
	}
	
	
	}
	}
	else if (typ==2) //shear
	{
	
	while (!feof(myfile))  //runs over all SPH particles in the event
        {
         fscanf(myfile,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %s %lf %lf\n"
		 ,&par2[i].n.x[0],&par2[i].n.x[1],&par2[i].n.x[2],&par2[i].u.x[0],&par2[i].u.x[1],&par2[i].u.x[2],&par2[i].vol
		 ,&par2[i].bulkpi,&par2[i].pi00,&par2[i].pi11,&par2[i].pi22,&par2[i].pi33,&par2[i].pi12,&par2[i].tau,&par2[i].r.x[0]
		 ,&par2[i].r.x[1],&par2[i].s,&par2[i].energy,&par2[i].T,&par2[i].muB,&par2[i].muS,&par2[i].muQ,&par2[i].eosname
		 ,&par2[i].enthalpy,&par2[i].cs2);
           par2[i].vol/=sc3;
           par2[i].pi33*=par2[i].tau*par2[i].tau;
           par2[i].nu=par2[i].n&par2[i].u;
	if (par2[i].nu<0)
	{
		
		j++;
	}
	else
	{
		
		avgnu+=par2[i].nu;
		i++;
	}
	
	
	}
	}
	else if (typ==3 ) //bulk+shear
	{
	
	while (!feof(myfile))   // runs over all the SPH particles in the event
	{
	fscanf(myfile,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %s %lf %lf\n"
		 ,&par2[i].n.x[0],&par2[i].n.x[1],&par2[i].n.x[2],&par2[i].u.x[0],&par2[i].u.x[1],&par2[i].u.x[2],&par2[i].vol
		 ,&par2[i].bulkpi,&par2[i].pi00,&par2[i].pi11,&par2[i].pi22,&par2[i].pi33,&par2[i].pi12,&par2[i].tau,&par2[i].r.x[0]
		 ,&par2[i].r.x[1],&par2[i].s,&par2[i].energy,&par2[i].T,&par2[i].muB,&par2[i].muS,&par2[i].muQ,&par2[i].eosname
		 ,&par2[i].enthalpy,&par2[i].cs2);; 
	par2[i].vol/=sc3;
	par2[i].pi33*=par2[i].tau*par2[i].tau;
	par2[i].nu=par2[i].n&par2[i].u;
	if (par2[i].nu<0)
	{
		
		j++;
	}
	else
	{
		
		avgnu+=par2[i].nu;
		i++;
	}
	
	
	}
	}
        else if (typ==4)//bulk+shear+idealBSQ
        {

        while (!feof(myfile))  //runs over all SPH particles in the event
        {
         fscanf(myfile,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %s %lf %lf\n"
		 ,&par2[i].n.x[0],&par2[i].n.x[1],&par2[i].n.x[2],&par2[i].u.x[0],&par2[i].u.x[1],&par2[i].u.x[2],&par2[i].vol
		 ,&par2[i].bulkpi,&par2[i].pi00,&par2[i].pi11,&par2[i].pi22,&par2[i].pi33,&par2[i].pi12,&par2[i].tau,&par2[i].r.x[0]
		 ,&par2[i].r.x[1],&par2[i].s,&par2[i].energy,&par2[i].T,&par2[i].muB,&par2[i].muS,&par2[i].muQ,&par2[i].eosname
		 ,&par2[i].enthalpy,&par2[i].cs2);
           par2[i].vol/=sc3;
           par2[i].pi33*=par2[i].tau*par2[i].tau;
           par2[i].nu=par2[i].n&par2[i].u;
		   //cout << par2[i].T << " " << par2[i].muB << " " << par2[i].muS << " " << par2[i].muQ << endl;
        if (par2[i].nu<0)
        {
          j++;
        }
        else
        {

                 avgnu+=par2[i].nu;
                 i++;
        }

        }
        } 
	else if (typ==5) //v2
	{
	while (!feof(myfile))   // runs over all the SPH particles in the event
	{
	fscanf(myfile,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n"
		 ,&par2[i].n.x[0],&par2[i].n.x[1],&par2[i].n.x[2],&par2[i].u.x[0],&par2[i].u.x[1],&par2[i].u.x[2],&par2[i].vol
		 ,&par2[i].bulkpi,&par2[i].pi00,&par2[i].pi11,&par2[i].pi22,&par2[i].pi33,&par2[i].pi12,&par2[i].tau,&par2[i].r.x[0]
		 ,&par2[i].r.x[1],&par2[i].s,&par2[i].energy,&par2[i].T,&par2[i].muB,&par2[i].muS,&par2[i].muQ
		 ,&par2[i].enthalpy,&par2[i].cs2);; 
	par2[i].vol/=sc3;
	par2[i].pi33*=par2[i].tau*par2[i].tau;
	par2[i].nu=par2[i].n&par2[i].u;
	if (par2[i].nu<0)
	{
		
		j++;
	}
	else
	{
		
		avgnu+=par2[i].nu;
		i++;
	}
	
	
	}
	}	
	fclose(myfile);
	
	avgnu/=i;
	delavg(avgnu,i,par2);
	delete [] par2;
	
	// not sure if it is needed?
	evn=i;
	evncor=cev;
	cout << evn << endl;
	cout << "% of negative SPH particles=" << (j*100.)/(i+j) << "%" << endl;
	cout << "Pions=" << totmul(0) <<  endl;
	//cout << "Pi0=" << totmul(0) <<  " Pi+=" << totmul(1) <<  " K+=" << totmul(2)<<  " P+=" << totmul(3) <<  endl;
	cout << "Event "<< cev << endl;
	//cout << &par2.T << " " << &par2.muB << " " << &par2.muS << " " << &par2.muQ << " " << endl;
	//cout << par[i].T << " " << par[i].muB << " " << par[i].muS << " " << par[i].muQ << endl;
    //cout << &par2[i].T << " " << &par2[i].muB << " " << &par2[i].muS << " " << &par2[i].muQ << endl;
	//getchar();
	
}

template <int D,int DD>
void SPH<D,DD>::delavg(double avg, int &n, PAR * par2)
{
	
	
	int *list=new int [n];
	
	int nl=0;
	double avgp=avg/100;
	for(int i=0;i<n;i++)
	{
		if (par2[i].nu<avgp)
		{
			list[i]=1;
			//cout << i << endl;
		
		}
		else if (par2[i].n.x[0]==1){
			list[i]=1;
		}
		else 
		{
			list[i]=0;
			nl++;
		}
	
	}
	
	par=new PAR [nl];
	
	int nc=0;
	for(int i=0;i<n;i++)
	{
		if (list[i]==0)
		{
			par[nc]=par2[i];
			nc++;
			//delete [] par;
		}
	}
	delete [] list;
	//delete [] par;
	//delete [] qv;
	
	cout << "SPH with too small nu= " << n-nl << endl;
	n=nl;
	
	
}

template <int D,int DD>
double SPH<D,DD>::totmul(int ht)
{

	double tot=0,tot2=0;
	double pre,inner;
	pre=had[ht].deg/(2*PI*PI)*had[ht].mass*had[ht].mass*T;
	
	//cout << had[ht].mass << " " << T << endl;
	//cout << "hmass= " <<  had[ht].mass << "T= " <<  T << "theta= " <<  had[ht].theta << endl;
	//getchar();
	
	inner=had[ht].mass/T;
	
	for(int j=0;j<10;j++)
	{
		double add=j+1;
		double sub=add*inner;
		Bessel bes;
		tot2+=pow(-had[ht].theta,j)/add*bes.Kn(2,sub);
	
	}
	tot2*=pre;
	//cout << "rho_pi=" << tot2 << endl;
	
	
	for(int i=0;i<evn;i++)
	{
		tot+=tot2*par[i].vol;
	
	}
	return tot;

}


template <int D,int DD>
void SPH<D,DD>::printrun(string ofolder){

	ofstream OUT;
	
	string vtyp;
	
	
	
	
	if (typ==0) vtyp="i"; 
	if (typ==1) vtyp="bvc"; 
	else if (typ==2) vtyp="svc"; 
	else if (typ==3) vtyp="sbvc"; 
        else if (typ==4) vtyp="bsqsbvc";	
	else if (typ==5) vtyp="v2";
	
	string pre=ofolder+"/ev";
	string post=vtyp+"_dNdphidpp.dat";
	
	
	
  	}




template <int D,int DD>
void SPH<D,DD>::flist()
{
//////////// CHECK hlist.T ????????????
	
	
	ifstream input(flist2.c_str());
	
	if (!input.is_open())
 	{
 	cout << "Can't open " << flist2 << endl;
 	exit(1);
 	}

	
	string line;
	
	h1tot=101;
	getline(input,line);
	
	hlist=new HLIST [h1tot];
	int i=0;
	while (input >> hlist[i].T >>  hlist[i].E0 >> hlist[i].D0 >> hlist[i].B0)   
	{ 	
	
	i++;
	if (i>h1tot) break;
	}
	input.close();
	
	h1tot=i;
	
	double Tnew=T*1000;
	for (int ru=0;ru<h1tot;ru++){
	
	if ((hlist[ru].T<(Tnew+0.01))&&(hlist[ru].T>(Tnew-0.01))) {
	h1tot=1;
	hlist[0].E0=hlist[ru].E0;
	hlist[0].D0=hlist[ru].D0;
	hlist[0].B0=hlist[ru].B0;
	cout <<  hlist[ru].T << " " << hlist[ru].B0 << endl;
	break;
	}
	
	}
}

template <int D,int DD>
void SPH<D,DD>::setcoef(int i)
{
	

	
  if (h1tot==1)
  {
	had[i].E0=hlist[0].E0;
  	had[i].D0=hlist[0].D0;
  	had[i].B0=hlist[0].B0;
  }
  else{
	double mT=had[i].mass/T;
	if (mT<hlist[0].mT) cout << "Error: mass is smaller than the pion!" << endl;

  	for (int j=1;j<h1tot;j++)
  	{
  	if (mT<hlist[j].mT) 
  	{
  		double mfac=(mT-hlist[j-1].mT)/(hlist[j].mT-hlist[j-1].mT);
  		had[i].E0=hlist[j-1].E0+mfac*(hlist[j].E0-hlist[j-1].E0);
  		had[i].D0=hlist[j-1].D0+mfac*(hlist[j].D0-hlist[j-1].D0);
  		had[i].B0=hlist[j-1].B0+mfac*(hlist[j].B0-hlist[j-1].B0);
  		cout << "constants" << endl;
  		cout << had[i].E0 << " " << had[i].D0 << " "  << had[i].B0 << endl;
  		break;
  	
  	}
  
 	}
  }
  
  
  
}

template <int D,int DD>
string SPH<D,DD>::convertInt(int number)
{
   stringstream ss;//create a stringstream
   ss << number;//add number to the stream
   return ss.str();//return a string with the contents of the stream
}





template <int D,int DD>
double SPH<D,DD>::dNdpdphi(double p, double phi, HAD cur)
{
	
	double vfac=cur.vfac;
	double out=0,outsc=0;
	outc=0;
	string negc="neg";
	
	
	
	for (int i=0;i<evn;i++)
	{
	
	double I1,I2;//added mu's
	
	Iout(I1,I2,p,phi,cur,i);//here is where we pass the mu's to be read in deltaf
	
	
	
	double qp=p*cos(phi)*qv[i].x[1]+p*sin(phi)*qv[i].x[2];
	
	double qtot=qv[i].x[0]+qp;
	
	double sub=qv[i].x[0]*I1+qp*I2;
	if (neg!=negc){ if ((sub<0)||qtot<0) sub=0;}
	out+=sub;  
	
	
	
        //Added or typ==4	 
	if ((typ==1)||(typ==3)||(typ==4) || (typ==5)) {double sub2=qv[i].x[0]*I1c+qp*I2c;
	  if (neg!=negc){if ((sub2<0)||qtot<0||isnan(sub2)) sub2=0;}
	outc+=sub2;}
	if (typ>1) {double sub3=(qv[i].x[0]*I1sc+qp*I2sc)/par[i].s; 
	  if (neg!=negc){ if ((sub3<0)||qtot<0||isnan(sub3)) sub3=0;}
	outsc+=sub3; }
	
	//cout << cur.id <<" sub " << sub << " outc=" << outc << " outsc=" << outsc << endl;
	}
	
	
	

	if (isnan(out)==1) cout << "negc =" <<out << endl;
	

	if (typ==1)  outc*=vfac;
	else if (typ==2) outc=vfac*out+cur.svfac*outsc;
	else if (typ==3) outc=vfac*outc+cur.svfac*outsc;
	else if (typ==4) outc=vfac*outc+cur.svfac*outsc;
	else if (typ==5) outc=vfac*out+cur.svfac*outsc; //added for v2
        //else if (typ==4) outc=like 3 	
	//cout << p << " " << phi << " " <<  out*vfac << " " << outc << " " <<  cur.svfac*outsc << endl;
	
	
	

	return out*=vfac;
	//cout << p << " " << phi << " " <<  out*vfac << " " << outc << " " <<  cur.svfac*outsc << endl;

}

template <int D,int DD>
void SPH<D,DD>::checknu( )
{
	qv=new Vector<double,DD> [evn];
	for (int i=0;i< evn;i++)
	{
	qv[i]=(par[i].vol/par[i].nu)*par[i].n;
	}
}



template <int D,int DD>
void SPH<D,DD>::Iout(double &I1, double &I2, double p, double phi, HAD cur,int nsph)// Added to read mu's as well
{
	double out1=0,out2=0;
	double out1c=0,out2c=0;
	double pd=pperp(p,phi,par[nsph].u);
	double b0,b1,b2,bsub,fac,pre;
	double g=par[nsph].u.x[0];
	double eperp=Eperp(p,cur.mass);
	double f0s,f1s,f2s;
	double F0c,F1c,F2c;
	double px=p*cos(phi);
	double py=p*sin(phi);
	double px2=px*px,py2=py*py,pxy=2*px*py;
	double ep2=eperp*eperp;
	double ep3=ep2*eperp/4.;
	double TG=par[nsph].T/g;
	double bfac=eperp/TG;
	
	
	
      //I added type==4 here
	if ((typ==1)||(typ==3)||(typ==4)||(typ==5)) calcF2(cur,nsph,pd,f0s,f1s,f2s);
	if (typ>1) {I1sc=0;I2sc=0;}
	
	double expT;
	
	 if ((abs(cur.baryon*par[nsph].muB*0.1973)+ abs(cur.strange*par[nsph].muS*0.1973)+ abs(cur.electric*par[nsph].muQ*0.1973)) > cur.mass)
	{expT=0;}  
	
	 /* if (((cur.baryon*par[nsph].muB*0.1973)+ (cur.strange*par[nsph].muS*0.1973)+ (cur.electric*par[nsph].muQ*0.1973)) > cur.mass)
	{expT=0;} */

	 else {expT=exp((pd/par[nsph].T)-((cur.baryon*par[nsph].muB*0.1973 
	              + cur.strange*par[nsph].muS*0.1973+ cur.electric*par[nsph].muQ*0.1973)/par[nsph].T));}  
	
	
	/* expT=exp(pd/par[nsph].T); */


	
	//cout << " T= " << (par[nsph].T) <<  " muB/T= " << par[nsph].muB/par[nsph].T << " muS/T= " << par[nsph].muS/par[nsph].T << " muQ/T= " << par[nsph].muQ/par[nsph].T << endl;
	//cout << " T= " << (par[nsph].T) <<  " muB/T= " << par[nsph].muB*0.197/par[nsph].T << " muS/T= " << par[nsph].muS*0.197/par[nsph].T << " muQ/T= " << par[nsph].muQ*0.197/par[nsph].T << endl;
	/* cout << " fwrong =  " << exp((pd/par[nsph].T)-((cur.baryon*par[nsph].muB 
	              + cur.strange*par[nsph].muS+ cur.electric*par[nsph].muQ)/par[nsph].T)) <<  " condition = " << (((cur.baryon*par[nsph].muB*0.1973 
	              + cur.strange*par[nsph].muS*0.1973+ cur.electric*par[nsph].muQ*0.1973))) << endl; */
	//getchar();
	//cout << "expT = " << (expT) << endl;
	for (int nn=0;nn<=Nmax;nn++)
	{
		
		Bessel bes;
		double add=(nn+1);
		bsub=add*bfac;
		b0=bes.K0(bsub);
		b1=bes.K1(bsub);
		
		
		
		pre=pow(-cur.theta,nn)*pow(expT,add);
		double preb1=pre*b1;
		out1+=preb1;
		double preb0=pre*b0;
		out2+=preb0;
		
		//cout << "out1 = " << out1 << endl;
		
		
		
	//add a routine down for type 4	
		if ((typ==1)||(typ==3)||(typ==4)||(typ==5))
		{
		
		
		fac=TG/add;
		b2=bes.Kn(2,bsub);
		F0c=1+add*f0s;
		F1c=add*f1s;
		F2c=add*f2s;
		
		double prep=preb1*eperp;
		double facF2=fac*F2c;
		double F0F2=F0c+F2c*ep2;
		
		
		out1c+=ep2*(F1c*preb0+pre*facF2*b2)+prep*(F0F2+fac*F1c);
        //cout << "out1c = " << out1c << endl;		
		out2c+=preb0*F0F2+prep*(F1c+facF2);
		}
		
		if (typ==2 ||typ==3 ||(typ==4)||(typ==5)) 
		b2=bes.Kn(2,bsub);
		double pred=pre*add;
		
		double spi1=par[nsph].pi00+par[nsph].pi33;
		double spi3=px2*par[nsph].pi11 +py2*par[nsph].pi22 +pxy* par[nsph].pi12;
		
		
		I1sc+=pred*( ep3*spi1*bes.Kn(3,bsub)+(ep3*(3*par[nsph].pi00- par[nsph].pi33)+eperp*spi3   )*b1   );
		I2sc+=pred*( 0.5*ep2*spi1*b2+(0.5*ep2*(par[nsph].pi00- par[nsph].pi33)+spi3)*b0);
		//cout << "I2sc = " << I2sc << endl;
		
	}
	


	I1=2*out1*eperp;
	I2=2*out2;
	//cout << "I1 = " << I1  << "out1 = " << out1 << "eperp = " << eperp << endl;
	//cout << "I2 = " << I2  << "out2 = " << out2 <<  endl;
	//cout << "I1sc = " << I1sc  << "I2sc = " << I2sc <<  endl;
	
	I1c=2*out1c;
	I2c=2*out2c;
	//cout << "I1c = " << I1c << "out1c = " << out1c<< endl;
	//cout << "I2c = " << I1c << "out2c = " << out2c<< endl;
	//cout << " outc = " << outc << endl;
	
	//getchar();
//	if ((pd/par[nsph].T)>64) {I1=0;
//	I2=0;
//	I1c=0;
//	I2c=0;
//	
//	}
	
	

	
}

template <int D,int DD>
void SPH<D,DD>::calcsPI(int h)
{
   
   m3=pow(had[h].mass,3);
   double bfac=had[h].mass/T;
   double bsum=0;
   for(int i=1;i<10;i++)
   {
   double bsub2=i*bfac;
   Bessel bes;
   bsum+=bes.Kn(3,bsub2)/i;
   }
   
   
   had[h].spar=had[h].deg*m3*facc*bsum/sc3;
 
   
}

template <int D,int DD>
void SPH<D,DD>::calcF2(HAD cur, int nsph, double pd,double &F0,double &F1, double &F2)
{
	F0=par[nsph].bulkpi*(cur.E0 - cur.D0*pd+cur.B0*pd*pd);
	F1=par[nsph].bulkpi*par[nsph].u.x[0]*(cur.D0-2*cur.B0*pd);
	F2=par[nsph].bulkpi*par[nsph].u.x[0]*par[nsph].u.x[0]*cur.B0;
}

template <int D,int DD>
void SPH<D,DD>::readweights(string ptpoints, string phipoints)
{

	//cout << "start" << endl;

	FILE * myfile2 = fopen (phipoints.c_str(),"r");
	
	cout << phipoints << endl;
	double sub;
	while (fscanf(myfile2,"%lf %*f",&sub)==1)
	{
	phi.push_back(sub);
	
	}
	
	fclose(myfile2);



	FILE * myfile = fopen (ptpoints.c_str(),"r");
	 
	
	while (fscanf(myfile,"%lf %*f",&sub)==1)
	{
	pt.push_back(sub); 
	
	}
	
	fclose(myfile);

	

}



#endif

#ifndef _INT_H_
#define _INT_H_


#include "SPH.h"
#include "spectra.h"

const double scale2=1./(2*PI);


template <int D,int DD>
//void printstart(string ofolder,string ev,string &out1, string &out2,  SPH<D,DD> &sph)
void printstart(string ofolder,string ev,string &out1, string &out2,  SPH<D,DD> &sph)
/* {
    std::string negc="neg";
    std::string viscosity_prefix;	
	
    // Define prefix depending on viscosity type
    if (sph.typ==0)      viscosity_prefix = "i";
    else if (sph.typ==1) viscosity_prefix = "bv"; 
	else if (sph.typ==2) viscosity_prefix = "sv"; 
	else if (sph.typ==3) viscosity_prefix = "sbv"; 
	else if (sph.typ==4) viscosity_prefix = "bsqsbv";
	
    std::string negative_suffix;
    // Append suffix
    if (sph.neg == negc) {
        negative_suffix = "_neg";
    }
    else {
        negative_suffix = "";
    }
    
    out2 = ofolder + "/event_"
                   + ev
                   + "/corrected_spectra/"
                   + viscosity_prefix
                   + "_dNdpdphi"
                   + negative_suffix
                   + ".dat";

	ofstream OUT2;
	OUT2.open(out2.c_str() );
	if (!OUT2.is_open()) {
		cout << "Error: cannot open corrected spectra file!" << endl;
		cout << "Filename: " << out2 << endl;
		exit(1);
	}
	
	OUT2.close();


    out1 = ofolder  + "/event_"
                    + ev
                    + "/uncorrected_spectra/"
                    + viscosity_prefix
                    + "_dNdpdphi"
                    + negative_suffix
                    + ".dat";

    ofstream OUT;
    OUT.open(out1.c_str() );
    if (!OUT.is_open()) {
        cout << "Error: cannot open spectra file!" << endl;
        exit(1);
    }

    OUT.close();
}
 */

{
	string negc="neg";
	
	string vtyp;	
	
	if (sph.typ==1) vtyp="bvc"; 
	else if (sph.typ==2) vtyp="svc"; 
	else if (sph.typ==3) vtyp="sbvc"; 
	else if (sph.typ==4) vtyp="bsqsbvc";
	else if (sph.typ==5) vtyp="v2c";
	
	out2=ofolder+"/ev"+vtyp+"_dNdphidpp.dat";
	if (sph.neg==negc) out2=ofolder+"/ev"+ev+vtyp+"_dNdphidpp_neg.dat";
	ofstream OUT2;
	OUT2.open(out2.c_str() );
	if (!OUT2.is_open())
	{
		cout << "Error: cannot open corrected spectra file!" << endl;
		cout << "Filename: " << out2 << endl;
		exit(1);
	}
	
	OUT2.close();



if (sph.typ==0) vtyp="i"; 
	  else if (sph.typ==1) vtyp="bv"; 
	  else if (sph.typ==2) vtyp="sv"; 
	  else if (sph.typ==3) vtyp="sbv";
      else if (sph.typ==4) vtyp="bsqsbv";	  
	  else if (sph.typ==5) vtyp="v2";
	  
	  
	  out1=ofolder+"/ev"+vtyp+"_dNdphidpp.dat";
	  if (sph.neg==negc) out1=ofolder+"/ev"+ev+vtyp+"_dNdphidpp_neg.dat";
	  
	  ofstream OUT;
	  OUT.open(out1.c_str() );
	  if (!OUT.is_open())
	    {
	      cout << "Error: cannot open spectra file!" << endl;
	      exit(1);
	    }
	  
	  OUT.close();
   }

template <int D,int DD>
void print(string out1, int h, list & l , SPH<D,DD> sph)
{
	
	ofstream OUT;
	OUT.open(out1.c_str(), ios::out | ios::app );
	if (!OUT.is_open())
	{
		cout << "Error: cannot open final spectra file!" << endl;
		exit(1);
	}
	
	OUT << sph.had[h].id << endl;
	for(int i=0;i<l.phimax;i++)
	{
	for(int ps=0;ps<l.pTmax;ps++)
	{	
	OUT <<  l.dNdpdphi.x[ps][i] << " " ;
	}
	OUT <<  endl;
	}
	
	OUT.close();

}

template <int D,int DD>
void printc(string out1, int h, list & l , SPH<D,DD> sph)
{
	
	ofstream OUT;
	OUT.open(out1.c_str(), ios::out | ios::app );
	if (!OUT.is_open())
	{
		cout << "Error: cannot open final corrected spectra file!" << endl;
		exit(1);
	}
	
	OUT << sph.had[h].id << endl;
	for(int i=0;i<l.phimax;i++)
	{
	for(int ps=0;ps<l.pTmax;ps++)
	{	
	OUT <<  l.dNdpdphic.x[ps][i] << " " ;
	}
	OUT <<  endl;
	}
	
	OUT.close();

}


template <int D,int DD>
  void nprint(string out1, int h, list & l , SPH<D,DD> sph,int nsub)
{
  
  ofstream OUT;
  OUT.open(out1.c_str(), ios::out | ios::app );
  if (!OUT.is_open())
    {
      cout << "Error: cannot open final spectra file!" << endl;
      exit(1);
    }
  
  OUT << nsub << endl;
  for(int i=0;i<l.phimax;i++)
    {
      for(int ps=0;ps<l.pTmax;ps++)
	{
	  OUT <<  l.dNdpdphi.x[ps][i] << " " ;
	}
      OUT <<  endl;
    }
  
  OUT.close();

}

template <int D,int DD>
  void nprintc(string out1, int h, list & l , SPH<D,DD> sph,int nsub)
{
  
  ofstream OUT;
  OUT.open(out1.c_str(), ios::out | ios::app );
  if (!OUT.is_open())
    {
      cout << "Error: cannot open final corrected spectra file!" << endl;
      exit(1);
    }
  
  OUT << nsub << endl;
  for(int i=0;i<l.phimax;i++)
    {
      for(int ps=0;ps<l.pTmax;ps++)
	{
	  OUT <<  l.dNdpdphic.x[ps][i] << " " ;
	}
      OUT <<  endl;
    }
  
  OUT.close();

}


#endif

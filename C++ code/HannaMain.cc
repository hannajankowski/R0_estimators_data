#include <math.h>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <list>
#include <algorithm>
#include <functional>
#include "rannumgen.h"
#include "HannaClass.h"

using namespace std;

long int seed = -1;

//Note that C++ compilers may not compile this code the same on all operating systems, or in C++ use programs. 

//Start with one initial infected, and 10000 susceptibles
int InitSuscept = 10000;
int InitExpose = 0;
int InitInfect = 1;
int InitSymp = 0;
int InitRecov = 0;

double pop_amp = 1.0;

//These need to be modified to run SIR, SEIR, or SEIR models. 
double EImean = 1.0; //E to I      
double IYmean = 2.0; //I to Y
double YRmean = 2.0; //I to R
double Infmean = 30000.0;
double infmean = Infmean/(InitSuscept + InitExpose + InitInfect + InitSymp + InitRecov);

#include "HannaFunction.h"

int main(int argc, char* argv[]){        
	int numsuscept = InitSuscept;
	int numexpose = InitExpose;
	int numinfect = InitInfect;
    	int numsymp = InitSymp;
	int numrecov = InitRecov;

	double pop = numsuscept + numexpose + numinfect + numsymp + numrecov;

	infmean = (Infmean)/pop;

	double currenttime = 0.0, prevtime = 0.0, printtime = 0.0;
	double printinc = 0.1;

	SList suscept;
	EList expose;
	IList infect;
    	YList symp;
	RList recov;

	SList::iterator spos;
	EList::iterator epos;
	IList::iterator ipos;
    	YList::iterator ypos;
	RList::iterator rpos;

	double etime = 0.0;
	double itime = 0.0;
    	double ytime = 0.0;
	double rtime = 0.0;
    	double mintime = 0.0;
    	int patnum = 0;

	seed = atoi(argv[1]);
	char fname[100];
	sprintf(fname, "Hannacomplete%dseed.m",-seed);

#ifdef USE_STD
	std::ofstream ofile(fname);
#else
	ofstream ofile(fname);
#endif

	initializesuscept(suscept,numsuscept);
	initializeexpose(expose,numexpose);
	initializeinfect(infect,numinfect,infmean);
    	initializesymp(symp,numsymp,infmean);
	initializerecov(recov,numrecov);

	numsuscept = suscept.size();
	numexpose = expose.size();
	numinfect = infect.size();
    	numsymp = symp.size();
	numrecov = recov.size();

	printtime = printtime + printinc;

    	ofile<<"dataarray = ["<<currenttime<<" "<<printtime<<" ";
    	ofile<<numsuscept<<" "<<suscept.size()<<" ";
    	ofile<<numexpose<<" "<<expose.size()<<" ";
    	ofile<<numinfect<<" "<<infect.size()<<" ";
    	ofile<<numsymp<<" "<<symp.size()<<" ";
    	ofile<<numrecov<<" "<<recov.size()<<endl;

	while(currenttime < 200.0 && (numexpose + numinfect + numsymp)>0){
		//cout<<currenttime<<" "<<endl;	
		if(currenttime >= printtime){
			ofile<<currenttime<<" "<<printtime<<" ";
            ofile<<numsuscept<<" "<<suscept.size()<<" ";
            ofile<<numexpose<<" "<<expose.size()<<" ";
            ofile<<numinfect<<" "<<infect.size()<<" ";
            ofile<<numsymp<<" "<<symp.size()<<" ";
            ofile<<numrecov<<" "<<recov.size()<<endl;
            
            cout<<currenttime<<" "<<printtime<<" ";
            cout<<numsuscept<<" "<<suscept.size()<<" ";
            cout<<numexpose<<" "<<expose.size()<<" ";
            cout<<numinfect<<" "<<infect.size()<<" ";
            cout<<numsymp<<" "<<symp.size()<<" ";
            cout<<numrecov<<" "<<recov.size()<<endl;
            
			printtime = printtime + printinc;
		}

		spos = suscept.begin();
		epos = expose.begin();
		ipos = infect.begin();
        	ypos = symp.begin();
		rpos = recov.begin();

		pop = numsuscept + numexpose + numinfect + numsymp + numrecov;
        	infmean = Infmean/pop;

		if(numexpose>0){
			etime = epos->geteventtime();
		}
		else{
			etime = currenttime + 10000.0;
		}
		if(numinfect>0){
			itime = ipos->geteventtime();
		}
		else{
			itime = currenttime + 10000.0;
		}
        	if(numsymp>0){
            	ytime = ypos->geteventtime();
        	}
        	else{
            	ytime = currenttime + 10000.0;
        	}
        
		mintime = getmintime(etime,itime,ytime);
/*        	cout<<etime<<" "<<itime<<" "<<ytime<<endl;
        	cout<<mintime<<" "<<currenttime<<" "<<printtime<<" ";
        	cout<<numsuscept<<" "<<suscept.size()<<" ";
        	cout<<numexpose<<" "<<expose.size()<<" ";
        	cout<<numinfect<<" "<<infect.size()<<" ";
        	cout<<numsymp<<" "<<symp.size()<<" ";
        	cout<<numrecov<<" "<<recov.size()<<endl;*/



		if(mintime == etime){
           		currenttime = etime;
           		insertinfect(infect, currenttime, infmean, etime, epos->getSEtime(), epos->getpatnum());
	           	numinfect++;
	           	numexpose--;
	           	expose.erase(epos);
	           	epos = expose.begin();
        	}
        	else if(mintime == itime){
      	     	currenttime = itime;
	           	if(ipos->geteventtype() == 'i'){
              		//cout<<"infect\n";
		            if(ran3(&seed)<(numsuscept+0.0)/(pop+0.0)){
                			if(numsuscept>0){
                    			patnum = takesuscept(suscept, numsuscept);
                    			numsuscept--;
                    			insertexpose(expose, currenttime, itime, patnum);
                    			numexpose++;
                		}
              	}
              	getnewInfecttime(infect, currenttime, infmean, ipos->getIYtime(), ipos->getEItime(), ipos->getSEtime(), ipos->getpatnum());
           }
           else{
              insertsymp(symp, currenttime, infmean, ipos->getIYtime(), ipos->getEItime(), ipos->getSEtime(), ipos->getpatnum());
              numsymp++;
              numinfect--;
           }
           infect.erase(ipos);
           ipos = infect.begin();
        }
        else if(mintime == ytime){
           currenttime = ytime;
           if(ypos->geteventtype() == 'i'){
           		if(ran3(&seed)<(numsuscept+0.0)/(pop+0.0)){
                		if(numsuscept>0){
                    		patnum = takesuscept(suscept, numsuscept);
			            numsuscept--;
                    		insertexpose(expose, currenttime, ytime, patnum);
		                  numexpose++;
                		}
              	}
              	getnewYnfecttime(symp, currenttime, infmean, ypos->getYRtime(), ypos->getIYtime(), ypos->getEItime(), ypos->getSEtime(), ypos->getpatnum());
           }
           else{
              	insertrecov(recov, ytime, ypos->getIYtime(), ypos->getEItime(), ypos->getSEtime(), ypos->getpatnum());
              	numrecov++;
              	numsymp--;
           }
           symp.erase(ypos);
           ypos = symp.begin();
        }
        else{
           cout<<"OOOOOOOOOO*******************************************************"<<endl;
        }
    }        

    ofile<<currenttime<<" "<<printtime<<" ";
    ofile<<numsuscept<<" "<<suscept.size()<<" ";
    ofile<<numexpose<<" "<<expose.size()<<" ";
    ofile<<numinfect<<" "<<infect.size()<<" ";
    ofile<<numsymp<<" "<<symp.size()<<" ";
    ofile<<numrecov<<" "<<recov.size()<<"]"<<endl;
    
    rpos = recov.begin();
    ofile<<"recovarray = [";
    for(int j=0;j<=recov.size();j++){
       ofile<<rpos->getpatnum()<<" ";
       ofile<<rpos->getYRtime()<<" ";
       ofile<<rpos->getIYtime()<<" ";
       ofile<<rpos->getEItime()<<" ";
       ofile<<rpos->getSEtime()<<endl;
       ++rpos;
    }
    ofile<<"]";

}

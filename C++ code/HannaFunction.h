void initializesuscept(SList& suscept, int numsuscept){

	for(int i=0;i<numsuscept;i++){
	   suscept.push_back(Suscept(i+1));
	}
	suscept.sort();
}

void initializeexpose(EList& expose, int numexpose){
	double etime = 0.0;
	
	for (int i=0;i<numexpose;i++){
		etime = expdev(&seed)*EImean;
		expose.push_back(Expose(etime,0,i));
	}	
	expose.sort();
}


void initializeinfect(IList& infect, int numinfect, double infmean){
	double inftime = 0.0;
	double iytime = 0.0;

	for(int i=0;i<numinfect;i++){
		inftime = expdev(&seed)*infmean;
		iytime = expdev(&seed)*IYmean;

		if(inftime < iytime){
			infect.push_back(Infect(inftime,iytime,0,0,i,'i'));
		}
		else{
			infect.push_back(Infect(iytime,iytime,0,0,i,'y'));
		}
    }
	infect.sort();
}

void initializesymp(YList& symp, int numsymp, double infmean){
	double inftime = 0.0;
	double yrtime = 0.0;

	for(int i=0;i<numsymp;i++){
		inftime = expdev(&seed)*infmean;
		yrtime = expdev(&seed)*YRmean;

		if(inftime < yrtime){
			symp.push_back(Symp(inftime,yrtime,0,0,0,i,'i'));
		}
		else{
			symp.push_back(Symp(yrtime,yrtime,0,0,0,i,'r'));
		}
    }
	symp.sort();
}

void initializerecov(RList& recov, int numrecov){

	for(int i=0;i<numrecov;i++){
		recov.push_back(Recov(0,0,0,0,0));
	}
	recov.sort();
}


double getmintime(double etime, double itime, double ytime){
	double mintime = 0;
	double timearray[] = {etime, itime, ytime};//, vtime};

	mintime = timearray[0];
	for(int i=1; i<3;i++){
		if(timearray[i]<mintime){
			mintime = timearray[i];
		}
	}
	return mintime;
}


void insertexpose(EList& expose, double currenttime, double setime, int patnum){
	double etime = expdev(&seed)*(EImean) + currenttime;

	EList::iterator epos=find_if(expose.begin(),expose.end(),std::bind2nd(std::greater<Expose>(),Expose(etime,setime,patnum)));
	expose.insert(epos,Expose(etime,setime,patnum));
}


void insertinfect(IList& infect, double currenttime, double infmean, double eitime, double setime, int patnum){
	double inftime = expdev(&seed)*infmean + currenttime;
	double iytime = expdev(&seed)*IYmean + currenttime;

	if(inftime < iytime){
		IList::iterator ipos=find_if(infect.begin(),infect.end(), std::bind2nd(std::greater<Infect>(), Infect(inftime,iytime,eitime,setime,patnum,'i')));
		infect.insert(ipos,Infect(inftime,iytime,eitime,setime,patnum,'i'));
	}
	else{
		IList::iterator ipos=find_if(infect.begin(),infect.end(), std::bind2nd(std::greater<Infect>(), Infect(iytime,iytime,eitime,setime,patnum,'i')));
		infect.insert(ipos,Infect(iytime,iytime,eitime,setime,patnum,'y'));
	}
}

void insertsymp(YList& symp, double currenttime, double infmean, double iytime, double eitime, double setime, int patnum){
    double inftime = expdev(&seed)*infmean + currenttime;
	double yrtime = expdev(&seed)*YRmean + currenttime;

	if(inftime < yrtime){
		YList::iterator ipos=find_if(symp.begin(),symp.end(), std::bind2nd(std::greater<Symp>(), Symp(inftime,yrtime,iytime,eitime,setime,patnum,'i')));
		symp.insert(ipos,Symp(inftime,yrtime,iytime,eitime,setime,patnum,'i'));
	}
	else{
		YList::iterator ipos=find_if(symp.begin(),symp.end(), std::bind2nd(std::greater<Symp>(), Symp(yrtime,iytime,iytime,eitime,setime,patnum,'y')));
		symp.insert(ipos,Symp(yrtime,iytime,iytime,eitime,setime,patnum,'r'));
	}
}


void insertrecov(RList& recov, double yrtime, double iytime, double eitime, double setime, int patnum){

	RList::iterator rpos=find_if(recov.begin(), recov.end(), std::bind2nd(std::greater<Recov>(), Recov(yrtime,iytime,eitime,setime,patnum)));
	recov.insert(rpos,Recov(yrtime,iytime,eitime,setime,patnum));
}


void getnewInfecttime(IList& infect, double currenttime, double infmean, double iytime, double eitime, double setime, int patnum){
	double inftime = expdev(&seed)*(infmean) + currenttime;

    if(inftime < iytime){
		IList::iterator ipos=find_if(infect.begin(),infect.end(), std::bind2nd(std::greater<Infect>(), Infect(inftime,iytime,eitime,setime,patnum,'i')));
		infect.insert(ipos,Infect(inftime,iytime,eitime,setime,patnum,'i'));
	}
	else{
		IList::iterator ipos=find_if(infect.begin(),infect.end(), std::bind2nd(std::greater<Infect>(), Infect(iytime,iytime,eitime,setime,patnum,'i')));
		infect.insert(ipos,Infect(iytime,iytime,eitime,setime,patnum,'y'));
	}
}


void getnewYnfecttime(YList& symp, double currenttime, double infmean, double yrtime, double iytime, double eitime, double setime, int patnum){
	double inftime = expdev(&seed)*(infmean) + currenttime;

    if(inftime < yrtime){
		YList::iterator ipos=find_if(symp.begin(),symp.end(), std::bind2nd(std::greater<Symp>(), Symp(inftime,yrtime,iytime,eitime,setime,patnum,'i')));
		symp.insert(ipos,Symp(inftime,yrtime,iytime,eitime,setime,patnum,'i'));
	}
	else{
		YList::iterator ipos=find_if(symp.begin(),symp.end(), std::bind2nd(std::greater<Symp>(), Symp(yrtime,iytime,iytime,eitime,setime,patnum,'i')));
		symp.insert(ipos,Symp(yrtime,yrtime,iytime,eitime,setime,patnum,'r'));
	}
}

int takesuscept(SList& suscept, int numsuscept){
	double whattake = ran3(&seed)*numsuscept;
    int patnum = 0;
	SList::iterator spos=suscept.begin();

	for(int i=1;i<whattake;i++){
		++spos;
	}          
    patnum = spos->getpatnum();      
	suscept.erase(spos);
    return patnum;
}




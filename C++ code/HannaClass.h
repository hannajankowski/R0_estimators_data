class Suscept{
public:
	Suscept(int);
	int getpatnum() const{return patnum;};
private:
	int patnum;
};
typedef std::list<Suscept> SList;

Suscept::Suscept(int pnum){
	patnum = pnum;
};
bool operator>(const Suscept& S1, const Suscept& S2)
{return S1.getpatnum() > S2.getpatnum();}
bool operator<(const Suscept& S1, const Suscept& S2)
{return S1.getpatnum() < S2.getpatnum();}



class Expose{
public:
	Expose(double,double,int);
	double geteventtime() const {return eventtime;};
    double getSEtime() const {return SEtime;};
    int getpatnum() const {return patnum;};
private:
	double eventtime;
    double SEtime;
	int patnum;
};
typedef std::list<Expose> EList;

Expose::Expose(double etime, double setime, int pnum){
	eventtime = etime;
    SEtime = setime;
	patnum = pnum;
};
bool operator>(const Expose& S1, const Expose& S2)
{return S1.geteventtime() > S2.geteventtime();}
bool operator<(const Expose& S1, const Expose& S2)
{return S1.geteventtime() < S2.geteventtime();}






class Infect{
public:
	Infect(double, double, double, double, int, char);
	double geteventtime() const {return eventtime;};
    double getIYtime() const {return IYtime;};
    double getEItime() const {return EItime;};
    double getSEtime() const {return SEtime;};
	int getpatnum() const {return patnum;};
    char geteventtype() const {return eventtype;};
private:
	double eventtime;
    double IYtime;
    double EItime;
	double SEtime;
    int patnum;
    char eventtype;
};
typedef std::list<Infect> IList;

Infect::Infect(double etime, double iytime, double eitime, double setime, int pnum, char etype){
	eventtime = etime;
    IYtime = iytime;
    EItime = eitime;
    SEtime = setime;
	patnum = pnum;
    eventtype = etype;
};

bool operator>(const Infect& S1, const Infect& S2)
{return S1.geteventtime() > S2.geteventtime();}
bool operator<(const Infect& S1, const Infect& S2)
{return S1.geteventtime() < S2.geteventtime();}






class Symp{
public:
	Symp(double,double,double,double,double,int,char);
	double geteventtime() const {return eventtime;};
    double getYRtime() const {return YRtime;};
    double getEItime() const {return EItime;};
    double getIYtime() const {return IYtime;};
    double getSEtime() const {return SEtime;};
	int getpatnum() const {return patnum;};
    char geteventtype() const {return eventtype;};
private:
	double eventtime;
    double YRtime;
    double EItime;
    double SEtime;
    double IYtime;
	int patnum;
    char eventtype;
};
typedef std::list<Symp> YList;

Symp::Symp(double ytime, double yrtime, double iytime, double eitime, double setime, int pnum, char etype){
	eventtime = ytime;
    YRtime  = yrtime;
    EItime = eitime;
    SEtime = setime;
    IYtime = iytime;
	patnum = pnum;
    eventtype = etype;
};

bool operator>(const Symp& S1, const Symp& S2)
{return S1.geteventtime() > S2.geteventtime();}
bool operator<(const Symp& S1, const Symp& S2)
{return S1.geteventtime() < S2.geteventtime();}



class Recov{
public:
	Recov(double, double, double, double, int);
	int getpatnum() const {return patnum;};
    double getYRtime() const {return YRtime;};
    double getEItime() const {return EItime;};
    double getIYtime() const {return IYtime;};
    double getSEtime() const {return SEtime;};
private:
	int patnum;
    double YRtime;
    double IYtime;
    double EItime;
    double SEtime;
};

typedef std::list<Recov> RList;

Recov::Recov(double yrtime, double iytime, double eitime, double setime, int pnum){
	patnum = pnum;
    YRtime = yrtime;
    IYtime = iytime;
    EItime = eitime;
    SEtime = setime;
};

bool operator>(const Recov& S1, const Recov& S2)
{return S1.getYRtime() > S2.getYRtime();}
bool operator<(const Recov& S1, const Recov& S2)
{return S1.getYRtime() < S2.getYRtime();}



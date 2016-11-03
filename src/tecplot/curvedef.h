struct node{
	double Xdata;
	double Ydata;
	double Zdata;
	double ph;
	double totht;
	node* next;
};

typedef node* NodePtr;

class curve{
	double xbounds[2];
	double ybounds[2];
	double zbounds[2];
	double phbounds[2];
	int length;

	void updatelength();
	void push(double,double,double,double,double );
	void push(double,double,double,double,double,int);
	public:

	NodePtr head;
	// appends node atthe end of list
	void added(double,double,double,double,double);    
	// inserts node in ascending order of 1st and 2nd arguments 
	void arrange(double,double,double,double,double);
	// brings together nodes having same Xcoords
	void reorder(void);		     
	void merge(curve);
	// removes node at specified location
	void remove(double,double,double);   
	// deletes nth location
	void deloc(int);		
	// writes file 'precessed.txt' with all node contents
	void writeout(void);		     

	void getlength(int & );
	// writes all nodes and returns length of list
	void getinfo(int & );		     
	// just prints out all info
	void getinfo(void);		     
	void getextents(double* ,double* ,double* ,double* );

	// returns points 
	friend curve getnearest(curve & megalist,double entry ,double error );
	// fills node array PH by interpolating points
	friend void SetPile(curve & megalist,node [][100],int,int);
	friend double GetDistance(double X1,double Y1,double,double);

	curve(void);
};

curve::curve(void){
	xbounds[0]=0.0;
	xbounds[1]=0.0;
	ybounds[0]=0.0;
	ybounds[1]=0.0;
	zbounds[0]=0.0;
	zbounds[1]=0.0;
	phbounds[0]=0.0;
	phbounds[1]=0.0;
	length=0;
	head=NULL;
}

double GetDistance(double X1,double Y1,double X2,double Y2);


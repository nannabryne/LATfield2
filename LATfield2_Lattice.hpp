#ifndef LATFIELD2_LATTICE_HPP
#define LATFIELD2_LATTICE_HPP


/*! \file LATfield2_Lattice.hpp
 \brief LATfield2_Lattice.hpp contains the class Lattice definition.
 \author David Daverio, Neil Bevis, edited by Nanna Bryne
 */ 






//CONSTANTS=====================

int Lattice::initialized = 1;        //Status flag for initialized

//CONSTRUCTORS===================

Lattice::Lattice()
{
	status_=0;
	arch_saved_=false;
}

Lattice::Lattice(int dim, const int* size, int halo)
{
	status_=0;
	arch_saved_=false;
	this->initialize(dim, size, halo);
}  


Lattice::Lattice(int dim, const int size, int halo)
{
	status_=0;
	arch_saved_=false;
	int* sizeArray=new int[dim];
	for(int i=0; i<dim; i++) { sizeArray[i]=size; }
	this->initialize(dim, sizeArray, halo);
	delete[] sizeArray;
}



//DESTRUCTOR=========================

Lattice::~Lattice() 
{ 
	if((status_ & initialized) > 0)
    {
		delete[] size_;
		delete[] sizeLocal_;
		delete[] jump_;
        delete[] sizeLocalAllProcDim0_;
        delete[] sizeLocalAllProcDim1_;
    }
}
//INITIALIZE=========================

void Lattice::initialize(int dim, const int size, int halo)
{
	int* sizeArray=new int[dim];
	for(int i=0; i<dim; i++) { sizeArray[i]=size; }
	this->initialize(dim, sizeArray, halo);
}

void Lattice::initialize(int dim, const int* size, int halo)
{
	int i;
	
	if((status_ & initialized) > 0)
    {
		delete[] size_;
		delete[] sizeLocal_;
		delete[] jump_;
        delete[] sizeLocalAllProcDim0_;
        delete[] sizeLocalAllProcDim1_;
    }
	//Store input lattice properties
	dim_ =dim;
	size_=new int[dim_];
	for(i=0;i<dim_;i++) size_[i]=size[i];
	halo_=halo;
	
	
	//Calculate local size
	sizeLocal_=new int[dim_];
	sizeLocal_[dim_-1]=int(ceil( (parallel.grid_size()[0]-parallel.grid_rank()[0])*size_[dim_-1]/float(parallel.grid_size()[0]) ));
	sizeLocal_[dim_-1]-=int(ceil((parallel.grid_size()[0]-parallel.grid_rank()[0]-1)*size_[dim_-1]/float(parallel.grid_size()[0]) ));
	sizeLocal_[dim_-2]=int(ceil( (parallel.grid_size()[1]-parallel.grid_rank()[1])*size_[dim_-2]/float(parallel.grid_size()[1]) ));
	sizeLocal_[dim_-2]-=int(ceil((parallel.grid_size()[1]-parallel.grid_rank()[1]-1)*size_[dim_-2]/float(parallel.grid_size()[1]) ));
	for(i=0;i<dim_-2;i++) sizeLocal_[i]=size_[i];
    
    sizeLocalAllProcDim0_ = new int[parallel.grid_size()[0]];
	sizeLocalAllProcDim1_ = new int[parallel.grid_size()[1]];
    
	for(i=0;i<parallel.grid_size()[0];i++)
    {
	    sizeLocalAllProcDim0_[i] = int(ceil( (parallel.grid_size()[0]-i)*size_[dim_-1]/float(parallel.grid_size()[0]) ));
	    sizeLocalAllProcDim0_[i] -= int(ceil((parallel.grid_size()[0]-i-1)*size_[dim_-1]/float(parallel.grid_size()[0]) ));	
    }
	for(i=0;i<parallel.grid_size()[1];i++)
    {
	    sizeLocalAllProcDim1_[i] = int(ceil( (parallel.grid_size()[1]-i)*size_[dim_-2]/float(parallel.grid_size()[1]) ));
	    sizeLocalAllProcDim1_[i] -= int(ceil((parallel.grid_size()[1]-i-1)*size_[dim_-2]/float(parallel.grid_size()[1]) ));	
    }
	
	//Calculate index jumps
	jump_=new long[dim_];
	jump_[0]=1;
	for(i=1;i<dim_;i++) jump_[i]=jump_[i-1]*(sizeLocal_[i-1]+2*halo_);
	
	//Calculate number of sites in lattice
	sitesLocal_=1;
	sitesLocalGross_=1;
	for(i=0;i<dim_;i++)
    {
		sitesLocal_*=sizeLocal_[i];
		sitesLocalGross_*=sizeLocal_[i]+(2*halo_);
    }
	sites_=sitesLocal_;
	sitesGross_=sitesLocalGross_;
	parallel.sum(sites_);
	parallel.sum(sitesGross_);
	
	//Calculate index of first and last local sites on lattice
	siteFirst_=0;
	siteLast_=sitesLocalGross_-1;
	for(i=0;i<dim_;i++)
    {
		siteFirst_+=jump_[i]*halo_;
		siteLast_-=jump_[i]*halo_;
    }
	
	////calculate coordSkip
	
	
	
	//Get each processor to tell the others in his dim0_group its local sizeLocal_[dim-1])
	int* sizes_dim0=new int[parallel.grid_size()[0]];
	for(i=0;i<parallel.grid_size()[0];i++)
	{
		if(i==parallel.grid_rank()[0]) { sizes_dim0[i]=sizeLocal_[dim_-1]; }
		parallel.broadcast_dim0(sizes_dim0[i],i);
	}
	//Sum up sizes for the processors of less than or equal rank
	coordSkip_[0]=0;
	for(i=0; i<parallel.grid_rank()[0]; i++)coordSkip_[0]+=sizes_dim0[i];
	
	
	
	//Get each processor to tell the others in his dim1_group its local sizeLocal_[dim-2])
	int* sizes_dim1=new int[parallel.grid_size()[1]];
	for(i=0;i<parallel.grid_size()[1];i++)
	{
		if(i==parallel.grid_rank()[1]) { sizes_dim1[i]=sizeLocal_[dim_-2]; }
		parallel.broadcast_dim1(sizes_dim1[i],i);
	}
	//Sum up sizes for the processors of less than or equal rank
	coordSkip_[1]=0;
	for(i=0; i<parallel.grid_rank()[1]; i++)coordSkip_[1]+=sizes_dim1[i];
	
	
	
	
	////calculate sitesSkip : sitesskip used for fastread , fastload (function witch need to be coded :-) )
	
	sitesSkip_=coordSkip_[0];
	for(i=0;i<dim_-1;i++)sitesSkip_*=size_[i];
	
	long siteSkiptemp= coordSkip_[1]*sizeLocal_[dim_-1];
	for(i=0;i<dim_-2;i++)siteSkiptemp*=size_[i];
	
	sitesSkip_+=siteSkiptemp;
	
	//calculate sitesSkip2d : 
	
	int* sizes1 = new int[parallel.grid_size()[1]];
	int* sizes0 = new int[parallel.grid_size()[0]];
	long* offset1 = new long[parallel.grid_size()[1]];
	long* offset0 = new long[parallel.grid_size()[0]];
	
	int b=1;
	int n;
	for(i=0;i<dim_-2;i++)b*=sizeLocal_[i];
	
	//calulate offset in dim-2
	for(i=0;i<parallel.grid_size()[1];i++)sizes1[i]=sizes_dim1[i] * b;
	for(n=0;n<parallel.grid_size()[1];n++)offset1[n]=0;
	
	for(n=1;n<parallel.grid_size()[1];n++)
	{
		for(i=0;i<n;i++)offset1[n]+=sizes1[i];
	}
	
	//calulate offset in dim-1
	for(i=0;i<parallel.grid_size()[0];i++)sizes0[i]=size_[dim_-2] * sizes_dim0[i] * b;
	for(n=0;n<parallel.grid_size()[0];n++)offset0[n]=0;
	
	for(n=1;n<parallel.grid_size()[0];n++)
	{
		for(i=0;i<n;i++)offset0[n]+=sizes0[i];
	}
	
	sitesSkip2d_ = offset0[parallel.grid_rank()[0]] + offset1[parallel.grid_rank()[1]];
	
	//Set status
	status_ = status_ | initialized;
	
	//Free memory
	delete[] sizes_dim0;
	delete[] sizes_dim1;
  
  // WV: Free more memory!
  delete[] sizes1;
  delete[] sizes0;
  delete[] offset1;
  delete[] offset0;

	
}
#ifdef FFT3D
void Lattice::initializeRealFFT(Lattice & lat_real, int halo)
{
	
	if(lat_real.dim()!=3)
	{
		if(parallel.isRoot())
		{
			cerr<<"Latfield2d::Lattice::initializeRealFFT : fft curently work only for 3d cubic lattice"<<endl;
			cerr<<"Latfield2d::Lattice::initializeRealFFT : coordinate lattice have not 3 dimensions"<<endl;
			cerr<<"Latfield2d : Abort Process Requested"<<endl;
			
		}
		parallel.abortForce();
	}
	
	if(lat_real.size(0)!=lat_real.size(1) | lat_real.size(2)!=lat_real.size(1))
	{
		if(parallel.isRoot())
		{
			cerr<<"Latfield2d::Lattice::initializeRealFFT : fft curently work only for 3d cubic lattice"<<endl;
			cerr<<"Latfield2d::Lattice::initializeRealFFT : coordinate lattice is not cubic"<<endl;
			cerr<<"Latfield2d : Abort Process Requested"<<endl;
			
		}
		parallel.abortForce();
	}
	
	int lat_size[3];
	
	lat_size[0]=lat_real.size(0)/2+1;
	lat_size[1]=lat_real.size(0);
	lat_size[2]=lat_real.size(0);
	
	this->initialize(3, lat_size, halo);
}
void Lattice::initializeComplexFFT(Lattice & lat_real, int halo)
{
	
	if(lat_real.dim()!=3)
	{
		if(parallel.isRoot())
		{
			cerr<<"Latfield2d::Lattice::initializeRealFFT : fft curently work only for 3d cubic lattice"<<endl;
			cerr<<"Latfield2d::Lattice::initializeRealFFT : coordinate lattice have not 3 dimensions"<<endl;
			cerr<<"Latfield2d : Abort Process Requested"<<endl;
			
		}
		parallel.abortForce();
	}
	
	if(lat_real.size(0)!=lat_real.size(1) | lat_real.size(2)!=lat_real.size(1))
	{
		if(parallel.isRoot())
		{
			cerr<<"Latfield2d::Lattice::initializeRealFFT : fft curently work only for 3d cubic lattice"<<endl;
			cerr<<"Latfield2d::Lattice::initializeRealFFT : coordinate lattice is not cubic"<<endl;
			cerr<<"Latfield2d : Abort Process Requested"<<endl;
			
		}
		parallel.abortForce();
	}
	
	int lat_size[3];
	
	lat_size[0]=lat_real.size(0);
	lat_size[1]=lat_real.size(0);
	lat_size[2]=lat_real.size(0);
	
	this->initialize(3, lat_size, halo);
}
#endif


void Lattice::save_arch(const string filename)
{
	fstream file;
	int p,i;
	
	for( p=0; p<parallel.size(); p++ ) 
	{
		if( parallel.rank()==p )
		{
			if( parallel.rank()==0)
			{
				//Truncate file if first process
				file.open(filename.c_str(), fstream::out | fstream::trunc);
			}
			else
			{
				//Append to file if another process
				file.open(filename.c_str(), fstream::out | fstream::app);
			}
			if(!file.is_open())
			{
				cerr<<"Latfield::Lattice::save_arch - Could not open file for writing"<<endl;
				cerr<<"Latfield::Lattice::save_arch - File: "<<filename<<endl;
				parallel.abortRequest();
			}
			
			if( parallel.rank()==0)
			{
				file<<"# Architerctur of the lattice"<<endl;
				file<<"# Number of dimension :"<<endl;
				file<<this->dim()<<endl;
				file<<"############################"<<endl;
			}
			file<<"############################"<<endl;
			file<<"#  Ranks of processor, world, dim-1, dim-2 :"<<endl;
			file<<parallel.rank()<<" "<<parallel.grid_rank()[0]<<" "<<parallel.grid_rank()[1] <<endl;
			file<<"#  Local size :"<<endl;
			for(i=0;i<this->dim();i++)
			{
				file<<this->sizeLocal()[i]<<endl;
			}
			file<<"############################"<<endl;
			
			
			file.close();
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
	arch_saved_=true;
	if(parallel.rank()==0)cout<<"Architecture saved in : "<<filename<<endl;
	
}

int Lattice::getRank(int* coord)
{
    int n,m;
    int temp;
    bool flag;
    
    temp=0;
    flag = true;
    for(n=0;flag;n++)
    {
        temp+= sizeLocalAllProcDim0_[n];
        if(temp>coord[dim_-1])flag=false;
    }
    n-=1;
    temp=0;
    flag = true;
    for(m=0;flag;m++)
    {
        temp+= sizeLocalAllProcDim1_[m];
        if(temp>coord[dim_-2])flag=false;
    }
    m-=1;
    
    
    return parallel.grid2world(n,m);
    
}
int Lattice::getRankDim0(int coord)
{
    int n;
    int temp;
    bool flag;
    
    temp=0;
    flag = true;
    for(n=0;flag;n++)
    {
        temp+= sizeLocalAllProcDim0_[n];
        if(temp>coord)flag=false;
    }
    return n-1;
}
int Lattice::getRankDim1(int coord)
{
    int m;
    int temp;
    bool flag;
    temp=0;
    flag = true;
    for(m=0;flag;m++)
    {
        temp+= sizeLocalAllProcDim1_[m];
        if(temp>coord)flag=false;
    }
    return  m-1;
}


//MISCELLANEOUS======================

bool Lattice::is_arch_saved() {return arch_saved_;};
int  Lattice::dim() { return dim_; };
int* Lattice::size() { return size_; };
int  Lattice::size(int i) { return size_[i]; }
long  Lattice::sites() { return sites_; }
long  Lattice::sitesGross() { return sitesGross_; }
int  Lattice::halo() { return halo_; }

int* Lattice::sizeLocal() { return sizeLocal_; };
int  Lattice::sizeLocal(int i) { return sizeLocal_[i]; }
long  Lattice::sitesLocal() { return sitesLocal_; }
long  Lattice::sitesLocalGross() { return sitesLocalGross_; }

long  Lattice::jump(int i) { return jump_[i]; }
long  Lattice::sitesSkip() { return sitesSkip_; }
long  Lattice::sitesSkip2d() { return sitesSkip2d_; }
long*  Lattice::coordSkip() { return coordSkip_; }
long Lattice::siteFirst() { return siteFirst_; }
long Lattice::siteLast() { return siteLast_; }

int * Lattice::sizeLocalAllProcDim0(){ return sizeLocalAllProcDim0_; }
int * Lattice::sizeLocalAllProcDim1(){ return sizeLocalAllProcDim1_; }



int Lattice::indexTransform(int* local_coord){

	const int halo = this->halo();
    const int dim = this->dim();

	int index_flat = 0;
    int jump = 1;   // lack of better name
    // for(int n=0; n<dim; n++){
    //     index_flat += jump * ( halo + local_coord[n] );
    //     jump *= ( this->sizeLocal(n) + 2*halo );
    // }

	index_flat = halo + local_coord[0];
	for(int n=1; n<dim; n++){
		jump *= ( this->sizeLocal(n-1) + 2*halo );
        index_flat += jump * ( halo + local_coord[n] );
    }

	return index_flat;
}





#ifndef _OPENMP  // non-OpenMP codes

void Lattice::for_each(std::function<void(Site&)> operation){
	/* (no OpenMP version) usual iterator method */
	Site x(*this);
	for(x.first(); x.test(); x.next()){
		operation(x);
	}
}


void Lattice::for_each(
        std::function<void(Site *)> onsite_operation, 
        Lattice ** other_lattices,
        const int num_lattices,
		string scheme){

	Site * sites = NULL;
	sites = new Site[num_lattices];
	Lattice * lats = NULL;
	lats = new Lattice[num_lattices];
	lats[0] = *this;
	sites[0].initialize(lats[0]);

	for(int l=1; l<num_lattices; l++){
		lats[l] = *other_lattices[l-1];
		sites[l].initialize(lats[l]);
		sites[l].first();
	}

	for(sites[0].first(); sites[0].test(); sites[0].next()){
		onsite_operation(sites);
		for(int l=1; l<num_lattices; l++)sites[l].next(); 
	}

	delete[] sites;

}


void Lattice::for_each(
        std::function<void(Site&, Site *)> onsite_operation, 
        Lattice ** other_lattices,
        const int num_other_lattices,
        string scheme){
		

	// Site * other_sites = NULL;
	// other_sites = new Site[num_other_lattices];
	Site other_sites[num_other_lattices];
	// Site * other_sites;

	for(int l=0; l<num_other_lattices; l++){
		other_sites[l].initialize(*other_lattices[l]);
		other_sites[l].first();
	}

	Site this_site(*this);

	for(this_site.first(); this_site.test(); this_site.next()){
		onsite_operation(this_site, other_sites);
		for(int l=0; l<num_other_lattices; l++)other_sites[l].next();
	}

	// delete[] other_sites;
	// free(other_sites);
		
}




template<typename F> // testing with some object that can point to a lattice
void Lattice::for_each(
        std::function<void(Site&, Site *)> onsite_operation, 
       	F** fields,
        const int num_fields,
        string scheme){
		

	Site other_sites[num_fields];

	for(int l=0; l<num_fields; l++){
		other_sites[l].initialize(fields[l]->lattice());
		other_sites[l].first();
	}

	Site this_site(*this);

	for(this_site.first(); this_site.test(); this_site.next()){
		onsite_operation(this_site, other_sites);
		for(int l=0; l<num_fields; l++)other_sites[l].next();
	}

		
}




#else // OpenMP parallelised loops:


void Lattice::for_each(std::function<void(Site&)> operation){

	/* Should modify to secure against infinite loop */

	/* currently only implemented for dim=3 */
	if(this->dim()==3)
	{	
		if(omp_in_parallel()){ /* if parallel region already exists */

			#pragma omp for collapse(2)
			for(int k=0; k<this->sizeLocal(2); k++)
				for(int j=0; j<this->sizeLocal(1); j++){
					
					int ijk[] = {0,j,k};
					int idx = this->indexTransform(ijk);
					
					Site x(*this, idx);
					for(int i=0; i<this->sizeLocal(0); i++){
						operation(x);
						x.indexAdvance(1);
					}

				}
		}
		else{ /* create new parallel region */
			#pragma omp parallel
			{
				for_each(operation);
			}
		}
	}
	else
	{
		Site x(*this);
		for(x.first(); x.test(); x.next()){
			operation(x);
		}

	}

}	





void Lattice::for_each(
        std::function<void(Site *)> onsite_operation, 
        Lattice ** other_lattices,
        const int num_lattices,
		string scheme){
	
	if(this->dim()==3){

		if(omp_in_parallel()){ // if parallel region exists

			Site * sites = NULL;		sites = new Site[num_lattices];
			Lattice * lats = NULL;		lats = new Lattice[num_lattices];
			
			lats[0] = *this;
	
			for(int l=1; l<num_lattices; l++)lats[l] = *other_lattices[l-1];


			auto loop = [&] (int j_start=0, int k_start=0, const int incr=1){
				#pragma omp for collapse(2)
				for(int k=k_start; k<this->sizeLocal(2); k+=incr)
					for(int j=j_start; j<this->sizeLocal(1); j+=incr){
						int ijk[] = {0,j,k};

						for(int l=0; l<num_lattices; l++)sites[l].initialize(lats[l], (long)lats[l].indexTransform(ijk));

						for(int i=0; i<this->sizeLocal(0); i++){
							onsite_operation(sites); 
							for(int l=0; l<num_lattices; l++)sites[l].indexAdvance(1);
						}
					}
				// implicit barrier
			};

			if(scheme=="arbitrary")loop(0,0,1);
			else if(scheme=="controlled"){ // fix name
				loop(0,0,2); 	/* even-even sites */
				loop(1,1,2); 	/* odd-odd sites */
				loop(1,0,2); 	/* even-odd sites */
				loop(0,1,2);	/* odd-even sites */
			}

			delete[] sites;
			// delete[] lats;
		}
		else{
			#pragma omp parallel
			{
				for_each(onsite_operation, other_lattices, num_lattices, scheme);
			}
		}

	}
	else{ // if dimension is not 3

		Site * sites = NULL;
		sites = new Site[num_lattices];
		Lattice * lats = NULL;
		lats = new Lattice[num_lattices];
		lats[0] = *this;
		sites[0].initialize(lats[0]);

		for(int l=1; l<num_lattices; l++){
			lats[l] = *other_lattices[l-1];
			sites[l].initialize(lats[l]);
			sites[l].first();
		}

		for(sites[0].first(); sites[0].test(); sites[0].next()){
			onsite_operation(sites);
			for(int l=1; l<num_lattices; l++)sites[l].next(); 
		}
	}

		
}




void Lattice::for_each(
        std::function<void(Site&, Site *)> onsite_operation, 
        Lattice ** other_lattices,
        const int num_other_lattices,
        string scheme){
		


	if(this->dim()==3){

		if(omp_in_parallel()){ // if parallel region exists

			Site * other_sites = NULL;
			other_sites = new Site[num_other_lattices];

			for(int l=0; l<num_other_lattices; l++){
				// other_sites[l].initialize(*other_lattices[l]);
				other_sites[l] = Site(*other_lattices[l]);
			}
			
			Site this_site(*this);

			auto loop = [&] (int j_start=0, int k_start=0, const int incr=1){
				#pragma omp for collapse(2)
				for(int k=k_start; k<this->sizeLocal(2); k+=incr)
					for(int j=j_start; j<this->sizeLocal(1); j+=incr){
						int ijk[] = {0,j,k};

						this_site.setIndex((long)this->indexTransform(ijk));

		
						for(int l=0; l<num_other_lattices; l++)other_sites[l].setIndex((long)other_lattices[l]->indexTransform(ijk));

						for(int i=0; i<this->sizeLocal(0); i++){
							COUT << "3\n";
							onsite_operation(this_site, other_sites); 

							this_site.indexAdvance(1);
							for(int l=0; l<num_other_lattices; l++)other_sites[l].indexAdvance(1);
						}
					}
				// implicit barrier
			};


			if(scheme=="arbitrary")loop(0,0,1);
			else if(scheme=="controlled"){ // fix name
				loop(0,0,2); 	/* even-even sites */
				loop(1,1,2); 	/* odd-odd sites */
				loop(1,0,2); 	/* even-odd sites */
				loop(0,1,2);	/* odd-even sites */
			}

			delete[] other_sites;
		}
		else{
			#pragma omp parallel
			{
				for_each(onsite_operation, other_lattices, num_other_lattices, scheme);
			}
		}

	}
	else{ // if dimension is not 3

		Site * other_sites = NULL;
		other_sites = new Site[num_other_lattices];

		for(int l=0; l<num_other_lattices; l++){
			other_sites[l].initialize(*other_lattices[l]);
			other_sites[l].first();
		}

		Site this_site(*this);

		for(this_site.first(); this_site.test(); this_site.next()){
			onsite_operation(this_site, other_sites);
			for(int l=0; l<num_other_lattices; l++)other_sites[l].next(); 
		}

		delete[] other_sites;
	}


}





template<typename F> // testing with some object that can point to a lattice
void Lattice::for_each(
        std::function<void(Site&, Site *)> onsite_operation, 
       	F** fields,
        const int num_fields,
        string scheme){


	if(this->dim()==3){

		if(omp_in_parallel()){ // if parallel region exists

			Site * other_sites = NULL;
			other_sites = new Site[num_fields];

			// Site other_sites[num_fields];

			for(int l=0; l<num_fields; l++){
				// other_sites[l].initialize(fields[l]->lattice());
				other_sites[l] = Site(fields[l]->lattice());
			}
			
			Site this_site(*this);

			auto loop = [&] (int j_start=0, int k_start=0, const int incr=1){
				#pragma omp for collapse(2)
				for(int k=k_start; k<this->sizeLocal(2); k+=incr)
					for(int j=j_start; j<this->sizeLocal(1); j+=incr){
						int ijk[] = {0,j,k};
						this_site.setIndex((long)this->indexTransform(ijk));

						for(int l=0; l<num_fields; l++)other_sites[l].setIndex((long)fields[l]->lattice().indexTransform(ijk));
						

						for(int i=0; i<this->sizeLocal(0); i++){
							onsite_operation(this_site, other_sites); 

							this_site.indexAdvance(1);
							for(int l=0; l<num_fields; l++)other_sites[l].indexAdvance(1);
						}
					}
				// implicit barrier
			};


			if(scheme=="arbitrary")loop(0,0,1);
			else if(scheme=="controlled"){ // fix name
				loop(0,0,2); 	/* even-even sites */
				loop(1,1,2); 	/* odd-odd sites */
				loop(1,0,2); 	/* even-odd sites */
				loop(0,1,2);	/* odd-even sites */
			}

			delete[] other_sites;

		}
		else{
			#pragma omp parallel
			{
				for_each(onsite_operation, fields, num_fields, scheme);
			}
		}

	}
	else{ // if dimension is not 3

		Site other_sites[num_fields];

		for(int l=0; l<num_fields; l++){
			other_sites[l].initialize(fields[l]->lattice());
			other_sites[l].first();
		}

		Site this_site(*this);

		for(this_site.first(); this_site.test(); this_site.next()){
			onsite_operation(this_site, other_sites);
			for(int l=0; l<num_fields; l++)other_sites[l].next();
		}
	}


}
		




#endif // OpenMP vs non-OpenMP codes



void Lattice::for_each(
        std::function<void(Site *)> onsite_operation, 
        Lattice *other_lattice,
        string scheme){

	Lattice * lat_list[1];
	lat_list[0] = other_lattice;
	for_each(onsite_operation, lat_list, 2, scheme);

}



void Lattice::for_each(
        std::function<void(Site&, Site&)> onsite_operation, 
        Lattice *other_lattice,
        string scheme){

	auto op = [&](Site& this_site, Site * other_sites){
		onsite_operation(this_site, other_sites[0]);
	};
	Lattice * lat_list[1];
	lat_list[0] = other_lattice;

	for_each(op, lat_list, 1, scheme);

}





#endif


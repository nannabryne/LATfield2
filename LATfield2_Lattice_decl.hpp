#ifndef LATFIELD2_LATTICE_DECL_HPP
#define LATFIELD2_LATTICE_DECL_HPP

/*! \class Lattice
 
 \brief The Lattice class describe a cartesian mesh (with 2 or more dimensions). The updateHalo method of the Field class generate the periodicity.
 
 
 It store the global and local geometry of the mesh. The last 2 dimension of the lattice are scattered into the MPI processes grid.
 
 
 */
class Lattice
{
public:
    //! Constructor.
    Lattice();
    
    /*!
     Constructor with initialization
     \sa initialize(int dim, const int* size, int halo);
     \param dim : number of dimension
     \param size : array containing the size of each dimension.
     \param halo : size of the halo (ghost cells, same for each dimension)
     */
    Lattice(int dim, const int* size, int halo);
    
    /*!
     Constructor with initialization
     \sa initialize(int dim, const int size, int halo);
     \param dim : number of dimension
     \param size : size of each dimension (same for each dimension)
     \param halo : size of the halo (same for each dimension)
     */
    Lattice(int dim, const int size, int halo);
    
    //! Destructor.
    ~Lattice();
    
    /*!
     Initialization of a dim-dimensional lattice, the size of each dimension is set by the second parameter: int *size. The ghost cell number (halo) is the same for each dimension.
     
     \param dim : number of dimension
     \param size : array containing the size of each dimension.
     \param halo : size of the halo (same for each dimension)
     */
    void initialize(int dim, const int* size, int halo);
    
    /*!
     Initialization of a dim-dimensional lattice, each dimension have the same size. The ghost cell number (halo) is the same for each dimension.
     
     
     \param dim : number of dimension
     \param size : size of each dimension (same for each dimension)
     \param halo : size of the halo (same for each dimension)
     */
    void initialize(int dim, const int size, int halo);
    
#ifdef FFT3D
    
    
    /*!
     Initialization of a lattice for Fourier space in case of real to complex transform. The Fourier space lattice size is defined according to the real space one. The fourier space lattice have "halo" ghost cells in each dimension (which can be different than the halo of the real space lattice).
     \param lat_real : pointer to a real space lattice.
     \param halo : size of the halo (same for each dimension)
     */
    void initializeRealFFT(Lattice & lat_real, int halo);
    
    /*!
     Initialization of a lattice for Fourier space in case of complex to complex transform. The Fourier space lattice size is defined according to the real space one.. The fourier space lattice have "halo" ghost cells in each dimension (which can be different than the halo of the real space lattice).
     \param lat_real : pointer to a real space lattice.
     \param halo : size of the halo (same for each dimension)
     */
    void initializeComplexFFT(Lattice & lat_real, int halo);
#endif
    
    
    
    /*!
     \return int. Number of dimensions of the lattice.
     */
    int  dim();
    /*!
     \return int. Size of the halo (ghost cells).
     */
    int  halo();
    
    /*!
     \return int*. Pointer to the array of the size of each dimension of the lattice.
     */
    int* size();
    
    /*!
     Function which returns the size of a given dimension of the lattice.
     \param direction : asked dimension.
     \return int. Global size of the lattice in the given dimension.
     */
    int  size(int direction);  //Size in a particular dimension
    
    /*!
     \return int*. Pointer to the array of the size of each dimension of the sublattice stored in this MPI process.
     */
    int* sizeLocal();               //Local version
    
    /*!
     Function which returns the size of a given dimension of the sublattice stored in this MPI process.
     \param direction : asked dimension.
     \return int. Global size of the sublattice (of this MPI process) in the given dimension.
     */
    int  sizeLocal(int direction);  //Local version
    
    
    /*!
     \return long. Number of sites on the lattice (excluding halo sites).
     */
    long  sites();              //Number of (global) sites
    /*!
     \return long. Number of sites on the lattice (including halo sites).
     */
    long  sitesGross();         //Number of (global) sites including halo
    
    /*!
     \return long. Number of sites (excluding halo sites) of the sublattice stored in this MPI process.
     */
    long  sitesLocal();              //Local version
    
    /*!
     \return long. Number of sites (including halo sites) of the sublattice stored in this MPI process.
     */
    long  sitesLocalGross();         //Local version
    
    /*!
     \return long. Array index of the first site which is not within the halo.
     */
    long siteFirst();
    
    /*!
     \return long. Array index of the last site which is not within the halo.
     */
    long siteLast();
    
    
    
    /*!
     Function which return the number of data_ array elements to jump to move to the next site in the given direction. (does not take into account the number of component of the fields, therefor should be multiplied by Field.components().) Should not be used by user.
     \param direction : asked direction.
     \return long. Number of array elements to jump.
     */
    long  jump(int direction);       //Number of sites jumped to move in direction
    
    /*!
     \return long. Number of sites before first local site in lattice. Should not be used by users.
     */
    long  sitesSkip();               //Number of sites before first local site in lattice (say in a file)
    
    /*!
     \return long. Number of sites before first local site in lattice. Should not be used by users.
     */
    long  sitesSkip2d();
    
    /*!
     \return *long. Pointer to an array which store the last 2 dimensions coordinate of the first local(in this MPI process) sites. Index 0 is for dim-1, index 1 is for dim-2/
     */
    long*  coordSkip();              //Number to add to coord[dim_-1] to get global value
    
    /*!
     Function which save in serial and in ASCII the global and local description of the Lattice. Usefull to read a file writen by fast_save or fast_write methods of the Field class.
     \param filename : filename of the architectur file.
     */
    void save_arch(const string filename);
    
    /*!
     \return return true if the description of the lattice has been written on disk.
     
     \sa save_arch(const string filename)
     */
    bool is_arch_saved();
    
    int getRank(int* coord) ; //return the world rank of the process who get the lattices site "coord"
    int getRankDim0(int coord) ;
    int getRankDim1(int coord) ;
    
    int * sizeLocalAllProcDim0();
    int * sizeLocalAllProcDim1();


    /* * * * * * * * * * * * * * * * * * * *
    [NEW!] OpenMP stuff 
    * * * * * * * * * * * * * * * * * * * */

    /*! \fn indexTransform(int* local_coord)
     \brief Transform from local coordinates (non-halo lattice points) to flattened (actual memory) index.

     \param local_coord local coordinate to transform (in 3D case: {i,j,k})

     \return int: flat index 
     */
    int indexTransform(int* local_coord);

   
   
    /*  - - - - functions 'for_each(...)' - - - -
        
        To be used when iterating over all of the lattice(s) sites

    */
    

    /*! \fn for_each(std::function<void(Site&)> operation)
     \brief Perform some operation on each lattice site.

     \details Can call this both inside an OpenMP parallel region, in which case it will use that environment, and outside. In the absence of an existing OpenMP parallel region, it will create and kill such an environment automatically iff the -fopenmp flag is given to the compiler, otherwise the function simply iterate through the lattice in the "old" way, using x.first() and so on, and execute the computation in the same way.

     \param operation lambda function of a site containing the computation to be performed in the loop
     */
    void for_each(std::function<void(Site&)> operation);


    /*! \fn for_each(std::function<void(Site*)> onsite_operation, Lattice ** other_lattices, const int num_lattices, string scheme)
     \brief Perform some operation on each lattice site.

     \details Can call this both inside an OpenMP parallel region, in which case it will use that environment, and outside. In the absence of an existing OpenMP parallel region, it will create and kill such an environment automatically iff the -fopenmp flag is given to the compiler, otherwise the function simply iterate through the lattice in the "old" way, using x.first() and so on, and execute the computation in the same way.

     \param onsite_operation lambda function of every site containing the computation to be performed in the loop
     \param other_lattices other lattices to consider
     \param num_lattices the TOTAL number of lattices to consider ( 1 + #{other lattices} )
     \param scheme key word for how to loop through the lattice ("controlled": loops through the grid with jumps of 2 and 2 in y and z dir)
     */
    void for_each(
        std::function<void(Site *)> onsite_operation, 
        Lattice ** other_lattices=NULL,
        const int num_lattices=1,
        string scheme="arbitrary");


    /*! \fn for_each(std::function<void(Site&, Site *)> onsite_operation, Lattice ** other_lattices, const int num_other_lattices, string scheme)
     \brief Perform some operation on each lattice site.

     \details Can call this both inside an OpenMP parallel region, in which case it will use that environment, and outside. In the absence of an existing OpenMP parallel region, it will create and kill such an environment automatically iff the -fopenmp flag is given to the compiler, otherwise the function simply iterate through the lattice in the "old" way, using x.first() and so on, and execute the computation in the same way.

     \param onsite_operation lambda function of every site containing the computation to be performed in the loop
     \param other_lattices other lattices to consider
     \param num_other_lattices the number of secondary lattices to consider ( #{other lattices} )
     \param scheme key word for how to loop through the lattice ("controlled": loops through the grid with jumps of 2 and 2 in y and z dir)
     */
    void for_each(
        std::function<void(Site&, Site *)> onsite_operation, 
        Lattice ** other_lattices=NULL,
        const int num_other_lattices=0,
        string scheme="arbitrary");




    /*! \fn for_each(std::function<void(Site&, Site *)> onsite_operation, Lattice ** fields, const int num_fields, string scheme)
     \brief Perform some operation on each lattice site.

     \details Can call this both inside an OpenMP parallel region, in which case it will use that environment, and outside. In the absence of an existing OpenMP parallel region, it will create and kill such an environment automatically iff the -fopenmp flag is given to the compiler, otherwise the function simply iterate through the lattice in the "old" way, using x.first() and so on, and execute the computation in the same way.

     \param onsite_operation lambda function of every site containing the computation to be performed in the loop
     \param fields fields from which to extract other lattice(s) [typically Field<Real>, but can be any class F that contains a member function Lattice& lattice()]
     \param num_fields number of fields from which to extract other lattice(s)
     \param scheme key word for how to loop through the lattice ("controlled": loops through the grid with jumps of 2 and 2 in y and z dir)
     */
    template<typename F>
    void for_each(
        std::function<void(Site&, Site *)> onsite_operation, 
        F ** fields=NULL,
        const int num_fields=0,
        string scheme="arbitrary");


    
    /*! \fn for_each(std::function<void(Site *)> onsite_operation, Lattice * other_lattice, string scheme)
     \brief Perform some operation on each lattice site on two lattices (this and other_lattice).

     \details Can call this both inside an OpenMP parallel region, in which case it will use that environment, and outside. In the absence of an existing OpenMP parallel region, it will create and kill such an environment automatically iff the -fopenmp flag is given to the compiler, otherwise the function simply iterate through the lattice in the "old" way, using x.first() and so on, and execute the computation in the same way.

     \param onsite_operation lambda function of every site (as pointer list) containing the computation to be performed in the loop
     \param other_lattice other lattice to consider
     \param scheme key word for how to loop through the lattice ("controlled": loops through the grid with jumps of 2 and 2 in y and z dir)
     */
    void for_each(
        std::function<void(Site *)> onsite_operation, 
        Lattice * other_lattice,
        string scheme="controlled");

    /*! \fn for_each(std::function<void(Site&, Site&)> onsite_operation, Lattice * other_lattice, string scheme)
     \brief Perform some operation on each lattice site on two lattices (this and other_lattice).

     \details Can call this both inside an OpenMP parallel region, in which case it will use that environment, and outside. In the absence of an existing OpenMP parallel region, it will create and kill such an environment automatically iff the -fopenmp flag is given to the compiler, otherwise the function simply iterate through the lattice in the "old" way, using x.first() and so on, and execute the computation in the same way.

     \param onsite_operation lambda function of both sites containing the computation to be performed in the loop
     \param other_lattice other lattice to consider
     \param scheme key word for how to loop through the lattice ("controlled": loops through the grid with jumps of 2 and 2 in y and z dir)
     */
    void for_each(
        std::function<void(Site&, Site&)> onsite_operation, 
        Lattice * other_lattice,
        string scheme="controlled");

    





private:
    int        status_;
    static int initialized;
    
    //Global variables==============
    int  dim_;              //ok//Number of dimensions
    int* size_;             //ok//Number of lattice sites in each direction
    long  sites_;            //ok//Number of sites in lattice
    long  sitesGross_;       //ok//Number of sites in lattice plus halos
    int  halo_;             //ok//Number of sites extra in each direction
    
    //Local variables===============
    int* sizeLocal_;       //ok//Number of local lattice sites in each direction
    int* sizeLocalAllProcDim0_;
    int* sizeLocalAllProcDim1_;
    long  sitesLocal_;      //ok//Number of local sites in lattice
    long  sitesLocalGross_; //ok//Number of local sites in lattice plus halo
    long* jump_;            //ok//Jumps needed to move up one in each direction
    
    
    long  siteFirst_;       //ok//Index of first local site in lattice
    
    long  siteLast_;        //ok//Index of last local site in lattice
    long  sitesSkip_;      //Number of global lattice sites before first local site (say in a file)
    long  sitesSkip2d_;
    long  coordSkip_[2];       //Number to add to coord[dim_-1] and coord[dim_-2] to get global value
    
    
    //save variable for fast save
    int arch_saved_;
    
};

#endif

!======================= 
!    type definition
!=======================

!-------------------------------
! element characteristics
!-------------------------------
TYPE elem_type

SEQUENCE

   ! normal vector - vnorm(edge index, cartesian coordinates)
   MYFLOAT vnorm(3,2)

   ! cell surface    
   MYFLOAT surface

   ! insphere radius    
   MYFLOAT inradius

   ! nb_degree of freedom 
   INTEGER(kind=4) nb_dof

   ! cpml index
   ! pointer to the table of memory variable pml 
   INTEGER(kind=4) pml_idx 

   ! type of face - face_type(face index)
   INTEGER(kind=4) face_type(3)

   ! nodes numbering - node_num(node index, face index)
   ! face index is the face number in face i  
   ! if node index = 1 --> local numbering in the cell i  
   ! if node index = 2 --> local numbering in the neighbour cell sharing this face  
   INTEGER(kind=2) node_perm_idx(2,3)

   ! neighbour cell number - neighbour(face index)
   INTEGER(kind=4) neighbour(3)    

   ! val idx
   ! pointer to the table of time solution val 
   INTEGER(kind=4) val_idx 

   ! dft idx
   ! pointer to the table of frequency solution dft 
   INTEGER(kind=4) dft_idx 

   ! pointer to the table of fault element
   INTEGER(kind=4) ft_idx

END TYPE elem_type

#ifdef DOUBLE_PRECISION
   INTEGER, PARAMETER :: ELEM_SIZE = 132
#else
   INTEGER, PARAMETER :: ELEM_SIZE = 100
#endif

!----------------------
! config parameters
!----------------------
TYPE config_type

SEQUENCE

   ! interpolation degree			  
   INTEGER inter

   ! =0 the mesh is read on file; =1: the mesh is constructed (see MESH)			  
   INTEGER imesh

   ! if imesh=1 then (iflag=1:mesh is uniform, iflag=2:mesh is an english flag) 
   INTEGER iflag

   ! level of output information
   INTEGER impre

   ! CFL condition
   MYFLOAT cfl
			  
   ! frequency of output information
   MYFLOAT ifre
 
   ! time max of seismogram
   MYFLOAT tmax

   ! snapshot parameters 
   INTEGER snapshot
   MYFLOAT min1, max1, min2, max2, snap_pixel, snap_time

   ! name of source signature file 
   CHARACTER(LEN=80) src_file

   ! dt used in the source file
   MYFLOAT src_dt

   ! sigma for smooth source 
   MYFLOAT src_sigma

   ! source tensor definition
   MYFLOAT mxz, myz
  
   ! length of PML layer in each direction
   MYFLOAT lpml_xmin, lpml_xmax, lpml_ymin, lpml_ymax

   ! maximum authorized elapse time 
   MYFLOAT tlimit 

   ! adapt interpolation order 
   MYFLOAT fmax, lim_P1, lim_P2
   INTEGER Pk_adapt

   ! CPML parameters
   INTEGER(kind=4) cpml_order   
   MYFLOAT cpml_npower, cpml_rcoef, cpml_freq, cpml_vmax, cpml_add_damp

   ! switch to activate DFT 
   LOGICAL dft_switch

   ! output stress component
   INTEGER output_stress

   ! test case
   INTEGER icas

   ! cell average for extraction of results
   INTEGER cell_aver	

   ! number of time step between 2 seismograms recording 
   INTEGER seismo_buf

   ! number of frequencies to be extracted
   INTEGER nb_freq  

   ! source position
   INTEGER src_loc

   ! nb of source computed in parallel
   INTEGER nb_src_paral

   ! source type (0 = point source, 1 = finite source)
   INTEGER src_type 

   ! fault simulation (0 = kinematic, 1 = dynamic)
   INTEGER fault_sim

   ! name of fault vectors file 
   CHARACTER(LEN=80) fault_vectors_file

   ! name of fault file 
   CHARACTER(LEN=80) fault_file

   ! Normal stress evolution (0 = fixed, 1 = variable)
   INTEGER ns_evol

   ! fault snapshot parameters 
   INTEGER fault_snapshot
   MYFLOAT fault_snap_time

! ------------- Inversion ----------------

   ! Mode (0 = Single FWI gradient, 1 = FWI)
   INTEGER mode

   ! Inversion algorithm (0 = Steepest descent, 1 = Non-linear CG, 2 = LBFGS)
   INTEGER inv_alg

   ! Convergence criteria
   MYFLOAT conv 

   ! Maximum number of iterations
   INTEGER max_iter

   ! Precondition choice
   INTEGER precond

   ! Debug option
   INTEGER debug_inv
			  
END TYPE config_type

!------------------------------
! xtab - communication table
!------------------------------
TYPE xtab_type

SEQUENCE

   ! process id
   INTEGER id_sub

   ! number of faces shared between subdomains
   INTEGER nb_face

   ! nb values to send 
   INTEGER nb_send

   ! nb values to receive
   INTEGER nb_rec

END TYPE xtab_type

INTEGER, PARAMETER :: XTAB_SIZE = 16

!----------------------
! subdomain data
!----------------------
TYPE sub_type

SEQUENCE

   ! extremum coordinates of subdomain 
   MYFLOAT xmin, ymin, xmax, ymax  

   ! number of internal elem 
   INTEGER nti

   ! number of external elem (physical boundary)  
   INTEGER nte

   ! number of artificial elem (intra subdomain boundary)
   INTEGER ntx 

   ! number of fault faces
   INTEGER nfault 

   ! number of elem in CPML layers
   INTEGER ntpml  

   ! number of cell per interpolation order
   INTEGER ntp2  

   ! number of nodes 
   INTEGER ns

   ! number of neighbour subdomain 
   INTEGER nsub

   ! number of elements over the fault
   ! (without the artificial elements)
   INTEGER nfault_elem

   ! number of border fault elements 
   INTEGER nfault_ext
  
END TYPE sub_type

!----------------------
! global domain data
!----------------------
TYPE glob_type

SEQUENCE

   ! extremum
   MYFLOAT xmin, ymin, xmax, ymax

   ! total surface
   MYFLOAT tot_sur

   ! number of grid points 
   INTEGER nxg, nyg
     
   ! number of nodes
   INTEGER nsg

   ! number of elements
   INTEGER ntg

   ! number of external faces
   INTEGER nfrg

   ! number of fault faces
   INTEGER nfault                 

   ! number of elem in CPML layers
   INTEGER ntpml  

   ! number of cell per interpolation order
   INTEGER ntp0, ntp1, ntp2  

   ! number of nodes in the fault
   INTEGER ns_fault

END TYPE glob_type

!----------------------
! mesh parameters
!----------------------
TYPE mesh_param_type

SEQUENCE

   ! size of medium (m) 
   MYFLOAT rlx, rly

   ! origin 
   MYFLOAT xor, yor

   ! number of points in the 2 direction of the numerical model 
   INTEGER nx, ny

   ! number of process per direction
   INTEGER ndimx, ndimy  

   ! boundary conditions at the 4 edges of the model
   INTEGER ifx1, ifx2, ify1, ify2

   ! flag to store the mesh on file
   INTEGER iecr

END TYPE mesh_param_type

!---------------------------
! physical files parameters
!---------------------------
TYPE physical_param_type

SEQUENCE

   ! space sampling used in the cartesian grid
   MYFLOAT dx

   ! Number of points in the physical model in the 2 axis
   INTEGER nx, ny

   ! name of vp, vs and rho files
   CHARACTER (LEN=80) vp_file, vs_file, rho_file    

END TYPE physical_param_type

!------------------------
! acquisition parameters
!------------------------
TYPE acqui_type

SEQUENCE

   ! amplitude of source (cell source index, node index)
   MYFLOAT amp_loc(MAX_SRC, MAX_DOF) 

   ! source coordinates (source index, cartesian coordinates)
   MYFLOAT shot(MAX_SHOT, 2)

   ! Number of source cells (local)
   INTEGER nb_src_loc 

   ! Indice of the source cells (local)
   INTEGER src_loc(MAX_SRC) 

   ! Number of shot
   INTEGER n_shot

   ! Number of receiver per shot - n_receiver(shot index)
   INTEGER n_receiver(MAX_SHOT)

   ! Number of receiver located in the subdomain 
   INTEGER n_rec_loc 

   ! Receiver coordinates per source (source indice, receiver indice, cartesian coordinates)
   MYFLOAT, DIMENSION(:,:,:), POINTER :: receiver

   ! Global receiver number and local cell number   
   ! rec_loc(receiver number in subdomain, 1=global receiver number, 2=local cell number)
   INTEGER, DIMENSION(:,:), POINTER :: rec_loc

   ! Coefficient to apply at each node for a receiver cell
   ! rec_coef(receiver number in domain, node index)  
   MYFLOAT, DIMENSION(:,:), POINTER :: rec_coef  

   ! Number of receivers (FOR SNAPSHOT)
   INTEGER n_receiver_snap 

   ! Number of receiver located in the subdomain (FOR SNAPSHOT) 
   INTEGER n_snap_loc      

   ! Receiver coordinates (receiver index, cartesian coordinates) (FOR SNAPSHOT)
   MYFLOAT, DIMENSION(:,:), POINTER :: receiver_snap(:,:)

   ! Global receiver number and local cell cell number (FOR SNAPSHOT)   
   ! snap_loc(receiver number in subdomain, 1=global receiver number, 2=local cell number)
   INTEGER, DIMENSION(:,:), POINTER :: snap_loc

   ! Coefficient to apply at each node for a receiver cell (FOR SNAPSHOT)
   ! snap_coef(receiver number in domain, node index)  
   MYFLOAT, DIMENSION(:,:), POINTER :: snap_coef   

!---------- For the adjoint sources ----------------

   ! Normalisation of the excitation
   ! rec_norm(receiver index in subdomain)
   MYFLOAT, DIMENSION(:), POINTER :: rec_norm 

   ! *** For smooth adjoint sources ***
   
   ! Number of smoothed adjoint sources located in the subdomain per receiver
   ! 1st index : Global receiver index
   INTEGER, DIMENSION(:), POINTER :: n_sas_loc

   ! Local cell number of the adjoint sources located in the subdomain per receiver
   ! 1st index : Global receiver index
   ! 2nd index : Number of smoothed adjoint sources located in the subdomain per receiver (n_sas_loc(irec)) 
   INTEGER, DIMENSION(:,:), POINTER :: sas_loc

   ! Coefficient to apply at each node for every adjoint sources (it's included the normalisation factor)
   ! 1st index : Global receiver index
   ! 2nd index : Number of smoothed adjoint sources located in the subdomain per receiver (n_sas_loc(irec)) 
   ! 3rd index : DOF per element
   MYFLOAT, DIMENSION(:,:,:), POINTER :: sas_amp

END TYPE acqui_type

!--------------
! PML tables
!--------------
TYPE pml_type

SEQUENCE

   ! 1st coefficient - coef1(x, y)
   MYFLOAT coef1(2)

   ! 2nd coefficient - coef2(x, y)
   MYFLOAT coef2(2)

   ! memory variable for velocity components
   MYFLOAT mem_v(2, MAX_DOF)

   ! memory variable for stress components
   MYFLOAT mem_sigmax(1, MAX_DOF)
   MYFLOAT mem_sigmay(1, MAX_DOF)

END TYPE pml_type

#ifdef DOUBLE_PRECISION
   INTEGER, PARAMETER :: PML_SIZE = 176
#else
   INTEGER, PARAMETER :: PML_SIZE = 88
#endif

!---------------------------------
! acquisition fault parameters
!---------------------------------
TYPE acqui_fault_type

SEQUENCE

   ! Number of receivers over the fault (global)
   INTEGER nrec 

   ! Number of receiver located in the subdomain 
   INTEGER nrec_loc 

   ! Coordinates of the receiver over the fault (global)
   MYFLOAT, DIMENSION(:,:), POINTER :: rec_coor 

   ! Global receiver number and local cell cell number   
   ! rec_loc(receiver number in subdomain, 1=global receiver number, 2=local cell number)
   INTEGER, DIMENSION(:,:), POINTER :: rec_loc

   ! Coefficient to apply at each node for a receiver cell
   ! rec_coef(receiver number in domain, node index)  
   MYFLOAT, DIMENSION(:,:), POINTER :: rec_coef

END TYPE acqui_fault_type

!---------------------------------
! element over the fault
!---------------------------------
TYPE elem_fault_type

SEQUENCE

   ! Interior element (0 - Exterior, 1 - Interior)
   INTEGER interior

   ! Nucleation element (0 - Outside, 1 - Nucleation)
   INTEGER nucleation

   ! Side element (0 - Minus side, 1 - Plus side)
   INTEGER side

   ! Local index of the element
   INTEGER elem_idx

   ! Local index of the face fault of elem_idx
   INTEGER face_fault

   ! Physical id fault segment
   INTEGER physical_id

   ! Slip magnitude
   MYFLOAT U(MAX_DOF)

   ! Slip rate
   MYFLOAT V(MAX_DOF)

   ! Slip rate history (For kinematic simulations)
   MYFLOAT, DIMENSION(:), POINTER :: V_history

   ! Shear stress
   MYFLOAT ss(MAX_DOF)

   ! Normal stress
   MYFLOAT ns(MAX_DOF)

   ! Rupture time
   MYFLOAT rt(MAX_DOF)

   ! Neighbour flux 
   MYFLOAT neighbour_flux(MAX_DOF*1)

   ! Neighbour density
   MYFLOAT neighbour_rho

   ! Initial normal stress
   MYFLOAT init_normal(MAX_DOF)

   ! Initial shear stress along strike
   MYFLOAT init_shear(MAX_DOF)

   ! Frictional parameters (Slip-weakening friction)
   MYFLOAT fault_friction(5,MAX_DOF)

END TYPE elem_fault_type

#ifdef DOUBLE_PRECISION
   INTEGER, PARAMETER :: ELEM_FAULT_SIZE = 856
#else
   INTEGER, PARAMETER :: ELEM_FAULT_SIZE = 440
#endif

!---------------------
! Full wave inversion 
!---------------------
TYPE fwi_type

SEQUENCE

   ! Adjoint sources
   MYFLOAT, DIMENSION(:,:), POINTER :: adj_sources

   ! Model ( slip-rate(x,t) of the forward problem)
   MYFLOAT, DIMENSION(:), POINTER :: model

   ! Gradient ( Traction(x,t) of the adjoint problem)
   MYFLOAT, DIMENSION(:), POINTER :: gradient

   ! Observations 
   MYFLOAT, DIMENSION(:,:), POINTER :: obs

   ! Cost function evaluation
   MYFLOAT fcost

   ! *** For smoothed adjoint sources ***
   
   ! Adjoint sources
   MYFLOAT, DIMENSION(:,:), POINTER :: sadj_sources

   ! Observations 
   MYFLOAT, DIMENSION(:,:), POINTER :: sobs
  
END TYPE fwi_type

!---------------------
! Butterworth filter
!---------------------
TYPE butterworth_filter_type

SEQUENCE

   ! Order of the butterworth filter
   INTEGER order

   ! Cutoff frequency
   MYFLOAT fc

   ! Filter coefficients
   MYFLOAT C(5,10)

   ! Coefficients defining the filter memory
   MYFLOAT D(2,10)

   ! Group delay in seconds
   MYFLOAT Tg

   ! Number of required 2nd order sections
   INTEGER NSections
  
END TYPE butterworth_filter_type

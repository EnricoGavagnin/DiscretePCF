Pair Correlation Funcitons for discrete domains

PCF=PCF_discrete(type,M,opt) returns the vector with the values of the
discrete pair correlation function (PCF) with metric specified by 'type'
and occupancy matrix M.

Ref: 
Gavagnin E., Owen J.P. and Yates C.A. "Pair correlation functions  
for identifying spatial correlation in discrete domains"  
doi: arXiv:1804.03452

-------------------------------------------------------------------------- 
INPUT:

type: string which specifies the type of PCF to use:
       PCF_discrete('taxicab',...): Square Taxicab PCF or Cube Taxicab PCF
       PCF_discrete('uniform',...): Square Uniform PCF or Cube Uniform PCF
       PCF_discrete('triangular',...): Triangular PCF
       PCF_discrete('hexagonal',...): Hexagonal PCF
       PCF_discrete('irregular',...): General PCF
      (see the paper in Ref. for detailed definitions).

   M: If type='taxicab', 'unifom', 'triangular' and 'hexagonal' M is 
      ccupancy matrix of the discerte domain. M(i,j)=1 if the site (i,j) 
      is occupied and M(i,j)=0 if site (i,j) is empty. 
      (For 'triangular' and 'hexagonal', sites are numbered as in 
      Figure 10 of the paper in Ref.)       

      If type='irregular', M is a vector with labels corresponding to
      occupied sites of the irregular domain.

 opt: If type='taxicab', 'unifom', 'triangular' or 'hexagonal'
      opt spcifies the type of boundary conditions in both the
      boundaries:

      PCF_discrete(...,'periodic'): periodic boundary conditions
                                    in both directions  
      PCF_discrete(...,'nonperiodic'): non periodic boundary conditions
                                       in both directions  
      
      If If type='irregular', opt reads the adjacency matrix of the
      irregular domain. 

OUTPUT:

PCF: vecotr with the values of the pair correlation funciton

-------------------------------------------------------------------------- 
EXAMPLES:

  PCF=PCF_discrete('taxicab',M,'periodic') returns the Square Taxicab or 
  Cube Taxicab PCF, depending on the dimension of M, and with periodic
  boundary conditions in both directions. 

  PCF=PCF_discrete('irregular',M,A) returns the General PCF with adjacency
  matrix A and occupied sites M.  

-------------------------------------------------------------------------- 

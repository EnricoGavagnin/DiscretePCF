%%--- help for PCF_discrete ---
%
% PCF=PCF_discrete(type,M,opt) returns the vector with the values of the
% discrete pair correlation function (PCF) with metric specified by 'type'
% and occupancy matrix M.
%
% Ref: 
% Gavagnin E., Owen J.P. and Yates C.A. "Pair correlation functions  
% for identifying spatial correlation in discrete domains"  
% doi: arXiv:1804.03452
%
% Created 07/02/2018
% By Enrico Gavagnin
% email: e.gavagnin@bath.ac.uk
%
%-------------------------------------------------------------------------- 
%INPUT:
%
% type: string which specifies the type of PCF to use:
%        PCF_discrete('taxicab',...): Square Taxicab PCF or Cube Taxicab PCF
%        PCF_discrete('uniform',...): Square Uniform PCF or Cube Uniform PCF
%        PCF_discrete('triangular',...): Triangular PCF
%        PCF_discrete('hexagonal',...): Hexagonal PCF
%        PCF_discrete('irregular',...): General PCF
%       (see the paper in Ref. for detailed definitions).
%
%    M: If type='taxicab', 'unifom', 'triangular' and 'hexagonal' M is 
%       ccupancy matrix of the discerte domain. M(i,j)=1 if the site (i,j) 
%       is occupied and M(i,j)=0 if site (i,j) is empty. 
%       (For 'triangular' and 'hexagonal', sites are numbered as in 
%       Figure 10 of the paper in Ref.)       
%
%       If type='irregular', M is a vector with labels corresponding to
%       occupied sites of the irregular domain.
% 
%  opt: If type='taxicab', 'unifom', 'triangular' or 'hexagonal'
%       opt spcifies the type of boundary conditions in both the
%       boundaries:
%
%       PCF_discrete(...,'periodic'): periodic boundary conditions
%                                     in both directions  
%       PCF_discrete(...,'nonperiodic'): non periodic boundary conditions
%                                        in both directions  
%       
%       If If type='irregular', opt reads the adjacency matrix of the
%       irregular domain. 
%
% OUTPUT:
% 
% PCF: vecotr with the values of the pair correlation funciton
%
%-------------------------------------------------------------------------- 
% EXAMPLES:
%
%   PCF=PCF_discrete('taxicab',M,'periodic') returns the Square Taxicab or 
%   Cube Taxicab PCF, depending on the dimension of M, and with periodic
%   boundary conditions in both directions. 
%
%   PCF=PCF_discrete('irregular',M,A) returns the General PCF with adjacency
%   matrix A and occupied sites M.  
%
%-------------------------------------------------------------------------- 



function [PCF] = PCF_discrete(type,M,opt)


%% Square Taxicab
if strcmp(type,'taxicab') && size(size(M),2)==2
    
    %Store the dimensions of the domain
    [Ly,Lx]=size(M);
    
    %Store the number of occupied sites
    N=nnz(M);
    
    %Store a vector with the coordinates of the occupied sites
    % C=¦y¦x¦
    C=zeros(N,2);
    [C(:,1),C(:,2)]=find(M);
    
    
    if strcmp(opt,'periodic')
        %% Square taxicab periodic BC
        
        %Define the maximum range of the PCF
        r_max=floor(min(Lx,Ly)/2);
        
        %Preallocation PCF
        PCF=zeros(1,r_max);
        
        % Loop through the occupied lattice sites
        for i=1:N-1
            for j=i+1:N
                
                %Compute the distance between sites i and j with periodic
                %taxicab metric
                r=min(abs(C(i,1)-C(j,1)),abs(-abs(C(i,1)-C(j,1))+Ly))+...
                    min(abs(C(i,2)-C(j,2)),abs(-abs(C(i,2)-C(j,2))+Lx));
                
                
                if r<=r_max
                    %Store the contribution for the PCF
                    PCF(r)=PCF(r)+1;
                end
            end
        end
        
        %Normalisation factor
        f_norm=(2*Lx*Ly)*(1:r_max);
        
        %Normalisation
        PCF=Lx*Ly*(Lx*Ly-1)*PCF./(N*(N-1)*f_norm);
        
    elseif strcmp(opt,'nonperiodic')
        %% Square taxicab non periodic BC
        
        %Define the maximum range of the PCF
        r_max=min(Lx,Ly)-1;
        
        %Preallocation PCF
        PCF=zeros(1,r_max);
        
        % Loop through the lattice sites
        for i=1:N-1
            for j=i+1:N
                
                %Compute the distance between sites i and j using non
                %periodic taxicab metric
                r=abs(C(i,1)-C(j,1))+abs(C(i,2)-C(j,2));
                
                if r<=r_max
                    %Store the contribution for the PCF
                    PCF(r)=PCF(r)+1;
                end
            end
        end
        
        %Normalisation factor
        f_norm=(2*Lx*Ly)*(1:r_max)-((1:r_max).^2*(Lx+Ly)-1/3*((1:r_max).^3-(1:r_max)));
        
        %Normalisation
        PCF=Lx*Ly*(Lx*Ly-1)*PCF./(N*(N-1)*f_norm);
        
        
    else
        error('Wrong boundary conditions mode')
    end
    
    %% Square uniform
elseif strcmp(type,'uniform') && size(size(M),2)==2
    
    %Store the dimensions of the domain
    [Ly,Lx]=size(M);
    
    %Store the number of occupied sites
    N=nnz(M);
    
    %Store a vector with the coordinates of the occupied sites
    % C=¦y¦x¦
    C=zeros(N,2);
    [C(:,1),C(:,2)]=find(M);
    
    
    if strcmp(opt,'periodic')
        %% Square uniform periodic BC
        
        %Define the maximum range of the PCF
        r_max=ceil(min(Lx,Ly)/2)-1;
        
        %Preallocation PCF
        PCF=zeros(1,r_max);
        
        % Loop through the occupied lattice sites
        for i=1:N-1
            for j=i+1:N
                
                %Compute the distance between sites i and j with periodic
                %uniform metric
                r=max(min(abs(C(i,1)-C(j,1)),abs(-abs(C(i,1)-...
                    C(j,1))+Ly)),min(abs(C(i,2)-C(j,2)),...
                    abs(-abs(C(i,2)-C(j,2))+Lx)));
                
                
                if r<=r_max
                    %Store the contribution for the PCF
                    PCF(r)=PCF(r)+1;
                end
            end
        end
        
        %Normalisation factor
        f_norm=(4*Lx*Ly)*(1:r_max);
        
        %Normalisation
        PCF=Lx*Ly*(Lx*Ly-1)*PCF./(N*(N-1)*f_norm);
        
    elseif strcmp(opt,'nonperiodic')
        %% Square uniform non periodic BC
        
        %Define the maximum range of the PCF
        r_max=min(Lx,Ly)-1;
        
        %Preallocation PCF
        PCF=zeros(1,r_max);
        
        % Loop through the lattice sites
        for i=1:N-1
            for j=i+1:N
                
                %Compute the distance between sites i and j using non
                %periodic uniform metric
                r=max(abs(C(i,1)-C(j,1)),abs(C(i,2)-C(j,2)));
                
                if r<=r_max
                    %Store the contribution for the PCF
                    PCF(r)=PCF(r)+1;
                end
            end
        end
        
        %Normalisation factor
        f_norm=(4*Lx*Ly)*(1:r_max)-((1:r_max).^2*(3*Lx+3*Ly)-2*(1:r_max).^3);
        
        %Normalisation
        PCF=Lx*Ly*(Lx*Ly-1)*PCF./(N*(N-1)*f_norm);
        
        
    else
        error('Wrong boundary conditions mode')
    end
    
    %% Cube Taxicab
elseif strcmp(type,'taxicab') && size(size(M),2)==3
    
    %Store the dimensions of the domain
    [Ly,Lx,Lz]=size(M);
    
    %Store the number of occupied sites
    N=nnz(M);
    
    %Store a vector with the coordinates of the occupied sites
    % C=¦y¦x¦
    C=zeros(N,3);
    %Store the dimensions of the domain
    ind = find(M);
    [C(:,1), C(:,2), C(:,3)] = ind2sub(size(M), ind);
    
    
    if strcmp(opt,'periodic')
        %% Cube taxicab periodic BC

       
        %Define the maximum range of the PCF
        r_max=floor(min(min(Lx,Ly),Lz)/2);
        
        %Preallocation PCF
        PCF=zeros(1,r_max);
        
        % Loop through the occupied lattice sites
        for i=1:N-1
            for j=i+1:N
                
                %Compute the distance between sites i and j with periodic
                %taxicab metric
                r=min(abs(C(i,1)-C(j,1)),abs(-abs(C(i,1)-C(j,1))+Ly))+...
                    min(abs(C(i,2)-C(j,2)),abs(-abs(C(i,2)-C(j,2))+Lx))+...
                    min(abs(C(i,3)-C(j,3)),abs(-abs(C(i,3)-C(j,3))+Lz));
                
                
                if r<=r_max
                    %Store the contribution for the PCF
                    PCF(r)=PCF(r)+1;
                end
            end
        end
        
        %Normalisation factor
        f_norm=(2*(1:r_max).^2+ones(1,r_max))*Lx*Ly*Lz;
        
        %Normalisation
        PCF=Lx*Ly*Lz*(Lx*Ly*Lz-1)*PCF./(N*(N-1)*f_norm);
        
    elseif strcmp(opt,'nonperiodic')
        %% Cube taxicab non periodic BC
        
        %Define the maximum range of the PCF
        r_max=min(Lx,min(Ly,Lz))-1;
        
        %Preallocation PCF
        PCF=zeros(1,r_max);
        
        % Loop through the lattice sites
        for i=1:N-1
            for j=i+1:N
                
                %Compute the distance between sites i and j using non
                %periodic taxicab metric
                r=abs(C(i,1)-C(j,1))+abs(C(i,2)-C(j,2))+abs(C(i,3)-C(j,3));
                
                if r<=r_max
                    %Store the contribution for the PCF
                    PCF(r)=PCF(r)+1;
                end
            end
        end
        
        %Normalisation factor
        f_norm=(2*(1:r_max).^2+ones(1,r_max))*Lx*Ly*Lz-1/3*(2*(1:r_max).^3 ...
            +(1:r_max))*(Lx*Ly+Ly*Lz+Lx*Lz)+(1:r_max).^2 .* ((1:r_max).^2 ...
            -ones(1,r_max))/6*(Lx+Ly+Lz)-(1:r_max).^5/30+(1:r_max).^3/6 ...
            -2/15*(1:r_max);
        
        %Normalisation
        PCF=Lx*Ly*Lz*(Lx*Ly*Lz-1)*PCF./(N*(N-1)*f_norm);
        
        
    else
        error('Wrong boundary conditions mode')
    end
    
    %% Cube uniform
elseif strcmp(type,'uniform') && size(size(M),2)==3
    
    %Store the dimensions of the domain
    [Ly,Lx,Lz]=size(M);
    
    %Store the number of occupied sites
    N=nnz(M);
    
    %Store a vector with the coordinates of the occupied sites
    % C=¦y¦x¦z¦
     %Store the dimensions of the domain
    ind = find(M);
    [C(:,1), C(:,2), C(:,3)] = ind2sub(size(M), ind);
    
    
    if strcmp(opt,'periodic')
        %% Square uniform periodic BC
        
        %Define the maximum range of the PCF
        r_max=ceil(min(min(Lx,Ly),Lz)/2)-1;
        
        %Preallocation PCF
        PCF=zeros(1,r_max);
        
        % Loop through the occupied lattice sites
        for i=1:N-1
            for j=i+1:N
                
                %Compute the distance between sites i and j with periodic
                %uniform metric
                r=max(min(abs(C(i,1)-C(j,1)),abs(-abs(C(i,1)-...
                    C(j,1))+Ly)),max(min(abs(C(i,2)-C(j,2)),...
                    abs(-abs(C(i,2)-C(j,2))+Lx)),min(abs(C(i,3)...
                    -C(j,3)),abs(-abs(C(i,3)-C(j,3))+Lz))));
                
                
                if r<=r_max
                    %Store the contribution for the PCF
                    PCF(r)=PCF(r)+1;
                end
            end
        end
        
        %Normalisation factor
        f_norm=(12*(1:r_max).^2+ones(1,r_max))*Lx*Ly*Lz;
        
        %Normalisation
        PCF=Lx*Ly*Lz*(Lx*Lz*Ly-1)*PCF./(N*(N-1)*f_norm);
        
    elseif strcmp(opt,'nonperiodic')
        %% Square uniform non periodic BC
        
        %Define the maximum range of the PCF
        r_max=min(Lx,min(Ly,Lz))-1;
        
        %Preallocation PCF
        PCF=zeros(1,r_max);
        
        % Loop through the lattice sites
        for i=1:N-1
            for j=i+1:N
                
                %Compute the distance between sites i and j using non
                %periodic uniform metric
                r=max(abs(C(i,1)-C(j,1)),max(abs(C(i,2)-C(j,2)),abs(C(i,3)-C(j,3))));
                
                if r<=r_max
                    %Store the contribution for the PCF
                    PCF(r)=PCF(r)+1;
                end
            end
        end
        
        %Normalisation factor
        f_norm=(12*(1:r_max).^2+ones(1,r_max))*Lx*Ly*Lz-(1:r_max).* ...
            (8*(1:r_max).^2+1)*(Lx*Ly+Ly*Lz+Lx*Lz)+(1:r_max).^2 .* ...
            (5*(1:r_max).^2+ones(1,r_max))*(Lx+Ly+Lz)-(1:r_max).^3 ...
            .* (3*(1:r_max).^2+ones(1,r_max));
        
        %Normalisation
        PCF=Lx*Ly*Lz*(Lx*Ly*Lz-1)*PCF./(N*(N-1)*f_norm);
        
        
    else
        error('Wrong boundary conditions mode')
    end
    
    
    %% Triangular
elseif strcmp(type,'triangular')
    
    %Store the dimensions of the domain
    [Ly,Lx]=size(M);
    
    %Check if the domain sizes are eligible for periodic BC
    if mod(Ly,2)~=0
        error('The y-dimension of the lattice has to be even in order to use periodic BC on a triangular tesselation. Type "help PCF_calculator" for more information');
    end
    
    %Store the number of occupied sites
    N=nnz(M);
    
    %Store a vector with the coordinates of the occupied sites
    % C=¦y¦x¦
    C=zeros(N,2);
    [C(:,1),C(:,2)]=find(M);
    
    
    if strcmp(opt,'periodic')
        %% Triangular with periodic BC 
        
        %Define the maximum range of the PCF
        r_max=floor(min(Lx,Ly)/2)-1;
        
        %Preallocation PCF
        PCF=zeros(1,r_max);
        
        % Loop through the occupied lattice sites
        for i=1:N-1
            for j=i+1:N
                
                %Compute the distance between sites i and j with periodic
                %taxicab metric
                
                dy=min(abs(C(i,1)-C(j,1)),abs(-abs(C(i,1)-...
                    C(j,1))+Ly));
                dx=min(abs(C(i,2)-C(j,2)),abs(-abs(C(i,2)-...
                    C(j,2))+Lx));
                M=ceil(dy/2);
                
                if mod(sum(C(i,:)),2)==0 && mod(sum(C(j,:)),2)==1
                    r=dx+dy+2*max(0,floor((dy-dx)/2));
                else
                    r=dx+dy+2*max(0,ceil((dy-dx)/2));
                end
                
                if r<=r_max
                    
                    %Store the contribution for the PCF
                    PCF(r)=PCF(r)+1;
                end
            end
        end
        
        %Normalisation factor
        f_norm=(3/2*Lx*Ly)*(1:r_max);
        
        %Normalisation
        PCF=Lx*Ly*(Lx*Ly-1)*PCF./(N*(N-1)*f_norm);
        
    elseif strcmp(opt,'nonperiodic')
        %% Triangular with non periodic BC 
        
        %Define the maximum range of the PCF
        r_max=floor(min(Lx,Ly))-1;
        
        %Preallocation PCF
        PCF=zeros(1,r_max);
        
        % Loop through the occupied lattice sites
        for i=1:N-1
            for j=i+1:N
                
                %Compute the distance between sites i and j with non periodic
                %taxicab metric
                
                dy=abs(C(i,1)-C(j,1));
                dx=abs(C(i,2)-C(j,2));
                
                M=ceil(dy/2);
                
                if mod(sum(C(i,:)),2)==0 && mod(sum(C(j,:)),2)==1
                    r=dx+dy+2*max(0,floor((dy-dx)/2));
                else
                    r=dx+dy+2*max(0,ceil((dy-dx)/2));
                end
                
                if r<=r_max
                    
                    %Store the contribution for the PCF
                    PCF(r)=PCF(r)+1;
                end
            end
        end
        
        
        %Normalisation factor
        f_norm=zeros(1,r_max);
        f_norm(1)=3/2*Lx*Ly-Lx/2-Ly;
        f_norm(2)=3*Lx*Ly-2*Lx-4*Ly+2;
        
        for m=3:r_max
            k_three=floor((m-3)/4);
            k_six=floor((m-6)/4);
            k_seven=floor((m-7)/4);
            f_norm(m)=3/2*m*Lx*Ly-m^2/2*Lx+(2*k_six*(k_six-2*k_three+1)...
                +k_three*(m-6)-m^2+m-2)*Ly+1/3*(m-1)*(m^2-2*m+6)...
                -1/3*(k_seven+1)*(20*k_seven^2+37*k_seven+12)...
                -(m-7-4*k_seven)*(k_seven+1)*(m+k_seven-2);
        end
        
        %Normalisation
        PCF=Lx*Ly*(Lx*Ly-1)*PCF./(N*(N-1)*f_norm);
        
        
    else
        error('Wrong boundary conditions mode')
    end
    
  %% Hexagonal
elseif strcmp(type,'hexagonal')
    
    %Store the dimensions of the domain
    [Ly,Lx]=size(M);
    
    %Check if the domain sizes are eligible for periodic BC
    if mod(Ly,2)~=0
        error('The y-dimension of the lattice has to be even in order to use periodic BC on a hexagonal tesselation. Type "help PCF_calculator" for more information');
    end
    
    %Store the number of occupied sites
    N=nnz(M);
    
    %Store a vector with the coordinates of the occupied sites
    % C=¦y¦x¦
    C=zeros(N,2);
    [C(:,1),C(:,2)]=find(M);
    
    
    if strcmp(opt,'periodic')
        %% Triangular with periodic BC 
        
        %Define the maximum range of the PCF
        r_max=floor(min(Lx,Ly)/2)-1;
        
        %Preallocation PCF
        PCF=zeros(1,r_max);
        
        % Loop through the occupied lattice sites
        for i=1:N-1
            for j=i+1:N
                
                %Compute the distance between sites i and j with periodic
                %taxicab metric
                
                dy=min(abs(C(i,1)-C(j,1)),abs(-abs(C(i,1)-...
                    C(j,1))+Ly));
                dx=min(abs(C(i,2)-C(j,2)),abs(-abs(C(i,2)-...
                    C(j,2))+Lx));
                
                
                if mod(C(i,1)+C(i,2),2)==0
                   r=max(dx,dy)+max(0,min(dy-ceil((dx-1)/2),floor((dx+1)/2)));
                  
                else
                  
                   r=max(dx,dy)+max(0,min(dy-ceil(dx/2),floor(dx/2)));
                  
                end
                  
                if r<=r_max
                    
                    %Store the contribution for the PCF
                    PCF(r)=PCF(r)+1;
                end
            end
        end
        
        %Normalisation factor
        f_norm=(3*Lx*Ly)*(1:r_max);
        
        %Normalisation
        PCF=Lx*Ly*(Lx*Ly-1)*PCF./(N*(N-1)*f_norm);
        
    elseif strcmp(opt,'nonperiodic')
        %% Hexagonal with non periodic BC 
        
        %Define the maximum range of the PCF
        r_max=floor(min(Lx,Ly))-1;
        
        %Preallocation PCF
        PCF=zeros(1,r_max);
        
        % Loop through the occupied lattice sites
        for i=1:N-1
            for j=i+1:N
                
                %Compute the distance between sites i and j with non periodic
                %taxicab metric              
                dy=abs(C(i,1)-C(j,1));
                dx=abs(C(i,2)-C(j,2));
              
                if mod(C(i,1)+C(i,2),2)==0
                   r=max(dx,dy)+max(0,min(dy-ceil((dx-1)/2),floor((dx+1)/2)));
                else
                   r=max(dx,dy)+max(0,min(dy-ceil(dx/2),floor(dx/2)));
                end
                  
                if r<=r_max
                    %Store the contribution for the PCF
                    PCF(r)=PCF(r)+1;
                end
               
            end
        end
        
        
        %Normalisation factor
        k=mod(1:r_max,2);
        
        f_norm=3*Lx*Ly*(1:r_max)-1/4*(7*(1:r_max).^2+k)*Lx ...
            -2*(1:r_max).^2*Ly++11/12*(1:r_max).^3 ...
            -(2-3*k).*(1:r_max)/12;
        
        %Normalisation
        PCF=Lx*Ly*(Lx*Ly-1)*PCF./(N*(N-1)*f_norm);
        
        
    else
        error('Wrong boundary conditions mode')
    end   
    
    
    %% Irregular domain (General PCF)
elseif strcmp(type,'irregular')
    if min(size(M))>1
        error('L should be a vecotr. Type "help PCF_calculator" for more details.')
    end
    
    %Define the maximum range of the PCF
    r_max=size(opt,1);
    
    %Preallocation PCF
    PCF=zeros(1,r_max);
    
    %Computation of the distance matrix and normalisation

    %Initialisation
    Dist=eye(size(opt));
    prevA=opt;    
    r=2;
    
    %Number of occupied sites
    N=length(M);
 
    %Loop that compute the powers of the A matrix by induction and
    %automatically update the matrix of the distances
    
    while size(opt,1)^2-nnz(Dist)>0
        
        %Compute the matrix indicating of the element which become non
        %zeros for the first time at the r-th power of the A
        Reached=double(prevA>0).*double(Dist==0);
  
        %Update the matrix that store the pairwise distances
        Dist=Dist+(r-1)*Reached;
        
        %Update the PCF
        PCF(r-1)=size(opt,1)*(size(opt,1)-1)*nnz(Reached(M,M))/(N*(N-1)*nnz(Reached));
        
        %Update the power of A
        prevA=prevA*opt;
        
        %Update the index r
        r=r+1;
    end
    
    %Cut the tail of the PCF
    PCF=PCF(1:r-3);
    
else
    error('Wrong metric mode')
end


end

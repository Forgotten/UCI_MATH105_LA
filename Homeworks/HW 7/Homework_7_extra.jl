# Functions to create FD matrices discretizing the Helmholtz equation in
# 2D with absorbing boundary conditions implemented via PML

# by Leonardo Zepeda-Nunez March 2015 @ MIT

function FDweights(z,x,m::Int)
#---------------------------------
# finite-difference weights
# (Fornberg algorithm)
#
# z:  expansion point
# x:  vector of evaluation points
# m:  order of derivative
#
# Example: cwei = FDweights(0,[0 1 2],1);
# gives    cwei = [-3/2  2  -1/2]
#
# h f'_0 = -3/2 f_0 + 2 f_1 - 1/2 f_2
#
#---------------------------------

  n  = length(x)-1;
  c1 = 1;
  c4 = x[1]-z;
  c = zeros(n+1,m+1);
  c[1,1] = 1;
  for i=1:n
    mn = min(i,m);
    c2 = 1;
    c5 = c4;
    c4 = x[i+1]-z;
    for j=0:i-1
      c3 = x[i+1]-x[j+1];
      c2 = c2*c3;
      if (j == (i-1))
        for k=mn:-1:1
          c[i+1,k+1] = c1*(k*c[i,k]-c5*c[i,k+1])/c2;
        end
        c[i+1,1] = -c1*c5*c[i,1]/c2;
      end
      for k=mn:-1:1
        c[j+1,k+1] = (c4*c[j+1,k+1]-k*c[j+1,k])/c3;
      end
      c[j+1,1] = c4*c[j+1,1]/c3;
    end
    c1 = c2;
  end
  cwei = c[:,end];
  return cwei
end

function FirstOrderDifferenceMatrix1d(nx::Int,h::Float64,order::Int)
  # function Dx = FirstOrderDifferenceMatrix1d(nx,h::Float64,order::Int)
  # fucntion to compute the first order Finiter different matrix using
  # the Fornberg algorithm to compute the stencil, with a descentered stencil
  # at the edges
  # input :   nx    size of the matrix to be generated
  #           h     discretization step used
  #           order order of the discretizatin used
  # output:   Dx    Finite difference Matrix

  #computing the FD matrix at the interior
  diagonals = repmat(FDweights(order/2,linspace(0,order,order+1),1)'/h,nx,1);
  Dx = spzeros(nx,nx);
  for ii = 1:order+1
    bound = abs(-round(Integer,order/2)-1+ii);
    Dx = Dx + spdiagm(diagonals[1:end-round(Integer,bound),ii] ,-round(Integer,order/2)-1+ii,nx,nx);
  end

  # modifing the matrix at the boundaries, using descentered stencils
   for ii = 1:(round(Integer,order/2)-1)
     weights = FDweights(ii,linspace(0,order+2,order+3),1)'/h;
     Dx[ii,1:order+2]=weights[2:end];

     weights = FDweights(order+2-ii,linspace(0,order+2,order+3),1)'/h;
     Dx[end-(ii-1),(end-(order+1)):end]=weights[1:end-1];
   end
  return Dx
end

function stiffnessMatrix(nx::Int, dx::Float64, order::Int)
  # function Dx = stiffnessMatrix(nx, dx::Float64, order::Int)
  # function to compute a 1D stiffness Matris using  finite differences
  # (Dirichlet boundary nodes are not on the grid)
  # input :   nx    size of the matrix to be generated
  #           dx    discretization step used
  #           order order of the discretizatin used
  # output:   Dxx   Stiffness Matrix

  # computing the FD matrix at the interior
  diagonals = repmat(FDweights(order/2,linspace(0,order,order+1),2)'/(dx^2),nx,1);
  Dxx = spzeros(nx,nx);
  for ii = 1:order+1
    bound = abs(-round(Integer,order/2)-1+ii);
    Dxx = Dxx + spdiagm(diagonals[1:end-round(Integer,bound),ii] ,-round(Integer,order/2)-1+ii,nx,nx);
  end

  # modifing the matrix at the boundaries to obtain an uniform accuracy
  # using descentered stencil at the boundaries
  for ii = 1:(round(Integer,order/2)-1)
     weights = FDweights(ii,linspace(0,order+2,order+3),2)'/(dx^2);
     Dxx[ii,1:order+2]=weights[2:end];

     weights = FDweights(order+2-ii,linspace(0,order+2,order+3),2)'/(dx^2);
     Dxx[end-(ii-1),(end-(order+1)):end]=weights[1:end-1];
   end

   return Dxx
end

function DistribPML(nx::Int,ny::Int, nPML::Int,fac)
  # function (sigmaX, sigmaY, sigmaZ) = DistribPML(nx,ny,nz, nPML,fac)
  # function to create the damping profile on the PML's

  sigmaX = [fac*sigma(i,nx,nPML)+0*j for i=1:nx,j=1:ny];
  sigmaY = [fac*sigma(j,ny,nPML)+0*i for i=1:nx,j=1:ny];

  return (sigmaX, sigmaY)
end

function DistribPMLDerivative(nx::Int,ny::Int,nPML::Int,fac)
  # function to create the derivative of the damping profile on the PML's

  DxsigmaX = [fac*Dxsigma(i,nx,nPML)+0*j for i=1:nx,j=1:ny ];
  DysigmaY = [fac*Dxsigma(j,ny,nPML)+0*i for i=1:nx,j=1:ny ];

  return (DxsigmaX, DysigmaY)
end

function sigma(i,n,nPML)
  # pml function, we start one point after the boundary to enforce the continuity
  res = (i.<nPML).*((i-nPML).^2)/(nPML-1).^2 + (i.> (n-nPML+1)).*((i-n+nPML-1).^2)/(nPML-1).^2;
end

function Dxsigma(i,n,nPML)
  # derivative of the pml function, we start one point after the boundary to enforce the continuity
  res = -2*(i.<nPML).*((i-nPML))/(nPML-1).^2 + 2*(i.> (n-nPML+1)).*((i-n+nPML-1))/(nPML-1).^2;
end

function HelmholtzMatrix(m,nx::Int,ny::Int,npml::Int,h,fac,order::Int,omega)
  #function HelmholtzMatrix(m,nx,ny,nz,npml,h,fac,order,omega)
  #
  # H = -(\triangle + \omega^2 m I)
  # total number of degrees of freedom
  n = nx*ny;

  (sx,sy)    = DistribPML(nx,ny,npml,fac);
  (dsx,dsy) = DistribPMLDerivative(nx,ny,npml,fac);
  # assembling the 1-dimensional stiffness matrices
  Dxx1d = stiffnessMatrix(nx,h,order);
  Dyy1d = stiffnessMatrix(ny,h,order);
  # assembking the 1-dimensional finite difference matrices
  Dx1d  = FirstOrderDifferenceMatrix1d(nx,h,order);
  Dy1d  = FirstOrderDifferenceMatrix1d(ny,h,order);

  # assembling the 2D matrices using Kronecker products
  Dx    = kron(speye(ny),Dx1d    );
  Dy    = kron(Dy1d     ,speye(nx));

  Dxx   = kron(speye(ny), Dxx1d    );
  Dyy   = kron(Dyy1d    ,speye(nx) );

  # assembling the slowness matrix
  M     = spdiagm(m[:],0,n,n);

# assembling the Helmholtz matrix
  H = - omega^2*M +
        spdiagm(-1.im/(omega*(npml-1)*h)*dsx[:]./(1-1.im/omega*sx[:]).^3,0,n,n)*Dx +
        spdiagm(-1.im/(omega*(npml-1)*h)*dsy[:]./(1-1.im/omega*sy[:]).^3,0,n,n)*Dy -
        spdiagm(1./(1-1.im/omega*sx[:]).^2,0,n,n)*Dxx-
        spdiagm(1./(1-1.im/omega*sy[:]).^2,0,n,n)*Dyy;
  return H;
end


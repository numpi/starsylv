program starsylv
  integer :: n, r, i, iargc
  character(LEN=32) :: arg

 if ( iargc() > 0 ) then
     call getarg(1, arg)
     read(arg, *) n
     call getarg(2, arg)
     read(arg, *) r
     call test_equation(n, r)
  else
     do i = 4, 12
        n = 2**i
        call test_equation(n, 32)
     end do
  end if

end program starsylv

subroutine test_equation(n, r)

  implicit none

  ! Declarations
  integer :: n, r
  real*8, allocatable :: A(:,:,:), B(:,:,:), C(:,:,:), D(:,:,:), E(:,:,:), X(:,:,:)
  real*8 :: res, start, end

  ! Generate random triangular matrices such that we can
  ! solve the problem. 
  print *, 'Generating random problem of size N = ', n, ' R =', r
  allocate(A(n,n,r), B(n,n,r), C(n,n,r), D(n,n,r), E(n,n,r), X(n,n,r))
  call generate_random_problem(A,B,C,D,E,n,r)

  ! Solve the equation
  call cpu_time(start)
  call solve_tsylvester(A, B, C, D, E, X, n, r)
  call cpu_time(end)
  print *, 'Time:', end - start

  ! Evaluate residual
  ! call evaluate_tsylvester(A, B, C, D, E, X, n, r, res)
  call relative_residual(A, B, C, D, E, X, n, r, res)

  print *, 'Residual:', res

  ! Deallocate all the storage
  deallocate(A,B,C,D,E,X)

end subroutine test_equation

subroutine relative_residual(A, B, C, D, E, X, n, r, res)
  implicit none

  integer :: n, r, k, i, j
  real*8 :: res, Mnorm, nrm1, nrm2
  real*8 :: A(n,n,r), B(n,n,r), C(n,n,r), D(n,n,r), E(n,n,r), X(n,n,r)

  call evaluate_tsylvester(A, B, C, D, E, X, n, r, res)
  
  ! Compute the Frobenius norm of the large matrix      
  Mnorm = 0.d0
  do k = 1, r
     ! Compute the Frobenius norm of A_k \otimes B_k
     nrm1 = 0.d0
     nrm2 = 0.d0
     
     do i = 1, n
        do j = 1, n
           nrm1 = nrm1 + abs(A(i,j,k))**2
           nrm2 = nrm2 + abs(B(i,j,k))**2
        end do
     end do
     
     Mnorm = Mnorm + nrm1 * nrm2

     do i = 1, n
        do j = 1, n
           nrm1 = nrm1 + abs(C(i,j,k))**2
           nrm2 = nrm2 + abs(D(i,j,k))**2
        end do
     end do
     
     Mnorm = Mnorm + nrm1 * nrm2
  end do
  
  ! Given a lower bound to the 2-norm of M by frobnorm(M) / n
  Mnorm = dsqrt(Mnorm) / ( n * dsqrt(r * 1.d0) )
  
  ! Evaluate the norm of the solution X
  nrm1 = 0.d0
  do i = 1, n
     do j = 1, n
        do k = 1, r
           nrm1 = nrm1 + abs(X(i,j,k))**2
        end do
     end do
  end do
  
  nrm1 = dsqrt(nrm1)
  
  res = res / (Mnorm * nrm1)
end subroutine relative_residual

subroutine evaluate_tsylvester(A, B, C, D, E, X, n, r, res)
  implicit none

  integer :: n, r, k, i, j
  real*8 :: res
  real*8 :: A(n,n,r), B(n,n,r), C(n,n,r), D(n,n,r), E(n,n,r), X(n,n,r)
  real*8, allocatable :: F(:,:,:)

  allocate(F(n,n,r))

  do k = 1, r - 1
     F(:,:,k) = MATMUL(A(:,:,k), MATMUL(X(:,:,k), B(:,:,k))) + & 
          MATMUL(C(:,:,k), MATMUL(X(:,:,k+1), D(:,:,k))) - E(:,:,k)
  end do

  F(:,:,r) = MATMUL(A(:,:,r), MATMUL(X(:,:,r), B(:,:,r))) + & 
       MATMUL(C(:,:,r), MATMUL(TRANSPOSE(X(:,:,1)), D(:,:,r))) - E(:,:,r)

  ! Evaluate the residual in Frobenius norm
  res = 0
  do k = 1, r
     do i = 1, n
        do j = 1, n
           res = res + F(i,j,k)**2
        end do
     end do
  end do
  
  res = dsqrt(res)    
  deallocate(F)
end subroutine evaluate_tsylvester

subroutine generate_random_problem(A, B, C, D, E, n, r)
  integer :: n, r, i, j, k
  real*8 :: A(n,n,r), B(n,n,r), C(n,n,r), D(n,n,r), E(n,n,r), x

  A = 0
  B = 0
  C = 0
  D = 0
  E = 0

  do k = 1, r
     do i = 1, n
        do j = i, n
           ! We need to fill entirely the matrix E, whilst all the other are triangular
           call random_number(x) 
           E(i,j,k) = x
           call random_number(x) 
           E(j,i,k) = x

           call random_number(x)
           A(i,j,k) = x
           call random_number(x)
           B(j,i,k) = x
           call random_number(x)
           C(i,j,k) = x
           call random_number(x)
           D(j,i,k) = x
        end do
     end do
  end do

  do i = 1, n
     do k = 1, r
        C(i,i,k) = C(i,i,k) + dsqrt(n * 1.0d0)
        D(i,i,k) = D(i,i,k) + dsqrt(n * 1.0d0)
     end do
  end do

end subroutine generate_random_problem
  	
subroutine solve_tsylvester(A, B, C, D, E, X, n, r)
  implicit none
  
  integer :: k, i, n, r, l
  real*8 :: A(n,n,r), B(n,n,r), C(n,n,r)
  real*8 :: D(n,n,r), E(n,n,r), X(n,n,r)
  real*8, allocatable :: XD(:,:,:), XB(:,:,:)
  real*8, allocatable :: Md(:), Ms(:), f(:), xx(:)
  
  allocate(XB(n,n,r), XD(n,n,r))
  allocate(Md(2*r), Ms(2*r), f(2*r), xx(2*r))
    
  ! Copy the input data into the newly allocated matrices
  X = 0
  
  do k = n, 1, -1
     ! Handle the diagonal element
     call BuildSystem(A, B, C, D, E, X, n, r, k, k, Md, Ms, f, XB, XD)    
     call SolveBidiag(Md, Ms, r, f, X(k,k,:))

     do i = 1, r
        XB(k,k,i) = XB(k,k,i) + X(k,k,i) * B(k,k,i)
        XD(k,k,i) = XD(k,k,i) + X(k,k,1+mod(i,r)) * D(k,k,i)
     end do

     do i = k - 1, 1, -1
       call BuildSystem(A, B, C, D, E, X, n, r, i, k, Md, Ms, f, XB, XD)
       call SolveBidiag(Md, Ms, 2 * r, f, xx)

       X(i,k,:) = xx(1:r)
       X(k,i,:) = xx(r+1:2*r)

       do l = 1, r - 1
          XB(i,k,l) = XB(i,k,l) + X(i,k,l) * B(k,k,l)
          XD(i,k,l) = XD(i,k,l) + X(i,k,l+1) * D(k,k,l)
          XB(k,i,l) = XB(k,i,l) + X(k,i,l) * B(i,i,l)
          XD(k,i,l) = XD(k,i,l) + X(k,i,l+1) * D(i,i,l)
       end do

       XB(i,k,r) = XB(i,k,r) + X(i,k,r) * B(k,k,r)
       XD(i,k,r) = XD(i,k,r) + X(k,i,1) * D(k,k,r)
       XB(k,i,r) = XB(k,i,r) + X(k,i,r) * B(i,i,r)
       XD(k,i,r) = XD(k,i,r) + X(i,k,1) * D(i,i,r)
     end do
  end do

  deallocate(Md, Ms, f, xx)
  deallocate(XB, XD)

end subroutine solve_tsylvester

subroutine SolveBidiag(Md, Ms, n, f, x)
  implicit none

  integer :: n, i
  real*8 :: Md(n), Ms(n), f(n), x(n), c, s, lc(n)

  ! lc contains the last column of the matrix
  lc = 0
  lc(n) = Md(n)
  lc(n-1) = Ms(n-1)

  ! Perform a QR decomposition of the bidiagonal matrix, 
  ! and accumulate the transformations on f. 
  do i = 1, n - 1
    call drotg(Md(i), Ms(n), c, s)
    Ms(n) = 0
    call drot(1, Ms(i), 1, Ms(n), 1, c, s)
    call drot(1, lc(i), 1, lc(n), 1, c, s)
    call drot(1, f(i), 1, f(n), 1, c, s)
  end do

  ! Solve for x
  x(n) = f(n) / lc(n)
  x(n-1) = ( f(n-1) - lc(n-1) * x(n) ) / Md(n-1)

  do i = n-2, 1, -1
    x(i) = ( f(i) - lc(i) * x(n) - Ms(i) * x(i+1) ) / Md(i)
  end do
end subroutine SolveBidiag

subroutine BuildSystem(A, B, C, D, E, X, n, r, i, j, Md, Ms, f, XB, XD)
  ! Build the linear system associated with the computation of the
  ! element in position (i,j) and (j,i). The bidiagonal matrix is stored
  ! in Md (diagonal elements) and Ms(superdiagonal elements plus the
  ! one at the bottom-left corner).
  !
  ! The right-handside is stored inside the vector f.

  implicit none
  
  ! Variables
  integer n, r, i, j, k, l
  real*8 :: A(n,n,r), B(n,n,r), C(n,n,r), D(n,n,r), X(n,n,r)
  real*8 :: E(n,n,r), Md(2*r), Ms(2*r), f(2*r)
  real*8 :: XB(n,n,r), XD(n,n,r), ddot, g

  ! Might be useful to debug the code
  ! character(len=5120) :: line

  ! Building the linear system
  Md = 0
  Ms = 0
  f  = 0

  if (i .eq. j) then
     do k = 1, r
        Md(k) = B(i,j,k) * A(i,j,k)
        Ms(k) = C(i,j,k) * D(i,j,k)
     end do

     do k = 1, r - 1
        XB(i,j,k) = ddot(n-i+1, X(i,i,k), n, B(i,i,k), 1)
        XD(i,j,k) = ddot(n-i+1, X(i,i,k+1), n, D(i,i,k), 1)
        f(k) = E(i,j,k) - &
             ddot(n - i, A(i,i+1,k), n, XB(i+1,j,k), 1) - & 
             A(i,i,k) * XB(i,j,k) - &
             ddot(n-i, C(i,i+1,k), n, XD(i+1,j,k), 1) - &
             C(i,i,k) * XD(i,j,k)
     end do

     XB(i,j,r) = ddot(n-i+1, X(i,i,r), n, B(i,i,r), 1)
     XD(i,j,r) = ddot(n-i+1, X(i,i,1), 1, D(i,i,r), 1)

     f(r) = E(i,j,r) - &
          ddot(n-i, A(i,i+1,r), n, XB(i+1,j,r), 1) - &
          A(i,i,r) * XB(i,j,r) - &
          ddot(n-i, C(i,i+1,r), n, XD(i+1,i,r), 1) - &
          C(i,i,r) * XD(i,j,r)
  else
     ! Handle the case i ~= j, which is a 2r x 2r system
     do k = 1, r
        Md(k) = B(j,j,k) * A(i,i,k)
        Ms(k) = C(i,i,k) * D(j,j,k)
     end do

     do k = 1, r
        Md(r+k) = B(i,i,k) * A(j,j,k)
        Ms(r+k) = C(j,j,k) * D(i,i,k)
     end do

     do k = 1, r - 1
        XB(i,j,k) = ddot(n-j+1, X(i,j,k), n, B(j,j,k), 1)
        XD(i,j,k) = ddot(n-j+1, X(i,j,k+1), n, D(j,j,k), 1)
        f(k) = E(i,j,k) - &
             ddot(n-i, A(i,i+1,k), n, XB(i+1,j,k), 1) - &
             A(i,i,k) *  XB(i,j,k) - &
             ddot(n-i, C(i,i+1,k), n, XD(i+1,j,k), 1) - &
             C(i,i,k) * XD(i,j,k)
     end do

     XB(i,j,r) = ddot(n-j+1, X(i,j,r), n, B(j,j,r), 1)
     XD(i,j,r) = ddot(n-j+1, X(j,i,1), 1, D(j,j,r), 1)

     f(r) = E(i,j,r) - &
          ddot(n-i, A(i,i+1,r), n, XB(i+1,j,r), 1) - &
          A(i,i,r) * XB(i,j,r) - &
          ddot(n-i, C(i,i+1,r), n, XD(i+1,j,r), 1) - &
          C(i,i,r) * XD(i,j,r)

     do k = 1, r - 1
        XB(j,i,k) = ddot(n-i+1, X(j,i,k), n, B(i,i,k), 1)
        XD(j,i,k) = ddot(n-i+1, X(j,i,k+1), n, D(i,i,k), 1)
        f(r+k) = E(j,i,k) - &
             ddot(n-j, A(j,j+1,k), n, XB(j+1,i,k), 1) - &
             A(j,j,k) *  XB(j,i,k) - &
             ddot(n-j, C(j,j+1,k), n, XD(j+1,i,k), 1) - &
             C(j,j,k) * XD(j,i,k)
     end do

     XB(j,i,r) = ddot(n-i+1, X(j,i,r), n, B(i,i,r), 1)
     XD(j,i,r) = ddot(n-i+1, X(i,j,1), 1, D(i,i,r), 1)
     
     f(r+k) = E(j,i,r) - &
          ddot(n-j, A(j,j+1,r), n, XB(j+1, i, r), 1) - &
          A(j,j,r) *  XB(j,i,r) - &
          ddot(n-j, C(j,j+1,r), n, XD(j+1,i,r), 1) - &
          C(j,j,r) * XD(j,i,r)
  end if
  
  
end subroutine BuildSystem

! This file contains an example that can be used to test the starsylv
! implementation.

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
  call starsylv_solve_tsylvester(A, B, C, D, E, X, n, r)
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
  

!-----------------------------------------------------------------------------
! main.f90
! Test of Cholesky decomposition based on Cholesky-Crout algorithm
! Author: Tianbai Xiao
! Licence: GPLv3
! > A is a n by n positive definite real symmetric matrix
! > A = L L^T, where L is a lower triangular matrix and L^T is transpose of L.
!-----------------------------------------------------------------------------

program main
  
implicit none
real*8 :: A(5,5)
real*8 :: L(5,5)
  
A = reshape((/ 5.0d0, 4.0d0, 3.0d0, 2.0d0, 1.0d0, &
     &         4.0d0, 4.0d0, 3.0d0, 2.0d0, 1.0d0, &
     &         3.0d0, 3.0d0, 3.0d0, 2.0d0, 1.0d0, &
     &         2.0d0, 2.0d0, 2.0d0, 2.0d0, 1.0d0, &
     &         1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0 /), (/5,5/))

call calc_cholesky(5, A, L)
  
write(6,'(5f6.2)') L

end program main


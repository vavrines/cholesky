!-----------------------------------------------------------------------------
! cholesky.f90
! Calculate Cholesky decomposition based on Cholesky-Crout algorithm
! Author: Tianbai Xiao
! Licence: GPLv3
! > A is a n by n positive definite real symmetric matrix
! > A = L L^T, where L is a lower triangular matrix and L^T is transpose of L.
!-----------------------------------------------------------------------------

subroutine calc_cholesky(n, A, L)

implicit none
integer, intent(in) :: n
real(kind=8), intent(in) :: A(n,n)
real(kind=8), intent(out) :: L(n,n)
integer :: i,j

! check of positive definite
do i = 1, n
    if (A(i,i).le.0.0d0) then
        return
    end if
end do

! Cholesky-Crout algorithm
L(:,:) = 0.d0
do j = 1, n
    L(j,j) = sqrt( A(j,j) - dot_product(L(j,1:j-1),L(j,1:j-1)) )
    do i = j+1, n
        L(i,j)  = ( A(i,j) - dot_product(L(i,1:j-1),L(j,1:j-1)) ) / L(j,j)
    end do
end do

end subroutine calc_cholesky

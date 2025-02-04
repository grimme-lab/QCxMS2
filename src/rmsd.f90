! This file is  taken from mctc-rmsd.
! SPDX-Identifier: LGPL-3.0-or-later
!
! Copyright (C) 2004, 2005 Chaok Seok, Evangelos Coutsias and Ken Dill
!      UCSF, Univeristy of New Mexico, Seoul National University
! Written by Chaok Seok and Evangelos Coutsias 2004.
!
! Modified and adapted for mctc-rmsd by Sebastian Ehlert.
!
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version.
!
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with this library; if not, write to the Free Software
! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

!> Implementation of the least-square RMSD fit of two structures
module rmsd_ls
   use xtb_mctc_accuracy, only : dp
   implicit none


contains
!> This subroutine calculates the least-square RMSD of two coordinate
!> sets coord1(3, n) and coord2(3, n) using a method based on quaternion.
pure subroutine get_rmsd_for_coord(coord1, coord2, rmsd, gradient, trafo, mask)

   !> Coordinates of the reference structure
   real(dp), intent(in) :: coord1(:, :)

   !> Coordinates of the structure to compare
   real(dp), intent(in) :: coord2(:, :)

   !> Root mean square deviation between the input structure
   real(dp), intent(out) :: rmsd

   !> Derivative of the RMSD w.r.t. the coordinates of the reference structure
   real(dp), intent(out), optional :: gradient(:, :)

   !> Transformation matrix to rotate into reference structure
   real(dp), intent(out), optional :: trafo(:, :)

   !> Atoms to include into RMSD calculation
   logical, intent(in), optional :: mask(:)

   integer :: n, m, i, j
   real(dp), allocatable :: x(:, :), y(:, :)
   real(dp), allocatable :: xi(:), yi(:)
   real(dp) :: x_center(3), y_center(3)
   real(dp) :: x_norm, y_norm, lambda
   real(dp) :: Rmatrix(3, 3), tmp(3)
   real(dp) :: S(4, 4), q(4)

   if (present(mask)) then
      m = min(size(coord1, 2), size(coord2, 2), size(mask))
      n = count(mask(1:m))
      allocate(x(3, n), y(3, n), xi(n), yi(n))
      ! make copies of the original coordinates
      j = 0
      do i = 1, m
         if (mask(i)) then
            j = j + 1
            x(:, j) = coord1(:, i)
            y(:, j) = coord2(:, i)
         end if
      end do
   else
      n = min(size(coord1, 2), size(coord2, 2))
      allocate(x(3, n), y(3, n), xi(n), yi(n))
      ! make copies of the original coordinates
      x(:, 1:n) = coord1(:, 1:n)
      y(:, 1:n) = coord2(:, 1:n)
   end if

   ! calculate the barycenters
   x_center(:) = 0.0_dp
   y_center(:) = 0.0_dp
   do i = 1, n
      x_center(:) = x_center(:) + x(:, i)/real(n, dp)
      y_center(:) = y_center(:) + y(:, i)/real(n, dp)
   end do

   ! calculate centroidal coordinates and the norms
   x_norm = 0.0_dp
   y_norm = 0.0_dp
   do i = 1, n
      x(:, i) = x(:, i) - x_center
      y(:, i) = y(:, i) - y_center
      x_norm = x_norm + dot_product(x(:, i), x(:, i))
      y_norm = y_norm + dot_product(y(:, i), y(:, i))
   end do

   ! calculate the R matrix
   Rmatrix(:, :) = matmul(x, transpose(y))

   ! S matrix
   S(1, 1) = Rmatrix(1, 1) + Rmatrix(2, 2) + Rmatrix(3, 3)
   S(2, 1) = Rmatrix(2, 3) - Rmatrix(3, 2)
   S(3, 1) = Rmatrix(3, 1) - Rmatrix(1, 3)
   S(4, 1) = Rmatrix(1, 2) - Rmatrix(2, 1)

   S(1, 2) = S(2, 1)
   S(2, 2) = Rmatrix(1, 1) - Rmatrix(2, 2) - Rmatrix(3, 3)
   S(3, 2) = Rmatrix(1, 2) + Rmatrix(2, 1)
   S(4, 2) = Rmatrix(1, 3) + Rmatrix(3, 1)

   S(1, 3) = S(3, 1)
   S(2, 3) = S(3, 2)
   S(3, 3) = -Rmatrix(1, 1) + Rmatrix(2, 2) - Rmatrix(3, 3)
   S(4, 3) = Rmatrix(2, 3) + Rmatrix(3, 2)

   S(1, 4) = S(4, 1)
   S(2, 4) = S(4, 2)
   S(3, 4) = S(4, 3)
   S(4, 4) = -Rmatrix(1, 1) - Rmatrix(2, 2) + Rmatrix(3, 3)

   ! Calculate eigenvalues and eigenvectors, and
   ! take the maximum eigenvalue lambda and the corresponding eigenvector q.
   call dstmev(S, lambda, q)

   if (present(trafo) .or. present(gradient)) then
      ! convert quaternion q to rotation matrix
      call rotation_matrix(q, Rmatrix)
      if (present(trafo)) trafo(:, :) = Rmatrix
   end if

   ! RMS Deviation, small number added to avoid NAN in gradient, SG 12/18
   rmsd = sqrt(max(0.0_dp, ((x_norm + y_norm) - 2.0_dp*lambda))/real(n, dp))

   if (present(gradient)) then
      if (present(mask)) then
         j = 0
         gradient(:, :) = 0.0_dp
         do i = 1, n
            if (mask(i)) then
               j = j + 1
               tmp(:) = matmul(transpose(Rmatrix), y(:, i))
               gradient(:, j) = (x(:, i) - tmp)/max(epsilon(0.0_dp), rmsd*real(n, dp))
            else
            end if
         end do
      else
         do i = 1, n
            tmp(:) = matmul(transpose(Rmatrix), y(:, i))
            gradient(:, i) = (x(:, i) - tmp)/max(epsilon(0.0_dp), rmsd*real(n, dp))
         end do
      end if
   end if

end subroutine get_rmsd_for_coord


!> This subroutine constructs rotation matrix U from quaternion q.
pure subroutine rotation_matrix(q, U)
   real(dp), intent(in) :: q(:)
   real(dp), intent(out) :: U(:, :)
   real(dp) :: q0, q1, q2, q3, b0, b1, b2, b3
   real(dp) :: q00, q01, q02, q03, q11, q12, q13, q22, q23, q33

   q0 = q(1)
   q1 = q(2)
   q2 = q(3)
   q3 = q(4)

   b0 = 2.0_dp*q0
   b1 = 2.0_dp*q1
   b2 = 2.0_dp*q2
   b3 = 2.0_dp*q3

   q00 = b0*q0 - 1.0_dp
   q01 = b0*q1
   q02 = b0*q2
   q03 = b0*q3

   q11 = b1*q1
   q12 = b1*q2
   q13 = b1*q3

   q22 = b2*q2
   q23 = b2*q3

   q33 = b3*q3

   U(1, 1) = q00 + q11
   U(1, 2) = q12 - q03
   U(1, 3) = q13 + q02

   U(2, 1) = q12 + q03
   U(2, 2) = q00 + q22
   U(2, 3) = q23 - q01

   U(3, 1) = q13 - q02
   U(3, 2) = q23 + q01
   U(3, 3) = q00 + q33

end subroutine rotation_matrix


!> A simple subroutine to compute the leading eigenvalue and eigenvector
!> of a symmetric, traceless 4x4 matrix A by an inverse power iteration:
!> (1) the matrix is converted to tridiagonal form by 3 Givens
!> rotations;  V*A*V' = T
!> (2) Gershgorin's theorem is used to estimate a lower
!> bound for the leading negative eigenvalue:
!> lambda_1 > g=min(T11-t12, -t21+T22-t23, -t32+T33-t34, -t43+T44)
!>          =
!> where tij=abs(Tij)
!> (3) Form the positive definite matrix
!>     B = T-gI
!> (4) Use svd (algorithm svdcmp from "Numerical Recipes")
!>     to compute eigenvalues and eigenvectors for SPD matrix B
!> (5) Shift spectrum back and keep leading singular vector
!>     and largest eigenvalue.
!> (6) Convert eigenvector to original matrix A, through
!>     multiplication by V'.
pure subroutine dstmev(A, lambda, evec)
   real(dp), intent(inout) :: A(4, 4)
   real(dp), intent(out) :: evec(4)
   real(dp), intent(out) :: lambda
   real(dp) :: T(4, 4), V(4, 4), SV(4, 4)
   integer :: i
   integer :: max_loc(1) ! must be an array
   real(dp) :: SW(4)
   real(dp) :: rv1(8)

   ! (I). Convert to tridiagonal form, keeping similarity transform
   ! (a product of 3 Givens rotations)
   call givens4(A, T, V)

   ! (II) Estimate lower bound of smallest eigenvalue by Gershgorin's theorem
   lambda = min(T(1, 1) - abs(T(1, 2)), -abs(T(2, 1)) + T(2, 2) - abs(T(2, 3)), &
      -abs(T(3, 2)) + T(3, 3) - abs(T(3, 4)), -abs(T(4, 3)) + T(4, 4))

   ! (III). Form positive definite matrix     T <== lambda*I - T
   do i = 1, 4
      T(i, i) = T(i, i) - lambda
   enddo

   ! (IV). Compute singular values/vectors of SPD matrix B
   call svdcmp(4, T, 4, 4, SW, SV, rv1)

    !(V). Shift spectrum back
   max_loc = maxloc(SW)
   lambda = SW(max_loc(1)) + lambda

   ! (VI). Convert eigenvector to original coordinates: (V is transposed!)
   evec = matmul(V, SV(:, max_loc(1)))

end subroutine dstmev


!> Performs givens rotations to reduce symmetric 4x4 matrix to tridiagonal
pure subroutine givens4(S, T, V)
   real(dp), intent(in) :: S(4, 4)
   real(dp), intent(out) :: T(4, 4)
   real(dp), intent(out) :: V(4, 4)
   real(dp)  :: c1, c2, c3, s1, s2, s3, r1, r2, r3, c1c2, s1c2

   T = S
   V = 0.0_dp

   ! Zero out entries T(4, 1) and T(1, 4)
   ! compute cos and sin of rotation angle in the 3-4 plane
   r1 = pythag(T(3, 1), T(4, 1))
   if (r1 .ne. 0.0_dp) then
      c1 = T(3, 1)/r1
      s1 = T(4, 1)/r1
      V(3, 3) = c1
      V(3, 4) = s1
      V(4, 3) = -s1
      V(4, 4) = c1
      T(3, 1) = r1
      T(4, 1) = 0.0_dp
      T(3:4, 2:4) = matmul(V(3:4, 3:4), T(3:4, 2:4))
      T(1:2, 3:4) = transpose(T(3:4, 1:2))
      T(3:4, 3:4) = matmul(T(3:4, 3:4), transpose(V(3:4, 3:4)))
   else
      c1 = 1.0_dp
      s1 = 0.0_dp
   endif

   ! Zero out entries T(3, 1) and T(1, 3)
   ! compute cos and sin of rotation angle in the 2-3 plane
   r2 = pythag(T(3, 1), T(2, 1))
   if (r2 .ne. 0.0_dp) then
      c2 = T(2, 1)/r2
      s2 = T(3, 1)/r2
      V(2, 2) = c2
      V(2, 3) = s2
      V(3, 2) = -s2
      V(3, 3) = c2
      T(2, 1) = r2
      T(3, 1) = 0.0_dp
      T(2:3, 2:4) = matmul(V(2:3, 2:3), T(2:3, 2:4))
      T(1, 2:3) = T(2:3, 1)
      T(4, 2:3) = T(2:3, 4)
      T(2:3, 2:3) = matmul(T(2:3, 2:3), transpose(V(2:3, 2:3)))
   else
      c2 = 1.0_dp
      s2 = 0.0_dp
   endif

   ! Zero out entries T(4, 2) and T(2, 4)
   ! compute cos and sin of rotation angle in the 3-4 plane
   r3 = pythag(T(4, 2), T(3, 2))
   if (r3 .ne. 0.0_dp) then
      c3 = T(3, 2)/r3
      s3 = T(4, 2)/r3
      V(3, 3) = c3
      V(3, 4) = s3
      V(4, 3) = -s3
      V(4, 4) = c3
      T(3, 2) = r3
      T(4, 2) = 0.0_dp
      T(3:4, 3:4) = matmul(V(3:4, 3:4), T(3:4, 3:4))
      T(1:2, 3:4) = transpose(T(3:4, 1:2))
      T(3:4, 3:4) = matmul(T(3:4, 3:4), transpose(V(3:4, 3:4)))
   else
      c3 = 1.0_dp
      s3 = 0.0_dp
   endif

   ! Compute net rotation matrix (accumulate similarity for evec. computation)
   ! To save transposing later, This is the transpose!
   V(1, 1) = 1.0_dp
   V(1, 2:4) = 0.0_dp
   V(2:4, 1) = 0.0_dp
   V(2, 2) = c2
   V(3, 2) = c1*s2
   V(4, 2) = s1*s2
   c1c2 = c1*c2
   s1c2 = s1*c2
   V(2, 3) = -s2*c3
   V(3, 3) = c1c2*c3 - s1*s3
   V(4, 3) = s1c2*c3 + c1*s3
   V(2, 4) = s2*s3
   V(3, 4) = -c1c2*s3 - s1*c3
   V(4, 4) = -s1c2*s3 + c1*c3

end subroutine givens4


pure subroutine svdcmp(mmax, a, m, n, w, v, rv1)

   integer, intent(in) :: mmax
   integer, intent(in) :: m
   integer, intent(in) :: n
   real(dp), intent(inout) :: a(mmax, *)
   real(dp), intent(inout) :: v(mmax, *)
   real(dp), intent(inout) :: w(*)
   real(dp), intent(inout) :: rv1(*)
   integer :: i, its, j, jj, k, l, nm
   real(dp) :: anorm, c, f, g, h, s, scale, x, y, z

   g = 0.0_dp
   scale = 0.0_dp
   anorm = 0.0_dp
   do i = 1, n
      l = i + 1
      rv1(i) = scale*g
      g = 0.0_dp
      s = 0.0_dp
      scale = 0.0_dp
      if (i .le. m) then
         do k = i, m
            scale = scale + abs(a(k, i))
         end do
         if (scale .ne. 0.0_dp) then
            do k = i, m
               a(k, i) = a(k, i)/scale
               s = s + a(k, i)*a(k, i)
            end do
            f = a(i, i)
            g = -sign(sqrt(s), f)
            h = f*g - s
            a(i, i) = f - g
            do j = l, n
               s = 0.0_dp
               do k = i, m
                  s = s + a(k, i)*a(k, j)
               end do
               f = s/h
               do k = i, m
                  a(k, j) = a(k, j) + f*a(k, i)
               end do
            end do
            do k = i, m
               a(k, i) = scale*a(k, i)
            end do
         endif
      endif
      w(i) = scale*g
      g = 0.0_dp
      s = 0.0_dp
      scale = 0.0_dp
      if ((i .le. m) .and. (i .ne. n)) then
         do k = l, n
            scale = scale + abs(a(i, k))
         end do
         if (scale .ne. 0.0_dp) then
            do k = l, n
               a(i, k) = a(i, k)/scale
               s = s + a(i, k)*a(i, k)
            end do
            f = a(i, l)
            g = -sign(sqrt(s), f)
            h = f*g - s
            a(i, l) = f - g
            do k = l, n
               rv1(k) = a(i, k)/h
            end do
            do j = l, m
               s = 0.0_dp
               do k = l, n
                  s = s + a(j, k)*a(i, k)
               end do
               do k = l, n
                  a(j, k) = a(j, k) + s*rv1(k)
               end do
            end do
            do k = l, n
               a(i, k) = scale*a(i, k)
            end do
         endif
      endif
      anorm = max(anorm, (abs(w(i)) + abs(rv1(i))))
   end do

   do i = n, 1, -1
      if (i .lt. n) then
         if (g .ne. 0.0_dp) then
            do j = l, n
               v(j, i) = (a(i, j)/a(i, l))/g
            end do
            do j = l, n
               s = 0.0_dp
               do k = l, n
                  s = s + a(i, k)*v(k, j)
               end do
               do k = l, n
                  v(k, j) = v(k, j) + s*v(k, i)
               end do
            end do
         endif
         do j = l, n
            v(i, j) = 0.0_dp
            v(j, i) = 0.0_dp
         end do
      endif
      v(i, i) = 1.0_dp
      g = rv1(i)
      l = i
   end do

   do i = min(m, n), 1, -1
      l = i + 1
      g = w(i)
      do j = l, n
         a(i, j) = 0.0_dp
      end do
      if (g .ne. 0.0_dp) then
         g = 1.0_dp/g
         do j = l, n
            s = 0.0_dp
            do k = l, m
               s = s + a(k, i)*a(k, j)
            end do
            f = (s/a(i, i))*g
            do k = i, m
               a(k, j) = a(k, j) + f*a(k, i)
            end do
         end do
         do j = i, m
            a(j, i) = a(j, i)*g
         end do
      else
         do j = i, m
            a(j, i) = 0.0_dp
         end do
      endif
      a(i, i) = a(i, i) + 1.0_dp
   end do

   do k = n, 1, -1
      do its = 1, 30
         do l = k, 1, -1
            nm = l - 1
            if ((abs(rv1(l)) + anorm) .eq. anorm) goto 2
            if ((abs(w(nm)) + anorm) .eq. anorm) goto 1
         end do
      1  c = 0.0_dp
         s = 1.0_dp
         do i = l, k
            f = s*rv1(i)
            rv1(i) = c*rv1(i)
            if ((abs(f) + anorm) .eq. anorm) goto 2
            g = w(i)
            h = pythag(f, g)
            w(i) = h
            h = 1.0_dp/h
            c = (g*h)
            s = -(f*h)
            do j = 1, m
               y = a(j, nm)
               z = a(j, i)
               a(j, nm) = (y*c) + (z*s)
               a(j, i) = -(y*s) + (z*c)
            end do
         end do
      2  z = w(k)
         if (l .eq. k) then
            if (z .lt. 0.0_dp) then
               w(k) = -z
               do j = 1, n
                  v(j, k) = -v(j, k)
               end do
            endif
            goto 3
         endif
         if (its .eq. 30) then
            error stop 'no convergence in svdcmp'
         endif
         x = w(l)
         nm = k - 1
         y = w(nm)
         g = rv1(nm)
         h = rv1(k)
         f = ((y - z)*(y + z) + (g - h)*(g + h))/(2.0_dp*h*y)
         g = pythag(f, 1.0_dp)
         f = ((x - z)*(x + z) + h*((y/(f + sign(g, f))) - h))/x
         c = 1.0_dp
         s = 1.0_dp
         do j = l, nm
            i = j + 1
            g = rv1(i)
            y = w(i)
            h = s*g
            g = c*g
            z = pythag(f, h)
            rv1(j) = z
            c = f/z
            s = h/z
            f = (x*c) + (g*s)
            g = -(x*s) + (g*c)
            h = y*s
            y = y*c
            do jj = 1, n
               x = v(jj, j)
               z = v(jj, i)
               v(jj, j) = (x*c) + (z*s)
               v(jj, i) = -(x*s) + (z*c)
            end do
            z = pythag(f, h)
            w(j) = z
            if (z .ne. 0.0_dp) then
               z = 1.0_dp/z
               c = f*z
               s = h*z
            endif
            f = (c*g) + (s*y)
            x = -(s*g) + (c*y)
            do jj = 1, m
               y = a(jj, j)
               z = a(jj, i)
               a(jj, j) = (y*c) + (z*s)
               a(jj, i) = -(y*s) + (z*c)
            end do
         end do
         rv1(l) = 0.0_dp
         rv1(k) = f
         w(k) = x
      end do
   3  continue
   end do

end subroutine svdcmp


elemental function pythag(a, b)
   real(dp), intent(in) :: a, b
   real(dp) :: pythag
   real(dp) :: absa, absb

   absa = abs(a)
   absb = abs(b)
   if (absa .gt. absb) then
      pythag = absa*dsqrt(1.0_dp + (absb/absa)**2)
   else
      if (absb .eq. 0.0_dp) then
         pythag = 0.0_dp
      else
         pythag = absb*dsqrt(1.0_dp + (absa/absb)**2)
      endif
   endif

end function pythag


end module rmsd_ls

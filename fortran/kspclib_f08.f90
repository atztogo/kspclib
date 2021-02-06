! Copyright (C) 2021 Atsushi Togo
! All rights reserved.

! This file is part of kspclib.

! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:

! * Redistributions of source code must retain the above copyright
!   notice, this list of conditions and the following disclaimer.

! * Redistributions in binary form must reproduce the above copyright
!   notice, this list of conditions and the following disclaimer in
!   the documentation and/or other materials provided with the
!   distribution.

! * Neither the name of the kspclib project nor the names of its
!   contributors may be used to endorse or promote products derived
!   from this software without specific prior written permission.

! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
! COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
! POSSIBILITY OF SUCH DAMAGE.

module kspclib_f08

  use iso_c_binding, only: c_char, c_int, c_long, c_double

  implicit none

  private

  interface

     function ksp_get_major_version() bind(c)
       import c_int
       integer(c_int) :: ksp_get_major_version
     end function ksp_get_major_version


     function ksp_get_minor_version() bind(c)
       import c_int
       integer(c_int) :: ksp_get_minor_version
     end function ksp_get_minor_version


     function ksp_get_micro_version() bind(c)
       import c_int
       integer(c_int) :: ksp_get_micro_version
     end function ksp_get_micro_version


     subroutine ksp_get_all_grid_addresses(grid_address, mesh) bind(c)
       import c_int
       integer(c_int), intent(in) :: mesh(3)
       integer(c_int), intent(inout) :: grid_address(3, *)
     end subroutine ksp_get_all_grid_addresses


     subroutine ksp_get_double_grid_address(address_double, address, mesh, &
          is_shift) bind(c)
       import c_int
       integer(c_int), intent(inout) :: address_double(3)
       integer(c_int), intent(in) :: address(3)
       integer(c_int), intent(in) :: mesh(3)
       integer(c_int), intent(in) :: is_shift(3)
     end subroutine ksp_get_double_grid_address


     function ksp_get_double_grid_index(address_double, mesh) bind(c)
       import c_int, c_long
       integer(c_int), intent(in) :: mesh(3)
       integer(c_int), intent(in) :: address_double(3)
       integer(c_long) :: ksp_get_double_grid_index
     end function ksp_get_double_grid_index


     subroutine ksp_get_thm_relative_grid_addresses(relative_grid_addresses, &
          rec_lattice) bind(c)
       import c_int, c_double
       integer(c_int), intent(inout) :: relative_grid_addresses(3, 4, 24)
       real(c_double), intent(in) :: rec_lattice(3, 3)
     end subroutine ksp_get_thm_relative_grid_addresses


     function ksp_get_thm_integration_weight(omega, &
          tetrahedra_omegas, function_char) bind(c)
       import c_char, c_double
       real(c_double), intent(in), value :: omega
       real(c_double), intent(in) :: tetrahedra_omegas(4, 24)
       character(kind=c_char), intent(in), value :: function_char
       real(c_double) :: ksp_get_thm_integration_weight
     end function ksp_get_thm_integration_weight


     function ksp_snf_transform_rotations(transformed_rots, &
          rotations, num_rot, D_diag, Q) bind(c)
       import c_long, c_int
       integer(c_long), intent(inout) :: transformed_rots(3, 3, *)
       integer(c_int), intent(in) :: rotations(3, 3, *)
       integer(c_int), intent(in), value :: num_rot
       integer(c_long), intent(in) :: D_diag(3)
       integer(c_long), intent(in) :: Q(3, 3)
     end function ksp_snf_transform_rotations


     function ksp_get_snf3x3(D_diag, P, Q, A) bind(c)
       import c_long, c_int
       integer(c_long), intent(inout) :: D_diag(3)
       integer(c_long), intent(inout) :: P(3, 3)
       integer(c_long), intent(inout) :: Q(3, 3)
       integer(c_long), intent(in) :: A(3, 3)
       integer(c_int) :: ksp_get_snf3x3
     end function ksp_get_snf3x3


     subroutine ksp_get_all_grgrid_addresses(grid_address, D_diag) bind(c)
       import c_long
       integer(c_long), intent(inout) :: grid_address(3, *)
       integer(c_long), intent(in) :: D_diag(3)
     end subroutine ksp_get_all_grgrid_addresses


     subroutine ksp_get_double_grgrid_address(address_double, address, &
          D_diag, PS) bind(c)
       import c_long
       integer(c_long), intent(inout) :: address_double(3)
       integer(c_long), intent(in) :: D_diag(3)
       integer(c_long), intent(in) :: address(3)
       integer(c_long), intent(in) :: PS(3)
     end subroutine ksp_get_double_grgrid_address


     function ksp_get_grgrid_index(address, D_diag) bind(c)
       import c_long
       integer(c_long), intent(in) :: address(3)
       integer(c_long), intent(in) :: D_diag(3)
       integer(c_long) :: ksp_get_grgrid_index
     end function ksp_get_grgrid_index


     function ksp_get_double_grgrid_index(address_double, D_diag, PS) bind(c)
       import c_long
       integer(c_long), intent(in) :: address_double(3)
       integer(c_long), intent(in) :: D_diag(3)
       integer(c_long), intent(in) :: PS(3)
       integer(c_long) :: ksp_get_double_grgrid_index
     end function ksp_get_double_grgrid_index


     subroutine ksp_get_grgrid_address_from_index(address, grid_index, D_diag) &
          bind(c)
       import c_long
       integer(c_long), intent(inout) :: address(3)
       integer(c_long), intent(in), value :: grid_index
       integer(c_long), intent(in) :: D_diag(3)
     end subroutine ksp_get_grgrid_address_from_index


     function ksp_rotate_grgrid_index(grid_index, rotation, D_diag, PS) bind(c)
       import c_long
       integer(c_long), intent(in), value :: grid_index
       integer(c_long), intent(in) :: rotation(3, 3)
       integer(c_long), intent(in) :: D_diag(3)
       integer(c_long), intent(in) :: PS(3)
       integer(c_long) :: ksp_rotate_grgrid_index
     end function ksp_rotate_grgrid_index


     subroutine ksp_get_ir_grgrid_map(ir_grid_indices, rotations, num_rot, &
          D_diag, PS) bind(c)
       import c_long, c_int
       integer(c_long), intent(inout) :: ir_grid_indices(*)
       integer(c_long), intent(in) :: rotations(3, 3, *)
       integer(c_int), intent(in), value :: num_rot
       integer(c_long), intent(in) :: D_diag(3)
       integer(c_long), intent(in) :: PS(3)
     end subroutine ksp_get_ir_grgrid_map

  end interface


  public :: ksp_get_major_version,  ksp_get_minor_version, &
       ksp_get_micro_version, ksp_get_all_grid_addresses, ksp_get_snf3x3, &
       ksp_get_double_grid_address, ksp_get_double_grid_index, &
       ksp_get_thm_relative_grid_addresses, ksp_get_thm_integration_weight, &
       ksp_snf_transform_rotations, ksp_get_all_grgrid_addresses, &
       ksp_get_double_grgrid_address, ksp_get_grgrid_index, &
       ksp_get_double_grgrid_index, ksp_get_grgrid_address_from_index, &
       ksp_rotate_grgrid_index, ksp_get_ir_grgrid_map


end module kspclib_f08

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

     ! function spg_get_symmetry( rotation, translation, max_size, lattice, &
     !      & position, types, num_atom, symprec) bind(c)
     !   import c_int, c_double
     !   integer(c_int), intent(inout) :: rotation(3,3,*)
     !   real(c_double), intent(inout) :: translation(3,*)
     !   integer(c_int), intent(in), value :: max_size
     !   real(c_double), intent(in) :: lattice(3,3), position(3,*)
     !   integer(c_int), intent(in) :: types(*)
     !   integer(c_int), intent(in), value :: num_atom
     !   real(c_double), intent(in), value :: symprec
     !   integer(c_int) :: spg_get_symmetry
     ! end function spg_get_symmetry


     ! function spgat_get_symmetry( rotation, translation, max_size, lattice, &
     !      & position, types, num_atom, symprec, angle_tolerance) bind(c)
     !   import c_int, c_double
     !   integer(c_int), intent(out) :: rotation(3,3,*)
     !   real(c_double), intent(out) :: translation(3,*)
     !   integer(c_int), intent(in), value :: max_size
     !   real(c_double), intent(in) :: lattice(3,3), position(3,*)
     !   integer(c_int), intent(in) :: types(*)
     !   integer(c_int), intent(in), value :: num_atom
     !   real(c_double), intent(in), value :: symprec, angle_tolerance
     !   integer(c_int) :: spgat_get_symmetry
     ! end function spgat_get_symmetry


     ! function spg_get_symmetry_with_collinear_spin( rotation, translation, &
     !      & equivalent_atoms, max_size, lattice, position, types, spins, &
     !      & num_atom, symprec) bind(c)
     !   import c_int, c_double
     !   integer(c_int), intent(out) :: rotation(3,3,*)
     !   real(c_double), intent(out) :: translation(3,*)
     !   integer(c_int), intent(out) :: equivalent_atoms(*)
     !   integer(c_int), intent(in), value :: max_size
     !   real(c_double), intent(in) :: lattice(3,3), position(3,*)
     !   integer(c_int), intent(in) :: types(*)
     !   real(c_double), intent(in) :: spins(*)
     !   integer(c_int), intent(in), value :: num_atom
     !   real(c_double), intent(in), value :: symprec
     !   integer(c_int) :: spg_get_symmetry_with_collinear_spin
     ! end function spg_get_symmetry_with_collinear_spin


     ! function spgat_get_symmetry_with_collinear_spin( rotation, translation, &
     !      & equivalent_atoms, max_size, lattice, position, types, spins, &
     !      & num_atom, symprec, angle_tolerance) bind(c)
     !   import c_int, c_double
     !   integer(c_int), intent(out) :: rotation(3,3,*)
     !   real(c_double), intent(out) :: translation(3,*)
     !   integer(c_int), intent(out) :: equivalent_atoms(*)
     !   integer(c_int), intent(in), value :: max_size
     !   real(c_double), intent(in) :: lattice(3,3), position(3,*)
     !   integer(c_int), intent(in) :: types(*)
     !   real(c_double), intent(in) :: spins(*)
     !   integer(c_int), intent(in), value :: num_atom
     !   real(c_double), intent(in), value :: symprec, angle_tolerance
     !   integer(c_int) :: spgat_get_symmetry_with_collinear_spin
     ! end function spgat_get_symmetry_with_collinear_spin


     ! function spg_get_symmetry_with_site_tensors( rotation, translation, &
     !      & equivalent_atoms, primitive_lattice, spin_flips, max_size, lattice, &
     !      & position, types, tensors, tensor_rank, num_atom, is_magnetic, &
     !      & symprec) bind(c)
     !   import c_int, c_double
     !   integer(c_int), intent(inout) :: rotation(3,3,*)
     !   real(c_double), intent(inout) :: translation(3,*)
     !   integer(c_int), intent(out) :: equivalent_atoms(*)
     !   real(c_double), intent(out) :: primitive_lattice(3,3)
     !   integer(c_int), intent(out) :: spin_flips(*)
     !   integer(c_int), intent(in), value :: max_size
     !   real(c_double), intent(in) :: lattice(3,3), position(3,*)
     !   integer(c_int), intent(in) :: types(*)
     !   real(c_double), intent(in) :: tensors(*)
     !   integer(c_int), intent(in), value :: tensor_rank
     !   integer(c_int), intent(in), value :: num_atom
     !   integer(c_int), intent(in), value :: is_magnetic
     !   real(c_double), intent(in), value :: symprec
     !   integer(c_int) :: spg_get_symmetry_with_site_tensors
     ! end function spg_get_symmetry_with_site_tensors


     ! function spgat_get_symmetry_with_site_tensors( rotation, translation, &
     !      & equivalent_atoms, primitive_lattice, spin_flips, max_size, lattice, &
     !      & position, types, tensors, tensor_rank, num_atom, is_magnetic, &
     !      & symprec, angle_tolerance) bind(c)
     !   import c_int, c_double
     !   integer(c_int), intent(inout) :: rotation(3,3,*)
     !   real(c_double), intent(inout) :: translation(3,*)
     !   integer(c_int), intent(out) :: equivalent_atoms(*)
     !   real(c_double), intent(out) :: primitive_lattice(3,3)
     !   integer(c_int), intent(out) :: spin_flips(*)
     !   integer(c_int), intent(in), value :: max_size
     !   real(c_double), intent(in) :: lattice(3,3), position(3,*)
     !   integer(c_int), intent(in) :: types(*)
     !   real(c_double), intent(in) :: tensors(*)
     !   integer(c_int), intent(in), value :: tensor_rank
     !   integer(c_int), intent(in), value :: num_atom
     !   integer(c_int), intent(in), value :: is_magnetic
     !   real(c_double), intent(in), value :: symprec, angle_tolerance
     !   integer(c_int) :: spgat_get_symmetry_with_site_tensors
     ! end function spgat_get_symmetry_with_site_tensors


     ! function spg_get_multiplicity( lattice, position, types, num_atom, &
     !      & symprec) bind(c)
     !   import c_int, c_double
     !   real(c_double), intent(in) :: lattice(3,3), position(3,*)
     !   integer(c_int), intent(in) :: types(*)
     !   integer(c_int), intent(in), value :: num_atom
     !   real(c_double), intent(in), value :: symprec
     !   integer(c_int) :: spg_get_multiplicity
     ! end function spg_get_multiplicity


     ! function spgat_get_multiplicity( lattice, position, types, num_atom, &
     !      & symprec, angle_tolerance) bind(c)
     !   import c_int, c_double
     !   real(c_double), intent(in) :: lattice(3,3), position(3,*)
     !   integer(c_int), intent(in) :: types(*)
     !   integer(c_int), intent(in), value :: num_atom
     !   real(c_double), intent(in), value :: symprec, angle_tolerance
     !   integer(c_int) :: spgat_get_multiplicity
     ! end function spgat_get_multiplicity


     ! function spg_delaunay_reduce( lattice, symprec) bind(c)
     !   import c_int, c_double
     !   real(c_double), intent(inout) :: lattice(3,3)
     !   real(c_double), intent(in), value :: symprec
     !   integer(c_int) :: spg_delaunay_reduce
     ! end function spg_delaunay_reduce


     ! function spg_niggli_reduce( lattice, symprec) bind(c)
     !   import c_int, c_double
     !   real(c_double), intent(inout) :: lattice(3,3)
     !   real(c_double), intent(in), value :: symprec
     !   integer(c_int) :: spg_niggli_reduce
     ! end function spg_niggli_reduce


     ! function spg_find_primitive(lattice, position, types, num_atom, symprec) &
     !   & bind(c)
     !   import c_int, c_double
     !   real(c_double), intent(inout) :: lattice(3,3), position(3,*)
     !   integer(c_int), intent(inout) :: types(*)
     !   integer(c_int), intent(in), value :: num_atom
     !   real(c_double), intent(in), value :: symprec
     !   integer(c_int) :: spg_find_primitive
     ! end function spg_find_primitive


     ! function spgat_find_primitive(lattice, position, types, num_atom, symprec, &
     !      & angle_tolerance) bind(c)
     !   import c_int, c_double
     !   real(c_double), intent(inout) :: lattice(3,3), position(3,*)
     !   integer(c_int), intent(inout) :: types(*)
     !   integer(c_int), intent(in), value :: num_atom
     !   real(c_double), intent(in), value :: symprec, angle_tolerance
     !   integer(c_int) :: spgat_find_primitive
     ! end function spgat_find_primitive


     ! function spg_get_international( symbol, lattice, position, types, num_atom, &
     !      & symprec) bind(c)
     !   import c_char, c_int, c_double
     !   character(kind=c_char), intent(out) :: symbol(11)
     !   real(c_double), intent(in) :: lattice(3,3), position(3, *)
     !   integer(c_int), intent(in) :: types(*)
     !   integer(c_int), intent(in), value :: num_atom
     !   real(c_double), intent(in), value :: symprec
     !   integer(c_int) :: spg_get_international ! the number corresponding to 'symbol'. 0 on failure
     ! end function spg_get_international


     ! function spgat_get_international( symbol, lattice, position, types, num_atom, &
     !      & symprec, angle_tolerance) bind(c)
     !   import c_char, c_int, c_double
     !   character(kind=c_char), intent(out) :: symbol(11)
     !   real(c_double), intent(in) :: lattice(3,3), position(3, *)
     !   integer(c_int), intent(in) :: types(*)
     !   integer(c_int), intent(in), value :: num_atom
     !   real(c_double), intent(in), value :: symprec, angle_tolerance
     !   integer(c_int) :: spgat_get_international ! the number corresponding to 'symbol'. 0 on failure
     ! end function spgat_get_international


     ! function spg_get_schoenflies( symbol, lattice, position, types, num_atom, &
     !      & symprec) bind(c)
     !   import c_char, c_int, c_double
     !   character(kind=c_char), intent(out) :: symbol(7)
     !   real(c_double), intent(in) :: lattice(3,3), position(3, *)
     !   integer(c_int), intent(in) :: types(*)
     !   integer(c_int), intent(in), value :: num_atom
     !   real(c_double), intent(in), value :: symprec
     !   integer(c_int) :: spg_get_schoenflies ! the number corresponding to 'symbol'. 0 on failure
     ! end function spg_get_schoenflies


     ! function spgat_get_schoenflies( symbol, lattice, position, types, num_atom, &
     !      & symprec, angle_tolerance) bind(c)
     !   import c_char, c_int, c_double
     !   character(kind=c_char), intent(out) :: symbol(7)
     !   real(c_double), intent(in) :: lattice(3,3), position(3, *)
     !   integer(c_int), intent(in) :: types(*)
     !   integer(c_int), intent(in), value :: num_atom
     !   real(c_double), intent(in), value :: symprec, angle_tolerance
     !   integer(c_int) :: spgat_get_schoenflies ! the number corresponding to 'symbol'. 0 on failure
     ! end function spgat_get_schoenflies


     ! function spg_get_pointgroup( symbol, trans_mat, rotations, num_rotations) &
     !      & bind(c)
     !   import c_char, c_int, c_double
     !   character(kind=c_char) :: symbol(6)
     !   integer(c_int), intent(inout) :: trans_mat(3,3)
     !   integer(c_int), intent(in) :: rotations(3,3,*)
     !   integer(c_int), intent(in), value :: num_rotations
     !   integer(c_int) :: spg_get_pointgroup
     ! end function spg_get_pointgroup


     ! function spg_refine_cell( lattice, position, types, num_atom, symprec) bind(c)
     !   import c_int, c_double
     !   real(c_double), intent(inout) :: lattice(3,3), position(3,*)
     !   integer(c_int), intent(inout) :: types(*)
     !   integer(c_int), intent(in), value :: num_atom
     !   real(c_double), intent(in), value :: symprec
     !   integer(c_int) :: spg_refine_cell
     ! end function spg_refine_cell


     ! function spgat_refine_cell( lattice, position, types, num_atom, symprec, &
     !      & angle_tolerance) bind(c)
     !   import c_int, c_double
     !   real(c_double), intent(inout) :: lattice(3,3), position(3,*)
     !   integer(c_int), intent(inout) :: types(*)
     !   integer(c_int), intent(in), value :: num_atom
     !   real(c_double), intent(in), value :: symprec, angle_tolerance
     !   integer(c_int) :: spgat_refine_cell
     ! end function spgat_refine_cell


     ! function spg_standardize_cell( lattice, position, types, num_atom, &
     !      & to_primitive, no_idealize, symprec) bind(c)
     !   import c_int, c_double
     !   real(c_double), intent(inout) :: lattice(3,3), position(3,*)
     !   integer(c_int), intent(inout) :: types(*)
     !   integer(c_int), intent(in), value :: num_atom
     !   integer(c_int), intent(in), value :: to_primitive, no_idealize
     !   real(c_double), intent(in), value :: symprec
     !   integer(c_int) :: spg_refine_cell
     ! end function spg_standardize_cell


     ! function spgat_standardize_cell( lattice, position, types, num_atom, &
     !      & to_primitive, no_idealize, symprec, angle_tolerance) bind(c)
     !   import c_int, c_double
     !   real(c_double), intent(inout) :: lattice(3,3), position(3,*)
     !   integer(c_int), intent(inout) :: types(*)
     !   integer(c_int), intent(in), value :: num_atom
     !   integer(c_int), intent(in), value :: to_primitive, no_idealize
     !   real(c_double), intent(in), value :: symprec, angle_tolerance
     !   integer(c_int) :: spgat_refine_cell
     ! end function spgat_standardize_cell


     ! function spg_get_ir_reciprocal_mesh(grid_point, map, mesh, is_shift, &
     !      & is_time_reversal, lattice, position, types, num_atom, symprec) bind(c)
     !   import c_int, c_double
     !   !   Beware the map refers to positions starting at 0
     !   integer(c_int), intent(out) :: grid_point(3, *), map(*) ! size is product(mesh)
     !   integer(c_int), intent(in) :: mesh(3), is_shift(3)
     !   integer(c_int), intent(in), value :: is_time_reversal
     !   real(c_double), intent(in) :: lattice(3,3), position(3, *)
     !   integer(c_int), intent(in) :: types(*)
     !   integer(c_int), intent(in), value :: num_atom
     !   real(c_double), intent(in), value :: symprec
     !   integer(c_int) :: spg_get_ir_reciprocal_mesh ! the number of points in the reduced mesh
     ! end function spg_get_ir_reciprocal_mesh


     ! function spg_get_stabilized_reciprocal_mesh(grid_point, map, mesh, is_shift, &
     !      & is_time_reversal, lattice, num_rot, rotations, num_q, qpoints) bind(c)
     !   import c_int, c_double
     !   !   Beware the map refers to positions starting at 0
     !   integer(c_int), intent(out) :: grid_point(3,*), map(*)
     !   integer(c_int), intent(in) :: mesh(3)
     !   integer(c_int), intent(in) :: is_shift(3)
     !   integer(c_int), intent(in), value :: is_time_reversal
     !   real(c_double), intent(in) :: lattice(3,3)
     !   integer(c_int), intent(in), value :: num_rot
     !   integer(c_int), intent(in) :: rotations(3,3,*)
     !   integer(c_int), intent(in), value :: num_q
     !   real(c_double), intent(in) :: qpoints(3,*)
     !   integer(c_int) :: spg_get_stabilized_reciprocal_mesh
     ! end function spg_get_stabilized_reciprocal_mesh

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
       integer(c_int), intent(in) :: mesh(3)
       integer(c_int), intent(in) :: address(3)
       integer(c_int), intent(in) :: is_shift(3)
       integer(c_int), intent(inout) :: address_double(3)
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


     function ksp_get_snf3x3(D_diag, P, Q, A) bind(c)
       import c_long, c_int
       integer(c_long), intent(inout) :: D_diag(3)
       integer(c_long), intent(inout) :: P(3, 3)
       integer(c_long), intent(inout) :: Q(3, 3)
       integer(c_long), intent(in) :: A(3, 3)
       integer(c_int) :: ksp_get_snf3x3
     end function ksp_get_snf3x3

  end interface

  public :: ksp_get_major_version,  ksp_get_minor_version, &
       ksp_get_micro_version, ksp_get_all_grid_addresses, ksp_get_snf3x3, &
       ksp_get_double_grid_address, ksp_get_double_grid_index, &
       ksp_get_thm_relative_grid_addresses

end module kspclib_f08

! Copyright (c) 2019-2022, The Neko Authors
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:
!
!   * Redistributions of source code must retain the above copyright
!     notice, this list of conditions and the following disclaimer.
!
!   * Redistributions in binary form must reproduce the above
!     copyright notice, this list of conditions and the following
!     disclaimer in the documentation and/or other materials provided
!     with the distribution.
!
!   * Neither the name of the authors nor the names of its
!     contributors may be used to endorse or promote products derived
!     from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
! COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
! POSSIBILITY OF SUCH DAMAGE.
!
!> NEKTON session data reader
!! @details This module is used to read NEKTON session data in ascii
module rea_file
  use generic_file
  use num_types
  use utils
  use mesh
  use point
  use map
  use rea
  use re2_file, only: re2_file_t
  use map_file
  use comm
  use datadist
  use htable
  use logger
  implicit none
  private

  ! Defines the conventions for conversion of rea labels to labeled zones.
  integer, public :: NEKO_W_BC_LABEL = -1
  integer, public :: NEKO_V_BC_LABEL = -1
  integer, public :: NEKO_O_BC_LABEL = -1
  integer, public :: NEKO_SYM_BC_LABEL = -1
  integer, public :: NEKO_ON_BC_LABEL = -1
  integer, public :: NEKO_SHL_BC_LABEL = -1

  !> Interface for NEKTON ascii files
  type, public, extends(generic_file_t) :: rea_file_t
   contains
     procedure :: read => rea_file_read
     procedure :: write => rea_file_write
  end type rea_file_t

contains

  !> Load NEKTON session data from an ascii file
  subroutine rea_file_read(this, data)
    class(rea_file_t) :: this
    class(*), target, intent(inout) :: data
    type(mesh_t), pointer :: msh
    real(kind=dp), pointer :: params(:)
    character(len=3), pointer :: cbc(:,:)
    integer, allocatable :: curve_type(:,:)
    logical, allocatable :: curve_element(:)
    character(len=1) :: chtemp
    integer :: ndim, nparam, nskip, nlogic, ncurve
    integer :: nelgs, nelgv, i, j, ierr, l
    integer :: el_idx, pt_idx
    logical :: read_param, read_bcs, read_map
    real(kind=dp) :: xc(8), yc(8), zc(8), curve(5)
    real(kind=dp), allocatable :: bc_data(:,:,:)
    real(kind=dp), allocatable :: curve_data(:,:,:)
    type(point_t) :: p(8)
    type(re2_file_t) :: re2_file
    type(map_file_t) :: map_file
    character(len=1024) :: re2_fname, map_fname
    integer :: start_el, end_el, nel, edge
    type(linear_dist_t) :: dist
    type(map_t) :: nm
    type(htable_pt_t) :: htp
    integer :: sym_facet, pids(4), p_el_idx, p_facet
    integer :: off
    integer, parameter, dimension(6) :: facet_map = (/3, 2, 4, 1, 5, 6/)
    logical :: curve_skip = .false.
    character(len=LOG_SIZE) :: log_buf
    ! ---- Offsets for conversion from internal zones to labeled zones
    integer :: user_labeled_zones(NEKO_MSH_MAX_ZLBLS)
    integer :: labeled_zone_offsets(NEKO_MSH_MAX_ZLBLS)
    integer :: total_labeled_zone_offset, current_internal_zone
    total_labeled_zone_offset = 0
    labeled_zone_offsets = 0
    user_labeled_zones = 0
    current_internal_zone = 1
    ! ----

    select type(data)
    type is (rea_t)
       call rea_free(data)
       msh => data%msh
       params => data%params
       cbc => data%cbc
       read_param = .true.
       read_bcs = .true.
    type is (mesh_t)
       msh => data
       read_param = .false.
       read_bcs = .false.
    class default
       call neko_error('Invalid output data')
    end select

    if (read_param .and. read_bcs .and. pe_size .gt. 1) then
       call neko_error('Reading NEKTON session data only implemented in serial')
    end if


    open(unit=9,file=trim(this%fname), status='old', iostat=ierr)
    call neko_log%message('Reading NEKTON file ' // this%fname)

    read(9, *)
    read(9, *)
    read(9, *) ndim
    read(9, *) nparam

    if (.not. read_param) then
       ! Skip parameters
       do i = 1, nparam
          read(9, *)
       end do
    else
       allocate(params(nparam))
       do i = 1, nparam
          read(9, *) params(i)
       end do
    end if

    ! Skip passive scalars
    read(9, *) nskip
    do i = 1, nskip
       read(9, *)
    end do

    ! Skip logic switches
    read(9, *) nlogic
    do i = 1, nlogic
       read(9, *)
    end do

    ! Read mesh info
    read(9, *)
    read(9, *)
    read(9, *) nelgs,ndim, nelgv
    if (nelgs .lt. 0) then
       re2_fname = trim(this%fname(1:scan(trim(this%fname), &
            '.', back=.true.)))//'re2'
       call re2_file%init(re2_fname)
       call re2_file%read(msh)
    else
       write(log_buf,1) ndim, nelgv
1      format('gdim = ', i1, ', nelements =', i7)
       call neko_log%message(log_Buf)

       call filename_chsuffix(this%fname, map_fname, 'map')
       inquire(file=map_fname, exist=read_map)
       if (read_map) then
          call map_init(nm, nelgv, 2**ndim)
          call map_file%init(map_fname)
          call map_file%read(nm)
       else
          call neko_log%warning('No NEKTON map file found')
       end if

       ! Use a load-balanced linear distribution
       dist = linear_dist_t(nelgv, pe_rank, pe_size, NEKO_COMM)
       nel = dist%num_local()
       start_el = dist%start_idx() + 1
       end_el = dist%end_idx() + 1

       call msh%init(ndim, dist)

       call htp%init(2**ndim * nel, ndim)

       el_idx = 1
       pt_idx = 0
       do i = 1, nelgv
          read(9, *)
          if (ndim .eq. 2) then
             read(9, *) (xc(j),j=1,4)
             read(9, *) (yc(j),j=1,4)
             if (i .ge. start_el .and. i .le. end_el) then
                do j = 1, 4
                   p(j) = point_t(real(xc(j),dp), real(yc(j),dp),real(0d0,dp))
                   call rea_file_add_point(htp, p(j), pt_idx)
                end do
                ! swap vertices to keep symmetric vertex numbering in neko
                call msh%add_element(el_idx, p(1), p(2), p(4), p(3))
             end if
          else if (ndim .eq. 3) then
             read(9, *) (xc(j),j=1,4)
             read(9, *) (yc(j),j=1,4)
             read(9, *) (zc(j),j=1,4)
             read(9, *) (xc(j),j=5,8)
             read(9, *) (yc(j),j=5,8)
             read(9, *) (zc(j),j=5,8)
             if (i .ge. start_el .and. i .le. end_el) then
                do j = 1, 8
                   p(j) = point_t(real(xc(j),dp), real(yc(j),dp), real(zc(j),dp))
                   call rea_file_add_point(htp, p(j), pt_idx)
                end do
                ! swap vertices to keep symmetric vertex numbering in neko
                call msh%add_element(el_idx, &
                     p(1), p(2), p(4), p(3), p(5), p(6), p(8), p(7))
             end if
          end if
          if (i .ge. start_el .and. i .le. end_el) then
             el_idx = el_idx + 1
          end if
       end do

       call htp%free()

       read(9, *)
       read(9, *) ncurve
       allocate(curve_data(5,8,nelgv))
       allocate(curve_element(nelgv))
       allocate(curve_type(8,nelgv))
       do i = 1, nelgv
          curve_element(i) = .false.
          do j = 1, 8
             curve_type(j,i) = 0
             do l = 1, 5
                curve_data(l,j,i) = 0d0
             end do
          end do
       end do
       do i = 1, ncurve
          read(9, *) edge, el_idx, (curve(j),j=1,5), chtemp
          do j = 1, 5
             curve_data(j,edge,el_idx) = curve(j)
          end do
          curve_element(el_idx) = .true.
          select case(trim(chtemp))
          case ('s')
             curve_type(edge,el_idx) = 1
             curve_skip = .true.
          case ('e')
             curve_type(edge,el_idx) = 2
             curve_skip = .true.
          case ('C')
             curve_type(edge,el_idx) = 3
          case ('m')
             curve_type(edge,el_idx) = 4
          end select
       end do
       if (curve_skip) then
          call neko_log%warning('Curve type: s, e are not supported, treating mesh as non-curved.')
       else
          do el_idx = 1, nelgv
             if (curve_element(el_idx)) then
                call msh%mark_curve_element(el_idx, &
                     curve_data(1,1,el_idx), curve_type(1,el_idx))
             end if
          end do
       end if
       deallocate(curve_data)
       deallocate(curve_element)
       deallocate(curve_type)

       ! Read fluid boundary conditions
       read(9,*)
       read(9,*)
       if (.not. read_bcs) then ! Mark zones in the mesh
          call neko_log%message("Reading boundary conditions", neko_log_debug)
          allocate(cbc(6,nelgv))
          allocate(bc_data(6,2*ndim,nelgv))
          off = 0
          !Fix for different horrible .rea periodic bc formats.
          if (nelgv .lt. 1000) off = 1
          do i = 1, nelgv
             if (i .ge. start_el .and. i .le. end_el) then
                el_idx = i - start_el + 1
                do j = 1, 2*ndim
                   read(9, *) cbc(j, i), (bc_data(l,j,i),l=1,6)
                   sym_facet = facet_map(j)

                   select case(trim(cbc(j,i)))
                   case ('W')
                      if (NEKO_W_BC_LABEL .eq. -1) then
                         NEKO_W_BC_LABEL = current_internal_zone
                         current_internal_zone = current_internal_zone + 1
                      end if

                      call rea_file_mark_labeled_bc(msh, el_idx, sym_facet, &
                           cbc(j,i), &
                           NEKO_W_BC_LABEL, total_labeled_zone_offset, &
                           labeled_zone_offsets(NEKO_W_BC_LABEL) .eq. 0)

                      labeled_zone_offsets(NEKO_W_BC_LABEL) = 1

                   case ('v', 'V')
                      if (NEKO_V_BC_LABEL .eq. -1) then
                         NEKO_V_BC_LABEL = current_internal_zone
                         current_internal_zone = current_internal_zone + 1
                      end if

                      call rea_file_mark_labeled_bc(msh, el_idx, sym_facet, &
                           cbc(j,i), &
                           NEKO_V_BC_LABEL, total_labeled_zone_offset, &
                           labeled_zone_offsets(NEKO_V_BC_LABEL) .eq. 0)

                      labeled_zone_offsets(NEKO_V_BC_LABEL) = 1
                   case ('O', 'o')
                      if (NEKO_O_BC_LABEL .eq. -1) then
                         NEKO_O_BC_LABEL = current_internal_zone
                         current_internal_zone = current_internal_zone + 1
                      end if

                      call rea_file_mark_labeled_bc(msh, el_idx, sym_facet, &
                           cbc(j,i), &
                           NEKO_O_BC_LABEL, total_labeled_zone_offset, &
                           labeled_zone_offsets(NEKO_O_BC_LABEL) .eq. 0)

                      labeled_zone_offsets(NEKO_O_BC_LABEL) = 1
                   case ('SYM')
                      if (NEKO_SYM_BC_LABEL .eq. -1) then
                         NEKO_SYM_BC_LABEL = current_internal_zone
                         current_internal_zone = current_internal_zone + 1
                      end if

                      call rea_file_mark_labeled_bc(msh, el_idx, sym_facet, &
                           cbc(j,i), &
                           NEKO_SYM_BC_LABEL, total_labeled_zone_offset, &
                           labeled_zone_offsets(NEKO_SYM_BC_LABEL) .eq. 0)

                      labeled_zone_offsets(NEKO_SYM_BC_LABEL) = 1
                   case ('ON', 'on')
                      if (NEKO_ON_BC_LABEL .eq. -1) then
                         NEKO_ON_BC_LABEL = current_internal_zone
                         current_internal_zone = current_internal_zone + 1
                      end if

                      call rea_file_mark_labeled_bc(msh, el_idx, sym_facet, &
                           cbc(j,i), &
                           NEKO_ON_BC_LABEL, total_labeled_zone_offset, &
                           labeled_zone_offsets(NEKO_ON_BC_LABEL) .eq. 0)

                      labeled_zone_offsets(NEKO_ON_BC_LABEL) = 1
                   case ('s', 'sl', 'sh', 'shl', 'S', 'SL', 'SH', 'SHL')
                      if (NEKO_SHL_BC_LABEL .eq. -1) then
                         NEKO_SHL_BC_LABEL = current_internal_zone
                         current_internal_zone = current_internal_zone + 1
                      end if

                      call rea_file_mark_labeled_bc(msh, el_idx, sym_facet, &
                           cbc(j,i), &
                           NEKO_SHL_BC_LABEL, total_labeled_zone_offset, &
                           labeled_zone_offsets(NEKO_SHL_BC_LABEL) .eq. 0)

                      labeled_zone_offsets(NEKO_SHL_BC_LABEL) = 1
                   case ('P')
                      p_el_idx = int(bc_data(2+off,j,i))
                      p_facet = facet_map(int(bc_data(3+off,j,i)))
                      call msh%get_facet_ids(sym_facet, el_idx, pids)
                      call msh%mark_periodic_facet(sym_facet, el_idx, &
                           p_facet, p_el_idx, pids)
                   end select
                end do
             end if
          end do

          do i = 1, nelgv
             if (i .ge. start_el .and. i .le. end_el) then
                el_idx = i - start_el + 1
                do j = 1, 2*ndim
                   sym_facet = facet_map(j)
                   select case(trim(cbc(j,i)))
                   case ('P')
                      p_el_idx = int(bc_data(2+off,j,i))
                      p_facet = facet_map(int(bc_data(3+off,j,i)))
                      call msh%create_periodic_ids(sym_facet, el_idx, &
                           p_facet, p_el_idx)
                   end select
                end do
             end if
          end do
          do i = 1, nelgv
             if (i .ge. start_el .and. i .le. end_el) then
                el_idx = i - start_el + 1
                do j = 1, 2*ndim
                   sym_facet = facet_map(j)
                   select case(trim(cbc(j,i)))
                   case ('P')
                      p_el_idx = int(bc_data(2+off,j,i))
                      p_facet = facet_map(int(bc_data(3+off,j,i)))
                      call msh%create_periodic_ids(sym_facet, el_idx, &
                           p_facet, p_el_idx)
                   end select
                end do
             end if
          end do
          do i = 1, nelgv
             if (i .ge. start_el .and. i .le. end_el) then
                el_idx = i - start_el + 1
                do j = 1, 2*ndim
                   sym_facet = facet_map(j)
                   select case(trim(cbc(j,i)))
                   case ('P')
                      p_el_idx = int(bc_data(2+off,j,i))
                      p_facet = facet_map(int(bc_data(3+off,j,i)))
                      call msh%create_periodic_ids(sym_facet, el_idx, &
                           p_facet, p_el_idx)
                   end select
                end do
             end if
          end do
          deallocate(cbc)
          deallocate(bc_data)
       else  ! Store bcs in a NEKTON session structure
          allocate(cbc(6,nelgv))
          do i = 1, nelgv
             do j = 1, 2*ndim
                read(9,'(a1, a3)') chtemp, cbc(j, i)
             end do
          end do
       end if

       call msh%finalize()

       call neko_log%message('Done')
       close(9)
    endif

  end subroutine rea_file_read

  subroutine rea_file_write(this, data, t)
    class(rea_file_t), intent(inout) :: this
    class(*), target, intent(in) :: data
    real(kind=rp), intent(in), optional :: t
  end subroutine rea_file_write

  subroutine rea_file_add_point(htp, p, idx)
    type(htable_pt_t), intent(inout) :: htp
    type(point_t), intent(inout) :: p
    integer, intent(inout) :: idx
    integer :: tmp

    if (htp%get(p, tmp) .gt. 0) then
       idx = idx + 1
       call htp%set(p, idx)
       call p%set_id(idx)
    else
       call p%set_id(tmp)
    end if

  end subroutine rea_file_add_point

  !> Mark a labeled zone based on an offset
  !! @param msh The mesh on which to mark the labeled zone.
  !! @param el_idx The index of the element on which to mark the labeled zone.
  !! @param facet The facet index to mark.
  !! @param type The rea label type (e.g. "W", "ON", "SYM", etc).
  !! @param label The integer label with which to mark the labeled zone.
  !! @param offset The offset with which to increment the label, in case there
  !! are any existing user labeled zones.
  !! @param print_info Wether or not to print information to the standard
  !! output.
  subroutine rea_file_mark_labeled_bc(msh, el_idx, facet, type, label, offset, print_info)
    type(mesh_t), intent(inout) :: msh
    integer, intent(in) :: el_idx
    integer, intent(in) :: facet
    character(len=*), intent(in) :: type
    integer, intent(in) :: label
    integer, intent(in) :: offset
    logical, intent(in) :: print_info

    integer :: mark_label
    character(len=LOG_SIZE) :: log_buf

    mark_label = offset + label

    if (mark_label .lt. 1 .or. mark_label .gt. NEKO_MSH_MAX_ZLBLS) then
       call neko_error("You have reached the maximum amount of allowed labeled&
& zones (max allowed: 20). This happened when converting re2 internal labels&
& like e.g. 'w', 'V' or 'o' to labeled zones. Please reduce the number of&
& labeled zones that you have defined or make sure that they are labeled&
& from [1,...,20].")
    end if

    if (print_info) then
       write (log_buf, "(A3,A,I2)") trim(type), " => Labeled index ", mark_label
       call neko_log%message(log_buf)
    end if
    call msh%mark_labeled_facet(facet, el_idx, mark_label)

  end subroutine rea_file_mark_labeled_bc

end module rea_file

module vars
  integer, parameter :: nlayer = 200
  real(8), parameter :: thick = 0.5d0 ! Ang, thickness of one atomic layer
  integer, save :: atoms, nbulk_Mo, nedge_Mo, nbulk_S, nedge_S, fixed(3)
  integer, save :: atoms_layer(nlayer)
  real(8), save :: side(3), thickness_edge
  integer, allocatable, save :: tipo(:), layer(:)
  real(8), allocatable, save :: mass(:), pos_eq(:,:), poteng_eq(:), distance(:,:)
  real(8), allocatable, save :: pos(:,:), vel(:,:), force(:,:), poteng(:), kineng(:,:)
  character(6), allocatable, save :: tag(:)
end module vars

program all
  use vars, only: atoms, mass, tag, pos_eq, fixed, pos, vel, force, poteng, poteng_eq, &
       kineng, distance, layer, nbulk_Mo, nbulk_S, nedge_Mo, nedge_S, thickness_edge
  implicit none
  real(8), parameter :: metal_to_eV = 1.d-3 / 6.02214076d0 / 1.602176634d0
  real(8), parameter :: tol = 1.d-3
  integer :: i, j, id, iconf, timestep, timestep0, ny, iperiod, min
  real(8) :: dt, pos_tmp(3), poteng_tmp, vel_tmp(3), force_tmp(3)
  real(8) :: err, displ(3), time, ly, velocity, displacement, offset, period, position, rmsd_bulk, rmsd_edge
  character(20) :: exclude

  read(*,*) velocity, dt, min, displacement, ly, ny, exclude, thickness_edge ! velocity in Ang/ps, dt in ps

  ! set parameters
  period = ly / dble(ny) ! Ang
  offset = dble(min) * displacement ! Ang
  if ( exclude == 'first' ) offset = offset + period
  
  call read_eq

  call read_lists

  allocate ( pos(atoms,3), stat = i )
  if ( i /= 0 ) stop 'Allocation POS'
  allocate ( vel(atoms,3), stat = i )
  if ( i /= 0 ) stop 'Allocation VEL'
  allocate ( force(atoms,3), stat = i )
  if ( i /= 0 ) stop 'Allocation FORCE'
  allocate ( poteng(atoms), stat = i )
  if ( i /= 0 ) stop 'Allocation PE'
  allocate ( kineng(atoms,3), stat = i )
  if ( i /= 0 ) stop 'Allocation KE'
  
  open ( unit=100, file='system.lammpstrj' )
  open ( unit=200, file='3-displacements.dat' )
  open ( unit=300, file='3-velocities.dat' )
  open ( unit=400, file='3-forces.dat' )
  open ( unit=500, file='3-potential_energies.dat' )
  open ( unit=501, file='3-potential_energies_abs.dat' )
  open ( unit=99, file='3-potential_energies_eq.dat' )
  do i = 1, atoms
     write(99,*) 0.d0, layer(i), minval(distance(i,:)), poteng_eq(i), i
  end do
  close(99)   
  open ( unit=600, file='3-kinetic_energies.dat' )
  open ( unit=700, file='3-rmsd.dat' )
  iconf = 0
  iperiod = 0
  do
     ! read
     read(100,*,iostat=i) ! ITEM: TIMESTEP
     if ( i /= 0 ) exit
     read(100,*) timestep
     if ( iperiod == 0 ) then
        timestep0 = timestep
        iperiod = 1
     end if
     time = dble(timestep-timestep0) * dt ! time in ps
     position = velocity * time
     if ( position < offset ) then ! starting from the first minimum
        do i = 1, atoms+7
           read(100,*)
        end do
        cycle
     else
        if ( position >= dble(iperiod) * period + offset ) then
           iperiod = iperiod + 1
           if ( dble(iperiod) * period + offset > 20.d0 ) exit
        end if
     end if
     iconf = iconf + 1
     if ( mod(iconf,100) == 0 ) write(0,*) iconf
     write(0,*) iconf
     do i = 1, 7
        read(100,*)
     end do
     do i = 1, atoms
        read(100,*) id, (pos_tmp(j),j=1,3), (vel_tmp(j),j=1,3), (force_tmp(j),j=1,3), poteng_tmp
        pos(id,:) = pos_tmp(:)
        vel(id,:) = vel_tmp(:)
        force(id,:) = force_tmp(:)
        poteng(id) = poteng_tmp
        kineng(id,:) = 0.5d0 * mass(id) * vel_tmp(:) * vel_tmp(:) * metal_to_eV
     end do

     ! shift
     do j = 1, 3
        pos(:,j) = pos(:,j) - ( pos(fixed(1),j) - pos_eq(fixed(1),j) )
     end do

     ! check   
     do i = 2, 3
        err = sqrt( dot_product( pos(fixed(i),:)-pos_eq(fixed(i),:), pos(fixed(i),:)-pos_eq(fixed(i),:) ) )
        if ( err > tol ) write(0,*) 'CAREFUL! ', i, err, (j, pos(fixed(i),j), pos_eq(fixed(i),j),j=1,3)
     end do

     ! printout
     rmsd_edge = 0.d0
     rmsd_bulk = 0.d0
     do i = 1, atoms
        if ( tag(i) == 'bulk' .or. tag(i) == 'edge' ) then
           ! displacements
           displ(:) = pos(i,:) - pos_eq(i,:)
           write(200,*) position, layer(i), minval(distance(i,:)), (displ(j),j=1,3), sqrt(dot_product(displ,displ)), i
           ! velocities
           write(300,*) position, layer(i), minval(distance(i,:)), (vel(i,j),j=1,3), i
           ! forces
           write(400,*) position, layer(i), minval(distance(i,:)), (force(i,j),j=1,3), sqrt(dot_product(force(i,:),force(i,:))), i
           ! potential energies
           write(500,*) position, layer(i), minval(distance(i,:)), poteng(i)-poteng_eq(i), i
           write(501,*) position, layer(i), minval(distance(i,:)), poteng(i), i
           ! kinetic energies
           write(600,*) position, layer(i), minval(distance(i,:)), (kineng(i,j),j=1,3), sum(kineng(i,:)), i
           if ( tag(i) == 'edge' ) then
              rmsd_edge = rmsd_edge + dot_product(displ,displ)
           else if ( tag(i) == 'bulk' ) then
              rmsd_bulk = rmsd_bulk + dot_product(displ,displ)
           end if
        end if
     end do
     rmsd_edge = sqrt( rmsd_edge / dble(nedge_Mo+nedge_S) )
     if ( nbulk_Mo+nbulk_S > 0 ) rmsd_bulk = sqrt( rmsd_bulk / dble(nbulk_Mo+nbulk_S) )
     write(700,*) position, rmsd_edge, rmsd_bulk
     
  end do
  close(100)
  close(200)
  close(300)
  close(400)
  close(500)
  close(600)
  close(700)
  
  stop
end program all

subroutine read_eq
  use vars, only: atoms, side, tipo, mass, pos_eq, distance, nlayer, &
       atoms_layer, layer, atoms, poteng_eq, thick
  implicit none
  integer :: i, j, id, tipo_tmp, atom_types, ix(3)
  integer :: iA, iB, iC, l, k, kk
  real(8) :: xlo, xhi, m, pos_tmp(3), xxx
  real(8) :: xmin, xmax, ymax, d_min, poteng_tmp!, y1, y2
  real(8), allocatable :: mass_pertype(:)
  real(8) :: slope(3), intercept(3)

  ! positions
  open(unit=10,file='equ.data')
  read(10,*)
  read(10,*)
  read(10,*) atoms
  read(10,*) atom_types
  read(10,*)
  do i = 1, 3
     read(10,*) xlo, xhi
     side(i) = xhi - xlo
  end do
  read(10,*)
  read(10,*) ! Masses
  read(10,*)
  allocate ( mass_pertype(atom_types), stat = i )
  if ( i /= 0 ) stop 'Allocation MASS_TMP'
  do i = 1, atom_types
     read(10,*) j, m
     mass_pertype(j) = m
  end do
  read(10,*)
  read(10,*) ! Atoms # atomic
  read(10,*)
  allocate ( tipo(atoms), stat = i )
  if ( i /= 0 ) stop 'Allocation TIPO'
  allocate ( pos_eq(atoms,3), stat = i )
  if ( i /= 0 ) stop 'Allocation POS_EQ'
  allocate ( mass(atoms), stat = i )
  if ( i /= 0 ) stop 'Allocation MASS'
  do i = 1, atoms
     read(10,*) id, tipo_tmp, (pos_tmp(j),j=1,3), (ix(j),j=1,3)
     mass(id) = mass_pertype(tipo_tmp)
     tipo(id) = tipo_tmp
     pos_eq(id,:) = pos_tmp(:) + dble(ix(:)) * side(:)
  end do
  close(10)
  deallocate ( mass_pertype )

  !    C
  !   / \
  !  A---B
  ! find vertexes
  xmin = 1.d99
  xmax = -1.d99
  ymax = -1.d99
  do i = 1, atoms
     if ( tipo(i) < 4 ) then
        if ( pos_eq(i,1) < xmin ) then
           xmin = pos_eq(i,1)
           iA = i
        end if
        if ( pos_eq(i,1) > xmax ) then
           xmax = pos_eq(i,1)
           iB = i
        end if
        if ( pos_eq(i,2) > ymax ) then
           ymax = pos_eq(i,2)
           iC = i
        end if
     end if
  end do

  ! equations of segments
  ! 1) A-B, 2) A-C, 3) B-C
  slope(1) = ( pos_eq(iB,2) - pos_eq(iA,2) ) / ( pos_eq(iB,1) - pos_eq(iA,1) )
  slope(2) = ( pos_eq(iC,2) - pos_eq(iA,2) ) / ( pos_eq(iC,1) - pos_eq(iA,1) )
  slope(3) = ( pos_eq(iC,2) - pos_eq(iB,2) ) / ( pos_eq(iC,1) - pos_eq(iB,1) )
  intercept(1) = pos_eq(iA,2) - slope(1) * pos_eq(iA,1)
  intercept(2) = pos_eq(iA,2) - slope(2) * pos_eq(iA,1)
  intercept(3) = pos_eq(iB,2) - slope(3) * pos_eq(iB,1)

  ! calculate distance from edge
  allocate ( distance(atoms,3), stat = i )
  if ( i /= 0 ) stop 'Allocation DISTANCE'
  do i = 1, atoms
     if ( tipo(i) < 4 ) then
        do j = 1, 3
           distance(i,j) = abs(pos_eq(i,2)-(slope(j)*pos_eq(i,1)+intercept(j))) / sqrt(1.d0+slope(j)**2)
        end do
     end if
  end do

  ! find layers
  allocate ( layer(atoms), stat = i )
  if ( i /= 0 ) stop 'Allocation LAYER'
  layer(:) = 0
  atoms_layer(:) = 0
  do l = 1, nlayer
     d_min = 1.d99
     do i = 1, atoms
        if ( tipo(i) < 4 ) then
           if ( layer(i) == 0 ) then
              if ( minval(distance(i,:)) < d_min ) d_min = minval(distance(i,:))
           end if
        end if
     end do
     do i = 1, atoms
        if ( tipo(i) < 4 ) then
           if ( layer(i) == 0 ) then
              if ( minval(distance(i,:)) < d_min+thick ) then
                 atoms_layer(l) = atoms_layer(l) + 1
                 layer(i) = l
              end if
           end if
        end if
     end do
  end do

  ! potential energy
  allocate ( poteng_eq(atoms), stat = i )
  if ( i /= 0 ) stop 'Allocation POTENG_EQ'
  open(unit=99,file='equ.lammpstrj')
  do i = 1, 9
     read(99,*)
  end do
  do i = 1, atoms
     read(99,*) id, (xxx,j=1,9), poteng_tmp
     poteng_eq(id) = poteng_tmp
  end do
  close(99)

kk = 0  
! ~ -0.8 eV/atom 
k = 0  
do i = 1, atoms
   if ( tipo(i) == 1 .or. tipo(i) == 2 .or. tipo(i) == 3 ) then
      if ( poteng_eq(i) < -0.7d0 .and. poteng_eq(i) > -0.9d0 ) k = k + 1
   end if
end do
kk = kk + k
open(unit=666,file='ref_-0.8.xyz')
write(666,*) k
write(666,*)
do i = 1, atoms
   if ( tipo(i) == 1 .or. tipo(i) == 2 .or. tipo(i) == 3 ) then
      if ( poteng_eq(i) < -0.7d0 .and. poteng_eq(i) > -0.9d0 ) then
         if ( tipo(i) == 2 ) write(666,*) 'Mo', (pos_eq(i,j),j=1,3)
         if ( tipo(i) == 1 .or. tipo(i) == 3 ) write(666,*) 'S', (pos_eq(i,j),j=1,3)
      end if
   end if
end do
close(666)

! ~ -1.5 eV/atom 
k = 0  
do i = 1, atoms
   if ( tipo(i) == 1 .or. tipo(i) == 2 .or. tipo(i) == 3 ) then
      if ( poteng_eq(i) < -1.5d0 .and. poteng_eq(i) > -1.6d0 ) k = k + 1
   end if
end do
kk = kk + k
open(unit=666,file='ref_-1.5.xyz')
write(666,*) k
write(666,*)
do i = 1, atoms
   if ( tipo(i) == 1 .or. tipo(i) == 2 .or. tipo(i) == 3 ) then
      if ( poteng_eq(i) < -1.5d0 .and. poteng_eq(i) > -1.6d0 ) then
         if ( tipo(i) == 2 ) write(666,*) 'Mo', (pos_eq(i,j),j=1,3)
         if ( tipo(i) == 1 .or. tipo(i) == 3 ) write(666,*) 'S', (pos_eq(i,j),j=1,3)
      end if
   end if
end do
close(666)

! ~ -2.3 eV/atom 
k = 0  
do i = 1, atoms
   if ( tipo(i) == 1 .or. tipo(i) == 2 .or. tipo(i) == 3 ) then
      if ( poteng_eq(i) < -2.2d0 .and. poteng_eq(i) > -2.4d0 ) k = k + 1
   end if
end do
kk = kk + k
open(unit=666,file='ref_-2.3.xyz')
write(666,*) k
write(666,*)
do i = 1, atoms
   if ( tipo(i) == 1 .or. tipo(i) == 2 .or. tipo(i) == 3 ) then
      if ( poteng_eq(i) < -2.2d0 .and. poteng_eq(i) > -2.4d0 ) then
         if ( tipo(i) == 2 ) write(666,*) 'Mo', (pos_eq(i,j),j=1,3)
         if ( tipo(i) == 1 .or. tipo(i) == 3 ) write(666,*) 'S', (pos_eq(i,j),j=1,3)
      end if
   end if
end do
close(666)

! ~ -5.2 eV/atom 
k = 0  
do i = 1, atoms
   if ( tipo(i) == 1 .or. tipo(i) == 2 .or. tipo(i) == 3 ) then
      if ( poteng_eq(i) < -5.1d0 .and. poteng_eq(i) > -5.3d0 ) k = k + 1
   end if
end do
kk = kk + k
open(unit=666,file='ref_-5.2.xyz')
write(666,*) k
write(666,*)
do i = 1, atoms
   if ( tipo(i) == 1 .or. tipo(i) == 2 .or. tipo(i) == 3 ) then
      if ( poteng_eq(i) < -5.1d0 .and. poteng_eq(i) > -5.3d0 ) then
         if ( tipo(i) == 2 ) write(666,*) 'Mo', (pos_eq(i,j),j=1,3)
         if ( tipo(i) == 1 .or. tipo(i) == 3 ) write(666,*) 'S', (pos_eq(i,j),j=1,3)
      end if
   end if
end do
close(666)
  
! ~ -6.7 eV/atom 
k = 0  
do i = 1, atoms
   if ( tipo(i) == 1 .or. tipo(i) == 2 .or. tipo(i) == 3 ) then
      if ( poteng_eq(i) < -6.7d0 .and. poteng_eq(i) > -6.8d0 ) k = k + 1
   end if
end do
kk = kk + k
open(unit=666,file='ref_-6.7.xyz')
write(666,*) k
write(666,*)
do i = 1, atoms
   if ( tipo(i) == 1 .or. tipo(i) == 2 .or. tipo(i) == 3 ) then
      if ( poteng_eq(i) < -6.7d0 .and. poteng_eq(i) > -6.8d0 ) then
         if ( tipo(i) == 2 ) write(666,*) 'Mo', (pos_eq(i,j),j=1,3)
         if ( tipo(i) == 1 .or. tipo(i) == 3 ) write(666,*) 'S', (pos_eq(i,j),j=1,3)
      end if
   end if
end do
close(666)

! ~ -7.0 eV/atom ! graphene
!k = 0  
!do i = 1, atoms
!   if ( tipo(i) == 1 .or. tipo(i) == 2 .or. tipo(i) == 3 ) then
!      if ( poteng_eq(i) < -7.0d0 .and. poteng_eq(i) > -7.1d0 ) k = k + 1
!   end if
!end do
!kk = kk + k
!open(unit=666,file='ref_-7.0.xyz')
!write(666,*) k
!write(666,*)
!do i = 1, atoms
!   if ( tipo(i) == 1 .or. tipo(i) == 2 .or. tipo(i) == 3 ) then
!      if ( poteng_eq(i) < -7.0d0 .and. poteng_eq(i) > -7.1d0 ) then
!         if ( tipo(i) == 2 ) write(666,*) 'Mo', (pos_eq(i,j),j=1,3)
!         if ( tipo(i) == 1 .or. tipo(i) == 3 ) write(666,*) 'S', (pos_eq(i,j),j=1,3)
!      end if
!   end if
!end do
!close(666)

! ~ -8.2 eV/atom 
k = 0  
do i = 1, atoms
   if ( tipo(i) == 1 .or. tipo(i) == 2 .or. tipo(i) == 3 ) then
      if ( poteng_eq(i) < -8.2d0 .and. poteng_eq(i) > -8.3d0 ) k = k + 1
   end if
end do
kk = kk + k
open(unit=666,file='ref_-8.2.xyz')
write(666,*) k
write(666,*)
do i = 1, atoms
   if ( tipo(i) == 1 .or. tipo(i) == 2 .or. tipo(i) == 3 ) then
      if ( poteng_eq(i) < -8.2d0 .and. poteng_eq(i) > -8.3d0 ) then
         if ( tipo(i) == 2 ) write(666,*) 'Mo', (pos_eq(i,j),j=1,3)
         if ( tipo(i) == 1 .or. tipo(i) == 3 ) write(666,*) 'S', (pos_eq(i,j),j=1,3)
      end if
   end if
end do
close(666)

k = 0
do i = 1, atoms
   if ( tipo(i) == 1 .or. tipo(i) == 2 .or. tipo(i) == 3 ) k = k + 1
end do
if ( k /= kk ) then
   write(0,*) 'tot', k
   write(0,*) 'xyz', kk
   stop 'aiudo!'
end if

  return
end subroutine read_eq

subroutine read_lists
  use vars, only: atoms, tipo, tag, nbulk_Mo, nedge_Mo, nbulk_S, nedge_S, fixed, pos_eq, distance, thickness_edge
  implicit none
  integer :: i, id, iA, iB, iC
  real(8) :: xmax, xmin, ymax
  
  allocate ( tag(atoms), stat = i )
  if ( i /= 0 ) stop 'Allocation TAG'

  tag(:) = 'none'
  
!  open(unit=20,file='bulk.txt')
!  do
!     read(20,*,iostat=i) id
!     if ( i /= 0 ) exit
!     tag(id) = 'bulk'
!  end do
!  close(20)
!  
!  open(unit=21,file='edge.txt')
!  do
!     read(21,*,iostat=i) id
!     if ( i /= 0 ) exit
!     tag(id) = 'edge'
!  end do
!  close(21)

  ! using new definition of bulk/edge
  do i = 1, atoms
     if ( tipo(i) < 4 ) then
        if ( minval(distance(i,:)) <= thickness_edge ) then
           tag(i) = 'edge'
        else
           tag(i) = 'bulk'
        end if
     end if
  end do
  
  open(unit=22,file='fixed_sulfurs.txt')
  read(22,*) (fixed(i),i=1,3)
  close(22)
  do i = 1, 3
     tag(fixed(i)) = 'fix'
  end do
  
  !    C
  !   / \
  !  A---B
  ! find corners
  xmin = 1.d99
  xmax = -1.d99
  ymax = -1.d99
  do i = 1, atoms
     if ( tipo(i) == 3 ) then
        if ( pos_eq(i,1) < xmin ) then
           xmin = pos_eq(i,1)
           iA = i
        end if
        if ( pos_eq(i,1) > xmax ) then
           xmax = pos_eq(i,1)
           iB = i
        end if
        if ( pos_eq(i,2) > ymax ) then
           ymax = pos_eq(i,2)
           iC = i
        end if
     end if
  end do
  tag(iA) = 'corner'
  tag(iB) = 'corner'
  tag(iC) = 'corner'
  xmin = 1.d99
  xmax = -1.d99
  ymax = -1.d99
  do i = 1, atoms
     if ( tipo(i) == 1 ) then
        if ( pos_eq(i,1) < xmin ) then
           xmin = pos_eq(i,1)
           iA = i
        end if
        if ( pos_eq(i,1) > xmax ) then
           xmax = pos_eq(i,1)
           iB = i
        end if
        if ( pos_eq(i,2) > ymax ) then
           ymax = pos_eq(i,2)
           iC = i
        end if
     end if
  end do
  tag(iA) = 'corner'
  tag(iB) = 'corner'
  tag(iC) = 'corner'

  ! count atoms
  nbulk_Mo = 0
  nbulk_S = 0
  nedge_Mo = 0
  nedge_S = 0
  do i = 1, atoms
     if ( tag(i) == 'bulk' ) then
        if ( tipo(i) == 2 ) nbulk_Mo = nbulk_Mo + 1
        if ( tipo(i) == 1 .or. tipo(i) == 3 ) nbulk_S = nbulk_S + 1
     else if ( tag(i) == 'edge' ) then
        if ( tipo(i) == 2 ) nedge_Mo = nedge_Mo + 1
        if ( tipo(i) == 1 .or. tipo(i) == 3 ) nedge_S = nedge_S + 1
     end if
  end do

  return
end subroutine read_lists


!open(unit=1000,file='x.lammpstrj')
!write(1000,'(a)') 'ITEM: TIMESTEP'
!write(1000,*) 0
!write(1000,'(a)') 'ITEM: NUMBER OF ATOMS'
!write(1000,*) sum( atoms_layer(:) )
!write(1000,'(a)') 'ITEM: BOX BOUNDS pp pp pp'
!write(1000,*) 0.d0, side(1)
!write(1000,*) 0.d0, side(2)
!write(1000,*) 0.d0, side(3)
!write(1000,'(a)') 'ITEM: ATOMS element xu yu zu'
!do i = 1, atoms
!   if ( layer(i) == 1 ) write(1000,*) 'H', (pos_eq(i,j),j=1,3)
!   if ( layer(i) == 2 ) write(1000,*) 'He', (pos_eq(i,j),j=1,3)
!   if ( layer(i) == 3 ) write(1000,*) 'Li', (pos_eq(i,j),j=1,3)
!   if ( layer(i) == 4 ) write(1000,*) 'Be', (pos_eq(i,j),j=1,3)
!   if ( layer(i) == 5 ) write(1000,*) 'B', (pos_eq(i,j),j=1,3)
!   if ( layer(i) == 6 ) write(1000,*) 'C', (pos_eq(i,j),j=1,3)
!   if ( layer(i) == 7 ) write(1000,*) 'N', (pos_eq(i,j),j=1,3)
!   if ( layer(i) == 8 ) write(1000,*) 'O', (pos_eq(i,j),j=1,3)
!   if ( layer(i) == 9 ) write(1000,*) 'F', (pos_eq(i,j),j=1,3)
!   if ( layer(i) == 10 ) write(1000,*) 'Ne', (pos_eq(i,j),j=1,3)
!   if ( layer(i) == 11 ) write(1000,*) 'Na', (pos_eq(i,j),j=1,3)
!   if ( layer(i) == 12 ) write(1000,*) 'Mg', (pos_eq(i,j),j=1,3)
!   if ( layer(i) == 13 ) write(1000,*) 'Al', (pos_eq(i,j),j=1,3)
!   if ( layer(i) == 14 ) write(1000,*) 'Si', (pos_eq(i,j),j=1,3)
!   if ( layer(i) == 15 ) write(1000,*) 'P', (pos_eq(i,j),j=1,3)
!   if ( layer(i) == 16 ) write(1000,*) 'S', (pos_eq(i,j),j=1,3)
!   if ( layer(i) == 17 ) write(1000,*) 'Cl', (pos_eq(i,j),j=1,3)
!   if ( layer(i) == 18 ) write(1000,*) 'Ar', (pos_eq(i,j),j=1,3)
!   if ( layer(i) == 19 ) write(1000,*) 'K', (pos_eq(i,j),j=1,3)
!   if ( layer(i) == 20 ) write(1000,*) 'Ca', (pos_eq(i,j),j=1,3)
!   if ( layer(i) == 21 ) write(1000,*) 'Sc', (pos_eq(i,j),j=1,3)
!   if ( layer(i) == 22 ) write(1000,*) 'Ti', (pos_eq(i,j),j=1,3)
!   if ( layer(i) == 23 ) write(1000,*) 'V', (pos_eq(i,j),j=1,3)
!   if ( layer(i) == 24 ) write(1000,*) 'Cr', (pos_eq(i,j),j=1,3)
!   if ( layer(i) == 25 ) write(1000,*) 'Mn', (pos_eq(i,j),j=1,3)
!   if ( layer(i) == 26 ) write(1000,*) 'Fe', (pos_eq(i,j),j=1,3)
!   if ( layer(i) == 27 ) write(1000,*) 'Co', (pos_eq(i,j),j=1,3)
!   if ( layer(i) == 28 ) write(1000,*) 'Ni', (pos_eq(i,j),j=1,3)
!   if ( layer(i) == 29 ) write(1000,*) 'Cu', (pos_eq(i,j),j=1,3)
!   if ( layer(i) == 30 ) write(1000,*) 'Zn', (pos_eq(i,j),j=1,3)
!   if ( layer(i) == 31 ) write(1000,*) 'Ga', (pos_eq(i,j),j=1,3)
!   if ( layer(i) == 32 ) write(1000,*) 'Ge', (pos_eq(i,j),j=1,3)
!   if ( layer(i) == 33 ) write(1000,*) 'As', (pos_eq(i,j),j=1,3)
!   if ( layer(i) == 34 ) write(1000,*) 'Se', (pos_eq(i,j),j=1,3)
!end do

  
!do l = 1, nlayer
!   write(1000+l,'(a)') 'ITEM: TIMESTEP'
!   write(1000+l,*) 0
!   write(1000+l,'(a)') 'ITEM: NUMBER OF ATOMS'
!   write(1000+l,*) atoms_layer(l)
!   write(1000+l,'(a)') 'ITEM: BOX BOUNDS pp pp pp'
!   write(1000+l,*) 0.d0, side(1)
!   write(1000+l,*) 0.d0, side(2)
!   write(1000+l,*) 0.d0, side(3)
!   write(1000+l,'(a)') 'ITEM: ATOMS element xu yu zu'
!   do i = 1, atoms
!      if ( layer(i) == l ) write(1000+l,*) 'X', (pos_eq(i,j),j=1,3)
!   end do
!end do
!stop

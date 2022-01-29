module vars
  integer, parameter :: nlayer = 200
  real(8), parameter :: thick = 0.5d0 ! Ang, thickness of one atomic layer
  integer, save :: atoms, nbulk_Mo, nedge_Mo, nbulk_S, nedge_S, fixed(3)
  integer, save :: atoms_layer(nlayer)
  real(8), save :: side(3)
  integer, allocatable, save :: tipo(:), layer(:)
  real(8), allocatable, save :: mass(:), pos_eq(:,:), poteng_eq(:), distance(:,:)
  real(8), allocatable, save :: pos(:,:), vel(:,:), force(:,:), poteng(:), kineng(:,:)
  character(6), allocatable, save :: tag(:)
end module vars

program maps
  use vars, only: atoms, mass, tag, pos_eq, fixed, pos, vel, force, poteng, poteng_eq, kineng, tipo
  implicit none
  real(8), parameter :: metal_to_eV = 1.d-3 / 6.02214076d0 / 1.602176634d0
  real(8), parameter :: tol = 1.d-3
  integer :: i, j, id, iconf, timestep, timestep0, ny, iperiod, min
  real(8) :: dt, pos_tmp(3), poteng_tmp, vel_tmp(3), force_tmp(3)
  real(8) :: err, displ(3), time, ly, velocity, displacement, offset, period, position
  character(20) :: exclude
  real(8), allocatable :: displ_ave(:,:), vel_ave(:,:), force_ave(:,:), poteng_ave(:), kineng_ave(:,:), rmsd_ave(:,:), rmsv_ave(:,:)

  read(*,*) velocity, dt, min, displacement, ly, ny, exclude ! velocity in Ang/ps, dt in ps

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

  allocate ( displ_ave(atoms,4), stat = i )
  if ( i /= 0 ) stop 'Allocation DISPL_AVE'
  allocate ( vel_ave(atoms,4), stat = i )
  if ( i /= 0 ) stop 'Allocation VEL_AVE'
  allocate ( force_ave(atoms,4), stat = i )
  if ( i /= 0 ) stop 'Allocation FORCE_AVE'
  allocate ( poteng_ave(atoms), stat = i )
  if ( i /= 0 ) stop 'Allocation PE_AVE'
  allocate ( kineng_ave(atoms,4), stat = i )
  if ( i /= 0 ) stop 'Allocation KE_AVE'
  allocate ( rmsd_ave(atoms,4), stat = i )
  if ( i /= 0 ) stop 'Allocation RMSD_AVE'
  allocate ( rmsv_ave(atoms,4), stat = i )
  if ( i /= 0 ) stop 'Allocation RMSV_AVE'
  
  open ( unit=100, file='system.lammpstrj' )
  iconf = 0
  iperiod = 0
  displ_ave(:,:) = 0.d0
  vel_ave(:,:) = 0.d0
  force_ave(:,:) = 0.d0
  poteng_ave(:) = 0.d0
  kineng_ave(:,:) = 0.d0
  rmsd_ave(:,:) = 0.d0
  rmsv_ave(:,:) = 0.d0
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
        if ( err > tol ) write(0,*) 'CAREFUL NOW! ', i, err!, (j, pos(fixed(i),j), pos_eq(fixed(i),j),j=1,3)
     end do

     ! averages
     do i = 1, atoms
        ! displacements
        displ(:) = pos(i,:) - pos_eq(i,:)
        displ_ave(i,1:3) = displ_ave(i,1:3) + abs(displ(1:3))
        displ_ave(i,4) = displ_ave(i,4) + sqrt(dot_product(displ,displ))
        ! velocities
        vel_ave(i,1:3) = vel_ave(i,1:3) + abs(vel(i,1:3))
        vel_ave(i,4) = vel_ave(i,4) + sqrt(dot_product(vel(i,:),vel(i,:)))
        ! forces
        force_ave(i,1:3) = force_ave(i,1:3) + abs(force(i,1:3))
        force_ave(i,4) = force_ave(i,4) + sqrt(dot_product(force(i,:),force(i,:)))
        ! potential energies
        poteng_ave(i) = poteng_ave(i) + poteng(i)-poteng_eq(i)
        ! kinetic energies
        kineng_ave(i,1:3) = kineng_ave(i,1:3) + kineng(i,1:3)
        kineng_ave(i,4) = kineng_ave(i,4) + sum(kineng(i,:))
        ! rmsd
        rmsd_ave(i,1:3) = rmsd_ave(i,1:3) + displ(1:3)*displ(1:3)
        rmsd_ave(i,4) = rmsd_ave(i,4) + dot_product(displ,displ)
        ! rmsv
        rmsv_ave(i,1:3) = rmsv_ave(i,1:3) + vel(i,1:3)*vel(i,1:3)
        rmsv_ave(i,4) = rmsv_ave(i,4) + dot_product(vel(i,:),vel(i,:))
     end do
     
  end do
  close(100)

  ! printout
  open ( unit=200, file='5-displ_Mo.dat' )
  open ( unit=201, file='5-displ_Sbottom.dat' )
  open ( unit=202, file='5-displ_Stop.dat' )
  open ( unit=300, file='5-vel_Mo.dat' )
  open ( unit=301, file='5-vel_Sbottom.dat' )
  open ( unit=302, file='5-vel_Stop.dat' )
  open ( unit=400, file='5-force_Mo.dat' )
  open ( unit=401, file='5-force_Sbottom.dat' )
  open ( unit=402, file='5-force_Stop.dat' )
  open ( unit=500, file='5-poteng_Mo.dat' )
  open ( unit=501, file='5-poteng_Sbottom.dat' )
  open ( unit=502, file='5-poteng_Stop.dat' )
  open ( unit=600, file='5-kineng_Mo.dat' )
  open ( unit=601, file='5-kineng_Sbottom.dat' )
  open ( unit=602, file='5-kineng_Stop.dat' )
  open ( unit=700, file='5-rmsd_Mo.dat' )
  open ( unit=701, file='5-rmsd_Sbottom.dat' )
  open ( unit=702, file='5-rmsd_Stop.dat' )
  open ( unit=800, file='5-rmsv_Mo.dat' )
  open ( unit=801, file='5-rmsv_Sbottom.dat' )
  open ( unit=802, file='5-rmsv_Stop.dat' )
  do i = 1, atoms
     if ( tag(i) /= 'none' ) then
        if ( tipo(i) == 2 ) then
           write(200,*) (pos_eq(i,j),j=1,3), (displ_ave(i,j)/dble(iconf),j=1,4)
           write(300,*) (pos_eq(i,j),j=1,3), (vel_ave(i,j)/dble(iconf),j=1,4)
           write(400,*) (pos_eq(i,j),j=1,3), (force_ave(i,j)/dble(iconf),j=1,4)
           write(500,*) (pos_eq(i,j),j=1,3), poteng_ave(i)/dble(iconf)
           write(600,*) (pos_eq(i,j),j=1,3), (kineng_ave(i,j)/dble(iconf),j=1,4)
           write(700,*) (pos_eq(i,j),j=1,3), (sqrt(rmsd_ave(i,j)/dble(iconf)),j=1,4)
           write(800,*) (pos_eq(i,j),j=1,3), (sqrt(rmsv_ave(i,j)/dble(iconf)),j=1,4)
        else if ( tipo(i) == 1 ) then
           write(201,*) (pos_eq(i,j),j=1,3), (displ_ave(i,j)/dble(iconf),j=1,4)
           write(301,*) (pos_eq(i,j),j=1,3), (vel_ave(i,j)/dble(iconf),j=1,4)
           write(401,*) (pos_eq(i,j),j=1,3), (force_ave(i,j)/dble(iconf),j=1,4)
           write(501,*) (pos_eq(i,j),j=1,3), poteng_ave(i)/dble(iconf)
           write(601,*) (pos_eq(i,j),j=1,3), (kineng_ave(i,j)/dble(iconf),j=1,4)
           write(701,*) (pos_eq(i,j),j=1,3), (sqrt(rmsd_ave(i,j)/dble(iconf)),j=1,4)
           write(801,*) (pos_eq(i,j),j=1,3), (sqrt(rmsv_ave(i,j)/dble(iconf)),j=1,4)
        else if ( tipo(i) == 3 ) then
           write(202,*) (pos_eq(i,j),j=1,3), (displ_ave(i,j)/dble(iconf),j=1,4)
           write(302,*) (pos_eq(i,j),j=1,3), (vel_ave(i,j)/dble(iconf),j=1,4)
           write(402,*) (pos_eq(i,j),j=1,3), (force_ave(i,j)/dble(iconf),j=1,4)
           write(502,*) (pos_eq(i,j),j=1,3), poteng_ave(i)/dble(iconf)
           write(602,*) (pos_eq(i,j),j=1,3), (kineng_ave(i,j)/dble(iconf),j=1,4)
           write(702,*) (pos_eq(i,j),j=1,3), (sqrt(rmsd_ave(i,j)/dble(iconf)),j=1,4)
           write(802,*) (pos_eq(i,j),j=1,3), (sqrt(rmsv_ave(i,j)/dble(iconf)),j=1,4)
        end if
     end if
  end do
  do i = 2, 8
     close(100*i)
     close(100*i+1)
     close(100*i+2)
  end do

  stop
end program maps

subroutine read_eq
  use vars, only: atoms, side, tipo, mass, pos_eq, distance, nlayer, &
       atoms_layer, layer, atoms, poteng_eq, thick
  implicit none
  integer :: i, j, id, tipo_tmp, atom_types, ix(3)
  integer :: iA, iB, iC, l
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

  return
end subroutine read_eq

subroutine read_lists
  use vars, only: atoms, tipo, tag, nbulk_Mo, nedge_Mo, nbulk_S, nedge_S, fixed, pos_eq
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

  ! fake definition of tags
  do i = 1, atoms
     if ( tipo(i) <= 3 ) tag(i) = 'scio'
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

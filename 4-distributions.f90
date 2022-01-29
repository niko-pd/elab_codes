module vars
  integer, parameter :: nlayer = 200
  real(8), parameter :: thick = 0.5d0 ! Ang, thickness of one atomic layer
  integer, save :: atoms, nedge_Mo, nedge_S, n2nd_Mo, n2nd_S, n3rd_Mo, n3rd_S, nbulk_Mo, nbulk_S, fixed(3)
  integer, save :: atoms_layer(nlayer)
  real(8), save :: side(3), thickness_edge
  integer, allocatable, save :: tipo(:), layer(:)
  real(8), allocatable, save :: mass(:), pos_eq(:,:), poteng_eq(:), distance(:,:)
  real(8), allocatable, save :: pos(:,:), vel(:,:), force(:,:), poteng(:), kineng(:,:)
  character(6), allocatable, save :: tag(:)
end module vars

program distributions
  use vars, only: atoms, mass, tag, pos_eq, fixed, pos, vel, force, poteng, poteng_eq, kineng, tipo, nbulk_Mo, nbulk_S, nedge_Mo, nedge_S, thickness_edge, n2nd_Mo, n2nd_S, n3rd_Mo, n3rd_S
  implicit none
  real(8), parameter :: metal_to_eV = 1.d-3 / 6.02214076d0 / 1.602176634d0
  real(8), parameter :: tol = 1.d-3
  ! displacements
  real(8), parameter :: displ_dbin=0.005d0, displ_max=5.d0-displ_dbin/2.d0, displ_min=-displ_max ! Ang
  integer, parameter :: displ_nbin=int((displ_max-displ_min)/displ_dbin)
  ! velocities
  real(8), parameter :: vel_dbin=0.1d0, vel_max=20.d0-vel_dbin/2.d0, vel_min=-vel_max ! Ang/ps
  integer, parameter :: vel_nbin=int((vel_max-vel_min)/vel_dbin)
  ! forces
  real(8), parameter :: force_dbin=0.05d0, force_max=10.d0-force_dbin/2.d0, force_min=-force_max ! eV/Ang
  integer, parameter :: force_nbin=int((force_max-force_min)/force_dbin)
  ! potential energies
  real(8), parameter :: poteng_dbin=0.0002d0, poteng_max=0.5d0-poteng_dbin/2.d0, poteng_min=-poteng_max ! eV
  integer, parameter :: poteng_nbin=int((poteng_max-poteng_min)/poteng_dbin)
  ! kinetic energies
  real(8), parameter :: kineng_dbin=0.001d0, kineng_max=0.6d0-kineng_dbin/2.d0, kineng_min=-kineng_max ! eV
  integer, parameter :: kineng_nbin=int((kineng_max-kineng_min)/kineng_dbin) 
  integer :: i, j, k, id, iconf, timestep, timestep0, ny, iperiod, min
  real(8) :: dt, pos_tmp(3), poteng_tmp, vel_tmp(3), force_tmp(3)
  real(8) :: err, displ(3), time, ly, velocity, displacement, offset, period, position
  character(20) :: exclude
  real(8) :: displ_distr_edge_Mo(0:displ_nbin,4), vel_distr_edge_Mo(0:vel_nbin,4), force_distr_edge_Mo(0:force_nbin,4)
  real(8) :: poteng_distr_edge_Mo(0:poteng_nbin), kineng_distr_edge_Mo(0:kineng_nbin,4)
  real(8) :: displ_distr_edge_S(0:displ_nbin,4), vel_distr_edge_S(0:vel_nbin,4), force_distr_edge_S(0:force_nbin,4)
  real(8) :: poteng_distr_edge_S(0:poteng_nbin), kineng_distr_edge_S(0:kineng_nbin,4)
  real(8) :: displ_distr_2nd_Mo(0:displ_nbin,4), vel_distr_2nd_Mo(0:vel_nbin,4), force_distr_2nd_Mo(0:force_nbin,4)
  real(8) :: poteng_distr_2nd_Mo(0:poteng_nbin), kineng_distr_2nd_Mo(0:kineng_nbin,4)
  real(8) :: displ_distr_2nd_S(0:displ_nbin,4), vel_distr_2nd_S(0:vel_nbin,4), force_distr_2nd_S(0:force_nbin,4)
  real(8) :: poteng_distr_2nd_S(0:poteng_nbin), kineng_distr_2nd_S(0:kineng_nbin,4)
  real(8) :: displ_distr_3rd_Mo(0:displ_nbin,4), vel_distr_3rd_Mo(0:vel_nbin,4), force_distr_3rd_Mo(0:force_nbin,4)
  real(8) :: poteng_distr_3rd_Mo(0:poteng_nbin), kineng_distr_3rd_Mo(0:kineng_nbin,4)
  real(8) :: displ_distr_3rd_S(0:displ_nbin,4), vel_distr_3rd_S(0:vel_nbin,4), force_distr_3rd_S(0:force_nbin,4)
  real(8) :: poteng_distr_3rd_S(0:poteng_nbin), kineng_distr_3rd_S(0:kineng_nbin,4)
  real(8) :: displ_distr_bulk_Mo(0:displ_nbin,4), vel_distr_bulk_Mo(0:vel_nbin,4), force_distr_bulk_Mo(0:force_nbin,4)
  real(8) :: poteng_distr_bulk_Mo(0:poteng_nbin), kineng_distr_bulk_Mo(0:kineng_nbin,4)
  real(8) :: displ_distr_bulk_S(0:displ_nbin,4), vel_distr_bulk_S(0:vel_nbin,4), force_distr_bulk_S(0:force_nbin,4)
  real(8) :: poteng_distr_bulk_S(0:poteng_nbin), kineng_distr_bulk_S(0:kineng_nbin,4)
  real(8) :: rmsd_distr_edge(0:displ_nbin), rmsd_distr_2nd(0:displ_nbin), rmsd_distr_3rd(0:displ_nbin), rmsd_distr_bulk(0:displ_nbin)
  real(8) :: rmsd_edge, rmsd_2nd, rmsd_3rd, rmsd_bulk
  real(8) :: eq_pe_edge_Mo, eq_pe_edge_S, eq_pe_2nd_Mo, eq_pe_2nd_S, eq_pe_3rd_Mo, eq_pe_3rd_S, eq_pe_bulk_Mo,  eq_pe_bulk_S
  real(8) :: pe_edge_Mo, pe_edge_S, pe_2nd_Mo, pe_2nd_S, pe_3rd_Mo, pe_3rd_S, pe_bulk_Mo, pe_bulk_S
  real(8) :: pe_distr_edge_Mo(0:poteng_nbin), pe_distr_2nd_Mo(0:poteng_nbin), pe_distr_3rd_Mo(0:poteng_nbin), pe_distr_bulk_Mo(0:poteng_nbin) 
  real(8) :: pe_distr_edge_S(0:poteng_nbin), pe_distr_2nd_S(0:poteng_nbin), pe_distr_3rd_S(0:poteng_nbin), pe_distr_bulk_S(0:poteng_nbin) 

  read(*,*) velocity, dt, min, displacement, ly, ny, exclude, thickness_edge ! velocity in Ang/ps, dt in ps

  ! set parameters
  period = ly / dble(ny) ! Ang
  offset = dble(min) * displacement ! Ang
  if ( exclude == 'first' ) offset = offset + period
  
  call read_eq

  call read_lists

  eq_pe_edge_Mo = 0.d0
  eq_pe_edge_S = 0.d0
  eq_pe_2nd_Mo = 0.d0
  eq_pe_2nd_S = 0.d0
  eq_pe_3rd_Mo = 0.d0
  eq_pe_3rd_S = 0.d0
  eq_pe_bulk_Mo = 0.d0
  eq_pe_bulk_S = 0.d0
  do i = 1, atoms
     if ( tipo(i) == 2 ) then
        if ( tag(i) == 'edge' ) then
           eq_pe_edge_Mo = eq_pe_edge_Mo + poteng_eq(i)
        else if ( tag(i) == '2nd' ) then
           eq_pe_2nd_Mo = eq_pe_2nd_Mo + poteng_eq(i)
        else if ( tag(i) == '3rd' ) then
           eq_pe_3rd_Mo = eq_pe_3rd_Mo + poteng_eq(i)
        else if ( tag(i) == 'bulk' ) then
           eq_pe_bulk_Mo = eq_pe_bulk_Mo + poteng_eq(i)
        end if
     elseif ( tipo(i) == 1 .or. tipo(i) == 3 ) then
        if ( tag(i) == 'edge' ) then
           eq_pe_edge_S = eq_pe_edge_S + poteng_eq(i)
        else if ( tag(i) == '2nd' ) then
           eq_pe_2nd_S = eq_pe_2nd_S + poteng_eq(i)
        else if ( tag(i) == '3rd' ) then
           eq_pe_3rd_S = eq_pe_3rd_S + poteng_eq(i)
        else if ( tag(i) == 'bulk' ) then
           eq_pe_bulk_S = eq_pe_bulk_S + poteng_eq(i)
        end if
     end if
  end do
  eq_pe_edge_Mo = eq_pe_edge_Mo / dble(nedge_Mo) 
  eq_pe_edge_S = eq_pe_edge_S / dble(nedge_S) 
  eq_pe_2nd_Mo = eq_pe_2nd_Mo / dble(n2nd_Mo) 
  eq_pe_2nd_S = eq_pe_2nd_S / dble(n2nd_S) 
  eq_pe_3rd_Mo = eq_pe_3rd_Mo / dble(n3rd_Mo) 
  eq_pe_3rd_S = eq_pe_3rd_S / dble(n3rd_S) 
  eq_pe_bulk_Mo = eq_pe_bulk_Mo / dble(nbulk_Mo) 
  eq_pe_bulk_S = eq_pe_bulk_S / dble(nbulk_S) 
  
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
  iconf = 0
  iperiod = 0
  displ_distr_edge_Mo(:,:) = 0.d0
  vel_distr_edge_Mo(:,:) = 0.d0
  force_distr_edge_Mo(:,:) = 0.d0
  poteng_distr_edge_Mo(:) = 0.d0
  kineng_distr_edge_Mo(:,:) = 0.d0
  displ_distr_edge_S(:,:) = 0.d0
  vel_distr_edge_S(:,:) = 0.d0
  force_distr_edge_S(:,:) = 0.d0
  poteng_distr_edge_S(:) = 0.d0
  kineng_distr_edge_S(:,:) = 0.d0
  rmsd_distr_edge(:) = 0.d0
  pe_distr_edge_Mo(:) = 0.d0
  pe_distr_edge_S(:) = 0.d0
  displ_distr_2nd_Mo(:,:) = 0.d0
  vel_distr_2nd_Mo(:,:) = 0.d0
  force_distr_2nd_Mo(:,:) = 0.d0
  poteng_distr_2nd_Mo(:) = 0.d0
  kineng_distr_2nd_Mo(:,:) = 0.d0
  displ_distr_2nd_S(:,:) = 0.d0
  vel_distr_2nd_S(:,:) = 0.d0
  force_distr_2nd_S(:,:) = 0.d0
  poteng_distr_2nd_S(:) = 0.d0
  kineng_distr_2nd_S(:,:) = 0.d0
  rmsd_distr_2nd(:) = 0.d0
  pe_distr_2nd_Mo(:) = 0.d0
  pe_distr_2nd_S(:) = 0.d0
  displ_distr_3rd_Mo(:,:) = 0.d0
  vel_distr_3rd_Mo(:,:) = 0.d0
  force_distr_3rd_Mo(:,:) = 0.d0
  poteng_distr_3rd_Mo(:) = 0.d0
  kineng_distr_3rd_Mo(:,:) = 0.d0
  displ_distr_3rd_S(:,:) = 0.d0
  vel_distr_3rd_S(:,:) = 0.d0
  force_distr_3rd_S(:,:) = 0.d0
  poteng_distr_3rd_S(:) = 0.d0
  kineng_distr_3rd_S(:,:) = 0.d0
  rmsd_distr_3rd(:) = 0.d0
  pe_distr_3rd_Mo(:) = 0.d0
  pe_distr_3rd_S(:) = 0.d0
  displ_distr_bulk_Mo(:,:) = 0.d0
  vel_distr_bulk_Mo(:,:) = 0.d0
  force_distr_bulk_Mo(:,:) = 0.d0
  poteng_distr_bulk_Mo(:) = 0.d0
  kineng_distr_bulk_Mo(:,:) = 0.d0
  displ_distr_bulk_S(:,:) = 0.d0
  vel_distr_bulk_S(:,:) = 0.d0
  force_distr_bulk_S(:,:) = 0.d0
  poteng_distr_bulk_S(:) = 0.d0
  kineng_distr_bulk_S(:,:) = 0.d0
  rmsd_distr_bulk(:) = 0.d0
  pe_distr_bulk_Mo(:) = 0.d0
  pe_distr_bulk_S(:) = 0.d0
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

     ! distributions
     rmsd_edge = 0.d0
     pe_edge_Mo = 0.d0
     pe_edge_S = 0.d0
     rmsd_2nd = 0.d0
     pe_2nd_Mo = 0.d0
     pe_2nd_S = 0.d0
     rmsd_3rd = 0.d0
     pe_3rd_Mo = 0.d0
     pe_3rd_S = 0.d0
     rmsd_bulk = 0.d0
     pe_bulk_Mo = 0.d0
     pe_bulk_S = 0.d0
     do i = 1, atoms
        if ( tag(i) == 'edge' ) then
           if ( tipo(i) == 2 ) then
              ! displacements
              displ(:) = pos(i,:) - pos_eq(i,:)
              do j = 1, 3
                 k = int( ( displ(j) - displ_min ) / displ_dbin )
                 displ_distr_edge_Mo(k,j) = displ_distr_edge_Mo(k,j) + 1.d0
              end do
              k = int( ( sqrt(dot_product(displ,displ)) - displ_min ) / displ_dbin )
              displ_distr_edge_Mo(k,4) = displ_distr_edge_Mo(k,4) + 1.d0
              ! velocities
              do j = 1, 3
                 k = int( ( vel(i,j) - vel_min ) / vel_dbin )
                 vel_distr_edge_Mo(k,j) = vel_distr_edge_Mo(k,j) + 1.d0
              end do
              k = int( ( sqrt(dot_product(vel(i,:),vel(i,:))) - vel_min ) / vel_dbin )
              vel_distr_edge_Mo(k,4) = vel_distr_edge_Mo(k,4) + 1.d0
              ! forces
              do j = 1, 3
                 k = int( ( force(i,j) - force_min ) / force_dbin )
                 force_distr_edge_Mo(k,j) = force_distr_edge_Mo(k,j) + 1.d0
              end do
              k = int( ( sqrt(dot_product(force(i,:),force(i,:))) - force_min ) / force_dbin )
              force_distr_edge_Mo(k,4) = force_distr_edge_Mo(k,4) + 1.d0
              ! potential energies
              k = int( ( poteng(i)-poteng_eq(i) - poteng_min ) / poteng_dbin )
              poteng_distr_edge_Mo(k) = poteng_distr_edge_Mo(k) + 1.d0
              pe_edge_Mo = pe_edge_Mo + poteng(i)
              ! kinetic energies
              do j = 1, 3
                 k = int( ( kineng(i,j) - kineng_min ) / kineng_dbin )
                 kineng_distr_edge_Mo(k,j) = kineng_distr_edge_Mo(k,j) + 1.d0
              end do
              k = int( ( sum(kineng(i,:)) - kineng_min ) / kineng_dbin )
              kineng_distr_edge_Mo(k,4) = kineng_distr_edge_Mo(k,4) + 1.d0
           else if ( tipo(i) == 1 .or. tipo(i) == 3 ) then
              ! displacements
              displ(:) = pos(i,:) - pos_eq(i,:)
              do j = 1, 3
                 k = int( ( displ(j) - displ_min ) / displ_dbin )
                 displ_distr_edge_S(k,j) = displ_distr_edge_S(k,j) + 1.d0
              end do
              k = int( ( sqrt(dot_product(displ,displ)) - displ_min ) / displ_dbin )
              displ_distr_edge_S(k,4) = displ_distr_edge_S(k,4) + 1.d0
              ! velocities
              do j = 1, 3
                 k = int( ( vel(i,j) - vel_min ) / vel_dbin )
                 vel_distr_edge_S(k,j) = vel_distr_edge_S(k,j) + 1.d0
              end do
              k = int( ( sqrt(dot_product(vel(i,:),vel(i,:))) - vel_min ) / vel_dbin )
              vel_distr_edge_S(k,4) = vel_distr_edge_S(k,4) + 1.d0
              ! forces
              do j = 1, 3
                 k = int( ( force(i,j) - force_min ) / force_dbin )
                 force_distr_edge_S(k,j) = force_distr_edge_S(k,j) + 1.d0
              end do
              k = int( ( sqrt(dot_product(force(i,:),force(i,:))) - force_min ) / force_dbin )
              force_distr_edge_S(k,4) = force_distr_edge_S(k,4) + 1.d0
              ! potential energies
              k = int( ( poteng(i)-poteng_eq(i) - poteng_min ) / poteng_dbin )
              poteng_distr_edge_S(k) = poteng_distr_edge_S(k) + 1.d0
              pe_edge_S = pe_edge_S + poteng(i)
              ! kinetic energies
              do j = 1, 3
                 k = int( ( kineng(i,j) - kineng_min ) / kineng_dbin )
                 kineng_distr_edge_S(k,j) = kineng_distr_edge_S(k,j) + 1.d0
              end do
              k = int( ( sum(kineng(i,:)) - kineng_min ) / kineng_dbin )
              kineng_distr_edge_S(k,4) = kineng_distr_edge_S(k,4) + 1.d0
           end if
           rmsd_edge = rmsd_edge + dot_product(displ,displ)
        else if ( tag(i) == '2nd' ) then
           if ( tipo(i) == 2 ) then
              ! displacements
              displ(:) = pos(i,:) - pos_eq(i,:)
              do j = 1, 3
                 k = int( ( displ(j) - displ_min ) / displ_dbin )
                 displ_distr_2nd_Mo(k,j) = displ_distr_2nd_Mo(k,j) + 1.d0
              end do
              k = int( ( sqrt(dot_product(displ,displ)) - displ_min ) / displ_dbin )
              displ_distr_2nd_Mo(k,4) = displ_distr_2nd_Mo(k,4) + 1.d0
              ! velocities
              do j = 1, 3
                 k = int( ( vel(i,j) - vel_min ) / vel_dbin )
                 vel_distr_2nd_Mo(k,j) = vel_distr_2nd_Mo(k,j) + 1.d0
              end do
              k = int( ( sqrt(dot_product(vel(i,:),vel(i,:))) - vel_min ) / vel_dbin )
              vel_distr_2nd_Mo(k,4) = vel_distr_2nd_Mo(k,4) + 1.d0
              ! forces
              do j = 1, 3
                 k = int( ( force(i,j) - force_min ) / force_dbin )
                 force_distr_2nd_Mo(k,j) = force_distr_2nd_Mo(k,j) + 1.d0
              end do
              k = int( ( sqrt(dot_product(force(i,:),force(i,:))) - force_min ) / force_dbin )
              force_distr_2nd_Mo(k,4) = force_distr_2nd_Mo(k,4) + 1.d0
              ! potential energies
              k = int( ( poteng(i)-poteng_eq(i) - poteng_min ) / poteng_dbin )
              poteng_distr_2nd_Mo(k) = poteng_distr_2nd_Mo(k) + 1.d0
              pe_2nd_Mo = pe_2nd_Mo + poteng(i)
              ! kinetic energies
              do j = 1, 3
                 k = int( ( kineng(i,j) - kineng_min ) / kineng_dbin )
                 kineng_distr_2nd_Mo(k,j) = kineng_distr_2nd_Mo(k,j) + 1.d0
              end do
              k = int( ( sum(kineng(i,:)) - kineng_min ) / kineng_dbin )
              kineng_distr_2nd_Mo(k,4) = kineng_distr_2nd_Mo(k,4) + 1.d0
           else if ( tipo(i) == 1 .or. tipo(i) == 3 ) then
              ! displacements
              displ(:) = pos(i,:) - pos_eq(i,:)
              do j = 1, 3
                 k = int( ( displ(j) - displ_min ) / displ_dbin )
                 displ_distr_2nd_S(k,j) = displ_distr_2nd_S(k,j) + 1.d0
              end do
              k = int( ( sqrt(dot_product(displ,displ)) - displ_min ) / displ_dbin )
              displ_distr_2nd_S(k,4) = displ_distr_2nd_S(k,4) + 1.d0
              ! velocities
              do j = 1, 3
                 k = int( ( vel(i,j) - vel_min ) / vel_dbin )
                 vel_distr_2nd_S(k,j) = vel_distr_2nd_S(k,j) + 1.d0
              end do
              k = int( ( sqrt(dot_product(vel(i,:),vel(i,:))) - vel_min ) / vel_dbin )
              vel_distr_2nd_S(k,4) = vel_distr_2nd_S(k,4) + 1.d0
              ! forces
              do j = 1, 3
                 k = int( ( force(i,j) - force_min ) / force_dbin )
                 force_distr_2nd_S(k,j) = force_distr_2nd_S(k,j) + 1.d0
              end do
              k = int( ( sqrt(dot_product(force(i,:),force(i,:))) - force_min ) / force_dbin )
              force_distr_2nd_S(k,4) = force_distr_2nd_S(k,4) + 1.d0
              ! potential energies
              k = int( ( poteng(i)-poteng_eq(i) - poteng_min ) / poteng_dbin )
              poteng_distr_2nd_S(k) = poteng_distr_2nd_S(k) + 1.d0
              pe_2nd_S = pe_2nd_S + poteng(i)
              ! kinetic energies
              do j = 1, 3
                 k = int( ( kineng(i,j) - kineng_min ) / kineng_dbin )
                 kineng_distr_2nd_S(k,j) = kineng_distr_2nd_S(k,j) + 1.d0
              end do
              k = int( ( sum(kineng(i,:)) - kineng_min ) / kineng_dbin )
              kineng_distr_2nd_S(k,4) = kineng_distr_2nd_S(k,4) + 1.d0
           end if
           rmsd_2nd = rmsd_2nd + dot_product(displ,displ)
        else if ( tag(i) == '3rd' ) then
           if ( tipo(i) == 2 ) then
              ! displacements
              displ(:) = pos(i,:) - pos_eq(i,:)
              do j = 1, 3
                 k = int( ( displ(j) - displ_min ) / displ_dbin )
                 displ_distr_3rd_Mo(k,j) = displ_distr_3rd_Mo(k,j) + 1.d0
              end do
              k = int( ( sqrt(dot_product(displ,displ)) - displ_min ) / displ_dbin )
              displ_distr_3rd_Mo(k,4) = displ_distr_3rd_Mo(k,4) + 1.d0
              ! velocities
              do j = 1, 3
                 k = int( ( vel(i,j) - vel_min ) / vel_dbin )
                 vel_distr_3rd_Mo(k,j) = vel_distr_3rd_Mo(k,j) + 1.d0
              end do
              k = int( ( sqrt(dot_product(vel(i,:),vel(i,:))) - vel_min ) / vel_dbin )
              vel_distr_3rd_Mo(k,4) = vel_distr_3rd_Mo(k,4) + 1.d0
              ! forces
              do j = 1, 3
                 k = int( ( force(i,j) - force_min ) / force_dbin )
                 force_distr_3rd_Mo(k,j) = force_distr_3rd_Mo(k,j) + 1.d0
              end do
              k = int( ( sqrt(dot_product(force(i,:),force(i,:))) - force_min ) / force_dbin )
              force_distr_3rd_Mo(k,4) = force_distr_3rd_Mo(k,4) + 1.d0
              ! potential energies
              k = int( ( poteng(i)-poteng_eq(i) - poteng_min ) / poteng_dbin )
              poteng_distr_3rd_Mo(k) = poteng_distr_3rd_Mo(k) + 1.d0
              pe_3rd_Mo = pe_3rd_Mo + poteng(i)
              ! kinetic energies
              do j = 1, 3
                 k = int( ( kineng(i,j) - kineng_min ) / kineng_dbin )
                 kineng_distr_3rd_Mo(k,j) = kineng_distr_3rd_Mo(k,j) + 1.d0
              end do
              k = int( ( sum(kineng(i,:)) - kineng_min ) / kineng_dbin )
              kineng_distr_3rd_Mo(k,4) = kineng_distr_3rd_Mo(k,4) + 1.d0
           else if ( tipo(i) == 1 .or. tipo(i) == 3 ) then
              ! displacements
              displ(:) = pos(i,:) - pos_eq(i,:)
              do j = 1, 3
                 k = int( ( displ(j) - displ_min ) / displ_dbin )
                 displ_distr_3rd_S(k,j) = displ_distr_3rd_S(k,j) + 1.d0
              end do
              k = int( ( sqrt(dot_product(displ,displ)) - displ_min ) / displ_dbin )
              displ_distr_3rd_S(k,4) = displ_distr_3rd_S(k,4) + 1.d0
              ! velocities
              do j = 1, 3
                 k = int( ( vel(i,j) - vel_min ) / vel_dbin )
                 vel_distr_3rd_S(k,j) = vel_distr_3rd_S(k,j) + 1.d0
              end do
              k = int( ( sqrt(dot_product(vel(i,:),vel(i,:))) - vel_min ) / vel_dbin )
              vel_distr_3rd_S(k,4) = vel_distr_3rd_S(k,4) + 1.d0
              ! forces
              do j = 1, 3
                 k = int( ( force(i,j) - force_min ) / force_dbin )
                 force_distr_3rd_S(k,j) = force_distr_3rd_S(k,j) + 1.d0
              end do
              k = int( ( sqrt(dot_product(force(i,:),force(i,:))) - force_min ) / force_dbin )
              force_distr_3rd_S(k,4) = force_distr_3rd_S(k,4) + 1.d0
              ! potential energies
              k = int( ( poteng(i)-poteng_eq(i) - poteng_min ) / poteng_dbin )
              poteng_distr_3rd_S(k) = poteng_distr_3rd_S(k) + 1.d0
              pe_3rd_S = pe_3rd_S + poteng(i)
              ! kinetic energies
              do j = 1, 3
                 k = int( ( kineng(i,j) - kineng_min ) / kineng_dbin )
                 kineng_distr_3rd_S(k,j) = kineng_distr_3rd_S(k,j) + 1.d0
              end do
              k = int( ( sum(kineng(i,:)) - kineng_min ) / kineng_dbin )
              kineng_distr_3rd_S(k,4) = kineng_distr_3rd_S(k,4) + 1.d0
           end if
           rmsd_3rd = rmsd_3rd + dot_product(displ,displ)
        else if ( tag(i) == 'bulk' ) then
           if ( tipo(i) == 2 ) then
              ! displacements
              displ(:) = pos(i,:) - pos_eq(i,:)
              do j = 1, 3
                 k = int( ( displ(j) - displ_min ) / displ_dbin )
                 displ_distr_bulk_Mo(k,j) = displ_distr_bulk_Mo(k,j) + 1.d0
              end do
              k = int( ( sqrt(dot_product(displ,displ)) - displ_min ) / displ_dbin )
              displ_distr_bulk_Mo(k,4) = displ_distr_bulk_Mo(k,4) + 1.d0
              ! velocities
              do j = 1, 3
                 k = int( ( vel(i,j) - vel_min ) / vel_dbin )
                 vel_distr_bulk_Mo(k,j) = vel_distr_bulk_Mo(k,j) + 1.d0
              end do
              k = int( ( sqrt(dot_product(vel(i,:),vel(i,:))) - vel_min ) / vel_dbin )
              vel_distr_bulk_Mo(k,4) = vel_distr_bulk_Mo(k,4) + 1.d0
              ! forces
              do j = 1, 3
                 k = int( ( force(i,j) - force_min ) / force_dbin )
                 force_distr_bulk_Mo(k,j) = force_distr_bulk_Mo(k,j) + 1.d0
              end do
              k = int( ( sqrt(dot_product(force(i,:),force(i,:))) - force_min ) / force_dbin )
              force_distr_bulk_Mo(k,4) = force_distr_bulk_Mo(k,4) + 1.d0
              ! potential energies
              k = int( ( poteng(i)-poteng_eq(i) - poteng_min ) / poteng_dbin )
              poteng_distr_bulk_Mo(k) = poteng_distr_bulk_Mo(k) + 1.d0
              pe_bulk_Mo = pe_bulk_Mo + poteng(i)
              ! kinetic energies
              do j = 1, 3
                 k = int( ( kineng(i,j) - kineng_min ) / kineng_dbin )
                 kineng_distr_bulk_Mo(k,j) = kineng_distr_bulk_Mo(k,j) + 1.d0
              end do
              k = int( ( sum(kineng(i,:)) - kineng_min ) / kineng_dbin )
              kineng_distr_bulk_Mo(k,4) = kineng_distr_bulk_Mo(k,4) + 1.d0
           else if ( tipo(i) == 1 .or. tipo(i) == 3 ) then
              ! displacements
              displ(:) = pos(i,:) - pos_eq(i,:)
              do j = 1, 3
                 k = int( ( displ(j) - displ_min ) / displ_dbin )
                 displ_distr_bulk_S(k,j) = displ_distr_bulk_S(k,j) + 1.d0
              end do
              k = int( ( sqrt(dot_product(displ,displ)) - displ_min ) / displ_dbin )
              displ_distr_bulk_S(k,4) = displ_distr_bulk_S(k,4) + 1.d0
              ! velocities
              do j = 1, 3
                 k = int( ( vel(i,j) - vel_min ) / vel_dbin )
                 vel_distr_bulk_S(k,j) = vel_distr_bulk_S(k,j) + 1.d0
              end do
              k = int( ( sqrt(dot_product(vel(i,:),vel(i,:))) - vel_min ) / vel_dbin )
              vel_distr_bulk_S(k,4) = vel_distr_bulk_S(k,4) + 1.d0
              ! forces
              do j = 1, 3
                 k = int( ( force(i,j) - force_min ) / force_dbin )
                 force_distr_bulk_S(k,j) = force_distr_bulk_S(k,j) + 1.d0
              end do
              k = int( ( sqrt(dot_product(force(i,:),force(i,:))) - force_min ) / force_dbin )
              force_distr_bulk_S(k,4) = force_distr_bulk_S(k,4) + 1.d0
              ! potential energies
              k = int( ( poteng(i)-poteng_eq(i) - poteng_min ) / poteng_dbin )
              poteng_distr_bulk_S(k) = poteng_distr_bulk_S(k) + 1.d0
              pe_bulk_S = pe_bulk_S + poteng(i)
              ! kinetic energies
              do j = 1, 3
                 k = int( ( kineng(i,j) - kineng_min ) / kineng_dbin )
                 kineng_distr_bulk_S(k,j) = kineng_distr_bulk_S(k,j) + 1.d0
              end do
              k = int( ( sum(kineng(i,:)) - kineng_min ) / kineng_dbin )
              kineng_distr_bulk_S(k,4) = kineng_distr_bulk_S(k,4) + 1.d0
           end if
           rmsd_bulk = rmsd_bulk + dot_product(displ,displ)
        end if
     end do
     rmsd_edge = sqrt( rmsd_edge / dble(nedge_Mo+nedge_S) )
     k = int( ( rmsd_edge - displ_min ) / displ_dbin )
     rmsd_distr_edge(k) = rmsd_distr_edge(k) + 1.d0
     pe_edge_Mo = pe_edge_Mo / dble(nedge_Mo) 
     k = int( ( pe_edge_Mo-eq_pe_edge_Mo - poteng_min ) / poteng_dbin )
     pe_distr_edge_Mo(k) = pe_distr_edge_Mo(k) + 1.d0
     pe_edge_S = pe_edge_S / dble(nedge_S) 
     k = int( ( pe_edge_S-eq_pe_edge_S - poteng_min ) / poteng_dbin )
     pe_distr_edge_S(k) = pe_distr_edge_S(k) + 1.d0
     rmsd_2nd = sqrt( rmsd_2nd / dble(n2nd_Mo+n2nd_S) )
     k = int( ( rmsd_2nd - displ_min ) / displ_dbin )
     rmsd_distr_2nd(k) = rmsd_distr_2nd(k) + 1.d0
     pe_2nd_Mo = pe_2nd_Mo / dble(n2nd_Mo) 
     k = int( ( pe_2nd_Mo-eq_pe_2nd_Mo - poteng_min ) / poteng_dbin )
     pe_distr_2nd_Mo(k) = pe_distr_2nd_Mo(k) + 1.d0
     pe_2nd_S = pe_2nd_S / dble(n2nd_S) 
     k = int( ( pe_2nd_S-eq_pe_2nd_S - poteng_min ) / poteng_dbin )
     pe_distr_2nd_S(k) = pe_distr_2nd_S(k) + 1.d0
     rmsd_3rd = sqrt( rmsd_3rd / dble(n3rd_Mo+n3rd_S) )
     k = int( ( rmsd_3rd - displ_min ) / displ_dbin )
     rmsd_distr_3rd(k) = rmsd_distr_3rd(k) + 1.d0
     pe_3rd_Mo = pe_3rd_Mo / dble(n3rd_Mo) 
     k = int( ( pe_3rd_Mo-eq_pe_3rd_Mo - poteng_min ) / poteng_dbin )
     pe_distr_3rd_Mo(k) = pe_distr_3rd_Mo(k) + 1.d0
     pe_3rd_S = pe_3rd_S / dble(n3rd_S) 
     k = int( ( pe_3rd_S-eq_pe_3rd_S - poteng_min ) / poteng_dbin )
     pe_distr_3rd_S(k) = pe_distr_3rd_S(k) + 1.d0
     if ( nbulk_Mo+nbulk_S > 0 ) rmsd_bulk = sqrt( rmsd_bulk / dble(nbulk_Mo+nbulk_S) )
     k = int( ( rmsd_bulk - displ_min ) / displ_dbin )
     rmsd_distr_bulk(k) = rmsd_distr_bulk(k) + 1.d0
     pe_bulk_Mo = pe_bulk_Mo / dble(nbulk_Mo) 
     k = int( ( pe_bulk_Mo-eq_pe_bulk_Mo - poteng_min ) / poteng_dbin )
     pe_distr_bulk_Mo(k) = pe_distr_bulk_Mo(k) + 1.d0
     pe_bulk_S = pe_bulk_S / dble(nbulk_S) 
     k = int( ( pe_bulk_S-eq_pe_bulk_S - poteng_min ) / poteng_dbin )
     pe_distr_bulk_S(k) = pe_distr_bulk_S(k) + 1.d0
     
  end do
  close(100)

  ! printout
  open ( unit=200, file='4-displacements_edge_Mo.dat' )
  open ( unit=201, file='4-displacements_edge_S.dat' )
  open ( unit=202, file='4-displacements_2nd_Mo.dat' )
  open ( unit=203, file='4-displacements_2nd_S.dat' )
  open ( unit=204, file='4-displacements_3rd_Mo.dat' )
  open ( unit=205, file='4-displacements_3rd_S.dat' )
  open ( unit=210, file='4-displacements_bulk_Mo.dat' )
  open ( unit=211, file='4-displacements_bulk_S.dat' )
  do i = 0, displ_nbin
     write(200,*) displ_min+(dble(i)+0.5d0)*displ_dbin, (displ_distr_edge_Mo(i,j)/dble(iconf*nedge_Mo),j=1,4)
     write(201,*) displ_min+(dble(i)+0.5d0)*displ_dbin, (displ_distr_edge_S(i,j) /dble(iconf*nedge_S) ,j=1,4)
     write(202,*) displ_min+(dble(i)+0.5d0)*displ_dbin, (displ_distr_2nd_Mo(i,j)/dble(iconf*n2nd_Mo),j=1,4)
     write(203,*) displ_min+(dble(i)+0.5d0)*displ_dbin, (displ_distr_2nd_S(i,j) /dble(iconf*n2nd_S) ,j=1,4)
     write(204,*) displ_min+(dble(i)+0.5d0)*displ_dbin, (displ_distr_3rd_Mo(i,j)/dble(iconf*n3rd_Mo),j=1,4)
     write(205,*) displ_min+(dble(i)+0.5d0)*displ_dbin, (displ_distr_3rd_S(i,j) /dble(iconf*n3rd_S) ,j=1,4)
     if ( nbulk_Mo > 0 ) write(210,*) displ_min+(dble(i)+0.5d0)*displ_dbin, (displ_distr_bulk_Mo(i,j)/dble(iconf*nbulk_Mo),j=1,4)
     if ( nbulk_S  > 0 ) write(211,*) displ_min+(dble(i)+0.5d0)*displ_dbin, (displ_distr_bulk_S(i,j) /dble(iconf*nbulk_S) ,j=1,4)
  end do
  do i = 0, 5
     close(200+i)
  end do
  close(210)
  close(211)
  open ( unit=200, file='4-velocities_edge_Mo.dat' )
  open ( unit=201, file='4-velocities_edge_S.dat' )
  open ( unit=202, file='4-velocities_2nd_Mo.dat' )
  open ( unit=203, file='4-velocities_2nd_S.dat' )
  open ( unit=204, file='4-velocities_3rd_Mo.dat' )
  open ( unit=205, file='4-velocities_3rd_S.dat' )
  open ( unit=210, file='4-velocities_bulk_Mo.dat' )
  open ( unit=211, file='4-velocities_bulk_S.dat' )
  do i = 0, vel_nbin
     write(200,*) vel_min+(dble(i)+0.5d0)*vel_dbin, (vel_distr_edge_Mo(i,j)/dble(iconf*nedge_Mo),j=1,4)
     write(201,*) vel_min+(dble(i)+0.5d0)*vel_dbin, (vel_distr_edge_S(i,j) /dble(iconf*nedge_S) ,j=1,4)
     write(202,*) vel_min+(dble(i)+0.5d0)*vel_dbin, (vel_distr_2nd_Mo(i,j)/dble(iconf*n2nd_Mo),j=1,4)
     write(203,*) vel_min+(dble(i)+0.5d0)*vel_dbin, (vel_distr_2nd_S(i,j) /dble(iconf*n2nd_S) ,j=1,4)
     write(204,*) vel_min+(dble(i)+0.5d0)*vel_dbin, (vel_distr_3rd_Mo(i,j)/dble(iconf*n3rd_Mo),j=1,4)
     write(205,*) vel_min+(dble(i)+0.5d0)*vel_dbin, (vel_distr_3rd_S(i,j) /dble(iconf*n3rd_S) ,j=1,4)
     if ( nbulk_Mo > 0 ) write(210,*) vel_min+(dble(i)+0.5d0)*vel_dbin, (vel_distr_bulk_Mo(i,j)/dble(iconf*nbulk_Mo),j=1,4)
     if ( nbulk_S  > 0 ) write(211,*) vel_min+(dble(i)+0.5d0)*vel_dbin, (vel_distr_bulk_S(i,j) /dble(iconf*nbulk_S) ,j=1,4)
  end do
  do i = 0, 5
     close(200+i)
  end do
  close(210)
  close(211)
  open ( unit=200, file='4-forces_edge_Mo.dat' )
  open ( unit=201, file='4-forces_edge_S.dat' )
  open ( unit=202, file='4-forces_2nd_Mo.dat' )
  open ( unit=203, file='4-forces_2nd_S.dat' )
  open ( unit=204, file='4-forces_3rd_Mo.dat' )
  open ( unit=205, file='4-forces_3rd_S.dat' )
  open ( unit=210, file='4-forces_bulk_Mo.dat' )
  open ( unit=211, file='4-forces_bulk_S.dat' )
  do i = 0, force_nbin
     write(200,*) force_min+(dble(i)+0.5d0)*force_dbin, (force_distr_edge_Mo(i,j)/dble(iconf*nedge_Mo),j=1,4)
     write(201,*) force_min+(dble(i)+0.5d0)*force_dbin, (force_distr_edge_S(i,j) /dble(iconf*nedge_S) ,j=1,4)
     write(202,*) force_min+(dble(i)+0.5d0)*force_dbin, (force_distr_2nd_Mo(i,j)/dble(iconf*n2nd_Mo),j=1,4)
     write(203,*) force_min+(dble(i)+0.5d0)*force_dbin, (force_distr_2nd_S(i,j) /dble(iconf*n2nd_S) ,j=1,4)
     write(204,*) force_min+(dble(i)+0.5d0)*force_dbin, (force_distr_3rd_Mo(i,j)/dble(iconf*n3rd_Mo),j=1,4)
     write(205,*) force_min+(dble(i)+0.5d0)*force_dbin, (force_distr_3rd_S(i,j) /dble(iconf*n3rd_S) ,j=1,4)
     if ( nbulk_Mo > 0 ) write(210,*) force_min+(dble(i)+0.5d0)*force_dbin, (force_distr_bulk_Mo(i,j)/dble(iconf*nbulk_Mo),j=1,4)
     if ( nbulk_S  > 0 ) write(211,*) force_min+(dble(i)+0.5d0)*force_dbin, (force_distr_bulk_S(i,j) /dble(iconf*nbulk_S) ,j=1,4)
  end do
  do i = 0, 5
     close(200+i)
  end do
  close(210)
  close(211)
  open ( unit=200, file='4-potential_energies_edge_Mo.dat' )
  open ( unit=201, file='4-potential_energies_edge_S.dat' )
  open ( unit=202, file='4-potential_energies_2nd_Mo.dat' )
  open ( unit=203, file='4-potential_energies_2nd_S.dat' )
  open ( unit=204, file='4-potential_energies_3rd_Mo.dat' )
  open ( unit=205, file='4-potential_energies_3rd_S.dat' )
  open ( unit=210, file='4-potential_energies_bulk_Mo.dat' )
  open ( unit=211, file='4-potential_energies_bulk_S.dat' )
  open ( unit=250, file='4-potential_energies_edge_Mo_old.dat' )
  open ( unit=251, file='4-potential_energies_edge_S_old.dat' )
  open ( unit=252, file='4-potential_energies_2nd_Mo_old.dat' )
  open ( unit=253, file='4-potential_energies_2nd_S_old.dat' )
  open ( unit=254, file='4-potential_energies_3rd_Mo_old.dat' )
  open ( unit=255, file='4-potential_energies_3rd_S_old.dat' )
  open ( unit=260, file='4-potential_energies_bulk_Mo_old.dat' )
  open ( unit=261, file='4-potential_energies_bulk_S_old.dat' )
  do i = 0, poteng_nbin
     write(200,*) poteng_min+(dble(i)+0.5d0)*poteng_dbin, poteng_distr_edge_Mo(i)/dble(iconf*nedge_Mo)
     write(201,*) poteng_min+(dble(i)+0.5d0)*poteng_dbin, poteng_distr_edge_S(i) /dble(iconf*nedge_S) 
     write(202,*) poteng_min+(dble(i)+0.5d0)*poteng_dbin, poteng_distr_2nd_Mo(i)/dble(iconf*n2nd_Mo)
     write(203,*) poteng_min+(dble(i)+0.5d0)*poteng_dbin, poteng_distr_2nd_S(i) /dble(iconf*n2nd_S) 
     write(204,*) poteng_min+(dble(i)+0.5d0)*poteng_dbin, poteng_distr_3rd_Mo(i)/dble(iconf*n3rd_Mo)
     write(205,*) poteng_min+(dble(i)+0.5d0)*poteng_dbin, poteng_distr_3rd_S(i) /dble(iconf*n3rd_S) 
     if ( nbulk_Mo > 0 ) write(210,*) poteng_min+(dble(i)+0.5d0)*poteng_dbin, poteng_distr_bulk_Mo(i)/dble(iconf*nbulk_Mo)
     if ( nbulk_S  > 0 ) write(211,*) poteng_min+(dble(i)+0.5d0)*poteng_dbin, poteng_distr_bulk_S(i) /dble(iconf*nbulk_S) 
     write(250,*) poteng_min+(dble(i)+0.5d0)*poteng_dbin, pe_distr_edge_Mo(i)/dble(iconf)
     write(251,*) poteng_min+(dble(i)+0.5d0)*poteng_dbin, pe_distr_edge_S(i) /dble(iconf) 
     write(252,*) poteng_min+(dble(i)+0.5d0)*poteng_dbin, pe_distr_2nd_Mo(i)/dble(iconf)
     write(253,*) poteng_min+(dble(i)+0.5d0)*poteng_dbin, pe_distr_2nd_S(i) /dble(iconf) 
     write(254,*) poteng_min+(dble(i)+0.5d0)*poteng_dbin, pe_distr_3rd_Mo(i)/dble(iconf)
     write(255,*) poteng_min+(dble(i)+0.5d0)*poteng_dbin, pe_distr_3rd_S(i) /dble(iconf) 
     if ( nbulk_Mo > 0 ) write(260,*) poteng_min+(dble(i)+0.5d0)*poteng_dbin, pe_distr_bulk_Mo(i)/dble(iconf)
     if ( nbulk_S  > 0 ) write(261,*) poteng_min+(dble(i)+0.5d0)*poteng_dbin, pe_distr_bulk_S(i) /dble(iconf) 
  end do
  do i = 0, 5
     close(200+i)
     close(250+i)
  end do
  close(210)
  close(211)
  close(260)
  close(261)
  open ( unit=200, file='4-kinetic_energies_edge_Mo.dat' )
  open ( unit=201, file='4-kinetic_energies_edge_S.dat' )
  open ( unit=202, file='4-kinetic_energies_2nd_Mo.dat' )
  open ( unit=203, file='4-kinetic_energies_2nd_S.dat' )
  open ( unit=204, file='4-kinetic_energies_3rd_Mo.dat' )
  open ( unit=205, file='4-kinetic_energies_3rd_S.dat' )
  open ( unit=210, file='4-kinetic_energies_bulk_Mo.dat' )
  open ( unit=211, file='4-kinetic_energies_bulk_S.dat' )
  do i = 0, kineng_nbin
     write(200,*) kineng_min+(dble(i)+0.5d0)*kineng_dbin, (kineng_distr_edge_Mo(i,j)/dble(iconf*nedge_Mo),j=1,4)
     write(201,*) kineng_min+(dble(i)+0.5d0)*kineng_dbin, (kineng_distr_edge_S(i,j) /dble(iconf*nedge_S) ,j=1,4)
     write(202,*) kineng_min+(dble(i)+0.5d0)*kineng_dbin, (kineng_distr_2nd_Mo(i,j)/dble(iconf*n2nd_Mo),j=1,4)
     write(203,*) kineng_min+(dble(i)+0.5d0)*kineng_dbin, (kineng_distr_2nd_S(i,j) /dble(iconf*n2nd_S) ,j=1,4)
     write(204,*) kineng_min+(dble(i)+0.5d0)*kineng_dbin, (kineng_distr_3rd_Mo(i,j)/dble(iconf*n3rd_Mo),j=1,4)
     write(205,*) kineng_min+(dble(i)+0.5d0)*kineng_dbin, (kineng_distr_3rd_S(i,j) /dble(iconf*n3rd_S) ,j=1,4)
     if ( nbulk_Mo > 0 ) write(210,*) kineng_min+(dble(i)+0.5d0)*kineng_dbin, (kineng_distr_bulk_Mo(i,j)/dble(iconf*nbulk_Mo),j=1,4)
     if ( nbulk_S  > 0 ) write(211,*) kineng_min+(dble(i)+0.5d0)*kineng_dbin, (kineng_distr_bulk_S(i,j) /dble(iconf*nbulk_S) ,j=1,4)
  end do
  do i = 0, 5
     close(200+i)
  end do
  close(210)
  close(211)
  open ( unit=200, file='4-rmsd_edge.dat' )
  open ( unit=201, file='4-rmsd_2nd.dat' )
  open ( unit=202, file='4-rmsd_3rd.dat' )
  open ( unit=210, file='4-rmsd_bulk.dat' )
  do i = 0, displ_nbin
     write(200,*) displ_min+(dble(i)+0.5d0)*displ_dbin, rmsd_distr_edge(i)/dble(iconf)
     write(201,*) displ_min+(dble(i)+0.5d0)*displ_dbin, rmsd_distr_2nd(i)/dble(iconf)
     write(202,*) displ_min+(dble(i)+0.5d0)*displ_dbin, rmsd_distr_3rd(i)/dble(iconf)
     if ( nbulk_Mo+nbulk_S > 0 ) write(210,*) displ_min+(dble(i)+0.5d0)*displ_dbin, rmsd_distr_bulk(i)/dble(iconf)
  end do
  close(200)
  close(201)
  close(202)
  close(210)

  stop
end program distributions

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
  use vars, only: atoms, tipo, tag, nbulk_Mo, nedge_Mo, nbulk_S, nedge_S, fixed, pos_eq, distance, thickness_edge, n2nd_Mo, n2nd_S, n3rd_Mo, n3rd_S
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
        else if ( minval(distance(i,:)) > thickness_edge .and. minval(distance(i,:)) <= 2.d0*thickness_edge ) then
           tag(i) = '2nd'
        else if ( minval(distance(i,:)) > 2.d0*thickness_edge .and. minval(distance(i,:)) <= 3.d0*thickness_edge ) then
           tag(i) = '3rd'
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
  nedge_Mo = 0
  nedge_S = 0
  n2nd_Mo = 0
  n2nd_S = 0
  n3rd_Mo = 0
  n3rd_S = 0
  nbulk_Mo = 0
  nbulk_S = 0
  do i = 1, atoms
     if ( tag(i) == 'edge' ) then
        if ( tipo(i) == 2 ) nedge_Mo = nedge_Mo + 1
        if ( tipo(i) == 1 .or. tipo(i) == 3 ) nedge_S = nedge_S + 1
     else if ( tag(i) == '2nd' ) then
        if ( tipo(i) == 2 ) n2nd_Mo = n2nd_Mo + 1
        if ( tipo(i) == 1 .or. tipo(i) == 3 ) n2nd_S = n2nd_S + 1
     else if ( tag(i) == '3rd' ) then
        if ( tipo(i) == 2 ) n3rd_Mo = n3rd_Mo + 1
        if ( tipo(i) == 1 .or. tipo(i) == 3 ) n3rd_S = n3rd_S + 1
     else if ( tag(i) == 'bulk' ) then
        if ( tipo(i) == 2 ) nbulk_Mo = nbulk_Mo + 1
        if ( tipo(i) == 1 .or. tipo(i) == 3 ) nbulk_S = nbulk_S + 1
     else 
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

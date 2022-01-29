program elab
  implicit none
  real(8), parameter :: eV_Ang_to_N = 1.602176634d-9 ! eV/Ang to N
  real(8), parameter :: Ang_to_m = 1.d-10            ! Ang to m
  real(8), parameter :: m2_to_nm2 = 1.d18            ! m^2 to nm^2
  real(8), parameter :: Pa_to_kPa = 1.d-3            ! Pa to kPa
  real(8), parameter :: m_to_um = 1.d6               ! m to um
  real(8), parameter :: N_m_to_nN_um = 1.d3          ! N/m to nN/um
  integer :: i, j, k, nsa, npe, imin
  real(8) :: fave, fstd, fmin, A, E, E_err, lx, p, S, S_err!, xxx
  character(20) :: excsa, excpe
  real(8), allocatable :: f(:,:), ff(:)
  
  read(*,*) nsa, npe, excsa, excpe
  if ( excsa /= 'none' .and. excsa /= 'min' ) stop 'error sample'
  if ( excpe /= 'none' .and. excpe /= 'first' ) stop 'error period'
  
  allocate ( f(nsa,npe) )
  
  open(unit=10,file='x.0')
  do i = 1, nsa
     do j = 1, npe
        read(10,*) k, f(i,j) ! eV/Ang
     end do
  end do
  close(10)

  allocate ( ff(nsa) )

  if ( excpe == 'none' ) then
     do i = 1, nsa
        ff(i) = 0.d0
        do j = 1, npe
           ff(i) = ff(i) + f(i,j)
        end do
        ff(i) = ff(i) / dble(npe)
     end do
  else if ( excpe == 'first' ) then
     do i = 1, nsa
        ff(i) = 0.d0
        do j = 2, npe
           ff(i) = ff(i) + f(i,j)
        end do
        ff(i) = ff(i) / dble(npe-1)
     end do
  end if

  deallocate ( f )
  
  fave = 0.d0
  fstd = 0.d0
  open(unit=90,file='2-considered.txt')
  if ( excsa == 'none' ) then
     do i = 1, nsa
        fave = fave + ff(i)
        fstd = fstd + ff(i)*ff(i)
        write(90,*) i
     end do
     fave = fave / dble(nsa)
     fstd = sqrt( fstd / dble(nsa) - fave*fave )
  else if ( excsa == 'min' ) then
     fmin = 1.d99
     do i = 1, nsa
        if ( ff(i) < fmin ) then
           imin = i
           fmin = ff(i)
        end if
     end do
     do i = 1, nsa
        if ( i /= imin ) then
           fave = fave + ff(i)
           fstd = fstd + ff(i)*ff(i)
           write(90,*) i
        end if
     end do
     fave = fave / dble(nsa-1)
     fstd = sqrt( fstd / dble(nsa-1) - fave*fave )
  end if
  close(90)
  fave = fave * eV_Ang_to_N     ! N
  fstd = fstd * eV_Ang_to_N     ! N

  open(unit=20,file='input')
  read(20,*) lx      ! Ang
  lx = lx * Ang_to_m ! m
  close(20)

  A = sqrt(3.d0) / 4.d0 * lx**2 ! m^2
  S = fave / A                     ! Pa
  S_err = fstd / A             ! Pa
  p = 3.d0 * lx                 ! m
  E = fave / p                     ! N/m
  E_err = fstd / p             ! N/m
  
  write(*,*) A*m2_to_nm2, S*Pa_to_kPa, S_err*Pa_to_kPa, p*m_to_um, E*N_m_to_nN_um, E_err*N_m_to_nN_um ! nm^2, kPa, um, nN/um
  
  stop
end program elab

program friction
  implicit none
  integer :: i, j, k, iconf, ny, nprint, iperiod, min
  real(8) :: t, f, f_old, xxx, W, ly, pos, pos_old, f_ave, f_std, vel, displ, offset, period, pos0
  
  read(*,*) vel, nprint, min, displ, ly, ny ! vel in Ang/ps

  ! set parameters
  period = ly / dble(ny) ! Ang
  offset = dble(min) * displ ! Ang

  ! open units
  open(unit=10,file='system.out') !time etotal pe ke temp f_tether_x f_tether_y f_tether_z f_tether_f c_f_s_x c_f_s_y c_f_s_z
  read(10,*) (xxx,j=1,10), f ! eV/Ang
  open(unit=100,file='1-work.dat')
  open(unit=300,file='1-res.txt')

  ! init
  iconf = 0
  k = 0
  iperiod = 1
  W = 0.d0
  pos = 0.d0
  f_ave = 0.d0
  f_std = 0.d0
  pos_old = pos
  f_old = -f
  pos0 = pos
  
  ! loop over steps
  do
     read(10,*,iostat=i) t, (xxx,j=1,9), f ! eV/Ang
     if ( i /= 0 ) exit
     pos = t * vel
     
     ! skip until reaching the first minimum
     if ( pos < offset ) then 
        pos_old = pos
        f_old = -f
        pos0 = pos
        cycle
     else
        
        ! reset counters and averages
        if ( pos >= dble(iperiod) * period + offset ) then
           write(300,*) iperiod, f_ave/dble(iconf), sqrt( f_std/dble(iconf) - (f_ave/dble(iconf))**2 )
           iconf = 0
           iperiod = iperiod + 1
           if ( dble(iperiod) * period + offset > 20.d0 ) exit
           f_ave = 0.d0
           f_std = 0.d0
        end if

        ! accumulate averages
        W = W + 0.5d0 * ( f_old - f ) * ( pos - pos_old ) ! eV
        iconf = iconf + 1
        k = k + 1
        f_ave = f_ave - f
        f_std = f_std + f*f
        if ( mod(k,nprint) == 0 ) write(100,*) pos-pos0, W
        
        pos_old = pos
        f_old = -f

     end if

  end do
  close(10)
  close(100)
  close(200)
  close(300)
  
  stop
end program friction

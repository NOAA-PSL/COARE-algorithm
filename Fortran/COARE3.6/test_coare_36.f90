program test_coare_36
  use machine , only : kind_phys
  use coare36_module, only: coare36flux, coare36flux_coolskin


  implicit none

  real :: zet
  integer :: ii
  ! input arguments for coare36_flux
  real(kind=kind_phys):: u, zu, t, zt, rh, zq, P, ts,    &
       sw_dn, lw_dn, grav, zi, rain, Ss, cp, latitude, dT_skinx, dz_skin
  ! input/output arguments for coare36_flux
  real(kind_phys) :: sigH
  ! output arguments for coare36_flux
  real(kind=kind_phys) :: tau, hsb, hlb, Le, &
       sw_net, lw_dn_abs, lw_up
  real(kind=kind_phys), dimension(21) :: output_vec, output_vec_python, output_vec_cskin

  open(unit=12,file='test_35_data.txt',form='FORMATTED',action='read')
  read(12,*)
  open(unit=14,file='../MATLABv3p6/test_36_output_matlab_072821.txt',form='FORMATTED',action='read')
  open(unit=15,file='test_36_output_f90coolskin.txt',form='FORMATTED',action='write')
  read(14,*)
  write(15,*) '# usr	tau	hsb	hlb	hlwebb	tsr	qsr	zot	zoq	Cd	Ch	Ce	L	zet	dter	dqer	tkt	RF	Cdn_10	Chn_10	Cen_10'
  do ii = 1,116
    write(*,*) 'Running case ', ii
    read(12,*) u, zu, t, zt, rh, zq, P, ts,    &
         sw_dn, lw_dn, latitude, zi, rain, cp, sigH
    read(14,*) output_vec_python(1:21)

    grav = 9.81
    Ss = 35.5
    call coare36flux_coolskin( &
         u,zu,t,zt,rh,zq,P,ts,sw_dn,lw_dn,grav,zi,rain,Ss,cp,sigH, &
         tau,hsb,hlb,Le,sw_net,lw_dn_abs,lw_up,dT_skinx,dz_skin,output_vec_cskin)
    write(15,'(21E14.6)') output_vec_cskin(1:21)


    ! convert units on fortran output
    output_vec_cskin(7) = 1.e3*output_vec_cskin(7)

    write(*,*) 'Values (matlab36, fortran36, fortran36CoolSkin)'
    write(*,*) '#      usr           tau           hsb           hlb           hlwebb           tsr           qsr           zot           zoq           '
    write(*,'(9E14.6)') output_vec_python(1:9)
    write(*,'(21E14.6)') output_vec_cskin(1:9)
    write(*,*)
    write(*,*) '#      Cd           Ch           Ce           L           zet           dter           dqer           tkt      '
    write(*,'(9E14.6)') output_vec_python(10:17)
    write(*,'(21E14.6)') output_vec_cskin(10:17)

    write(*,*)
    write(*,*) '#     RF           Cdn_10           Chn_10           Cen_10'
    write(*,'(9E14.6)') output_vec_python(18:21)
    write(*,'(21E14.6)') output_vec_cskin(18:21)

    write(*,*)
    write(*,*) 'Relative Error (fort36-py35, fort36CoolSkin-py35)'
    write(*,*) '#      usr           tau           hsb           hlb           hlwebb           tsr           qsr           zot           zoq           '
    write(*,'(21E14.3)') (output_vec_cskin(1:9)-output_vec_python(1:9))/output_vec_python(1:9)

    write(*,*)
    write(*,*) '#      Cd           Ch           Ce           L           zet           dter           dqer           tkt      '
    write(*,'(21E14.3)') (output_vec_cskin(10:17)-output_vec_python(10:17))/output_vec_python(10:17)

    write(*,*)
    write(*,*) '#     RF           Cdn_10           Chn_10           Cen_10'
    write(*,'(21E14.3)') (output_vec_cskin(18:21)-output_vec_python(18:21))/output_vec_python(18:21)

  end do

  close(unit=12)
  close(unit=14)
  close(unit=15)

end program test_coare_36

program main

#ifdef __NVCOMPILER
#ifdef ACC
use openacc_curand
#else
use curand
#endif ACC

#endif __NVCOMPILER

implicit none

INTEGER,PARAMETER :: rp = KIND(0.d0)
REAL(rp),PARAMETER :: C_PI = 4.0_rp*ATAN(1.0_rp) !< Definition of @f$\pi@f$
REAL(rp),PARAMETER :: C_E = 1.602176E-19_rp !< Absolute value of electron charge in Coulombs (C).
REAL(rp),PARAMETER :: C_ME = 9.109382E-31_rp !< Electron mass in kg
REAL(rp),PARAMETER :: C_C = 299792458.0_rp !< Light speed in m/s
REAL(rp), PARAMETER :: C_MU = 4.0_rp*C_PI*1E-7_rp !< Vacuum permeability in N/A^2
REAL(rp), PARAMETER :: C_E0 = 1.0_rp/(C_MU*C_C**2) !< Vacuum permittivity in C^2/(N*m^2)
INTEGER  :: c1,c2,cr
REAL(rp) :: rate
CHARACTER(100) :: path_to_inputs,path_to_outputs
INTEGER,PARAMETER 	:: default_unit_open = 101
INTEGER,PARAMETER 	:: output_write = 202,data_write = 102
INTEGER :: argn,read_stat,flag
REAL(rp)	:: dt,simulation_time,output_time
INTEGER :: nRE,iout,it,pp,num_outputs,t_steps
REAL(rp),ALLOCATABLE,DIMENSION(:) :: KE,eta
INTEGER,ALLOCATABLE,DIMENSION(:) :: flagCol
REAL(rp) :: rnd1,rnd2
REAL(rp) :: KE0,eta0,gam0,v0
REAL(rp) :: Epar,ne,Te,Zeff,KEmin
REAL(rp) :: Clog0,gammac0,vthe0,gammin,vmin,pmin,vratmin,psivratmin,CBmin,nuD
LOGICAL :: slowing_down,pitch_diffusion,energy_diffusion
REAL(rp) :: dW1,dW2,gam,vmag,pmag,xi,vrat,psivrat,CB,CA,dCA,CF,dp,dxi

#ifdef __NVCOMPILER
#ifdef ACC
type(curandStateXORWOW) :: h
integer(8) :: seed, seq, offset
#else
type(curandGenerator) :: g
INTEGER :: istat
#endif ACC
#endif __NVCOMPILER

NAMELIST /input_parameters/ nRE,simulation_time,output_time,KE0,eta0,Epar,ne,Te,Zeff,KEmin, &
slowing_down,pitch_diffusion,energy_diffusion

!! initialize system_clock
CALL system_clock(count_rate=cr)
rate = REAL(cr)

CALL system_clock(c1)

!! find input/output file
argn = command_argument_count()

call get_command_argument(1,path_to_inputs)
call get_command_argument(2,path_to_outputs)

!! open log and data output files
OPEN(UNIT=data_write, &
FILE=TRIM(path_to_outputs)//"data.korc", &
STATUS='UNKNOWN',FORM='FORMATTED',POSITION='REWIND')

OPEN(UNIT=output_write, &
FILE=TRIM(path_to_outputs)//"output.korc", &
STATUS='UNKNOWN',FORM='FORMATTED',POSITION='REWIND')

!! set defaults for inputs, open input file, read from input file

nRE=1
simulation_time=0._rp
output_time=0._rp
KE0=1.E7
eta0=10._rp
KEmin=1.E5
Epar=1._rp
ne=1.E20
Te=1.5
Zeff=1.
slowing_down=.true.
pitch_diffusion=.true.
energy_diffusion=.true.

OPEN(UNIT=default_unit_open,FILE=TRIM(path_to_inputs), &
STATUS='OLD',FORM='formatted',POSITION='REWIND')
read(default_unit_open,nml=input_parameters,IOSTAT=read_stat)
close(default_unit_open)

write(output_write,'("* * * Particle IC * * *")')

ALLOCATE(KE(nRE))
ALLOCATE(eta(nRE))
ALLOCATE(flagCol(nRE))

eta0=eta0*C_PI/180._rp
KE0=KE0*C_E
KEmin=KEmin*C_E
Te=Te*C_E

gam0=1._rp+KE0/(C_ME*C_C**2)
v0=C_C*sqrt(1._rp-1/gam0**2)

#ifdef ACC
!$acc  parallel loop
#endif ACC
do pp=1,nRE
  KE(pp)=KE0
  eta(pp)=eta0
  flagCol(pp)=1
end do
#ifdef ACC  
!$acc end parallel loop
#endif ACC

write(output_write,'("Number of electrons: ",I16)') nRE
write(output_write,'("KE0 (ev),eta0 (deg): ",E17.10,E17.10)') KE0/C_E,eta0*180._rp/C_PI
write(data_write,'("KE0 (ev),eta0 (deg): ",E17.10,E17.10)') KE0/C_E,eta0*180._rp/C_PI
flush(output_write)

write(output_write,'("* * * Timings * * *")')

Clog0=14.9_rp - LOG(1E-20_rp*ne)/2._rp + LOG(1E-3_rp*Te/C_E)
gammac0=ne*Clog0*C_E**4/(4._rp*C_PI*C_E0**2)
vthe0=sqrt(2*Te/C_ME)

gammin=1._rp+KEmin/(C_ME*C_C**2)
vmin=C_C*sqrt(1._rp-1/gammin**2)
pmin=C_ME*C_C*sqrt(gammin**2-1._rp)

vratmin=vmin/vthe0
psivratmin=0.5_rp*(ERF(vratmin) - 2.0_rp*vratmin*EXP(-vratmin**2)/SQRT(C_PI))/vratmin**2

CBmin=gammac0/(2*vmin)*(Zeff+erf(vratmin)-psivratmin+0.5_rp*(vthe0*vmin/C_C**2)**2)
nuD=2/pmin**2*CBmin

!! Set timestep to resolve relativistic gyrofrequency, simulation time and
!! number of time steps
dt = 0.05_rp/nuD

write(output_write,'("minimum inverse pitch angle scattering frequency: ",E17.10)') 1._rp/nuD

num_outputs=ceiling(simulation_time/output_time)
t_steps=ceiling(output_time/dt)
if (t_steps.eq.0) then
   dt=0._rp
else
   dt=output_time/real(t_steps)
endif

write(output_write,'("simulation time: ",E17.10)') simulation_time
write(output_write,'("number of outputs: ",I16)') num_outputs
write(output_write,'("dt: ",E17.10)') dt
write(output_write,'("t_steps: ",I16)') t_steps*num_outputs
flush(output_write)


call system_clock(c2)

write(output_write,'("Setup time: ",E17.10)') (c2-c1)/rate

write(output_write,'("* * * * * * * * * Begin FP * * * * * * * * *")')

#ifndef ACC
#ifdef __NVCOMPILER
istat = curandCreateGeneratorHost(g,CURAND_RNG_PSEUDO_XORWOW)
#else
call random_seed()
#endif __NVCOMPILER
#endif ACC

do iout=1,num_outputs

#ifdef ACC
  !$acc  parallel loop

  seed = 12345
  seq = 0
  offset = 0
  call curand_init(seed, seq, offset, h)
#endif ACC
  do pp=1,nRE

    do it=1,t_steps

#ifdef ACC
      rnd1=curand_uniform(h)
      rnd2=curand_uniform(h)
#else
#ifdef __NVCOMPILER
      istat = curandGenerateUniformDouble(g, rnd1, 1)
      istat = curandGenerateUniformDouble(g, rnd2, 1)
#else      
      call random_number(rnd1)
      call random_number(rnd2)
#endif __NVCOMPILER
#endif ACC

      flag=flagCol(pp)

      dW1 = SQRT(3*dt)*(-1+2*rnd1)
      dW2 = SQRT(3*dt)*(-1+2*rnd2)

      gam=1._rp+KE(pp)/(C_ME*C_C**2)
      vmag=C_C*sqrt(1._rp-1/gam**2)
      pmag=C_ME*C_C*sqrt(gam**2-1._rp)
      xi=cos(eta(pp))

      vrat=vmag/vthe0
      psivrat=0.5_rp*(ERF(vrat) - 2.0_rp*vrat*EXP(-vrat**2)/SQRT(C_PI))/vrat**2

      CB=gammac0/(2*vmag)*(Zeff+erf(vrat)-psivrat+0.5_rp*(vthe0*vmag/C_C**2)**2)
      CA=gammac0/vmag*psivrat
      dCA=gammac0*((2*(gam*vmag/C_C)**2-1)*psivrat+2.0_rp*vrat*EXP(-vrat**2)/SQRT(C_PI))/(gam**3*C_ME*vmag**2)
      CF=gammac0*psivrat/Te

      if (.not.slowing_down) CF=0._rp
      if (.not.pitch_diffusion) CB=0._rp
      if (.not.energy_diffusion) THEN
         CA=0._rp
         dCA=0._rp
      ENDIF

      dp=REAL(flag)*((-CF+dCA+Epar*xi)*dt+sqrt(2.0_rp*CA)*dW1)
      dxi=REAL(flag)*((-2*xi*CB/(pmag*pmag)+Epar*(1-xi*xi)/pmag)*dt- &
      sqrt(2.0_rp*CB*(1-xi*xi))/pmag*dW1)

      pmag=pmag+dp
      xi=xi+dxi

      if (xi>1) then
        xi=1-mod(xi,1._rp)
      else if (xi<-1) then
        xi=-1-mod(xi,-1._rp)
      endif

      KE(pp)=C_ME*C_C**2*(sqrt(1+(pmag/(C_ME*C_C))**2)-1._rp)
      eta(pp)=acos(xi)

      if ((KE(pp).lt.KEmin).and.(flag.eq.1)) then
        flagCol(pp)=0
      end if

    end do

    write(data_write,'("pp,KE (ev),eta (deg): ",I16,E17.10,E17.10)') pp,KE(pp)/C_E,eta(pp)*180._rp/C_PI
  end do 
end do

#ifdef __NVCOMPILER
istat = curandDestroyGenerator(g)
#endif __NVCOMPILER

write(output_write,'("* * * * * * * * * Final Conditions * * * * * * * * *")')
write(output_write,'("KE (ev),eta (deg): ",E17.10,E17.10)') KE(1)/C_E,eta(1)*180._rp/C_PI

call system_clock(c1)

write(output_write,'("FP time: ",E17.10)') (c1-c2)/rate

! * * * FINALIZING SIMULATION * * *

write(output_write,'("FP ran suppessfully!")')
close(output_write)
close(data_write)

end program main
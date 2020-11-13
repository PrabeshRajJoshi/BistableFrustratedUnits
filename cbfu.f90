! program is compiled with the command >gfortran cbfu.f90 −o cbfu.out
! executed with > ./cbfu.out options
! options - option_symbol option_value
! default values are used if symbol not called in the execution command
! symbol    description             default value
! nx        number of rows          10
! ny        number of columns       10
! bR        beta_R                  0.1d0
! sA        seed for initial        1
!           condt. of protein A 
! sB        seed for intiail        2
!           condt. of protein b
! rt        real time               50000.d0
! dt        integration step        0.01d0
! (resulting number of iterations is rt/dt = 5000000 by default)
! ws        writing step            10.d0
!           (in real time units,
!           not iteration time) 
! a0        alpha at initial time   10.d0
! a1        value by which alpha    15.d0
!           is increased
! sa        "speed" of alpha        4000.d0
! fn        name of data file       'nx'x'ny'_yyyymmddhhmm

! alpha changes liner in time by formula 
!   alpha = alpha_0 + alpha_1 * (t/alpha_speed)
! t is the iteration step and alpha_speed = sa/dt, so every
!   alpha_speed steps alpha is increased by the value of alpha_1
! for adiabatic change a1 should be maximum value of alpha we want to
!   reach and 'sa' should be equal to 'dt'.
! for constant value of alpha a0 should be equal to that value and
!   'sa' should be equal to 0.

! only nx, ny, sA, and sB are integers, all other values should be
!   written as xxxx.xxxd0

module calculations

    implicit none
    !! system size
    integer :: nx, ny, latt_size
    !! integration parameters
    real(8) :: beta_R, alpha1, alpha0, speedOFalpha, dt, real_time, s
    !! model parameters
    real(8), parameter  :: sixth = 1.d0/6.d0, K = 0.02d0, b = 0.01d0, gamma = 0.01d0
    !! initial conditions
    !!  r - random initial conditions
    !!  d - diagonal initial conditions
    !!  e - equal initial conditions
    character(1) :: ic = 'r'
    !! seeds for random number generator
    integer :: idum1 = -1, idum2 = 2
    character(22) :: fname ! fname is now a global variable
    public :: f, ran1

    contains
        subroutine calc

            real(8), allocatable, dimension(:) :: A, B, A1, B1, x, y
            real(8), allocatable, dimension(:) :: m1, m2, m3, m4, l1, l2, l3, l4
            integer :: time, step, t, check_step, idum1, idum2
            integer :: site, i, j, alpha_step
            character :: ft*5, fm*20
            real(8) :: alpha

            allocate( A(1:latt_size), B(1:latt_size), A1(1:latt_size), B1(1:latt_size))
            allocate( y(1:latt_size), x(1:latt_size))
            allocate( m1(1:latt_size), m2(1:latt_size), m3(1:latt_size), m4(1:latt_size))
            allocate( l1(1:latt_size), l2(1:latt_size), l3(1:latt_size), l4(1:latt_size))

            call initial_conditioins(A,B)

            ! integrate in time
            time = int(real_time/dt)
            step = int(s/dt)
            check_step = step
            alpha_step = int(speedOFalpha/dt)

            ! writing all parameter values in a file
            call writing

            ! write intial conditions
            open(unit=2, file=fname//"_A.dat", access = "sequential")
            !! open(unit=2, file=fname//"_B.dat", access = "sequential")

            write(ft, '(i5.5)') latt_size
            fm = '(' // ft // '(F8.4, 2X))'

            write(2, fm) A
            !! write(3, fm) B


            !! start Runge Kutta 4
            do t=0, time-1
                if( mod(t,alpha_step) == 0 ) then
                    alpha = alpha0 + alpha1*dfloat(t)/dfloat(alpha_step)
                endif

                x = A
                y = B
                do i = 1, nx
                    do j = 1, ny
                        site = (i-1)*ny + j
                        m1(site) = dt * f(x, y, 'a', i, j, site, alpha)
                        l1(site) = dt * f(x, y, 'b', i, j, site, alpha)
                    enddo
                enddo

                x = A + 0.5*m1
                y = B + 0.5*l1
                do i = 1, nx
                    do j = 1, ny
                        site = (i-1)*ny + j
                        m2(site) = dt * f(x, y, 'a', i, j, site, alpha)
                        l2(site) = dt * f(x, y, 'b', i, j, site, alpha)
                    enddo
                enddo

                x = A + 0.5*m2
                y = B + 0.5*l2
                do i = 1, nx
                    do j = 1, ny
                        site = (i-1)*ny + j
                        m3(site) = dt * f(x, y, 'a', i, j, site, alpha)
                        l3(site) = dt * f(x, y, 'b', i, j, site, alpha)
                    enddo
                enddo

                x = A + 0.5*m3
                y = B + 0.5*l3
                do i = 1, nx
                    do j = 1, ny
                        site = (i-1)*ny + j
                        m4(site) = dt * f(x, y, 'a', i, j, site, alpha)
                        l4(site) = dt * f(x, y, 'b', i, j, site, alpha)
                    enddo
                enddo

                do i = 1, latt_size
                    A1(i) = A(i) + sixth * ( m1(i) + 2.d0*m2(i) + 2.d0*m3(i) + m4(i) )
                    B1(i) = B(i) + sixth * ( l1(i) + 2.d0*l2(i) + 2.d0*l3(i) + l4(i) )
                enddo

                A = A1
                B = B1

                if (t == check_step) then
                    write(2,fm) A
                    ! write(3,fm) B
                    check_step = check_step + step
                endif

            enddo

            close(2)
            !close(3)
            
            deallocate(A, A1, B1, m1, m2, m3, m4, l1, l2, l3, l4, x, y)

        end subroutine calc

        ! --------------------------------------------------------

        subroutine initial_conditioins(A,B)
            integer :: i, j, site, r1, r2
            real(8), dimension(1:latt_size), intent(out) :: A,B

            select case(ic)
            case('r')
                r1 = idum1
                r2 = idum2
                do i = 1, latt_size
                    A(i) = ran1(r1)
                    B(i) = ran1(r2)
                enddo

            case('d')
                do i = 1, nx-1, 2
                    do j = 1, ny-1, 2
                        site = (i-1)*ny + j
                        A(site) = 1.d0
                        B(site) = 1.d0
                        site = i*ny + j + 1
                        A(site) = 1.d0
                        B(site) = 1.d0
                    enddo
                enddo

            case('e')
                do i = 1, latt_size
                    A(i) = 0.4d0
                    B(i) = 0.4d0
                enddo
                ! one oscillator kicked out
                A(3) = 0.5
            end select
        end subroutine
        ! --------------------------------------------------------

        ! write parameters values
        subroutine writing
            open(unit = 1, file = fname // "_param.dat", access = "sequential")
            write(1,*) "parameter values for the data in the files called cx_cy_*.dat"
            write(1, *) "************************************************************"
            write( 1 , '(2 (A4, 1X, I3 , 6X) )' ) "nx=", nx , "ny=", ny
            write( 1 , '(A4, F5.3 ,3X, A6, F8.0 )' ) " dt=", dt , "time=", real_time
            write( 1 , '(A14, 1X, F5.0, 1X, A10 )' ) " writing every ", s , " time units "
            write( 1 , * ) " initial conditions = " , ic
            write( 1 , '( 2 (A7, 1X, I3 , 6X))' ) " seedA=", idum1 , " seedB=" , idum2
            write( 1 , * ) " ******************* model parameters ************************** "
            write( 1 , '( 3 (A7, 1X, F5.3, 3X))' ) "gamma=", gamma, " K=" ,K, " b=" , b
            write( 1 , '(A8, 1X, F6.3) ' ) " beta_R=" , beta_R
            write( 1 , * ) " alpha = ", alpha0 , " + " , alpha1, " * ( time / " , speedOFalpha , " ) "
            close( 1 )
        end subroutine
        ! --------------------------------------------------------

        real(8) function f(x, y, var, i, j, site, alpha)
            integer, intent(in) :: site, i, j
            real(8), intent(in) :: alpha
            character, intent(in):: var*1
            real(8), dimension(1:latt_size), intent(in) :: x,y

            real(8) :: sum_repr
            integer :: site_r, lattice_site

            select case(var)
            case('a')
                if( mod(i,2)==1 .and. mod(j,2)==1 ) then
                    lattice_site = 1
                elseif( mod(i,2)==1 .and. mod(j,2)==0 ) then
                    lattice_site = 2
                elseif( mod(i,2)==0 .and. mod(j,2)==1 ) then
                    lattice_site = 2
                elseif( mod(i,2)==0 .and. mod(j,2)==0 ) then
                    lattice_site = 1
                else
                    print*, 'there has been a mistake with choosing the partiy of a lattice site'
                    stop
                endif
                
                sum_repr = 0.d0
                select case(lattice_site)
                case(1) ! odd row - odd column and even row - even column
                    site_r = mod( (nx - mod( nx - i + 2, nx )), nx ) * ny + j ! site up
                    sum_repr = 1.d0 / ( 1.d0 + ( x(site_r)/K )**2 ) + sum_repr
                    site_r = mod( i , nx )*ny + j !site down
                    sum_repr = 1.d0 / ( 1.d0 + ( x(site_r)/K )**2 ) + sum_repr

                case(2) ! odd row -even column
                    site_r = ( i - 1)*ny + ny - mod( ny - j + 1,ny ) ! site left
                    sum_repr = 1.d0 / ( 1.d0 + ( x( site_r ) /K)**2 ) + sum_repr
                    site_r = ( i - 1)*ny + mod( j , ny ) + 1 ! site right
                    sum_repr = 1.d0 / ( 1.d0 + ( x ( site_r ) /K)**2 ) + sum_repr
                endselect

                f = alpha / ( 1.d0 + ( y (site) / K) ) * ( b + x (site)**2 ) / ( 1.d0 + x (site)**2 ) &
                    & - x(site) + beta_R * sum_repr !+ beta_A * sum_active

            case('b')
                f = gamma * ( x(site) - y(site) )

            end select

        end function f
        ! --------------------------------------------------------

        real(8) function ran1(idum)
            integer :: idum, ia, im, iq, ir, ntab, ndiv
            REAL( 8 ) :: am, eps, rnmx
            PARAMETER ( ia =16807 ,im=2147483647 ,am=1./im , iq =127773 , ir=2836, &
                & ntab=32, ndiv=1+(im-1)/ntab, eps=1.2e-7,rnmx=1.- eps )
            INTEGER :: j, k, iv(ntab), iy
            SAVE iv, iy
            DATA iv /ntab*0/, iy /0/
        
            if ( ( idum == 0 ) .OR. ( iy==0) ) then
                idum = MAX(-idum , 1 )
                DO j = ntab +8,1,-1
                    k = idum/ iq
                    idum = ia *(idum - k* iq )-ir *k
                    IF ( idum < 0 ) idum = idum+im
                    IF ( j <= ntab ) iv(j) = idum
                ENDDO
                iy = iv( 1 )
            endif
            k = idum/iq
            idum = ia *(idum - k*iq ) - ir*k
            if ( idum < 0 ) idum = idum+IM
            j = 1+iy / ndiv
            iy = iv ( j )
            iv( j ) = idum
            ran1 = MIN(am* iy , rnmx )
            
        end function ran1

end module
! ***********************************************************

program cbfu
    use calculations
    implicit none
    integer :: beginning, rate, en, cn, i
    character(len=22) :: arg
    character :: cx*3 , cy*3 , broj*12, da*8 ,ti*10 ! new
    call system_clock(beginning, rate)
    ! set the default values for the program variables

    nx = 10
    ny = 10
    beta_R = 0.1d0
    idum1 = 1
    idum2 = 2
    real_time = 50000.d0
    dt = 0.01d0
    s = 10.d0
    alpha0 = 10.d0
    alpha1 = 15.d0
    speedOFalpha = 4000.d0

    ! this block will give the name to the fila as before, if you give it in the  command line it will be overwritten
    ! *************************************
    write ( cx , ' ( i3.3 ) ' ) nx
    write ( cy , ' ( i3.3 ) ' ) ny
    call date_and_time ( da , ti )
    write ( broj, ' (A8, A4) ' ) da , ti
    fname = cx // "x" // cy // "_" // broj
    ! *************************************

    cn = command_argument_count()

    open( unit = 1, file = "input.dat", access = "sequential")
        do i=2, cn, 2
            call getarg(i,arg)
            write ( 1 , ' (a22) ' ) arg
        enddo
    close(1)
    
    open( unit=2, file="input.dat" , status="old" , action="read" , access="sequential" )
        do i=1 , cn-1, 2
            call getarg(i, arg)
            
            if(arg == "nx") then
                read(2, '(i3)') nx
            elseif(arg == "ny") then
                read(2, '(i3)') ny
            elseif(arg == "bR") then
                read(2, '(f10.5)') beta_R
            elseif(arg == "sA") then
                read(2, '(i5)') idum1
            elseif(arg == "sB") then
                read(2, '(i5)') idum2
            elseif(arg == "rt") then
                read(2, '(f10.5)') real_time
            elseif(arg == "dt") then
                read(2, '(f10.5)') dt
            elseif(arg == "ws") then
                read(2, '(f10.5)') s
            elseif(arg == "a0") then
                read(2, '(f10.5)') alpha0
            elseif(arg == "a1") then
                read(2, '(f10.5)') alpha1
            elseif(arg == "sa") then
                read(2, '(f10.5)') speedOFalpha
            elseif(arg == "fn") then
                read(2, '(a22)') fname
            endif
        enddo
    close(2)

    latt_size = nx*ny

    call calc

    call system_clock(en)
    print*, real(en-beginning)/real(rate)

    end program



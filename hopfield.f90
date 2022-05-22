! Simulates a classical ('84) continuous time (Cohen-Grossberg-)Hopfield
! network, using the simple additive model formulation described in equation 8
! of Grossberg (1988). The energy function used is that from Cohen and
! Grossberg (1983).
!
! Sources:
! Cohen, M. A., & Grossberg, S. (1983). Absolute stability of global pattern
!   formation and parallel memory storage by competitive neural networks. IEEE
!   transactions on systems, man, and cybernetics, (5), 815-826.
! Grossberg, S. (1988). Nonlinear neural networks: Principles, mechanisms and
!   architectures. Neural networks, 1(1), 17-61.
! Hopfield, J. J. (1984). Neurons with graded response have collective
!   computational properties like those of two-state neurons. Proceedings of the
!   national academy of sciences, 81(10), 3088-3092.

program hopfield
    implicit none

    ! Constants and hyperparameters
    logical, parameter :: debug    = .true.  ! debug flag
    integer, parameter :: n        = 20      ! number of neurons
    real, parameter    :: duration = 5.      ! (fixed) simulation duration
    real, parameter    :: dt       = 0.1     ! (fixed) step size
    real, parameter    :: lambda   = 3.      ! steepness parameter

    ! Model parameters
    real               :: x(n)               ! short-term memory (STM) traces
    real               :: A(n)               ! exponential decay rates
    real               :: z(n, n)            ! long-term memory (LTM) traces

    ! Misc variables
    real               :: p1(n), p2(n)       ! patterns to be stored
    integer            :: i, j               ! loop variables

    ! Set the patterns to store
    p1 = (/ -1., -1., -1., -1., -1.,  1., -1., -1., -1., -1.,  1.,  1.,  1., &
             1., -1., -1.,  1., -1., -1.,  1. /)
    p2 = (/ -1., -1.,  1.,  1.,  1., -1., -1.,  1.,  1., -1., -1., -1., -1., &
             1.,  1., -1., -1.,  1., -1., -1. /)

    ! Hardcoded query / initial state (Hamming distance of 2 compared to p2)
    x  = (/  1., -1.,  1.,  1.,  1., -1., -1.,  1.,  1., -1., -1., -1., -1., &
             1.,  1., -1., -1.,  1., -1.,  1. /)

    ! Initialize model parameters
    A = 1.
    do concurrent (i = 1:n)
        do concurrent (j = 1:n)
            z(i, j) = p1(i) * p1(j) + p2(i) * p2(j)
        end do
    end do

    ! Print info 
    print *, new_line('a'), 'Input pattern:'
    write (*,'(20f5.1)') x
    if (debug) then
        print *, new_line('a'), 'Debug log (set debug to .false. to turn off):'
    end if

    ! Start simulation
    call simulation

    ! Print result
    print *, new_line('a'), 'Recovered pattern:'
    write (*,'(20f5.1)') tanh(x)

contains

    ! Driver for Runge-Kutta using a fixed stepsize (dt)
    subroutine simulation
        real    :: t              ! time variable
        integer :: nsteps         ! number of steps
        integer :: j              ! loop variable

        ! Initialization
        nsteps = ceiling(duration / dt)
        t = 0
        j = 0

        ! Run simulation
        do
            if (debug) then
                ! Print diagnostics
                print *, j, t, energy(), x
            end if

            ! Check no. of iterations
            if (j == nsteps) then
                exit
            end if

            ! Update state for next Runge-Kutta step
            t = t + dt
            j = j + 1

            ! Take next Runge-Kutta step
            call rk
        end do
    end subroutine simulation

    ! Computes a single Runge-Kutta step for a fixed stepsize dt
    subroutine rk
        real :: k(n, 4) ! slopes

        call dxdt(k(:, 1), x)
        call dxdt(k(:, 2), x + 0.5 * dt * k(:, 1))
        call dxdt(k(:, 3), x + 0.5 * dt * k(:, 2))
        call dxdt(k(:, 4), x + dt * k(:, 3))
        x = x + (dt / 6) * (k(:, 1) + 2 * k(:, 2) + 2 * k(:, 3) + k(:, 4))
    end subroutine rk

    subroutine dxdt(k, x)
        real, intent(out) :: k(n)    ! slopes
        real, intent(in)  :: x(n)    ! STM traces
        integer           :: i       ! loop variable

        do concurrent (i = 1:n)
            k(i) = -A(i) * x(i) + sum(tanh(lambda * x) * z(:, i))
        end do
    end subroutine dxdt

    ! Computes the energy (Cohen & Grossberg (1983)) using Simpson's 1/3 rule
    function energy()
        real    :: fx(n)     ! mean firing rates / activations
        real    :: energy    ! output
        integer :: i         ! loop variable

        ! Precalculate activations and derivatives
        fx = tanh(lambda * x)

        ! Calculate energy
        energy = 0.0
        do i = 1, n
            energy = energy + fx(i) * sum(fx * z(i,:))
        end do
        energy = -0.5 * energy + (lambda / 6.) &
            * sum(A * x**2 * (2. - tanh(lambda*0.5*x)**2 - tanh(lambda*x)**2))
    end function energy

end program hopfield

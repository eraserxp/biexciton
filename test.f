      program test
        ! using lapack95
        use mkl95_precision, only: WP => DP
        use mkl95_lapack !, only: GGEV
        use mkl95_blas !, only: DOT
        use f77_subroutine
        use unit_conversion
        use exciton
        use find_exciton_pairs
 
        implicit none
        interface 
          subroutine probability_dis2(time, eig_vector_matrix, eig_value, initial_coefficients, 
     &                                coefficients_squared_k, K_n)
            double precision :: time
            double precision :: eig_vector_matrix(:,:)
            double precision :: eig_value(:)
            double precision :: initial_coefficients(:)
            double precision :: coefficients_squared_k(:)
            integer :: K_n
          end subroutine
        end interface
        double precision, parameter :: Pi=3.141592653589793115997963468d0
        character(len=30) :: file_name, file_name2
        character(len=4) :: number_string       
        double precision ::  time
        double precision :: electric_field
        double precision :: delta_E
        ! exclude N=2 exciton, only N=1 exciton pairs are included
        double precision, allocatable :: ham(:,:) 
        double precision, allocatable :: eig_vector_n1(:,:)
        double precision, allocatable :: eig_vector_full_n1(:,:)
        double precision, allocatable :: B(:,:)
!        double precision, allocatable :: probability(:)
        double precision, allocatable :: xx_energy(:) ! the energy of two exciton (including biexciton state)
        double precision, allocatable :: E_continuum(:) ! the energy of continuum states
        double precision, allocatable :: alphar(:)
        double precision, allocatable :: alphai(:)
        double precision, allocatable :: beta(:)

        double precision, allocatable :: evolution_matrix(:,:) 
        double precision, allocatable :: eig_vector_matrix(:,:)
        double precision, allocatable :: eig_value(:) 
        double precision, allocatable :: initial_coefficients(:)
        double precision, allocatable :: coefficients_squared_k(:)
        integer ::  i, j, K, matrixdim
        integer :: ndim
        double precision :: wave_vector
        double precision :: theta ! the angle between electric field and the intermolecular axis
        double precision :: k_biexciton
        integer :: xxindex, imin, imax
        integer :: info
        common /energy_difference/ delta_E

        call system("rm -rf *.dat")
        theta = pi/2 !0.D0

        electric_field = 687629.834D0 
        open(unit=18, file='electric_field.dat')
        write(18,*) "electric field = ", electric_field/1.D5, "kV/cm"
        close(18)
        call Field_in_atomic_unit(electric_field)

        Number_of_molecule = 501 

          ! the energy difference: 2*|N=1,M=0> - |N=2,M=0>
          ! define the energy of dressed state |N=2,M=0> as zero
          delta_E = delta_energy_ex_angle(electric_field, number_of_molecule,theta)

        do K_n = -Number_of_molecule/2, Number_of_molecule/2 ! the definition of K_n comes from find_exciton_pairs 
          ! the wavevector of N=2 exciton is K_n
!          K_n = 0
!          write(17,*) "Number of molecule = ", Number_of_molecule
          k_biexciton = K_n*pi/( Number_of_molecule/2) ! the quantum number of biexciton (k1 + k2)
          call two_exciton_pairs ! find all possible N=1 excitons pairs
          call assign_values ! store wavevector information about N=1 exciton pairs in k1_array, k2_array
          write(*,*) "K_n:", K_n, ",       ", "Number of exciton pairs:", pair_number !it's defined in find_exciton_pairs 

 
          matrixdim = pair_number 
          ndim = pair_number -1
          allocate(ham(ndim,ndim))
          allocate(eig_vector_n1(ndim,ndim))
          allocate(eig_vector_full_n1(ndim + 1,ndim + 1)) 
          allocate(xx_energy(ndim) )
          allocate(E_continuum(ndim-1))
                   
          allocate(B(ndim,ndim))
          allocate(alphar(ndim))
          allocate(alphai(ndim))
          allocate(beta(ndim))

          call  biexciton_ham_K(Number_of_molecule, electric_field, theta, ham, pair_number, k1_array, k2_array,K_n) 
          B = 0.D0
          do i = 1, ndim
            B(i,i) = 1.D0
          end do
          call ggev(ham, B, alphar, alphai, beta, VR=eig_vector_n1,info=info)
          if (info/=0) then
            write(*,*) "Subroutine ggev fails!"
          endif
          do i = 1, ndim
            xx_energy(i) = alphar(i)/beta(i)
          end do
          eig_vector_full_n1(1:ndim, 1:ndim) = eig_vector_n1(1:ndim, 1:ndim)
          eig_vector_full_n1(1:ndim+1, ndim + 1) = 1.D0
          do i = 1, ndim
            eig_vector_full_n1(ndim+1, i) = -SUM( eig_vector_n1(:,i) )
          end do 
          call normalize_eigvector_matrix(eig_vector_full_n1, ndim + 1)         
          xx_energy = xx_energy + delta_E

          ! find out which state corresponds to biexciton
          if (theta > DACOS(1.D0/SQRT(3.D0)) ) then
            xxindex = MAXLOC(xx_energy,dim=1)
            call kick_one_out(xx_energy, ndim, xxindex, E_continuum)
            xxindex = xxindex + 1
          else
            xxindex = MINLOC(xx_energy,dim=1) 
            call kick_one_out(xx_energy, ndim, xxindex, E_continuum)
            xxindex = xxindex + 1
          endif
          write(*,*) "xxindex = ", xxindex
!----------record the biexciton spectrum----------------------------
          open(unit=23, file='biexciton_spectrum.dat', position='append')
          write(23,*) k_biexciton, xx_energy(xxindex-1)
          close(23)
!-------------------------------------------------------------------
!----------record the exciton pair spectrum----------------------------
          open(unit=33, file='exciton_pair_spectrum.dat', position='append')
          do i=1, ndim
            write(33,*) k_biexciton, xx_energy(i)
          end do
          write(33,*) "#====================================="
          close(33)
!-------------------------------------------------------------------

!-----------record the max and min energies of continuum states----------
          open(unit=44, file='Emax_continuum.dat', position='append')
          write(44,*) k_biexciton, MAXVAL(E_continuum)
          close(44)
          open(unit=55, file='Emin_continuum.dat', position='append')
          write(55,*) k_biexciton, MINVAL(E_continuum)
          close(55)
!---------------------------------------------------------------------------
          deallocate(E_continuum)


          deallocate(ham, B, alphar, alphai, beta, eig_vector_n1) 


          allocate( evolution_matrix(matrixdim,matrixdim) )
          allocate( eig_vector_matrix(matrixdim,matrixdim) )
          allocate( eig_value(matrixdim) )
          allocate( initial_coefficients(matrixdim) )
          allocate( coefficients_squared_k(matrixdim) )

          call form_evolution_matrix_kinematic(electric_field, theta, 
     &                            xx_energy,eig_vector_full_n1, evolution_matrix) 
          deallocate(eig_vector_full_n1, xx_energy)

          ! record the diagonal terms (regard M(1,1) as 0)
          open(unit=97, file='energy_difference.dat', position='append')
          write(97,*) K_n, ( evolution_matrix(xxindex,xxindex)-evolution_matrix(1,1) )/(2000*pi)           
          close(97)
          
          ! record the matrix element for transition from N=2 exciton to biexciton
          open(unit=98, file='tran_mat_N2_biex.dat', position='append')
          write(98,*) K_n, evolution_matrix(1,xxindex)/(2000*pi)           
          close(98)

          call  lapack_eig2(evolution_matrix, matrixdim, eig_value, eig_vector_matrix)
          ! after the above call, the Upper triangle part (including the diagonal elements) of evolution_matrix are destoried
          deallocate(evolution_matrix)
          write(88,*) eig_value
          initial_coefficients=0.D0
          initial_coefficients(1)=1.D0 ! start from a N=2 exciton


          do i = 0, 10000
            time = i*1.D-7
            write(*,*) i
            call probability_dis2(time, eig_vector_matrix, eig_value, initial_coefficients, coefficients_squared_k, K_n) 
            write(*,*) "===================="

            write(number_string, "(i4)") K_n 
! ADJUSTL(STRING) will left adjust a string by removing leading spaces. 
! Spaces are inserted at the end of the string as needed. 
            number_string = adjustl(number_string)
!-----------------record probability of N=2 exciton as a function of time------
            file_name = "n2ex" // trim(number_string) // ".dat"
            file_name = trim(file_name)
            open(unit=10,file=file_name,position="append")   
            write(10,*)  time, coefficients_squared_k(1)
            close(10)
!-----------------record probability of biexciton as a function of time------
            file_name2 = "biex" // trim(number_string) // ".dat"
            file_name2 = trim(file_name2)
            open(unit=20,file=file_name2,position="append")   
            write(20,*)  time, coefficients_squared_k(xxindex)
            close(20)           
          end do

          call deallocation ! deallocate k1_array and k2_array
          deallocate( eig_vector_matrix )
          deallocate( eig_value )
          deallocate( initial_coefficients )
          deallocate( coefficients_squared_k )


       end do


      end program

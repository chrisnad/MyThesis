!----------------------------------------------------------------------------------------------------------
! MODULE: mumps
!
!> @author JL Tailhan
!> @author J. Foulliaron (version 1.0 - 2011)
!
!> @brief
!> Module permettant l'utilisation de la librairie MUMPS
!----------------------------------------------------------------------------------------------------------
module mumps


contains


!------------------------------------------------------------------------------------------------------
!> @author JL Tailhan
!> @author J. Foulliaron (version 1.0 - 2011)
!
!> @brief
!> Procedure de resolution de deux système à l'aide de la librairie MUMPS \n
!>  matg*sol1 = res1 et matg*sol2 = res2
!
!> @details
!> #### DESCRIPTION:
!>  Entrees :\n
!>   - matg : matrice du systeme lineaire\n
!>   - res1 : second membre du premier systeme\n
!>   - res2 : second membre du second systeme\n\n
!
!>  Sorties :\n
!>   - sol1 : solution du premier systeme\n
!>   - sol2 : solution du secon systeme\n
!------------------------------------------------------------------------------------------------------
subroutine Solve_mumps(matg, res1,res2, sol1, sol2)

    !**************************************************************************************
    !* DEFINITIONS
    !**************************************************************************************

    !--------------------------------------------------------------------------------------
    !--- Modules generaux
    use variables, only : ndlt
    use sparse_matrix

    implicit none

    !include 'mpif.h'
    include 'dmumps_struc.h'

    !--------------------------------------------------------------------------------------
    !--- Variables in/out
    type(sparse), intent(in) :: matg
    real*8, dimension(:), intent(in) :: res1, res2
    real*8, dimension(:), intent(out) :: sol1, sol2

    !--------------------------------------------------------------------------------------
    !--- Variables locales
    type (DMUMPS_STRUC) :: id
    integer :: i, j

    !**************************************************************************************
    !* CORPS DE PROCEDURE
    !**************************************************************************************

    !--------------------------------------------------------------------------------------
    !--- Definition du communicateur
    !id%COMM = mpi_comm_world

    !--------------------------------------------------------------------------------------
    !--- Symmetrie
    id%SYM = 0

    !--------------------------------------------------------------------------------------
    !--- Travail du processus maitre
    id%PAR = 1

    id%ICNTL(5) = 0 ! La matrice est assemblee

    !--------------------------------------------------------------------------------------
    !--- Initialisation
    id%JOB = -1
    call DMUMPS(id)

    !--------------------------------------------------------------------------------------
    !--- Parametre de control
    id%ICNTL(2)=-1
    id%ICNTL(3)=-1
    id%ICNTL(4) = 1
    id%ICNTL(5) = 0 ! La matrice est assemble
    id%ICNTL(7) = 4 ! Modif Josquin 09/06/2011

    !--------------------------------------------------------------------------------------
    !--- Definition du probleme dans le processus 0
    if(id%MYID == 0) then
        id%N = ndlt
        id%NZ = matg%ia(matg%n+1)-1
        id%NRHS = 2
        id%LRHS = id%N
        allocate(id%IRN(id%NZ)) ; id%IRN = 0
        allocate(id%JCN(id%NZ)) ; id%JCN = 0
        allocate(id%A(id%NZ))   ; id%A = 0.D0
        allocate(id%RHS(id%N*2)); id%RHS = 0.D0

        id%JCN = (/(matg%ja(i), i =1,id%NZ)/)
        id%A = (/(matg%a(i), i =1,id%NZ)/)
        id%RHS(1:id%N)= (/(res1(i), i =1,id%N)/)
        id%RHS(id%N+1:id%N*2) = (/(res2(i), i =1,id%N)/)


        do i = 1,matg%n
              do j = matg%ia(i), matg%ia(i+1)-1
                id%IRN(j) = i
            enddo
        enddo

    endif

    !--------------------------------------------------------------------------------------
    !--- Resolution du systeme
    id%JOB = 6

    call DMUMPS(id)

    !--------------------------------------------------------------------------------------
    !--- Deallocation utilisateur
    if(id%MYID == 0) then
        sol1 = (/(id%RHS(i), i =1, id%N)/)
        sol2 = (/(id%RHS(i), i =id%N+1, id%N*2)/)
        deallocate (id%IRN)
        deallocate (id%JCN)
        deallocate (id%A)
        deallocate (id%RHS)
    endif

    !--------------------------------------------------------------------------------------
    !--- Deallocation interne
    id%JOB = -2
    call DMUMPS(id)

end subroutine Solve_mumps


!------------------------------------------------------------------------------------------------------
!> @author JL Tailhan
!> @author J. Foulliaron (version 1.0 - 2011)
!
!> @brief
!> Procedure de resolution de deux système à l'aide de la librairie MUMPS \n
!>  matg*sol1 = res1 et matg*sol2 = res2 \n\n
!
!> => Meme routine que Solve_mumps, mais adaptée au format CSRC \n
!> NB : Cette routine considere la matrice du systeme comme symetrique
!
!> @details
!> #### DESCRIPTION:
!>  Entrees :\n
!>   - matg : matrice du systeme lineaire\n
!>   - res1 : second membre du premier systeme\n
!>   - res2 : second membre du second systeme\n\n
!
!>  Sorties :\n
!>   - sol1 : solution du premier systeme\n
!>   - sol2 : solution du secon systeme\n
!------------------------------------------------------------------------------------------------------
subroutine MUMPS_CSRC(matg, res1, res2, sol1, sol2)

    !**************************************************************************************
    !* DEFINITIONS
    !**************************************************************************************

    !--------------------------------------------------------------------------------------
    !--- Modules generaux
    use variables, only : ndlt
    use sparse_matrix
    use formatCSRC

    implicit none

    !include 'mpif.h'
    include 'dmumps_struc.h'

    !--------------------------------------------------------------------------------------
    !--- Variables in/out
    type(matrice_CSRC), intent(in) :: matg
    real*8, dimension(:), intent(in) :: res1, res2
    real*8, dimension(:), intent(out) :: sol1, sol2

    !--------------------------------------------------------------------------------------
    !--- Variables locales
    type (DMUMPS_STRUC) :: id
    integer :: i, k
    integer(kind=4)::pos

    !real*8 :: time1, time2

    !**************************************************************************************
    !* CORPS DE PROCEDURE
    !**************************************************************************************

    !--------------------------------------------------------------------------------------
    !--- Definition du communicateur
    !id%COMM = mpi_comm_world
    id%COMM = 0

    !--------------------------------------------------------------------------------------
    !--- Symmetrie  (JLT 05/2014) remarque : id%SYM = 1 pour pb 3D est tres couteux en temps de calcul)
    id%SYM = 1

    !--------------------------------------------------------------------------------------
    !--- Travail du processus maitre
    id%PAR = 1

    !--------------------------------------------------------------------------------------
    !--- Initialisation
    id%JOB = -1
    call DMUMPS(id)

    !--------------------------------------------------------------------------------------
    !--- Parametres de controle

    !--------------------------------------------------------------------------------------
    !--- Type d'envoi de la matrice : assemblee + centralisee
    id%ICNTL(5) = 0    ! matrice assemblee
    id%ICNTL(18) = 0   ! matrice non distribuee (centralisee)

    !--------------------------------------------------------------------------------------
    !--- Output streams
    id%ICNTL(2)=-1
    id%ICNTL(3)=-1

    !--------------------------------------------------------------------------------------
    !--- Level of printing for errors
    id%ICNTL(4) = 1   ! =1 : seuls les messages d'erreur sont imprimes

    !--------------------------------------------------------------------------------------
    !--- Permutations / Scaling de la matrice
    id%ICNTL(6) = 7   ! 7=choix automatique fait par le package

    !--------------------------------------------------------------------------------------
    !--- Reordonnancement de la matrice
    id%ICNTL(7) = 7   ! 7=choix automatique fait par le package (tout sauf la 1 - Josk1)

    !--------------------------------------------------------------------------------------
    !--- Definition du probleme dans le processus 0 (phase peu couteuse)
    if (id%MYID == 0) then

        !call CPU_TIME(time1)

        id%N = ndlt
        !id%NZ = (matg%ia(matg%nl+1)-1)*2+matg%nl
        id%NZ = (matg%ia(matg%nl+1)-1)+matg%nl
        id%NRHS = 2
        id%LRHS = id%N
        allocate(id%IRN(id%NZ)) ; id%IRN = 0
        allocate(id%JCN(id%NZ)) ; id%JCN = 0
        allocate(id%A(id%NZ))   ; id%A = 0.D0
        allocate(id%RHS(id%N*2)); id%RHS = 0.D0

        !id%A(1:matg%nl)=matg%Da
        pos=1
        do i=1,matg%nl
            do k=matg%ia(i),matg%ia(i+1)-1

                id%IRN(pos)=i
                id%JCN(pos)=matg%ja(k)                  ! id%JCN = (/(matg%ja(i), i =1,id%NZ)/)
                id%A(pos)=matg%La(k)

                !id%IRN(pos)=matg%ja(k)
                !id%JCN(pos)=i                         ! id%JCN = (/(matg%ja(i), i =1,id%NZ)/)
                !id%A(pos)=matg%Ua(k)

                pos=pos+1

                !id%IRN(pos+1)=matg%ja(k)
                !id%JCN(pos+1)=i                         ! id%JCN = (/(matg%ja(i), i =1,id%NZ)/)
                !id%A(pos+1)=matg%Ua(k)
                !pos=pos+2
            enddo

            id%A(pos)=matg%Da(i)
            id%JCN(pos)=i
            id%IRN(pos)=i
            pos=pos+1         !id%A = (/(matg%a(i), i =1,id%NZ)/)
        enddo

        id%RHS(1:id%N)= (/(res1(i), i =1,id%N)/)
        id%RHS(id%N+1:id%N*2) = (/(res2(i), i =1,id%N)/)

    endif

    !--------------------------------------------------------------------------------------
    !--- Analyse, Factorisation, Resolution du systeme
    ! id%JOB = 1,2,3 correspond aux phases (Analyse, Facto, Resol)
    ! id%JOB = 6 regroupe ces trois phases directement
    id%JOB = 6
    call DMUMPS(id)

    !--------------------------------------------------------------------------------------
    !--- Deallocation utilisateur
    if(id%MYID == 0) then
        sol1 = (/(id%RHS(i), i =1, id%N)/)
        sol2 = (/(id%RHS(i), i =id%N+1, id%N*2)/)
        deallocate (id%IRN)
        deallocate (id%JCN)
        deallocate (id%A)
        deallocate (id%RHS)
    endif

    !--------------------------------------------------------------------------------------
    !--- Deallocation interne
    id%JOB = -2
    call DMUMPS(id)

end subroutine MUMPS_CSRC


!------------------------------------------------------------------------------------------------------
!> @author JL Tailhan
!> @author J. Foulliaron (version 1.0 - 2011)
!
!> @brief
!> Procedure de resolution de deux système à l'aide de la librairie MUMPS \n
!>  matg*sol1 = res1 et matg*sol2 = res2 \n\n
!
!> => Meme routine que Solve_mumps, mais adaptée au format CSRC
!> (mais censee elle, etre optimpisee et plus rapide.)\n
!> NB : Cette routine considere la matrice du systeme comme non symetrique
!
!> @details
!> #### DESCRIPTION:
!>  Entrees :\n
!>   - matg : matrice du systeme lineaire\n
!>   - res1 : second membre du premier systeme\n
!>   - res2 : second membre du second systeme\n\n
!
!>  Sorties :\n
!>   - sol1 : solution du premier systeme\n
!>   - sol2 : solution du secon systeme\n
!------------------------------------------------------------------------------------------------------
subroutine MUMPS_CSRC222(matg, res1,res2, sol1, sol2)

    !**************************************************************************************
    !* DEFINITIONS
    !**************************************************************************************

    !--------------------------------------------------------------------------------------
    !--- Modules generaux
    use variables, only : ndlt
    use sparse_matrix
    use formatCSRC

    implicit none

    !include 'mpif.h'
    include 'dmumps_struc.h'

    !--------------------------------------------------------------------------------------
    !--- Variables in/out
    type(matrice_CSRC), intent(in) :: matg
    real*8, dimension(:), intent(in) :: res1, res2
    real*8, dimension(:), intent(out) :: sol1, sol2

    !--------------------------------------------------------------------------------------
    !--- Variables locales
    real*8, dimension(:), allocatable :: A
    integer (kind=4), dimension(:), allocatable :: ia, ja
    type (DMUMPS_STRUC) :: id
    integer :: i, k
    integer(kind=4)::pos
    pos=1

    !**************************************************************************************
    !* CORPS DE PROCEDURE
    !**************************************************************************************

    !--------------------------------------------------------------------------------------
    !--- Definition du communicateur
    !id%COMM = mpi_comm_world
    id%COMM = 0

    !--------------------------------------------------------------------------------------
    !--- Symmetrie
    id%SYM = 0

    !--------------------------------------------------------------------------------------
    !--- Travail du processus maitre
    id%PAR = 1

    !--------------------------------------------------------------------------------------
    !--- Initialisation
    id%JOB = -1
    call DMUMPS(id)

    !--------------------------------------------------------------------------------------
    !--- Parametres de controle
    ! Type d'envoi de la matrice : assemblee + centralisee
    id%ICNTL(5) = 0    ! matrice assemblee
    id%ICNTL(18) = 0   ! matrice non distribuee (centralisee)

    !--------------------------------------------------------------------------------------
    !--- Output streams
    id%ICNTL(2)=-1
    id%ICNTL(3)=-1

    !--------------------------------------------------------------------------------------
    !--- NIveau d'impression des erreurs
    id%ICNTL(4) = 1   ! =1 : seuls les messages d'erreur sont imprimes

    !--------------------------------------------------------------------------------------
    !--- Permutations / Scaling de la matrice
    id%ICNTL(6) = 7   ! 7=choix automatique fait par le package

    !--------------------------------------------------------------------------------------
    !--- Reordonnancement de la matrice
    id%ICNTL(7) = 7   ! 7=choix automatique fait par le package (tout sauf la 1 - Josk1)

    !--------------------------------------------------------------------------------------
    !--- Definition du probleme dans le processus 0
    if(id%MYID == 0) then

        id%N = ndlt
        id%NZ = (matg%ia(matg%nl+1)-1)*2+matg%nl
        id%NRHS = 2
        id%LRHS = id%N

        allocate(id%RHS(id%N*2)) ; id%RHS =0.d0
        allocate(A(id%NZ))       ; A = 0.d0
        allocate(ia(id%NZ))      ; ia = 0
        allocate(ja(id%NZ))      ; ja = 0
        pos=1

        do i=1,matg%nl
            do k=matg%ia(i),matg%ia(i+1)-1

                if (matg%La(k)/=0) then
                    ia(pos)=i
                    ja(pos)=matg%ja(k)                ! id%JCN = (/(matg%ja(i), i =1,id%NZ)/)
                    A(pos)=matg%La(k)
                    pos=pos+1
                endif

                if(matg%Ua(k)/=0)then
                    ia(pos)=matg%ja(k)
                    ja(pos)=i                         ! id%JCN = (/(matg%ja(i), i =1,id%NZ)/)
                    A(pos)=matg%Ua(k)
                    pos=pos+1
                endif
            enddo

            A(pos)=matg%Da(i)
            ja(pos)=i
            ia(pos)=i
            pos=pos+1         !id%A = (/(matg%a(i), i =1,id%NZ)/)
        enddo

        id%NZ=pos-1
        allocate(id%IRN(id%NZ)) ; id%IRN = 0
        allocate(id%JCN(id%NZ)) ; id%JCN = 0
        allocate(id%A(id%NZ))   ; id%A = 0.d0

        id%IRN=ia(1:pos-1)
        id%JCN=ja(1:pos-1)
        id%A=A(1:pos-1)
        deallocate(A,ia,ja)

        id%RHS(1:id%N)= (/(res1(i), i =1,id%N)/)
        id%RHS(id%N+1:id%N*2) = (/(res2(i), i =1,id%N)/)

    endif

    !--------------------------------------------------------------------------------------
    !--- Analyse, Factorisation, Resolution du systeme
    ! id%JOB = 1,2,3 correspond aux phases (Analyse, Facto, Resol)
    ! id%JOB = 6 regroupe ces trois phases directement
    id%JOB = 6

    call DMUMPS(id)

    !--------------------------------------------------------------------------------------
    !--- Deallocation utilisateur
    if(id%MYID == 0) then
        sol1 = (/(id%RHS(i), i =1, id%N)/)
        sol2 = (/(id%RHS(i), i =id%N+1, id%N*2)/)
        deallocate (id%IRN)
        deallocate (id%JCN)
        deallocate (id%A)
        deallocate (id%RHS)
    endif

    !--------------------------------------------------------------------------------------
    !--- Deallocation interne
    id%JOB = -2
    call DMUMPS(id)

end subroutine MUMPS_CSRC222

!------------------------------------------------------------------------------------------------------

end module mumps

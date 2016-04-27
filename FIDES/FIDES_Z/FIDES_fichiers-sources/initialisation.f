!----------------------------------------------------------------------------------------------------------
! MODULE: initialisation
!
!> @author JL Tailhan
!> @author J. Goncalvez (version 1.0 - 2010)
!
!> @brief
!> Routines d'initialisation.
!----------------------------------------------------------------------------------------------------------
module initialisation

    implicit none

    interface init_vec
        module procedure initint_vec, initreal_vec, initreal8_vec, initchaine_vec, initlogi_vec
    end interface

    interface init_mat
        module procedure initint_mat, initreal_mat, initreal8_mat, initlogi_mat
    end interface

contains


!------------------------------------------------------------------------------------------------------
!> @author JL Tailhan
!> @author J. Goncalvez (version 1.0 - 2010)
!
!> @brief
!> Allocation et initialisation (mise a zero) des differents vecteurs de champs (par element) 
!
!> @details
!> #### DESCRIPTION:
!>  A faire
!------------------------------------------------------------------------------------------------------
subroutine initva (va, kva, nomcomp)

    !**************************************************************************************
    !* DEFINITIONS
    !**************************************************************************************

    !--------------------------------------------------------------------------------------
    !--- Modules generaux
    use variables, only : varint, nelt, infele, ktypel, calco

    implicit none

    !--------------------------------------------------------------------------------------
    !--- Variables in/out
    real*8, dimension(:), allocatable, intent(out) :: va
    integer, dimension(:), allocatable, intent(out) :: kva
    character(len=*), intent(in) :: nomcomp

    !--------------------------------------------------------------------------------------
    !--- Variables locales
    integer :: idime, ie, nlieu, ncmp

    !**************************************************************************************
    !* CORPS DE PROCEDURE
    !**************************************************************************************

    idime = 0
    allocate(kva(nelt))
    kva = 0

    do ie = 1, nelt
        nlieu = size((infele(ktypel(ie))%Q),1)
        if(calco == 'NOEUD') nlieu = infele(ktypel(ie))%nnel
        ncmp = infele(ktypel(ie))%ncgr
        if(nomcomp == 'NOLI') ncmp = ncmp + 1
        if(nomcomp == 'VINT') ncmp = varint
        kva(ie) = idime + 1
        idime = idime + ncmp*nlieu
    enddo

    allocate(va(idime))
    va=0.d0

end subroutine initva


!------------------------------------------------------------------------------------------------------
!> @author JL Tailhan
!> @author J. Goncalvez (version 1.0 - 2010)
!
!> @brief
!> Procedure d'allocation et d'initialisation de la structure de donnees relative aux elements
!
!> @details
!> #### DESCRIPTION:
!>  Entree :\n
!>  - n : nombre de type d'element\n
!------------------------------------------------------------------------------------------------------
subroutine init_infele(n)

    !**************************************************************************************
    !* DEFINITIONS
    !**************************************************************************************

    !--------------------------------------------------------------------------------------
    !--- Modules generaux
    use variables

    implicit none

    !--------------------------------------------------------------------------------------
    !--- Variables in/out
    integer, intent(in) :: n

    !--------------------------------------------------------------------------------------
    !--- Variables locales
    integer :: i

    !**************************************************************************************
    !* CORPS DE PROCEDURE
    !**************************************************************************************

    allocate(infele(n))
    do i = 1, n
        infele(i)%nnel = 0
        infele(i)%ndln = 0
        infele(i)%order= 0
        infele(i)%ncgr = 0
    enddo
end subroutine init_infele


!----------------------------------------------------------------------------------------------------------
!--- INITIALISATION DES VECTEURS
!----------------------------------------------------------------------------------------------------------


!------------------------------------------------------------------------------------------------------
!> @author JL Tailhan
!> @author J. Goncalvez (version 1.0 - 2010)
!
!> @brief
!> Procedure d'initialisation de vecteurs d'entiers
!
!> @details
!> #### DESCRIPTION:
!>  A faire
!------------------------------------------------------------------------------------------------------
subroutine initint_vec(vec, taille , valeur)

    !**************************************************************************************
    !* DEFINITIONS
    !**************************************************************************************

    implicit none

    !--------------------------------------------------------------------------------------
    !--- Variables in/out
    integer, dimension(:), allocatable, intent(inout) :: vec
    integer, intent(in) :: taille
    integer, intent(in), optional :: valeur

    !**************************************************************************************
    !* CORPS DE PROCEDURE
    !**************************************************************************************

    allocate(vec(taille))
    vec = 0
    if (present(valeur)) vec = valeur

end subroutine initint_vec


!------------------------------------------------------------------------------------------------------
!> @author JL Tailhan
!> @author J. Goncalvez (version 1.0 - 2010)
!
!> @brief
!> Procedure d'initialisation de vecteurs de reels
!
!> @details
!> #### DESCRIPTION:
!>  A faire
!------------------------------------------------------------------------------------------------------
subroutine initreal_vec(vec, taille, valeur)

    !**************************************************************************************
    !* DEFINITIONS
    !**************************************************************************************

    implicit none

    !--------------------------------------------------------------------------------------
    !--- Variables in/out
    real, dimension(:), allocatable, intent(inout) :: vec
    integer, intent(in) :: taille
    real, intent(in), optional :: valeur

    !**************************************************************************************
    !* CORPS DE PROCEDURE
    !**************************************************************************************

    allocate(vec(taille))
    vec = 0.
    if (present(valeur)) vec = 0.

end subroutine initreal_vec


!------------------------------------------------------------------------------------------------------
!> @author JL Tailhan
!> @author J. Goncalvez (version 1.0 - 2010)
!
!> @brief
!> Procedure d'initialisation de vecteurs de reels (double precision)
!
!> @details
!> #### DESCRIPTION:
!>  A faire
!------------------------------------------------------------------------------------------------------
subroutine initreal8_vec(vec, taille, valeur)

    !**************************************************************************************
    !* DEFINITIONS
    !**************************************************************************************

    implicit none

    !--------------------------------------------------------------------------------------
    !--- Variables in/out
    real*8, dimension(:), allocatable, intent(inout) :: vec
    integer, intent(in) :: taille
    real*8, intent(in), optional :: valeur

    !**************************************************************************************
    !* CORPS DE PROCEDURE
    !**************************************************************************************

    allocate(vec(taille))
    vec = 0.d0
    if (present(valeur)) vec = valeur

end subroutine initreal8_vec


!------------------------------------------------------------------------------------------------------
!> @author JL Tailhan
!> @author J. Goncalvez (version 1.0 - 2010)
!
!> @brief
!> Procedure d'initialisation de vecteurs de chaines de caracteres
!
!> @details
!> #### DESCRIPTION:
!>  A faire
!------------------------------------------------------------------------------------------------------
subroutine initchaine_vec(vec, taille, L, valeur)

    !**************************************************************************************
    !* DEFINITIONS
    !**************************************************************************************

    implicit none

    !--------------------------------------------------------------------------------------
    !--- Variables in/out
    integer, intent(in) :: L
    character(len=L), dimension(:), allocatable, intent(inout) :: vec
    integer, intent(in) :: taille
    character(len=L), intent(in), optional :: valeur

    !**************************************************************************************
    !* CORPS DE PROCEDURE
    !**************************************************************************************

    allocate(vec(taille))
    vec = ''
    if (present(valeur)) vec = valeur

end subroutine initchaine_vec


!------------------------------------------------------------------------------------------------------
!> @author JL Tailhan
!> @author J. Goncalvez (version 1.0 - 2010)
!
!> @brief
!> Procedure d'initialisation de vecteurs de type booleen
!
!> @details
!> #### DESCRIPTION:
!>  A faire
!------------------------------------------------------------------------------------------------------
subroutine initlogi_vec(vec, taille, valeur)

    !**************************************************************************************
    !* DEFINITIONS
    !**************************************************************************************

    implicit none

    !--------------------------------------------------------------------------------------
    !--- Variables in/out
    logical, dimension(:), allocatable, intent(inout) :: vec
    integer, intent(in) :: taille
    logical, intent(in), optional :: valeur

    !**************************************************************************************
    !* CORPS DE PROCEDURE
    !**************************************************************************************

    allocate(vec(taille))
    vec = .false.
    if (present(valeur)) vec = valeur

end subroutine initlogi_vec

!----------------------------------------------------------------------------------------------------------
!--- INITIALISATION DES MATRICES
!----------------------------------------------------------------------------------------------------------


!------------------------------------------------------------------------------------------------------
!> @author JL Tailhan
!> @author J. Goncalvez (version 1.0 - 2010)
!
!> @brief
!> Procedure d'initialisation de matrice d'entiers
!
!> @details
!> #### DESCRIPTION:
!>  A faire
!------------------------------------------------------------------------------------------------------
subroutine initint_mat(mat, nbl, nbc, valeur)

    !**************************************************************************************
    !* DEFINITIONS
    !**************************************************************************************

    implicit none

    !--------------------------------------------------------------------------------------
    !--- Variables in/out
    integer, dimension(:,:), allocatable, intent(inout) :: mat
    integer, intent(in) :: nbl, nbc
    integer, intent(in), optional :: valeur

    !**************************************************************************************
    !* CORPS DE PROCEDURE
    !**************************************************************************************

    allocate(mat(nbl,nbc))
    mat = 0
    if (present(valeur)) mat = valeur

end subroutine initint_mat


!------------------------------------------------------------------------------------------------------
!> @author JL Tailhan
!> @author J. Goncalvez (version 1.0 - 2010)
!
!> @brief
!> Procedure d'initialisation de matrice de reels
!
!> @details
!> #### DESCRIPTION:
!>  A faire
!------------------------------------------------------------------------------------------------------
subroutine initreal_mat(mat, nbl, nbc, valeur)

    !**************************************************************************************
    !* DEFINITIONS
    !**************************************************************************************

    implicit none

    !--------------------------------------------------------------------------------------
    !--- Variables in/out
    real, dimension(:,:), allocatable, intent(inout) :: mat
    integer, intent(in) :: nbl, nbc
    real, intent(in), optional :: valeur

    !**************************************************************************************
    !* CORPS DE PROCEDURE
    !**************************************************************************************

    allocate(mat(nbl,nbc))
    mat = 0.
    if (present(valeur)) mat = valeur

end subroutine initreal_mat


!------------------------------------------------------------------------------------------------------
!> @author JL Tailhan
!> @author J. Goncalvez (version 1.0 - 2010)
!
!> @brief
!> Procedure d'initialisation de matrice de reels (double precision)
!
!> @details
!> #### DESCRIPTION:
!>  A faire
!------------------------------------------------------------------------------------------------------
subroutine initreal8_mat(mat, nbl, nbc, valeur)

    !**************************************************************************************
    !* DEFINITIONS
    !**************************************************************************************

    implicit none

    !--------------------------------------------------------------------------------------
    !--- Variables in/out
    real*8, dimension(:,:), allocatable, intent(inout) :: mat
    integer, intent(in) :: nbl, nbc
    real*8, intent(in), optional :: valeur

    !**************************************************************************************
    !* CORPS DE PROCEDURE
    !**************************************************************************************

    allocate(mat(nbl,nbc))
    mat = 0.D0
    if (present(valeur)) mat = valeur

end subroutine initreal8_mat


!------------------------------------------------------------------------------------------------------
!> @author JL Tailhan
!> @author J. Goncalvez (version 1.0 - 2010)
!
!> @brief
!> Procedure d'initialisation de matrice de booleens
!
!> @details
!> #### DESCRIPTION:
!>  A faire
!------------------------------------------------------------------------------------------------------
subroutine initlogi_mat(mat, nbl, nbc, valeur)

    !**************************************************************************************
    !* DEFINITIONS
    !**************************************************************************************

    implicit none

    !--------------------------------------------------------------------------------------
    !--- Variables in/out
    logical, dimension(:,:), allocatable, intent(inout) :: mat
    integer, intent(in) :: nbl, nbc
    logical, intent(in), optional :: valeur

    !**************************************************************************************
    !* CORPS DE PROCEDURE
    !**************************************************************************************

    allocate(mat(nbl,nbc))
    mat = .false.
    if (present(valeur)) mat = valeur

end subroutine initlogi_mat

!------------------------------------------------------------------------------------------------------

end module initialisation

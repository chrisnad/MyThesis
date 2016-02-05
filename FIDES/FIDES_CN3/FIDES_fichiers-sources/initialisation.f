!------------------------------------------------------------------------------
! MODULE: initialisation
!
!> @author JL Tailhan
!
!> @brief
!> Routines d'initialisation.
!>
!------------------------------------------------------------------------------
module initialisation

    implicit none

    interface init_vec
       module procedure initint_vec, initreal_vec, initreal8_vec, initchaine_vec, initlogi_vec
    end interface

    interface init_mat
       module procedure initint_mat, initreal_mat, initreal8_mat, initlogi_mat
    end interface

contains

!------------------------------------------------------------------------!
!       Allocation et initialisation (mise a zero) des differents        !
!                 vecteurs de champs (par element)                       !
!------------------------------------------------------------------------!
   subroutine initva (va, kva, nomcomp)

      use variables, only : varint, nelt, infele, ktypel, calco

      implicit none

      real*8, dimension(:), allocatable, intent(out) :: va
      integer, dimension(:), allocatable, intent(out) :: kva
      character(len=*), intent(in) :: nomcomp
      integer :: idime, ie, nlieu, ncmp

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
      end do

      allocate(va(idime))
      va=0.

   end subroutine initva

!--------------------------------------------------------------------------------!
!   Procedure d'allocation et d'initialisation de la structure de donnees        !
!                          relative aux elements                                 !
!  Entree :                                                                      !
!  - n : nombre de type d'element                                                !
!--------------------------------------------------------------------------------!

   subroutine init_infele(n)

       use variables

       implicit none

       integer, intent(in) :: n
       integer :: i

       allocate(infele(n))
       do i = 1, n
           infele(i)%nnel = 0
           infele(i)%ndln = 0
           infele(i)%order= 0
           infele(i)%ncgr = 0
       end do
   end subroutine init_infele


!*********************************************************************************
!                                Initialisation des vecteurs                     *
!*********************************************************************************

!--------------------------------------------------------------------------------!
!                                   Vecteurs d'entiers                           !
!--------------------------------------------------------------------------------!

   subroutine initint_vec(vec, taille , valeur)

      implicit none

      integer, dimension(:), allocatable, intent(inout) :: vec
      integer, intent(in) :: taille
      integer, intent(in), optional :: valeur

      allocate(vec(taille))
      vec = 0
      if (present(valeur)) vec = valeur

   end subroutine initint_vec


!---------------------------------------------------------------------------------!
!                             Vecteurs de reels                                   !
!---------------------------------------------------------------------------------!

   subroutine initreal_vec(vec, taille, valeur)

      implicit none

      real, dimension(:), allocatable, intent(inout) :: vec
      integer, intent(in) :: taille
      real, intent(in), optional :: valeur

      allocate(vec(taille))
      vec = 0.
      if (present(valeur)) vec = 0.

   end subroutine initreal_vec

!---------------------------------------------------------------------------------!
!                     Vecteurs de reels (double precision)                        !
!---------------------------------------------------------------------------------!

   subroutine initreal8_vec(vec, taille, valeur)

      implicit none

      real*8, dimension(:), allocatable, intent(inout) :: vec
      integer, intent(in) :: taille
      real*8, intent(in), optional :: valeur

      allocate(vec(taille))
      vec = 0.d0
      if (present(valeur)) vec = valeur

   end subroutine initreal8_vec

!---------------------------------------------------------------------------------!
!                           Vecteurs de chaines de caracteres                     !
!---------------------------------------------------------------------------------!

   subroutine initchaine_vec(vec, taille, L, valeur)

      implicit none

      integer :: L
      character(len=L), dimension(:), allocatable, intent(inout) :: vec
      integer, intent(in) :: taille
      character(len=L), intent(in), optional :: valeur

      allocate(vec(taille))
      vec = ''
      if (present(valeur)) vec = valeur

   end subroutine initchaine_vec

!---------------------------------------------------------------------------------!
!                                  Vecteurs de type logique                       !
!---------------------------------------------------------------------------------!

   subroutine initlogi_vec(vec, taille, valeur)

      implicit none

      logical, dimension(:), allocatable, intent(inout) :: vec
      integer, intent(in) :: taille
      logical, intent(in), optional :: valeur

      allocate(vec(taille))
      vec = .false.
      if (present(valeur)) vec = valeur

   end subroutine initlogi_vec



!**********************************************************************************
!                                Initialisation des matrices                      *
!**********************************************************************************

!---------------------------------------------------------------------------------!
!                                   Matrices d'entiers                            !
!---------------------------------------------------------------------------------!

   subroutine initint_mat(mat, nbl, nbc, valeur)

      implicit none

      integer, dimension(:,:), allocatable, intent(inout) :: mat
      integer, intent(in) :: nbl, nbc
      integer, intent(in), optional :: valeur

      allocate(mat(nbl,nbc))
      mat = 0
      if (present(valeur)) mat = valeur

   end subroutine initint_mat

!---------------------------------------------------------------------------------!
!                                  Matrices de reels                              !
!---------------------------------------------------------------------------------!

   subroutine initreal_mat(mat, nbl, nbc, valeur)

      implicit none

      real, dimension(:,:), allocatable, intent(inout) :: mat
      integer, intent(in) :: nbl, nbc
      real, intent(in), optional :: valeur

      allocate(mat(nbl,nbc))
      mat = 0.
      if (present(valeur)) mat = valeur

   end subroutine initreal_mat

!---------------------------------------------------------------------------------!
!                      Matrices de reels (double precision)                       !
!---------------------------------------------------------------------------------!

   subroutine initreal8_mat(mat, nbl, nbc, valeur)

      implicit none

      real*8, dimension(:,:), allocatable, intent(inout) :: mat
      integer, intent(in) :: nbl, nbc
      real*8, intent(in), optional :: valeur

      allocate(mat(nbl,nbc))
      mat = 0.D0
      if (present(valeur)) mat = valeur

   end subroutine initreal8_mat

!---------------------------------------------------------------------------------!
!                                  Matrices de type logique                       !
!---------------------------------------------------------------------------------!

   subroutine initlogi_mat(mat, nbl, nbc, valeur)

      implicit none

      logical, dimension(:,:), allocatable, intent(inout) :: mat
      integer, intent(in) :: nbl, nbc
      logical, intent(in), optional :: valeur

      allocate(mat(nbl,nbc))
      mat = .false.
      if (present(valeur)) mat = valeur

   end subroutine initlogi_mat


end module initialisation

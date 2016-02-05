module Assemblage_CSRC

!------------------------------------------------------------------------------------------------!
! Module effectuant l'assemblage lors de l'utilisation du format CSRC
!------------------------------------------------------------------------------------------------!

contains

 subroutine Assem_CSRC(matg,vres)

    use variables, only : nelt, vprelg, idmax, kprop, vsol, explicite, kloce
    use formatCSRC
    use lib_elem, only : elem_kloce2
    use calc_elem, only : elem_ke
    use math

    implicit none

    !include 'mpif.h'

    type(matrice_CSRC), intent(inout) :: matg   ! Matrice globale
    real*8, dimension(:), intent(inout) :: vres  ! Vecteur residu global

    real*8, dimension(:,:), allocatable :: vme, vke  ! Matrice elementaire
    real*8, dimension(:), allocatable :: vre  ! Vecteur residu elementaire
    real*8, dimension(:), allocatable :: vdle
    real*8, dimension(idmax) :: vprel ! proprietes elementaires

    !real*8 :: T1, T2, tdeb, tfin ! temps cpu pour mesurer dans le detail temps phase d'assemblage
    integer :: ie, ndle

    ! Mise a zeros des valeurs de matg
    matg%La = 0.0d0
    matg%Ua = 0.0d0
    matg%Da = 0.0d0

    !--------------------- Boucle sur les elements ---------------------
    do ie = 1, nelt

       !--------------------- Deplacement total ------------------------
       allocate(vdle(ndle))
       call elem_kloce2(ie,ndle)
       vdle = vsol(kloce(1:ndle))

       !-- Calcul de vke et vre : matrice tangente et residu elementaires --
       vprel = vprelg(kprop(ie),1:idmax)  ! Proprietes elementaires

       if (explicite) then
           !------- Assemblage de vme --------
           ! call elem_me(vme,vre,ie,kloce,vprel,vdle)
           ! todo ...
       else
           !------- Assemblage de vke ---------
           allocate(vke(ndle,ndle)) ;   allocate(vre(ndle)) !-- vke et vre initialises dans elem_ke !
           call elem_ke(vke,vre,ie,ndle,vprel,vdle)

           !============================================================================!

           call Assem_rapide22(matg,vke,kloce(1:ndle))


       end if !fin if explicite

       !------------ Assemblage de vre ----------
       vres(kloce(1:ndle)) = vres(kloce(1:ndle)) + vre

       !------------ Desallocation --------------
       deallocate(vre,vdle)

       if (explicite) then
           deallocate(vme)
       else
           deallocate(vke)
       end if

    end do

  end subroutine Assem_CSRC

end module Assemblage_CSRC

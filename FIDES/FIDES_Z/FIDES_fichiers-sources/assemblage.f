!----------------------------------------------------------------------------------------------------------
! MODULE: assemblage
!
!> @author JL Tailhan
!> @author J. Goncalvez (version 1.0 - 2010)
!
!> @brief
!> Routines d'assemblage de la matrice de rigidite globale et du vecteur
!> residu global
!----------------------------------------------------------------------------------------------------------
module assemblage


contains


!------------------------------------------------------------------------------------------------------
!> @authors J. Goncalvez (version 1.0 - 2010)
!
!> @brief
!> Calcul et assemblage de la rigidite globale Kg et du vecteur residu global vres
!
!> @details
!> #### DESCRIPTION:
!>  Entree :\n
!>  - ipas : pas de calcul\n\n
!>
!>  Sorties :\n
!>  - matg : matrice globale\n
!>  - vres : vecteur residu global\n
!------------------------------------------------------------------------------------------------------
subroutine assem(matg,vres)  !,ipas)

    !**************************************************************************************
    !* DEFINITIONS
    !**************************************************************************************

    !--------------------------------------------------------------------------------------
    !--- Modules generaux
    use variables, only : nelt, vprelg, idmax, kprop, vsol, explicite, kloce
    use lib_elem, only : elem_kloce2
    use calc_elem, only : elem_ke
    use sparse_matrix

    implicit none

    !--------------------------------------------------------------------------------------
    !--- Variables in/out
    type(sparse), intent(inout) :: matg   ! Matrice globale
    real*8, dimension(:), intent(inout) :: vres  ! Vecteur residu global
    !integer, intent(in) :: ipas

    !--------------------------------------------------------------------------------------
    !--- Variables locales
    real*8, dimension(:,:), allocatable :: vme, vke  ! Matrice elementaire
    real*8, dimension(:), allocatable :: vre  ! Vecteur residu elementaire
    real*8, dimension(:), allocatable :: vdle
    real*8, dimension(idmax) :: vprel ! proprietes elementaires
    integer :: ie, ndle

    ! Mise a zeros des valeurs de matg
    matg%a = 0.

    !**************************************************************************************
    !* CORPS DE PROCEDURE
    !**************************************************************************************

    !--------------------------------------------------------------------------------------
    !--- Boucle sur les elements
    do ie = 1, nelt
        call elem_kloce2(ie, ndle)

        !----------------------------------------------------------------------------------
        !--- Deplacement total
        allocate(vdle(ndle))  ;
        vdle = vsol(kloce(1:ndle))

        !----------------------------------------------------------------------------------
        !--- Calcul de vke et vre : matrice tangente et residu elementaires
        vprel = vprelg(kprop(ie),1:idmax)  ! Proprietes elementaires

        !----------------------------------------------------------------------------------
        !--- Assemblage de vme
        if(explicite) then
            ! call elem_me(vme,vre,ie,ndle,vprel,vdle)
            ! todo ...

        !----------------------------------------------------------------------------------
        !--- Assemblage de vke
        else
            allocate(vke(ndle,ndle)) ;   allocate(vre(ndle)) ! vke et vre initialises dans elem_ke
            call elem_ke(vke,vre,ie,ndle,vprel,vdle)
            call assem_sparse(matg, vke, kloce(1:ndle), kloce(1:ndle))
        endif

        !----------------------------------------------------------------------------------
        !--- Assemblage de vre
        vres(kloce(1:ndle)) = vres(kloce(1:ndle)) + vre

        !----------------------------------------------------------------------------------
        !--- Desallocation
        deallocate(vdle,vre)
        if(explicite) then
            deallocate(vme)
        else
            deallocate(vke)
        endif

    enddo

    !--------------------------------------------------------------------------------------
    !--- Economie de la place memoire
    !if(ipas == 1) call Sparse_reduc_taille(matg)
    !  On elimine les valeurs nulles de la matrice vkg
    call del_zeros(matg)

end subroutine assem

!------------------------------------------------------------------------------------------------------

end module assemblage

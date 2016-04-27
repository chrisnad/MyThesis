!----------------------------------------------------------------------------------------------------------
! MODULE: aiguillage
!
!> @author JL Tailhan
!> @author J. Goncalvez (version 1.0 - 2010)
!
!> @brief
!> Routines d'aiguillage entre modèles de fissuration (explicite et semi-explicite)
!----------------------------------------------------------------------------------------------------------
module aiguillage


contains

!------------------------------------------------------------------------------------------------------
!> @authors J. Goncalvez (version 1.0 - 2010)
!
!> @brief
!> Fonction d'aiguillage pour le pilotage du calcul
!
!> @details
!> #### DESCRIPTION:
!> A faire
!------------------------------------------------------------------------------------------------------
subroutine pilot(imetpilo,vsol,vduI,vduII,dlam)

    !**************************************************************************************
    !* DEFINITIONS
    !**************************************************************************************

    !--------------------------------------------------------------------------------------
    !--- Modules generaux
    use variables, only :   interf, fiss, iedngi, iedngm, iedng, &
    &                       nelcritic, mrtrav, dime, nelt, inorm, irloc, &
    !    Ajout pour element macro BA
    &                       fissBA, iedngba, &
    !    Ajout pour pilotage indirect (npilot=3)
    & npilot, &
    !    Ajout pour algo pb evolutif (ievtp)
    & alp_Dt

    use element_interface, only : interf_pilot
    use fissuration,    only : fiss_pilot
    use fissuration_BA, only : fissBA_pilot
    use initialisation, only : init_mat

    implicit none

    !--------------------------------------------------------------------------------------
    !--- Variables in/out
    real*8, dimension(:), intent(in)    :: vduII, vsol
    real*8, dimension(:), intent(inout) :: vduI
    real*8, intent(inout)   :: dlam
    integer, intent(in)     :: imetpilo

    !--------------------------------------------------------------------------------------
    !--- Variables locales
    real*8 :: alpha, alphai, alpham
    integer :: elefissi, elefissm
    character(len=8) :: MOTi, MOTm, MOT
    !    Ajout pour element macro BA
    real*8 :: alphaba
    integer :: elefissba
    character(len=8) :: MOTba

    !**************************************************************************************
    !* CORPS DE PROCEDURE
    !**************************************************************************************

    alpha   = 1.d0
    alphai  = 1.d0
    alpham  = 1.d0
    alphaba = 1.d0

    MOTi  = '    '
    MOTm  = '    '
    MOT   = '    '
    iedng = 0

    elefissi = 0
    elefissm = 0
    elefissba = 0

    call init_mat(mrtrav,nelt,dime*dime)

    if (interf == 1) call interf_pilot(imetpilo,vsol,vduI,vduII,dlam,alphai,MOTi,elefissi)
    if (fiss == 1)   call fiss_pilot(imetpilo,vsol,vduI,vduII,dlam,alpham,MOTm,elefissm)
    if (fissBA == 1) call fissBA_pilot(imetpilo,vsol,vduI,vduII,dlam,alphaba,MOTba,elefissba)

    !--------------------------------------------------------------------------------------
    !--- Choix du minimum de pas de chargement entre les trois pilotages
    alpha = minval((/alphai,alpham,alphaba/))

    if (alpha /= 1.d0) then

        !----------------------------------------------------------------------------------
        !--- Element d'interface
        if (alpha == alphai) then
            iedng = iedngi
            MOT = MOTi
        endif

        !----------------------------------------------------------------------------------
        !--- Element fissurant
        if (alpha == alpham) then
            iedng = iedngm
            MOT = MOTm
            !----- Et on stocke son repere principal
            !irloc(iedngm,:) = mrtrav(iedngm,:)
            !inorm(iedngm,:) = mrtrav(iedngm,1:dime)
            !inorm(iedngm,:) = normfissm
        endif
	
        !----------------------------------------------------------------------------------
        !--- Element BA fissurant
        if (alpha == alphaba) then
            iedng = iedngba
            MOT = MOTba
            !----- Et on stocke son repere principal
            !irloc(iedngba,:) = mrtrav(iedngba,:)
            !inorm(iedngba,:) = mrtrav(iedngba,1:dime)
            !inorm(iedngba,:) = normfissba
        end if
        !----------------------------------------------------------------------------------
        !--- Pour tous les modes de pilotage sauf pilotage indirect (npilot=3)
        if (npilot/=3) then
            if (alpha <= 0.d0) alpha = 1.d-10
        endif
        if (alpha >= 1.d0) alpha = 1.d0

        print'(a11,e11.6,a11,e11.6,a13,i7,a3,a8,a21,i7,2x,a4)','Pilotage : ', &
        &      alpha,' dlam : ',dlam,'   element : ',iedng,'  ',MOT,' [ elts critiques : ',elefissi+elefissm,' ]'

        nelcritic = elefissi+elefissm+elefissba

        !----------------------------------------------------------------------------------
        !--- Si nelcritic est diff de zeros alors on enleve 1 pour tenir compte de iedng dans le nombre
        if (nelcritic/=0) nelcritic = nelcritic
    else
        nelcritic = 0
    endif

    !--------------------------------------------------------------------------------------
    !--- Ajout pour algo pb evolutif (option: ievtp)
    alp_Dt = alpha


    !--------------------------------------------------------------------------------------
    !--- On reactualise le champ de deplacements (vduI et facteur de charge)
    vduI = alpha*vduI
    dlam = alpha*dlam

    !-----------------------------------------------------------------
    !---  Recalcul du repere de fissuration ou du repere principal (fiss=1)
    if (iedng /= 0) then
        if ((fiss == 1).and.(iedng==iedngm)) then
            call fiss_pilot(imetpilo,vsol,vduI,vduII,dlam,alpham,MOTm,elefissm)
            if (.not.allocated(irloc)) call init_mat(irloc,nelt,dime*dime)
            !----- Et on stocke son repere principal
            irloc(iedng,:) = mrtrav(iedng,:)
            inorm(iedng,:) = mrtrav(iedng,1:dime)
         endif
        if ((fissba == 1).and.(iedng==iedngba)) then
            call fissBA_pilot(imetpilo,vsol,vduI,vduII,dlam,alpham,MOTm,elefissm,iedng)
            if (.not.allocated(irloc)) call init_mat(irloc,nelt,dime*dime)
            !----- Et on stocke son repere principal
            irloc(iedng,:) = mrtrav(iedng,:)
            inorm(iedng,:) = mrtrav(iedng,1:dime)
         endif
    endif
    !-----------------------------------------------------------------

    deallocate(mrtrav)

end subroutine pilot


!------------------------------------------------------------------------------------------------------
!> @authors J. Goncalvez (version 1.0 - 2010)
!
!> @brief
!> Fonction d'aiguillage pour la gestion du changement d'etat
!
!> @details
!> #### DESCRIPTION:
!> A faire
!------------------------------------------------------------------------------------------------------
subroutine change_etat()

    !**************************************************************************************
    !* DEFINITIONS
    !**************************************************************************************

    !--------------------------------------------------------------------------------------
    !--- Modules generaux
    use variables, only : interf, fiss
    use element_interface, only : interf_change_etat
    use fissuration, only : fiss_change_etat

    implicit none

    !**************************************************************************************
    !* CORPS DE PROCEDURE
    !**************************************************************************************

    if (interf == 1) call interf_change_etat
    if (fiss == 1)   call fiss_change_etat

end subroutine change_etat


!------------------------------------------------------------------------------------------------------
!> @authors J. Goncalvez (version 1.0 - 2010)
!
!> @brief
!> Fonction d'aiguillage pour la rupture
!
!> @details
!> #### DESCRIPTION:
!> A faire
!------------------------------------------------------------------------------------------------------
subroutine rupture(ie,vsig,vnle,wplael,vprel)

    !**************************************************************************************
    !* DEFINITIONS
    !**************************************************************************************

    !--------------------------------------------------------------------------------------
    !--- Modules generaux
    use variables, only : interf, fiss
    use element_interface, only : interf_rupture
    use fissuration, only : fiss_rupture

    implicit none

    !--------------------------------------------------------------------------------------
    !--- Variables in/out
    real*8, dimension(:,:), intent(inout) :: vsig, vnle
    integer, intent(in) :: ie
    real*8, dimension(:), intent(in), optional :: vprel
    real*8, intent(in), optional :: wplael

    !**************************************************************************************
    !* CORPS DE PROCEDURE
    !**************************************************************************************

    if (interf == 1) call interf_rupture(ie,vsig,vnle)
    if (fiss == 1)   call fiss_rupture(ie,vsig,vnle,wplael,vprel)

end subroutine rupture


!------------------------------------------------------------------------------------------------------
!> @authors J. Goncalvez (version 1.0 - 2010)
!
!> @brief
!> Fonction d'aiguillage pour la gestion de la raideur
!
!> @details
!> #### DESCRIPTION:
!> A faire
!------------------------------------------------------------------------------------------------------
subroutine modul(ie,ipg,vh,vprel)

    !**************************************************************************************
    !* DEFINITIONS
    !**************************************************************************************

    !--------------------------------------------------------------------------------------
    !--- Modules generaux
    use variables, only : interf, fiss
    use element_interface, only : interf_modul
    use fissuration, only : fiss_modul

    implicit none

    !--------------------------------------------------------------------------------------
    !--- Variables in/out
    real*8, dimension(:), intent(in) :: vprel
    real*8, dimension(:,:), intent(inout) :: vh
    integer, intent(in) :: ie, ipg

    !**************************************************************************************
    !* CORPS DE PROCEDURE
    !**************************************************************************************

    if (interf == 1) call interf_modul(ie,ipg,vh,vprel)
    if (fiss == 1)   call fiss_modul(ie,vh,vprel)

end subroutine modul


!------------------------------------------------------------------------------------------------------
!> @authors J. Goncalvez (version 1.0 - 2010)
!
!> @brief
!> Fonction d'aiguillage pour la gestion du comportement
!
!> @details
!> #### DESCRIPTION:
!> A faire
!------------------------------------------------------------------------------------------------------
subroutine loi(iloi,icomp,ie)!,ipg)

    !**************************************************************************************
    !* DEFINITIONS
    !**************************************************************************************

    !--------------------------------------------------------------------------------------
    !--- Modules generaux
    use variables, only : interf, fiss, fissBA
    use element_interface, only :interf_loi
    use fissuration, only :fiss_loi
    use fissuration_BA, only :fissBA_loi

    !--------------------------------------------------------------------------------------
    !--- Variables in/out
    integer, intent(in) :: icomp, ie!, ipg
    integer, intent(out) :: iloi


    !**************************************************************************************
    !* CORPS DE PROCEDURE
    !**************************************************************************************

    iloi = icomp

    if (fiss == 1) call fiss_loi(iloi,icomp,ie)
    if (interf==1) call interf_loi(iloi,icomp,ie)!,ipg)
    if (fissBA==1) call fissBA_loi(iloi,icomp)

end subroutine loi


!------------------------------------------------------------------------------------------------------
!> @authors J. Goncalvez (version 1.0 - 2010)
!
!> @brief
!> Fonction pour le stockage des resultats a chaque pas de calcul
!
!> @details
!> #### DESCRIPTION:
!> A faire
!------------------------------------------------------------------------------------------------------
subroutine stock()

    !**************************************************************************************
    !* DEFINITIONS
    !**************************************************************************************

    !--------------------------------------------------------------------------------------
    !--- Modules generaux
    use variables, only : interf, fiss
    use element_interface, only : interf_stock
    use fissuration, only : fiss_stock

    implicit none

    !**************************************************************************************
    !* CORPS DE PROCEDURE
    !**************************************************************************************

    if (interf == 1) call interf_stock()
    if (fiss == 1)   call fiss_stock()

end subroutine stock


!------------------------------------------------------------------------------------------------------
!> @authors J. Goncalvez (version 1.0 - 2010)
!
!> @brief
!> Fonction d'aiguillage pour la definition des proprietes mecaniques aleatoires
!
!> @details
!> #### DESCRIPTION:
!> A faire
!------------------------------------------------------------------------------------------------------
subroutine distal()

    !**************************************************************************************
    !* DEFINITIONS
    !**************************************************************************************

    !--------------------------------------------------------------------------------------
    !--- Modules generaux
    use variables, only : interf, fiss, fissBA
    use element_interface, only : interf_distal
    use fissuration, only : fiss_distal
    use fissuration_BA, only : fissBA_distal

    implicit none

    !**************************************************************************************
    !* CORPS DE PROCEDURE
    !**************************************************************************************

    if (interf == 1) call interf_distal
    if (fiss == 1)   call fiss_distal
    if (fissBA == 1) call fissBA_distal

end subroutine distal


!------------------------------------------------------------------------------------------------------


end module aiguillage

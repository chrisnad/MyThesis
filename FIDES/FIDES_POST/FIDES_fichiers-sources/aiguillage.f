!------------------------------------------------------------------------------
! MODULE: aiguillage
!
!> @author JL Tailhan
!
!> @brief
!> Routines d'aiguillage entre modèles de fissuration (explicite et semi-explicite)
!>
!------------------------------------------------------------------------------
module aiguillage


contains

!*****************************************************************!

subroutine pilot(imetpilo,vsol,vduI,vduII,dlam)

!-----------------------------------------------------------------!
! Fonction d'aiguillage pour le pilotage du calcul                !
!-----------------------------------------------------------------!

    use variables, only : interf, fiss, iedngi, iedngm, iedng, &
                 & iecritic, mrtrav, dime, nelt, inorm, irloc, &
                 !--- Ajout pour algo pb evolutif (ievtp) ---!
                 & alp_Dt
                 
    use element_interface, only : interf_pilot
    use fissuration, only : fiss_pilot
    use initialisation, only : init_mat

    implicit none

    real*8, dimension(:), intent(in) :: vduII, vsol
    real*8, dimension(:), intent(inout) :: vduI
    real*8, intent(inout) :: dlam
    integer, intent(in) :: imetpilo
    real*8 :: alpha, alphai, alpham
        integer :: elefissi, elefissm
    character(len=8) :: MOTi, MOTm, MOT

     alpha   = 1.d0
     alphai  = 1.d0
     alpham  = 1.d0
     MOTi  = '    '
     MOTm  = '    '
     MOT   = '    '
     iedng = 0

     elefissi = 0
     elefissm = 0

     call init_mat(mrtrav,nelt,dime*dime)

     if (interf == 1) call interf_pilot(imetpilo,vsol,vduI,vduII,dlam,alphai,MOTi,elefissi)
     if (fiss == 1) call fiss_pilot(imetpilo,vsol,vduI,vduII,dlam,alpham,MOTm,elefissm)

     !----- Choix du minimum de pas de chargement entre les deux pilotages
     alpha = minval((/alphai,alpham/))
     
     if (alpha /= 1.d0) then

        ! Element d'interface
        if (alpha == alphai) then
            iedng = iedngi
            MOT = MOTi
        end if

        ! Element fissurant
        if (alpha == alpham) then
            iedng = iedngm
            MOT = MOTm
            !----- Et on stocke son repere principal
            irloc(iedngm,:) = mrtrav(iedngm,:)
            inorm(iedngm,:) = mrtrav(iedngm,1:dime)
        end if

        if (alpha <= 0.d0) alpha = 1.d-10
        if (alpha >= 1.d0) alpha = 1.d0

        print'(a11,e11.6,a11,e11.6,a13,i7,a3,a8,a21,i7,2x,a4)','Pilotage : ', &
             &      alpha,' dlam : ',dlam,'   element : ',iedng,'  ',MOT,' [ elts critiques : ',elefissi+elefissm,' ]'
        
        iecritic = elefissi+elefissm
    end if
    
    !----- Ajout pour algo pb evolutif (option: ievtp)
    alp_Dt = alpha


    !----- On reactualise le champ de deplacements (vduI et facteur de charge)
    vduI = alpha*vduI
    dlam = alpha*dlam

    deallocate(mrtrav)

end subroutine pilot

!-------------------------------------------------------!

subroutine change_etat()

!-----------------------------------------------------------------!
! Fonction d'aiguillage pour la gestion du changement d'etat      !
!-----------------------------------------------------------------!

     use variables, only : interf, fiss
     use element_interface, only : interf_change_etat
     use fissuration, only : fiss_change_etat

     implicit none

     if (interf == 1) call interf_change_etat
     if (fiss == 1)   call fiss_change_etat

end subroutine change_etat

!-------------------------------------------------------!

subroutine rupture(ie,vsig,vnle,wplael,vprel)

!-----------------------------------------------------------------!
! Fonction d'aiguillage pour la gestion du changement d'etat      !
!-----------------------------------------------------------------!

     use variables, only : interf, fiss
     use element_interface, only : interf_rupture
     use fissuration, only : fiss_rupture

     implicit none

     real*8, dimension(:,:), intent(inout) :: vsig, vnle
     integer, intent(in) :: ie
     real*8, dimension(:), intent(in), optional :: vprel
     real*8, intent(in), optional :: wplael

     if (interf == 1) call interf_rupture(ie,vsig,vnle)
     if (fiss == 1)   call fiss_rupture(ie,vsig,vnle,wplael,vprel)

end subroutine rupture

!-------------------------------------------------------!

subroutine modul(ie,ipg,vh,vprel)

!-----------------------------------------------------------------!
! Fonction d'aiguillage pour la gestion de la raideur             !
!-----------------------------------------------------------------!

     use variables, only : interf, fiss
     use element_interface, only : interf_modul
     use fissuration, only : fiss_modul

     implicit none

     real*8, dimension(:), intent(in) :: vprel
     real*8, dimension(:,:), intent(inout) :: vh
     integer, intent(in) :: ie, ipg

     if (interf == 1) call interf_modul(ie,ipg,vh,vprel)
     if (fiss == 1)   call fiss_modul(ie,vh,vprel)

end subroutine modul

!-------------------------------------------------------!


subroutine loi(iloi,icomp,ie)!,ipg)

!-----------------------------------------------------------------!
! Fonction d'aiguillage pour la gestion du comportement           !
!-----------------------------------------------------------------!

     use variables, only : interf, fiss, fissBA
     use element_interface, only :interf_loi
     use fissuration, only :fiss_loi
     use fissuration_BA, only :fissBA_loi

     integer, intent(in) :: icomp, ie!, ipg
     integer, intent(out) :: iloi
     
     if (icomp /= 25) then
        if (fiss == 1) call fiss_loi(iloi,icomp,ie)
        if (interf==1) call interf_loi(iloi,icomp,ie)!,ipg)
     else
        if (fissBA==1) call fissBA_loi(iloi,icomp,ie)

     endif

end subroutine loi

!-------------------------------------------------------!

subroutine stock()

!-----------------------------------------------------------------!
!            Fonction pour le stockage des resultats a chaque     !
!                              pas de calcul                      !
!-----------------------------------------------------------------!

     use variables, only : interf, fiss
     use element_interface, only : interf_stock
     use fissuration, only : fiss_stock

     implicit none

     if (interf == 1) call interf_stock()
     if (fiss == 1)   call fiss_stock()

end subroutine stock

!-------------------------------------------------------!

subroutine distal()

!-----------------------------------------------------------------!
! Fonction d'aiguillage pour la definition des proprietes         !
!                     mecaniques aleatoires                       !
!-----------------------------------------------------------------!

     use variables, only : interf, fiss, fissBA
     use element_interface, only : interf_distal
     use fissuration, only : fiss_distal
     use fissuration_BA, only : fissBA_distal

     implicit none

     if (interf == 1) call interf_distal
     if (fiss == 1)   call fiss_distal
     if (fissBA == 1) call fissBA_distal

end subroutine distal

!-------------------------------------------------------!

end module aiguillage

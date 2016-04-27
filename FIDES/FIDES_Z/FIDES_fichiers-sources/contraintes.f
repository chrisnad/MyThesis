!----------------------------------------------------------------------------------------------------------
! MODULE: contraintes
!
!> @author JL Tailhan
!> @author J. Goncalvez (version 1.0 - 2010)
!
!> @brief
!> Fonction d'aiguillage pour le calcul des contraintes.
!----------------------------------------------------------------------------------------------------------
module contraintes


contains


!------------------------------------------------------------------------------------------------------
!> @author JL Tailhan
!> @authors J. Goncalvez (version 1.0 - 2010)
!> @authors T.S. Song    (version 1.0 - 2012)
!> @authors C. Nader     (version 1.0 - 2014)
!
!> @brief
!> Fonction d'aiguillage pour le calcul des contraintes
!
!> @details
!> #### DESCRIPTION:
!>  A faire
!------------------------------------------------------------------------------------------------------
subroutine elem_sig(vsi,vnl,vin,ie,vhep,vb,vnl0,vprel,vdle,ipg)

    !**************************************************************************************
    !* DEFINITIONS
    !**************************************************************************************

    !--------------------------------------------------------------------------------------
    !--- Modules generaux
    use initialisation, only : init_mat, init_vec
    use aiguillage, only : loi, modul
    use element_interface, only : interf_sig, interf_sig_endo
    use fissuration, only : fiss_sig
    use fissuration, only : fiss_sig_endo
    use acier, only : acier_sig
    use fissuration_BA, only : fissBA_sig_endo_ortho
    use lib_elem, only : elem_hooke
    use variables, only : ktypel, nomtype, ietatpg

    implicit none

    !--------------------------------------------------------------------------------------
    !--- Variables in/out
    real*8, dimension(:,:), intent(in) :: vb
    integer, intent(in) :: ie, ipg
    real*8, dimension(:), intent(in) :: vnl0, vdle
    real*8, dimension(:), intent(in) :: vprel
    real*8, dimension(:), intent(inout) :: vsi
    real*8, dimension(:), intent(out) :: vnl, vin
    real*8, dimension(:,:), allocatable, intent(out) :: vhep

    !--------------------------------------------------------------------------------------
    !--- Variables locales
    real*8, dimension(:,:), allocatable :: vh
    integer :: icomp, iloi

    !--- Variables permettant l'utilisation de lois complexes (combinant des lois simples)
    !    loi 158 = 105 + 108 (108 ne pouvant fonctionner seule)
    real*8, dimension(size(vb,1),size(vb,1)) :: vhep2
    real*8, dimension(size(vsi,1)) :: vsi2
    real*8, dimension(size(vin,1)) :: vin2
    integer :: iloi1,iloi2

    !**************************************************************************************
    !* CORPS DE PROCEDURE
    !**************************************************************************************

    !--------------------------------------------------------------------------------------
    !--- On recupere le type de comportement
    icomp = int(vprel(1))   ! icomp pour chaque elements

    !--------------------------------------------------------------------------------------
    !--- On verifie la compatibilite : Comportement / Loi
    if (icomp /= 1) then
        call loi(iloi,icomp,ie)!,ipg)
    else
        iloi = 1
    endif

    !--------------------------------------------------------------------------------------
    !--- matrice d'elasticite vh, en fonction du type d'element
    call elem_hooke(vh,nomtype(ktypel(ie)),vprel)
    call modul(ie,ipg,vh,vprel)

    call init_mat(vhep,size(vh,1),size(vh,2))

    !--------------------------------------------------------------------------------------
    !--- Aiguillage en fonction de la loi de comportement
    select case (iloi)

        !----------------------------------------------------------------------------------
        !--- Cas de l'elasticite lineaire
        case(1, 2)
            vsi = matmul(vh,matmul(vb,vdle))
            vhep = vh

        !----------------------------------------------------------------------------------
        !--- Cas du modele de Von Mises pour les aciers
        case(11)
            !2016-03-18 JLT MODIF call acier_sig(vhep,vsi,vnl,vh,vb,vnl0,vprel,vdle,iloi,ie,ipg)
            call acier_sig(vhep,vsi,vnl,vh,vb,vnl0,vdle,ie)

        !----------------------------------------------------------------------------------
        !--- Cas des modèles non linéaires (type plasticité ou fissuration)
        case(12, 13, 14)
            call fiss_sig(vsi,vnl,vh,vb,vnl0,vprel,vdle,iloi,ie)
            vhep = vh

!           case(12)
!               call fiss_sig_plas(vsi,vnl,vh,vb,vnl0,vprel,vdle,iloi,ie)
!               vhep = vh

        !----------------------------------------------------------------------------------
        !--- Cas du modele macro de fissuration (endomagement)
        case(15)
            call fiss_sig_endo(vsi,vnl,vin,vh,vb,vprel,vdle,iloi,ie)
            vhep = vh

        !----------------------------------------------------------------------------------
        !--- Cas du modele macro de fissuration (endomagement) BA orthotrope
        case(25, 26, 27)
            call fissBA_sig_endo_ortho(vsi,vin,vh,vb,vprel,vdle,iloi,ie)
            vhep = vh

        !----------------------------------------------------------------------------------
        !--- Cas du modele de fissuration d'interface elastique fragile
        case(100, 101)
            call interf_sig(vhep,vsi,vnl,vh,vb,vnl0,vprel,vdle,iloi,ie,ipg)

        !----------------------------------------------------------------------------------
        !--- Cas des modeles simples d'interface endommageable puis frottante (pour beton, BA et BRF)
        case(102, 105, 106, 108)
            !if (all(ietatpg(ie,:)==0).or.all(ietatpg(ie,:)==3)) then
            if (any(ietatpg(ie,:)==0).or.any(ietatpg(ie,:)==3)) then
                call interf_sig_endo(vhep,vsi,vnl,vin,vh,vb,vnl0,vprel,vdle,iloi,ie,ipg)
            else
                call interf_sig(vhep,vsi,vnl,vh,vb,vnl0,vprel,vdle,iloi,ie,ipg)
            endif

        !----------------------------------------------------------------------------------
        !--- Cas du modele d'interface endommageable (pour interf acier/beton)
        case(107)
            call interf_sig_endo(vhep,vsi,vnl,vin,vh,vb,vnl0,vprel,vdle,iloi,ie,ipg)

        !----------------------------------------------------------------------------------
        !--- Cas du modele d'interface endommageable puis frottante (pour BRF(elem line))
        !    Combinaison des lois 105 et 108
        case(158)
            if ((ietatpg(ie,1)==4).and.(ietatpg(ie,ipg)==5)) then
                ietatpg(ie,ipg) = 4
                vin(1:3) = 0.d0
            endif

            if (all(ietatpg(ie,:)==0).or.all(ietatpg(ie,:)==3)) then

                !--------------------------------------------------------------------------
                !--- Initialisations
                vsi2  = 0.D0
                vin2  = 0.D0
                vin2(4:5) = vin(4:5)
                vhep2 = 0.D0

                !--------------------------------------------------------------------------
                !--- Comportement de la matrice cimentaire
                iloi1 = 105
                call interf_sig_endo(vhep,vsi,vnl,vin,vh,vb,vnl0,vprel,vdle,iloi1,ie,ipg)

                if (all(ietatpg(ie,:)==1).or.all(ietatpg(ie,:)==2)) then
                    ietatpg(ie,:) = 5
                    ietatpg(ie,1) = 4
                    vin(1:3) = 0.d0
                endif

                !--------------------------------------------------------------------------
                !--- Comportement du pontage par fibres
                iloi2 = 108
                call interf_sig_endo(vhep2,vsi2,vnl,vin2,vh,vb,vnl0,vprel,vdle,iloi2,ie,ipg)

                vsi = vsi+vsi2
                vhep = vhep+vhep2

            !elseif (all(ietatpg(ie,:)==4)) then
            elseif (ietatpg(ie,ipg)==4) then
                iloi2 = 108
                call interf_sig_endo(vhep,vsi,vnl,vin,vh,vb,vnl0,vprel,vdle,iloi2,ie,ipg)
            else
                call interf_sig(vhep,vsi,vnl,vh,vb,vnl0,vprel,vdle,iloi,ie,ipg)
            endif

        case default
            stop "Elem_sig : cas non encore implante"

    end select

    deallocate(vh)

end subroutine elem_sig

!------------------------------------------------------------------------------------------------------

end module contraintes

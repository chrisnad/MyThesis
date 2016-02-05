!------------------------------------------------------------------------------
! MODULE: contraintes
!
!> @author JL Tailhan
!
!> @brief
!> Fonction d'aiguillage pour le calcul des contraintes.
!>
!------------------------------------------------------------------------------
module contraintes


contains

!---------------------------------------------------------!
!   Fonction d'aiguillage pour le calcul des contraintes  !
!---------------------------------------------------------!
   subroutine elem_sig(vsi,vnl,vin,ie,vhep,vb,vnl0,vprel,vdle,ipg)

    use initialisation, only : init_mat, init_vec
    use aiguillage, only : loi
    use element_interface, only : interf_sig, interf_sig_endo
    use fissuration, only : fiss_sig, fiss_sig_elas, fiss_sig_endo
    use acier, only : acier_sig
    use fissuration_BA, only : fiss_sig_endo_ortho
    use variables, only : ktypel, nomtype, ietatpg
    use lib_elem, only : elem_hooke
    use aiguillage, only : modul

    implicit none

    ! Variables IN
    real*8, dimension(:,:), intent(in) :: vb
    integer, intent(in) :: ie, ipg
    real*8, dimension(:), intent(in) :: vnl0, vdle
    real*8, dimension(:), intent(in) :: vprel

    ! Variables OUT
    real*8, dimension(:), intent(inout) :: vsi
    real*8, dimension(:), intent(out) :: vnl, vin
    real*8, dimension(:,:), allocatable, intent(out) :: vhep

    real*8, dimension(:,:), allocatable :: vh
    integer :: icomp, iloi

    ! Variables permettant l'utilisation de lois complexes (combinant des lois simples)
    ! loi 158 = 105 + 108 (108 ne pouvant fonctionner seule)
    real*8, dimension(size(vb,1),size(vb,1)) :: vhep2
    real*8, dimension(size(vsi,1)) :: vsi2
    real*8, dimension(size(vin,1)) :: vin2
    integer :: iloi1,iloi2

    ! ----- On recupere le type de comportement --------
    icomp = int(vprel(1))   ! icomp pour chaque elements

    ! ----- On verifie la compatibilite : Comportement / Loi ---
    if (icomp /= 1) then
       call loi(iloi,icomp,ie)!,ipg)
    else
       iloi = 1
    end if

    !----- matrice d'elasticite vh, en fonction du type d'element -----
    call elem_hooke(vh,nomtype(ktypel(ie)),vprel)

    call modul(ie,ipg,vh,vprel)

    call init_mat(vhep,size(vh,1),size(vh,2))

    select case (iloi)

        case(1,2)
            ! Cas de l'elasticite lineaire (isotrope=1, orthotrope=2)
            vsi = matmul(vh,matmul(vb,vdle))
            vhep = vh

        case(5)
            ! Cas de l'elasticite lineaire (iloi=1,2) et elasticite fragile (iloi=5)
            call fiss_sig_elas(vsi,vh,vb,vdle,iloi,ie)
            vhep = vh

        case(100,101)
            ! Modele de fissuration d'interface elastique fragile
            call interf_sig(vhep,vsi,vnl,vh,vb,vnl0,vprel,vdle,iloi,ie,ipg)

        case(107)
            ! Modele d'interface endommageable (pour interf acier/beton)
            call interf_sig_endo(vhep,vsi,vnl,vin,vh,vb,vnl0,vprel,vdle,iloi,ie,ipg)

        case(102, 105, 106, 108)
            ! Modeles simples d'interface endommageable puis frottante (pour beton et BRF)
            if (any(ietatpg(ie,:)==0).or.any(ietatpg(ie,:)==3)) then
            !if (all(ietatpg(ie,:)==0).or.all(ietatpg(ie,:)==3)) then
                call interf_sig_endo(vhep,vsi,vnl,vin,vh,vb,vnl0,vprel,vdle,iloi,ie,ipg)
            else
                call interf_sig(vhep,vsi,vnl,vh,vb,vnl0,vprel,vdle,iloi,ie,ipg)
            end if

        case(158)
            ! Modele complexe d'interface endommageable puis frottante (pour BRF(elem line))
            ! Combinaison des lois 105 et 108
            if ((ietatpg(ie,1)==4).and.(ietatpg(ie,ipg)==5)) then
                ietatpg(ie,ipg) = 4
                vin(1:3) = 0.d0
            end if

            if (all(ietatpg(ie,:)==0).or.all(ietatpg(ie,:)==3)) then
                ! Initialisations
                vsi2  = 0.D0
                vin2  = 0.D0
                vin2(4:5) = vin(4:5)
                vhep2 = 0.D0

                ! Comportement de la matrice cimentaire
                iloi1 = 105
                call interf_sig_endo(vhep,vsi,vnl,vin,vh,vb,vnl0,vprel,vdle,iloi1,ie,ipg)

                if (all(ietatpg(ie,:)==1).or.all(ietatpg(ie,:)==2)) then
                    ietatpg(ie,:) = 5
                    ietatpg(ie,1) = 4
                    vin(1:3) = 0.d0
                end if

                ! Comportement du pontage par fibres
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
            end if

        case(11)
            call acier_sig(vhep,vsi,vnl,vh,vb,vnl0,vprel,vdle,iloi,ie,ipg)

        case(12, 13, 14)
            call fiss_sig(vsi,vnl,vh,vb,vnl0,vprel,vdle,iloi,ie)
            vhep = vh

!        case(12)
!            call fiss_sig_plas(vsi,vnl,vh,vb,vnl0,vprel,vdle,iloi,ie)
!            vhep = vh

        case(15)
            ! Modele macro de fissuration (endomagement)
            call fiss_sig_endo(vsi,vnl,vin,vh,vb,vprel,vdle,iloi,ie)
            vhep = vh

        case(25, 26, 27)
            ! Modele macro de fissuration (endomagement) BA orthotrope
            call fiss_sig_endo_ortho(vsi,vnl,vin,vh,vb,vprel,vdle,iloi,ie)
            vhep = vh

    case default
        stop "Elem_sig : cas non encore implante"

    end select

    deallocate(vh)

   end subroutine elem_sig

end module contraintes

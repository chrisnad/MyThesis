!------------------------------------------------------------------------------
! MODULE: fissuration
!
!> @author JL Tailhan
!
!> @brief
!> Ensemble des routines relative au(x) modèle(s) macroscopique(s)
!> de fissuration.
!>
!------------------------------------------------------------------------------

module fissuration


contains


!---------------------------------------------------------------------------
!> @author JL Tailhan
!> @brief Routine de calcul des contraintes (en semi implicite).
!
!> @details
!> ### DESCRIPTION:
!> Routine de calcul des contraintes basée sur un algorithme
!> semi implicite.
!>
!>--------------------------------------------------------------------------
subroutine fiss_sig(vsig,vnle,vh,vb,vnl0,vprel,vdle,iloi,ie)

        use variables, only : dime, inict, irloc, iedng, inorm
        use initialisation, only : init_vec, init_mat
        use utilitaire, only : princ

        implicit none

        !--- Variables IN
        real*8, dimension(:,:), intent(in) :: vh, vb
        real*8, dimension(:), intent(in) :: vnl0, vdle, vprel
        integer, intent(in) :: iloi, ie

        !--- Variables OUT
        real*8, dimension(:), intent(inout) :: vsig
        real*8, dimension(:), intent(out) :: vnle

        !--- Quantites globales
        real*8, dimension(size(vh,1),size(vh,2)) ::  P, Q
        real*8, dimension(dime,dime) :: V
        real*8, dimension(dime) :: vsp
        integer :: it, nc1
        character(len=5) :: calcu

        !--- Quantites principales
        real*8, dimension(size(vh,1),size(vh,2)) :: vhloc
        real*8 :: vcrit, dlam
        logical :: nconv

        ! Quantites principales
        real*8, dimension(size(vb,1)) :: vepspl, vsigm, vsigi
        real*8, dimension(size(vh,1)) :: vsigmloc, vsigiloc, vepsplloc
        real*8, dimension(size(vh,1)) :: vdfdsigloc, vdgdsigloc

!********************************************************!

        !----- Recuperation des grandeurs de non-linearite
        vepspl = vnl0(1:size(vb,1))            ! deformations plastiques initiales

        !----- Calcul de sigma + dsigma
        vsigm = matmul(vh,(matmul(vb,vdle)-vepspl))

        !----- Calcul et rangement par ordre decroissant des contraintes principales
        if (inict) then
            vsigm = vsig + vsigm
        end if

        !----- Adaptation pour le pilotage par l'element le plus dangereux
        if (ie==iedng) then
            V=reshape(irloc(ie,:),(/dime,dime/))
        else
            vsp = 0.d0
            call princ(vsigm,vsp,V)
        end if

        !----- Recuperation de la matrice de changement de repere global --> local principal
        nc1 = size(vb,1)
        P = fiss_changement_repere(V,nc1,1,1) !--- matrice de changement de repere pour les contraintes
        Q = fiss_changement_repere(V,nc1,2,2) !--- matrice inverse de changement de repere pour les deformations
        vsigi = 0.d0

        ! Dans le repere principal...
        vhloc = matmul(P,matmul(vh,Q))
        vsigmloc = matmul(P,vsigm)
        vsigiloc = matmul(P,vsigi)
        vepsplloc = matmul(P,vepspl)

        !----- Calcul des contraintes verifiant le critere
        calcu = 'D1F'
        call fiss_crit(iloi,vprel,vsigmloc,ie,vcrit,calcu,vdfdsigloc,vdgdsigloc)
        calcu = 'D0F'

        !----- Calcul du multiplicateur plastique (algo semi-implicite)
        it = 0
        vsigiloc = vsigmloc
        dlam = 0.0d0
        nconv = .true.

        if (vcrit>0.d0) then


        ! ajouter la modification de vh dans le cas de cisaillement

            do while (nconv)

                dlam = dlam + vcrit/dot_product(vdfdsigloc,matmul(vhloc,vdgdsigloc))
                vsigmloc = vsigiloc - dlam*matmul(vhloc,vdgdsigloc)
                call fiss_crit(iloi,vprel,vsigmloc,ie,vcrit,calcu)
                it = it + 1
                if (abs(vcrit) < 1.e-10) nconv = .false.
                if (it>50) stop 'fissuration: non convergence dans le calcul des contraintes'

            end do

            !if (iloi==12) inorm(ie,:) = v(:,1)

        end if

        vepsplloc = vepsplloc + dlam*vdgdsigloc

        !----- On retourne dans le repere global...
        vsigm = matmul(Q,vsigmloc)
        vepspl = matmul(Q,vepsplloc)

        !----- Sorties :

        vsig = vsigm
        vnle = vepspl

end subroutine fiss_sig


!---------------------------------------------------------------------------
!> @author JL Tailhan
!> @brief Routine de calcul des contraintes (en semi implicite / plast).
!
!> @details
!> ### DESCRIPTION:
!> Routine de calcul des contraintes elastiques (y compris elastique-fragile)
!>
!>--------------------------------------------------------------------------
subroutine fiss_sig_elas(vsig,vh,vb,vdle,iloi,ie)


    use variables, only : iedng, ietat, dime, irloc
    use math

    implicit none

    !--- Variables IN
    real*8, dimension(:,:), intent(in) :: vb
    real*8, dimension(:,:), intent(inout) :: vh
    real*8, dimension(:), intent(in) :: vdle
    integer, intent(in) :: iloi, ie

    !--- Variables OUT
    real*8, dimension(:), intent(inout) :: vsig

    !****************************************************************************
    !--- Variables locales
    real*8, dimension(size(vh,1),size(vh,2)) :: P1, P2, vhloc, vhpetit
    real*8, dimension(dime,dime) :: V
    real*8, dimension(size(vsig,1)) :: vsiloc
    real*8  :: Dn, muG
    integer :: nc1
    logical :: vcrit

    !****************************************************************************
    vhpetit = vh/1.d15

    !----- Cas de(s) loi(s) fragile(s)
    select case(iloi)

        case(5)    !----- Loi elastique fragile (beton fissurant)

            !----- Test pour activation de la rupture fragile
            vcrit = .false.
            if (ie == iedng) then
                !----- Pour l'element le plus dangereux
                vcrit = .true.
            else
                if (ietat(ie) /= 0) then
                    !----- Pour les elements deja rompus
                    vcrit = .true.
                end if
            end if

            !----- Declanchement de la rupture fragile
            if (vcrit .eqv. .true.) then

                !--    Constitution du repere de la fissure
                V = reshape(irloc(ie,:),(/dime,dime/))

                !--    Changement de repere pour revenir dans le rep global
                nc1 = size(vb,1)
                P1 = fiss_changement_repere(V,nc1,1,1) !--- matrice de changt de rep pr les contraintes
                P2 = fiss_changement_repere(V,nc1,2,2) !--- matrice inverse de changt de rep pr les deformations
                vhloc = matmul(P1,matmul(vh,P2))

                !--    Calcul de vhloc pour prendre en compte la fissure (repere local)
                if (ie/=iedng) then
                    muG = 0.3d0
                    if (dime==2) then
                        vhloc(3,3) = muG*vhloc(3,3)
                    elseif (dime==3) then
                        vhloc(4,4) = muG*vhloc(4,4)
                        vhloc(5,5) = muG*vhloc(5,5)
                        vhloc(6,6) = muG*vhloc(6,6)
                    else
                        print *, 'Cas non prevu'
                        stop
                    end if
                endif

                !-- Prise en compte de l'etat de l'element
                Dn = 0.d0
                if (ietat(ie)>0) ietat(ie)=2 ! <= pour traiter la refermeture
                if (ietat(ie)==-2) then
                !elseif (ietat(ie)==1) then
                    Dn = 1.d0
                endif
                vhloc = (1.d0 - Dn)*vhloc

                !-- Calcul de l'etat de contrainte ds rep fissure
                P2 = fiss_changement_repere(V,nc1,2,1) !--- matrice de changt de rep pr les deformations
                vsiloc = matmul(vhloc,matmul(P2,matmul(vb,vdle)))

!if (iedng==1736) affiche=.true.
!if ((affiche)) then !.and.(ie==2217)) then
!    print*,ie,ietat(ie),vsiloc(1),Dn
!endif

                !-- Correction des modules en fonction de l'etat de l'element
                if (ietat(ie)/=-2) then
                    if (vsiloc(1) >= 0.d0) then
                        Dn = 1.d0 ! <- l'element est ouvert !
                        ietat(ie) = 1
                    else
                        ietat(ie) = 2
                    endif
                endif
                vhloc = (1.d0 - Dn)*vhloc

                !-- Calcul de vh dans le repere global
                nc1 = size(vb,1)
                P1 = fiss_changement_repere(V,nc1,1,2) !--- matrice inverse de changt de rep pr les contraintes
                P2 = fiss_changement_repere(V,nc1,2,1) !--- matrice de changt de rep pr les deformations
                vh = matmul(P1,matmul(vhloc,P2))

            end if

        case default

            stop 'fiss_sig_elas : loi non implantee'

    end select

    !----- Calcul des contraintes elastiques
    vsig = matmul(vh,matmul(vb,vdle))
    vh   = vh + vhpetit


!if ((affiche).and.(ietat(ie)/=0)) then !.and.(ie==2217)) then
!                vsiloc = matmul(vhloc,matmul(P2,matmul(vb,vdle)))
!    print*,ie,ietat(ie),vsiloc(1),Dn
!    !pause
!endif

end subroutine fiss_sig_elas





!---------------------------------------------------------------------------
!> @author JL Tailhan
!> @brief Routine de calcul des contraintes (en semi implicite / plast).
!
!> @details
!> ### DESCRIPTION:
!> Routine de calcul des contraintes basée sur un algorithme
!> semi implicite (modèle plastique).
!>
!>--------------------------------------------------------------------------
subroutine fiss_sig_plas(vsig,vnle,vh,vb,vnl0,vprel,vdle,iloi,ie)

        use variables, only : dime, inict, irloc, iedng, inorm, ietat
        use initialisation, only : init_vec, init_mat
        use utilitaire, only : princ

        implicit none

        !--- Variables IN
        real*8, dimension(:,:), intent(in) :: vh, vb
        real*8, dimension(:), intent(in) :: vnl0, vdle, vprel
        integer, intent(in) :: iloi, ie

        !--- Variables OUT
        real*8, dimension(:), intent(inout) :: vsig
        real*8, dimension(:), intent(out) :: vnle

        !--- Quantites globales
        real*8, dimension(size(vh,1),size(vh,2)) ::  P, Q
        real*8, dimension(dime,dime) :: V
        real*8, dimension(dime) :: vsp
        integer :: it, nc1
        character(len=5) :: calcu

        !--- Quantites principales
        real*8, dimension(size(vh,1),size(vh,2)) :: vhloc
        real*8 :: vcrit, dlam
        logical :: nconv

        ! Quantites principales
        real*8, dimension(size(vb,1)) :: vepspl, vsigm, vsigi
        real*8, dimension(size(vh,1)) :: vsigmloc, vsigiloc, vepsplloc
        real*8, dimension(size(vh,1)) :: vdfdsigloc, vdgdsigloc

!********************************************************!

        !----- Recuperation des grandeurs de non-linearite
        vepspl = vnl0(1:size(vb,1))            ! deformations plastiques initiales

        !----- Calcul de sigma + dsigma
        vsigm = matmul(vh,(matmul(vb,vdle)-vepspl))

        !----- Calcul et rangement par ordre decroissant des contraintes principales
        if (inict) then
            vsigm = vsig + vsigm
        end if

        !----- Adaptation pour le pilotage par l'element le plus dangereux
        if (ie==iedng) then
            V=reshape(irloc(ie,:),(/dime,dime/))
        else
            vsp = 0.d0
            call princ(vsigm,vsp,V)
        end if

        !----- Recuperation de la matrice de changement de repere global --> local principal
        nc1 = size(vb,1)
        P = fiss_changement_repere(V,nc1,1,1) !--- matrice de changement de repere pour les contraintes
        Q = fiss_changement_repere(V,nc1,2,2) !--- matrice inverse de changement de repere pour les deformations
        vsigi = 0.d0

        ! Dans le repere principal...
        vhloc = matmul(P,matmul(vh,Q))
        vsigmloc = matmul(P,vsigm)
        vsigiloc = matmul(P,vsigi)
        vepsplloc = matmul(P,vepspl)

        !----- Calcul des contraintes verifiant le critere
        calcu = 'D1F'
        call fiss_crit(iloi,vprel,vsigmloc,ie,vcrit,calcu,vdfdsigloc,vdgdsigloc)
        calcu = 'D0F'

        if ((ie==iedng).or.(ietat(ie)==3)) then
            !----- Calcul du multiplicateur plastique (algo semi-implicite)
            it = 0
            vsigiloc = vsigmloc
            dlam = 0.0d0
            nconv = .true.

            if (vcrit>0.d0) then


            ! ajouter la modification de vh dans le cas de cisaillement

                do while (nconv)

                    dlam = dlam + vcrit/dot_product(vdfdsigloc,matmul(vhloc,vdgdsigloc))
                    vsigmloc = vsigiloc - dlam*matmul(vhloc,vdgdsigloc)
                    call fiss_crit(iloi,vprel,vsigmloc,ie,vcrit,calcu)
                    it = it + 1
                    if (abs(vcrit) < 1.e-10) nconv = .false.
                    if (it>50) stop 'fissuration: non convergence dans le calcul des contraintes'

                end do

                !if (iloi==12) inorm(ie,:) = v(:,1)

            end if

            vepsplloc = vepsplloc + dlam*vdgdsigloc

        end if

        !----- On retourne dans le repere global...
        vsigm = matmul(Q,vsigmloc)
        vepspl = matmul(Q,vepsplloc)

        !----- Sorties :

        vsig = vsigm
        vnle = vepspl

end subroutine fiss_sig_plas


!!---------------------------------------------------------------------------
!!> @author JL Tailhan
!!> @brief Routine de calcul des contraintes pour modèles d'endommagement.
!!
!!> @details
!!> ### DESCRIPTION:
!!> Routine de calcul des contraintes pour modèles d'endommagement
!!>
!!>--------------------------------------------------------------------------
!subroutine fiss_sig_endo2(vsig,vnle,vin,vh,vb,vnl0,vprel,vdle,iloi,ie,ipg)

!        use variables, only : dime, inict, irloc, iedng, ietat, inorm
!        use initialisation, only : init_vec, init_mat
!        use utilitaire, only : princ

!        implicit none

!        !--- Variables IN
!        real*8, dimension(:,:), intent(in) :: vb
!        real*8, dimension(:,:), intent(inout) :: vh
!        real*8, dimension(:), intent(in) :: vnl0, vdle, vprel
!        integer, intent(in) :: iloi, ie, ipg

!        !--- Variables OUT
!        real*8, dimension(:), intent(inout) :: vsig, vin
!        real*8, dimension(:), intent(out) :: vnle

! !****************************************************************************

!        !--- Quantites globales
!        real*8, dimension(size(vh,1),size(vh,2)) ::  P, Q
!        real*8, dimension(dime,dime) :: V  !, V0, VVsig
!        real*8, dimension(dime) :: vsp  !, vsp0
!        integer :: nc1, id  !, NROT

!        !--- Quantites principales
!        real*8, dimension(size(vh,1),size(vh,2)) :: vhloc
!        logical :: vcrit
!        !integer, dimension(dime) :: ilocp

!        ! Quantites principales
!        real*8, dimension(size(vb,1)) :: vsigm, vsigi, veps
!        real*8, dimension(size(vh,1)) :: vsigmloc, vsigiloc, vepsloc

!        ! Variables pour le modele d'endommagement
!        real*8 :: Dini, var0, var, vars, varl, D, sigma1, epsilon1
!        real*8 :: vcrit1, varcrit, varc, F_var, RT, W, E, nu

!!********************************************************!

!        !---- Initialisations
!        Dini = vin(1) ! endommagement initial
!        var0 = vin(2) ! seuil initial
!        vars = vin(3) ! seuil courant
!        varl = 1000. ! deformation limite


!        ! Endommagement initial
!        D = Dini

!        !----- Recuperation des parametres pour le critere
!        id = 6
!        if (dime == 3) id = id-1
!        RT = vprel(id)  ! Contrainte limite
!        W = vprel(id+1) ! Energie de post-fissuration
!        E = vprel(id-2)
!        nu = vprel(id-1)

!        ! Etat de contrainte et deformation
!        vsigm = matmul((1.d0-D)*vh,(matmul(vb,vdle)))
!        veps = (matmul(vb,vdle))

!        !----- Calcul et rangement par ordre decroissant des contraintes principales
!        if (inict) then
!            vsigm = vsig + vsigm
!        end if

!        !----- Adaptation pour le pilotage par l'element le plus dangereux
!        if (ie==iedng) then
!            V=reshape(irloc(ie,:),(/dime,dime/))
!        else
!            vsp = 0.d0
!            call princ(vsigm,vsp,V)
!        end if

!        !----- Recuperation de la matrice de changement de repere global --> local principal
!        nc1 = size(vb,1)
!        P = fiss_changement_repere(V,nc1,1,1) !--- matrice de changement de repere pour les contraintes
!        Q = fiss_changement_repere(V,nc1,2,2) !--- matrice inverse de changement de repere pour les deformations
!        vsigi = 0.d0

!        ! Dans le repere principal...
!        vhloc = matmul(P,matmul(vh,Q))
!        vsigmloc = matmul(P,vsigm)
!        vsigiloc = matmul(P,vsigi)
!        vepsloc = matmul(P,veps)

!        !----- Recuperation de la contrainte principale et de la deformation correspondante
!        sigma1 = vsigmloc(1)
!        epsilon1 = vepsloc(1)

!        !---- Variqble de pilotqge de l'endommagement
!        var = epsilon1

!        !----- Difference entre les contraintes pour verification du critere
!        !vcrit1 = sigma1 - 1.0001*RT   ! en ouverture pure

!        vcrit = .false.
!        !if (vcrit1 >= 0.d0 .and. ietat(ie) == 0) then
!        if (ie == iedng) then
!            !----- Changement "vcrit"
!            vcrit = .true.
!            !----- Calcul du seuil de deformation initial
!            if(dime == 2) then
!                var0 = (sigma1 - nu*vsigmloc(2))/E !contraintes planes uniquement
!            else if(dime == 3) then
!                var0 = (sigma1 - nu*(vsigmloc(2)+vsigmloc(3)))/E
!            end if
!            !----- Initialisation du seuil de deformation courant
!            vars = var0
!            F_var = 0.d0
!            !----- Changement etat de l'element
!            ietat(ie) = 3

!        else
!            if(ietat(ie) == 3) then
!                vcrit = .true.
!            end if
!        end if

!        !----- Declanchement et evolution de l'endommagement
!        if(vcrit .eqv. .true.) then

!            inorm(ie,:) = v(:,1)

!            !---- Valeur critique de la variable d'endommagement
!            varc = 2.*W/RT

!            if(varc <= var0) then
!                W = 0.5*var0*RT
!                varc = var0
!                varl = var0
!            else
!                varl = sqrt(2*W/E) + 1.*(varc-sqrt(2*W/E))
!                !varl=1.75d-03
!            end if

!            varcrit = varc
!            if(varc > varl) then
!                varl = max(varl,sqrt(2*W/E))
!                varc = var0 - ((varl-var0)**2)/(varc-2*varl+var0)
!                varcrit = varl
!            end if

!            if(var .ge. vars) then

!                !---- Actualisation du seuil d'endommagement
!                vars = var
!                !---- Loi d'evolution d'endommagement
!                !if(varc < var0) then
!                !     varc = var0
!                if(var0 == varc) then
!                    F_var = 1.d0
!                else
!                    F_var = (var - var0)/(varc - var0 + 1.d-20)
!                    if(F_var > 1.d0) then
!                        F_var = 1.d0
!                        !else if(F_var < 0.d0) then
!                        !F_var = 0.d0
!                    end if
!                end if
!                !F_var = floor(10000.d0 * F_var) / 10000.d0

!            end if

!            !---- Calcul de l'endommagement (updated)
!            D = 1.d0 - var0/var * (1.d0 - F_var)

!            if(D<0.d0) D=0.d0
!            D = max(Dini,D)
!            D = min(D,.9999d0)
!            if((var>varcrit) .or. (D .ge. 0.9999)) then
!                D=.9999d0
!                ietat(ie) = 1
!                inorm(ie,:) = v(:,1)
!            end if

!        end if

!        D = floor(1.d4 * D) * 1.d-4

!        ! Rigidité
!        vsig = matmul((1.d0-D)*vh,(matmul(vb,vdle)))
!        vh = (1.d0-D)*vh

!        vin(1) = D
!        vin(2) = var0
!        vin(3) = vars

!end subroutine fiss_sig_endo2


!---------------------------------------------------------------------------
!> @author JL Tailhan
!> @brief Routine de calcul des contraintes pour modèles d'endommagement.
!
!> @details
!> ### DESCRIPTION:
!> Routine de calcul des contraintes pour modèles d'endommagement
!>
!>--------------------------------------------------------------------------
subroutine fiss_sig_endo(vsig,vnle,vin,vh,vb,vprel,vdle,iloi,ie)

    use variables, only : dime, iedng, ietat, inorm, irloc
    use initialisation, only : init_vec, init_mat
    use utilitaire, only : princ

    implicit none

    !--- Variables IN
    real*8, dimension(:,:), intent(in) :: vb
    real*8, dimension(:,:), intent(inout) :: vh
    real*8, dimension(:), intent(in) :: vdle, vprel
    integer, intent(in) :: iloi, ie

    !--- Variables OUT
    real*8, dimension(:), intent(inout) :: vsig, vin
    real*8, dimension(:), intent(out) :: vnle


 !****************************************************************************

    !--- Quantites globales
    real*8, dimension(size(vh,1),size(vh,2)) ::  P, Q, vhpetit
    real*8, dimension(dime,dime) :: V  !, V0, VVsig
    real*8, dimension(dime) :: vsp  !, vsp0
    integer :: nc1, id  !, NROT

    !--- Quantites principales
    logical :: vcrit
    !integer, dimension(dime) :: ilocp

    ! Quantites principales
    real*8, dimension(size(vb,1)) :: veps, vepsloc
    real*8, dimension(size(vh,1)) :: vsigm, vsigmloc

    ! Variables pour le modele d'endommagement
    real*8 :: Dini, var0, var, vars, varl, D, sigma1, epsilon1
    real*8 :: varcrit, varc, F_var, RT, W, E, nu
    real*8 :: ouv, Dnul

!********************************************************!

    !----- Calcul de l'ouverture de fissure initiale (pour non-interpenetration)
    ouv = 0.D0 ; Dnul = 1.d0
    if (ietat(ie) == 3) ouv = fiss_ouv(ie) ! Calcul de l'ouverture fait selon la normale au pas precedent
    if (ouv < 0.D0) Dnul = 0.d0

    vhpetit = vh/1.d15

    select case(iloi)

        case(15)  !----- Loi d'endommagement (beton fissurant)

            !---- Initialisations
            Dini = vin(1) ! endommagement initial
            var0 = vin(2) ! seuil initial
            vars = vin(3) ! seuil courant
            varl = 1000. ! deformation limite

            !----- Etat de contrainte et de deformation
            veps  = matmul(vb,vdle)
            vsigm = matmul((1.d0-Dnul*Dini)*vh,veps)

            !----- Parametres specifiques à la loi
            id = 6
            if (dime == 3) id = id-1
            E = vprel(id-2) ! Module d'Young
            nu = vprel(id-1)! Coefficient de Poisson
            RT = vprel(id)  ! Resistance a la traction
            W = vprel(id+1) ! Energie post-fissuration

            !----- Recuperation de la matrice de changement de repere global --> local principal
            !----- pour l'element le plus dangereux on recupere le repere local a la fissure : V
            !if (ie==iedng) then
            !    V=reshape(irloc(ie,:),(/dime,dime/))
            !else
            !----- sinon on recalcule ce repere
                vsp = 0.d0
                call princ(vsigm,vsp,V)
            !end if
            !----- puis on calcule les matrices de changement de repere (Global->principal)
            nc1 = size(vb,1)
            P = fiss_changement_repere(V,nc1,1,1) !--- matrice de changement de repere pour les contraintes
            !Q = fiss_changement_repere(V,nc1,2,2) !--- matrice inverse de changement de repere pour les deformations

            !----- Dans le repere principal...
            vsigmloc = matmul(P,vsigm)
            vepsloc = matmul(P,veps)

            !----- Recuperation de la contrainte principale et de la deformation correspondante
            sigma1 = vsigmloc(1)
            epsilon1 = vepsloc(1)

            !----- Test pour activation du calcul de l'endommagement
            vcrit = .false.
            if (ie == iedng) then
                !----- Changement "vcrit"
                vcrit = .true.
                !----- Calcul du seuil de deformation initial
                if(dime == 2) then
                    var0 = (sigma1 - nu*vsigmloc(2))/E !contraintes planes uniquement
                else if(dime == 3) then
                    var0 = (sigma1 - nu*(vsigmloc(2)+vsigmloc(3)))/E
                end if
                !----- Initialisation du seuil de deformation courant
                vars = var0
                F_var = 0.d0
                !----- Changement etat de l'element
                ietat(ie) = 3

            else
                if(ietat(ie) == 3) then
                    vcrit = .true.
                end if
            end if

            !---- Valeur critique de la variable d'endommagement
            varc = 2.*W/RT

            if(varc <= var0) then
                W = 0.5*var0*RT
                varc = var0
                varl = var0
            else
                varl = sqrt(2*W/E) + 1.*(varc-sqrt(2*W/E)) !! a modifier !!
                !varl=1.75d-03
            end if

            varcrit = varc
            if(varc > varl) then
                varl = max(varl,sqrt(2*W/E))
                varc = var0 - ((varl-var0)**2)/(varc-2*varl+var0)
                varcrit = varl
            end if

            !----- Declanchement et evolution de l'endommagement
            if(vcrit .eqv. .true.) then

                !---- Recuperation de la normale a la fissure
                !inorm(ie,:) = V(:,1)

                !---- Variable de pilotage de l'endommagement
                var = Dnul * epsilon1

                if(var .ge. vars) then

                    !---- Actualisation du seuil d'endommagement
                    vars = var

                    !---- Loi d'evolution d'endommagement
                    !if(varc < var0) then
                    !     varc = var0
                    if(var0 == varc) then
                        F_var = 1.d0
                    else
                        F_var = (var - var0)/(varc - var0 + 1.d-20)
                        if(F_var > 1.d0) then
                            F_var = 1.d0
                            !else if(F_var < 0.d0) then
                            !F_var = 0.d0
                        end if
                    end if
                    !F_var = floor(10000.d0 * F_var) / 10000.d0

                    !---- Calcul de l'endommagement (updated)
                    D = 1.d0 - var0/var * (1.d0 - F_var)

                else

                    D = Dini

                end if

            else

                D = Dini

            end if

            if(D < 0.d0) D = 0.d0
            D = max(Dini,D)
            D = min(D,1.d0)

            if((var>varcrit) .or. (D == .9999d0)) then
                D = .9999d0
                ietat(ie) = 1
                inorm(ie,:) = V(:,1) ! la valeur est déjà stockee dans inorm (cf. plus haut)
            end if

            ! Contraintes et Rigidité
            vsig = matmul((1.d0-Dnul*D)*vh,(matmul(vb,vdle)))
            vnle(1:size(veps,1)) = (-Dnul*D*matmul(vb,vdle))

            vh = (1.d0-Dnul*D)*vh
            vh = vh + vhpetit

            vin(1) = D
            vin(2) = var0
            vin(3) = vars

        case default

            stop 'fiss_sig_endo : cas non implante'

    end select


end subroutine fiss_sig_endo


!---------------------------------------------------------------------------
!> @author JL Tailhan
!> @brief Routine de calcul des critères et dérivées des critères.
!
!> @details
!> ### DESCRIPTION:
!> Routine de calcul des critères de rupture et des dérivées des critères.
!> Les critères sont exprimés dans le repère principal des contraintes.
!> N.B. : DANS LE REPERE PRINCIPAL (sigma1, sigma2, sigma3)
!> => sigma1 > sigma2 > sigma3
!>
!>--------------------------------------------------------------------------
subroutine fiss_crit(iloi,vprel,vsig,ie,vcrit,calcu,vdfdsig,vdgdsig)


        use variables, only : dime, pi, &
                           &  ietat, iedng
        use math, only : Jacobi, INV, norme, find_vec
        use initialisation, only : init_vec, init_mat

        implicit none

        !--- Variables IN
        real*8, dimension(:), intent(in) :: vsig
        real*8, dimension(:), intent(in) :: vprel
        integer, intent(in) :: iloi, ie
        character(len=5), intent(in), optional :: calcu

        !--- Variables OUT
        real*8, intent(out) :: vcrit
        real*8, dimension(:), intent(out), optional :: vdfdsig, vdgdsig

        !--- Variables locales procedure
        real*8 :: s1, s3, vcrit1
        integer :: id
        logical :: ideriv

        !--- Loi de Rankine (case 12)
        real*8 :: RT

        !--- Loi de Von Mises (case 11)
        !real*8, dimension(:), allocatable :: vsigd, vun
        !real*8 :: cnu, a, trvsig, vmis


!********************************************************!

        !----- Test(s) de compatibilite des arguments d'entree
        ideriv = .false.
        if (calcu == 'D1F') then
            ideriv = .true.
            vdfdsig = 0.0d0
            vdgdsig = 0.0d0
        end if

!********************************************************!

 select case(iloi)

    case(12)     !----- Loi de Tresca "orientee"

        !----- Recuperation des parametres
        id = 6
        if (dime == 3) id = id-1
        RT = vprel(id)                        ! Contrainte limite

        !----- Recuperation des contraintes principales 1 et 3
        s1 = vsig(1)
        s3 = vsig(3)

        !----- Criteres de rupture de l'element
        !---
        !--- dans le cas du calcul iteratif des contraintes
        if (calcu=='D0F')then
                vcrit = s1 - RT
            return
        end if
        !---
        !--- dans le cas de la detection de la rupture
        vcrit1 = s1 - RT                      ! ... en ouverture pure
        !vcrit2 = s1 - s3 - 2.d0 * C           ! ... en cisaillement

        vcrit = -1.d0

        !----- Calcul des criteres et des derivees
        if ((iedng==ie) .or. (ietat(ie) /=0 )) then

            !----- Rupture en ouverture
            if (vcrit1 >= 0.d0) then

                if (ie==iedng) ietat(ie) = 3 ! ... en plasticite l'element n'est pas encore fissure

                !!!!! ATTENTION : PARTIE A ADAPTER POUR PRENDRE EN COMPTE LE CISAILLEMENT !!!!
                vcrit = vcrit1
                if (ideriv) then
                    if (dime == 2) then
                        vdfdsig = (/ 1.d0, 0.d0, 0.d0 /)
                        vdgdsig = (/ 1.d0, 0.d0, 0.d0 /)
                    elseif (dime == 3) then
                        vdfdsig = (/ 1.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0 /)
                        vdgdsig = (/ 1.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0 /)
                    end if
                end if
            end if

        end if

!********************************************************!

!         case(11)     !----- Loi de Von Mises
!
!         !----- Recuperation des parametres de la loi
!         id = 6
!         if (dime == 3) id = id-1
!         RT = vprel(id)                        ! Contrainte limite
!         cnu = vprel(id-1)                     ! Coefficient de Poisson
!
!        !----- Calcul de la contrainte de Von Mises
!        call init_vec(vun,size(vsig))
!        call init_vec(vsigd,size(vsig))
!
!        if (dime==2) then
!                if (vprel(2)==3) then       ! Calcul en CP (a=0)
!                        a = 0.d0
!                else
!                        a = 1.d0            ! Calcul en DP (a=1)
!                end if
!                trvsig = (1.+a*cnu)*(vsig(1)+vsig(2))
!                vun = (/ 1.d0, 1.d0, 0.d0 /)
!                vsigd = vsig - 1.d0/3.d0*trvsig*vun
!                vmis = sqrt(3.*(vsigd(1)**2. + vsigd(2)**2. + vsigd(3)**2. + vsigd(1)*vsigd(2)))
!        elseif (dime==3) then
!                trvsig = vsig(1)+vsig(2)+vsig(3)
!                vun = (/ 1.d0, 1.d0, 1.d0, 0.d0, 0.d0, 0.d0 /)
!                vsigd = vsig - 1.d0/3.d0*trvsig*vun
!                vmis = sqrt(3.d0*(vsigd(1)**2. + vsigd(2)**2. + vsigd(3)**2. + vsigd(1)*vsigd(2) + vsig(2)*vsig(3) + vsigd(1)*vsigd(3)))
!        else
!                print*,'Loi de Von Mises non implantee en 1D'
!                stop
!        end if
!
!        deallocate (vun, vsigd)
!
!        !----- Critere de Von Mises
!        vcrit = vmis - RT
!
!        !----- Calcul des derivees du critere
!        if (ideriv .eqv. .true.) then
!                !vdfdsig = 0.d0
!                vdfdsigloc = 0.d0
!                if (vmis /= 0.d0) then
!                        !vdfdsig(1,:) = 3.d0/(2.d0*vmis)*vsigd
!                        vdfdsigloc(1,:) = 3.d0/(2.d0*vmis)*vsigd         ! a corriger
!                end if
!           end if
!
!        if (ideriv .eqv. .true.)  vdgdsigloc = transpose(vdfdsigloc) !vdgdsig = transpose(vdfdsig)

!********************************************************!

!        case(13)     !----- Loi de Coulomb

!        !----- Recuperation des parametres de la loi
!        id = 6
!        if (dime == 3) id = id-1
!        C  = vprel(id+2)                      ! Cohesion
!        phi = pi*vprel(id+3)/180.             ! Angle de frottement
!        psi = pi*vprel(id+4)/180.             ! Angle de dilatance
!        KP = (1.+sin(phi))/(1.-sin(phi))      ! Poussee
!        KPS = (1.+sin(psi))/(1.-sin(psi))     ! Poussee
!        RP = 2.*C*cos(phi)/(1.-sin(phi))      ! Limite elastique en compression
!
!        !----- Critere de Coulomb
!        !vcrit = KP*vspmax - vspmin - RP
!        vcrit = KP*vsp(1) - vsp(3) - RP
!
!        !----- Calcul des derivees du critere
!        if (ideriv .eqv. .true.) then
!                if (dime == 2) then
!                           !vdfdsig(1,:) = (/ KPS*nf(1)**2.-nf(2)**2., KPS*nf(2)**2.-nf(1)**2., 2.*(KPS+1.d0)*nf(1)*nf(2)/)
!                        vdfdsigloc(1,:) = (/ KPS*nf(1)**2.-nf(2)**2., KPS*nf(2)**2.-nf(1)**2., 2.*(KPS+1.d0)*nf(1)*nf(2)/)     ! a corriger
!                elseif (dime == 3) then
!                        !vdfdsig(1,:) = (/ nf(1)**2., nf(2)**2., nf(3)**2., 2.*nf(1)*nf(2), 2.*nf(2)*nf(3), 2.*nf(3)*nf(2)/)
!                        vdfdsigloc(1,:) = (/ nf(1)**2., nf(2)**2., nf(3)**2., 2.*nf(1)*nf(2), 2.*nf(2)*nf(3), 2.*nf(3)*nf(2)/)   ! a corriger
!                end if
!        end if
!
!        if (ideriv .eqv. .true.)  vdgdsigloc = transpose(vdfdsigloc) ! vdgdsig = transpose(vdfdsig)

!********************************************************!

!         case(14)     !----- Loi de Coulomb "orientee"

!        !----- Recuperation des parametres de la loi
!        id = 6
!        if (dime == 3) id = id-1
!        C  = vprel(id+2)                      ! Cohesion
!        phi = pi*vprel(id+3)/180.             ! Angle de frottement
!        psi = pi*vprel(id+4)/180.             ! Angle de dilatance
!
!        !----- Calcul de la (des) direction(s) tangentielle(s) a la fissure
!        if (dime == 2) then
!                t1(:) = (/ -nf(2), nf(1)/)
!                t1(:) = t1(:)/norme(t1)
!        elseif (dime == 3) then
!                t1(:) = (/ -nf(2), nf(1), 0.d0/)
!                t2(:) = (/ -nf(3)*nf(1), -nf(3)*nf(2), nf(1)*nf(1)+nf(2)*nf(2) /)
!                t2(:) = t2(:)/norme(t2)
!        end if
!
!        !----- Calcul du vecteur des contraintes
!        T = matmul(VVsig,nf)
!
!        !----- Decomposition de T en deux vecteurs :
!        !le premier dirige selon la normale a la fissure (sigma*nf) et le deuxieme appartenant au plan de la fissure (Tt)
!        sigma =  dot_product(T,nf)
!        Tt = T - (sigma*nf)
!
!        !----- Decomposition de T selon la (les) direction(s) tangentielle(s) a la fissure
!        Tt1 = dot_product(Tt,t1)
!        if (dime == 3) then
!                Tt2 = dot_product(Tt,t1)
!        end if
!
!        !----- Definition de tau
!        tau = norme(Tt)
!
!        !----- Critere de Mohr-Coulomb oriente
!        vcrit = abs(tau) + sigma*tan(phi)
!
!        !----- Calcul des derivees du critere dans le repere principal
!        if (ideriv .eqv. .true.) then
!                if (tau /= 0.d0) then
!                        if (dime == 2) then  ! a verifier
!                                vdfdsigloc(1,:) = (/ tan(phi), 0.d0, Tt1/abs(tau) /)
!                                vdgdsigloctr(1,:) = (/ tan(psi), 0.d0, Tt1/abs(tau) /)
!                        elseif (dime == 3) then   ! a verifier
!                                vdfdsigloc(1,:) = (/ tan(phi), 0.d0, 0.d0, 0.d0, Tt2/abs(tau), Tt1/abs(tau) /)
!                                vdgdsigloctr(1,:) = (/ tan(psi), 0.d0, 0.d0, 0.d0, Tt2/abs(tau), Tt1/abs(tau) /)
!                        end if
!                else
!                        if (dime == 2) then   ! a verifier
!                                vdfdsigloc(1,:) = (/ tan(phi), 0.d0, 0.d0 /)
!                                vdgdsigloctr(1,:) = (/ tan(psi), 0.d0, 0.d0 /)
!                        elseif (dime == 3) then   ! a verifier
!                                vdfdsigloc(1,:) = (/ tan(phi), 0.d0, 0.d0, 0.d0, 0.d0, 0.d0 /)
!                                vdgdsigloctr(1,:) = (/ tan(psi), 0.d0, 0.d0, 0.d0, 0.d0, 0.d0 /)
!                        end if
!                end if
!                vdgdsigloc = transpose(vdgdsigloctr)
!
!           end if


 end select


end subroutine fiss_crit

!********************************************************!

!---------------------------------------------------------------------------
!> @author JL Tailhan
!> @brief Changement de repère des contraintes (global/principal)
!
!> @details
!> ### DESCRIPTION:
!> Definition de la matrice de changement de repere :
!> du repere global au repere principal
!>
!>--------------------------------------------------------------------------
function fiss_changement_repere_old(V,n,ind) result (M)

        use variables, only : dime
        use math, only : norme

        implicit none

        real*8, dimension(dime,dime), intent(in) :: V
        integer, intent(in) :: n, ind

        real*8, dimension(n,n) :: M
        real*8 :: n1, n2, n3, a, b, c, d, e, f, g, h, i

!********************************************************!
        n1 = norme(V(:,1))
        n2 = norme(v(:,2))
        a = V(1,1)/n1
        b = V(2,1)/n1

        d = V(1,2)/n2
        e = V(2,2)/n2

        if (dime == 3) then
            n3 = norme(V(:,3))
            c = V(3,1)/n1
            f = V(3,2)/n2
            g = V(1,3)/n3
            h = V(2,3)/n3
            i = V(3,3)/n3
        endif

        !----- Calcul de la matrice M=P : passage repere global => repere principal
        if (ind==1) then

            if (dime == 2) then
                M(1,:) = (/ (a**2), (b**2), (2.d0*a*b) /)
                M(2,:) = (/ (d**2), (e**2), (2.d0*d*e) /)
                M(3,:) = (/ (a*d) , (b*e) , (b*d+a*e)  /)
            elseif (dime == 3) then
                M(1,:) = (/ (a**2), (b**2), (c**2), (2.d0*a*b), (2.d0*b*c), (2.d0*a*c) /)
                M(2,:) = (/ (d**2), (e**2), (f**2), (2.d0*d*e), (2.d0*e*f), (2.d0*f*d) /)
                M(3,:) = (/ (g**2), (h**2), (i**2), (2.d0*g*h), (2.d0*h*i), (2.d0*i*g) /)
                M(4,:) = (/ (a*d) , (b*e) , (c*f) , (b*d+a*e) , (c*e+b*f) , (c*d+a*f)  /)
                M(5,:) = (/ (d*g) , (e*h) , (f*i) , (e*g+d*h) , (f*h+e*i) , (f*g+d*i)  /)
                M(6,:) = (/ (a*g) , (b*h) , (c*i) , (b*g+a*h) , (c*h+b*i) , (c*g+a*i)  /)
            endif

        !----- Calcul de la matrice M=Q : passage repere principal => repere global
        elseif (ind==2) then

            if (dime == 2) then
                M(1,:) = (/ (a**2), (d**2), (2.d0*a*d) /)
                M(2,:) = (/ (b**2), (e**2), (2.d0*b*e) /)
                M(3,:) = (/ (a*b) , (d*e) , (b*d+a*e)  /)
            elseif (dime == 3) then
                M(1,:) = (/ (a**2), (d**2), (g**2), (2.d0*a*d), (2.d0*d*g), (2.d0*g*a) /)
                M(2,:) = (/ (b**2), (e**2), (h**2), (2.d0*b*e), (2.d0*e*h), (2.d0*h*b) /)
                M(3,:) = (/ (c**2), (f**2), (i**2), (2.d0*c*f), (2.d0*f*i), (2.d0*i*c) /)
                M(4,:) = (/ (a*b) , (d*e) , (g*h),  (b*d+a*e) , (e*g+d*h) , (g*b+a*h)  /)
                M(5,:) = (/ (b*c) , (e*f) , (h*i),  (c*e+b*f) , (h*f+e*i) , (h*c+b*i)  /)
                M(6,:) = (/ (a*c) , (d*f) , (i*g),  (c*d+a*f) , (g*f+d*i) , (g*c+a*i)  /)
            endif

        endif

end function fiss_changement_repere_old


!---------------------------------------------------------------------------
!> @author JL Tailhan
!> @brief Changement de repère des contraintes (global/principal)
!
!> @details
!> ### DESCRIPTION:
!> Definition de la matrice de changement de repere :
!> du repere global au repere principal
!>
!>--------------------------------------------------------------------------
!   ind = 1 matrice de changement de rep global -> local
!   ind = 2 matrice de changement inverse (transposee de la prec.)
!   iopt = 1 travail sur les contraintes
!   iopt = 2 travail sur les deformations
!---------------------------------------------------------------------------

function fiss_changement_repere(V,n,iopt,typind) result (M)

        use variables, only : dime
        use math, only : norme

        implicit none

        !----- Input variables
        real*8, dimension(dime,dime), intent(in) :: V
        integer, intent(in) :: n, iopt
        integer, intent(in), optional :: typind

        !----- Output result
        real*8, dimension(n,n) :: M

        !----- Local variables
        real*8 :: n1, n2, n3, a, b, c, d, e, f, g, h, i
        integer:: ind

!********************************************************!
        n1 = norme(V(:,1))
        n2 = norme(v(:,2))
        
        a = V(1,1)/n1
        b = V(2,1)/n1

        d = V(1,2)/n2
        e = V(2,2)/n2

        if (dime == 3) then
            n3 = norme(V(:,3))
            c = V(3,1)/n1
            f = V(3,2)/n2
            g = V(1,3)/n3
            h = V(2,3)/n3
            i = V(3,3)/n3
        endif
!********************************************************!

        ind = 1
        if (present(typind)) then
            ind = typind
        endif

        !----- Pour les contraintes
        if (iopt==1) then

            if (ind==1) then ! Matrice de passage : Sig_p = Ps Sig
                if (dime == 2) then
                    M(1,:) = (/ (a**2), (b**2), (2.d0*a*b) /)
                    M(2,:) = (/ (d**2), (e**2), (2.d0*d*e) /)
                    M(3,:) = (/ (a*d) , (b*e) , (b*d+a*e)  /)
                elseif (dime == 3) then
                    M(1,:) = (/ (a**2), (b**2), (c**2), (2.d0*a*b), (2.d0*b*c), (2.d0*a*c) /)
                    M(2,:) = (/ (d**2), (e**2), (f**2), (2.d0*d*e), (2.d0*e*f), (2.d0*f*d) /)
                    M(3,:) = (/ (g**2), (h**2), (i**2), (2.d0*g*h), (2.d0*h*i), (2.d0*i*g) /)
                    M(4,:) = (/ (a*d) , (b*e) , (c*f) , (b*d+a*e) , (c*e+b*f) , (c*d+a*f)  /)
                    M(5,:) = (/ (d*g) , (e*h) , (f*i) , (e*g+d*h) , (f*h+e*i) , (f*g+d*i)  /)
                    M(6,:) = (/ (a*g) , (b*h) , (c*i) , (b*g+a*h) , (c*h+b*i) , (c*g+a*i)  /)
                endif
            elseif (ind==2) then ! Matrice de passage inverse : Sig = Transpose(Ps) Sig_p
                if (dime == 2) then
                    M(1,:) = (/ (a**2)     , (d**2)     , (2.d0*a*d)     /)
                    M(2,:) = (/ (b**2)     , (e**2)     , (2.d0*b*e)     /)
                    M(3,:) = (/ (a*b) , (d*e) , (b*d+a*e) /)
                elseif (dime == 3) then
                    M(1,:) = (/ (a**2), (d**2), (g**2), (2.d0*a*d), (2.d0*d*g), (2.d0*a*g)/)
                    M(2,:) = (/ (b**2), (e**2), (h**2), (2.d0*b*e), (2.d0*e*h), (2.d0*b*h)/)
                    M(3,:) = (/ (c**2), (f**2), (i**2), (2.d0*c*f), (2.d0*f*i), (2.d0*c*i)/)
                    M(4,:) = (/ (a*b) , (d*e) , (g*h) , (b*d+a*e) , (e*g+d*h) , (b*g+a*h) /)
                    M(5,:) = (/ (b*c) , (e*f) , (h*i) , (c*e+b*f) , (f*h+e*i) , (c*h+b*i) /)
                    M(6,:) = (/ (a*c) , (f*d) , (i*g) , (c*d+a*f) , (f*g+d*i) , (c*g+a*i) /)
                endif
            endif

        !----- Pour les deformations
        elseif (iopt==2) then

            if (ind==1) then ! Matrice de passage Eps_p = Pe Eps : 
                if (dime == 2) then
                    M(1,:) = (/ (a**2)     , (b**2)     , (a*b)     /)
                    M(2,:) = (/ (d**2)     , (e**2)     , (d*e)     /)
                    M(3,:) = (/ (2.d0*a*d) , (2.d0*b*e) , (b*d+a*e) /)
                elseif (dime == 3) then
                    M(1,:) = (/ (a**2)     , (b**2)     , (c**2)     , (a*b)     , (b*c)     , (a*c)     /)
                    M(2,:) = (/ (d**2)     , (e**2)     , (f**2)     , (d*e)     , (e*f)     , (f*d)     /)
                    M(3,:) = (/ (g**2)     , (h**2)     , (i**2)     , (g*h)     , (h*i)     , (i*g)     /)
                    M(4,:) = (/ (2.d0*a*d) , (2.d0*b*e) , (2.d0*c*f) , (b*d+a*e) , (c*e+b*f) , (c*d+a*f) /)
                    M(5,:) = (/ (2.d0*d*g) , (2.d0*e*h) , (2.d0*f*i) , (e*g+d*h) , (f*h+e*i) , (f*g+d*i) /)
                    M(6,:) = (/ (2.d0*a*g) , (2.d0*b*h) , (2.d0*c*i) , (b*g+a*h) , (c*h+b*i) , (c*g+a*i) /)
                end if
            elseif (ind==2) then ! Matrice de passage inverse : Eps = Transpose(Pe) Eps_p
                if (dime == 2) then
                    M(1,:) = (/ (a**2) , (d**2) , (a*d) /)
                    M(2,:) = (/ (b**2) , (e**2) , (b*e) /)
                    M(3,:) = (/ (2.d0*a*b)  , (2.d0*d*e)  , (b*d+a*e)  /)
                elseif (dime == 3) then
                    M(1,:) = (/ (a**2)    , (d**2)    , (g**2)    , (a*d)    , (d*g)    , (a*g)    /)
                    M(2,:) = (/ (b**2)    , (e**2)    , (h**2)    , (b*e)    , (e*h)    , (b*h)    /)
                    M(3,:) = (/ (c**2)    , (f**2)    , (i**2)    , (c*f)    , (f*i)    , (c*i)    /)
                    M(4,:) = (/ (2.d0*a*b), (2.d0*d*e), (2.d0*g*h), (b*d+a*e), (e*g+d*h), (b*g+a*h)/)
                    M(5,:) = (/ (2.d0*b*c), (2.d0*e*f), (2.d0*h*i), (c*e+b*f), (f*h+e*i), (c*h+b*i)/)
                    M(6,:) = (/ (2.d0*a*c), (2.d0*f*d), (2.d0*i*g), (c*d+a*f), (f*g+d*i), (c*g+a*i)/)
                endif
            endif

        end if

end function fiss_changement_repere



!---------------------------------------------------------------------------
!> @author JL Tailhan
!> @brief Pilotage du calcul sur l'avancement de la fissuration
!
!> @details
!> ### DESCRIPTION:
!> Pilotage du calcul sur l'avancement de la fissuration pourle modèle de
!> fissuration semi-explicite
!>
!>--------------------------------------------------------------------------
        subroutine fiss_pilot(imetpilo,vsol,vduI,vduII,dlam,alpham,MOTm,elefiss,iedg)

        !********************************************************!
        !  Pilotage du calcul : gestion du facteur de chargement !
        !********************************************************!

            use variables, only : dime, nelt, idmax, ktypel, nomtype, infele, &
                        & vprelg, kprop, ietat, iedngm, fiss, pi,     &
                        & mrtrav, vrtrav, kloce, fc
            use initialisation, only : init_vec, init_mat
            use lib_elem, only : elem_B, elem_kloce2, elem_hooke
            use math
            use utilitaire, only : princ

            implicit none

            ! Variables IN
            real*8, dimension(:), intent(in) :: vduII, vsol
            integer, intent(in) :: imetpilo
            integer, intent(in), optional :: iedg

            ! Variables IN-OUT
            real*8, dimension(:), intent(in) :: vduI
            real*8, intent(in) :: dlam
            real*8, intent(out) :: alpham
            character(len=8), intent(inout) :: MOTm
            integer, intent(inout) :: elefiss

            !real*8, dimension(dime) :: normfissm

            ! Variables locales
            real*8, dimension(:), allocatable :: vdle, vdl0, vsi0, vdsi, vsi, vep, &
                                                & vsi0loc, vdsiloc, veploc, alpg, vdlI, vdlII
            real*8, dimension(:,:), allocatable :: vn, vb, vh, ksig, Ps, Pe
            
            real*8, dimension(dime, dime) :: rotmatm
            real*8 :: vprel(idmax), alph(nelt), vbid(dime), vbid2(1,dime*dime), V(dime,dime)
            real*8 :: a, b, c, RT, RC, &
                        & vcrit1, vcrit2, detj, &
                        & s1, s3, alpt, alpc
            character(len=5) :: typel
            character(len=8) :: MOT
            character(len=8), dimension(nelt) :: rupt
            integer :: ie, ipg, npg, id, iloi, ndle, imin, imax
            integer :: binf, bsup, iedang
            
        !----- Switch en fonction de l'option de calcul
        !      iedang => fiss_pilot applique a iedang seulement (pour obtenir le repere local)
            binf = 1
            bsup = nelt
            if (present(iedg)) then
                iedang = iedg
                binf = iedang
                bsup = iedang
            else
                iedang = 0
            endif

        !----- Switch en fonction de la methode de pilotage

            select case(imetpilo)

                case(1)
                    !----- Pilotage sur element le plus dangereux par recalcul de vduI

                    !- Pour le cas beton fissurant
                    if (fiss==1) then

                        !- Initialisations
                        alph = 1.d0
                        iedngm = 0

                        !- Boucle sur les elements
                        do ie = binf, bsup

                            MOT = '  '
                            alpt = 1.d0
                            alpc = 1.d0
                            typel = nomtype(ktypel(ie))

                            !----- Pour les elements massifs vierges
                            if ((typel == 'MBT3' .or. typel == 'MBQ4' .or. typel == 'MTT4' .or. typel == 'MBT6')  &
                                & .and. (ietat(ie)==0)) then

                                !----- Proprietes elementaires
                                vprel = vprelg(kprop(ie),1:idmax)

                                !----- Recuperation des parametres de la loi
                                iloi = int(vprel(1))
                                
                                !if (iloi==1) return

                                !-----  Recuperation des informations sur les elements
                                allocate(ksig(size(infele(ktypel(ie))%Q,1), size(infele(ktypel(ie))%Q,2)))
                                ksig = infele(ktypel(ie))%Q
                                npg = size(ksig,1)

                                !----- Initialisations relatives aux elements
                                call init_vec(vrtrav,npg*dime*dime)
                                allocate(alpg(npg))
                                alpg = 1.d0

                                !- Detection de la premiere fissuration (ouverture)
                                ipg = 1

                                !----- Vecteur kloce de localisation dans vecteur global
                                call elem_kloce2( ie, ndle)
                                allocate(vdl0(ndle));  allocate(vdle(ndle))
                                allocate(vdlI(ndle));  allocate(vdlII(ndle))

                                !----- Deplacements
                                vdl0  = vsol(kloce(1:ndle))
                                vdlI  = vduI(kloce(1:ndle))
                                vdlII = vduII(kloce(1:ndle))
                                vdle  = vdlI + dlam*vdlII

                                !----- Matrice d'elasticite vh, en fonction du type d'element
                                call elem_hooke(vh,nomtype(ktypel(ie)),vprel)

                                !----- Calcul des contraintes au centre de gravite de l'element
                                call elem_B(vn,vb,detj,ksig(ipg,:),ie)

                                ! dans le repere global
                                allocate(vsi0(size(vb,1)))
                                allocate(vdsi(size(vb,1)))
                                vsi0  = matmul(vh,matmul(vb,vdl0))
                                vdsi  = matmul(vh,matmul(vb,vdle))

                                !----- Calcul et rangement par ordre decroissant des contraintes principales
                                allocate(vsi(size(vsi0,1)))
                                allocate(vep(size(vdl0,1)))
                                vsi = vsi0 + vdsi
                                vep = matmul(vb, vdl0 + vdle)
                                call princ(vsi,vbid,V)
                                imax = 1
                                imin = 2
                                if (dime==3) imin = 3

                                !----- Recuperation de la matrice de changement de repere global --> local principal
                                call init_mat(Ps,size(vh,1),size(vh,2))
                                call init_mat(Pe,size(vh,1),size(vh,2))                                
                                Ps = fiss_changement_repere(V,size(vb,1),1,1) !--- changement de repere pour les contraintes
                                Pe = fiss_changement_repere(V,size(vb,1),2,1) !--- changement de repere pour les deformations

                                ! Dans le repere principal...
                                allocate(vdsiloc(size(vh,1)))
                                allocate(vsi0loc(size(vh,1)))
                                allocate(veploc(size(vh,1)))
                                vdsiloc = matmul(Ps,vdsi)
                                vsi0loc = matmul(Ps,vsi0)
                                veploc  = matmul(Pe,vep)

                                deallocate(vdl0,vdlI,vdlII,vdle,vep,vh,vb,vn)
                                deallocate(vsi0,vdsi,vsi,Ps,Pe)
                                
                                if (iloi==5) then

                                !----- Recuperation des parametres de la loi
                                    id = 6
                                    if (dime == 3) id = id-1

                                    RT = 1.001d0 * vprel(id)     ! Contrainte limite en traction pure
                                    RC = 1.001d0 * vprel(id+1)   ! Contrainte limite en cisaillement
                                    RC = 10.d0*RT

                                    !----- Recuperation des contraintes principales
                                    s1 = vsi0loc(imax)+vdsiloc(imax)
                                    s3 = vsi0loc(imin)+vdsiloc(imin)

                                    !----- Critere de rupture en ouverture pure
                                    vcrit1 = s1 - RT
                                   
                                    !----- Critere de rupture en cisaillement
                                    vcrit2 = s1 - s3 - RC

                                    if (ie==iedang) then
                                        if (vcrit1>vcrit2) then 
                                            vcrit1 = 1.d0
                                        else
                                            vcrit2 = 1.d0
                                        endif
                                    endif

                                    if (vcrit1 > 0.d0 .or. vcrit2 > 0.d0) then
                                        if (vcrit1 > 0.d0 .and. vcrit2 > 0.d0) then
                                            !- Critere en traction pure
                                            alpt = (RT-vsi0loc(imax))/vdsiloc(imax)
                                            !- Critere en cisaillement
                                            alpc = (RC-(vsi0loc(imax)-vsi0loc(imin)))/(vdsiloc(imax)-vdsiloc(imin))
                                            
                                            alpg(ipg) = min(alpt, alpc)
                                            if (alpc < alpt) then
                                                MOT = 'cisaille'
                                            else
                                                MOT = 'ouvert'
                                            end if
                                        elseif (vcrit1 > 0.d0) then
                                            MOT = 'ouvert'
                                            alpg(ipg) = (RT-vsi0loc(imax))/vdsiloc(imax)
                                        elseif (vcrit2 > 0.d0) then
                                            MOT = 'cisaille'
                                            alpg(ipg) = (RC-(vsi0loc(imax)-vsi0loc(imin)))/(vdsiloc(imax)-vdsiloc(imin))
                                        end if
                                    end if
                                    
                                    if (MOT == 'cisaille') then
                                        if (veploc(size(veploc)) > 0.d0) then
                                            if (dime == 2) then
                                                rotmatm = rotmat(-pi/4.d0, dime)
                                            else
                                                rotmatm = rotmat(-pi/4.d0, dime, V(:, 2))
                                            end if
                                        elseif (veploc(size(veploc)) < 0.d0) then
                                            if (dime == 2) then
                                                rotmatm = rotmat(+pi/4.d0, dime)
                                            else
                                                rotmatm = rotmat(+pi/4.d0, dime, V(:, 2))
                                            end if
                                        else
                                            if (veploc(1) < veploc(dime)) then
                                                if (dime == 2) then
                                                    rotmatm = rotmat(pi/4.d0, dime)
                                                else
                                                    rotmatm = rotmat(pi/4.d0, dime, V(:, 2))
                                                end if
                                            else
                                                if (dime == 2) then
                                                    rotmatm = rotmat(-pi/4.d0, dime)
                                                else
                                                    rotmatm = rotmat(-pi/4.d0, dime, V(:, 2))
                                                end if
                                            end if
                                        end if
                                        !normfissm = matmul(rotmatm, V(:, 1))
                                        if (dime==2) then
                                            V(:, 1) = matmul(rotmatm, V(:, 1))
                                            V(:, 1) = V(:, 1) / norme(V(:, 1))
                                            V(:, 2) = matmul(rotmatm, V(:, 2))
                                            V(:, 2) = V(:, 2) / norme(V(:, 2))
                                        elseif (dime==3) then
                                            V(:, 1) = matmul(rotmatm, V(:, 1))
                                            V(:, 1) = V(:, 1) / norme(V(:, 1))
                                            V(:, 2) = V(:, 2) / norme(V(:, 2))
                                            V(:, 3) = matmul(rotmatm, V(:, 3))
                                            V(:, 3) = V(:, 3) / norme(V(:, 3))
                                        endif
                                    elseif (MOT == 'ouvert') then
                                        !normfissm = V(:, 1)
                                        V(:, 1) = V(:, 1) / norme(V(:, 1))
                                        V(:, 2) = V(:, 2) / norme(V(:, 2))
                                        if (dime==3) V(:, 3) = V(:, 3) / norme(V(:, 3))
                                    end if
                                ! Fin du cas loi 5

                                elseif (iloi==12) then

                                !----- Recuperation des parametres de la loi
                                    id = 6
                                    if (dime == 3) id = id-1

                                    RT = 1.001d0 * vprel(id)                        ! Contrainte limite en traction pure
                                    RC = 1.001d0 * fc                               ! Contrainte limite en cisaillement

                                    !----- Recuperation des contraintes principales
                                    s1 = maxval(vsi0loc+vdsiloc)
                                    s3 = minval(vsi0loc+vdsiloc)

                                    !----- Critere de rupture en ouverture pure
                                    vcrit1 = s1 - RT
                                   
                                    !----- Critere de rupture en cisaillement
                                    vcrit2 = s1 - s3 - RC

                                    if (vcrit1 > 0.d0 .or. vcrit2 > 0.d0) then
                                        if (vcrit1 > 0.d0 .and. vcrit2 > 0.d0) then
                                            !- Critere en traction pure
                                            alpt = (RT-maxval(vsi0loc))/maxval(vdsiloc)
                                            !- Critere en cisaillement
                                            alpc = (RC-(maxval(vsi0loc)-minval(vsi0loc)))/(maxval(vdsiloc)-minval(vdsiloc))
                                            
                                            alpg(ipg) = min(alpt, alpc)
                                            if (alpc < alpt) then
                                                MOT = 'cisaille'
                                            else
                                                MOT = 'ouvert'
                                            end if
                                        elseif (vcrit1 > 0.d0) then
                                            MOT = 'ouvert'
                                            alpg(ipg) = (RT-maxval(vsi0loc))/maxval(vdsiloc)
                                        elseif (vcrit2 > 0.d0) then
                                            MOT = 'cisaille'
                                            alpg(ipg) = (RC-(maxval(vsi0loc)-minval(vsi0loc)))/(maxval(vdsiloc)-minval(vdsiloc))
                                        end if
                                    end if
                                    
                                    if (MOT == 'cisaille') then
                                        if (veploc(size(veploc)) > 0.d0) then
                                            if (dime == 2) then
                                                rotmatm = rotmat(-pi/4.d0, dime)
                                            else
                                                rotmatm = rotmat(-pi/4.d0, dime, V(:, 2))
                                            end if
                                        elseif (veploc(size(veploc)) < 0.d0) then
                                            if (dime == 2) then
                                                rotmatm = rotmat(+pi/4.d0, dime)
                                            else
                                                rotmatm = rotmat(+pi/4.d0, dime, V(:, 2))
                                            end if
                                        else
                                            if (veploc(1) < veploc(dime)) then
                                                if (dime == 2) then
                                                    rotmatm = rotmat(pi/4.d0, dime)
                                                else
                                                    rotmatm = rotmat(pi/4.d0, dime, V(:, 2))
                                                end if
                                            else
                                                if (dime == 2) then
                                                    rotmatm = rotmat(-pi/4.d0, dime)
                                                else
                                                    rotmatm = rotmat(-pi/4.d0, dime, V(:, 2))
                                                end if
                                            end if
                                        end if
                                        !normfissm = matmul(rotmatm, V(:, 1))
                                    elseif (MOT == 'ouvert') then
                                        !normfissm = V(:, 1)
                                    else
                                    end if
                                ! Fin du cas loi 12

                                elseif (iloi==15) then

                                    !----- Recuperation des parametres de la loi
                                    id = 6
                                    if (dime == 3) id = id-1

                                    RT = 1.0001d0 * vprel(id)                        ! Contrainte limite

                                    !----- Recuperation de la contrainte principale 1
                                    s1 = vsi0loc(1)+vdsiloc(1)

                                    !----- Critere de rupture en ouverture pure
                                    vcrit1 = s1 - RT

                                    if (vcrit1 > 0.d0) then

                                        !- Critere en traction pure
                                        MOT = 'ouvert'
                                        alpg(ipg) = (RT-vsi0loc(1))/vdsiloc(1)
                                        !normfissm = V(:, 1)
                                    else
                                        !- Critere en cisaillement
                                        ! MOT = 'cisaille'
                                    end if
                                    ! Fin du cas loi 15

                                elseif ((iloi==1).or.(iloi==11).or.(iloi==27)) then

                                else
                                    stop 'FIDES_fiss_pilot : loi non encore programmee'
                                end if

                                deallocate(vsi0loc,vdsiloc,veploc)

                                !----- On conserve le repere local sur le point de gauss considere
                                vbid2=reshape(V,(/1,dime*dime/))
                                vrtrav((ipg-1)*dime*dime+1:ipg*dime*dime)=vbid2(1,:)

                                !----- On conserve aussi
                                rupt(ie) = MOT
                                alph(ie) = minval(alpg)

                                !----- On conserve le repere local du point de Gauss le plus dangereux
                                ipg = find_num(alpg,minval(alpg))
                                mrtrav(ie,:)=vrtrav((ipg-1)*dime*dime+1:(ipg*dime*dime))

                                deallocate(vrtrav)
                                deallocate(alpg,ksig)
                            end if

                        end do
                        ! Fin de la boucle sur les elements

                        alpham = minval(alph)

                        if (alpham < 1.d0) then
                            !----- On stocke le numero de l'element le plus dangereux
                            iedngm = find_num(alph,alpham)
                            elefiss = count(alph<=0)
                            MOTm = rupt(iedngm)
                        else
                            alpham = 1.d0
                        end if
                    end if

                case(2)
                    !----- Pilotage sur l'increment de contrainte avec recalcul de "dlam"

                case default
                    stop 'FIDES_fiss_pilot : cas non encore programme'
            end select

        end subroutine fiss_pilot


!---------------------------------------------------------------------------
!> @author JL Tailhan
!> @brief Routine de gestion du changement d'état de l'élement
!> (modèle semi-explicite)
!
!> @details
!> ### DESCRIPTION:
!> Routine de gestion du changement d'état de l'élement
!> dans le cas du modèle de fissuration semi-explicite
!>
!>--------------------------------------------------------------------------
subroutine fiss_change_etat()


        use variables, only : nelt, vsol, vdu,ietat, detoscill, deboscill, detoscillp, comptoscill, &
                  & vprelg, kprop,iouver1, iouver2, histetat1, histetat2, fiss, kloce, iter
        use initialisation, only : init_vec
        use lib_elem, only : elem_kloce2

        implicit none

        real*8, dimension(:), allocatable :: vdle,vdle0,vddle
        real*8 :: ouv, ouv0, freqoscill
        integer :: ie, ietat0, ietat1, ietat2, ndle
        integer :: limite ! nombre de retour a partir duquel on considere qu'il y a oscillation '

!********************************************************!

    limite=5

        if (fiss==1) then
            do ie = 1, nelt                       ! Boucle sur les elements

                freqoscill = 0.d0

                if (ietat(ie) /= 0) then      ! Element fissure

                    ! Ouverture
                    ouv = fiss_ouv(ie)

                    if ((ietat(ie)==1) .or. (ietat(ie)==2)) then

                        !iouver2(ie) = iouver1(ie)
                        !iouver1(ie) = ouv

                        ietat2 = histetat2(ie)
                        ietat1 = histetat1(ie)
                        ietat0 = ietat(ie)

                        ! La detection des oscillations se fait sur l'estimation
                        ! d'une frequence d'oscillation estimee toutes les
                        ! n(=limite) iterations
                        if ((ietat0/=ietat1)) then
                            deboscill(ie)=deboscill(ie)+1
                            detoscill(ie)=detoscill(ie)+1
                            if (deboscill(ie) == 1) then
                                if (detoscill(ie)/=1) then
                                    freqoscill = real(detoscill(ie))-real(detoscillp(ie))
                                    freqoscill = freqoscill / real(limite)
                                endif
                                detoscillp(ie) = detoscill(ie)
                            endif
                            if (deboscill(ie) >= limite-1) then
                                deboscill(ie) = 0
                            endif
                        end if
      !if (detoscill(ie)/=0) print*,ie,detoscill(ie),deboscill(ie),freqoscill

                        histetat1(ie) = ietat0
                        histetat2(ie) = ietat1
                    end if

                    if ((freqoscill>0.5d0).and.(ietat(ie)/=-2)) then
                        ietat(ie) = -2
                        print*, 'Oscillation de l''element',ie,' iloi',int(vprelg(kprop(ie),1))
                        detoscill(ie) = 0
                        !detoscill(ie) = detoscill(ie)+1
                        !comptoscill = comptoscill+1

                        !----- Vecteur kloce de localisation pour assemblage
                        !call elem_kloce2(ie,ndle)
                        !allocate(vdle(ndle))
                        !allocate(vdle0(ndle))
                        !allocate(vddle(ndle))

                        !----- Modification du deplacement total
                        !vdle  = vsol(kloce(1:ndle))
                        !vddle = vdu(kloce(1:ndle))
                        !vdle0 = vdle - vddle 

                        !vsol(kloce(1:ndle)) = vdle0
                        !ouv0 = fiss_ouv(ie)

                        !vsol(kloce(1:ndle)) = (abs(ouv)*vdle0 + abs(ouv0)*vdle)/(abs(ouv0)+abs(ouv))

                        !deallocate(vdle,vddle,vdle0)

                    end if
!if (iter>100) pause
                end if
            end do
        end if

end subroutine fiss_change_etat


!---------------------------------------------------------------------------
!> @author JL Tailhan
!> @brief Calcul de l'ouverture moyenne de fissure de l'élément ie
!
!> @details
!> ### DESCRIPTION:
!> Calcul de l'ouverture moyenne de fissure de l'élément ie.
!>
!> L'ouverture est calculée par le déplacement relatif des noeuds de
!> l'élément situés de part et d'autre de la normale à la fissure.
!> Le déplacement relatif normal moyen est retenu comme ouverture sur
!> l'élément.
!>
!>--------------------------------------------------------------------------
function fiss_ouv(ie) result (ouv)


        use variables, only : dime, inorm, vcor, kconec, infele, ktypel, vsol, infnoe
        use math, only : norme
        use initialisation, only : init_vec, init_mat

        implicit none

        integer, intent(in) :: ie

        real*8 :: XG, YG, ZG, AN, BN, CN, DN, DENO, dep1, dep2, ouv
        real*8, dimension(dime) :: nf, NN, PG
        real*8, dimension(:), allocatable :: depNorm
        real*8, dimension(infele(ktypel(ie))%nnel) :: dist
        !real*8, dimension(dime,infele(ktypel(ie))%nnel*dime) :: depNds
        real*8, dimension(dime,infele(ktypel(ie))%nnel) :: depNds
        real*8, dimension(dime,infele(ktypel(ie))%nnel) :: vcore
        integer :: i, nnzp, nnzn, noel, ino, ndln
        integer, dimension(infele(ktypel(ie))%nnel) :: vec
        integer, dimension(infele(ktypel(ie))%nnel,infele(ktypel(ie))%ndln) :: ddl

        logical, dimension(infele(ktypel(ie))%nnel) :: cotep

!*********************************************************

        !----- On recupere la normale a la fissure (vecteur ligne)
        nf = inorm(ie,:)
        NN = nf/norme(nf)

        !----- On recupere quelques information sur l'element
        noel = infele(ktypel(ie))%nnel
        ndln = infele(ktypel(ie))%ndln
        vec = kconec(ie,1:noel)

        !----- On recupere les coordonnees des noeuds de l'element
        do i = 1, noel
            vcore(:,i) = vcor(:,vec(i))
        end do

        if (dime==1) then
            stop 'FIDES_fiss_ouv : cas dime = 1 non prevu !'
        else
            !----- On calcule les deplacements des noeuds de l'element
            do ino = 1, noel
                ddl(ino,:) = infnoe(vec(ino))%dln(1:ndln)
            end do

            depNds(1,:) = vsol(ddl(:,1))
            depNds(2,:) = vsol(ddl(:,2))
            if (dime==3) then
                depNds(3,:) = vsol(ddl(:,3))
            end if
        end if

        !----- On calcule la projection normale a la fissure de ces deplacements
        allocate(depNorm(size(depNds,2)))
        depNorm = matmul(NN,depNds)

        !----- On calcule les distances des sommets au plan de fissuration pour
        !----- detecter la position de ces points par rapport au plan de
        !----- fissuration (le plan de fissuration est defini par PG et N)

        PG = sum(vcore,2)/size(vcore,2)

        XG = PG(1)
        YG = PG(2)
        AN = NN(1)
        BN = NN(2)

        if (dime==2) then
           DN = -(AN*XG + BN*YG)
           DENO = sqrt(AN**2. + BN**2.)
           dist=(matmul(NN,vcore) + DN)/DENO
        elseif (dime==3) then
           ZG = PG(3)
           CN = NN(3)
           DN=-(AN*XG + BN*YG + CN*ZG)
           DENO = sqrt(AN**2. + BN**2. + CN**2.)
           dist=(matmul(NN,vcore) + DN)/DENO
        else
           stop 'FIDES_fiss_ouv : cas dime = 1 non prevu !'
        end if

        !----- Calcul de l'ouverture moyenne de la fissure
        cotep = (dist>=0.d0)

        nnzp = count(dist>=0.d0)
        nnzn = count(dist<0.d0)

        dep1 = 0.d0
        dep2 = 0.d0

        do i = 1, noel
            if (cotep(i)) then
                dep1 = dep1 + (depNorm(i))/nnzp
            else
                dep2 = dep2 + (depNorm(i))/nnzn
            end if
        end do

        ouv = dep1 - dep2

        deallocate(depNorm)

end function fiss_ouv


!---------------------------------------------------------------------------
!> @author JL Tailhan
!> @brief Annulation des contraintes et vi si rupture (en ouverture)
!
!> @details
!> ### DESCRIPTION:
!> Procédure d'annumation des contraintes et variables internes lors de
!> la détection de la rupture de l'élément (en cas d'ouverture seulement)
!>
!>--------------------------------------------------------------------------
subroutine fiss_rupture(ie,vsig,vnle,wplael,vprel)


        use math, only : norme, find_vec
        use variables, only : calco, dime, ktypel, infele, fiss, ietat !, inorm, &
                            !& iedng !, icompfiss
        use initialisation, only : init_vec, init_mat
        use utilitaire, only : princ
        use lib_elem, only : elem_B

        implicit none

        real*8, dimension(:,:), intent(inout) :: vnle, vsig
        integer, intent(in) :: ie
        real*8, dimension(:), intent(in) :: vprel
        real*8, intent(in) :: wplael

        real*8, dimension(:), allocatable :: vpg
        real*8, dimension(:,:), allocatable :: ksin
        real*8, dimension(:,:), allocatable :: vn
        real*8, dimension(:,:), allocatable :: vb
        real*8, dimension(:,:), allocatable :: VVsigPG
        real*8, dimension(:), allocatable :: vsigPG
        integer, dimension(:), allocatable :: icol

        real*8, dimension(size(infele(ktypel(ie))%ksin,2)) :: ksiPG
        real*8, dimension(dime,dime) :: V
        real*8, dimension(dime) :: vsp
        real*8, dimension(dime) :: nf
        real*8 :: wplamax, detj, vspmax

        integer :: i, id, nnel, iloi, itype

!********************************************************!

        itype = ktypel(ie)

        ! Recuperation des informations sur les elements
        call init_mat(ksin, size(infele(itype)%ksin,1), size(infele(itype)%ksin,2))
        ksin = infele(itype)%ksin
        nnel = infele(itype)%nnel

        ! Recuperation des points d'integration
        call init_vec(vpg, size(infele(itype)%W,1))
        vpg = infele(itype)%W

        iloi=int(vprel(1))

        !if ((fiss==1) .and. (iloi==icompfiss)) then
        if (fiss==1) then
            if (iloi==12) then
                id=7
                if (dime==3) id=id-1
                !wplmax=vprel(id)   ! ... Modifie par JLT pour coherence entre lois 12 et 15
                wplamax = 0.d0
                if ((ietat(ie)==0).or.(ietat(ie)==3)) wplamax=vprel(id)
                ! ... Fin modif JLT
            else if (iloi==15) then
                return
            else
                !print *, 'Calcul impossible avec cette loi de comportement'
                !stop
            end if

            !if ((wplael>wplamax) .and. (ietat(ie)==0)) then ! ... modif JLT element plastifie: ietat=3
            if ((wplael>wplamax) .and. (ietat(ie)==3)) then
            !if ((ie==iedng) .or. (ietat(ie)/=0))then   ! (ietat(ie)==1)
            !if (ie==iedng) then

                ! Energie plastique maximale atteinte : l'element est fissure
                ietat(ie)=1

                !--- Calcul de la normale au plan de fissuration
                !
                ! On calcule les contraintes au centre de gravite par
                ! interpolation des contraintes calculees soit aux noeuds,
                ! soit aux points de Gauss
                call init_mat(VVsigPG, size(vsig,1), size(vsig,2))
                call init_vec(vsigPG, size(vsig,1))

                if (calco == 'NOEUD') then
                    ! Cas des contraintes deja calculees aux noeuds :
                    ! on se place au centre de gravite de l'element
                    ksiPG = sum(ksin,1)/nnel
                    call elem_B(vn,vb,detj,ksiPG,ie)
                    VVsigPG = matmul(vsig,transpose(vn))
                    vsigPG = sum(VVsigPG,2)
                    deallocate(vn,vb)
                else
                    ! Cas des contraintes calculees aux points de Gauss :
                    ! on calcule la moyenne (ponderee) des contraintes au
                    ! barycentre des points de Gauss
                    do i=1,size(vsig,2)
                        VVsigPG(:,i) = vpg(i)*vsig(:,i)/sum(vpg)
                    end do
                    vsigPG = sum(VVsigPG,2)
                end if

                call princ(vsigPG,vsp,V)

                ! Calcul de la normale a la fissure
                vspmax = maxval(vsp)                   ! Contrainte principale majeure
                call init_vec(icol,count(vsp==vspmax))
                icol = find_vec(vsp,vspmax)
                nf = sum(V(:,icol),2)
                nf = nf/norme(nf)
                !inorm(ie,:) = nf   ! TEMPORAIREMENT !!

                vnle=0.0d0*vnle    ! Attention !!! UNIQUEMENT POUR ELASTIQUE FRAGILE
                vsig=0.0d0*vsig

                deallocate(icol, VVsigPG, vsigPG)           ! AM
            end if
        end if

deallocate(ksin,vpg)

end subroutine fiss_rupture


!---------------------------------------------------------------------------
!> @author JL Tailhan
!> @brief Gestion de la raideur de l'élément
!
!> @details
!> ### DESCRIPTION:
!> Aiguillage pour le calcul de la raideur de l'élément en fonction
!> de son état (vierge, fissuré ouvert, fissuré refermé
!>
!>--------------------------------------------------------------------------
subroutine fiss_modul(ie,vh,vprel)

!********************************************************!
!                Gestion de la raideur                   !
!********************************************************!

        use variables, only : ietat, smear, depen!, icompfiss

        implicit none

        real*8, dimension(:), intent(in) :: vprel
        real*8, dimension(:,:), intent(inout) :: vh
        integer, intent(in) :: ie
        integer :: iloi

!********************************************************!

        iloi=int(vprel(1))

        !if (iloi==icompfiss) then
        if (iloi==12) then

        !***** Modification de la matrice pour element fissure
            if ((ietat(ie) == 1) .or. (ietat(ie) == -2)) then
                if (smear == 1) then   ! Smeared actif => degradation de vh active
                        call fiss_modul_smeared(ie,vh)
                elseif (depen == 1) then  ! Depenalisation active => on annule la raideur
                         vh=vh/1.d12
                end if
            end if
        end if
end subroutine fiss_modul


!---------------------------------------------------------------------------
!> @author JL Tailhan
!> @brief Calcul de la raideur de l'élément
!
!> @details
!> ### DESCRIPTION:
!> Calcul de la raideur de l'élément en fonction de son état
!> (vierge, fissuré ouvert, fissuré refermé
!>
!>--------------------------------------------------------------------------
subroutine fiss_modul_smeared(ie,vh)

!********************************************************!
!  Fonction donnant acces aux matrices de comportement   !
!           relatives aux elements fissures              !
!********************************************************!

        use variables, only : dime, inorm
        use math, only : INV

        implicit none

        real*8, dimension(:,:), intent(inout) :: vh
        integer, intent(in) :: ie

        real*8 :: Ef, Gf, lx, ly, lz, mx, my, mz, nx, ny, nz
        real*8, dimension(dime,dime) :: Dcr, M0, M1
        real*8, dimension(size(vh,1),size(vh,2)) :: Def, M2
        real*8, dimension(dime) :: nf
        real*8, dimension(size(vh,1),dime) :: M

!********************************************************!

        ! Parametres du comportement fissure

        Ef = 1.0d-10
        Gf = 1.0d-10

        ! Cas de l'element rompu
        if (dime==2) then
                Dcr(1,:) = (/ Ef, 0.d0 /)
                Dcr(2,:) = (/ 0.d0, Gf /)
        elseif (dime==3) then
                Dcr(1,:) = (/ Ef, 0.d0, 0.d0 /)
                Dcr(2,:) = (/ 0.d0, Gf, 0.d0 /)
                Dcr(3,:) = (/ 0.d0, 0.d0, Gf /)
        else
                print *, 'Cas non prevu'
                stop
        end if

        nf = inorm(ie,:)
        lx = nf(1)
        ly = nf(2)
        mx = -ly
        my = lx
        if (dime==2) then
                M(1,:) = (/ lx**2.,     lx*mx /)
                M(2,:) = (/ ly**2.,     ly*my /)
                M(3,:) = (/ 2.d0*lx*ly, lx*my + ly*mx /)
        elseif (dime==3) then
                lz = nf(3)
                mz = 0.0d0
                nx = ly*mz - lz*my
                ny = lz*mx - lx*mz
                nz = lx*my - ly*mx
                M(1,:) = (/ lx**2., lx*mx, lx*nx /)
                M(2,:) = (/ ly**2., ly*my, ly*my /)
                M(3,:) = (/ lz**2., lz*mz, lz*mz /)
                M(4,:) = (/ 2.*lx*ly, lx*my + ly*mx, nx*ly + ny*lx /)
                M(5,:) = (/ 2.*ly*lz, ly*mz + lz*my, ny*lz + nz*ly /)
                M(6,:) = (/ 2.*lz*lx, lz*mx + lx*mz + nz*lx + nx*lz /)
        else
                print *, 'Cas non prevu'
                stop
        end if

        ! Matrice Def = vh * N*inv(Dcr + N'*vh*N)*N' * vh;
        M0 = matmul(transpose(M),matmul(vh,M))
        M1 = INV(Dcr + M0)
        M2 = matmul(M,matmul(M1,transpose(M)))
        Def = matmul(vh,matmul(M2,vh))
        vh = vh - Def

end subroutine fiss_modul_smeared


!---------------------------------------------------------------------------
!> @author JL Tailhan
!> @brief Aiguillage sur les lois de comp (modèle de fissuration semi-explicite)
!
!> @details
!> ### DESCRIPTION:
!> Aiguillage sur les lois de comportement pour le modèle de fissuration
!> semi-explicite
!>
!>--------------------------------------------------------------------------
subroutine fiss_loi(iloi,icomp,ie)


    use variables, only : ietat !, icompfiss

    implicit none

    integer, intent(in) :: icomp, ie
    integer, intent(inout) :: iloi

!********************************************************!

    if (iloi==12) then 
        if ((ietat(ie)==1) .or. (ietat(ie) == -1) .or. (ietat(ie) == -2)) then
        ! Element fissure => aiguillage vers calcul elastique avec raideur nulle
            iloi=1
        elseif (ietat(ie)==2) then
            ! Element (fissure) referme => aiguillage vers calcul loi de Coulomb
            !iloi = 14
            iloi = 1                ! Loi elastique
        else
            iloi = icomp
        endif
    else
        ! Element vierge ou en cours de fissuration => aiguillage vers calcul elasto-plastique
        iloi = icomp
    endif

end subroutine fiss_loi


!---------------------------------------------------------------------------
!> @author JL Tailhan
!> @brief Initialisations (modèle de fissuration semi-explicite)
!
!> @details
!> ### DESCRIPTION:
!> Initialisation des vecteurs pour modèle de fissuration semi-explicite
!> états des éléments, ouvertures et orientations des fissures, historique des
!> états, etc. ...
!>
!>--------------------------------------------------------------------------
subroutine fiss_init()

!********************************************************!
!  Initialisations liees au calcul fissure               !
!********************************************************!

    use variables, only : vprelg, nelt, dime, ietat, inorm, irloc, iouver1, iouver2, &
                          & histetat1, histetat2, listcompfiss, & ! icompfiss, &
                          & fiss, smear, depen, iprint, iedngm
    use initialisation, only : init_mat, init_vec

    implicit none

    integer, dimension(3) :: test_ind
    integer i
    logical ismember
    character reply

!********************************************************!

    fiss = 0
    smear = 0
    depen = 0

    !listcompfiss = (/12, 15/)   ! Numero de lois de comportement pour beton probabiliste

    ismember = .false.
    reply = ''

    do i = 1, size(vprelg,1)
        ! Cas particulier de la loi traduisant la fissuration du beton
        if (any(listcompfiss==vprelg(i,1))) then
            fiss=1
            ismember = .true.
            exit
        end if
    end do

    !ismember = .false.
    if (ismember) then
        ! Choix de la representation de la fissuration (smeared ou depenalisation)
        do while ((reply/='d') .and. (reply/='s'))
            !if(iprint>0) then
            !        print*, ' Bonjour ! Qu''est-ce qu''on va faire aujourd''hui pour s''amuser ?? Smeared (s) ou depenalisation (d) ? (Enter = d) '
            !        read*, reply
            !end if
            !if (reply=='') reply = 'd'

            ! JLT deb : Petite modification pour éviter ce choix interactif
            reply='d'
            ! JLT Fin : Fin de modification

            select case(reply)
               case ('s')
                   !fiss = 1
                   smear = 1
                   if (iprint>0) print*, ' Ok, on va faire du SMEARED ! '
               case ('d')
                   !fiss = 1
                   depen = 1
                   if (iprint>0) print*, ' Ok, on va faire de la DEPENALISATION ! '
            end select
        end do
    end if

    test_ind = (/ fiss, smear, depen /)
    if (all(test_ind /= 0)) stop 'Erreur, verifier les index fiss, depen !!'

    if (fiss == 1 ) then
        iedngm = 0
        if (.not.allocated(ietat)) call init_vec(ietat,nelt)
        call init_vec(iouver1,nelt)
        call init_vec(iouver2,nelt)
        call init_vec(histetat1,nelt)
        call init_vec(histetat2,nelt)
        if (.not.allocated(inorm)) call init_mat(inorm,nelt,dime)
        if (.not.allocated(irloc)) call init_mat(irloc,nelt,dime*dime) ! JLT pour tenir compte de tout le repere
    end if

end subroutine fiss_init


!---------------------------------------------------------------------------
!> @author JL Tailhan
!> @brief Gestion du passage ietat=-2 à ietat=2
!
!> @details
!> ### DESCRIPTION:
!> Passage ietat=-2 à ietat=2
!>
!>--------------------------------------------------------------------------
subroutine fiss_stock()


    use variables, only :  fiss, ietat, nelt
    implicit none

    integer :: ie

! *******************************************************!

    if (fiss==1) then
        where (ietat==-2) ietat = 1
    end if
end subroutine fiss_stock


!---------------------------------------------------------------------------
!> @author JL Tailhan
!> @brief Définition des propriétés aléatoires
!
!> @details
!> ### DESCRIPTION:
!> Fonction de definition des proprietes mecaniques aleatoires
!> @param[in] :
!> - kprop : numerotation des groupes des elements
!> - vprelg : proprietes mecaniques (fonction du modele)
!>
!> @param[out] :
!> - kprop : actualise en fonction du numero de l'element
!>      (1 element prob = 1 groupe)
!> - vprelg : modifie pour prendre en compte la propriete
!>      aleatoire consideree - Attention depend du modele choisi
!>
!>--------------------------------------------------------------------------
subroutine fiss_distal()


    use variables, only : nelt, dime, ktypel, &
            & infele, vprelg, kprop, kprop0, idmax, &
            !& prob, pi, fiss, icompfiss, Dg, fc, &
            & pi, fiss, Dg, fc, vitrav, &
            !--- gestion de l'alea (deb)
            !& alea,iloip, param, iparam
            & alea, young, resist, enerj, enerjpf
            !--- gestion de l'alea (fin)
    use lib_elem, only : elem_B
    use initialisation
    use proprietes_aleatoires
    use math

    implicit none

    real*8, dimension(:,:), allocatable :: vprelg0, ksig, vn, vb
    real*8, dimension(:), allocatable   :: vpg, RR, EE, WW, WWP
    integer, dimension(:), allocatable  :: IGRR, IGEE, IGWW, IGWWP
    real*8, dimension(nelt) :: Vel
    real*8, dimension(1) :: Rt, E, W, WP
    real*8, dimension(5) :: loi1, loi2, loi3, loi4
    real*8 :: detj, epais, rp, &
              & poids, Ve, Vg, r !, &
              !& alpha   ! alpha est utilise pour la determination des parametres du modele elasto-plastique. AM
    integer :: iebe, i, ie, ie0, nbg, ipg, npg, id, iloi, dPlay, igrp, igmp, ine
    logical :: resp, modp, nrjp, iprobVO, iprobFI, iprobEN, iprobMY, iprobWP, nrjppf

!********************************************************!

    !----- Preparation de la renumerotation de kprop et vprelg  !
    ! Principe : on conserve la numerotation pour tous les      !
    ! groupes dont l'ancien numero est inferieur a nga. Puis on !
    ! "decale" d'un cran dans la numerotation tous les groupes  !
    ! qui ont un ancien numero plus grand que nga. Enfin, tous  !
    ! les elements du groupe nga sont renumerote en nga+i.      !

    if ((fiss==1) .and. (alea)) then !---Pour le modèle de fissuration aleatoire
        rp = 1.e10
        ine = 0
        nbg = maxval(kprop)
        call init_mat(vprelg0,(size(kprop,1)+nbg),size(vprelg,2))
        vprelg0(1:nbg,:) = vprelg(1:nbg,:)

        Vel = 0.0
        call init_vec(RR,nelt)
        call init_vec(IGRR,nelt)
        call init_vec(EE,nelt)
        call init_vec(IGEE,nelt)
        call init_vec(WW,nelt)
        call init_vec(IGWW,nelt)
        call init_vec(WWP,nelt)
        call init_vec(IGWWP,nelt)

        !--- compteur pour le décalage de numerotation des groupes dans vprelg
        ie0 = 0
        !--- compteur sur les groupes effectivement probabilises
        iebe = 0

        !----- Indices pour trace de distributions de proprietes aleatoires
        iprobVO = .false.
        iprobFI = .false.
        iprobEN = .false.
        iprobMY = .false.
        iprobWP = .false.


        !-- Boucle sur les elements
        do ie = 1, nelt

            !-- detection d'un element probabiliste (module et/ou resistance et/ou energie)
            resp = .false.
            modp = .false.
            nrjp = .false.
            nrjppf = .false.
            if (allocated(resist)) then
                if (count(resist%num==kprop0(ie))==1) resp=.true.
            end if
            if (allocated(young)) then
                if (count(young%num==kprop0(ie))==1) modp=.true.
            end if
            if (allocated(enerj)) then
                if (count(enerj%num==kprop0(ie))==1) nrjp=.true.
            end if
            if (allocated(enerjpf)) then
                if (count(enerjpf%num==kprop0(ie))==1) nrjppf=.true.
            end if

            !-- si l'element appartient a un groupe probabiliste (module et/ou resistance et/ou energie)
            if (resp .or. modp .or. nrjp .or. nrjppf) then

                !-- on verifie la compatibilite de la loi de comp avec sa formulation probabiliste
                iloi=int(vprelg0(kprop(ie),1))
                !if (iloi==icompfiss) then

                    !-- Incrementation des compteurs
                    ie0 = ie0 + 1
                    iebe = iebe + 1

                    !-- Calcul du rapport des volumes
                    ! recuperation des poids et points d'integration de l'element
                    allocate(vpg(size(infele(ktypel(ie))%W)))
                    allocate(ksig(size(infele(ktypel(ie))%Q,1),size(infele(ktypel(ie))%Q,2)))

                    vpg = infele(ktypel(ie))%W
                    ksig = infele(ktypel(ie))%Q
                    npg = size(ksig,1)

                    epais = 1.d0

                    ! recuperation de l'epaisseur si 2D
                    if ((dime == 2).and.(vprelg0(kprop(ie),2)==3)) epais = vprelg0(kprop(ie),idmax-2) ! En contraintes planes (CP)

                    Ve = 0.d0

                    do ipg = 1, npg
                        poids = vpg(ipg)
                        ! calcul des fonctions d'interpolation et des derivees
                        call elem_B(vn,vb,detj,ksig(ipg,:),ie)
                        deallocate(vb,vn)
                        Ve = Ve + detj*poids*epais
                    end do

                    Vel(ie0) = Ve
                    Vg = (4./3.)*pi*(Dg/2.)**3.
                    r = Ve/Vg

                    if (Ve < 0) then
                    print*, ie
                    r = 1
                    Ve = 0.d0
                    end if

                    rp = minval((/rp,r/))

                    iprobVO = .true.

                    deallocate(ksig,vpg)

                    !-- Distributions aleatoires de proprietes mecaniques

                    !-- Resistance a la traction
                    if (resp) then

                        ! on recupere le numero de groupe probabiliste correspondant a l'element
                        allocate(vitrav(1))
                        vitrav= find_vec(resist%num,kprop0(ie))
                        igrp = vitrav(1)
                        deallocate(vitrav)

                        !loi1(1) = 0.        ! b
                        !loi1(2) = 0.        ! c
                        !loi1(3) = 0.        ! moy
                        !loi1(4) = 0.        ! ect
                        loi1(5)=resist(igrp)%loi      ! num

                        if (loi1(5)==1) then
                            ! Loi normale
                            if (resist(igrp)%ipa/=0) then
                                loi1(3) = resist(igrp)%param(1)
                                loi1(4) = resist(igrp)%param(2)
                            else
                                call rossi_coef(loi1(3),loi1(4),r,fc,'RESI')
                            end if
                        elseif (loi1(5)==2) then
                            ! Loi de Weibull
                            if (resist(igrp)%ipa/=0) then
                                loi1(1) = resist(igrp)%param(1)
                                loi1(2) = resist(igrp)%param(2)
                            else
                                !if ((r<.004d0).or.(r>1000.d0)) stop 'Sortie du domaine de validité du modele'
                                if (dime==2) then
                                   loi1(1) = 1.532d0*dexp(-1.222d0*dlog10(r))+ &
                                   &         6.073d0*dexp(-0.216d0*dlog10(r))
                                   loi1(2) = .3475d0*(dlog10(r))**2 +2.182d0*dlog10(r) +4.249d0
                                elseif (dime==3) then
                                    loi1(1) = 1.574d0*dexp(.2044d0*dlog10(r))+ &
                                    &         6.633d0*dexp(-0.8092d0*dlog10(r))
                                    loi1(2) = 131.2d0*dexp(.3415d0*dlog10(r))- &
                                    &         128.7d0*dexp(.335d0*dlog10(r))
                                end if
                            end if

                        elseif (loi1(5)==3) then
                            ! Loi log-normale
                            if (resist(igrp)%ipa/=0) then
                                loi1(3) = resist(igrp)%param(1)
                                loi1(4) = resist(igrp)%param(2)
                            else
                                call rossi_coef(loi1(3),loi1(4),r,fc,'RESI')
                            end if

                        end if

                        !====> Tirage aleatoire de la resistance
                        Rt = distr_alea(loi1,1,ie0)
                        ! preparation pour trace
                        iprobFI = .true.

                        if (Rt(1) <= 0) then
                             print*, Rt(1)
                             stop ' Resistance negative'
                        end if
                        RR(iebe) = Rt(1)
                        IGRR(iebe) = resist(igrp)%num
                    end if

                    !-- Module d'Young
                    if (modp) then

                        ! on recupere le numero de groupe probabiliste correspondant a l'element
                        allocate(vitrav(1))
                        vitrav= find_vec(young%num,kprop0(ie))
                        igmp = vitrav(1)
                        deallocate(vitrav)

                        loi3(1) = 0.       ! b
                        loi3(2) = 0.       ! c
                        loi3(3) = 0.       ! moy
                        loi3(4) = 0.       ! ect
                        loi3(5) = young(igmp)%loi        ! num

                        if (loi3(5)==1) then
                            ! Loi normale
                            if (young(igmp)%ipa/=0) then
                                loi3(3) = young(igmp)%param(1)
                                loi3(4) = young(igmp)%param(2)
                            else
                                call rossi_coef(loi3(3),loi3(4),r,fc,'MODU',vprelg(kprop(ie),id))
                            end if
                        elseif (loi3(5) == 2) then
                            ! Loi de Weibull
                            stop 'FIDES_interf_distal : loi probabiliste non programmee pour les modules'
                        elseif (loi3(5) == 3) then
                            ! Loi log-normale
                            if (young(igmp)%ipa/=0) then
                                loi3(3) = young(igmp)%param(1)
                                loi3(4) = young(igmp)%param(2)
                            else
                                call rossi_coef(loi3(3),loi3(4),r,fc,'MODU',vprelg(kprop(ie),id))
                            end if
                        end if
                        !====> Tirage aleatoire du module
                        E = distr_alea(loi3,1,ie0)

                        ! preparation pour trace
                        iprobMY = .true.
                        EE(iebe) = E(1)
                        IGEE(iebe) = young(igmp)%num

                    end if! fin if modp

                    !-- Energie plastique
                    if (nrjp) then
                        ! on recupere le numero de groupe probabiliste correspondant a l'element
                        allocate(vitrav(1))
                        vitrav= find_vec(enerj%num,kprop0(ie))
                        igrp = vitrav(1)
                        deallocate(vitrav)

                        loi2(1) = 0.       ! b
                        loi2(2) = 0.       ! c
                        loi2(3) = 0.       ! moy
                        loi2(4) = 0.       ! ect
                        loi2(5) = enerj(igrp)%loi        ! num

                        if (loi2(5)==1) then
                            ! Loi normale
                            loi2(3) = enerj(igrp)%param(1)
                            loi2(4) = enerj(igrp)%param(2)
                        elseif (loi2(5)==2) then
                            ! Loi de Weibull
                            loi2(1) = enerj(igrp)%param(1)
                            loi2(2) = enerj(igrp)%param(2)
                        elseif (loi2(5)==3) then
                            ! Loi log-normale
                            loi2(3) = enerj(igrp)%param(1)
                            loi2(4) = enerj(igrp)%param(2)
                        end if
                        !====> Tirage aleatoire de l'energie
                        W = distr_alea(loi2,1,ie0)
                        ! Modifs pour modele elasto-plastique. AM
                        !W(1) = (Rt(1)**2)*(1.d0 - alpha**2)/2.d0/E(1)
                        ! AM

                        ! preparation pour trace
                        iprobEN = .true.
                        !WW(iebe) = W(1)*Ve
                        WW(iebe) = W(1) ! On enleve le Ve pour la visualisation de W
                        IGWW(iebe) = enerj(igrp)%num

                    end if!fin if nrjp

                    !-- Energie post-fissuration
                    if (nrjppf) then

                        ! on recupere le numero de groupe probabiliste correspondant a l'element
                        allocate(vitrav(1))
                        vitrav= find_vec(enerjpf%num,kprop0(ie))
                        igrp = vitrav(1)
                        deallocate(vitrav)

                        loi4(1) = 0.       ! b
                        loi4(2) = 0.       ! c
                        loi4(3) = 0.       ! moy
                        loi4(4) = 0.       ! ect
                        loi4(5) = enerjpf(igrp)%loi        ! num

                        if (loi4(5)==1) then
                            ! Loi normale
                            loi4(3) = enerjpf(igrp)%param(1)
                            loi4(4) = enerjpf(igrp)%param(2)
                        elseif (loi4(5)==2) then
                            ! Loi de Weibull
                            print*, 'error: loi non adaptee'
                            stop
                        elseif (loi4(5)==3) then
                            ! Loi log-normale
                            loi4(3) = enerjpf(igrp)%param(1)
                            loi4(4) = enerjpf(igrp)%param(2)
                        end if
                        !====> Tirage aleatoire de l'energie
                        WP = distr_alea(loi4,1,ie0)

                        ! preparation pour trace
                        iprobWP = .true.
                        WWP(iebe) = WP(1)
                        IGWWP(iebe) = enerjpf(igrp)%num

                    end if

                    !---- Stockage des valeurs (dans vprelg0 "decale")
                    vprelg0(nbg + ie0,:) = vprelg(kprop(ie),:)
                    kprop(ie) = nbg + ie0

                    !----- Stockage resistance a la traction
                    id = 4
                    if (dime==3) id = id-1

                    if (resp) vprelg0(kprop(ie),id+2) = Rt(1)

                    !----- Stockage modules elastiques
                    if (modp) vprelg0(kprop(ie),id) = E(1)

                    !----- Stockage de l'energie
                    if (nrjp) vprelg0(kprop(ie),id+3) = W(1) * Ve   ! Energie de deformation plastique

                    !----- Stockage de l'energie post-pic
                    if (nrjppf) then

                        if (iloi==12) then
                            WP(1) =  WP(1) - (RT(1)**2)/(2.*E(1))*((Ve/epais)**(1./2.))
                            WP(1) = max(0., WP(1)*Ve)
                        end if

                        if(dime == 2) then
                             vprelg0(kprop(ie),id+3) = WP(1) / ((Ve/epais)**(1./2.))   ! Energie post-pic
                        else if(dime == 3) then
                             vprelg0(kprop(ie),id+3) = WP(1) / (Ve**(1./3.))   ! Energie post-pic
                        end if
                    end if
                !end if
            end if
        end do


        !-- Sorties : on reecrit vprelg a partir de vprelg0 (puis on detruit vprelg0)
        deallocate(vprelg)
        call init_mat(vprelg,ie0+nbg,size(vprelg,2))
        vprelg(1:ie0+nbg,:) = vprelg0(1:ie0+nbg,:)

        !-- Traces graphiques eventuels

        ! Visualisation de la distribution des volumes
        if (iprobVO) then
            dPlay = 1
            if (dplay>=1) then
                print*,' '
                print*,'Maillage beton probabiliste :'
                print*,'-----------------------------'
            end if
            call distrib_volumes(Vel,Dg,50,dPlay)
        end if

        ! Visualisation de la distribution de resistance
        if (iprobFI) then
            dPlay = 1
            if (dplay>=1) then
                print*,'Fissuration probabiliste :'
                print*,'--------------------------'
                print'(a10,e12.5,a16,e12.5)','- Dmax :',Dg,' - Resistance :',fc
                print'(a10,e12.5)','- Vg :',Vg
            end if
            do i=1,size(resist%num)
                loi1(5)=resist(i)%loi
                call alea_trac(RR,loi1,200,dPlay,'resistance',IGRR,resist(i)%num)
            end do
        end if

        ! Visualisation de la distribution de module elastique
        if (iprobMY) then
            dPlay = 1
            do i=1,size(young%num)
                loi3(5)=young(i)%loi
                call alea_trac(EE,loi3,200,dPlay,'module',IGEE,young(i)%num)
            end do
        end if

        ! Visualisation de la distribution de l'energie
        if (iprobEN) then
            dPlay = 1
            do i=1,size(enerj%num)
                loi3(5)=enerj(i)%loi
                call alea_trac(WW,loi2,200,dPlay,'energie',IGWW,enerj(i)%num)
            end do
        end if

        if (iprobWP) then
            dPlay = 1
            do i=1,size(enerjpf%num)
                loi3(5)=enerjpf(i)%loi
                call alea_trac(WWP,loi4,200,dPlay,'energie post-fissuration',IGWWP,enerjpf(i)%num)
            end do
        end if

        deallocate(vprelg0,RR,EE,WW,WWP,IGRR,IGEE,IGWW,IGWWP)
    end if

    print*, 'Le rapport du volume le plus petit : ', rp

end subroutine fiss_distal
!-------------------------------------------------------!

end module fissuration

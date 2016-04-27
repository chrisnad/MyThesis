!------------------------------------------------------------------------------------------------------
! MODULE: fissuration_BA
!
!> @author C Nader      (version 1.0 - 2016)
!
!> @brief
!> Ensemble des routines relatives au(x) modele(s) macroscopique(s) orthotrope(s)
!> de fissuration.
!------------------------------------------------------------------------------------------------------
module fissuration_BA


contains


!------------------------------------------------------------------------------------------------------
!> @author C Nader      (version 1.0 - 2016)
!
!> @brief Routine de calcul des contraintes pour modËles d'endommagement orthotrope.
!
!> @details
!> ### DESCRIPTION:
!> Routine de calcul des contraintes pour modËles d'endommagement
!------------------------------------------------------------------------------------------------------
subroutine fissBA_sig_endo_ortho(vsig,vin,vh,vb,vprel,vdle,iloi,ie)

    !**************************************************************************************
    !* DEFINITIONS
    !**************************************************************************************

    !--------------------------------------------------------------------------------------
    !--- Modules generaux
    use variables  , only : dime, ietat, iedng

    implicit none

    !--------------------------------------------------------------------------------------
    !--- Variables in/out
    real*8, dimension(:,:), intent(in) :: vb
    real*8, dimension(:,:), intent(inout) :: vh
    real*8, dimension(:), intent(in) :: vdle, vprel
    integer, intent(in) :: iloi, ie
    real*8, dimension(:), intent(inout) :: vsig, vin

    !--------------------------------------------------------------------------------------
    !--- Variables locales
    !--- Quantites globales
    integer :: id

    !--- Quantites principales
    real*8, dimension(size(vb,1)) :: veps, vepsloc, vepslim

    ! Variables pour le modele d'endommagement
    real*8 :: Dini, vars, D, epsilon1, E1
    real*8 :: ouv, Dnul1, Dnul2

    !--- Variables JLT
    real*8, dimension(size(vh,1),size(vh,2)) ::  P, Pe !, Q
    real*8, dimension(dime,dime) :: V
    real*8, dimension(size(vh,1)) :: vsigm
    logical :: vcrit, knull
    integer :: nc1, i, j
    real*8 :: var0, var
    real*8 :: F_var  !varcrit, varc, F_var, RT, W, E, nu
    real*8 :: sigpa, Ea, eps0, epsc, epsl
    real*8 :: D0, D0ini, Dk, RBA

    !**************************************************************************************
    !* CORPS DE PROCEDURE
    !**************************************************************************************

    !--------------------------------------------------------------------------------------
    !--- Calcul de l'ouverture de fissure initiale (pour non-interpenetration)
    ouv = 0.D0 ; Dnul1 = 1.d0 ; Dnul2 = 1.d0 ; D = 0. ; Dk = 1.d0

    !--------------------------------------------------------------------------------------
    !--- Selection en fonction de la loi
    select case(iloi)

        !----------------------------------------------------------------------------------
        !--- Loi d'endommagement isotrope beton + armatures
        case(27)

            !------------------------------------------------------------------------------
            !--- Initialisations
            Dini = vin(1) ! endommagement initial
            D0ini= vin(2) ! seuil initial d'endommagement
            vars = vin(3) ! seuil courant d'endommagement

            !------------------------------------------------------------------------------
            !--- Etat de contrainte et de deformation
            veps  = matmul(vb,vdle)
            vsigm = matmul(vh,veps)
            vepslim = 0.d0

            !------------------------------------------------------------------------------
            !--- Parametres specifiques a la loi
            if (dime == 2) then
                id = 4
            elseif (dime == 3) then
                id = 3
            else
                STOP '1D NON IMPLEMENTE pour loi 27'
            end if

            E1 = vprel(id)

            Rba   = vprel(id+dime**2)            !    \
            eps0  = Rba/E1                       !     | Parametres
            Ea    = vprel(id+dime**2+1)          !     | de la loi
            epsl  = eps0 + vprel(id+dime**2+2)   !    /

            sigpa = (epsl - eps0) * Ea   ! debut de non-linearite du comportement des aciers
            epsc = -Ea*eps0/(E1-Ea)

            !------------------------------------------------------------------------------
            !--- Recuperation du repere principal de ferraillage
            do i = 1, dime
                do j = 1, dime
                    V(i, j) = 0.d0
                end do
                V(i, i) = 1.d0
            end do

            !------------------------------------------------------------------------------
            !--- Calcul des matrices de changement de repere (Global->ferraillage)
            P = fissBA_changement_repere(V,nc1,1)
            Pe = fissBA_changement_repere(V,nc1,2,1)
            
            !------------------------------------------------------------------------------
            !--- Dans le repere du ferraillage...
            vepsloc = matmul(P,veps)

            !------------------------------------------------------------------------------
            !--- Recuperation de la deformation correspondante dans l'axe principal de ferraillage
            epsilon1 = vepsloc(1)

            !------------------------------------------------------------------------------
            !--- Test pour activation du calcul de l'endommagement
            vcrit = .false.
            knull = .false.
            if (ie == iedng) then
                !----- Changement "vcrit"
                vcrit = .true.
                var0  = epsl-eps0
                !----- Initialisation du seuil de deformation courant
                vars  = var0
                F_var = 0.d0
                !----- Changement etat de l'element
                ietat(ie) = 3
            else
                if(ietat(ie) == 3) then
                    vcrit = .true.
                    var0  = epsl-eps0
                end if
            end if

            !------------------------------------------------------------------------------
            !--- Initialisation de l'endommagement
            D0 = D0ini
            D  = Dini

            !------------------------------------------------------------------------------
            !--- Declanchement et evolution de l'endommagement
            if(vcrit .eqv. .true.) then

                !--------------------------------------------------------------------------
                !--- Actualisation du seuil de rupture
                if (epsilon1 < epsc) then
                    D0 = 0.d0
                    Dnul1 = 0.d0
                    Dnul2 = 0.d0

                elseif (epsc <= epsilon1 .and. epsilon1 < eps0) then
                    D0 = 1.d0 - (Ea/E1)
                    Dnul1 = 1.d0
                    Dnul2 = 0.d0
                    vepslim(1) = eps0

                elseif (eps0 <= epsilon1 .and. epsilon1 < epsl) then
                    D0 = 1.d0 - (Ea/E1)
                    Dnul1 = 1.d0
                    Dnul2 = 0.d0
                    if (vars >= epsl) Dnul2 = 1.d0
                    vepslim(1) = eps0

                !--------------------------------------------------------------------------
                !--- Loi d'evolution d'endommagement
                elseif (epsl <= epsilon1) then
                    D0 = 1.d0 - (Ea/E1)
                    Dnul1 = 1.d0
                    Dnul2 = 1.d0

                    !---- Variable de pilotage de l'endommagement
                    var = epsilon1-eps0
                    vepslim(1) = eps0

                    if (var >= vars) then
                        !---- Actualisation du seuil de rupture
                        vars = var

                        !---- Calcul de l'endommagement
                        F_var = 1.d0 - (sigpa/(Ea*var0))
                        D = 1.d0 - var0/var * (1.d0 - F_var)

                    end if

                end if

            elseif(knull .eqv. .true.) then

                Dk = 0.d0

            end if

            if (D < 0.d0) D = 0.d0
            D = max(Dini,D)
            D = min(D,1.d0)
            !------------------------------------------------------------------------------
            !--- Contraintes et Rigidité
            vh(1,1) = (1.d0-Dnul2*D)*vh(1,1)
            vh(1,2) = (1.d0-Dnul2*D)*vh(1,2)
            vh(2,1) = (1.d0-Dnul2*D)*vh(2,1)
            if (dime == 3) then
                vh(1,3) = (1.d0-Dnul2*D)*vh(1,3)
                vh(3,1) = (1.d0-Dnul2*D)*vh(3,1)
            end if

            vh = (1.d0-Dnul1*D0)*Dk*vh
            vsig = matmul(vh,veps-vepslim)
            vin(1) = D
            vin(2) = D0
            vin(3) = vars

        case default

            stop 'fissBA_sig_endo_ortho : cas non implante'

    end select

end subroutine fissBA_sig_endo_ortho


!------------------------------------------------------------------------------------------------------
!> @author JL Tailhan
!> @author C Nader      (version 1.0 - 2016)
!
!> @brief Changement de repère des contraintes (global/principal)
!
!> @details
!> ### DESCRIPTION:
!> Definition de la matrice de changement de repere :
!> du repere global au repere principal
!------------------------------------------------------------------------------------------------------
!   ind = 1 matrice de changement de rep global -> local
!   ind = 2 matrice de changement inverse (transposee de la prec.)
!   iopt = 1 travail sur les contraintes
!   iopt = 2 travail sur les deformations
!------------------------------------------------------------------------------------------------------
function fissBA_changement_repere(V,n,iopt,typind) result (M)

    !**************************************************************************************
    !* DEFINITIONS
    !**************************************************************************************

    !--------------------------------------------------------------------------------------
    !--- Modules generaux
    use variables, only : dime
    use math, only : norme

    implicit none

    !--------------------------------------------------------------------------------------
    !--- Variables in/out
    real*8, dimension(dime,dime), intent(in) :: V
    integer, intent(in) :: n, iopt
    integer, intent(in), optional :: typind
    !--- Output result
    real*8, dimension(n,n) :: M

    !--------------------------------------------------------------------------------------
    !--- Variables locales
    real*8 :: n1, n2, n3, a, b, c, d, e, f, g, h, i
    integer:: ind

    !**************************************************************************************
    !* CORPS DE PROCEDURE
    !**************************************************************************************

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

    ind = 1
    if (present(typind)) then
        ind = typind
    endif

    !--------------------------------------------------------------------------------------
    !--- Pour les contraintes
    if (iopt==1) then

        !----------------------------------------------------------------------------------
        !--- Matrice de passage : Sig_p = Ps Sig
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

        !----------------------------------------------------------------------------------
        !--- Matrice de passage inverse : Sig = Transpose(Ps) Sig_p
        elseif (ind==2) then
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

    !--------------------------------------------------------------------------------------
    !--- Pour les deformations
    elseif (iopt==2) then

        !----------------------------------------------------------------------------------
        !--- Matrice de passage Eps_p = Pe Eps :
        if (ind==1) then
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

        !----------------------------------------------------------------------------------
        !--- Matrice de passage inverse : Eps = Transpose(Pe) Eps_p
        elseif (ind==2) then
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

end function fissBA_changement_repere


!---------------------------------------------------------------------------
!> @author C Nader      (version 1.0 - 2016)
!
!> @brief Affecter des valeurs aleatoires aux parametres de la loi BA
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
!>--------------------------------------------------------------------------
subroutine fissBA_distal()

    !**************************************************************************************
    !* DEFINITIONS
    !**************************************************************************************

    !--------------------------------------------------------------------------------------
    !--- Modules generaux
    use variables, only :   nelt, dime, &
    &                       vprelg, kprop, kprop0, fissBA, vitrav, &
    !    gestion de l'alea (deb)
    &                       alea, enerjba, enerjba2
    use initialisation
    use proprietes_aleatoires
    use math

    implicit none

    !--------------------------------------------------------------------------------------
    !--- Variables locales
    real*8, dimension(:,:), allocatable :: vprelg0
    real*8, dimension(1) :: PBA1, PBA2
    real*8, dimension(5) :: loi1, loi2

    integer :: ie, ie0, nbg, id, iloi, igrp
    logical :: nrjbap, nrjbap2

    !**************************************************************************************
    !* CORPS DE PROCEDURE
    !**************************************************************************************

    !--------------------------------------------------------------------------------------
    !--- Preparation de la renumerotation de kprop et vprelg
    !    Principe : on conserve la numerotation pour tous les
    !    groupes dont l'ancien numero est inferieur a nga. Puis on
    !    "decale" d'un cran dans la numerotation tous les groupes
    !    qui ont un ancien numero plus grand que nga. Enfin, tous
    !    les elements du groupe nga sont renumerote en nga+i.

    !--------------------------------------------------------------------------------------
    !--- Pour le modèle de fissuration aleatoire
    if ((fissBA==1) .and. (alea)) then

        nbg = maxval(kprop)
        call init_mat(vprelg0,(size(kprop,1)+nbg),size(vprelg,2))
        vprelg0(1:nbg,:) = vprelg(1:nbg,:)

        ie0 = 0

        !----------------------------------------------------------------------------------
        !--- Boucle sur les elements
        do ie = 1, nelt

            !------------------------------------------------------------------------------
            !--- detection d'un element probabiliste (module et/ou resistance et/ou energie)
            nrjbap = .false.
            if (allocated(enerjba)) then
                if (count(enerjba%num==kprop0(ie))==1) nrjbap=.true.
            end if

            nrjbap2 = .false.
            if (allocated(enerjba2)) then
                if (count(enerjba2%num==kprop0(ie))==1) nrjbap2=.true.
            end if

            !------------------------------------------------------------------------------
            !--- si l'element appartient a un groupe probabiliste (energie)
            if (nrjbap .or. nrjbap2) then

                !--------------------------------------------------------------------------
                !--- on verifie la compatibilite de la loi de comp avec sa formulation probabiliste
                iloi=int(vprelg0(kprop(ie),1))

                    ie0 = ie0 + 1

                    !----------------------------------------------------------------------
                    !--- Distributions aleatoires de proprietes mecaniques

                    !----------------------------------------------------------------------
                    !--- Energie plastique (ce test a ete conserve au cas ou...)
                    if (nrjbap) then

                        ! on recupere le numero de groupe probabiliste correspondant a l'element
                        allocate(vitrav(1))
                        vitrav= find_vec(enerjba%num,kprop0(ie))
                        igrp = vitrav(1)
                        deallocate(vitrav)

                        loi1(1) = 0.d0     ! b
                        loi1(2) = 0.d0     ! c
                        loi1(3) = 0.d0     ! moy
                        loi1(4) = 0.d0     ! ect
                        loi1(5) = enerjba(igrp)%loi        ! num

                        if (loi1(5)==1) then
                            ! Loi normale
                            loi1(3) = enerjba(igrp)%param(1)
                            loi1(4) = enerjba(igrp)%param(2)
                        elseif (loi1(5)==2) then
                            ! Loi de Weibull
                            loi1(1) = enerjba(igrp)%param(1)
                            loi1(2) = enerjba(igrp)%param(2)
                        elseif (loi1(5)==3) then
                            ! Loi log-normale
                            loi1(3) = enerjba(igrp)%param(1)
                            loi1(4) = enerjba(igrp)%param(2)
                        end if
                        !====> Tirage aleatoire de l'energie
                        PBA1 = distr_alea(loi1,1,ie0)

                    end if !fin if nrjp

                    !----------------------------------------------------------------------
                    !--- Energie plastique 2
                    if (nrjbap2) then
                        ! on recupere le numero de groupe probabiliste correspondant a l'element
                        allocate(vitrav(1))
                        vitrav= find_vec(enerjba2%num,kprop0(ie))
                        igrp = vitrav(1)
                        deallocate(vitrav)

                        loi2(1) = 0.d0     ! b
                        loi2(2) = 0.d0     ! c
                        loi2(3) = 0.d0     ! moy
                        loi2(4) = 0.d0     ! ect
                        loi2(5) = enerjba2(igrp)%loi        ! num

                        if (loi2(5)==1) then
                            ! Loi normale
                            loi2(3) = enerjba2(igrp)%param(1)
                            loi2(4) = enerjba2(igrp)%param(2)
                        elseif (loi2(5)==2) then
                            ! Loi de Weibull
                            loi2(1) = enerjba(igrp)%param(1)
                            loi2(2) = enerjba(igrp)%param(2)
                        elseif (loi2(5)==3) then
                            ! Loi log-normale
                            loi2(3) = enerjba2(igrp)%param(1)
                            loi2(4) = enerjba2(igrp)%param(2)
                        end if

                        !====> Tirage aleatoire de l'energie
                        PBA2 = distr_alea(loi2,1,ie0)

                    end if !fin if nrjp

                    !----------------------------------------------------------------------
                    !--- Stockage des valeurs
                    vprelg0(nbg + ie0,:) = vprelg(kprop(ie),:)
                    kprop(ie) = nbg + ie0

                    id = 1
                    if (dime==3) id = 5

                    !----------------------------------------------------------------------
                    !--- Stockage de l'energie BA
                    if (nrjbap .or. nrjbap2) then

                        if (iloi==26.or.iloi==27) then  ! Loi orthotrope endommageable
                            vprelg0(kprop(ie),id+7) = PBA1(1)
                        end if

                    end if

                    !----------------------------------------------------------------------
                    !--- Stockage de l'energie BA2
                    if (nrjbap2) then

                        if (iloi==26.or.iloi==27) then  ! Loi orthotrope endommageable
                            vprelg0(kprop(ie),id+8) = PBA2(1)
                        end if

                    end if

            end if
        end do

        !----------------------------------------------------------------------------------
        !--- Sorties :
        deallocate(vprelg)
        call init_mat(vprelg,ie0+nbg,size(vprelg,2))
        vprelg(1:ie0+nbg,:) = vprelg0(1:ie0+nbg,:)
        deallocate(vprelg0)

    end if

end subroutine fissBA_distal


!---------------------------------------------------------------------------
!> @author JL Tailhan
!> @author C Nader      (version 1.0 - 2016)
!
!> @brief Pilotage du calcul sur l'avancement de la fissuration
!
!> @details
!> ### DESCRIPTION:
!> Pilotage du calcul sur l'avancement de la fissuration pourle modèle de
!> fissuration semi-explicite
!>--------------------------------------------------------------------------
subroutine fissBA_pilot(imetpilo,vsol,vduI,vduII,dlam,alphaba,MOTba,elefiss,iedg)

    !**************************************************************************************
    !* DEFINITIONS
    !**************************************************************************************

    !--------------------------------------------------------------------------------------
    !--- Modules generaux
    use variables, only :   dime, nelt, idmax, ktypel, nomtype, infele, &
    &                       vprelg, kprop, ietat, iedngba, fissBA, pi,  &
    &                       mrtrav, vrtrav, kloce, fc
    use initialisation, only : init_vec, init_mat
    use lib_elem, only : elem_B, elem_kloce2, elem_hooke
    use math
    use utilitaire, only : princ

    implicit none

    !--------------------------------------------------------------------------------------
    !--- Variables in/out
    real*8, dimension(:), intent(in) :: vduII, vsol
    integer, intent(in) :: imetpilo
    integer, intent(in), optional :: iedg
    real*8, dimension(:), intent(in) :: vduI
    !real*8, dimension(dime), intent(out) :: normfissba
    real*8, intent(in) :: dlam
    real*8, intent(out) :: alphaba
    character(len=8), intent(inout) :: MOTba
    integer, intent(inout) :: elefiss

    !--------------------------------------------------------------------------------------
    !--- Variables locales
    real*8, dimension(:), allocatable ::    vdle, vdl0, vsi0, vdsi, vsi, vep, &
    &                                       vsi0loc, vdsiloc, veploc, alpg, vdlI, vdlII
    real*8, dimension(:,:), allocatable ::  vn, vb, vh, ksig, Ps, Pe

    real*8, dimension(dime, dime) :: rotmatba
    real*8 ::   vprel(idmax), alph(nelt), vbid(dime), vbid2(1,dime*dime), V(dime,dime)
    real*8 ::   RT, RC, vcrit1, vcrit2, detj, &
    &           s1, s3, alpt, alpc
    character(len=5) :: typel
    character(len=8) :: MOT
    character(len=8), dimension(nelt) :: rupt
    integer :: ie, ipg, npg, id, iloi, ndle, imin, imax
    integer :: binf, bsup, iedang

    !**************************************************************************************
    !* CORPS DE PROCEDURE
    !**************************************************************************************

    id = 3
    if (dime == 2) then
        id = 4
    elseif (dime == 3) then
        id = 3
    end if

    !--------------------------------------------------------------------------------------
    !--- Switch en fonction de l'option de calcul
    !    iedang => fiss_pilot applique a iedang seulement (pour obtenir le repere local)
    binf = 1
    bsup = nelt
    if (present(iedg)) then
        iedang = iedg
        binf = iedang
        bsup = iedang
    else
        iedang = 0
    endif


    !--------------------------------------------------------------------------------------
    !--- Switch en fonction de la methode de pilotage
    select case(imetpilo)

        !----------------------------------------------------------------------------------
        !--- Pilotage sur element le plus dangereux par recalcul de vduI
        case(1)

            !------------------------------------------------------------------------------
            !--- Pour le cas beton fissurant
            if (fissBA==1) then

                !--------------------------------------------------------------------------
                !--- Initialisations
                alph = 1.d0
                iedngba = 0

                !--------------------------------------------------------------------------
                !--- Boucle sur les elements
                do ie = binf, bsup

                    MOT = '  '
                    alpt = 1.d0
                    alpc = 1.d0
                    typel = nomtype(ktypel(ie))

                    !----------------------------------------------------------------------
                    !--- Pour les elements BA vierges
                    if ((typel == 'MBQ4' .or. typel == 'MBT3' .or. typel == 'MTH8')  &
                        & .and. (ietat(ie)==0)) then

                        !------------------------------------------------------------------
                        !--- Proprietes elementaires
                        vprel = vprelg(kprop(ie),1:idmax)

                        !------------------------------------------------------------------
                        !--- Recuperation des parametres de la loi
                        iloi = int(vprel(1))

                        !if (iloi==1) return

                        !------------------------------------------------------------------
                        !--- Recuperation des informations sur les elements
                        allocate(ksig(size(infele(ktypel(ie))%Q,1), size(infele(ktypel(ie))%Q,2)))
                        ksig = infele(ktypel(ie))%Q
                        npg = size(ksig,1)

                        !------------------------------------------------------------------
                        !--- Initialisations relatives aux elements
                        call init_vec(vrtrav,npg*dime**2)
                        allocate(alpg(npg))
                        alpg = 1.d0

                        !------------------------------------------------------------------
                        !--- Detection de la premiere fissuration (ouverture)
                        ipg = 1

                        !------------------------------------------------------------------
                        !--- Vecteur kloce de localisation dans vecteur global
                        call elem_kloce2( ie, ndle)
                        allocate(vdl0(ndle));  allocate(vdle(ndle))
                        allocate(vdlI(ndle));  allocate(vdlII(ndle))

                        !------------------------------------------------------------------
                        !--- Deplacements
                        vdl0  = vsol(kloce(1:ndle))
                        vdlI  = vduI(kloce(1:ndle))
                        vdlII = vduII(kloce(1:ndle))
                        vdle  = vdlI + dlam*vdlII

                        !------------------------------------------------------------------
                        !--- Matrice d'elasticite vh, en fonction du type d'element
                        call elem_hooke(vh,nomtype(ktypel(ie)),vprel)

                        !------------------------------------------------------------------
                        !--- Calcul des contraintes au centre de gravite de l'element
                        call elem_B(vn,vb,detj,ksig(ipg,:),ie)

                        !--- Dans le repere global...
                        allocate(vsi0(size(vb,1)))
                        allocate(vdsi(size(vb,1)))
                        vsi0  = matmul(vh,matmul(vb,vdl0))
                        vdsi  = matmul(vh,matmul(vb,vdle))

                        !------------------------------------------------------------------
                        !--- Calcul et rangement par ordre decroissant des contraintes principales
                        allocate(vsi(size(vsi0,1)))
                        allocate(vep(size(vdl0,1)))
                        vsi = vsi0 + vdsi
                        vep = matmul(vb, vdl0 + vdle)
                        call princ(vsi,vbid,V)
                        imax = 1
                        imin = 2
                        if (dime==3) imin = 3

                        !------------------------------------------------------------------
                        !--- Recuperation de la matrice de changement de repere global --> local principal
                        call init_mat(Ps,size(vh,1),size(vh,2))
                        call init_mat(Pe,size(vh,1),size(vh,2))
                        Ps = fissBA_changement_repere(V,size(vb,1),1,1) !--- changement de repere pour les contraintes
                        Pe = fissBA_changement_repere(V,size(vb,1),2,1) !--- changement de repere pour les deformations

                        !--- Dans le repere principal...
                        allocate(vdsiloc(size(vh,1)))
                        allocate(vsi0loc(size(vh,1)))
                        allocate(veploc(size(vh,1)))
                        vdsiloc = matmul(Ps,vdsi)
                        vsi0loc = matmul(Ps,vsi0)
                        veploc  = matmul(Pe,vep)

                        deallocate(vdl0,vdlI,vdlII,vdle,vep,vh,vb,vn)
                        deallocate(vsi0,vdsi,vsi,Ps,Pe)

                        !------------------------------------------------------------------
                        !--- Cas de la loi 27
                        if (iloi==27) then

                            !--------------------------------------------------------------
                            !--- Recuperation des parametres de la loi
                            RT = 1.001d0 * vprel(id+dime**2)     ! Contrainte limite en traction pure
                            RC = 1.001d0 * fc                    ! Contrainte limite en cisaillement

                            !--------------------------------------------------------------
                            !--- Recuperation des contraintes principales
                            s1 = vsi0loc(imax)+vdsiloc(imax)
                            s3 = vsi0loc(imin)+vdsiloc(imin)

                            !--------------------------------------------------------------
                            !--- Critere de rupture en ouverture pure
                            vcrit1 = s1 - RT

                            !--------------------------------------------------------------
                            !--- Critere de rupture en cisaillement
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
                                            rotmatba = rotmat(-pi/4.d0, dime)
                                        else
                                            rotmatba = rotmat(-pi/4.d0, dime, V(:, 2))
                                        end if
                                    elseif (veploc(size(veploc)) < 0.d0) then
                                        if (dime == 2) then
                                            rotmatba = rotmat(+pi/4.d0, dime)
                                        else
                                            rotmatba = rotmat(+pi/4.d0, dime, V(:, 2))
                                        end if
                                else
                                    if (veploc(1) < veploc(dime)) then
                                        if (dime == 2) then
                                            rotmatba = rotmat(pi/4.d0, dime)
                                        else
                                            rotmatba = rotmat(pi/4.d0, dime, V(:, 2))
                                        end if
                                    else
                                        if (dime == 2) then
                                            rotmatba = rotmat(-pi/4.d0, dime)
                                        else
                                            rotmatba = rotmat(-pi/4.d0, dime, V(:, 2))
                                        end if
                                    end if
                                end if
                                !normfissba = matmul(rotmatba, V(:, 1))
                                if (dime==2) then
                                    V(:, 1) = matmul(rotmatba, V(:, 1))
                                    V(:, 1) = V(:, 1) / norme(V(:, 1))
                                    V(:, 2) = matmul(rotmatba, V(:, 2))
                                    V(:, 2) = V(:, 2) / norme(V(:, 2))
                                elseif (dime==3) then
                                    V(:, 1) = matmul(rotmatba, V(:, 1))
                                    V(:, 1) = V(:, 1) / norme(V(:, 1))
                                    V(:, 2) = matmul(rotmatba, V(:, 2))
                                    V(:, 2) = V(:, 2) / norme(V(:, 2))
                                    V(:, 3) = matmul(rotmatba, V(:, 3))
                                    V(:, 3) = V(:, 3) / norme(V(:, 3))
                                endif
                            elseif (MOT == 'ouvert') then
                                !normfissba = V(:, 1)
                                V(:, 1) = V(:, 1) / norme(V(:, 1))
                                V(:, 2) = V(:, 2) / norme(V(:, 2))
                                if (dime==3) V(:, 3) = V(:, 3) / norme(V(:, 3))
                            else
                            end if
                        ! Fin du cas loi 27

                        elseif (iloi==1 .or. iloi==12) then

                        else
                            stop 'FIDES_fissBA_pilot : loi non encore programmee'
                        end if

                        deallocate(vsi0loc,vdsiloc,veploc)

                        !------------------------------------------------------------------
                        !--- On conserve le repere local sur le point de gauss considere
                        vbid2=reshape(V,(/1,dime*dime/))
                        vrtrav((ipg-1)*dime*dime+1:ipg*dime*dime)=vbid2(1,:)

                        !------------------------------------------------------------------
                        !--- On conserve
                        rupt(ie) = MOT
                        alph(ie) = minval(alpg)

                        !------------------------------------------------------------------
                        !--- On conserve le repere local du point de Gauss le plus dangereux
                        ipg = find_num(alpg,minval(alpg))
                        mrtrav(ie,:)=vrtrav((ipg-1)*dime*dime+1:(ipg*dime*dime))

                        deallocate(vrtrav)
                        deallocate(alpg,ksig)
                    end if
                end do
                !--- Fin de la boucle sur les elements
                !--------------------------------------------------------------------------

                alphaba = minval(alph)

                !--------------------------------------------------------------------------
                !--- Stockage du numero de l'element le plus dangereux
                if (alphaba < 1.d0) then
                    iedngba = find_num(alph,alphaba)
                    elefiss = count(alph<=0)
                    MOTba = rupt(iedngba)
                else
                    alphaba = 1.d0
                end if
            end if

        !----------------------------------------------------------------------------------
        !--- Pilotage sur l'increment de contrainte avec recalcul de "dlam"
        case(2)

            !--- A faire ...

        case default
            stop 'FIDES_fissBA_pilot : cas non encore programme'
    end select

end subroutine fissBA_pilot


!---------------------------------------------------------------------------
!> @author JL Tailhan
!> @author C Nader      (version 1.0 - 2016)
!
!> @brief Aiguillage sur les lois de comp (modèle de fissuration semi-explicite)
!
!> @details
!> ### DESCRIPTION:
!> Aiguillage sur les lois de comportement pour le modèle de fissuration
!> semi-explicite
!>--------------------------------------------------------------------------
subroutine fissBA_loi(iloi,icomp)

    !**************************************************************************************
    !* DEFINITIONS
    !**************************************************************************************

    implicit none

    !--------------------------------------------------------------------------------------
    !--- Variables in/out
    integer, intent(in) :: icomp
    integer, intent(inout) :: iloi

    !**************************************************************************************
    !* CORPS DE PROCEDURE
    !**************************************************************************************

    iloi = icomp

end subroutine fissBA_loi


!------------------------------------------------------------------------------------------------------
!> @author JL Tailhan
!> @author C Nader      (version 1.0 - 2016)
!
!> @brief Initialisations (modèle de fissuration element macro BA)
!
!> @details
!> ### DESCRIPTION:
!> A faire...
!------------------------------------------------------------------------------------------------------
subroutine fissBA_init()

    !**************************************************************************************
    !* DEFINITIONS
    !**************************************************************************************

    !--------------------------------------------------------------------------------------
    !--- Modules generaux
    use variables,      only : vprelg, fissBA, nbrgrp, nelt, dime, ietat, inorm, listcompfissBA
    use initialisation, only : init_mat, init_vec

    implicit none

    !--------------------------------------------------------------------------------------
    !--- Variables locales
    integer :: i

    !**************************************************************************************
    !* CORPS DE PROCEDURE
    !**************************************************************************************

    fissBA = 0

    do i = 1, nbrgrp

        !----------------------------------------------------------------------------------
        !--- Cas particulier de la loi traduisant la fissuration du beton
        if (any(listcompfissBA==vprelg(i,1))) then
            fissBA = 1
            if (.not.allocated(ietat)) call init_vec(ietat,nelt)
            if (.not.allocated(inorm)) call init_mat(inorm,nelt,dime)
            exit
        endif

    enddo

end subroutine fissBA_init

!------------------------------------------------------------------------------------------------------

end module fissuration_BA

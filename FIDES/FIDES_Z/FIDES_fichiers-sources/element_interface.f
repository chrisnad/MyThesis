!----------------------------------------------------------------------------------------------------------
! MODULE: element_interface
!
!> @author JL Tailhan
!> @author J. Goncalvez (version 1.0 - 2010)
!
!> @brief
!> Ensemble des routines relative aux éléments d'interface pour modéliser la fissuration.
!----------------------------------------------------------------------------------------------------------
module element_interface


contains


!------------------------------------------------------------------------------------------------------
!> @author JL Tailhan
!> @authors J. Goncalvez (version 1.0 - 2010)
!
!> @brief
!> Calcul des contraintes pour la loi de contact
!
!> @details
!> #### DESCRIPTION:
!>  Calcul des contraintes normale et tangente, dans le respect de la loi de Coulomb
!------------------------------------------------------------------------------------------------------
subroutine interf_sig(vhep,vsig,vnle,vh,vb,vnl0,vprel,vdle,iloi,ie,ipg)

    !**************************************************************************************
    !* DEFINITIONS
    !**************************************************************************************

    !--------------------------------------------------------------------------------------
    !--- Modules generaux
    use variables, only :   ietatpg, ktypel, nomtype, imet
    use lib_elem, only : elem_hooke
    use math, only : vecpro

    implicit none

    !--------------------------------------------------------------------------------------
    !--- Variables in/out
    real*8, dimension(:), intent(in) :: vnl0, vdle, vprel
    real*8, dimension(:,:), intent(inout), allocatable :: vh
    real*8, dimension(:,:), intent(in) :: vb
    integer, intent(in) :: ie, ipg, iloi
    real*8, dimension(:), intent(inout) :: vsig, vnle
    real*8, dimension(:,:), intent(out) :: vhep

    !--------------------------------------------------------------------------------------
    !--- Variables locales
    real*8, dimension(size(vh,1)) :: vdfdsig, vdgdsig
    real*8, dimension(size(vb,1)) :: vepspl, vsigm, vsigi
    real*8 :: vcrit, dlam
    logical :: nconv
    integer :: it
    character(len=5) :: calcu

    !**************************************************************************************
    !* CORPS DE PROCEDURE
    !**************************************************************************************

    !--------------------------------------------------------------------------------------
    !--- Recuperation des grandeurs de non-linearite
    vepspl = vnl0(1:size(vb,1))            ! deformations plastiques initiales

    !--------------------------------------------------------------------------------------
    !--- Calcul de sigma + dsigma
    vsigm = matmul(vh,(matmul(vb,vdle)-vepspl))

    !--------------------------------------------------------------------------------------
    !--- Calcul des contraintes verifiant le critere
    calcu = 'D1F'
    call interf_crit(iloi,vprel,vsigm,ie,ipg,vcrit,calcu,vdfdsig,vdgdsig)
    calcu = 'D0F'

    !--------------------------------------------------------------------------------------
    !--- Si l'élément est ouvert
    if (ietatpg(ie,ipg)==1) then

        !----------------------------------------------------------------------------------
        !--- Sorties :
        vsig = 0.d0
        vnle = 0.D0
        vhep = 0.d0
        return

    endif

    !--------------------------------------------------------------------------------------
    !--- Calcul du multiplicateur plastique (algo semi-implicite)
    it = 0
    vsigi = vsigm
    dlam = 0.d0
    nconv = .true.

    if (vcrit>0.d0) then

        do while (nconv)

            dlam = dlam + vcrit/dot_product(vdfdsig,matmul(vh,vdgdsig))
            vsigm = vsigi - dlam*matmul(vh,vdgdsig)
            call interf_crit(iloi,vprel,vsigm,ie,ipg,vcrit,calcu)
            it = it + 1
            if (abs(vcrit) < 1.d-10) nconv = .false.
            if (it>50) then
                print*,ie,ipg,ietatpg(ie,ipg)
                stop 'elem_interf: non convergence dans le calcul des contraintes'
            endif
        enddo

    endif

    !--------------------------------------------------------------------------------------
    !--- Calcul du module elastoplastique
    if (imet == 1) then
        call elem_hooke(vh,nomtype(ktypel(ie)),vprel)
        call interf_modul(ie,ipg,vh,vprel)
        vhep = vh
    elseif (imet == 2) then
        if (ietatpg(ie,ipg)==1 .or. ietatpg(ie,ipg)==-1) then
            vhep = 0.0d0

        elseif (ietatpg(ie,ipg)==2 .or. ietatpg(ie,ipg)==-2 .or. ietatpg(ie,ipg)==3) then
            call elem_hooke(vh,nomtype(ktypel(ie)),vprel)
            vhep = vh - vecpro(matmul(vh,vdgdsig),matmul(vdfdsig,vh))/dot_product(vdfdsig,matmul(vh,vdgdsig))

        endif
    else
        vhep = vh
    endif

    vepspl = vepspl + dlam*vdgdsig

    !--------------------------------------------------------------------------------------
    !--- Sorties :
    vsig = vsigm
    vnle = vepspl

end subroutine interf_sig


!------------------------------------------------------------------------------------------------------
!> @author JL Tailhan
!> @authors J. Goncalvez (version 1.0 - 2010)
!
!> @brief
!> Calcul des criteres de rupture et des derivees pour l'element d'interface
!
!> @details
!> #### DESCRIPTION:
!>  sorties :  vcrit, vdfdsig, vdgdsig
!------------------------------------------------------------------------------------------------------
subroutine interf_crit(iloi,vprel,vsigm,ie,ipg,vcrit,calcu,vdfdsig,vdgdsig)

    !**************************************************************************************
    !* DEFINITIONS
    !**************************************************************************************

    !--------------------------------------------------------------------------------------
    !--- Modules generaux
    use variables, only :   dime, kprop0, ietatpg, &
    &                       iedng, pi, ktypel, nomtype
    use lib_elem, only :    elem_hooke

    implicit none

    !--------------------------------------------------------------------------------------
    !--- Variables in/out
    real*8, dimension(:), intent(in) :: vsigm, vprel
    integer, intent(in) :: ie, iloi, ipg
    character(len=5), intent(in), optional :: calcu
    real*8, dimension(:), intent(out), optional :: vdfdsig, vdgdsig
    real*8, intent(out) :: vcrit

    !--------------------------------------------------------------------------------------
    !--- Variables locales
    character(len=5) :: typel
    real*8 :: sigma, tau, tau1, tau2, RT, C, phi, psi
    real*8 :: vcrit1, vcrit2, vcrit2a, vcrit2b, vcrit3, vcrit4, kn, kt, kt1, kt2, A
    integer :: id, ipc
    logical :: ipremf, ideriv
    real*8, dimension(:,:), allocatable :: vh

    !**************************************************************************************
    !* CORPS DE PROCEDURE
    !**************************************************************************************

    !--------------------------------------------------------------------------------------
    !--- Test(s) de compatibilite des arguments d'entree
    ideriv = .false.
    typel = nomtype(ktypel(ie))

    if (calcu == 'D1F') then
        ideriv = .true.
        vdfdsig = 0.0d0
        vdgdsig = 0.0d0

        call elem_hooke(vh,typel,vprel)
        call interf_modul(ie,ipg,vh,vprel)
    endif

    !--------------------------------------------------------------------------------------
    !--- Recuperation des parametres de la loi
    id = 6
    if (dime == 3) id = id-1
    C  = vprel(id+1)
    phi = pi*vprel(id+2)/180.           ! Angle de frottement
    psi = pi*vprel(id+3)/180.           ! Angle de dilatance
    RT = vprel(id)                      ! RT

    !--------------------------------------------------------------------------------------
    !--- Recuperation des contraintes normale(s) et tangentielle(s)
    if (dime == 2) then
        sigma = vsigm(1)                  ! Contrainte normale
        tau   = vsigm(2)                  ! Contrainte tangente

    elseif (dime == 3) then
        sigma  = vsigm(1)                ! Contrainte normale
        tau1   = vsigm(2)                ! Contrainte tangente 1
        tau2   = vsigm(3)                ! Contrainte tangente 2
        tau = sqrt(tau1**2+tau2**2)
    endif

    !--------------------------------------------------------------------------------------
    !--- Gestion des etats de l'interface
    if (typel(1:4)=='EJQ4'.or.typel(1:4)=='EJQ6'.or.typel(1:4)=='EJT6') then
        ipc = 1
    else
        print*,'type',typel(1:4),'numero de groupe',kprop0(ie)
        print*,'===>interf_crit : element non implante'
    endif

    !--------------------------------------------------------------------------------------
    !--- Detection premiere fissuration
    ipremf = .false.
    if (all(ietatpg(ie,:)==0)) ipremf = .true.

    !--------------------------------------------------------------------------------------
    !--- Aiguillage en fonction de la loi de comportement
    select case(iloi)

        !----------------------------------------------------------------------------------
        !--- Lois complexes
        case(100,102,105,106,108,158)

            if (RT > (C/tan(phi))) RT = C/tan(phi)

            if ((ietatpg(ie,ipg)/=0).or.(ie==iedng)) then  !<---element fissure ou element dangereux
                RT=0.0d0
                C=0.0d0
            endif

            !------------------------------------------------------------------------------
            !--- Criteres de rupture de l'element d'interface
            !    dans le cas du calcul iteratif des contraintes
            if (calcu=='D0F') then
                if((ietatpg(ie,ipg)==1).or.(ietatpg(ie,ipg)<0)) then
                    vcrit = sqrt(sigma**2 + tau**2)
                elseif (ietatpg(ie,ipg)==2) then
                    vcrit = abs(tau) + sigma * tan(phi)
                endif

                return
            endif

            !------------------------------------------------------------------------------
            !--- dans le cas de la detection de la rupture
            vcrit1 = sigma - RT                               ! Tension cut-off
            vcrit2 = abs(tau) + sigma * tan(phi) - C          ! Mohr-Coulomb
            !vcrit3 = abs(tau) + sigma * tan(psi-pi/2.)
            vcrit3 = -sigma

            vcrit = -1.D0

            !------------------------------------------------------------------------------
            !--- Pour l'element dangereux ou un element deja fissure (ietatpg/=0)
            if ((iedng==ie) .or. (ietatpg(ie,ipg) /=0 )) then                 ! strategie 1
            !if ((iedng==ie) .or. (ipremf.eqv. .false.)) then                 ! strategie 2

                !--------------------------------------------------------------------------
                !--- Contrainte normale positive
                if ((vcrit3 < 0.d0) .or. (ietatpg(ie,ipg) < 0.d0)) then

                    !----------------------------------------------------------------------
                    !--- Critere en ouverture pure
                    if ((vcrit1> 0.d0) .or. (vcrit2> 0.d0)) then
                        vcrit = sqrt(sigma**2 + tau**2)

                        !if (ideriv .eqv. .true.) then
                        !if (ietatpg(ie,ipg) >= 0)  ietatpg(ie,ipg) = 1     ! strategies 2 4
                        if (ietatpg(ie,ipg) > 0)  ietatpg(ie,ipg) = 1       ! strategie 1
                        !if (iedng==ie) ietatpg(ie,:) = 1                   ! strategie 1

                        kn = vh(1,1)

                        if (dime == 2) then
                            kt = vh(2,2)
                            ! Dans le repere local...
                            vdfdsig = (/ sigma, tau /)
                            vdfdsig = vdfdsig/vcrit
                            vdgdsig = (/ (sigma/kn), (tau/kt) /)

                        elseif (dime == 3) then
                            kt1 = vh(2,2)
                            kt2 = vh(3,3)
                            vdfdsig = (/ sigma, tau1, tau2 /)
                            vdfdsig = vdfdsig/vcrit
                            vdgdsig = (/ (sigma/kn), tau1/kt1, tau2/kt2 /)
                        endif
                        !endif
                    endif

                    if (iedng==ie) ietatpg(ie,ipg) = 1

                else

                    !----------------------------------------------------------------------
                    !--- Critere en cisaillement
                    if ((vcrit2 > 0.d0) .or. (vcrit1 > 0.d0)) then
                        vcrit = abs(tau) + sigma * tan(phi)

                        !if (ideriv) then
                        ietatpg(ie,ipg) = 2                             ! strategies 1 2 3 4
                        if (iedng==ie) ietatpg(ie,:) = 2                ! strategie 1

                        ! dans le repere local...
                        if (dime == 2) then
                            vdfdsig = (/ tan(phi), tau/abs(tau) /)
                            vdgdsig = (/ tan(psi), tau/abs(tau) /)

                        elseif (dime == 3) then
                            vdfdsig = (/ tan(phi), tau1/abs(tau), tau2/abs(tau) /)
                            vdgdsig = (/ tan(psi), tau1/abs(tau), tau2/abs(tau) /)
                        endif
                        !endif
                    endif

                    if (iedng==ie) ietatpg(ie,ipg) = 2

                endif
            endif                                                          ! strategies 1 et 2

        !----------------------------------------------------------------------------------
        !--- Loi complexe : tension cut-off + Tresca puis Coulomb
        case(101)

            if (RT > (C/tan(phi))) RT = C/tan(phi)

            !psi = atan(RT/(C+1.d-20)) dans la loi 101 probabilisee psi est calcule dans distal
            !if ((ietatpg(ie,ipg) /=0 ).or.(ie==iedng)) then     ! Element fissure
            if ((ietatpg(ie,ipg) /=0 )) then     ! Element fissure
                RT=0.0d0
                !C=0.d0
            endif

            !------------------------------------------------------------------------------
            !--- Criteres de rupture de l'element d'interface
            !    dans le cas du calcul iteratif des contraintes
            if (calcu=='D0F')then
                if ((ietatpg(ie,ipg)==1).or.(ietatpg(ie,ipg)<0)) then
                    vcrit = sqrt(sigma**2 + tau**2)
                elseif (ietatpg(ie,ipg)==2) then
                    A = C*(1-(tan(psi-pi/2.)/tan(phi)))
                    vcrit4  = abs(tau) + sigma * tan(psi-pi/2.) - A
                    if (vcrit4 < 0.d0) then
                        vcrit = abs(tau) + sigma * tan(phi)
                    else
                        vcrit = abs(tau)  - C
                    endif
                endif

                return
            endif

            !------------------------------------------------------------------------------
            !--- dans le cas de la detection de la rupture
            vcrit1  = sigma - RT                          ! Tension cut-off
            vcrit2a = abs(tau)  - C                       ! Tresca
            vcrit2b = abs(tau) + sigma * tan(phi)         ! Mohr-Coulomb
            !vcrit3  = abs(tau) + sigma * tan(psi-pi/2.)
            vcrit3  = -sigma
            A = C*(1-(tan(psi-pi/2.)/tan(phi)))
            vcrit4  = vcrit3 - A

            vcrit = -1.D0

            if ((iedng==ie) .or. (ietatpg(ie,ipg) /=0 )) then       ! strateg 1
            !if ((iedng==ie) .or. (ietatpg(ie,ipg) /=0 ) .or. (.not. ipremf)) then       ! strateg 1
            !if ((iedng==ie) .or. ( ipremf .eqv. .false. )) then    ! strateg 2

                if ((vcrit3 < 0.d0) .or. (ietatpg(ie,ipg) < 0.d0)) then

                    !vcrit2=vcrit2a
                    !if (ietatpg(ie,ipg) /=0 ) vcrit2 = vcrit2b

                    !--------------------------------------------------------------------------
                    !--- Critere en ouverture pure
                    !if ((vcrit1> 0.d0) .or. (vcrit2> 0.d0) .or. (ietatpg(ie,ipg) < 0)) then
                    if ((vcrit1> 0.d0) .or. (ie==iedng) .or. (ietatpg(ie,ipg) < 0)) then
                        vcrit = sqrt(sigma**2 + tau**2)

                        !if (ideriv) then

                        !if (ietatpg(ie,ipg) >= 0)  ietatpg(ie,ipg) = 1  ! strateg 2 4
                        if (ietatpg(ie,ipg) > 0)  ietatpg(ie,ipg) = 1    ! strateg 1 -deb
                        if (iedng==ie) ietatpg(ie,ipg) = 1                 ! strateg 1 -fin

                        kn = vh(1,1)

                        if (dime == 2) then
                            kt = vh(2,2)
                            ! Dans le repere local...
                            vdfdsig = (/ sigma, tau /)
                            vdfdsig = vdfdsig/vcrit
                            vdgdsig = (/ (sigma/kn), (tau/kt) /)

                        elseif (dime == 3) then
                            kt1 = vh(2,2)
                            kt2 = vh(3,3)
                            vdfdsig = (/ sigma, tau1, tau2 /)
                            vdfdsig = vdfdsig/vcrit
                            vdgdsig = (/ (sigma/kn), tau1/kt1, tau2/kt2 /)
                        endif
                        !endif
                    endif

                elseif ((vcrit3 > 0.d0) .and. (vcrit4 < 0.d0)) then

                    !vcrit2=vcrit2a
                    !if (ietatpg(ie,ipg) /=0 ) vcrit2 = vcrit2b
                    vcrit2 = abs(tau) + sigma * tan(phi)

                    !--------------------------------------------------------------------------
                    !--- Critere en cisaillement
                    if ((vcrit2 > 0.d0) .or. (vcrit1 > 0.d0) .or. (ie==iedng)) then

                        vcrit = abs(tau) + sigma * tan(phi)

                        !if (ideriv) then

                        ietatpg(ie,ipg) = 2    ! strateg 1 2 3 4
                        !if (iedng==ie) ietatpg(ie,ipg) = 2     ! strateg 1

                        if (dime == 2) then
                               ! dans le repere local...
                               vdfdsig = (/ tan(phi), tau/abs(tau) /)
                               if (vcrit2==vcrit2a) vdfdsig = (/ 0d0, tau/abs(tau) /)
                               vdgdsig = (/ tan(psi), tau/abs(tau) /)

                        elseif (dime == 3) then
                               vdfdsig = (/ tan(phi), tau1/abs(tau), tau2/abs(tau) /)
                               if (vcrit2==vcrit2a) vdfdsig = (/ 0d0, tau1/abs(tau), tau2/abs(tau) /)
                               vdgdsig = (/ tan(psi), tau1/abs(tau), tau2/abs(tau) /)
                        endif
                        !endif
                    endif

                elseif ((vcrit3 > 0.d0) .and. (vcrit4 > 0.d0)) then

                    vcrit2=vcrit2a
                    !--------------------------------------------------------------------------
                    !--- Critere de Tresca
                    if ((vcrit2 > 0.d0) .or. (vcrit1 > 0.d0) .or. (ie==iedng)) then
                        vcrit = abs(tau) - C

                        !if (ideriv) then

                        ietatpg(ie,ipg) = 2    ! strateg 1 2 3 4
                        !if (iedng==ie) ietatpg(ie,ipg) = 2      ! strateg 1 -fin

                        if (dime == 2) then
                            ! dans le repere local...
                            vdfdsig = (/ 0.D0    , tau/abs(tau) /)
                            vdgdsig = (/ tan(psi), tau/abs(tau) /)

                        elseif (dime == 3) then
                            vdfdsig = (/ 0.D0    , tau1/abs(tau), tau2/abs(tau) /)
                            vdgdsig = (/ tan(psi), tau1/abs(tau), tau2/abs(tau) /)

                        endif
                        !endif
                    endif

                endif
            endif ! strateg 1 2

        !----------------------------------------------------------------------------------
        !--- Elastoplasticite Coulomb--Tresca
        case(103)

            psi = atan(RT/C)
            if (ietatpg(ie,ipg) /=0 ) then     ! Element fissure ou element referme
                RT=0.0d0
            endif

            !------------------------------------------------------------------------------
            !--- Criteres de rupture de l'element d'interface
            vcrit1  = sigma - RT                          ! Tension cut-off
            vcrit2a = abs(tau)  - C                       ! Tresca
            vcrit2b = abs(tau) + sigma * tan(phi)         ! Mohr-Coulomb
            vcrit3  = abs(tau) + sigma * tan(psi-pi/2.)
            A = C*(1-(tan(psi-pi/2.)/tan(phi)))
            vcrit4  = vcrit3 - A

            vcrit = -1.D0

            if ((iedng==ie) .or. (ietatpg(ie,ipg) /=0 )) then
                if ((vcrit3 < 0.d0) .or. (ietatpg(ie,ipg) < 0.d0)) then
                    vcrit2=vcrit2a
                    if (ietatpg(ie,ipg) /=0 ) vcrit2 = vcrit2b

                    !----------------------------------------------------------------------
                    !--- Critere en ouverture pure
                    if ((vcrit1> 0.d0) .or. (vcrit2> 0.d0)) then
                        vcrit = sqrt(sigma**2 + tau**2)

                        if (ideriv) then
                            if (ietatpg(ie,ipg) > 0)  ietatpg(ie,ipg) = 1

                            if ((iedng==ie) .and. (ipremf .eqv. .true.)) then
                                ietatpg(ie,:) = 1
                            endif
                            kn = vh(1,1)

                            if (dime == 2) then
                                kt = vh(2,2)
                                ! Dans le repere local...
                                vdfdsig = (/ sigma, tau /)
                                vdfdsig = vdfdsig/vcrit
                                vdgdsig = (/ (sigma/kn), (tau/kt) /)

                            elseif (dime == 3) then
                                kt1 = vh(2,2)
                                kt2 = vh(3,3)
                                vdfdsig = (/ sigma, tau1, tau2 /)
                                vdfdsig = vdfdsig/vcrit
                                vdgdsig = (/ (sigma/kn), tau1/kt1, tau2/kt2 /)
                            endif
                        endif
                    endif

                elseif ((vcrit3 > 0.d0) .and. (vcrit4 < 0.d0)) then

                    if ((ietatpg(ie,ipg)/=0) .and. (ietatpg(ie,ipg)/=3)) then

                        vcrit2 = vcrit2b

                        !------------------------------------------------------------------
                        !--- Critere en cisaillement
                        if ((vcrit2 > 0.d0) .or. (vcrit1 > 0.d0)) then
                            vcrit = abs(tau) + sigma * tan(phi)

                            if (ideriv .eqv. .true.) then
                                ietatpg(ie,ipg) = 2

                                if ((iedng==ie) .and. (ipremf .eqv. .true.)) then
                                    ietatpg(ie,:) = 2
                                endif

                                if (dime == 2) then
                                    ! dans le repere local...
                                    vdfdsig = (/ tan(phi), tau/abs(tau) /)
                                    vdgdsig = (/ tan(psi), tau/abs(tau) /)

                                elseif (dime == 3) then
                                    vdfdsig = (/tan(phi), tau1/abs(tau), tau2/abs(tau) /)
                                    vdgdsig = (/tan(psi), tau1/abs(tau), tau2/abs(tau) /)
                                endif
                            endif
                        endif

                    !----------------------------------------------------------------------
                    !--- Elastoplastique
                    !    Si l'element est vierge et en elastoplastique
                    elseif ((ietatpg(ie,ipg)==0) .or. (ietatpg(ie,ipg)==3)) then
                        vcrit2 = vcrit2a

                        if ((vcrit2 > 0.d0) .or. (vcrit1 > 0.d0)) then
                            vcrit = abs(tau) - C

                            if (ideriv .eqv. .true.) then
                                ietatpg(ie,ipg) = 3
                                if ((iedng==ie) .and. (ipremf .eqv. .true.)) then
                                    ietatpg(ie,:) = 3
                                endif

                                if (dime == 2) then
                                    vdfdsig = (/0.d0 ,  tau/abs(tau)/)
                                    vdgdsig = (/0.d0 ,   1.d0/)

                                elseif (dime == 3) then
                                    vdfdsig = (/0.d0 ,   tau1/abs(tau) , tau2/abs(tau)/)
                                    vdgdsig= (/0.d0 ,   1.d0    ,  1.d0/)
                                endif

                            endif
                        endif
                    endif

                !--------------------------------------------------
                else
                    vcrit2 = vcrit2a
                    ! Critere de Tresca
                    if ((vcrit2 > 0.d0) .or. (vcrit1 > 0.d0)) then
                        vcrit = abs(tau) - C

                        if (ideriv .eqv. .true.) then

                            if ((ietatpg(ie,ipg)/=0) .and. (ietatpg(ie,ipg)/=3)) then

                                ietatpg(ie,ipg) = 2
                                if ((iedng==ie) .and. (ipremf .eqv. .true.)) then
                                    ietatpg(ie,:) = 2
                                endif

                                if (dime == 2) then
                                    ! dans le repere local...
                                    vdfdsig = (/ 0.D0    , tau/abs(tau) /)
                                    vdgdsig = (/ tan(psi), tau/abs(tau) /)

                                elseif (dime == 3) then
                                    vdfdsig = (/0.D0    , tau1/abs(tau), tau2/abs(tau) /)
                                    vdgdsig = (/tan(psi), tau1/abs(tau), tau2/abs(tau) /)
                                endif

                            ! Elastoplastique
                            elseif ((ietatpg(ie,ipg)==0) .or. (ietatpg(ie,ipg)==3)) then
                                ietatpg(ie,ipg) = 3
                                if ((iedng==ie) .and. (ipremf .eqv. .true.)) then
                                    ietatpg(ie,:) = 3
                                endif

                                if (dime == 2) then
                                    ! dans le repere local...
                                    vdfdsig = (/0.d0 ,   tau/abs(tau) /)
                                    vdgdsig = (/0.d0 ,   1.d0/)

                                elseif (dime == 3) then
                                    vdfdsig = (/0.d0 ,   tau1/abs(tau) , tau2/abs(tau)/)
                                    vdgdsig = (/0.d0 ,   1.d0    ,  1.d0/)
                                endif

                            endif
                        endif
                    endif
                endif
            endif

        case default
            print*, 'interf_crit : loi non implantee pour l''element d''interface'
            stop
    end select

end subroutine interf_crit


!------------------------------------------------------------------------------------------------------
!> @author JL Tailhan
!> @authors J. Goncalvez (version 1.0 - 2010)
!
!> @brief
!> Calcul des contraintes pour l'élément de contact endommageable
!
!> @details
!> #### DESCRIPTION:
!>  A faire
!------------------------------------------------------------------------------------------------------
subroutine interf_sig_endo(vhep,vsig,vnle,vin,vh,vb,vnl0,vprel,vdle,iloi,ie,ipg)

    !**************************************************************************************
    !* DEFINITIONS
    !**************************************************************************************

    !--------------------------------------------------------------------------------------
    !--- Modules generaux
    use variables, only :   dime, ietatpg, ktypel, nomtype, endo, &
    &                       iedng
    use lib_elem, only :    elem_hooke

    implicit none

    !--------------------------------------------------------------------------------------
    !--- Variables in/out
    real*8, dimension(:), intent(in) :: vnl0, vdle
    real*8, dimension(:), intent(in) :: vprel
    real*8, dimension(:,:), intent(in) :: vh
    real*8, dimension(:,:), intent(in) :: vb
    integer, intent(in) :: ie, ipg, iloi
    real*8, dimension(:), intent(inout) :: vsig, vnle, vin
    real*8, dimension(:,:), intent(out) :: vhep

    !--------------------------------------------------------------------------------------
    !--- Variables locales
    character(len=5) :: typel
    real*8, dimension(size(vb,1)) :: vepspl, vsigm, deprel
    real*8, dimension(size(vh,1),size(vh,2)) :: vhpen, vhpetit
    real*8 :: depreln, deprelt, depreln0, deprelt0
    real*8 :: RT, C, D, D0, penal, alpha , Dn, Dt
    integer :: id, ipc

    !--- Variables pour modèle 102
    real*8 :: depcrin, depcrit

    !--- Variables pour modèles 105, 106 et 107
    !real*8, dimension(:,:), allocatable :: vhendo
    real*8, dimension(dime,dime) :: vhendo
    real*8 :: sigma, tau, tau1, tau2
    real*8 :: var, var0, varc, vars, F_var
    real*8 :: Dini, vari, deprelni, penali, deprelti, F_vari, dF_vari
    real*8 :: W, beta, xc, xp, xi, x0, x
    real*8 :: alphar, coef1, coef2
    integer :: imeth
    logical :: vcrit

    !--- Variables pour la loi 108
    real*8, dimension(dime,dime) :: vhep2
    real*8 :: Wc, Dini2, var02, vars2, varc2, D2, var2

    !**************************************************************************************
    !* CORPS DE PROCEDURE
    !**************************************************************************************

    !--------------------------------------------------------------------------------------
    !--- Recuperation des grandeurs de non-linearite
    vepspl = vnl0(1:size(vb,1))            ! deformations anelastiques initiales
    vhep = vh
    vhpen = 0.D0

    !--------------------------------------------------------------------------------------
    !--- Recuperation du point d'integration central en fonction du type d'element
    typel = nomtype(ktypel(ie))
    if ((typel(1:4)=='EJQ4').or.(typel(1:4)=='EJT6').or.(typel(1:4)=='EJQ6')) then
        ipc = 1
    else
        stop 'interf_crit : element non implante'
    endif

    !--------------------------------------------------------------------------------------
    !--- Recuperation des parametres de la loi
    id = 6
    if (dime == 3) id = id-1
    C  = vprel(id+1)
    RT = vprel(id)

    !--------------------------------------------------------------------------------------
    !--- Calcul des deplacements relatifs
    if (dime==2) then
        deprel = matmul(vb,vdle)
        depreln = deprel(1)                          ! Deplacement normal
        deprelt = deprel(2)                          ! Deplacement tangentiel

    elseif (dime==3) then
        deprel = matmul(vb,vdle)
        depreln = deprel(1)                          ! Deplacement normal
        deprelt = sqrt(deprel(2)**2 + deprel(3)**2)  ! Deplacement tangentiel

    endif

    penal = 0.d0
    if (depreln <= 0.d0) penal = 1.d0

    !--------------------------------------------------------------------------------------
    !--- Aiguillage en fonction de la moi de comportement
    select case(iloi)

        !----------------------------------------------------------------------------------
        !--- Loi d'endommagement (interface acier/beton)
        case(102)

            if (depreln < 0.d0) depreln = 0.d0

            !------------------------------------------------------------------------------
            !--- Deplacements relatifs critiques en traction et en cisaillement
            depcrin = vprel(id+4)
            depcrit = vprel(id+5)

            Dn = 0.d0
            Dt = 0.d0

            !------------------------------------------------------------------------------
            !--- Seuils d'endommagement
            depreln0 = RT/vh(1,1)
            deprelt0 = C/vh(2,2)
            if (depcrin < depreln0) depcrin=depreln0
            if (depcrit < deprelt0) depcrit=deprelt0

            !------------------------------------------------------------------------------
            !--- Calcul de l'endommagement
            if ((depreln >= depreln0) .and. (depreln < depcrin)) Dn = 1.d0 - depreln0/depreln
            if ((abs(deprelt) >= deprelt0) .and. (abs(deprelt) < depcrit)) Dt = 1.d0 - deprelt0/abs(deprelt)

            if (iedng==ie .or. ietatpg(ie,ipg) == 1 .or. ietatpg(ie,ipg) == 2) then
                if (abs(deprelt) >= depcrit) Dt = 1.d0
                if (depreln >= depcrin) Dn = 1.d0
            endif

            alpha = 0.d0
            D  = alpha * Dn + (1.-alpha) * Dt
            D0 = endo(ie,ipg)
            D  = max(D,D0)
            D  = min(1.D0,D)

            !if ((ipg==ipc) .and. (all(ietatpg(ie,:)==0)) .and. (D > 0)) ietatpg(ie,:) = 3
            !if ((ipg==ipc) .and. (all(ietatpg(ie,:)==3)) .and. (Dt == 1.d0)) ietatpg(ie,:) = 2
            !if ((ipg==ipc) .and. (all(ietatpg(ie,:)==3)) .and. (Dn == 1.d0)) ietatpg(ie,:) = 1

            if (dime==2) then
                if ((ipg==ipc).and.(ietatpg(ie,1)==0).and.(ietatpg(ie,2)==0).and.(ietatpg(ie,3)==0)) then
                    if (D > 0) ietatpg(ie,:) = 3
                endif

                if ((ipg==ipc).and.(ietatpg(ie,1)==3).and.(ietatpg(ie,2)==3).and.(ietatpg(ie,3)==3)) then
                    if (Dt == 1.d0) ietatpg(ie,:) = 2
                endif

                if ((ipg==ipc).and.(ietatpg(ie,1)==3).and.(ietatpg(ie,2)==3).and.(ietatpg(ie,3)==3)) then
                    if (Dn == 1.d0) ietatpg(ie,:) = 1
                endif

            elseif (dime==3) then

                if ((ipg==ipc) .and. all(ietatpg(ie,:)==0)) then
                    if (D > 0) ietatpg(ie,:) = 3
                endif

                if ((ipg==ipc) .and. all(ietatpg(ie,:)==3)) then
                    if (Dt == 1.d0) ietatpg(ie,:) = 2
                endif

                if ((ipg==ipc) .and. all(ietatpg(ie,:)==3)) then
                    if (Dn == 1.d0) ietatpg(ie,:) = 1
                endif

            endif

            if (ietatpg(ie,ipg)==0) D = 0.d0
            if ((ietatpg(ie,ipg)==1) .or. (ietatpg(ie,ipg)==2)) D = 1.d0

            if (dime==2) then
                vhpen(1,:) = D*penal*(/vh(1,1), 0.d0/)
                vhpen(2,:) = D*penal*(/0.d0   , 0.d0/)
            elseif (dime==3) then
                vhpen(1,:) = D*penal*(/vh(1,1), 0.d0, 0.d0/)
                vhpen(2,:) = D*penal*(/0.d0   , 0.d0, 0.d0/)
                vhpen(3,:) = D*penal*(/0.d0   , 0.d0, 0.d0/)
            endif

            vhpetit = vh/1.e15
            vhep = (1-D)*vh + vhpen + vhpetit

            !------------------------------------------------------------------------------
            !--- Calcul des contraintes
            vsigm = matmul(vhep,matmul(vb,vdle))
            vepspl = (-D*matmul(vb,vdle))

        !----------------------------------------------------------------------------------
        !--- Loi d'endommagement simple (beton)
        case(105)

            !------------------------------------------------------------------------------
            !--- Initialisations
            Dini = vin(1)   ! endommagement initial
            var0 = vin(2)   ! seuil initial
            vars = vin(3)   ! seuil courant
            deprelni = vin(4)
            deprelti = vin(5)

            D = Dini
            penali = 0.d0
            if (deprelni <= 0.d0) penali = 1.d0

            !------------------------------------------------------------------------------
            !--- Parametres du modele
            W      = vprel(id+4)   ! energie

            !------------------------------------------------------------------------------
            !--- Methode de calcul de la rigidite (imeth=1 secante, imeth=2 tangente)
            imeth = 2

            !------------------------------------------------------------------------------
            !--- Preparation pour calcul du module tangent
            vhpetit = vh/1.d50

            !------------------------------------------------------------------------------
            !--- Test d'element endommage
            vcrit = .false.

            if (ie==iedng) then

                !--------------------------------------------------------------------------
                !--- Cas de l element le plus dangereux
                vcrit = .true.
                ! NB : pour l'element le plus dangeureux, on favorise ici la rupture des points de Gauss
                !      autres que celui du centre, meme si le critere n'est pas verifie en ces points.

                !--------------------------------------------------------------------------
                !--- Calcul de l'etat de contrainte
                if (dime==2) then
                    tau1  = vh(2,2)*deprel(2)
                    tau   = sqrt(tau1**2)            ! Contrainte tangente elastique ( abs(tau) )
                    sigma = vh(1,1)*deprel(1)        ! Contrainte normale elastique

                elseif (dime==3) then
                    tau1  = vh(2,2)*deprel(2)
                    tau2  = vh(3,3)*deprel(3)
                    tau   = sqrt(tau1**2 + tau2**2) ! Contrainte tangente elastique ( abs(tau) )
                    sigma = vh(1,1)*depreln         ! Contrainte normale elastique
                endif

                !--------------------------------------------------------------------------
                !--- Mecanisme de rupture en ouverture pure
                if ((tau/C - sigma/RT) <= 0.D0) then
                    depreln0 = RT/vh(1,1)
                    !===> DEB AJOUT JLT
                    if(sigma<RT) then
                        depreln0 = depreln !sigma/vh(1,1)
                    endif
                    !<=== FIN AJOUT JLT
                    deprelt0 = depreln0/depreln*deprelt

                !--------------------------------------------------------------------------
                !--- Mecanisme de rupture en cisaillement
                elseif ((tau/C - sigma/RT) > 0.D0) then
                    deprelt0 = C/vh(2,2)
                    !===> DEB AJOUT JLT
                    if (tau < C) then
                        deprelt0 = deprelt !tau/vh(2,2)
                    endif
                    !<=== FIN AJOUT JLT
                    depreln0 = deprelt0/deprelt*depreln

                endif

                ietatpg(ie,:)=3 ! Declaration de l'element endommage

                var0 = sqrt(depreln0**2 + deprelt0**2) + 1.D-20
                if (depreln0 <= 0.d0) var0 = sqrt(deprelt0**2) + 1.D-20

                vars = var0
                varc = var0 + 2.*W/(vh(1,1)*var0)

                deprelni = depreln0
                deprelti = deprelt0

            else
                !--------------------------------------------------------------------------
                !--- Cas des elements deja endommages
                !if (ietatpg(ie,ipg)==3) then
                if (ietatpg(ie,ipg)/=0) then
                    vcrit=.true.
                    varc = var0 + 2.*W/(vh(1,1)*var0)
                endif
            endif

            !------------------------------------------------------------------------------
            !--- Parametres complementaires du modele
            x0 = var0
            xc = varc
            if (W == 0.d0) xc = 1.001d0 * x0

            !------------------------------------------------------------------------------
            !--- Si element endommage
            if (vcrit) then

                !--------------------------------------------------------------------------
                !--- variable pilotant l'endommagement
                var = dsqrt((1.-penal)*depreln**2 + deprelt**2)         ! Deplacement relatif total

                !--------------------------------------------------------------------------
                !--- Cas de l'evolution de l'endommagement
                if (var>=vars) then

                    !----------------------------------------------------------------------
                    !--- Actualisation du seuil d'endommagement
                    vars = var

                    !----------------------------------------------------------------------
                    !--- loi d'evolution d'endommagement
                    x = var
                    F_var = 0.d0
                    if ((x>=x0).and.(x<xc)) then
                        F_var  = 1.d0 - (xc-x)/(xc-x0)
                    else
                        F_var  = 1.d0
                    endif

                    !----------------------------------------------------------------------
                    !--- Calcul de l'endommagement
                    D = 1.d0 - var0/var * (1.d0 - F_var)

                    if (D<0.d0) D=0.d0
                    D = max(Dini,D)
                    D = min(D,1.d0)

                    !----------------------------------------------------------------------
                    !--- Calcul du module vhep (imet==1 et imet==2)
                    if (imeth == 1) then

                        !------------------------------------------------------------------
                        !--- Calcul du module secant (endommage)
                        vhep = (1.-D)*vh

                    elseif (imeth == 2) then

                        !------------------------------------------------------------------
                        !--- Calcul du module tangent
                        vari = dsqrt((1.-penali)*deprelni**2 + deprelti**2) + 1.D-20

                        x = vari
                        F_vari = 0.d0
                        dF_vari = 0.d0
                        if ((x>=x0).and.(x<xc)) then
                            F_vari  = 1.d0 - (xc-x)/(xc-x0)
                            dF_vari = 1.d0 / (xc-x0)
                        else
                            F_vari  = 1.d0
                            dF_vari  = 0.d0
                        endif

                        if (dime==2) then
                            coef1 = (deprelni/vari)**2
                            coef2 = (deprelni*deprelti)/(vari**2)
                            if (coef1>.99d0) then
                                coef1 = 1.d0
                                coef2 = 0.d0
                            endif
                            if (coef1<.01d0) then
                                coef1 = 0.d0
                                coef2 = 0.d0
                            endif
                            vhep(1,1) =  vh(1,1) * ((1.-Dini)- var0*((1.-F_vari)/vari - dF_vari)*coef1)
                            vhep(2,2) =  vh(2,2) * ((1.-Dini)- var0*((1.-F_vari)/vari - dF_vari)*(1.-coef1))
                            vhep(1,2) =  vh(1,2) * (-var0*((1.-F_vari)/vari - dF_vari)*coef2)
                            vhep(1,2) =  vh(2,1) * (-var0*((1.-F_vari)/vari - dF_vari)*coef2)

                        else if (dime==3) then
                            stop 'interf_sig_endo (105) : cas 3D non encore implante'
                        endif
                    endif

                    if ((var>varc) .or. (D==1.d0)) then
                        if (var>varc) D=1.d0
                        !if (ipg==ipc) then
                            ietatpg(ie,ipg) = 1           ! Element ouvert
                            if (depreln <= 0.d0) then
                                ietatpg(ie,ipg) = 2       ! Element referme
                            endif
                        !endif
                    endif

                !--------------------------------------------------------------------------
                !--- Cas ou l'endommagement n'evolue pas
                else
                    vhep = (1.-D) * vh

                endif

                !--------------------------------------------------------------------------
                !--- Cas particulier du traitement de la refermeture de fissure
                !    (on referme sur un comportement elastique partiellement endommage)
                if (depreln <= 0.d0) then
                    vhep(1,1) = vh(1,1)
                endif

            endif

            !------------------------------------------------------------------------------
            !--- Calcul de vhep
            vhep = vhep + vhpetit

            !------------------------------------------------------------------------------
            !--- Calcul des contraintes
            vsigm(1) = (1.-D*(1.-penal)) * vh(1,1) * deprel(1)
            vsigm(2) = (1.-D) * vh(2,2) * deprel(2)
            if (dime==3) vsigm(3) = (1.-D) * vh(3,3) * deprel(3)

            !------------------------------------------------------------------------------
            !--- Sorties : variables internes du modèle
            vin(1) = D
            vin(2) = var0
            vin(3) = vars
            vin(4) = depreln
            vin(5) = deprelt


        !----------------------------------------------------------------------------------
        !--- Loi d'endommagement (BRF)
        case(106)

            !------------------------------------------------------------------------------
            !--- Initialisations
            Dini = vin(1)   ! endommagement initial
            var0 = vin(2)   ! seuil initial
            vars = vin(3)   ! seuil courant
            deprelni = vin(4)
            deprelti = vin(5)

            D = Dini
            penali = 0.d0
            if (deprelni <= 0.d0) penali = 1.d0

            !------------------------------------------------------------------------------
            !--- Parametres du modele
            W      = vprel(id+4)   ! energie totale post fissuration
            xp     = vprel(id+5)   ! seuil d’endommagement du pontage
            xc     = vprel(id+6)   ! depl relat critique pour le pontage
            alpha  = vprel(id+7)   ! rapport KT/KN pour l'ouverture
            alphar = vprel(id+8)   ! rapport KT/KN pour la refermeture

            if (xc <= xp) stop 'element_interface.f : loi 106, Erreur de donnees (xc<xp)'
            beta   = 2.*(W)/(xp*xc)/vh(1,1)  !xp constant, beta variable

            !------------------------------------------------------------------------------
            !--- Methode de calcul de la rigidite (imeth=1 secante, imeth=2 tangente)
            imeth = 1

            !------------------------------------------------------------------------------
            !--- Preparation pour calcul du module tangent
            !call elem_hooke(vhendo,nomtype(ktypel(ie)),vprel)
            vhendo = vh
            vhpetit = vh/1.d50

            !------------------------------------------------------------------------------
            !--- Test d'element endommage
            vcrit = .false.

            !------------------------------------------------------------------------------
            !--- Cas de l element le plus dangereux
            if (ie==iedng) then

                vcrit = .true.
                ! NB : pour l'element le plus dangeureux, on favorise ici la rupture des points de Gauss
                !      autres que celui du centre, meme si le critere n'est pas verifie en ces points.

                !--------------------------------------------------------------------------
                !--- Calcul de l'etat de contrainte (alpha tient compte de Kt~Kn)
                if (dime==2) then
                    tau1  = vh(2,2)*deprel(2)
                    tau   = sqrt(tau1**2)            ! Contrainte tangente elastique ( abs(tau) )
                    sigma = vh(1,1)*deprel(1)        ! Contrainte normale elastique

                elseif (dime==3) then
                    tau1  = vh(2,2)*deprel(2)
                    tau2  = vh(3,3)*deprel(3)
                    tau   = sqrt(tau1**2 + tau2**2) ! Contrainte tangente elastique ( abs(tau) )
                    sigma = vh(1,1)*depreln         ! Contrainte normale elastique
                endif

                !--------------------------------------------------------------------------
                !--- Mecanisme de rupture en ouverture pure
                if ((tau/C - sigma/RT) <= 0.D0) then
                    depreln0 = RT/vh(1,1)
                    !===> DEB AJOUT JLT
                    if(sigma<RT) then
                        depreln0 = depreln !sigma/vh(1,1)
                    endif
                    !<=== FIN AJOUT JLT
                    deprelt0 = depreln0/depreln*deprelt

                !--------------------------------------------------------------------------
                !--- Mecanisme de rupture en cisaillement
                elseif ((tau/C - sigma/RT) > 0.D0) then
                    deprelt0 = C/vh(2,2)
                    !===> DEB AJOUT JLT
                    if (tau < C) then
                        deprelt0 = deprelt !tau/vh(2,2)
                    endif
                    !<=== FIN AJOUT JLT
                    depreln0 = deprelt0/deprelt*depreln

                endif

                ietatpg(ie,:)=3 ! Declaration de l'element endommage
                var0 = xp
                vars = var0
                varc = xc

                deprelni = depreln0
                deprelti = deprelt0

            else

                !--------------------------------------------------------------------------
                !--- Cas des elements deja endommages
                if (ietatpg(ie,ipg)/=0) then
                    vcrit=.true.
                    varc = xc
                endif
            endif

            !------------------------------------------------------------------------------
            !--- Parametres complementaires du modele
            x0 = var0

            !------------------------------------------------------------------------------
            !--- Si element endommage
            if (vcrit) then

                !--------------------------------------------------------------------------
                !--- Complements pour calcul du module tangent
                vhendo = beta*vh
                vhendo(2,2) = alpha*vhendo(1,1)
                if (dime==3) vhendo(3,3) = alpha*vhendo(1,1)
                !vhendo(1,1) = beta*vhendo(1,1)
                !vhendo(2,2) = beta*alpha*vhendo(2,2)
                !if (dime==3) vhendo(3,3) = beta*alpha*vhendo(3,3)

                !--------------------------------------------------------------------------
                !--- variable pilotant l'endommagement
                var = (1.-penal)*dsqrt(depreln**2 + deprelt**2)

                !--------------------------------------------------------------------------
                !--- Cas de l'evolution de l'endommagement
                if (var>=vars) then

                    !----------------------------------------------------------------------
                    !--- Actualisation du seuil d'endommagement
                    vars = var

                    !----------------------------------------------------------------------
                    !--- loi d'evolution d'endommagement
                    x = var
                    F_var = 0.d0
                    if ((x>=x0).and.(x<xc)) then
                        F_var  = 1.d0 - (xc-x)/(xc-x0)
                    else
                        F_var  = 1.d0
                    endif

                    !----------------------------------------------------------------------
                    !--- Calcul de l'endommagement
                    D = 1.d0 - var0/var * (1.d0 - F_var)

                    if (D<0.d0) D=0.d0
                    D = max(Dini,D)
                    D = min(D,1.d0)

                    !----------------------------------------------------------------------
                    !--- Calcul du module vhep (imet==1 et imet==2)
                    if (imeth == 1) then
                        !----- Calcul du module secant (endommage)
                        vhep = (1.-D)*vhendo

                    elseif (imeth == 2) then
                        !----- Calcul du module tangent
                        vari = (1.-penali)*dsqrt(deprelni**2 + deprelti**2) + 1.D-20

                        x = vari
                        F_vari = 0.d0
                        dF_vari = 0.d0
                        if ((x>=x0).and.(x<xc)) then
                            F_vari  = 1.d0 - (xc-x)/(xc-x0)
                            dF_vari = 1.d0/(xc-x0)
                        else
                            F_vari  = 1.d0
                            dF_vari  = 0.d0
                        endif

                        if (dime==2) then
                            coef1 = (deprelni/vari)**2
                            coef2 = (deprelni*deprelti)/(vari**2)
                            if (coef1>.99d0) then
                                coef1 = 1.d0
                                coef2 = 0.d0
                            endif
                            if (coef1<.01d0) then
                                coef1 = 0.d0
                                coef2 = 0.d0
                            endif
                            vhep(1,1) =  vhendo(1,1) * ((1.-Dini)- var0*((1.-F_vari)/vari - dF_vari)*coef1)
                            vhep(2,2) =  vhendo(2,2) * ((1.-Dini)- var0*((1.-F_vari)/vari - dF_vari)*(1.-coef1))
                            vhep(1,2) =  vhendo(1,2) * (-var0*((1.-F_vari)/vari - dF_vari)*coef2)
                            vhep(1,2) =  vhendo(2,1) * (-var0*((1.-F_vari)/vari - dF_vari)*coef2)

                        else if (dime==3) then
                            stop 'interf_sig_endo (108) : cas 3D non encore implante'
                        endif
                    endif

                    if ((var>varc) .or. (D==1.d0)) then
                        if (var>varc) D=1.d0
                        !if (ipg==ipc) then
                            ietatpg(ie,ipg) = 1           ! Element ouvert
                            if (depreln <= 0.d0) then
                                ietatpg(ie,ipg) = 2       ! Element referme
                            endif
                        !endif
                    endif


                !--------------------------------------------------------------------------
                !--- Cas ou l'endommagement n'evolue pas
                else
                    vhep = (1.-D)*vhendo

                endif

                !--------------------------------------------------------------------------
                !--- Cas particulier du traitement de la refermeture de fissure
                !    (on referme sur un comportement elastique partiellement endommage)
                !    ( var = 0.d0   et D = Dini)
                if (depreln <= 0.d0) then
                    !----- Correction sur vhendo
                    vhendo(1,1) = vh(1,1)
                    vhendo(2,2) = (alphar/alpha)*vhendo(2,2)
                    if (dime==3) vhendo(3,3) = (alphar/alpha)*vhendo(3,3)
                    !----- Correction sur vhep
                    vhep = (1.-D)*vhendo
                    vhep(1,1) = vh(1,1)
                endif

            endif

            !------------------------------------------------------------------------------
            !--- Calcul de vhep
            vhep = vhep + vhpetit

            !------------------------------------------------------------------------------
            !--- Calcul des contraintes
            vsigm(1) = (1.-(1.-penal)*D) * vhendo(1,1) * deprel(1)
            vsigm(2) = (1.-D) * vhendo(2,2) * deprel(2)
            if (dime==3) vsigm(3) = (1.-D) * vhendo(3,3) * deprel(3)

            !------------------------------------------------------------------------------
            !--- Sorties : variables internes du modèle
            vin(1) = D
            vin(2) = var0
            vin(3) = vars
            vin(4) = depreln
            vin(5) = deprelt

            !------------------------------------------------------------------------------
            !--- Desallocation memoire
            !deallocate(vhendo)



        !----------------------------------------------------------------------------------
        !--- Loi d'interface simple (BA fissuré - acier elastique)
        case(107)

            !------------------------------------------------------------------------------
            !--- Initialisations
            D = 0.D0

            !------------------------------------------------------------------------------
            !--- Parametres du modele
            beta = vprel(id+4)   ! rapport KT/KN pour l'ouverture

            !------------------------------------------------------------------------------
            !--- Preparation pour calcul du module tangent
            !call elem_hooke(vhendo,nomtype(ktypel(ie)),vprel)
            vhendo = vh
            vhpetit = vh/1.d50

            !------------------------------------------------------------------------------
            !--- Test d'element endommage
            vcrit = .false.

            !------------------------------------------------------------------------------
            !--- Cas de l'element le plus dangereux
            if (ie==iedng) then
                vcrit=.true.
                ietatpg(ie,:)=3 ! Declaration de l'element endommage

            !------------------------------------------------------------------------------
            !--- Cas des autres elements
            else

                !--------------------------------------------------------------------------
                !--- Cas des elements deja endommages
                !if (ietatpg(ie,ipg)==3) then
                if (ietatpg(ie,ipg)/=0) then
                    vcrit=.true.
                endif
            endif

            !------------------------------------------------------------------------------
            !--- Calcul du module vhep
            if (vcrit) then

                vhep = beta*vhendo

                !--------------------------------------------------------------------------
                !--- Cas particulier du traitement de la refermeture de fissure
                !    (on referme sur un comportement elastique partiellement endommage)
                if (depreln <= 0.d0) then
                    vhep(1,1) = vh(1,1)
                endif
            endif

            vhep = vhep + vhpetit

            !------------------------------------------------------------------------------
            !--- Calcul des contraintes
            vsigm(1) = vhep(1,1)*deprel(1)
            vsigm(2) = vhep(2,2)*deprel(2)
            if (dime==3) vsigm(3) = vhep(3,3) * deprel(3)

            !------------------------------------------------------------------------------
            !--- Sorties : variables internes du modèle
            vin = 0.D0

            !------------------------------------------------------------------------------
            !--- Desallocation memoire
            !deallocate(vhendo)


        !----------------------------------------------------------------------------------
        !--- Loi d'endommagement complexe (beton)
        case(108)
            print*,'attention, il y a un pb avec la variable xi'

            !------------------------------------------------------------------------------
            !--- Initialisations
            Dini  = vin(1)   ! endommagement initial
            var0  = vin(2)   ! seuil initial
            vars  = vin(3)   ! seuil courant
            Dini2 = vin(4)   ! endommagement initial
            vars2 = vin(5)   ! seuil courant
            deprelni = vin(6)
            deprelti = vin(7)

            D = Dini
            D2 = Dini2
            penali = 0.d0
            if (deprelni <= 0.d0) penali = 1.d0

            !------------------------------------------------------------------------------
            !--- Parametres du modele
            W      = vprel(id+4)   ! energie totale post fissuration
            Wc     = vprel(id+5)   ! energie post fissuration matrice
            xp     = vprel(id+6)   ! seuil d’endommagement du pontage
            xc     = vprel(id+7)   ! depl relat critique pour le pontage
            alpha  = vprel(id+8)   ! rapport KT/KN pour l'ouverture
            alphar = vprel(id+9)   ! rapport KT/KN pour la refermeture

            if (xc <= xp) stop 'element_interface.f : loi 108, Erreur de donnees (xc<xp)'
            if (W<=Wc) W=Wc
            beta   = 2.*(W-Wc)/(xp*xc)/vh(1,1)  !xp constant, beta variable

            !------------------------------------------------------------------------------
            !--- Methode de calcul de la rigidite (imeth=1 secante, imeth=2 tangente)
            imeth = 2

            !------------------------------------------------------------------------------
            !--- Preparation pour calcul du module tangent
            vhendo = beta*vh
            vhendo(2,2) = alpha*vhendo(1,1)
            if (dime==3) vhendo(3,3) = alpha*vhendo(1,1)
            vhep2 = vhendo
            vhep  = vh
            vhpetit = vh/1.d50

            !------------------------------------------------------------------------------
            !--- Test d'element endommage
            vcrit = .false.

            !------------------------------------------------------------------------------
            !--- Cas de l element le plus dangereux
            if (ie==iedng) then
                vcrit = .true.
                ! NB : pour l'element le plus dangeureux, on favorise ici la rupture des points de Gauss
                !      autres que celui du centre, meme si le critere n'est pas verifie en ces points.

                !--------------------------------------------------------------------------
                !--- Calcul de l'etat de contrainte (alpha tient compte de Kt~Kn)
                if (dime==2) then
                    tau1  = vh(2,2)*deprel(2)
                    tau   = sqrt(tau1**2)            ! Contrainte tangente elastique ( abs(tau) )
                    sigma = vh(1,1)*deprel(1)        ! Contrainte normale elastique

                elseif (dime==3) then
                    tau1  = vh(2,2)*deprel(2)
                    tau2  = vh(3,3)*deprel(3)
                    tau   = sqrt(tau1**2 + tau2**2) ! Contrainte tangente elastique ( abs(tau) )
                    sigma = vh(1,1)*depreln         ! Contrainte normale elastique
                endif

                !--------------------------------------------------------------------------
                !--- Mecanisme de rupture en ouverture pure
                if ((tau/C - sigma/RT) <= 0.D0) then
                    depreln0 = RT/vh(1,1)
                    !===> DEB AJOUT JLT
                    if(sigma<RT) then
                        depreln0 = depreln !sigma/vh(1,1)
                    endif
                    !<=== FIN AJOUT JLT
                    deprelt0 = depreln0/depreln*deprelt

                !--------------------------------------------------------------------------
                !--- Mecanisme de rupture en cisaillement
                elseif ((tau/C - sigma/RT) > 0.D0) then
                    deprelt0 = C/vh(2,2)
                    !===> DEB AJOUT JLT
                    if (tau < C) then
                        deprelt0 = deprelt !tau/vh(2,2)
                    endif
                    !<=== FIN AJOUT JLT
                    depreln0 = deprelt0/deprelt*depreln

                endif

                ietatpg(ie,:)=3 ! Declaration de l'element endommage

                var0 = sqrt(depreln0**2 + deprelt0**2) + 1.D-20
                if (depreln0 <= 0.d0) var0 = sqrt(deprelt0**2) + 1.D-20
                vars  = var0
                varc = xi

                var02 = xp
                vars2 = var02
                varc2 = xc

                deprelni = depreln0
                deprelti = deprelt0

            else

                !--------------------------------------------------------------------------
                !--- Cas des elements deja endommages
                if (ietatpg(ie,ipg)/=0) then
                    vcrit=.true.
                    var02 = xp
                    varc = xi
                    varc2 = xc
                endif
            endif

            !------------------------------------------------------------------------------
            !--- Parametres complementaires du modele
            x0 = var0
            xi = 2.* Wc/(vh(1,1)*var0)
            if ((Wc == 0.d0) .or. (xi<=x0)) xi = x0


            !------------------------------------------------------------------------------
            !--- Si element endommage
            if (vcrit) then

                !--------------------------------------------------------------------------
                !--- variable pilotant l'endommagement D1
                var  = dsqrt((1.-penal)*depreln**2 + deprelt**2)         ! Deplacement relatif total

                !--------------------------------------------------------------------------
                !--- Cas ou l'endommagement D1 evolue
                if (var>=vars) then

                    !----------------------------------------------------------------------
                    !--- Cas de l'evolution de l'endommagement D1
                    !    Actualisation du seuil d'endommagement
                    vars = var

                    !----------------------------------------------------------------------
                    !--- loi d'evolution d'endommagement
                    x = var
                    F_var = 0.d0
                    if ((x>=x0).and.(x<xi)) then
                        F_var  = 1.d0 - (xi-x)/(xi-x0)
                    else
                        F_var  = 1.d0
                    endif

                    !----------------------------------------------------------------------
                    !--- Calcul de l'endommagement
                    D = 1.d0 - var0/var * (1.d0 - F_var)

                    if (D<0.d0) D=0.d0
                    D = max(Dini,D)
                    D = min(D,1.d0)

                    !----------------------------------------------------------------------
                    !--- Calcul du module vhep (imet==1 et imet==2)
                    if (imeth == 1) then
                        !----- Calcul du module secant (endommage)
                        vhep = (1.-D)*vh

                    elseif (imeth == 2) then
                        !----- Calcul du module tangent
                        vari = dsqrt((1.-penali)*deprelni**2 + deprelti**2) + 1.D-20

                        x = vari
                        F_vari = 0.d0
                        dF_vari = 0.d0
                        if ((x>=x0).and.(x<xi)) then
                            F_vari  = 1.d0 - (xi-x)/(xi-x0)
                            dF_vari = 1.d0 / (xi-x0)
                        else
                            F_vari  = 1.d0
                            dF_vari  = 0.d0
                        endif

                        if (dime==2) then
                            coef1 = (deprelni/vari)**2
                            coef2 = (deprelni*deprelti)/(vari**2)
                            if (coef1>.99d0) then
                                coef1 = 1.d0
                                coef2 = 0.d0
                            endif
                            if (coef1<.01d0) then
                                coef1 = 0.d0
                                coef2 = 0.d0
                            endif
                            vhep(1,1) =  vh(1,1) * ((1.-Dini)- var0*((1.-F_vari)/vari - dF_vari)*coef1)
                            vhep(2,2) =  vh(2,2) * ((1.-Dini)- var0*((1.-F_vari)/vari - dF_vari)*(1.-coef1))
                            vhep(1,2) =  vh(1,2) * (-var0*((1.-F_vari)/vari - dF_vari)*coef2)
                            vhep(1,2) =  vh(2,1) * (-var0*((1.-F_vari)/vari - dF_vari)*coef2)

                        else if (dime==3) then
                            stop 'interf_sig_endo (105) : cas 3D non encore implante'
                        endif

                    endif

                !--------------------------------------------------------------------------
                !--- Cas ou l'endommagement n'evolue pas
                else
                    vhep = (1.-D)*vh

                endif

                !--------------------------------------------------------------------------
                !--- variable pilotant l'endommagement D2
                var2 = (1.-penal)*dsqrt(depreln**2 + deprelt**2)         ! Deplacement relatif total

                !--------------------------------------------------------------------------
                !--- Cas ou l'endommagement D2 evolue
                if (var2>=vars2) then

                    !----------------------------------------------------------------------
                    !--- Actualisation du seuil d'endommagement
                    vars2 = var2

                    !----------------------------------------------------------------------
                    !--- loi d'evolution d'endommagement
                    x = var2
                    F_var = 0.d0
                    if ((x>=xp).and.(x<xc)) then
                        F_var  = 1.d0 - (xc-x)/(xc-xp)
                    else
                        F_var  = 1.d0
                    endif

                    !----------------------------------------------------------------------
                    !--- Calcul de l'endommagement
                    D2 = 1.d0 - var02/var2 * (1.d0 - F_var)

                    if (D2<0.d0) D2=0.d0
                    D2 = max(Dini2,D2)
                    D2 = min(D2,1.d0)

                    !----------------------------------------------------------------------
                    !--- Calcul du module vhep (imet==1 et imet==2)
                    if (imeth == 1) then
                        !----- Calcul du module secant (endommage)
                        vhep2 = (1.-D2)*vhendo

                    elseif (imeth == 2) then
                        !----- Calcul du module tangent
                        vari = (1.-penali)*dsqrt(deprelni**2 + deprelti**2) + 1.D-20

                        x = vari
                        F_vari = 0.d0
                        dF_vari = 0.d0
                        if ((x>=xp).and.(x<xc)) then
                            F_vari  = 1.d0 - (xc-x)/(xc-xp)
                            dF_vari = 1.d0/(xc-xp)
                        else
                            F_vari  = 1.d0
                            dF_vari  = 0.d0
                        endif

                        if (dime==2) then
                            coef1 = (deprelni/vari)**2
                            coef2 = (deprelni*deprelti)/(vari**2)
                            if (coef1>.99d0) then
                                coef1 = 1.d0
                                coef2 = 0.d0
                            endif
                            if (coef1<.01d0) then
                                coef1 = 0.d0
                                coef2 = 0.d0
                            endif
                            vhep2(1,1) =  vhendo(1,1) * ((1.-Dini2)- var02*((1.-F_vari)/vari - dF_vari)*coef1)
                            vhep2(2,2) =  vhendo(2,2) * ((1.-Dini2)- var02*((1.-F_vari)/vari - dF_vari)*(1.-coef1))
                            vhep2(1,2) =  vhendo(1,2) * (-var02*((1.-F_vari)/vari - dF_vari)*coef2)
                            vhep2(1,2) =  vhendo(2,1) * (-var02*((1.-F_vari)/vari - dF_vari)*coef2)

                        else if (dime==3) then
                            stop 'interf_sig_endo (108) : cas 3D non encore implante'
                        endif

                    endif

                    if ((var2>varc2) .or. (D2==1.d0)) then
                        if (var2>varc2) D2=1.d0
                        !if (ipg==ipc) then
                            ietatpg(ie,ipg) = 1           ! Element ouvert
                            if (depreln <= 0.d0) then
                                ietatpg(ie,ipg) = 2       ! Element referme
                            endif
                        !endif
                    endif

                !--------------------------------------------------------------------------
                !--- Cas ou l'endommagement D2 n'evolue pas
                else
                    vhep2 = (1.-D2)*vhendo

                endif

                !--------------------------------------------------------------------------
                !--- Cas particulier du traitement de la refermeture de fissure
                !    (on referme sur un comportement elastique partiellement endommage)
                if (depreln <= 0.d0) then
                    !----- Correction sur vhendo
                    vhendo(1,1) = 0.d0
                    vhendo(2,2) = (alphar/alpha)*vhendo(2,2)
                    if (dime==3)  vhendo(3,3) = (alphar/alpha)*vhendo(3,3)
                    !----- Correction sur vhep2
                    vhep2 = (1.-D2)*vhendo
                    !----- Correction sur vhep
                    vhep(1,1) = vh(1,1)
                endif

                vhep = vhep + vhep2

            endif

            !------------------------------------------------------------------------------
            !--- Calcul de vhep
            vhep = vhep + vhpetit

            !------------------------------------------------------------------------------
            !--- Calcul des contraintes (partie matrice)
            vsigm(1) = (1.-D*(1.-penal)) * vh(1,1) * deprel(1)
            vsigm(2) = (1.-D) * vh(2,2) * deprel(2)
            if (dime==3) vsigm(3) = (1.-D) * vh(3,3) * deprel(3)

            !------------------------------------------------------------------------------
            !--- Calcul des contraintes (partie fibres)
            if (vcrit) then
                vsigm(1) = vsigm(1) + (1.-D2) * vhendo(1,1) * deprel(1)
                vsigm(2) = vsigm(2) + (1.-D2) * vhendo(2,2) * deprel(2)
                if (dime==3) vsigm(3) = vsigm(3) + (1.-D2) * vhendo(3,3) * deprel(3)
            endif

            !------------------------------------------------------------------------------
            !--- Sorties : variables internes du modèle
            vin(1) = D
            vin(2) = var0
            vin(3) = vars
            vin(4) = D2
            vin(5) = vars2
            vin(6) = depreln
            vin(7) = deprelt


        case default
            stop 'interf_endo : cas non implante'

    end select

    !--------------------------------------------------------------------------------------
    !--- Sorties :
    vsig = vsigm
    vnle = vepspl
    endo(ie,ipg) = D

end subroutine interf_sig_endo


!------------------------------------------------------------------------------------------------------
!> @author JL Tailhan    (version 1.0 - 2014)
!
!> @brief
!> Pilotage du calcul : gestion du facteur de chargement
!
!> @details
!> #### DESCRIPTION:
!>  A faire
!------------------------------------------------------------------------------------------------------
subroutine interf_pilot(imetpilo, vsol, vduI, vduII, dlam, alphai, MOTi,elemfiss)

    !**************************************************************************************
    !* DEFINITIONS
    !**************************************************************************************

    !--------------------------------------------------------------------------------------
    !--- Modules generaux
    use variables, only :   dime, nelt, idmax, ktypel, nomtype, infele, &
    &                       vprelg, kprop, ietatpg, iedngi, interf, pi, kloce
    use lib_elem, only : elem_B, elem_kloce2, elem_hooke
    use math

    implicit none

    !--------------------------------------------------------------------------------------
    !--- Varibales in/out
    real*8, dimension(:), intent(in) :: vduI, vduII, vsol
    integer, intent(in) :: imetpilo
    real*8, intent(in) :: dlam
    real*8, intent(out) :: alphai
    integer, intent(inout) :: elemfiss
    character(len=8), intent(inout) :: MOTi

    !--------------------------------------------------------------------------------------
    !--- Varibales locales
    real*8, dimension(:), allocatable ::    vdle, vdl0, vsi0, vdsi, &
    &                                       alpg, vdlI, vdlII
    real*8, dimension(:,:), allocatable :: vn, vb, vh, ksig

    real*8 :: deprel(dime), deprel0(dime), vprel(idmax), alph(nelt)
    real*8 ::   alpc, C, RT, phi, psi, sigma, tau, dtau, &
    &           vcrit1, vcrit2, vcrit3, detj, signe, coef, alpgo, alpgc
    real*8 ::   depcrin, depcrit, depreln, deprelt, deprelt1, deprelt2, depreln0, deprelt0
    real*8 ::   vala, valb, valc, delta, alp1, alp2 !, alp0, zero=0.D0
    !real*8 :: ktkn

    character(len=5) :: typel
    character(len=8) :: MOT, MOT1, MOT2
    character(len=8), dimension(nelt) :: rupt
    integer :: ie, ipg, npg, id, iloi, ipc, ndle

    !**************************************************************************************
    !* CORPS DE PROCEDURE
    !**************************************************************************************

    !--------------------------------------------------------------------------------------
    !--- Switch en fonction de la methode de pilotage
    select case(imetpilo)

        !--------------------------------------------------------------------------------------
        !--- Pilotage sur element le plus dangereux par recalcul de vduI
        case(1)

            !----------------------------------------------------------------------------------
            !--- Pour le cas des elements d'interface seulement
            if (interf==1) then

                !------------------------------------------------------------------------------
                !--- Initialisation
                alph = 1.d0
                iedngi = 0

                !------------------------------------------------------------------------------
                !--- Boucle sur les elements
                do ie = 1, nelt

                    MOT = '  '
                    MOT1 = '  '
                    MOT2 = '  '
                    typel = nomtype(ktypel(ie))

                    !----- Pour les elements d'interface seulement
                    if ((typel=='EJQ4'.or.typel=='EJQ6'.or.typel=='EJT6'.or.typel=='EJT8') &
                    &   .and. (all(ietatpg(ie,:)==0) .or. all(ietatpg(ie,:)==3))) then

                        !----- Proprietes elementaires
                        vprel = vprelg(kprop(ie),1:idmax)

                        !----- Recuperation des parametres de la loi
                        iloi = int(vprel(1))

                        !-----  Recuperation des informations sur les elements
                        allocate(ksig(size(infele(ktypel(ie))%Q,1), size(infele(ktypel(ie))%Q,2)))
                        ksig = infele(ktypel(ie))%Q
                        npg = size(ksig,1)
                        allocate(alpg(npg))
                        alpg = 1.d0

                        !----- Vecteur kloce de localisation dans vecteur global
                        call elem_kloce2(ie, ndle)
                        allocate(vdl0(ndle));  allocate(vdle(ndle))
                        allocate(vdlI(ndle));  allocate(vdlII(ndle))

                        !----- Deplacements
                        vdl0  = vsol(kloce(1:ndle))
                        vdlI  = vduI(kloce(1:ndle))
                        vdlII = vduII(kloce(1:ndle))
                        vdle  = vdlI + dlam*vdlII

                        !----- Matrice d'elasticite vh, en fonction du type d'element
                        call elem_hooke(vh,nomtype(ktypel(ie)),vprel)

                        !---- Recuperation du point d'integration central en fonction du type d'element
                        if (typel(1:4)=='EJQ4'.or.typel(1:4)=='EJQ6'.or.typel(1:4)=='EJT6') then
                            ipc = 1
                        else
                            stop 'interf_crit : element non implante'
                        endif

                        !----- Pour toutes les lois on verifie au centre de l'element
                        npg = 1

                        !----- Boucle sur les points de Gauss
                        do ipg = 1, npg

                            !ipg = ipc

                            !----- Calcul des contraintes au centre de gravite de l'element
                            call elem_B(vn,vb,detj,ksig(ipg,:),ie)

                            allocate(vsi0(size(vb,1)))  ;  allocate(vdsi(size(vb,1)))
                            vsi0 = matmul(vh,matmul(vb,vdl0))
                            vdsi = matmul(vh,matmul(vb,vdle))

                            !----- Calcul du deplacement relatif normal au centre de gravite de l'element
                            deprel0 = matmul(vb,vdl0)
                            deprel = matmul(vb,vdle)

                            deallocate(vb,vn)

                            !-------------------------------------------------------------!
                            !----- Cas de la loi 1   !!!
                            !-------------------------------------------------------------!
                            if (iloi==1) then

                            !-------------------------------------------------------------!
                            !----- Cas de la loi 100 !!!
                            !-------------------------------------------------------------!
                            elseif (iloi==100) then
                                id = 6
                                if (dime == 3) id=id-1
                                coef = 1.001
                                RT = coef*vprel(id)
                                C  = coef*vprel(id+1)
                                phi = pi*vprel(id+2)/180.
                                psi = pi*vprel(id+3)/180.
                                if (RT > (C/tan(phi)))  RT = C/tan(phi)

                                sigma = vsi0(1)+vdsi(1)

                                if (dime == 2) then
                                    tau  = vsi0(2)+vdsi(2)
                                    dtau = vdsi(2)
                                else if (dime == 3) then
                                    tau   = sqrt((vsi0(2)+vdsi(2))**2+(vsi0(3)+vdsi(3))**2)
                                    dtau  = sqrt((vdsi(2))**2+(vdsi(3))**2)
                                endif

                                vcrit1 = sigma - RT
                                vcrit2 = abs(tau) + sigma*tan(phi) - C
                                vcrit3 = abs(tau) + sigma*tan(psi-pi/2.)

                                if (vcrit3 < 0.d0) then
                                    if (vcrit1 > 0.d0) then
                                        alpg(ipg) = (RT-vsi0(1))/vdsi(1)
                                        !--- si vsi0 > RT => l'element doit être rompu (ievtp)
                                        if (vsi0(1)-RT >= 0.d0) alpg(ipg) = 0.d0
                                        MOT = 'ouvert'
                                    elseif (vcrit2 >0.d0) then
                                        if (dime==2) then
                                            signe = tau/abs(tau)
                                            alpc = C - signe*vsi0(2) - vsi0(1)*tan(phi)
                                            alpc = alpc/(signe*vdsi(2) + vdsi(1)*tan(phi))
                                            alpg(ipg) = alpc
                                        elseif (dime==3) then
                                            vala =  vdsi(2)*vdsi(2) + vdsi(3)*vdsi(3) - &
                                            &       (vdsi(1)*tan(phi))*(vdsi(1)*tan(phi))
                                            valb =  2.d0*( vsi0(2)*vdsi(2) + vsi0(3)*vdsi(3) + &
                                            &       (C-vsi0(1)*tan(phi))*vdsi(1)*tan(phi) )
                                            valc =  vsi0(2)*vsi0(2) + vsi0(3)*vsi0(3) - &
                                            &       (C-vsi0(1)*tan(phi))*(C-vsi0(1)*tan(phi))
                                            !--- Si valc = 0 alors alpha est pris nul car valc represente le critere a l'iteration
                                            !     precedente. Donc s'il est positif c'est qu'il est deja depasse
                                            if (abs(valc) > epsilon(valc)) then
                                                alpc = 0.d0
                                            else
                                                if (abs(vala) < epsilon(vala)) then
                                                    if ((abs(valb) < epsilon(valb)) .or. (valb < 0.d0)) then
                                                        alpc = 0.d0
                                                    else
                                                        alpc = minval((/ 1.d0 , -valc/valb /))
                                                    endif
                                                else
                                                    delta = (valb * valb) - 4.d0 * vala * valc
                                                    alp2  = -(valb + dsqrt(delta)) / (2.d0 * vala)
                                                    alpc = minval((/ 1.d0 ,alp2 /))
                                                endif
                                            endif
                                            alpg(ipg) = alpc
                                        endif
                                        MOT = 'cisaille'
                                    endif

                                else
                                    if (vcrit1 > 0.d0) then
                                        alpg(ipg) = (RT-vsi0(1))/vdsi(1)
                                        !--- si vsi0 > RT => l'element doit être rompu (ievtp)
                                        if (vsi0(1)-RT >= 0.d0) alpg(ipg) = 0.d0
                                        MOT = 'ouvert'
                                    elseif (vcrit2 > 0.d0) then
                                        if (dime==2) then
                                        signe = tau/abs(tau)
                                        alpc = C - signe*vsi0(2) - vsi0(1)*tan(phi)
                                        alpc = alpc/(signe*vdsi(2) + vdsi(1)*tan(phi))
                                        alpg(ipg) = alpc
                                    elseif (dime==3) then
                                        vala = vdsi(2)*vdsi(2) + vdsi(3)*vdsi(3) - (vdsi(1)*tan(phi))*(vdsi(1)*tan(phi))
                                        valb = 2.d0*( vsi0(2)*vdsi(2) + vsi0(3)*vdsi(3) + (C-vsi0(1)*tan(phi))*vdsi(1)*tan(phi) )
                                        valc = vsi0(2)*vsi0(2) + vsi0(3)*vsi0(3) - (C-vsi0(1)*tan(phi))*(C-vsi0(1)*tan(phi))
                                        !--- Si valc = 0 alors alpha est pris nul car valc represente le critere a l'iteration
                                        !     precedente. Donc s'il est positif c'est qu'il est deja depasse
                                        if (abs(valc) > epsilon(valc)) then
                                            alpc = 0.d0
                                        else
                                            if (abs(vala) < epsilon(vala)) then
                                                if ((abs(valb) < epsilon(valb)) .or. (valb < 0.d0)) then
                                                    alpc = 0.d0
                                                else
                                                    alpc = minval((/ 1.d0 , -valc/valb /))
                                                endif
                                            else
                                                delta = (valb * valb) - 4.d0 * vala * valc
                                                alp2  = -(valb + dsqrt(delta)) / (2.d0 * vala)
                                                alpc = minval((/ 1.d0 ,alp2 /))
                                            endif
                                        endif
                                        alpg(ipg) = alpc
                                    endif
                                    MOT = 'cisaille'
                                endif
                            endif

                            !-------------------------------------------------------------!
                            !----- Cas de la loi 101 !!!
                            !-------------------------------------------------------------!
                            elseif (iloi==101) then

                                id = 6
                                if (dime == 3) id=id-1
                                coef = 1.001
                                RT = coef*vprel(id)
                                C  = coef*vprel(id+1)
                                phi = pi*vprel(id+2)/180.
                                psi = atan(RT/C)

                                sigma = vsi0(1)+vdsi(1)

                                if (dime==2) then
                                    tau   = vsi0(2)+vdsi(2)
                                elseif (dime==3) then
                                    tau   = sqrt((vsi0(2)+vdsi(2))**2+(vsi0(3)+vdsi(3))**2)
                                endif

                                vcrit1 = sigma - RT
                                vcrit2 = abs(tau) - C
                                vcrit3 = abs(tau) + sigma*tan(psi-pi/2.)

                                if (vcrit3 < 0.d0) then

                                    if (vcrit1 > 0.d0) then
                                        alpg(ipg) = (RT-vsi0(1))/vdsi(1)
                                        !--- si vsi0 > RT => l'element doit être rompu (ievtp)
                                        if (vsi0(1)-RT >= 0.d0) alpg(ipg) = 0.d0
                                        MOT = 'ouvert'

                                    elseif (vcrit2 > 0.d0) then
                                        if (dime==2) then
                                            signe = tau/abs(tau)
                                            alpc = C - signe*vsi0(2)
                                            alpc = alpc/(signe*vdsi(2))
                                            alpg(ipg) = alpc
                                        elseif (dime==3) then
                                            vala = vdsi(2)**2 + vdsi(3)**2
                                            valb = 2.d0*(vsi0(2)*vdsi(2) + vsi0(3)*vdsi(3))
                                            valc = vsi0(2)**2 + vsi0(3)**2 - C**2
                                            ! On doit avoir Tau_0 - C >0 (sinon alpc=0)
                                            if (valc < 0.d0) then
                                                delta = valb**2 - 4.*vala*valc
                                                if (delta >= 0) then
                                                    alp1 = (-valb + sqrt(delta))/(2.*vala)
                                                    alp2 = (-valb - sqrt(delta))/(2.*vala)
                                                    alpc = alp1 ! maxval((/alp1 , alp2/))
                                                else
                                                    stop 'interf_pilot (101 3D): cas impossible (1)'
                                                endif
                                            else
                                                alpc = 0.d0
                                            endif
                                            alpg(ipg) = alpc
                                        endif
                                        MOT = 'cisaille'
                                    endif

                                else

                                    if (vcrit1 > 0.d0) then
                                        alpg(ipg) = (RT-vsi0(1))/vdsi(1)
                                        !--- si vsi0 > RT => l'element doit être rompu (ievtp)
                                        if (vsi0(1)-RT >= 0.d0) alpg(ipg) = 0.d0
                                        MOT = 'ouvert'
                                    elseif (vcrit2 > 0.d0) then
                                        if (dime==2) then
                                            signe = tau/abs(tau)
                                            alpc = C - signe*vsi0(2)
                                            alpc = alpc/(signe*vdsi(2))
                                            alpg(ipg) = alpc
                                        elseif (dime==3) then
                                            vala = vdsi(2)**2 + vdsi(3)**2
                                            valb = 2.d0*(vsi0(2)*vdsi(2) + vsi0(3)*vdsi(3))
                                            valc = vsi0(2)**2 + vsi0(3)**2 - C**2
                                            ! On doit avoir Tau_0 - C >0 (sinon alpc=0)
                                            if (valc < 0.d0) then
                                                delta = valb**2 - 4.*vala*valc
                                                if (delta >= 0) then
                                                    alp1 = (-valb + sqrt(delta))/(2.*vala)
                                                    alp2 = (-valb - sqrt(delta))/(2.*vala)
                                                    alpc = alp1 ! maxval((/alp1 , alp2/))
                                                else
                                                    stop 'interf_pilot (101 3D): cas impossible (2)'
                                                endif
                                            else
                                                alpc = 0.d0
                                            endif
                                            alpg(ipg) = alpc
                                        endif
                                        MOT = 'cisaille'
                                    endif
                                endif

                            !-------------------------------------------------------------!
                            !----- Cas de la loi 102 !!!
                            !-------------------------------------------------------------!
                            elseif (iloi==102) then

                                alpgc = 1.d0
                                alpgo = 1.d0
                                id = 6
                                if (dime == 3) id=id-1
                                coef = 1.00001
                                depcrin = coef*vprel(id+4)
                                depcrit = coef*vprel(id+5)
                                C  = coef*vprel(id+1)
                                RT = coef*vprel(id)

                                if (depreln < 0.d0) depreln = 0.d0
                                depreln0 = RT/vh(1,1)
                                deprelt0 = C/vh(2,2)
                                if (depcrin < depreln0) depcrin = depreln0
                                if (depcrit < deprelt0) depcrit = deprelt0

                                ! Calcul des deplacements relatif
                                depreln  = deprel0(1) + deprel(1)
                                deprelt1 = deprel0(2) + deprel(2)
                                if (dime == 2) deprelt2 = 0.d0
                                if (dime == 3) deprelt2 = deprel0(3) + deprel(3)

                                ! Deplacement tangent relatif total
                                deprelt = sqrt(deprelt1**2. + deprelt2**2.)

                                if (abs(deprelt) >= depcrit) then
                                    MOT1 = 'cisaille'
                                    if (dime==2) then
                                        alpgc = (depcrit - abs(deprel0(2)))/abs(deprel(2))

                                    elseif (dime==3) then

                                        ! (ddt1^2+ddt2^2)alpha^2 + 2(ddt1.dt01+ddt2.dt02)alpha + dt01^2+dt02^2-dc^2 = 0
                                        !        a    .  alpha^2 +           b     .     alpha +         c          = 0

                                        vala = deprel(2)**2. + deprel(3)**2.
                                        valb = 2.*(deprel(2)*deprel0(2) + deprel(3)*deprel0(3))
                                        valc = deprel0(2)**2. + deprel0(3)**2. - depcrit**2.

                                        delta = valb**2. - 4.*vala*valc

                                        if (delta >= 0) then
                                            alp1 = (-valb + sqrt(delta))/(2.*vala)
                                            alp2 = (-valb - sqrt(delta))/(2.*vala)

                                            ! si alp = 0 !!!!!

                                            alpgc = maxval((/alp1 , alp2/))
                                        else
                                            !print*, 'interf_pilot - attention : discriminent négatif '
                                            alpgc = (depcrit - abs(deprelt))/depcrit
                                        endif
                                    endif
                                endif

                                if (depreln >= depcrin) then
                                    MOT2 = 'ouvert'
                                    alpgo = (depcrin - deprel0(1))/deprel(1)
                                endif

                                alpg(ipg) = minval((/alpgo, alpgc/))
                                if (alpg(ipg)==alpgo) MOT = MOT2
                                if (alpg(ipg)==alpgc) MOT = MOT1

                            !-------------------------------------------------------------!
                            !----- Cas des lois 105, 106, 107, 108 et 158 !!!
                            !-------------------------------------------------------------!
                            elseif ((iloi==105).or.(iloi==106).or.(iloi==107) &
                            &    .or. (iloi==108).or.(iloi==158)) then

                                if (all(ietatpg(ie,:)==0)) then
                                    id = 6
                                    if (dime == 3) id=id-1
                                    coef = 1.0001
                                    RT = coef*vprel(id)
                                    C  = coef*vprel(id+1)
!                                       !----- Petit traitement particulier pour loi 108
!                                        !   (prise en compte du rapport kt/kn
!                                        if (iloi==108) then
!                                           ktkn = vprel(id+8)
!                                           vsi0(2) = ktkn*vsi0(2)
!                                           vdsi(2) = ktkn*vdsi(2)
!                                           if (dime==3) then
!                                               vsi0(2) = ktkn*vsi0(2)
!                                               vdsi(2) = ktkn*vdsi(2)
!                                           endif
!                                         endif

                                    !----- Contrainte normale
                                    sigma = vsi0(1)+vdsi(1)

                                    !----- Contrainte tangentielle ( Attention ici : tau = abs(tau) )
                                    if (dime==2) then
                                        tau = sqrt((vsi0(2)+vdsi(2))**2)
                                    else if (dime==3) then
                                        tau = sqrt((vsi0(2)+vdsi(2))**2+(vsi0(3)+vdsi(3))**2)
                                        !stop 'interf_pilot (lois 105-107) : Cas 3D non encore implanté'
                                    endif

                                    !----- Criteres
                                    vcrit1 = sigma - RT
                                    vcrit2 = tau - C  ! Attention ici : tau = abs(tau) !
                                    vcrit3 = abs(tau)/C - sigma/RT

                                    if (vcrit3 < 0.d0) then
                                        !----- Critere en traction pure
                                        if (vcrit1 > 0.d0) then
                                            alpg(ipg) = (RT-vsi0(1))/vdsi(1)
                                            MOT = 'ouvert'
                                        endif
                                    else
                                        !----- Critere en cisaillement
                                        if (vcrit2 > 0.d0) then
                                            vala = vdsi(2)**2
                                            valb = 2.d0*(vsi0(2)*vdsi(2))
                                            valc = vsi0(2)**2 - C
                                            if (dime == 3) then
                                                vala = vala + vdsi(3)**2
                                                valb = valb + 2.d0*(vsi0(3)*vdsi(3))
                                                valc = valc + vsi0(3)**2
                                            endif
                                            delta = valb**2. - 4.*vala*valc
!                                               if (valc > 0.d0) then
!                                                 alpc = 0.d0     ! il n'est pas necessaire d'augmenter le deplacement !
!                                               else
                                            ! le discriminent est alors toujours positif
                                            if (delta < 0.d0) then
                                                print*, 'interf_pilot (105): cas theoriquement impossible (delta=',delta,'<0)'
                                                stop
                                            endif
                                            alp1 = (-valb - sqrt(delta))/(2.*vala)
                                            alp2 = (-valb + sqrt(delta))/(2.*vala)
                                            alpc = alp2 ! maxval((/alp1 , alp2/))
                                            if ((alp1 < 0.d0) .and. (alp2 < 0.d0)) alpc = 0.d0
                                            if ((alp1 > 0.d0) .and. (alp2 > 0.d0)) alpc = minval((/alp1 , alp2/))
                                            if ((alp1 < 0.d0) .and. (alp2 > 0.d0)) alpc = alp2
                                            if (alpc>1.d0) alpc=1.d0
                                            if (alpc<0.d0) alpc=0.d0
!                                              endif
                                            alpg(ipg) = alpc
                                            MOT = 'cisaille'

                                        endif
                                    endif
                                endif
                            else
                                stop "FIDES_interf_pilot : loi non encore programmee"
                            endif

                            deallocate(vsi0,vdsi)

                        enddo
                        !----- Fin de boucle sur les points de Gauss

                        deallocate(vdl0,vdlI,vdlII,vdle)
                        deallocate(vh)

                        rupt(ie) = MOT
                        alph(ie) = minval(alpg)

!if (abs(Rt)<epsilon(Rt)) print*,'pilot fin->',ie,Rt,vcrit1, rupt(ie), alph(ie)

                        deallocate(alpg,ksig)
                    endif
                enddo
                !------------------------------------------------------------------------------
                !--- Fin de la boucle d'element

                alphai = minval(alph)

                if (alphai < 1.d0) then
                    iedngi = find_num(alph,alphai)
                    print*,'dans interf_pilot',iedngi,alphai
                    MOTi = rupt(iedngi)
                    elemfiss = count(alph<=0)
                else
                    alphai = 1.d0
                endif
            endif

case default
    stop "FIDES_interf_pilot : loi non encore programmee"

end select

end subroutine interf_pilot


!------------------------------------------------------------------------------------------------------
!> @author JL Tailhan    (version 1.0 - 2014)
!
!> @brief
!> Gestion de l'interpenetration et des oscillations dans l'etat d'un element fissure
!> Traite aussi le cas d'un element totalement detache
!
!> @details
!> #### DESCRIPTION:
!>  A faire
!------------------------------------------------------------------------------------------------------
subroutine interf_change_etat()

    !**************************************************************************************
    !* DEFINITIONS
    !**************************************************************************************

    !--------------------------------------------------------------------------------------
    !--- Modules generaux
    use variables, only :   nelt, ktypel, nomtype, infele, ietatpg, osc_interf,&
    &                       histetatpg1, interf, dime, kloce,& !, comptoscill, ouvermax
    !    les deplacements et increments de deplacements
    &                       vsol,& ! vdu, &
    !    et les variables necessaires a la detection des elements libres
    &                       listelemautour, kconec, elemlibre
    use lib_elem, only : elem_B, elem_kloce2

    implicit none

    !--------------------------------------------------------------------------------------
    !--- Variables locales
    !    Quantites globales
    real*8, dimension(:),   allocatable :: vdle
    real*8, dimension(:,:), allocatable :: vn, vb, ksig

    !    Quantites principales : Deplacement, ouverture, etat
    real*8 :: detj, depreln, deprel(dime)
    integer :: ie, ipg, npg, ietat0, ietat1, ndle
    character(len=5) :: typel
    integer:: limite !< nombre de chagement d'etats a partir duquel on considere qu'il y a oscillation

    !    Detection des elements massifs libres
    integer :: nodetachel(nelt)
    integer :: ino,nbino,iel,ele,nnele,numno !,nnoface
    character(len=5) :: typele

    !**************************************************************************************
    !* CORPS DE PROCEDURE
    !**************************************************************************************

    limite=5

    if (interf==1) then

        nodetachel = 0

        !--------------------------------------------------------------------------------------
        !--- Boucle sur les elements
        do ie = 1, nelt

            !----------------------------------------------------------------------------------
            !--- Recuperation du type de l'element
            typel = nomtype(ktypel(ie))

            !----------------------------------------------------------------------------------
            !--- Pour les elements d'interface seulement
            if ((typel=='EJQ4') .or. (typel=='EJQ6') .or. (typel=='EJT6')) then

                !------------------------------------------------------------------------------
                !--- Detection des elements massifs libres
                ! (a partir des elements d'interface ouverts qui le bordent)
                if (all(ietatpg(ie,:)/=0).and.all(ietatpg(ie,:)/=3)) then

                    ! la boucle n'est faite que sur les noeuds extremes
                    nbino = infele(ktypel(ie))%nnel
                    if (typel=='EJQ6') nbino = nbino-2

                    do ino=1, nbino
                        numno = kconec(ie,ino)

                        do iel = 1,size(listelemautour(numno)%el,1)

                                ele = listelemautour(numno)%el(iel)

                                typele = nomtype(ktypel(ele))
                                nnele = infele(ktypel(ele))%nnel

                                if (typele(1:3)=='MBT') then
                                    nodetachel(ele) = nodetachel(ele) + 1

                                    ! un noeud connecte à deux elements (1EJ et 1 MB) est forcement
                                    ! sur un bord libre
                                    if (size(listelemautour(numno)%el)==2) nodetachel(ele) = nodetachel(ele) + 1

                                endif

                        enddo
                    enddo
                endif
                !------------------------------------------------------------------------------

                !------------------------------------------------------------------------------
                !--- Traitement du changement d'etat des elements
                !   - Comptabilisation des oscillations (changements d'etat recurrents)
                !   - Comptabilisation d'un changement d'etat lors du calcul des contraintes
                !   - Détection de l'interpenetration d'un element
                !----------------------------------------------------------------------------

                !------------------------------------------------------------------------------
                !--- Recuperation des informations sur les elements
                allocate(ksig(size(infele(ktypel(ie))%Q,1), size(infele(ktypel(ie))%Q,2)))
                ksig = infele(ktypel(ie))%Q
                npg = size(ksig,1)

                !------------------------------------------------------------------------------
                !--- Boucle sur les points d'integration
                do ipg = 1, npg

                    !--------------------------------------------------------------------------
                    !--- On considere les elements fissures
                    if (ietatpg(ie,ipg) /= 0) then

                        !================================================================================!
                        !             NOMBRE MAXIMAL D'OSCILLATIONS ATTEINT POUR UN ELEMENT              !
                        !                (on bloque l'etat de l'element en ouverture)                    !
                        !----- on bloque l'element dans une position refermee jusqu'a convergence totale !
                        !JLT - 2016-04_05 if ((osc_interf(ie)%el(ipg) == limite)) then
                        if ((osc_interf(ie)%el(ipg) == limite).and.(ietatpg(ie,ipg)==2)) then
                            ! On bloque l'etat de l'element a : -2 (ouvert)
                            ietatpg(ie,ipg)=-2
                            ! la ligne suivante est pour eviter de repasser dans le test la prochaine fois
                            !JLT - 2016-04_05 osc_interf(ie)%el(ipg)=osc_interf(ie)%el(ipg)+1
                            osc_interf(ie)%el(ipg)=0
                        endif
                        !================================================================================!

                        !----------------------------------------------------------------------
                        !--- Vecteur kloce de localisation pour assemblage
                        call elem_kloce2(ie, ndle)
                        allocate(vdle(ndle))

                        !----------------------------------------------------------------------
                        !--- Deplacement total
                        vdle  = vsol(kloce(1:ndle))

                        !----------------------------------------------------------------------
                        !--- Calcul du deplacement relatif normal au centre de gravite de l'element
                        call elem_B(vn,vb,detj,ksig(ipg,:),ie)
                        deprel = matmul(vb,vdle)
                        depreln = deprel(1)

                        !================================================================================!
                        !             POUR LES ELEMENTS MAINTENUS A (ietatpg=-2)
                        !----- liberation de l'etat bloque à -2
                        if (ietatpg(ie,ipg)==-2) then
                            osc_interf(ie)%el(ipg)=osc_interf(ie)%el(ipg)+1
                            !JLT - 2016-04_05 if (osc_interf(ie)%el(ipg)==limite+5) then
                            if (osc_interf(ie)%el(ipg)==limite) then
                                ! On libere element
                                osc_interf(ie)%el(ipg)=0
                                ietatpg(ie,ipg) = 1
                            endif
                            histetatpg1(ie,ipg) = ietatpg(ie,ipg)!<------- Etat actuel
                        endif
                        !================================================================================!

                        !----------------------------------------------------------------------
                        !--- Pour les elements en ouverture ou en refermeture
                        if ((ietatpg(ie,ipg)==1) .or. (ietatpg(ie,ipg)==2)) then

                            ietat1 = histetatpg1(ie,ipg) !<-- Etat precedent
                            ietat0 = ietatpg(ie,ipg)     !<-- Etat actuel

                            !------------------------------------------------------------------
                            !--- On teste d'abord un changement d'etat survenu lors du calcul des contraintes
                            !      (avant actualisation du deplacement) et on comptatbilise une oscillation
                            if ((ietat0/=ietat1) .and. (ietat1 /= 0)) then
                                osc_interf(ie)%el(ipg)=osc_interf(ie)%el(ipg)+1
                            endif

                            !------------------------------------------------------------------
                            !--- On teste ensuite une refermeture liée a la correction du deplacement a
                            !      l'iteration: on change l'etat de l'element et on comptatbilise une oscillation
                            if ((ietatpg(ie,ipg)==1) .and. (depreln < 0.d0)) then
                                ietatpg(ie,ipg) = 2
                                !----- on comptabilise une oscillation
                                osc_interf(ie)%el(ipg)=osc_interf(ie)%el(ipg)+1
                            endif

                            histetatpg1(ie,ipg) = ietat0

                        endif

                        deallocate(vdle,vb,vn)

                    endif
                enddo
                !----- Fin de boucle sur les points d'integration

                deallocate(ksig)

            endif
        enddo

        !--------------------------------------------------------------------------------------
        !--- Finalisation du controle d'elements libres
        do iel=1,nelt
            ! On detecte un element libre
            if(nodetachel(iel)==size(infele(ktypel(iel))%face)) then
                elemlibre(iel)=1
                !print*,'detection de l''element detache', iel
                call elem_kloce2(iel, ndle)
                vsol(kloce(1:ndle))=0.d0
                !pause
            endif
        enddo

    endif

end subroutine interf_change_etat


!------------------------------------------------------------------------------------------------------
!> @author JL Tailhan    (version 1.0 - 2014)
!
!> @brief
!> Fonction de detection de la rupture des elements d'interface
!
!> @details
!> #### DESCRIPTION:
!>  A faire
!------------------------------------------------------------------------------------------------------
subroutine interf_rupture(ie,vsig,vnle)

    !**************************************************************************************
    !* DEFINITIONS
    !**************************************************************************************

    !--------------------------------------------------------------------------------------
    !--- Modules generaux
    use variables, only : ktypel, nomtype, ietatpg

    implicit none

    !--------------------------------------------------------------------------------------
    !--- Variables in/out
    real*8, dimension(:,:), intent(inout) :: vsig, vnle
    integer, intent(in) :: ie

    !--------------------------------------------------------------------------------------
    !--- Variables locales
    character(len=5) :: typel

    !**************************************************************************************
    !* CORPS DE PROCEDURE
    !**************************************************************************************

    !--------------------------------------------------------------------------------------
    !--- Mise a zero des contraintes et deformations anelastiques
    typel = nomtype(ktypel(ie))

    if ((typel=='EJQ4').or.(typel=='EJQ6').or.(typel=='EJT6')) then

        if (all(ietatpg(ie,:) == 1) .or. all(ietatpg(ie,:) == -2)) then
            vnle = 0.d0
            vsig = 0.d0

        elseif (all(ietatpg(ie,:) == 2)) then
            vnle = 0.d0
        endif
    endif

end subroutine interf_rupture


!------------------------------------------------------------------------------------------------------
!> @author JL Tailhan    (version 1.0 - 2014)
!
!> @brief
!> Gestion de la raideur
!
!> @details
!> #### DESCRIPTION:
!>  A faire
!------------------------------------------------------------------------------------------------------
subroutine interf_modul(ie,ipg,vh,vprel)

    !**************************************************************************************
    !* DEFINITIONS
    !**************************************************************************************

    !--------------------------------------------------------------------------------------
    !--- Modules generaux
    use variables, only : dime, ktypel, nomtype, ietatpg

    implicit none

    !--------------------------------------------------------------------------------------
    !--- Variables in/out
    real*8, dimension(:), intent(in) :: vprel
    real*8, dimension(:,:), intent(inout) :: vh
    integer, intent(in) :: ie, ipg

    !--------------------------------------------------------------------------------------
    !--- Variables locales
    integer :: iloi
    character(len=5) :: typel

    !**************************************************************************************
    !* CORPS DE PROCEDURE
    !**************************************************************************************

    iloi=int(vprel(1))
    typel = nomtype(ktypel(ie))

    if (((typel=='EJQ4').or.(typel=='EJQ6').or.(typel=='EJT6')) &
    &   .and.  ((iloi==100) .or. (iloi==101) .or. (iloi==102) &
    &   .or.    (iloi==105) .or. (iloi==106) .or. (iloi==107) &
    &   .or.    (iloi==108) .or. (iloi==158) )) then

        !--------------------------------------------------------------------------------------
        !--- Modification de la matrice pour element casse
        if ((ietatpg(ie,ipg)==1).or.(ietatpg(ie,ipg)== -1)) then
            vh = vh/1.0d50
            !vh(1,1) = 1.d-20
            !vh(2,2) = 1.d-20

        !--------------------------------------------------------------------------------------
        !--- Modification de la matrice pour element cisaille
        elseif ((ietatpg(ie,ipg)==2).or.(ietatpg(ie,ipg)==-2)) then
            if (dime == 2) then
                vh(2,2) = vh(2,2)/1.0d50
                !vh(2,2) = 1.d-20
            elseif (dime == 3) then
                vh(2,2) = vh(2,2)/1.0d50
                vh(3,3) = vh(3,3)/1.0d50
            endif

        !--------------------------------------------------------------------------------------
        !--- Elastoplastique
        !elseif (ietatpg(ie,ipg)==3) then

        endif
    endif

end subroutine interf_modul


!------------------------------------------------------------------------------------------------------
!> @author JL Tailhan    (version 1.0 - 2014)
!
!> @brief
!> Verification de la compatibilite loi-element
!
!> @details
!> #### DESCRIPTION:
!>  A faire
!------------------------------------------------------------------------------------------------------
subroutine interf_loi(iloi,icomp,ie)!,ipg)

    !**************************************************************************************
    !* DEFINITIONS
    !**************************************************************************************

    !--------------------------------------------------------------------------------------
    !--- Modules generaux
    use variables, only : ktypel, nomtype

    implicit none

    !--------------------------------------------------------------------------------------
    !--- Variables in/out
    integer, intent(in) :: icomp, ie!, ipg
    integer, intent(out) :: iloi

    !--------------------------------------------------------------------------------------
    !--- Variables locales
    character(len=5) :: typel

    !**************************************************************************************
    !* CORPS DE PROCEDURE
    !**************************************************************************************

    typel = nomtype(ktypel(ie))

    if ((typel=='EJQ4').or.(typel=='EJQ6').or.(typel=='EJT6')) then
        if      ((icomp/=100).and.(icomp/=101).and.(icomp/=102) &
        &   .and.(icomp/=105).and.(icomp/=106) &
        &   .and.(icomp/=107).and.(icomp/=108).and.(icomp/=158)) then
            stop 'Discordance de la loi de l"element d"interface : Veuillez verifier !!! '
        !else
            !if (iloi==102) then
                !if (ietatpg(ie,ipg)==1 .or. ietatpg(ie,ipg)==2) iloi = 101
            !endif
        endif
    else
        iloi = icomp
    endif

end subroutine interf_loi


!------------------------------------------------------------------------------------------------------
!> @author JL Tailhan    (version 1.0 - 2014)
!
!> @brief
!> Initialisations liees a la gestion des etats de l'elements d'interface
!
!> @details
!> #### DESCRIPTION:
!>  A faire
!------------------------------------------------------------------------------------------------------
subroutine interf_init()

    !**************************************************************************************
    !* DEFINITIONS
    !**************************************************************************************

    !--------------------------------------------------------------------------------------
    !--- Modules generaux
    use variables, only :   nelt, nomtype, infele, interf, ietatpg, &
    &                       histetatpg1, histetatpg2, iedngi, endo
    use initialisation, only : init_mat

    implicit none

    !--------------------------------------------------------------------------------------
    !--- Variables locales
    integer, dimension(size(nomtype,1)) :: npg
    integer :: itype, ntype, npgm
    character(len=5) :: typel

    !**************************************************************************************
    !* CORPS DE PROCEDURE
    !**************************************************************************************

    interf = 0
    ntype = size(nomtype,1)
    npg = 0

    do itype = 1, ntype
        typel = nomtype(itype)
        if ((typel=='EJQ4').or.(typel=='EJQ6').or.(typel=='EJT6')) then

            npg(itype) = size(infele(itype)%W,1)
            interf = 1

        endif
    enddo

    if (interf==1) then
        iedngi = 0
        npgm = maxval(npg)
        call init_mat(ietatpg,nelt,npgm)
        call init_mat(endo,nelt,npgm)
        call init_mat(histetatpg1,nelt,npgm)
        call init_mat(histetatpg2,nelt,npgm)

    endif

end subroutine interf_init


!------------------------------------------------------------------------------------------------------
!> @author JL Tailhan    (version 1.0 - 2014)
!
!> @brief
!> Stockage des resultats à chaque pas de temps pour le modele de fissuration
!
!> @details
!> #### DESCRIPTION:
!>  Avec  gestion  des transitions etat(-1)--> etat(1), etat(-2)--> etat(2)
!------------------------------------------------------------------------------------------------------
subroutine interf_stock()

    !**************************************************************************************
    !* DEFINITIONS
    !**************************************************************************************

    !--------------------------------------------------------------------------------------
    !--- Modules generaux
    use variables, only :  interf, ietatpg, histetatpg1

    implicit none

    !**************************************************************************************
    !* CORPS DE PROCEDURE
    !**************************************************************************************

    if (interf==1) then
        where (ietatpg==-1) ietatpg = 1
        where (ietatpg==-2) ietatpg = 1
        histetatpg1 = ietatpg
    endif

end subroutine interf_stock


!------------------------------------------------------------------------------------------------------
!> @author JL Tailhan    (version 1.0 - 2014)
!
!> @brief
!> Fonction de definition des proprietes mecaniques aleatoires
!
!> @details
!> #### DESCRIPTION:
!>   - kprop : numerotation des groupes des elements\n
!>   - vprelg : proprietes mecaniques (fonction du modele)\n\n
!
!>   Valeurs en sortie :\n
!>   - kprop : actualise en fonction du numero de l'element\n
!>             (1 element probabiliste = 1 groupe)\n
!>   - vprelg : modifie pour prendre en compte la propriete\n
!>      aleatoire consideree - Attention depend du modele choisi\n
!------------------------------------------------------------------------------------------------------
subroutine interf_distal()

    !**************************************************************************************
    !* DEFINITIONS
    !**************************************************************************************

    !--------------------------------------------------------------------------------------
    !--- Modules generaux
    use variables, only :   nelt, dime, kconec, ktypel, nomtype, &
    &                       infele, vprelg, kprop, kprop0, interf, fiss, idmax, pi, &
    &                       Dg, fc, &
    &                       vitrav, &
    !    Gestion de l'alea
    &                       alea, young, resist, nrjbeto, nrjfibr
    use lib_elem, only :    elem_B
    use initialisation, only :  init_vec, init_mat
    use proprietes_aleatoires
    use math, only : find_vec

    implicit none

    !--------------------------------------------------------------------------------------
    !--- Variables locales
    real*8, dimension(:,:), allocatable ::  vprelg0, vn, vb, ksig
    real*8, dimension(:), allocatable ::    vpg, RR, EE, CC, WFF, WBB
    real*8 ::   Rt(1), E(1), C(1), Vel(nelt), WF(1), WB(1)
    real*8 ::   loi1(5), loi2(5), loi3(5), loi4(5)
    real*8 ::   detj, detj1, detj2, epais, poids, Ve, Se, Vg, Sg,&
    &           r, Vem, alphab, betab, gammab, &
    &           alphac, betac, gammac

    integer, dimension(:), allocatable :: ielm, ielf
    integer, dimension(:,:), allocatable :: cote
    integer ::  iebe, iefi, ienb, ienf
    integer ::  i, ie, ie0, nbg, ncote, np, ipg, npg, id, id1, nelf, &
    &           iloi, icote, dPlay !, cal
    integer ::  igfp, igma, igbt, igfb, j, k , l, m, n

    character(len=5) :: typea, typeb
    logical ::  iprobVO, iprobFI, iprobBT, iprobFB, iprobMY, &
    &           resp, modp, betp, fibp

    !**************************************************************************************
    !* CORPS DE PROCEDURE
    !**************************************************************************************

    !--------------------------------------------------------------------------------------
    !--- Dans le cas d'un calcul avec elements d'interface uniquement
    if ((fiss==0) .and. (interf == 1) .and. (alea)) then

        !----------------------------------------------------------------------------------
        !--- preparation de la renumerotation de kprop et vprelg
        ! Principe : On conserve la numerotation pour tous les
        ! groupes dont l'ancien numero est inferieur a nga. Puis on
        ! "decale" d'un cran dans la numerotation tous les groupes
        ! qui ont un ancien numero plus grand que nga. Enfin, tous
        ! les elements du groupe nga sont renumerote en nga+i.
        ! nga=nombre de groupes anciens ou initial (variable NG)

        !----------------------------------------------------------------------------------
        !--- Elements massifs probabilistes pour beton sain
        iebe  = 0

        !----------------------------------------------------------------------------------
        !--- Elements d'interface a resistance probabiliste
        iefi  = 0

        !----------------------------------------------------------------------------------
        !--- Elements d'interface a energie de fissuration probabiliste
        ienb  = 0

        !----------------------------------------------------------------------------------
        !--- Elements d'interface a energie de pontage fibre/matrice
        ienf  = 0

        !----------------------------------------------------------------------------------
        !--- Elements probabilistes totaux
        ie0 = 0

        nbg = maxval(kprop)
        call init_mat(vprelg0,(size(kprop,1)+nbg),size(vprelg,2))
        vprelg0(1:nbg,:) = vprelg(1:nbg,:)

        Vel = 0.d0
        call init_vec(RR,nelt)
        call init_vec(EE,nelt)
        call init_vec(CC,nelt)
        call init_vec(WBB,nelt)
        call init_vec(WFF,nelt)

        !----------------------------------------------------------------------------------
        !--- Indice pour trace de distributions de proprietes aleatoires
        iprobVO  = .false.
        iprobFI  = .false.
        iprobBT  = .false.
        iprobFB  = .false.
        iprobMY  = .false.

        !----------------------------------------------------------------------------------
        !--- Traitement des elements d'interface probabilistes d'abord
        do ie = 1, nelt

            !------------------------------------------------------------------------------
            !--- detection d'un element probabiliste (resistance et/ou energie)
            resp = .false.
            if (allocated(resist)) then
                if (count(resist%num==kprop0(ie))==1) resp=.true.
            endif

            !------------------------------------------------------------------------------
            !--- detection d'un element d'interface a energie de fissuration de beton
            betp = .false.
            if (allocated(nrjbeto)) then
                if (count(nrjbeto%num==kprop0(ie))==1) betp=.true.
            endif

            !------------------------------------------------------------------------------
            !--- detection d'un element d'interface a energie de pontage matrice/fibre
            fibp = .false.
            if (allocated(nrjfibr)) then
                if (count(nrjfibr%num==kprop0(ie))==1) fibp=.true.
            endif

            typea = nomtype(ktypel(ie))

            !------------------------------------------------------------------------------
            !--- si l'element appartient a un groupe probabiliste (resistance)
            if ((resp) .and. (typea=='EJQ4' .or. typea=='EJQ6' .or. typea=='EJT6')) then

                allocate(vitrav(1))
                vitrav= find_vec(resist%num,kprop0(ie))
                igfp = vitrav(1)
                deallocate(vitrav)

                !--------------------------------------------------------------------------
                !--- deb fissure probabiliste
                ie0 = ie0 + 1
                iefi = iefi + 1

                !--------------------------------------------------------------------------
                !--- Recherche des elements massifs voisins
                if (typea == 'EJQ4') then
                    allocate(cote(2,2))
                elseif (typea == 'EJQ6' .or. typea == 'EJT6') then
                    allocate(cote(2,3))
                else
                    stop 'interf_distal : element non encore implante'
                endif

                cote = infele(ktypel(ie))%face
                ncote = size(cote,1)
                np = size(cote,2)

                if (ncote/=2) then
                    stop 'interf_distal : erreur de programmation. Il faut 2 cotes'
                endif

                allocate(ielm(ncote))
                ielm = 0

                do icote = 1, ncote
                    do i = 1, nelt
                        if (i/=ie) then
                            k = 0 ; l = 0 ; m = 0 ; n = 0
                            do j = 1, size(kconec,2)
                                if (kconec(i,j)==kconec(ie,cote(icote,1))) k = 1
                                if (kconec(i,j)==kconec(ie,cote(icote,2))) l = 1

                                if (np==3) then
                                    if (kconec(i,j)==kconec(ie,cote(icote,3))) m = 1
                                endif

                                if (np==4) then
                                    if (kconec(i,j)==kconec(ie,cote(icote,3))) m = 1
                                    if (kconec(i,j)==kconec(ie,cote(icote,4))) n = 1
                                endif
                            enddo

                            if ((k==1.and.l==1.and.np==2).or.(k==1.and.l==1.and.m==1.and.np==3) &
                                & .or. (k==1.and.l==1.and.m==1.and.n==1.and.np==4)) then
                                ielm(icote) = i
                                exit
                            endif
                        endif
                    enddo

                    if (icote > 1) then
                        if (ktypel(ielm(1)) /= ktypel(ielm(icote))) &
                        &   stop 'interf_distal : elements massifs incompatibles'
                    endif
                enddo

                !--------------------------------------------------------------------------
                !--- Recuperation des poids et points d'integration de l'element
                allocate(vpg(size(infele(ktypel(ielm(1)))%W)))
                allocate(ksig(size(infele(ktypel(ielm(1)))%Q,1),size(infele(ktypel(ielm(1)))%Q,2)))

                vpg = infele(ktypel(ielm(1)))%W
                ksig = infele(ktypel(ielm(1)))%Q
                npg = size(ksig,1)

                !--------------------------------------------------------------------------
                !--- Recuperation de l'epaisseur si 2D
                epais = 1.
                iloi = int(vprelg(kprop(ielm(1)),1))

                if (iloi == 1) then
                    if ((dime==2).and.(int(vprelg(kprop0(ielm(1)),2))==3)) &
                    &   epais = vprelg(kprop0(ielm(1)),idmax-2) !En contraintes planes (CP) (position de l'epaisseur)
                else
                    print*,'Non encore programme pour la loi ', iloi
                    stop 'interf_distal : valable pour l''elasticite seulement'
                endif

                !--------------------------------------------------------------------------
                !--- Calcul du volume fissure (Ve(ie1)+Ve(ie2))
                Ve = 0.d0
                do ipg = 1, npg
                    poids = vpg(ipg)

                    !---------------------------------------------------------------------
                    !--- Calcul des fonctions d'interpolation et des derivees
                    call elem_B(vn,vb,detj,ksig(ipg,:),ielm(1))
                    detj1 = detj
                    deallocate(vb,vn)

                    call elem_B(vn,vb,detj,ksig(ipg,:),ielm(2))
                    detj2 = detj
                    deallocate(vb,vn)

                    !----------------------------------------------------------------------
                    !--- integration matrice vke
                    Ve = Ve + (detj1+detj2)*poids*epais
                enddo

                Vel(ie) = Ve

                Vg = (4./3.)*pi*(Dg/2.)**3.
                r = Ve/Vg

                deallocate(ksig,vpg,ielm,cote)

                !--------------------------------------------------------------------------
                !--- Resistance a la traction et cohesion
                if (resp) then

                    loi1(1) = 0.d0  ! b
                    loi1(2) = 0.d0  ! c
                    loi1(3) = 0.d0  ! moy
                    loi1(4) = 0.d0  ! ect
                    loi1(5) = resist(igfp)%loi

                    !----------------------------------------------------------------------
                    !--- Loi normale
                    if (loi1(5)==1) then

                        if (resist(igfp)%ipa/=0) then
                           loi1(3) = resist(igfp)%param(1)
                           loi1(4) = resist(igfp)%param(2)
                        else
                           call rossi_coef(loi1(3),loi1(4),r,fc,'RESI')
                        endif

                    !----------------------------------------------------------------------
                    !--- Loi de Weibull
                    elseif (loi1(5)==2) then

                        if (resist(igfp)%ipa/=0) then
                            loi1(1) = resist(igfp)%param(1)
                            loi1(2) = resist(igfp)%param(2)
                        else
                            !call interf_wblcoef(loi1(1),loi1(2),r,fc)

                            !Initialisations des paramètres
                            alphab = 0.d0
                            betab  = 0.d0
                            gammab = 0.d0
                            alphac = 0.d0
                            betac  = 0.d0
                            gammac = 0.d0

                            if ((typea=='EJQ4').or.(typea=='EJT6')) then
                                !Définition du paramètre b
                                alphab = 1.598e-04*fc**2-3.714e-02*fc+3.119e+00
                                betab  = 0.15
                                gammab = -9.782e-05*fc**2+3.156e-02*fc-0.8103e+00

                                !Définition du paramètre c
                                alphac = 1.605e-02*fc-2.79e+00
                                betac  = 0.11
                                gammac = 3.45

                            elseif (typea=='EJQ6') then
                                !Définition du paramètre b
                                alphab = 1.283e-04*fc**2-2.558e-02*fc+1.96e+00
                                betab  = 0.2
                                gammab = -1.214e-04*fc**2+2.811e-02*fc-0.1095e+00

                                !Définition du paramètre c
                                alphac = -7.569e-05*fc**2+2.04e-02*fc-4.989e+00
                                betac  = 0.06
                                gammac = 5.4

                            endif

                            loi1(1) = exp((alphab)*exp(-(betab)*log(r))+(gammab))
                            loi1(2) = exp((alphac)*exp(-(betac)*log(r))+(gammac))

                        endif

                    !----------------------------------------------------------------------
                    !--- Loi log-normale
                    elseif (loi1(5)==3) then

                       if (resist(igfp)%ipa/=0) then
                           loi1(3) = resist(igfp)%param(1)
                           loi1(4) = resist(igfp)%param(2)
                       else
                           call rossi_coef(loi1(3),loi1(4),r,fc,'RESI')
                       endif

                    endif

                    !----------------------------------------------------------------------
                    !--- Calcul de la resistance
                    Rt = distr_alea(loi1,1,ie0)

                    !----------------------------------------------------------------------
                    !--- Cas particulier des lois 105 et 108 (resistance prise égale a sa valeur moyenne)
                    if ((int(vprelg(kprop(ie),1))==105) .or. &
                    &    (int(vprelg(kprop(ie),1))==108))    then
                        if (loi1(5)==1) then
                            Rt(1) = loi1(3)
                        elseif (loi1(5)==2) then
                            Rt(1) = loi1(1) * gamma(1.d0 + (1.d0/loi1(2)))
                        elseif (loi1(5)==3) then
                            Rt(1) = loi1(3)
                        endif
                    endif

                    !----------------------------------------------------------------------
                    !--- Calcul de la cohesion
                    !C = 5.*Rt
                    id = 6
                    if (dime == 3) id = id-1
                    !----- Traitement particulier pour la cohesion (égale à fc/2, non aleatoire)
                    loi1(5) = 1
                    loi1(3) = vprelg(kprop0(ie),id+1)
                    !loi1(3) = fc/2.D0
                    !loi1(4) = 1.D-01 * loi1(3)
                    loi1(4) = 0.D0
                    C = distr_alea(loi1,1,ie0)
                    !----- Preparation pour trace graphique eventuel
                    RR(iefi) = Rt(1)
                    CC(iefi) = C(1)
                    iprobFI = .true.
                    iprobVO = .true.

                endif

                !--------------------------------------------------------------------------
                !--- Stockage resistance a la traction et cohesion
                vprelg0(nbg + ie0,:) = vprelg(kprop(ie),:)
                iloi = int(vprelg0(nbg + ie0,1))

                if (iloi == 100 .or. iloi == 101 .or. iloi == 102 .or. &
                &   iloi == 105 .or. iloi == 106 .or. iloi == 107 .or. &
                &   iloi == 108 .or. iloi == 158 ) then
                    id = 6
                    if (dime == 3) id = id-1
                    vprelg0(nbg + ie0,id)   = Rt(1)
                    vprelg0(nbg + ie0,id+1) = C(1)
                    !--- petite modification
                    if (iloi == 101) vprelg0(kprop(ie),id+3) = atan(Rt(1)/(C(1)+1.d-20))
                else
                    print*,'interf_distal (resi) : loi ',iloi ,' non probabilisee'
                    stop
                endif

                if ((.not.betp).and.(.not.fibp)) then
                    kprop(ie) = nbg + ie0
                endif

            endif

            !------------------------------------------------------------------------------
            !--- si l'element appartient a un groupe probabiliste (energie de fissuration du beton)
            if ((betp) .and. (typea=='EJQ4' .or. typea=='EJQ6' .or. typea=='EJT6')) then

                !--------------------------------------------------------------------------
                !--- Calcul de la section elementaire fissuree
                allocate(vpg(size(infele(ktypel(ie))%W)))
                allocate(ksig(size(infele(ktypel(ie))%Q,1),size(infele(ktypel(ie))%Q,2)))

                vpg = infele(ktypel(ie))%W
                ksig = infele(ktypel(ie))%Q
                npg = size(ksig,1)

                !--------------------------------------------------------------------------
                !--- Recuperation de l'epaisseur si 2D
                epais = 1.
                if (dime==2) then
                    !  En contraintes planes (CP) (position de l'epaisseur)
                    if (int(vprelg(kprop0(ie),2))==3) epais = vprelg(kprop0(ie),idmax-2)
                endif

                !--------------------------------------------------------------------------
                !--- Calcul de la surface elementaire Se
                Se = 0.d0
                do ipg = 1, npg
                    poids = vpg(ipg)
                    call elem_B(vn,vb,detj,ksig(ipg,:),ie)
                    deallocate(vb,vn)
                    Se = Se + detj*poids*epais
                enddo

                !--------------------------------------------------------------------------
                !--- Calcul du rapport de surfaces
                Sg = pi*(Dg/2.)**2.
                r = Se/Sg

                !--------------------------------------------------------------------------
                !--- Desallocations
                deallocate(ksig,vpg)

                allocate(vitrav(1))
                vitrav= find_vec(nrjbeto%num,kprop0(ie))
                igbt = vitrav(1)
                deallocate(vitrav)

                !--------------------------------------------------------------------------
                !--- distribution aleatoire d'energie de fissuration du beton
                if (.not.(resp)) then
                    ie0 = ie0 + 1
                endif
                ienb = ienb + 1

                loi2(1) = 0.d0  ! b
                loi2(2) = 0.d0  ! c
                loi2(3) = 0.d0  ! moy
                loi2(4) = 0.d0  ! ect
                loi2(5) = nrjbeto(igbt)%loi

                !--------------------------------------------------------------------------
                !--- verification compatibilite avec la loi
                id = 6
                if (dime == 3) id = id-1

                if ((iloi == 105).or.(iloi==107)) then  ! Attention test temporaire (JLT - 25/06/2014)
                    id1=4
                elseif ((iloi == 108).or.(iloi == 158)) then
                    id1 = 5
                else
                    print*,'interf_distal (energ fiss) : loi ',iloi ,' non probabilisee'
                    stop
                endif

                !--------------------------------------------------------------------------
                !--- Loi normale
                if (loi2(5)==1) then
                    if (nrjbeto(igbt)%ipa/=0) then
                        loi2(3) = nrjbeto(igbt)%param(1)
                        loi2(4) = nrjbeto(igbt)%param(2)
                    else
                            stop ' interf_distal (energ fiss) : possibilite non encore programmee pour loi normale '
                    endif

                !--------------------------------------------------------------------------
                !--- Loi log-normale
                elseif (loi2(5)==3) then
                    if (nrjbeto(igbt)%ipa/=0) then
                        loi2(3) = nrjbeto(igbt)%param(1)
                        loi2(4) = nrjbeto(igbt)%param(2)
                    else
                        loi2(3) = vprelg(kprop(ie),id+id1) ! Valeur moyenne de l'energie
                        loi2(4) = 101.0d0*dexp(-1.4d0*dlog(r+14.d0)) ! Ecart type sur l'energie
                        loi2(4) = loi2(4)*loi2(3)
                    endif
                else
                    stop ' interf_distal (energ fiss) : loi de probabilite non programmee'
                endif

                !--------------------------------------------------------------------------
                !--- Tirage aleatoire sur l'energie
                WB = distr_alea(loi2,1,ie0)

                !--------------------------------------------------------------------------
                !--- Preparation pour trace graphique eventuel
                WBB(ienb) = WB(1)
                iprobBT = .true.

                if (.not.(resp)) then
                    vprelg0(nbg + ie0,:) = vprelg(kprop(ie),:)
                endif
                iloi = int(vprelg0(nbg + ie0,1))

                !--------------------------------------------------------------------------
                !--- stockage de la valeur
                vprelg0(nbg + ie0,id+id1) = WB(1)

                if (.not.fibp) then
                    kprop(ie) = nbg + ie0
                endif

            endif

            !------------------------------------------------------------------------------
            !--- si l'element appartient a un groupe probabiliste (energie de pontage fibre/matrice)
            if ((fibp) .and. (typea=='EJQ4' .or. typea=='EJQ6' .or. typea=='EJT6')) then

                allocate(vitrav(1))
                vitrav= find_vec(nrjfibr%num,kprop0(ie))
                igfb = vitrav(1)
                deallocate(vitrav)

                if ((.not.(resp)).and.(.not.betp)) then
                    ie0 = ie0 + 1
                endif
                ienf = ienf + 1

                !--------------------------------------------------------------------------
                !--- Recuperation des parametres de la loi de distribution
                loi4(1) = 0.d0  ! b
                loi4(2) = 0.d0  ! c
                loi4(3) = 0.d0  ! moy
                loi4(4) = 0.d0  ! ect
                loi4(5) = nrjfibr(igfb)%loi

                !--------------------------------------------------------------------------
                !--- Loi normale
                if (loi4(5)==1) then
                    if (nrjfibr(igfb)%ipa/=0) then
                        loi4(3) = nrjfibr(igfb)%param(1)
                        loi4(4) = nrjfibr(igfb)%param(2)
                    else
                        stop 'interf_distal (pontage) : possibilite non encore programmee pour loi normale '
                    endif

                !--------------------------------------------------------------------------
                !--- Loi log-normale
                elseif (loi4(5)==3) then
                    if (nrjfibr(igfb)%ipa/=0) then
                        loi4(3) = nrjfibr(igfb)%param(1)
                        loi4(4) = nrjfibr(igfb)%param(2)
                    else
                        stop 'interf_distal (pontage) : possibilite non encore programmee pour loi log-normale '
                    endif
                else
                    stop 'interf_distal (pontage) : loi de distribution non programmee '
                endif

                !--------------------------------------------------------------------------
                !--- Tirage aleatoire sur l'energie
                WF = distr_alea(loi4,1,ie0)

                !--------------------------------------------------------------------------
                !--- Preparation pour trace graphique eventuel
                WFF(ienf) = WF(1)
                iprobFB = .true.

                if (.not.(resp).and.(.not.betp)) then
                    vprelg0(nbg + ie0,:) = vprelg(kprop(ie),:)
                endif
                iloi = int(vprelg0(nbg + ie0,1))

                if ((iloi==106).or.(iloi == 108).or.(iloi == 158)) then
                    id = 6
                    if (dime == 3) id = id-1
                    vprelg0(nbg + ie0,id+4) = WF(1)

                else
                    print*,'interf_distal (pontage fibre/matrice) : loi ',iloi ,' non probabilisee'
                    stop
                endif

                kprop(ie) = nbg + ie0

            endif

        enddo

        !----------------------------------------------------------------------------------
        !--- Traitement des elements massifs probabilistes pour beton ensuite
        do ie = 1, nelt

            !------------------------------------------------------------------------------
            !--- detection d'un element massif probabiliste (module d'young)
            modp = .false.
            if (allocated(young)) then
               if (count(young%num==kprop0(ie))==1) modp=.true.
            endif

            typeb = nomtype(ktypel(ie))

            if ((modp) .and. (typeb=='MBT3' .or. typeb=='MBT6' .or. typeb=='MTT4')) then

                allocate(vitrav(1))
                vitrav= find_vec(young%num,kprop0(ie))
                igma = vitrav(1)
                deallocate(vitrav)

                !--------------------------------------------------------------------------
                !--- deb fissure probabiliste
                ie0 = ie0 + 1
                iebe = iebe + 1

                !--------------------------------------------------------------------------
                !--- Recherche des elements d'interface voisins
                allocate(cote(size(infele(ktypel(ie))%face,1),size(infele(ktypel(ie))%face,2)))
                cote = infele(ktypel(ie))%face
                ncote = size(cote,1)
                np = size(cote,2)

                allocate(ielf(ncote))
                ielf = 0

                do icote = 1, ncote
                    do i = 1, nelt
                        if (i/=ie) then
                            k = 0 ; l = 0 ; m = 0 ; n = 0
                            do j = 1, size(kconec,2)
                                if (kconec(i,j)==kconec(ie,cote(icote,1))) k = 1
                                if (kconec(i,j)==kconec(ie,cote(icote,2))) l = 1

                                if (np==3) then
                                    if (kconec(i,j)==kconec(ie,cote(icote,3))) m = 1
                                endif

                                if (np==4) then
                                    if (kconec(i,j)==kconec(ie,cote(icote,3))) m = 1
                                    if (kconec(i,j)==kconec(ie,cote(icote,4))) n = 1
                                endif
                            enddo

                            if  ((k==1.and.l==1.and.np==2).or.(k==1.and.l==1.and.m==1.and.np==3) &
                            &.or.(k==1.and.l==1.and.m==1.and.n==1.and.np==4)) then

                                if (count(resist%num==kprop0(i))==1) ielf(icote) = i
                                exit

                            endif
                        endif
                    enddo
                enddo

                !--------------------------------------------------------------------------
                !--- Calcul du volume moyen massif : Vem = 1/n(Ve1 +...+ Ven)
                Vem = 0.d0
                nelf = count(ielf /= 0)

                do icote = 1, ncote
                    if (ielf(icote) /= 0) Vem = Vem + Vel(ielf(icote))
                enddo

                if (nelf /= 0) Vem = Vem/nelf
                Vg = (4./3.)*pi*(Dg/2.)**3
                r = Vem/Vg

                deallocate(ielf,cote)

                !--------------------------------------------------------------------------
                !--- Module d'Young (loi normale ou log-normale)
                loi3(1) = 0.d0  ! b
                loi3(2) = 0.d0  ! c
                loi3(3) = 0.d0  ! moy
                loi3(4) = 0.d0  ! ect
                loi3(5) = young(igma)%loi  ! num

                id=4
                if (dime==3) id=id-1

                !--------------------------------------------------------------------------
                !--- Loi normale
                if (loi3(5) == 1) then
                    if (young(igma)%ipa/=0) then
                        loi3(3) = young(igma)%param(1)
                        loi3(4) = young(igma)%param(2)
                    else
                        call rossi_coef(loi3(3),loi3(4),r,fc,'MODU',vprelg(kprop(ie),id))
                    endif
                    if (nelf==0) loi3(4)=0.d0

                !--------------------------------------------------------------------------
                !--- Loi de Weibull
                elseif (loi3(5) == 2) then
                    stop 'interf_distal : loi probabiliste non programmee pour les modules'

                !--------------------------------------------------------------------------
                !--- Loi log-normale
                elseif (loi3(5) == 3) then
                    if (young(igma)%ipa/=0) then
                        loi3(3) = young(igma)%param(1)
                        loi3(4) = young(igma)%param(2)
                    else
                        call rossi_coef(loi3(3),loi3(4),r,fc,'MODU',vprelg(kprop(ie),id))
                    endif

                    if (nelf==0) loi3(4)=0.d0
                endif

                !--------------------------------------------------------------------------
                !--- Calcul du module elastique
                E = distr_alea(loi3,1,ie0)

                !--------------------------------------------------------------------------
                !--- Preparation pour trace graphique eventuel
                EE(iebe) = E(1)
                iprobMY = .true.

                !--------------------------------------------------------------------------
                !--- Stockage modules elastiques
                vprelg0(nbg + ie0,:) = vprelg(kprop(ie),:)
                kprop(ie) = nbg + ie0

                id = 4
                if (dime == 3) id = id - 1
                vprelg0(kprop(ie),id) = E(1)

            endif
            !--- fin fissure probabiliste
            !------------------------------------------------------------------------------

        enddo

        deallocate(vprelg)

        call init_mat(vprelg,ie0+nbg,size(vprelg,2))

        !----------------------------------------------------------------------------------
        !--- Sorties
        vprelg(1:ie0+nbg,:) = vprelg0(1:ie0+nbg,:)

        !----------------------------------------------------------------------------------
        !--- Visualisation de la distribution des volumes
        if (iprobVO) then
            dPlay = 1
            if (dplay>=1) then
                print*,' '
                print*,'Maillage beton probabiliste :'
                print*,'-----------------------------'
            endif
            call distrib_volumes(Vel,Dg,50,dPlay)
        endif

        !----------------------------------------------------------------------------------
        !--- Visualisation des distributions de resistance et de cohesion
        if (iprobFI) then
            dPlay = 1
            if (dplay>=1) then
                print*,'Fissuration probabiliste :'
                print*,'--------------------------'
                print'(a10,e12.5,a16,e12.5)','- Dmax :',Dg,' - Resistance :',fc
                print'(a10,e12.5)','- Vg :',Vg
            endif
            call alea_trac(RR,loi1,200,dPlay,'resistance')
            call alea_trac(CC,loi1,200,dPlay,'cohesion beton')
        endif

        !----------------------------------------------------------------------------------
        !--- Visualisation des distributions de module elastique
        if (iprobMY) then
            dPlay = 1
            if (dplay>=1) then
                call alea_trac(EE,loi3,200,dPlay,'module')
            endif
        endif

        !----------------------------------------------------------------------------------
        !--- Visualisation des distributions d energie de fissuration
        if (iprobBT) then
            dPlay = 1
            if (dplay>=1) then
                print*,'Distribution d''energie (fissuration beton) :'
                print*,'---------------------------------------------'
            endif
            call alea_trac(WBB,loi2,200,dPlay,'energie (fissuration beton)')
        endif

        !----------------------------------------------------------------------------------
        !--- Visualisation des distributions d energie de pontage des fibres
        if (iprobFB) then
            dPlay = 1
            if (dplay>=1) then
                print*,'Distribution d''energie (pontage fibres/matrice) :'
                print*,'--------------------------------------------------'
            endif
            call alea_trac(WFF,loi4,200,dPlay,'energie (pontage fibres/matrice)')
        endif

        !call ecriture_distribution(RR,iefi,EE,iebe)   ! utilitaire.f

        !----------------------------------------------------------------------------------
        !--- Liberation memoire JLT
        deallocate(vprelg0,EE,RR,CC,WBB,WFF)

    endif

end subroutine interf_distal


!------------------------------------------------------------------------------------------------------
!> @author JL Tailhan    (version 1.0 - 2014)
!
!> @brief
!> Fonction pemettant le calcul de la chute de resistance pour la modélisation
!> de la propagation de fissure sous charge
!
!> @details
!> #### DESCRIPTION:
!>   Donnees en entree :\n
!>   - ipas    : numero du pas de temps\n\n
!
!>   Valeurs en sortie :\n
!>   - vprelg  : proprietes mecaniques au temps (n+1) apres\n
!>               prise en compte de la chute de resistance\n
!>               (prise en compte du fluage)\n
!------------------------------------------------------------------------------------------------------
subroutine interf_resichute(ipas) !,icalc_iter,iter)

    !**************************************************************************************
    !* DEFINITIONS
    !**************************************************************************************

    !--------------------------------------------------------------------------------------
    !--- Modules generaux
    use variables, only :   interf, nelt, nomtype, ktypel, kprop0, kprop, &
    &                       infele, vcont, kcont, dime, &
    !--- proprietes
    &                       vprelg, & !vprelg0, vprelg1, &
    !--- prise en compte de l'alea
    &                       resist, ietatpg, &
    !--- prise en compte de l'evolution temporelle de la resistance (propagation de fissure)
    &                       vdu, idmax, kloce, tps, dtps
    !&                      Dtm1, Dtm2, Dtrest, alp_Dt

    use lib_elem, only : elem_B, elem_kloce2, elem_hooke
    use initialisation, only : init_vec, init_mat
    use proprietes_aleatoires
    use math

    implicit none

    !--------------------------------------------------------------------------------------
    !--- Variables in/out
    integer, intent(in) :: ipas !, icalc_iter, iter

    !--------------------------------------------------------------------------------------
    !--- Variables locales
    real*8, dimension(:), allocatable   :: vpg, vdle, vdsi
    real*8, dimension(:,:), allocatable :: ksig, vsig, vn, vb, vh

    integer :: ie, nc, npg, idim, iloc, ipg, iloi, id, ndle, iedng
    real*8 :: vprel(idmax), temps(nelt), xi(nelt)
    real*8  :: detj, sigm, dRt, Rt   !, poids, ptot, Rt0, Rt1, dRc, Rc0, Rc1, Rc
    real*8  :: Dt, tmin, ximax, coef !, t, Dt0
    real*8  :: xa,xb                 !, ft, Lc, L0
    character(len=5) :: typel

    logical :: expl

    !**************************************************************************************
    !* CORPS DE PROCEDURE
    !**************************************************************************************

    !--------------------------------------------------------------------------------------
    coef = 1.d0
    !--------------------------------------------------------------------------------------
    expl = .true.

    temps = dtps(ipas)
    tmin = dtps(ipas)
    iedng = 0
    xi = 0.d0
    ximax = 0.d0

    xa = log(10.d0)/(-0.02552D0)
    xb = 0.9308D0

    !--------------------------------------------------------------------------------------
    !--- Dans le cas d'un calcul avec elements d'interface uniquement-!
    if (interf == 1) then

        !----------------------------------------------------------------------------------
        !--- Calcul de l'intervalle de temps (depend du pilotage par element le plus dangereux)
        !    Dt0 = dtps(ipas)
        !    if (icalc_iter==1) then
        !        if (iter==1) Dtrest = Dt0
        !        Dt = alp_Dt * Dtrest
        !    else
        !        if (iter==1) Dtrest = Dtm2 - Dtm1
        !        Dt = alp_DT * Dtrest
        !    endif
        !    if (iter==1) then
        !        Dtm2 = Dtrest
        !        Dtm1 = Dt
        !    endif
        !    if (ipas==1) then
        !        tps(ipas) = 0
        !    else
        !        tps(ipas) = tps(ipas-1) + Dt
        !    endif

        !----------------------------------------------------------------------------------
        !--- Traitement des elements d'interface béton
        do ie = 1, nelt
            typel = nomtype(ktypel(ie))

            if (allocated(resist)) then
                if ((count(resist%num==kprop0(ie))==1).and.&
                &   (typel=='EJQ4'.or.typel=='EJQ6'.or.typel=='EJT6').and.&
                &   (all(ietatpg(ie,:)==0))) then

                    !----------------------------------------------------------------------
                    !--- Recuperation des informations elementaires
                    call init_vec(vpg, size(infele(ktypel(ie))%W))
                    call init_mat(ksig, size(infele(ktypel(ie))%Q,1), size(infele(ktypel(ie))%Q,2))
                    nc   = infele(ktypel(ie))%ncgr
                    vpg  = infele(ktypel(ie))%W
                    ksig = infele(ktypel(ie))%Q
                    npg  = size(ksig,1)

                    !----------------------------------------------------------------------
                    !--- Recuperation des proprietes de l'element
                    vprel = vprelg(kprop(ie),1:idmax)

                    !----------------------------------------------------------------------
                    !--- localisation des contraintes elementaires dans vcont
                    idim = npg*nc;
                    iloc = kcont(ie)

                    !----------------------------------------------------------------------
                    !--- Recuperation de l'increment de deplacement (pour calcul de l'increment de contrainte)
                    allocate(vdle(ndle))
                    call elem_kloce2(ie, ndle)
                    vdle  = vdu(kloce(1:ndle))

                    !----------------------------------------------------------------------
                    !--- Recuperation des contraintes elementaires dans vcont (iteration precedante)
                    call init_mat(vsig, nc, npg)
                    call init_vec(vdsi, nc)
                    vsig  = reshape(vcont(iloc : (iloc+ idim - 1)),(/ nc, npg/))

                    !----------------------------------------------------------------------
                    !--- Calcul de la contrainte au centre de gravite de l'element
                    sigm = 0.d0
                    ! ptot = 0.d0
                    !
                    !-----------------------------------------------------------
                    !do ipg = 1, npg
                    !     poids = vpg(ipg)
                    !     !----- Calcul des fonctions d'interpolation et des derivees
                    !     call elem_B(vn,vb,detj,ksig(ipg,:),ie)

                    !     !----- calcul de la contrainte moyenne
                    !     !sigm = sigm + detj*poids*vsig(1,ipg)
                    !     ptot = ptot + detj*poids
                    !     deallocate(vb,vn)
                    !enddo
                    !-----------------------------------------------------------
                    !
                    !!sigm = sigm / ptot

                    ipg = 1
                    call elem_hooke(vh,nomtype(ktypel(ie)),vprel)
                    call elem_B(vn,vb,detj,ksig(ipg,:),ie)
                    vdsi = matmul(vh,matmul(vb,vdle))

                    sigm = vsig(1,1) + vdsi(1)

                    !L0 = 110.d0
                    !Lc = DSQRT(ptot)

                    deallocate(ksig,vpg,vsig,vdsi,vdle)
                    deallocate(vb,vn,vh)

                    !----------------------------------------------------------------------
                    !--- Recuperation des resistances au temps t(n)
                    iloi = int(vprelg(kprop(ie),1))
                    if (iloi == 100 .or. iloi == 101 .or. iloi == 106) then
                        id = 6
                        if (dime == 3) id = id-1
                        !Rt0 = vprelg0(kprop(ie),id)
                        Rt  = vprelg(kprop(ie),id)
                    else
                        print*,'FIDES_interf_resichu : la loi ',iloi ,' n''est pas une loi beton !'
                        stop
                    endif

                    !----------------------------------------------------------------------
                    !--- CALCUL EXPLICITE
                    if (expl) then

                        !------------------------------------------------------------------
                        !--- Calcul de la vitesse de chute de contrainte
                        !----- Loi numero 1
                        !xa = log(10.d0)/(0.2072D0)
                        !xb = exp(-(-2.414D0)*xa)
                        !----- Loi numero 2
                        if (sigm>=0) then
                            xi(ie) = min(sigm / RT , 1.d0)
                            !----- Loi numero 1
                            !dRt = -1.D0 * (RT0/xb) * ((xi/(1-xi))**xa)
                            !----- Loi numero 2
                            ! (longueur element / (hauteur - longueur de fissure))^2
                            dRt = -1.D0 * RT * exp(-xa*(xi(ie)-xb)) * coef ! 1/(Lc/L0)**2

                            Dt = (sigm - 1.001*RT) / dRt
                            Dt = max(Dt, 0.d0)
                            temps(ie) = Dt

                            if (tmin > Dt) then
                                tmin = Dt
                                iedng = ie
                            endif
                        endif

                    endif
                    !--- FIN CALCUL EXPLICITE
                    !----------------------------------------------------------------------

                endif
            endif
        enddo

        if (tmin>dtps(ipas)) tmin=dtps(ipas)

        if (ipas==1) then
            tps(ipas) = 0
        else
            tps(ipas) = tps(ipas-1) + tmin
        endif

        !----------------------------------------------------------------------------------
        !--- Chute de resistance
        !    do ie = 1, nelt
        ie = iedng
        if (allocated(resist)) then
            if ((count(resist%num==kprop0(ie))==1).and.&
            &(typel=='EJQ4'.or.typel=='EJQ6'.or.typel=='EJT6').and.&
            &(all(ietatpg(ie,:)==0))) then

                id = 6
                if (dime == 3) id = id-1
                !Rt0 = vprelg0(kprop(ie),id)
                Rt  =  vprelg(kprop(ie),id)

                dRt = -1.D0 * RT * exp(-xa*(xi(ie)-xb)) * coef
                Rt = Rt + tmin * dRt
                !if (ie==iedng) Rt = 0.d0

                if (sigm > Rt) Rt=0.d0

                iloi = int(vprelg(kprop(ie),1))
                if (iloi == 100 .or. iloi == 101.or. iloi == 106) then
                    id = 6
                    if (dime == 3) id = id-1
                    vprelg(kprop(ie),id)   = Rt
                else
                    print*,'FIDES_interf_resichu : la loi ',iloi ,' n''est pas une loi beton !'
                    stop
                endif

            endif
        endif
    !enddo

    endif

end subroutine interf_resichute


!------------------------------------------------------------------------------------------------------
!> @author JL Tailhan    (version 1.0 - 2014)
!
!> @brief
!> Calcul des contraintes pour la loi de contact
!
!> @details
!> #### DESCRIPTION:
!>  Calcul des contrainte normale et contrainte tangente dans le respect de la loi de Coulomb
!------------------------------------------------------------------------------------------------------
subroutine interf_penalModul(vh,vb,vdle,etat,ie)

    !**************************************************************************************
    !* DEFINITIONS
    !**************************************************************************************

    !--------------------------------------------------------------------------------------
    !--- Modules generaux
    use variables, only : dime

    implicit none

    !--------------------------------------------------------------------------------------
    !--- Variables in/out
    real*8, dimension(:), intent(in) :: vdle
    real*8, dimension(:,:), intent(inout) :: vh
    real*8, dimension(:,:), intent(in) :: vb
    integer, intent(in) :: etat,ie

    !--------------------------------------------------------------------------------------
    !--- Variables locales
    real*8, dimension(size(vb,1)) :: deprel
    real*8 :: depreln,deprelt,dep,vhmax,ouvermax

    !**************************************************************************************
    !* CORPS DE PROCEDURE
    !**************************************************************************************

    ouvermax=1.D11

    !--------------------------------------------------------------------------------------
    !--- Calcul des deplacements relatifs
    if (dime==2) then
        deprel = matmul(vb,vdle)
        depreln = deprel(1)                          ! Deplacement normal
        deprelt = deprel(2)                          ! Deplacement tangentiel

    elseif (dime==3) then
        deprel = matmul(vb,vdle)
        depreln = deprel(1)                          ! Deplacement normal
        deprelt = sqrt(deprel(2)**2 + deprel(3)**2)  ! Deplacement tangentiel

    endif

    dep = sqrt(depreln**2+deprelt**2)
    vhmax = 1.d11

    if (etat==1) then
        if (dep > ouvermax) then
            print*,'ie',ie,'element ouvert dans l''espace',dep,ouvermax
            vh(1,1) = vhmax*(dep-ouvermax)/(10.*ouvermax)
            if (vh(1,1)>vhmax) vh(1,1) = vhmax
            vh(2,2) = vhmax*(dep-ouvermax)/(10.*ouvermax)
            if (vh(2,2)>vhmax) vh(2,2) = vhmax
            print*,vh(1,1),vh(2,2)
        endif
    else if (etat==2) then
        if (dep > ouvermax) then
            print*,'ie',ie,'element cisaille dans l''espace',dep,ouvermax
            vh(2,2) = vhmax*(dep-ouvermax)/(10.*ouvermax)
            if (vh(2,2)>vhmax) vh(2,2) = vhmax
            print*,vh(1,1),vh(2,2)
        endif
    endif

end subroutine interf_penalModul

!------------------------------------------------------------------------------------------------------

end module element_interface

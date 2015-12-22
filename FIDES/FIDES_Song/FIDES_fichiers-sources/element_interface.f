module element_interface

!********************************************************!
!     Gestion du calcul des elements de l'interface      !
!********************************************************!

contains

!********************************************************!

subroutine interf_sig(vhep,vsig,vnle,vh,vb,vnl0,vprel,vdle,iloi,ie,ipg)

!********************************************************!
!     Calcul des contraintes pour la loi de contact      !
!      (contrainte normale et contrainte tangente )      !
!           Respect du critere Mohr-Coulomb              !
!********************************************************!
    use variables, only : vsol, dime, ietatpg, ktypel, nomtype, imet, vdle0
    use lib_elem, only : elem_hooke
    use initialisation, only : init_vec, init_mat
    use math, only : vecpro
    implicit none
	
    ! Variables IN
    real*8, dimension(:), intent(in) :: vnl0, vdle, vprel
    real*8, dimension(:,:), intent(inout), allocatable :: vh
    real*8, dimension(:,:), intent(in) :: vb
    integer, intent(in) :: ie, ipg, iloi

    ! Variables OUT
    real*8, dimension(:), intent(inout) :: vsig, vnle
    real*8, dimension(:,:), intent(out) :: vhep

    ! Quantites principales
    real*8, dimension(size(vh,1)) :: vdfdsig, vdgdsig
    real*8, dimension(size(vb,1)) :: vepspl, vsigm, vsigi, vdsig
    real*8, dimension(size(vh,1),size(vh,2)) :: vh1, vh2
    real*8 :: vcrit, dlam
    integer :: it
    character(len=5) :: calcu

!********************************************************!

    !----- Recuperation des grandeurs de non-linearite
    vepspl = vnl0(1:size(vb,1))            ! deformations plastiques initiales

    !----- Calcul de sigma + dsigma
    vsigm = matmul(vh,(matmul(vb,vdle)-vepspl))
    vdsig = matmul(vh,(matmul(vb,vdle0)-vepspl))   

    !----- Calcul des contraintes verifiant le critere
    calcu = 'D1F'
    call interf_crit(iloi,vprel,vsigm,vdsig,ie,ipg,vcrit,calcu,vdle,vb,vdfdsig,vdgdsig)
    calcu = 'D0F'

    !----- Calcul du multiplicateur plastique (algo semi-implicite)
    it = 0
    vsigi = vsigm
    dlam = 0.

    do while (vcrit > 1.d-30)

	    dlam = dlam + vcrit/dot_product(vdfdsig,matmul(vh,vdgdsig))
	    vsigm = vsigi - dlam*matmul(vh,vdgdsig)
        call interf_crit(iloi,vprel,vsigm,vdsig,ie,ipg,vcrit,calcu,vdle,vb)
	    it = it + 1
	    if (it > 50) exit

    end do
    
    !----- Calcul du module elastoplastique
    if (imet == 1) then
        call elem_hooke(vh,nomtype(ktypel(ie)),vprel)
        call interf_modul(ie,ipg,vh,vprel)
        vhep = vh
    elseif (imet == 2) then
        if (ietatpg(ie,ipg)==1) then
            vhep = 0.0d0
                 
        elseif (ietatpg(ie,ipg)==2 .or. ietatpg(ie,ipg)==3) then
            call elem_hooke(vh,nomtype(ktypel(ie)),vprel)
            vhep = vh - vecpro(matmul(vh,vdgdsig),matmul(vdfdsig,vh))/dot_product(vdfdsig,matmul(vh,vdgdsig))
                                     
        end if
    else
        vhep = vh
    end if

    vepspl = vepspl + dlam*vdgdsig

    !----- Sorties :
    vsig = vsigm
    vnle = vepspl

end subroutine interf_sig

!********************************************************!

subroutine interf_crit(iloi,vprel,vsigm,vdsig,ie,ipg,vcrit,calcu,vdle,vb,vdfdsig,vdgdsig)

!********************************************************!
!     Calcul des criteres de rupture et des derivees     !
!               pour l'element d'interface               !
!     sorties :  vcrit, vdfdsig, vdgdsig                 !
!********************************************************!

    use variables, only : dime, kprop0, ietatpg, &
            & iedng, nelt, pi, ktypel, infele, nomtype, ipas
    use initialisation, only : init_vec, init_mat
    use lib_elem, only : elem_hooke
    implicit none

    ! Variables IN
    real*8, dimension(:), intent(in) :: vsigm, vdsig, vprel, vdle
    real*8, dimension(:,:), intent(in) :: vb
    
    integer, intent(in) :: ie, iloi, ipg
    character(len=5), intent(in), optional :: calcu

    ! Variables OUT
    real*8, dimension(:), intent(out), optional :: vdfdsig, vdgdsig
    real*8, intent(out) :: vcrit

    ! Quantites globales
    character(len=5) :: typ
    real*8 :: sigma, tau, tau1, tau2, tau0, tau01, tau02, dtau1, dtau2, dtau, RT, C, phi, psi
    real*8 :: vcrit1, vcrit2, vcrit2a, vcrit2b, vcrit3, vcrit4, kn, kt, kt1, kt2, A
    integer :: id, ipc, i, npg
    logical :: ipremf, ideriv, iendo
    real*8, dimension(:,:), allocatable :: vh
    real*8, dimension(size(vb,1)) :: deprel
    real*8 :: depcrin, depcrit, depreln, deprelt, deprelt1, deprelt2, depreln0, deprelt0

!********************************************************!
    
    !----- Test(s) de compatibilite des arguments d'entree
    ideriv = .false.
    typ = nomtype(ktypel(ie))
    
    if (calcu == 'D1F') then
        ideriv = .true.
        vdfdsig = 0.0d0
        vdgdsig = 0.0d0
        
        call elem_hooke(vh,typ,vprel)
        call interf_modul(ie,ipg,vh,vprel)
    end if
    
    !----- Recuperation des parametres de la loi
    id = 6
    if (dime == 3) id = id-1
    C  = vprel(id+1)
    phi = pi*vprel(id+2)/180.           ! Angle de frottement
    psi = pi*vprel(id+3)/180.           ! Angle de dilatance
    RT = vprel(id)                      ! RT

    !----- Recuperation des contraintes normale(s) et tangentielle(s)
    if (dime == 2) then
        sigma = vsigm(1)                  ! Contrainte normale
        tau   = vsigm(2)                  ! Contrainte tangente

    elseif (dime == 3) then 
        sigma  = vsigm(1)                ! Contrainte normale
        tau1   = vsigm(2)                ! Contrainte tangente 1
        tau2   = vsigm(3)                ! Contrainte tangente 2
        tau = sqrt(tau1**2+tau2**2)

        if (iedng==ie) then
	        if (ideriv) then
	        ! Contraintes du pas precedent		
		        tau01  = vsigm(2)-vdsig(2)   ! Contrainte tangente 1
		        tau02  = vsigm(3)-vdsig(3)   ! Contrainte tangente 2
		
		        ! Increment de contraintes
		        dtau1 = vdsig(2)
		        dtau2 = vdsig(3)
				
		        tau0  = sqrt(tau01**2 + tau02**2)
		        dtau  = sqrt(dtau1**2 + dtau2**2)

               if ((abs(tau1)>=abs(tau01)).or.(abs(tau2)>=abs(tau02))) then
                  tau = tau0 + dtau     ! charge
               elseif ((abs(tau1)<abs(tau01)).or.(abs(tau2)<abs(tau02))) then
                  tau = tau0 - dtau    ! decharge
               end if                  
            end if
        end if           
    end if

    !----- Gestion des etats de l'interface
    if (typ=='EJQ4'.or.typ=='EJQ6') then
        ipc = 1
    elseif (typ=='EJT6') then
        ipc = 4
    elseif (typ=='EJT8') then
        ipc = 5
    end if
    
    !----- Detection premiere fissuration
    ipremf = .false.
    npg = size(infele(ktypel(ie))%Q,1)        
    if (ipg==ipc .and. count(ietatpg(ie,:)==0)==npg) ipremf = .true.   
   
    select case(iloi)
    case(100)          !----- Loi complexe : tension cut-off et Coulomb
        
        if (RT > (C/tan(phi))) RT = C/tan(phi)
                
        if (ietatpg(ie,ipg) /=0 ) then     ! Element fissure ou element referme
            RT=0.0d0
            C=0.0d0
        end if
        
        !----- Criteres de rupture de l'element d'interface
        vcrit1 = sigma - RT                               ! Tension cut-off
        vcrit2 = abs(tau) + sigma * tan(phi) - C          ! Mohr-Coulomb
        vcrit3 = abs(tau) + sigma * tan(psi-pi/2.)
        
        vcrit = -1.D0

        if ((iedng==ie) .or. (ietatpg(ie,ipg) /=0 )) then                 ! strategie 1
        !if ((iedng==ie) .or. (ipremf.eqv. .false.)) then                 ! strategie 2
        
            if ((vcrit3 < 0.d0) .or. (ietatpg(ie,ipg) < 0.d0)) then

                ! Critere en ouverture pure
                if ((vcrit1> 0.d0) .or. (vcrit2> 0.d0)) then
                    vcrit = sqrt(sigma**2 + tau**2)

                    if (ideriv .eqv. .true.) then
                        !if (ietatpg(ie,ipg) >= 0)  ietatpg(ie,ipg) = 1     ! strategies 2 4
                        if (ietatpg(ie,ipg) > 0)  ietatpg(ie,ipg) = 1       ! strategie 1
                        !if (iedng==ie) ietatpg(ie,:) = 1                    ! strategie 1
                        if ((iedng==ie) .and. (ipremf)) ietatpg(ie,:) = 1   ! strategie 1

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
                        end if
                    end if
                end if
            else

                ! Critere en cisaillement
                if ((vcrit2 > 0.d0) .or. (vcrit1 > 0.d0)) then
                    vcrit = abs(tau) + sigma * tan(phi)

                    if (ideriv) then
                        ietatpg(ie,ipg) = 2                                 ! strategies 1 2 3 4
                        !if (iedng==ie) ietatpg(ie,:) = 2                ! strategie 1
                        if ((iedng==ie) .and. (ipremf)) ietatpg(ie,:) = 2   ! strategie 1

                        if (dime == 2) then
                            ! dans le repere local...
                            vdfdsig = (/ tan(phi), tau/abs(tau) /)
                            vdgdsig = (/ tan(psi), tau/abs(tau) /)

                        elseif (dime == 3) then
                            vdfdsig = (/ tan(phi), tau1/abs(tau), tau2/abs(tau) /)
                            vdgdsig = (/ tan(psi), tau1/abs(tau), tau2/abs(tau) /)
                        end if
                    end if
                end if
            end if
        end if                                                          ! strategies 1 et 2

    case(101)     !----- Loi complexe : tension cut-off + Tresca puis Coulomb

        if (RT > (C/tan(phi))) RT = C/tan(phi)
        
        !psi = atan(RT/(C+1.d-20)) dans la loi 101 probabilisee psi est calcule dans distal
        if (ietatpg(ie,ipg) /=0 ) then     ! Element fissure ou element referme
            RT=0.0d0
            !C=0.d0
        end if

        !----- Criteres de rupture de l'element d'interface
        vcrit1  = sigma - RT                          ! Tension cut-off
        vcrit2a = abs(tau)  - C                       ! Tresca
        vcrit2b = abs(tau) + sigma * tan(phi)         ! Mohr-Coulomb
        vcrit3  = abs(tau) + sigma * tan(psi-pi/2.)
        A = C*(1-(tan(psi-pi/2.)/tan(phi)))
        vcrit4  = vcrit3 - A
        
        vcrit = -1.D0

        if ((iedng==ie) .or. (ietatpg(ie,ipg) /=0 )) then       ! strateg 1
        !if ((iedng==ie) .or. ( ipremf .eqv. .false. )) then    ! strateg 2
        
           if ((vcrit3 < 0.d0) .or. (ietatpg(ie,ipg) < 0.d0)) then
           
               vcrit2=vcrit2a
               if (ietatpg(ie,ipg) /=0 ) vcrit2 = vcrit2b

               ! Critere en ouverture pure
               if ((vcrit1> 0.d0) .or. (vcrit2> 0.d0)) then
                   vcrit = sqrt(sigma**2 + tau**2)

                   if (ideriv) then

                       !if (ietatpg(ie,ipg) >= 0)  ietatpg(ie,ipg) = 1   ! strateg 2 4
                       if (ietatpg(ie,ipg) > 0)  ietatpg(ie,ipg) = 1     ! strateg 1 -deb
                       !if (iedng==ie) ietatpg(ie,:) = 1                 ! strateg 1 -fin
                       if ((iedng==ie) .and. (ipremf)) ietatpg(ie,:) = 1 ! strateg 1 -fin
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
                       end if
                   end if
               end if

           elseif ((vcrit3 > 0.d0) .and. (vcrit4 < 0.d0)) then
           
               vcrit2=vcrit2a
               if (ietatpg(ie,ipg) /=0 ) vcrit2 = vcrit2b

                   ! Critere en cisaillement
                   if ((vcrit2 > 0.d0) .or. (vcrit1 > 0.d0)) then
                   
                       vcrit = abs(tau) + sigma * tan(phi)

                       if (ideriv) then
                           ietatpg(ie,ipg) = 2    ! strateg 1 2 3 4

                           !if (iedng==ie) ietatpg(ie,:) = 2     ! strateg 1
                           if ((iedng==ie) .and. (ipremf))  ietatpg(ie,:) = 2   ! strateg 1

                           if (dime == 2) then
                               ! dans le repere local...
                               vdfdsig = (/ tan(phi), tau/abs(tau) /)
                               if (vcrit2==vcrit2a) vdfdsig = (/ 0d0, tau/abs(tau) /)
                               vdgdsig = (/ tan(psi), tau/abs(tau) /)

                           elseif (dime == 3) then
                               vdfdsig = (/ tan(phi), tau1/abs(tau), tau2/abs(tau) /)
                               if (vcrit2==vcrit2a) vdfdsig = (/ 0d0, tau1/abs(tau), tau2/abs(tau) /)
                               vdgdsig = (/ tan(psi), tau1/abs(tau), tau2/abs(tau) /)    
                           end if
                       end if
                   end if

               else

                   vcrit2=vcrit2a
                   ! Critere de Tresca
                   if ((vcrit2 > 0.d0) .or. (vcrit1 > 0.d0)) then  
                       vcrit = abs(tau) - C

                       if (ideriv) then
                           ietatpg(ie,ipg) = 2    ! strateg 1 2 3 4
                           !if (iedng==ie) ietatpg(ie,:) = 2      ! strateg 1 -fin
                           if ((iedng==ie) .and. (ipremf)) ietatpg(ie,:) = 2    ! strateg 1 -fin

                           if (dime == 2) then
                               ! dans le repere local...
                               vdfdsig = (/ 0.D0    , tau/abs(tau) /)
                               vdgdsig = (/ tan(psi), tau/abs(tau) /)

                           elseif (dime == 3) then
                               vdfdsig = (/ 0.D0    , tau1/abs(tau), tau2/abs(tau) /)
                               vdgdsig = (/ tan(psi), tau1/abs(tau), tau2/abs(tau) /)

                           end if
                       end if
                   end if
               end if
           end if ! strateg 1 2

        case(103)   ! Elastoplasticite Coulomb--Tresca

        psi = atan(RT/C)
        if (ietatpg(ie,ipg) /=0 ) then     ! Element fissure ou element referme
            RT=0.0d0
        end if

        !----- Criteres de rupture de l'element d'interface
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

                ! Critere en ouverture pure
                if ((vcrit1> 0.d0) .or. (vcrit2> 0.d0)) then
                    vcrit = sqrt(sigma**2 + tau**2)

                    if (ideriv) then
                        if (ietatpg(ie,ipg) > 0)  ietatpg(ie,ipg) = 1

                        if ((iedng==ie) .and. (ipremf .eqv. .true.)) then
                            ietatpg(ie,:) = 1
                        end if
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
                        end if
                    end if
                end if

            elseif ((vcrit3 > 0.d0) .and. (vcrit4 < 0.d0)) then

                if ((ietatpg(ie,ipg)/=0) .and. (ietatpg(ie,ipg)/=3)) then

                   vcrit2 = vcrit2b

                   ! Critere en cisaillement
                   if ((vcrit2 > 0.d0) .or. (vcrit1 > 0.d0)) then
                       vcrit = abs(tau) + sigma * tan(phi)

                       if (ideriv .eqv. .true.) then
                           ietatpg(ie,ipg) = 2

                           if ((iedng==ie) .and. (ipremf .eqv. .true.)) then
                               ietatpg(ie,:) = 2
                           end if

                           if (dime == 2) then
                               ! dans le repere local...
                               vdfdsig = (/ tan(phi), tau/abs(tau) /)
                               vdgdsig = (/ tan(psi), tau/abs(tau) /)
 
                           elseif (dime == 3) then
                               vdfdsig = (/tan(phi), tau1/abs(tau), tau2/abs(tau) /)
                               vdgdsig = (/tan(psi), tau1/abs(tau), tau2/abs(tau) /)
                           end if
                       end if
                   end if

            ! Elastoplastique
            !--------------------------------------------------
            ! Si l'element est vierge et en elastoplastique
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
                            end if

                            if (dime == 2) then
                                ! dans le repere local...
                                vdfdsig = (/ 0.D0    , tau/abs(tau) /)
                                vdgdsig = (/ tan(psi), tau/abs(tau) /)

                            elseif (dime == 3) then
                                vdfdsig = (/0.D0    , tau1/abs(tau), tau2/abs(tau) /)
                                vdgdsig = (/tan(psi), tau1/abs(tau), tau2/abs(tau) /)
                            end if

                 ! Elastoplastique
                     
                        elseif ((ietatpg(ie,ipg)==0) .or. (ietatpg(ie,ipg)==3)) then
                            ietatpg(ie,ipg) = 3
                            if ((iedng==ie) .and. (ipremf .eqv. .true.)) then
                                ietatpg(ie,:) = 3
                            end if
                            
                            if (dime == 2) then
                                    ! dans le repere local...
                                    vdfdsig = (/0.d0 ,   tau/abs(tau) /)
                                    vdgdsig = (/0.d0 ,   1.d0/)

                            elseif (dime == 3) then
                                    vdfdsig = (/0.d0 ,   tau1/abs(tau) , tau2/abs(tau)/)
                                    vdgdsig = (/0.d0 ,   1.d0    ,  1.d0/)
                            end if

                         endif

                         !-----------------------------

                    end if
                end if
            end if
        end if
   
    case default
        stop 'interf_crit : loi non implantee pour l''element d''interface'
    end select

end subroutine interf_crit

!********************************************************!

subroutine interf_sig_endo(vhep,vsig,vnle,vh,vb,vnl0,vprel,vdle,iloi,ie,ipg)

!********************************************************!
!     Calcul des contraintes pour la loi de contact      !
!      (contrainte normale et contrainte tangente )      !
!           Respect du critere Mohr-Coulomb              !
!********************************************************!
    use variables, only : vsol, dime, ietatpg, ktypel, infele, nomtype, endo, &
                        & ieendo, iedng, ipas, vdle0
    use lib_elem, only : elem_hooke
    implicit none
    
    ! Variables IN
    real*8, dimension(:), intent(in) :: vnl0, vdle, vprel
    real*8, dimension(:,:), intent(inout), allocatable :: vh
    real*8, dimension(:,:), intent(in) :: vb
    integer, intent(in) :: ie, ipg, iloi

    ! Variables OUT
    real*8, dimension(:), intent(inout) :: vsig, vnle
    real*8, dimension(:,:), intent(out) :: vhep

    ! Quantites principales
    character(len=5) :: typ
    real*8, dimension(size(vb,1)) :: vepspl, vsigm, deprel, deprel0, ddeprel
    real*8, dimension(size(vh,1),size(vh,2)) :: vhpen, vhpetit, vhep0
    real*8 :: depcrin, depcrit, depreln, deprelt, deprelt0, depreln0, depret0, ddeprelt
    real*8 :: RT, C, D, D0, penal, alpha , Dn, Dt
    integer :: id, ipc, npg

!********************************************************!

    !----- Recuperation des grandeurs de non-linearite
    vepspl = vnl0(1:size(vb,1))            ! deformations plastiques initiales
    vhep = vh
    
    !---- Recuperation du point d'integration central en fonction du type d'element
    typ = nomtype(ktypel(ie))
    if (typ=='EJQ4'.or.typ=='EJQ6') then
        ipc = 1
    elseif (typ=='EJT6') then
        ipc = 4
    elseif (typ=='EJT8') then
        ipc = 5
    end if

    select case(iloi)
        case(102)          !----- Loi complexe : tension cut-off et Coulomb
            !----- Recuperation des parametres de la loi
            id = 6
            if (dime == 3) id = id-1
            C  = vprel(id+1)
            RT = vprel(id)                      ! RT
                
            !----- Deplacement relatif critique en traction et en cisaillement
            depcrin = vprel(id+4)
            depcrit = vprel(id+5)
        
            !----- Calcul du deplacement relatif
            if (dime==2) then
               deprel = matmul(vb,vdle)
               depreln = deprel(1)
               deprelt = deprel(2)        ! Deplacement tangentiel 1
               
            elseif (dime==3) then
               
               deprel = matmul(vb,vdle)
               
               ! Deplacement normal
               depreln = deprel(1)   
               
               ! Deplacement tangentiel 1, 2
               deprelt = sqrt(deprel(2)**2 + deprel(3)**2)
               
            end if

            penal = 0.d0
            if (depreln <= 0.d0) penal = 1.d0
            if (depreln < 0.d0) depreln = 0.d0
        
            Dn = 0.d0
            Dt = 0.d0
        
            alpha = 0.d0
                     
            call elem_hooke(vh,nomtype(ktypel(ie)),vprel)
                             
            depreln0 = RT/vh(1,1)
            deprelt0 = C/vh(2,2)
            if (depcrin < depreln0) depcrin=depreln0
            if (depcrit < deprelt0) depcrit=deprelt0

            if ((depreln >= depreln0) .and. (depreln < depcrin)) Dn = 1.d0 - depreln0/depreln
            if ((abs(deprelt) >= deprelt0) .and. (abs(deprelt) < depcrit)) Dt = 1.d0 - deprelt0/abs(deprelt)  
                       
            if (iedng==ie .or. ietatpg(ie,ipg) == 1 .or. ietatpg(ie,ipg) == 2) then
                if (abs(deprelt) >= depcrit) Dt = 1.d0
                if (depreln >= depcrin) Dn = 1.d0
            end if
        
            D  = alpha * Dn + (1.-alpha) * Dt
            D0 = endo(ie,ipg)
            D  = max(D,D0)
            D  = min(1.d0,D)

            !if ((ipg==ipc) .and. (all(ietatpg(ie,:)==0)) .and. (D > 0)) ietatpg(ie,:) = 3
            !if ((ipg==ipc) .and. (all(ietatpg(ie,:)==3)) .and. (Dt == 1.d0)) ietatpg(ie,:) = 2               
            !if ((ipg==ipc) .and. (all(ietatpg(ie,:)==3)) .and. (Dn == 1.d0)) ietatpg(ie,:) = 1

           if (dime==2) then
                  if ((ipg==ipc).and.(ietatpg(ie,1)==0).and.(ietatpg(ie,2)==0).and.(ietatpg(ie,3)==0)) then
                    if (D > 0) ietatpg(ie,:) = 3            
                  end if
                
                  if ((ipg==ipc).and.(ietatpg(ie,1)==3).and.(ietatpg(ie,2)==3).and.(ietatpg(ie,3)==3)) then
                    if (Dt == 1.d0) ietatpg(ie,:) = 2
                  end if
                
                  if ((ipg==ipc).and.(ietatpg(ie,1)==3).and.(ietatpg(ie,2)==3).and.(ietatpg(ie,3)==3)) then
                    if (Dn == 1.d0) ietatpg(ie,:) = 1
                  end if
                  
            elseif (dime==3) then

                  if ((ipg==ipc) .and. all(ietatpg(ie,:)==0)) then
                    if (D > 0) ietatpg(ie,:) = 3                                       
                  end if

                  if ((ipg==ipc) .and. all(ietatpg(ie,:)==3)) then
                    if (Dt == 1.d0) ietatpg(ie,:) = 2
                  end if

                  if ((ipg==ipc) .and. all(ietatpg(ie,:)==3)) then
                    if (Dn == 1.d0) ietatpg(ie,:) = 1
                  end if
                             
            end if

            if (ietatpg(ie,ipg)==0) D = 0.d0
            if ((ietatpg(ie,ipg)==1) .or. (ietatpg(ie,ipg)==2)) D = 1.d0

            if (dime==2) then
               vhpen(1,:) = D*penal*(/vh(1,1), 0.d0/)
               vhpen(2,:) = D*penal*(/0.d0   , 0.d0/)
            elseif (dime==3) then
               vhpen(1,:) = D*penal*(/vh(1,1), 0.d0, 0.d0/)
               vhpen(2,:) = D*penal*(/0.d0   , 0.d0, 0.d0/)
               vhpen(3,:) = D*penal*(/0.d0   , 0.d0, 0.d0/)
            end if

            vhpetit = vh/1.e15
            vhep = (1-D)*vh + vhpen + vhpetit

            !----- Calcul des contraintes
            vsigm = matmul(vhep,matmul(vb,vdle))
            vepspl = (-D*matmul(vb,vdle))

            !----- Sorties :
            vsig = vsigm
            vnle = vepspl
            endo(ie,ipg) = D
           
        case default
            stop 'interf_endo : cas non implante'

    end select

end subroutine interf_sig_endo

!********************************************************!

subroutine interf_pilot(imetpilo, vsol, vduI, vduII, dlam, alphai, MOTi,elemfiss)

!********************************************************!
!  Pilotage du calcul : gestion du facteur de chargement !
!********************************************************!

    use variables, only : dime, nelt, idmax, ktypel, nomtype, infele, &
        & vprelg, kprop, kprop0, ietatpg, iedngi, ieendo, interf, pi, ipas, kloce
    use lib_elem, only : elem_B, elem_kloce2, elem_hooke
    use math
	implicit none

	! Variables IN
	real*8, dimension(:), intent(in) :: vduI, vduII, vsol
	integer, intent(in) :: imetpilo
	real*8, intent(in) :: dlam	

    ! Variables IN-OUT
    real*8, intent(out) :: alphai
    integer, intent(inout) :: elemfiss
    character(len=8), intent(inout) :: MOTi 

    ! Quantites principales
    real*8, dimension(:), allocatable :: vdle, vdl0, vsi0, vdsi, &
                    & alpg, vdlI, vdlII
    real*8, dimension(:,:), allocatable :: vn, vb, vh, ksig

    real*8 :: deprel(dime), deprel0(dime), vprel(idmax), alph(nelt)
    real*8 :: alpc, sign0, C, RT, phi, psi, sigma, tau, tau0, dtau, vdsig, vsig0, &
        & vcrit1, vcrit2, vcrit3, detj, signe, coef, alpgo, alpgc
    real*8 :: depcrin, depcrit, depreln, deprelt, deprelt1, deprelt2, depreln0, deprelt0
    real*8 :: vala, valb, valc, delta, alp1, alp2

    character(len=5) :: typ
    character(len=8) :: MOT, MOT1, MOT2
    character(len=8), dimension(nelt) :: rupt
    integer :: i, ie, ipg, npg, id, iloi, ipc, ndle
    
!----- Switch en fonction de la methode de pilotage

select case(imetpilo)
    
    case(1)
        !----- Pilotage sur element le plus dangereux par recalcul de vduI

        !- Pour le cas des elements d'interface seulement
        if (interf==1) then
            
            !- Initialisation
            alph = 1.d0
            iedngi = 0
            ieendo = 0

            !- Boucle sur les elements
            do ie = 1, nelt

                MOT = '  '
                MOT1 = '  '
                MOT2 = '  '
                typ = nomtype(ktypel(ie))
                
                !----- Pour les elements d'interface seulement
                if ((typ=='EJQ4'.or.typ=='EJQ6'.or.typ=='EJT6'.or.typ=='EJT8') &
                   & .and. (all(ietatpg(ie,:)==0) .or. all(ietatpg(ie,:)==3))) then

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
                    if (typ=='EJQ4'.or.typ=='EJQ6') then
                       ipc = 1
                    elseif (typ=='EJT6') then
                       ipc = 4
                    elseif (typ=='EJT8') then
                       ipc = 5
                    end if

                    ipg = ipc

                    !----- Calcul des contraintes au centre de gravite de l'element
                    call elem_B(vn,vb,detj,ksig(ipg,:),ie)

                    allocate(vsi0(size(vb,1)))  ;  allocate(vdsi(size(vb,1)))
                    vsi0 = matmul(vh,matmul(vb,vdl0))
                    vdsi = matmul(vh,matmul(vb,vdle))

                    !----- Calcul du deplacement relatif normal au centre de gravite de l'element	
                    deprel0 = matmul(vb,vdl0)
                    deprel = matmul(vb,vdle)

                    deallocate(vdl0,vdlI,vdlII,vdle,vb,vn)

                    !----- Cas de la loi 100 !!!
                    if (iloi==1) then

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
                      tau   = vsi0(2)+vdsi(2)

                      vcrit1 = sigma - RT
                      vcrit2 = abs(tau) + sigma*tan(phi) - C
                      vcrit3 = abs(tau) + sigma*tan(psi-pi/2.)             

                      if (vcrit3 < 0.d0) then
                        !- Critere en traction pure
                        MOT = 'ouvert'
                        if (vcrit1 > 0.d0) then
                           alpg(ipg) = (RT-vsi0(1))/vdsi(1)

                        elseif (vcrit2 >0.d0) then
                           signe = tau/abs(tau)
                           alpc = C - signe*vsi0(2) - vsi0(1)*tan(phi)
                           alpc = alpc/(signe*vdsi(2) + vdsi(1)*tan(phi))
                           alpg(ipg) = alpc
                        end if
                                  
                      else

                        !- Critere en cisaillement
                        MOT = 'cisaille'
                        if (vcrit1 > 0.d0) then
                           alpg(ipg) = (RT-vsi0(1))/vdsi(1)

                        elseif (vcrit2 > 0.d0) then
                           signe = tau/abs(tau)
                           alpc = C - signe*vsi0(2) - vsi0(1)*tan(phi)
                           alpc = alpc/(signe*vdsi(2) + vdsi(1)*tan(phi))
                           alpg(ipg) = alpc
                        end if
                      end if

                    !----- Cas de la loi 101 !!!
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
                         tau0 = sqrt(vsi0(2)**2 + vsi0(3)**2)
                         dtau = sqrt(vdsi(2)**2 + vdsi(3)**2)
                         
                         tau  = tau0 + dtau

                         if ((abs(vsi0(2)+vdsi(2))>=abs(vsi0(2))).and.(abs(vsi0(3)+vdsi(3))>=abs(vsi0(3))))then
                           tau  = tau0 + dtau
                         elseif ((abs(vsi0(2)+vdsi(2))<abs(vsi0(2))).and.(abs(vsi0(3)+vdsi(3))<abs(vsi0(3))))then
                           tau  = tau0 - dtau
                         end if
                      end if

                      vcrit1 = sigma - RT
                      vcrit2 = abs(tau) - C
                      vcrit3 = abs(tau) + sigma*tan(psi-pi/2.)

                      if (vcrit3 < 0.d0) then
                        !- Critere en traction pure

                        if (vcrit1 > 0.d0) then
                          alpg(ipg) = (RT-vsi0(1))/vdsi(1)
                          MOT = 'ouvert'
                        elseif (vcrit2 > 0.d0) then
                           if (dime==2) then
                              signe = tau/abs(tau)
                              alpc = C - signe*vsi0(2)
                              alpc = alpc/(signe*vdsi(2))
                              alpg(ipg) = alpc

                           elseif (dime==3) then   
                              alpg(ipg) = (C - tau0)/dtau
                           end if                                 
                           MOT = 'cisaille'
                         end if
                                    
                      else
                        !- Critere en cisaillement
                        MOT = 'cisaille'
                        if (vcrit1 > 0.d0) then
                          alpg(ipg) = (RT-vsi0(1))/vdsi(1)

                        elseif (vcrit2 > 0.d0) then
                          if (dime==2) then
                             signe = tau/abs(tau)
                             alpc = C - signe*vsi0(2)
                             alpc = alpc/(signe*vdsi(2))
                             alpg(ipg) = alpc

                          elseif (dime==3) then                                 
                                alpg(ipg) = (C - tau0)/dtau
                             end if     

                         end if
                      end if

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
                             end if       
                          end if
                      end if
                                
                      if (depreln >= depcrin) then
                          MOT2 = 'ouvert'
                          alpgo = (depcrin - deprel0(1))/deprel(1)
                      end if
                              
                      alpg(ipg) = minval((/alpgo, alpgc/))
                      if (alpg(ipg)==alpgo) MOT = MOT2
                      if (alpg(ipg)==alpgc) MOT = MOT1
                    else
                      stop "FIDES_interf_pilot : loi non encore programmee"
                    end if

                    deallocate(vsi0,vdsi,vh)     
    
                    rupt(ie) = MOT
                    alph(ie) = minval(alpg)

                    deallocate(alpg,ksig)
                end if 
            end do   
            ! Fin de la boucle d'element

            alphai = minval(alph) 
               
            if (alphai < 1.d0) then
               iedngi = find_num(alph,alphai)          
               MOTi = rupt(iedngi)        
               elemfiss = count(alph<=0)
            else
               alphai = 1.d0
            end if
       end if

    case(2)
        !----- Pilotage sur l'increment de contrainte avec recalcul de "dlam"
        
    case default
		stop "FIDES_interf_pilot : loi non encore programmee"

    end select

end subroutine interf_pilot

!********************************************************!

subroutine interf_change_etat()

!********************************************************!
!             Gestion de l'interpenetration et           !
!           de l'element le plus dangereux               !
!********************************************************!

    use variables, only : nelt, ktypel, nomtype, infele, vsol, ietatpg, &
                 &  histetatpg1, histetatpg2, interf, detoscill, comptoscill, dime, kloce
    use initialisation, only : init_vec, init_mat
    use lib_elem, only : elem_B, elem_kloce2
    implicit none

    ! Quantites globales
    real*8, dimension(:), allocatable :: vdle
    real*8, dimension(:,:), allocatable :: vn, vb, ksig

    ! Quantites principales : Deplacement, ouverture, etat
    real*8 :: depreln, depreln1, depreln2, sign1, sign2, detj, deprel(dime)
    integer :: ie, ipg, npg, oscill, ietat0, ietat1, ietat2, ndle
    character(len=5) :: typ

!********************************************************!

if (interf==1) then
    do ie = 1, nelt
    typ = nomtype(ktypel(ie))
    !----- Pour les elements d'interface seulement
	if ((typ=='EJQ4') .or. (typ=='EJQ6') .or. (typ=='EJT6') .or. (typ=='EJT8')) then
        
           !-----  Recuperation des informations sur les elements
           allocate(ksig(size(infele(ktypel(ie))%Q,1), size(infele(ktypel(ie))%Q,2)))
           ksig = infele(ktypel(ie))%Q
           npg = size(ksig,1)
       
           do ipg = 1, npg
            
           if (ietatpg(ie,ipg) /= 0) then
            
              if(detoscill(ie)==6) then
                  !print*,'Oscillation detectee sur l''element',ie
                  comptoscill=comptoscill+1 
                  oscill=1
                  ietatpg(ie,ipg)=-2 
                  detoscill(ie)=detoscill(ie)+1
               end if
               
               !----- Vecteur kloce de localisation pour assemblage
               call elem_kloce2(ie, ndle)
               allocate(vdle(ndle))
   
               !----- Deplacement total
               vdle = vsol(kloce(1:ndle))
           
               !----- Calcul du deplacement relatif normal au centre de gravite de l'element 
 
               call elem_B(vn,vb,detj,ksig(ipg,:),ie)
               deprel = matmul(vb,vdle)
               depreln = deprel(1)
 
               if (ietatpg(ie,ipg) == -1) then

                  if (depreln <= 0.d0) ietatpg(ie,ipg) = 0  !la fissure s'est refermee apres ouverture 
                     
               elseif ((ietatpg(ie,ipg)==1) .or. (ietatpg(ie,ipg)==2)) then                     
                  ietat2 = histetatpg2(ie,ipg)
                  ietat1 = histetatpg1(ie,ipg)
                  ietat0 = ietatpg(ie,ipg)

                  oscill = 0
                  
                  if ((ietat0==1).and.(ietat1==2)) detoscill(ie)=detoscill(ie)+1
                  if ((ietat0==2).and.(ietat1==1)) detoscill(ie)=detoscill(ie)+1
                  
                  if ((ietat0==ietat2) .and. (ietat0/=ietat1)) then
                     oscill=1
                     ietatpg(ie,ipg)=-2
                     comptoscill=comptoscill+1 
                  end if
                  
                  histetatpg2(ie,ipg) = histetatpg1(ie,ipg)
                  histetatpg1(ie,ipg) = ietat0
                   
                  if ((ietatpg(ie,ipg)==1) .and. (depreln < 0.d0) .and. (oscill == 0)) then
                 	ietatpg(ie,ipg) = 2
                 	detoscill(ie)=detoscill(ie)+1  
                  !elseif ((ietatpg(ie,ipg)==2) .and. (depreln > 0.d0) .and. (oscill == 0)) then
                  !	ietatpg(ie,ipg) = 1
                  !	detoscill(ie)=detoscill(ie)+1
                  end if                  
               end if

               deallocate(vdle,vb,vn)

            else             
                ! IL FAUT PREVOIR ICI DE TRAITER LA DETECTION DE L'ELEMENT 
                ! LE PLUS DANGEREUX (I.E. CELUI QUI FISSURE EN PREMIER)
            end if
          end do

         deallocate(ksig)
    
    end if
    end do
end if

end subroutine interf_change_etat

!********************************************************!

subroutine interf_rupture(ie,ipg,vsig,vnle)

!********************************************************!
!       Fonction de detection de la rupture              !
!               des elements d'interface                 !
!********************************************************!

    use variables, only : ktypel, nomtype, ietatpg
    implicit none

    real*8, dimension(:,:), intent(inout) :: vsig, vnle
    integer, intent(in) :: ie, ipg
    character(len=5) :: typ

!********************************************************!

    typ = nomtype(ktypel(ie))

    if ((typ=='EJQ4').or.(typ=='EJQ6').or.(typ=='EJT6').or.(typ=='EJT8')) then
        if (ietatpg(ie,ipg) == 1) then
            vnle(:,ipg) = vnle(:,ipg)/1.0d40
            vsig(:,ipg) = vsig(:,ipg)/1.0d40

        elseif (ietatpg(ie,ipg) == 2) then
            vnle(:,ipg) = vnle(2,ipg)/1.0d40
        end if
    end if

end subroutine interf_rupture
!********************************************************!

subroutine interf_modul(ie,ipg,vh,vprel)

!********************************************************!
!                Gestion de la raideur                   !
!********************************************************!

    use variables, only : dime, ktypel, nomtype, ietatpg
    implicit none

    real*8, dimension(:), intent(in) :: vprel    
    real*8, dimension(:,:), intent(inout) :: vh
    integer, intent(in) :: ie, ipg    
    integer :: iloi
    character(len=5) :: typ    
    
!********************************************************!

    iloi=vprel(1)
    typ = nomtype(ktypel(ie))

    if (((typ=='EJQ4').or.(typ=='EJQ6').or.(typ=='EJT6').or.(typ=='EJT8')) &
        .and. ((iloi==100) .or. (iloi ==101) .or. (iloi ==102))) then

            !***** Modification de la matrice pour element casse
            if ((ietatpg(ie,ipg) == 1) .or. (ietatpg(ie,ipg) == -1) .or. (ietatpg(ie,ipg) == -2)) then
               vh = vh/1.0d50

            !***** Modification de la matrice pour element cisaille
            elseif (ietatpg(ie,ipg)==2) then
               if (dime == 2) then
                   vh(2,2) = vh(2,2)/1.0d50

               elseif (dime == 3) then
                   vh(2,2) = vh(2,2)/1.0d50
                   vh(3,3) = vh(3,3)/1.0d50
               end if

       !***** Elastoplastique
            elseif (ietatpg(ie,ipg)==3) then   

            end if
     end if

end subroutine interf_modul

!********************************************************!

subroutine interf_loi(iloi,icomp,ie,ipg)

!========================================================!
!     Calcul des contraintes pour element d'interface    !
!========================================================!

    use variables, only : ktypel, nomtype, ietatpg
    implicit none

    character(len=5) :: typ
    integer, intent(in) :: icomp, ie, ipg
    integer, intent(out) :: iloi

!********************************************************!
    
    typ = nomtype(ktypel(ie))

    if ((typ=='EJQ4').or.(typ=='EJQ6').or.(typ=='EJT6').or.(typ=='EJT8')) then
        if ((icomp /= 100) .and. (icomp /= 101) .and. (icomp /= 102)) then
            stop 'Discordance de la loi de l"element d"interface : Veuillez verifier !!! '
        else                        
            if (iloi==102) then               
              ! if (ietatpg(ie,ipg)==1 .or. ietatpg(ie,ipg)==2) iloi = 101                        
            end if
        end if
    else
       iloi = icomp
    end if

end subroutine interf_loi

!********************************************************!

subroutine interf_init()

!********************************************************!
!  Initialisations liees a la gestion des etats de       !
!               elements d'interface                     !
!********************************************************!
    
    use variables, only : nelt, nnt, nomtype, infele, interf, ietatpg, &
           & histetatpg1, histetatpg2, iedngi, ieendo, endo, interf_macro, mode
    use initialisation, only : init_mat, init_vec
    implicit none

    integer, dimension(size(nomtype,1)) :: npg
    integer :: itype, ntype, npgm
    character(len=5) :: typ

!********************************************************!
    
    interf = 0
    ntype = size(nomtype,1)
    npg = 0

    do itype = 1, ntype
       typ = nomtype(itype)
       if ((typ=='EJQ4').or.(typ=='EJQ6').or.(typ=='EJT6').or.(typ=='EJT8')) then

         npg(itype) = size(infele(itype)%W,1)
         interf = 1

       end if
    end do

    if (interf==1 .and. interf_macro==0)then
        ieendo = 0
        iedngi = 0
        npgm = maxval(npg)
        call init_mat(ietatpg,nelt,npgm)
        call init_mat(endo,nelt,npgm)
        call init_mat(histetatpg1,nelt,npgm)
        call init_mat(histetatpg2,nelt,npgm)
    end if

end subroutine interf_init

!********************************************************!

subroutine interf_stock()

!********************************************************!
!   Stockage des resultats à chaque pas de temps pour    !
!   le  modele  de  fissuration  avec  gestion  des      !
!   transitions etat(-1)--> etat(1), etat(-2)--> etat(2) !
!********************************************************!
    
    use variables, only :  interf, ietatpg, interf_macro
    implicit none

! *******************************************************!

	if ((interf==1).or.(interf_macro==1)) then
	    where (ietatpg==-1 .or. ietatpg==-2) ietatpg = 1
	end if

end subroutine interf_stock

!--------------------------------------------------------!

subroutine interf_distal()

!---------------------------------------------------------------!
! Fonction de definition des proprietes mecaniques aleatoires   !
!   Donnees en entree :                                         !
!   - kprop : numerotation des groupes des elements             !
!   - vprelg : proprietes mecaniques (fonction du modele)       !
!                                                               !
!   Valeurs en sortie :                                         !
!   - kprop : actualise en fonction du numero de l'element      !
!             (1 element probabiliste = 1 groupe)               !
!   - vprelg : modifie pour prendre en compte la propriete      !
!      aleatoire consideree - Attention depend du modele choisi !
!---------------------------------------------------------------!
    use variables, only : nelt, dime, kconec, ktypel, nomtype, infele, &
            & vprelg, kprop, kprop0, interf, interf_macro, fiss, idmax, pi, &
            & Dg, fc, &
            & vitrav, &
            !--- gestion de l'alea (deb)
            & alea, young, resist, interfa
            !--- gestion de l'alea (fin)
    use lib_elem, only : elem_B
    use initialisation
    use proprietes_aleatoires    
    use math, only : find_vec

    implicit none
    
    real*8, dimension(:,:), allocatable :: vprelg0, vn, vb, ksig
    real*8, dimension(:), allocatable :: vpg, RR, EE, CCb,CCi
    real*8 :: Rt(1), C(1), E(1), Vel(nelt), loi1(5), loi2(5), loi3(5)
    real*8 :: detj, detj1, detj2, epais, poids, Ve, Vg, &
            &  r, Vem, A, B, CD, b0, c0

    integer, dimension(:), allocatable :: ielm, ielf
    integer, dimension(:,:), allocatable :: cote
    integer :: igbep, iebe, igfip, iefi, igabp, ieab
    integer :: i, ie, ie0, nbg, ncote, np, ipg, npg, id, nelf, ele, &
            & iloi, icote, dPlay, cal            
    integer :: igfp, igma, igin , j, k , l, m, n

    character(len=5) :: typea, typeb
    logical :: iprobVO, iprobFI, iprobAB, iprobMY, resp, modp, intp

!---------------------------------------------------------------!

!- Dans le cas d'un calcul avec elements d'interface uniquement-!
if (fiss==0) then
  if ((interf==1 .or. interf_macro==1) .and. (alea)) then

    !----- preparation de la renumerotation de kprop et vprelg  !
    ! Principe : On conserve la numerotation pour tous les      !
    ! groupes dont l'ancien numero est inferieur a nga. Puis on !
    ! "decale" d'un cran dans la numerotation tous les groupes  !
    ! qui ont un ancien numero plus grand que nga. Enfin, tous  !
    ! les elements du groupe nga sont renumerote en nga+i.      !
    
    !-----  Elements massifs probabilistes pour beton sain  ----!
    iebe  = 0
    
    !-- Elements d'interface probabilistes pour beton fissure --!
    iefi  = 0
    
    !---- Elements d'interface probabilistes pour acier beton --!
    ieab  = 0
    
    !---- Elements probabilistes totaux --------------------- --! 
    ie0 = 0
        
    nbg = maxval(kprop)    
    call init_mat(vprelg0,(size(kprop,1)+nbg),size(vprelg,2))
    vprelg0(1:nbg,:) = vprelg(1:nbg,:)

    Vel = 0.d0
    call init_vec(RR,nelt)
    call init_vec(EE,nelt)
    call init_vec(CCb,nelt)
    call init_vec(CCi,nelt)

    !----- Indice pour trace de distributions de proprietes aleatoires
    iprobVO = .false.
    iprobFI = .false.
    iprobAB = .false.
    iprobMY = .false.

    !-- Traitement des elements d'interface probabilistes d'abord

    do ie = 1, nelt
        
       !-- detection d'un element probabiliste (module et/ou resistance et/ou energie)
       resp = .false.
       if (allocated(resist)) then
           if (count(resist%num==kprop0(ie))==1) resp=.true.
       end if

       typea = nomtype(ktypel(ie))

       !-- si l'element appartient a un groupe probabiliste (module et/ou resistance et/ou energie)
       if ((resp) .and. (typea=='EJQ4'.or.typea=='EJQ6'.or.typea=='EJT6'.or.typea=='EJT8')) then
        
               allocate(vitrav(1))
               vitrav= find_vec(resist%num,kprop0(ie))
               igfp = vitrav(1)
               deallocate(vitrav)
               
               !------------------ deb fissure probabiliste -----------------
               ie0 = ie0 + 1
               iefi = iefi + 1

               !----- Recherche des elements massifs voisins
               if (typea == 'EJQ4') then                   
                   allocate(cote(2,2))
               elseif (typea == 'EJQ6' .or. typea == 'EJT6') then
                   allocate(cote(2,3))
               elseif (typea == 'EJT8') then
                   allocate(cote(2,4))
               end if

               cote = infele(ktypel(ie))%face
               ncote = size(cote,1)
               np = size(cote,2)   

               if (ncote/=2) stop 'FIDES_interf_distal : erreur de programmation. Il faut 2 cotes'
              
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
                           end if
             
                           if (np==4) then
                             if (kconec(i,j)==kconec(ie,cote(icote,3))) m = 1
                             if (kconec(i,j)==kconec(ie,cote(icote,4))) n = 1
                           end if             
                       end do
         
                       if ((k==1.and.l==1.and.np==2).or.(k==1.and.l==1.and.m==1.and.np==3) &
                          & .or. (k==1.and.l==1.and.m==1.and.n==1.and.np==4)) then
                         ielm(icote) = i
                         exit
                       end if
                     end if
                   end do

                   if (icote > 1) then
                       if (ktypel(ielm(1)) /= ktypel(ielm(icote))) stop 'FIDES_interf_distal : elements massifs incompatibles'
                   end if
               end do

               !----- Recuperation des poids et points d'integration de l'element               
               allocate(vpg(size(infele(ktypel(ielm(1)))%W)))
               allocate(ksig(size(infele(ktypel(ielm(1)))%Q,1),size(infele(ktypel(ielm(1)))%Q,2)))			   
 			   
               vpg = infele(ktypel(ielm(1)))%W
               ksig = infele(ktypel(ielm(1)))%Q
               npg = size(ksig,1)
               
               !----- Recuperation de l'epaisseur si 2D
               epais = 1.
               iloi = vprelg(kprop(ielm(1)),1)
               
               if (iloi == 1) then                 
                   if (dime==2) epais = vprelg(kprop0(ielm(1)),idmax-2)   ! (position de l'epaisseur)
               else
                   print*,'Non encore programme pour la loi ', iloi                   
                   stop 'FIDES_interf_distal : valable pour l''elasticite seulement'
               end if
               
               !----- Calcul du volume fissure (Ve(ie1)+Ve(ie2))
               Ve = 0.d0
               do ipg = 1, npg
                   poids = vpg(ipg)
                   !----- Calcul des fonctions d'interpolation et des derivees
                   
                   call elem_B(vn,vb,detj,ksig(ipg,:),ielm(1))
                   detj1 = detj
                   deallocate(vb,vn)
                   
                   call elem_B(vn,vb,detj,ksig(ipg,:),ielm(2))
                   detj2 = detj                   
                   deallocate(vb,vn) 
   
                   !----- integration matrice vke
                   Ve = Ve + (detj1+detj2)*poids*epais                   
               end do
    
               Vel(ie) = Ve
               Vg = (4./3.)*pi*(Dg/2.)**3.
               r = Ve/Vg                

               deallocate(ksig,vpg,ielm,cote)
                              
               !----- Resistance a la traction  
               loi1(1) = 0.    ! b
               loi1(2) = 0.    ! c
               loi1(3) = 0.    ! moy
               loi1(4) = 0.    ! ect
               loi1(5) = resist(igfp)%loi
               
               if (loi1(5)==1) then
                   ! Loi normale
                   if (resist(igfp)%ipa/=0) then
                       loi1(3) = resist(igfp)%param(1)
                       loi1(4) = resist(igfp)%param(2)
                   else
                       call rossi_coef(loi1(3),loi1(4),r,fc,'RESI')
                   end if
               elseif (loi1(5)==2) then
                   ! Loi de Weibull
                   if (resist(igfp)%ipa/=0) then
                       loi1(1) = resist(igfp)%param(1)
                       loi1(2) = resist(igfp)%param(2)
                   else
                       !call interf_wblcoef(loi1(1),loi1(2),r,fc)
                       loi1(1) = 4.
                       loi1(2) = 1.
                   end if
               elseif (loi1(5)==3) then
                   ! Loi log-normale
                   if (resist(igfp)%ipa/=0) then
                       loi1(1) = resist(igfp)%param(1)
                       loi1(2) = resist(igfp)%param(2)
                   else
                       call rossi_coef(loi1(3),loi1(4),r,fc,'RESI')
                   end if
               end if
   
               !----- Calcul de la resistance et de la cohesion              
               Rt = distr_alea(loi1,1,ie0)
               C = 5.*Rt
               
               !----- Preparation pour trace graphique eventuel
               RR(iefi) = Rt(1)
               CCb(iefi) = C(1)
               iprobFI = .true.
               iprobVO = .true.
               
               !----- Stockage resistance a la traction et cohesion
               vprelg0(nbg + ie0,:) = vprelg(kprop(ie),:)
               kprop(ie) = nbg + ie0
               iloi = vprelg0(kprop(ie),1)

               if (iloi == 100 .or. iloi == 101 .or. iloi == 102) then
                   id = 6
                   if (dime == 3) id = id-1
                   vprelg0(kprop(ie),id) = Rt(1)
                   vprelg0(kprop(ie),id+1) = C(1)
                   
                   !--- petite modification
                   if (iloi == 101) vprelg0(kprop(ie),id+3) = atan(Rt(1)/(C(1)+1.d-20))
               else
                   print*,'FIDES_interf_distal : loi ',iloi ,' non probabilisee'
                   stop
               end if
               !------------------ fin fissure probabiliste -----------------
       end if
       
       !-- detection d'un element d'interface acier/beton probabiliste

       intp = .false.
       if (allocated(interfa)) then
           if (count(interfa%num==kprop0(ie))==1) intp=.true.
       end if

       if ((intp) .and. (typea=='EJQ4'.or.typea=='EJQ6'.or.typea=='EJT6'.or.typea=='EJT8')) then
       
               allocate(vitrav(1))
               vitrav= find_vec(interfa%num,kprop0(ie))
               igin = vitrav(1)
               deallocate(vitrav)

               !------------------ deb ac/beton probabiliste ----------------
               ! cal = 1 : distribution aleatoire sur la cohesion de l'interface acier beton
               cal = 0
               if (cal == 1) then
                   ie0 = ie0 + 1
                   ieab = ieab + 1        
        
                   loi2(1) = 0.    ! b
                   loi2(2) = 0.    ! c
                   loi2(3) = 0.    ! moy
                   loi2(4) = 0.    ! ect
                   loi2(5) = interfa(igin)%loi

                   if (loi2(5)==1) then
                       ! Loi normale
                       if (interfa(igin)%ipa/=0) then
                           loi2(3) = interfa(igin)%param(1)
                           loi2(4) = interfa(igin)%param(2)
                       else
                           loi2(3) = 25.0
                           loi2(4) = 12.5
                       end if
                   end if

                   !----- Calcul de la resistance et de la cohesion
                   C = distr_alea(loi2,1,ie0)

                   !----- Preparation pour trace graphique eventuel
                   CCi(ieab) = C(1)
                   iprobAB = .true.

                   vprelg0(nbg + ie0,:) = vprelg(kprop(ie),:)
                   kprop(ie) = nbg + ie0
                   iloi = vprelg0(kprop(ie),1)

                   if (iloi == 100 .or. iloi == 101 .or. iloi == 102) then
                       id = 6
                       if (dime == 3) id = id-1
                       vprelg0(kprop(ie),id+1) = C(1)
                   else
                       print*,'MEF_interf_distal : loi ',iloi ,' non probabilisee'
                       stop
                   end if              
               end if
       end if
   end do

   !----- Traitement des elements massifs probabilistes pour beton ensuite
   do ie = 1, nelt
       
       !-- detection d'un element massif probabiliste (module d'young)

       modp = .false.
       if (allocated(young)) then
           if (count(young%num==kprop0(ie))==1) modp=.true.
       end if

       typeb = nomtype(ktypel(ie))
 
       if ((modp).and.(iefi/=0).and.(typeb=='MBT3'.or.typeb=='MBT6'.or.typeb=='MBQ4'.or. &   ! element 2D
                                   & typeb=='MTT4'.or.typeb=='MTP6'.or.typeb=='MTH8')) then  !  element 3D

               allocate(vitrav(1))
               vitrav= find_vec(young%num,kprop0(ie))
               igma = vitrav(1)
               deallocate(vitrav)

               !------------------ deb fissure probabiliste -----------------
               ie0 = ie0 + 1
               iebe = iebe + 1
                
               !----- Recherche des elements d'interface voisins
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
                           end if
             
                           if (np==4) then
                             if (kconec(i,j)==kconec(ie,cote(icote,3))) m = 1
                             if (kconec(i,j)==kconec(ie,cote(icote,4))) n = 1
                           end if             
                       end do
         
                       if ((k==1.and.l==1.and.np==2).or.(k==1.and.l==1.and.m==1.and.np==3) &
                           &.or.(k==1.and.l==1.and.m==1.and.n==1.and.np==4)) then                          
                          
                           if (count(resist%num==kprop0(i))==1) ielf(icote) = i
                           exit
                         
                       end if
                     end if
                   end do
               end do
                       
               !----- Calcul du volume moyen massif : Vem = 1/n(Ve1 +...+ Ven)
               Vem = 0.
               nelf = count(ielf /= 0)

               do icote = 1, ncote
                   if (ielf(icote) /= 0) Vem = Vem + Vel(ielf(icote))
               end do

               if (nelf /= 0) Vem = Vem/nelf
               Vg = (4./3.)*pi*(Dg/2.)**3
               r = Vem/Vg        

               deallocate(ielf,cote)
               
               !----- Module d'Young (loi normale ou log-normale)

               loi3(1) = 0.    ! b
               loi3(2) = 0.    ! c
               loi3(3) = 0.    ! moy
               loi3(4) = 0.    ! ect   
               loi3(5) = young(igma)%loi  ! num

               id=4
               if (dime==3) id=id-1
               
               if (loi3(5) == 1) then
                   ! Loi normale
                   if (young(igma)%ipa/=0) then
                       loi3(3) = young(igma)%param(1)
                       loi3(4) = young(igma)%param(2)
                   else
                       call rossi_coef(loi3(3),loi3(4),r,fc,'MODU',vprelg(kprop(ie),id))
                   end if
                   if (nelf==0) loi3(4)=0.d0
               elseif (loi3(5) == 2) then
                   ! Loi de Weibull
                   stop 'FIDES_interf_distal : loi probabiliste non programmee pour les modules'
               elseif (loi3(5) == 3) then
                   ! Loi log-normale
                   if (young(igma)%ipa/=0) then
                       loi3(3) = young(igma)%param(1)
                       loi3(4) = young(igma)%param(2)
                   else
                       call rossi_coef(loi3(3),loi3(4),r,fc,'MODU',vprelg(kprop(ie),id))
                   end if
                   
                   if (nelf==0) loi3(4)=0.d0
               end if
                
               !----- Calcul du module elastique
               E = distr_alea(loi3,1,ie0)
    
               !----- Preparation pour trace graphique eventuel
               EE(iebe) = E(1)
               iprobMY = .true.

               !----- Stockage modules elastiques
               vprelg0(nbg + ie0,:) = vprelg(kprop(ie),:)
               kprop(ie) = nbg + ie0

               id = 4
               if (dime == 3) id = id - 1
               vprelg0(kprop(ie),id) = E(1)
                
               !------------------ fin fissure probabiliste -----------------
       end if
   end do

   deallocate(vprelg)
   call init_mat(vprelg,ie0+nbg,size(vprelg,2))

   ! Sorties :
   vprelg(1:ie0+nbg,:) = vprelg0(1:ie0+nbg,:)

   if (iprobVO) then
        dPlay = 1
        if (dplay>=1) then
            print*,' '
            print*,'Maillage beton probabiliste :'
            print*,'-----------------------------'
        end if
        ! Visualisation de la distribution des volumes
        call distrib_volumes(Vel,Dg,50,dPlay)
   end if
   
   if (iprobFI) then
        dPlay = 1
        if (dplay>=1) then
            print*,'Fissuration probabiliste :'
            print*,'--------------------------'
            print'(a10,e12.5,a16,e12.5)','- Dmax :',Dg,' - Resistance :',fc
            print'(a10,e12.5)','- Vg :',Vg
        end if
        ! Visualisation de la distribution de resistance
        call alea_trac(RR,loi1,200,dPlay,'resistance')
        
        ! Visualisation de la distribution de cohesion
        call alea_trac(CCb,loi1,200,dPlay,'cohesion beton')
   end if
   
   if (iprobMY) then  
        dPlay = 1
        ! Visualisation de la distribution de module elastique
        call alea_trac(EE,loi3,200,dPlay,'module')
   end if
    
   if (iprobAB) then  
        dPlay = 1
        if (dplay>=1) then
            print*,'Interface probabiliste :'
            print*,'------------------------'
        end if
        ! Visualisation de la distribution de cohesion
        call alea_trac(CCi,loi2,200,dPlay,'cohesion - interf. A-B')
   end if

   !call ecriture_distribution(RR,iefi,EE,iebe)   ! utilitaire.f

   !----- Liberation memoire JLT
   deallocate(vprelg0,EE,RR,CCb,CCi)
       
 end if
end if
end subroutine interf_distal

!-------------------------------------------------------!

subroutine interf_resichute(ipas)

!---------------------------------------------------------------!
! Fonction pemettant le calcul de la chute de resistance pour   !
! la modélisation de la propagation de fissure sous charge      !
!   Donnees en entree :                                         !
!   - ipas    : numero du pas de temps                          !
!   Valeurs en sortie :                                         !
!   - vprelg  : proprietes mecaniques au temps (n+1) apres      !
!               prise en compte de la chute de resistance       !
!               (prise en compte du fluage)                     !
!---------------------------------------------------------------!
    
    use variables, only : interf, nelt, nomtype, ktypel, kprop0, kprop, &
                        & infele, vcont, kcont, dtps, dime, vprelg, vprelg0, &
                        !--- prise en compte de l'alea
                        & resist, &
                        !--- prise en compte de l'alea
                        & tps, vprelg1
    use lib_elem, only : elem_B
    use initialisation
    use proprietes_aleatoires    
    use math

    implicit none

    integer :: ipas

    real*8, dimension(:), allocatable   :: vpg
    real*8, dimension(:,:), allocatable :: ksig, vsig, vn, vb

    integer :: ie, nc, npg, idim, iloc, ipg, iloi, id
    real*8  :: poids, detj, sigm, ptot, t, Dt, dRt, Rt0, Rt1, Rt, dRc, Rc0, Rc1, Rc
    real*8  :: xa,xb,xc,xi,ft
    character(len=5) :: typ

!---------------------------------------------------------------!

!- Dans le cas d'un calcul avec elements d'interface uniquement-!
if (interf == 1) then

    !-- Traitement des elements d'interface béton

    do ie = 1, nelt
       typ = nomtype(ktypel(ie))

       if (allocated(resist)) then
           if ((count(resist%num==kprop0(ie))==1).and.(typ=='EJQ4'.or.typ=='EJQ6'.or.typ=='EJT6'.or.typ=='EJT8')) then

                !----- Recuperation des informations elementaires -----
                call init_vec(vpg, size(infele(ktypel(ie))%W))
                call init_mat(ksig, size(infele(ktypel(ie))%Q,1), size(infele(ktypel(ie))%Q,2))
                nc   = infele(ktypel(ie))%ncgr
                vpg  = infele(ktypel(ie))%W
                ksig = infele(ktypel(ie))%Q
                npg  = size(ksig,1)

                !----- localisation des contraintes elementaires dans vcont -----
                idim = npg*nc;
                iloc = kcont(ie)

                !----- Recuperation des contraintes elementaires dans vcont -----
                call init_mat(vsig, nc, npg)
                vsig  = reshape(vcont(iloc : (iloc+ idim - 1)),(/ nc, npg/))     

                !----- Calcul de la contrainte moyenne sur l'élément (contrainte au centre de gravite) -----
                sigm = 0.d0
                ptot = 0.d0
                do ipg = 1, npg
                    poids = vpg(ipg)
                    !----- Calcul des fonctions d'interpolation et des derivees
                    call elem_B(vn,vb,detj,ksig(ipg,:),ie)
 
                    !----- calcul de la contrainte moyenne
                    sigm = sigm + detj*poids*vsig(1,ipg)
                    ptot = ptot + detj*poids
                    deallocate(vb,vn)
                end do
                deallocate(ksig,vpg,vsig)
                sigm = sigm / ptot

                !----- Recuperation des resistances au temps t(n)

                iloi = vprelg(kprop(ie),1)
                if (iloi == 100 .or. iloi == 101 .or. iloi == 102) then
                    id = 6
                    if (dime == 3) id = id-1
                    Rt0 = vprelg0(kprop(ie),id)
                    Rc0 = vprelg0(kprop(ie),id+1)
                    Rt1 = vprelg1(kprop(ie),id)
                    Rc1 = vprelg1(kprop(ie),id+1)
                else
                    print*,'FIDES_interf_resichu : la loi ',iloi ,' n''est pas une loi beton !'
                    stop
                end if               

                !----- Calcul de la vitesse de chute de contrainte

                !xa = log(10.d0)/(-0.016D0)
                !xb = -1.d0 * xa * 0.976D0
                xa = log(10.d0)/(-0.05D0)
                xb = -1.d0 * xa * 0.95D0
                if (sigm>=0) then
                    xi = sigm / RT0
                    if (xi >1.D0) xi = 1.d0
                    dRt = RT0*(xi-1)/exp(xa*xi+xb)
                else
                    dRt = 0.d0
                end if
                dRc = 5 * dRt
                
                t  = tps(ipas)
                Dt = dtps(ipas)

                !xc  = exp(xa*xi+xb)
                !ft  = 1.d0 - exp(-(2.d0*Dt/xc)**2)
                !ft  = 1.d0/(1+2000.d0*exp(-Dt/200))
                ft = 1.d0

                !----- Calcul de la resistance au temps t(n+1)
                Rt = Rt1 + Dt * ft * dRt
                Rc = Rc1 + Dt * ft * dRc
                if (Rt<=0.d0) Rt=0.d0
                if (Rc<=0.d0) Rc=0.d0

                iloi = vprelg(kprop(ie),1)
                if (iloi == 100 .or. iloi == 101 .or. iloi == 102) then
                    id = 6
                    if (dime == 3) id = id-1
                    vprelg(kprop(ie),id)   = Rt
                    vprelg(kprop(ie),id+1) = Rc
                else
                    print*,'FIDES_interf_resichu : la loi ',iloi ,' n''est pas une loi beton !'
                    stop
                end if
           end if
       end if              
    end do
end if

end subroutine interf_resichute

!-------------------------------------------------------!

end module element_interface

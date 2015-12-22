module aiguillage

!********************************************************************************!
!	 Gestion de l'aiguillage element d'interface et calcul fissure           ! 
!********************************************************************************!

contains

!*****************************************************************!

subroutine pilot(imetpilo,vsol,vduI,vduII,dlam)

!-----------------------------------------------------------------!
! Fonction d'aiguillage pour le pilotage du calcul                !
!-----------------------------------------------------------------!

     use variables, only : interf, fiss, interf_macro, iedngi, iedngm, iedngma, iedng, & 
                          & mrtrav, dime, nelt, inorm, irloc
     use element_interface, only : interf_pilot
     use fissuration, only : fiss_pilot
     use interface_macro, only : interf_macro_pilot
     use initialisation, only : init_mat

     real*8, dimension(:), intent(in) :: vduII, vsol
     real*8, dimension(:), intent(inout) :: vduI
     real*8, intent(inout) :: dlam
     integer, intent(in) :: imetpilo
     real*8 :: alpha, alphai, alpham, alphama
     integer :: elefissi, elefissm, elefissma
     character(len=8) :: MOTi, MOTm, MOTma, MOT 

     alpha   = 1.d0
     alphai  = 1.d0
     alpham  = 1.d0
     alphama = 1.d0
     MOTi  = '    '
     MOTm  = '    '
     MOTma = '    '
     MOT   = '    '
     iedng = 0

     elefissi = 0
     elefissm = 0
     elefissma = 0

     call init_mat(mrtrav,nelt,dime*dime)

     if (interf == 1) call interf_pilot(imetpilo,vsol,vduI,vduII,dlam,alphai,MOTi,elefissi)
     if (fiss == 1) call fiss_pilot(imetpilo,vsol,vduI,vduII,dlam,alpham,MOTm,elefissm)
     if (interf_macro == 1) call interf_macro_pilot(imetpilo,vsol,vduI,vduII,dlam,alphama,MOTma,elefissma)

     !----- Choix du minimum de pas de chargement entre les deux pilotages
     alpha = minval((/alphai,alpham,alphama/))
     
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

        ! Element d'interface macro
        if (alpha == alphama) then
            iedng = iedngma
            MOT = MOTma
        end if

        if (alpha <= 0.d0) alpha = 1.d-10
        if (alpha >= 1.d0) alpha = 1.d0

        print'(a11,e11.6,a11,e11.6,a13,i7,a3,a8,a21,i7,2x,a4)','Pilotage : ', &
             &      alpha,' dlam : ',dlam,'   element : ',iedng,'  ',MOT,' [ elts critiques : ',elefissi+elefissm+elefissma,' ]'
    end if

    !----- On reactualise le champ de deplacements (vduI et facteur de charge)
    vduI = alpha*vduI
    dlam = alpha*dlam

    deallocate(mrtrav)

end subroutine pilot

!-------------------------------------------------------!
subroutine pilot2(imetpilo,vsol,vduI,vduII,dlam)

    use variables, only : dime, nelt, idmax, ktypel, nomtype, infele, &
        & vprelg, kprop, kprop0, ietatpg, iedngi, ieendo, interf, pi, ipas, &
        & ietat, iedngm, fiss, mrtrav, vrtrav, inorm, irloc, iedng, kloce
    use utilitaire, only : princ
    use fissuration, only : fiss_changement_repere
    use initialisation, only : init_vec, init_mat
    use lib_elem, only : elem_B, elem_kloce2, elem_hooke
    use math
    implicit none

    ! Variables IN
    real*8, dimension(:), intent(in) :: vduII, vsol
    integer, intent(in) :: imetpilo

    ! Variables IN-OUT
    real*8, dimension(:), intent(inout) :: vduI
    real*8, intent(inout) :: dlam

    ! Quantites principales
    real*8, dimension(:), allocatable :: vdle, vdl0, vsi0, vdsi, &
                     & alpg, vdlI, vdlII, vsi, vsi0loc,vdsiloc
                     
    real*8, dimension(:,:), allocatable :: vn, vb, vh, ksig, P
    
    real*8 :: alpc, sign0, C, RT, phi, psi, sigma, tau, tau0, dtau, vdsig, vsig0, &
        & vcrit1, vcrit2, vcrit3, detj, signe, coef, alpgo, alpgc, s1
    real*8 :: depcrin, depcrit, depreln, deprelt,  deprelt1, deprelt2, depreln0, deprelt0
    real*8 :: vala, valb, valc, delta, alp1, alp2, alpha
    real*8 :: deprel(dime), deprel0(dime), vprel(idmax), alph(nelt), vbid(dime), vbid2(1,dime*dime), V(dime,dime)
    
    character(len=5) :: typ
    character(len=8) :: MOT, MOT1, MOT2
    character(len=8), dimension(nelt) :: rupt
    integer :: i, j, k, ie, ipg, npg, id, ipgd, iloi, ipc, ndle, elefiss
    logical :: ipremf, ivierge, iendo, irupt

    !---------------------------

    !- Initialisation
    alpha  = 1.d0
    iedng = 0
    call init_mat(mrtrav,nelt,dime*dime)
    alph = 1.d0

    !- Boucle sur les elements
    do ie = 1,nelt

        MOT = '  '
        MOT1 = '  '
        MOT2 = '  '
        typ = nomtype(ktypel(ie))

        !----- Proprietes elementaires
        vprel = vprelg(kprop(ie),1:idmax)

        !----- Recuperation des parametres de la loi
        iloi = int(vprel(1))
             
        !- Pour le cas des elements d'interface seulement
        if ((typ == 'EJQ4' .or. typ == 'EJQ6' .or. typ == 'EJT6' .or. typ == 'EJT8') &
           & .and. (iloi==100 .or. iloi==101 .or. iloi==102 .or. iloi==103)) then

             if (all(ietatpg(ie,:)==0) .or. all(ietatpg(ie,:)==3)) then
                    !-----  Recuperation des informations sur les elements            
                    call init_mat(ksig, size(infele(ktypel(ie))%Q,1), size(infele(ktypel(ie))%Q,2))
                    ksig = infele(ktypel(ie))%Q
                    npg = size(ksig,1)                   
  
                    allocate(alpg(npg))
                    alpg = 1.d0

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
                       
                    !---- Recuperation du point d'integration central en fonction du type d'element
                    if (typ=='EJQ4'.or.typ=='EJQ6') then
                       ipc = 1
                    elseif (typ=='EJT6') then
                       ipc = 4
                    elseif (typ=='EJT8') then
                       ipc = 5
                    end if

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
                    if (iloi==100) then
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
                        MOT = 'ouvert'
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

               end if   ! end if iloi
            end if   ! end if interf ==1

            !- Pour le cas beton fissurant
            if ((typ=='MBT3'.or.typ=='MBQ4'.or.typ=='MBT6'.or.typ=='MTT4'.or. &
                  & typ=='MTP6'.or.typ=='MTH8') .and. iloi==12) then

               if (ietat(ie)==0) then
                    !-----  Recuperation des informations sur les elements
                    call init_mat(ksig, size(infele(ktypel(ie))%Q,1), size(infele(ktypel(ie))%Q,2))
                    ksig = infele(ktypel(ie))%Q
                    npg = size(ksig,1)

                    !----- Initialisations relatives aux elements
                    call init_vec(vrtrav,npg*dime*dime)
                    allocate(alpg(npg))
                    alpg = 1.d0

                    do ipg = 1,1
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
                       vsi = vsi0+vdsi
                       call princ(vsi,vbid,V)

                       !----- Recuperation de la matrice de changement de repere global --> local principal
                       call init_mat(P,size(vh,1),size(vh,2))
                       P = fiss_changement_repere(V,size(vb,1),1)
                       vbid2=reshape(V,(/1,dime*dime/))
                       vrtrav((ipg-1)*dime*dime+1:ipg*dime*dime)=vbid2(1,:)

                       ! Dans le repere principal...
                       allocate(vdsiloc(size(vh,1)))
                       allocate(vsi0loc(size(vh,1)))
                       vdsiloc = matmul(P,vdsi)
                       vsi0loc = matmul(P,vsi0)

                       deallocate(vdl0,vdlI,vdlII,vdle,vh,vb,vn)
                       deallocate(vsi0,vdsi,vsi,P)

                       !----- Cas de la loi 12
                       if (iloi==12) then

                           !----- Recuperation des parametres de la loi
                           id = 6
                           if (dime == 3) id = id-1

                           RT = 1.001 * vprel(id)                        ! Contrainte limite
                           C  = 1.001 * vprel(id+2)                      ! Cohesion                            

                           !----- Recuperation de la contrainte principale 1 
                           s1 = vsi0loc(1)+vdsiloc(1)

                           !----- Critere de rupture en ouverture pure
                           vcrit1 = s1 - RT

                           if (vcrit1 > 0.d0) then
                               !- Critere en traction pure
                               MOT = 'ouvert'
                               alpg(ipg) = (RT-vsi0loc(1))/vdsiloc(1)
                           else
                               !- Critere en cisaillement
                               ! MOT = 'cisaille'
                           end if ! Fin du cas loi 12
                       else
                           stop 'FIDES_fiss_pilot : loi non encore programmee'
                       end if

                       deallocate(vsi0loc,vdsiloc)

                    end do

                    !----- On conserve 
                    rupt(ie) = MOT
                    alph(ie) = minval(alpg)
                    !----- On conserve le repere local du point de Gauss le plus dangereux
                    ipg = find_num(alpg,minval(alpg))
                    mrtrav(ie,:)=vrtrav((ipg-1)*dime*dime+1:(ipg*dime*dime))
                    
                    deallocate(vrtrav)
                    deallocate(alpg,ksig)
            end if  ! end if iloi
        end if  ! end if fiss == 1
    end do   
    ! Fin de la boucle d'element

    alpha = minval(alph) 

    if (alpha < 1.d0) then
        iedng = find_num(alph,alpha)
        MOT = rupt(iedng)
        elefiss = count(alph<=0)

        typ = nomtype(ktypel(iedng))

        if (fiss==1 .and. (typ=='MBT3') .or. (typ=='MTT4')) then
           !----- Et on stocke son repère principal
           irloc(iedng,:) = mrtrav(iedng,:)
           inorm(iedng,:) = mrtrav(iedng,1:dime)
        end if
        
        if (alpha <= 0.d0) alpha = 1.d-10
        if (alpha >= 1.d0) alpha = 1.d0

        print'(a11,e11.6,a11,e11.6,a13,i7,a3,a8,a21,i7,2x,a4)','Pilotage : ', &
             &      alpha,' dlam : ',dlam,'   element : ',iedng,'  ',MOT,' [ elts critiques : ',elefiss,' ]'
    else
        alpha = 1.d0
    end if

    !----- On reactualise le champ de deplacements (vduI et facteur de charge)
    vduI = alpha*vduI
    dlam = alpha*dlam

    deallocate(mrtrav)

end subroutine pilot2


subroutine change_etat()

!-----------------------------------------------------------------!
! Fonction d'aiguillage pour la gestion du changement d'etat      !
!-----------------------------------------------------------------!

     use variables, only : interf, fiss, interf_macro
     use element_interface, only : interf_change_etat
     use fissuration, only : fiss_change_etat 

     if (interf == 1) call interf_change_etat
     if (fiss == 1)   call fiss_change_etat 

end subroutine change_etat

!-------------------------------------------------------!

subroutine rupture(ie,ipg,vsig,vnle,wplael,vprel)

!-----------------------------------------------------------------!
! Fonction d'aiguillage pour la gestion du changement d'etat      !
!-----------------------------------------------------------------!

     use variables, only : interf, fiss, interf_macro
     use element_interface, only : interf_rupture
     use fissuration, only : fiss_rupture
     use acier, only : acier_rupture

     real*8, dimension(:,:), intent(inout) :: vsig, vnle
     integer, intent(in) :: ie, ipg
     real*8, dimension(:), intent(in), optional :: vprel
     real*8, intent(in), optional :: wplael

     if (interf == 1) call interf_rupture(ie,ipg,vsig,vnle)
     if (fiss == 1)   call fiss_rupture(ie,ipg,vsig,vnle,wplael,vprel)
     
     call acier_rupture(ie,ipg,vsig,vnle,vprel)
        
end subroutine rupture

!-------------------------------------------------------!

subroutine modul(ie,ipg,vh,vprel)

!-----------------------------------------------------------------!
! Fonction d'aiguillage pour la gestion de la raideur             !
!-----------------------------------------------------------------!

     use variables, only : interf, fiss, interf_macro
     use element_interface, only : interf_modul
     use fissuration, only : fiss_modul

     real*8, dimension(:), intent(in) :: vprel
     real*8, dimension(:,:), intent(inout) :: vh
     integer, intent(in) :: ie, ipg

     if (interf == 1) call interf_modul(ie,ipg,vh,vprel)
     if (fiss == 1)   call fiss_modul(ie,ipg,vh,vprel)

end subroutine modul

!-------------------------------------------------------!

subroutine loi(iloi,icomp,ie,ipg)

!-----------------------------------------------------------------!
! Fonction d'aiguillage pour la gestion du comportement           !
!-----------------------------------------------------------------!

     use variables, only : interf, fiss, interf_macro
     use element_interface, only :interf_loi
     use fissuration, only :fiss_loi

     integer, intent(in) :: icomp, ie, ipg
     integer, intent(out) :: iloi

     iloi = icomp

     if (fiss == 1) call fiss_loi(iloi,icomp,ie,ipg)
     if (interf==1) call interf_loi(iloi,icomp,ie,ipg)

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

     if (interf == 1) call interf_stock()
     if (fiss == 1)   call fiss_stock()

end subroutine stock

!-------------------------------------------------------!

subroutine distal()

!-----------------------------------------------------------------!
! Fonction d'aiguillage pour la definition des proprietes         !
!                     mecaniques aleatoires                       !
!-----------------------------------------------------------------!

     use variables, only : interf, fiss
     use element_interface, only : interf_distal
     use fissuration, only : fiss_distal

     if (interf == 1) call interf_distal        
     if (fiss == 1)   call fiss_distal

end subroutine distal

!-------------------------------------------------------!

end module aiguillage

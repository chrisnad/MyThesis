module fissuration

!********************************************************!
!         Gestion du calcul fissure                      ! 
!********************************************************!

contains

!********************************************************!

subroutine fiss_sig(vsig,vnle,vh,vb,vnl0,vprel,vdle,iloi,ie,ipg)

!********************************************************!
!     Calcul des contraintes pour le calcul fissure      !
!********************************************************!

        use variables, only : dime, inict, irloc, iedng
        use utilitaire, only : princ
        implicit none

        !--- Variables IN
        real*8, dimension(:,:), intent(in) :: vh, vb
        real*8, dimension(:), intent(in) :: vnl0, vdle, vprel
        integer, intent(in) :: iloi, ie, ipg
        
        !--- Variables OUT
        real*8, dimension(:), intent(inout) :: vsig
        real*8, dimension(:), intent(out) :: vnle
        
        !--- Quantites globales
        real*8, dimension(size(vh,1),size(vh,2)) ::  P, Q, vhloc
        real*8 :: vcrit, dlam, V(dime,dime), vsp(dime)
        real*8, dimension(size(vb,1)) :: vepspl, vsigm, vsigi, vsigmloc, vsigiloc, vepsplloc
        real*8, dimension(size(vh,1)) :: vdfdsigloc, vdgdsigloc

        integer :: it, nc1
        character(len=5) :: calcu
                            
!********************************************************!

        vdfdsigloc = 0.d0
        vdgdsigloc = 0.d0

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
        P = fiss_changement_repere(V,nc1,1)
        Q = fiss_changement_repere(V,nc1,2)
        vsigi = 0.d0
                  
        ! Dans le repere principal...
        vhloc = matmul(P,matmul(vh,Q))
        vsigmloc = matmul(P,vsigm)
        vsigiloc = matmul(P,vsigi)
        vepsplloc = matmul(P,vepspl)

        !----- Calcul des contraintes verifiant le critere         
        calcu = 'D1F'      
        call fiss_crit(iloi,vprel,vsigmloc,ie,ipg,vcrit,calcu,vdfdsigloc,vdgdsigloc)
        calcu = 'D0F'
        
        !----- Calcul du multiplicateur plastique (algo semi-implicite)
        it = 0
        vsigiloc = vsigmloc
        dlam = 0.0d0

        ! ajouter la modification de vh dans le cas de cisaillement

         do while (vcrit > 1.d-20)       
            dlam = dlam + vcrit/dot_product(vdfdsigloc,matmul(vhloc,vdgdsigloc))
            vsigmloc = vsigiloc - dlam*matmul(vhloc,vdgdsigloc)      
            call fiss_crit(iloi,vprel,vsigmloc,ie,ipg,vcrit,calcu)
            it = it + 1
            if (it > 50) exit   
        end do

        vepsplloc = vepsplloc + dlam*vdgdsigloc

        !----- On retourne dans le repere global...
        vsigm = matmul(Q,vsigmloc)
        vepspl = matmul(Q,vepsplloc)
   
        !----- Sorties :

        vsig = vsigm
        vnle = vepspl

end subroutine fiss_sig

!********************************************************!

subroutine fiss_crit(iloi,vprel,vsig,ie,ipg,vcrit,calcu,vdfdsig,vdgdsig)

!*****************************************************************!
!     Calcul des criteres de rupture et des derivees              !
!               pour le calcul fissure                            !
!         sorties :  vcrit, vdfdsig, vdgdsig                      !
!                                                                 ! 
! N.B. : ON EST DANS LE REPERE PRINCIPAL (sigma1, sigma2, sigma3) !
! sigma1 > sigma2 > sigma3 !     pas sur !!!!                     !
!*****************************************************************!

        use variables, only : dime, ktypel, nomtype, nelt, pi, &
                           &  ietat, iedng
        use math, only : Jacobi, INV, norme, find_vec
        use initialisation, only : init_vec, init_mat
        implicit none

        !--- Variables IN
        real*8, dimension(:), intent(in) :: vsig
        real*8, dimension(:), intent(in) :: vprel
        integer, intent(in) :: iloi, ie, ipg
        character(len=5), intent(in), optional :: calcu

        !--- Variables OUT
        real*8, intent(out) :: vcrit
        real*8, dimension(:), intent(out), optional :: vdfdsig, vdgdsig
 
        !--- Variables locales procedure
        real*8 :: s1, s3, vcrit1, vcrit2, vcrit2a, vcrit2b, vcrit3, vcrit4, A
        integer :: id, i, j, k
        logical :: ideriv
        character(len=5) :: typ

        !--- Loi de Tresca (case 12)
        real*8 :: RT, C

        !--- Loi de Von Mises (case 11)
        !real*8, dimension(:), allocatable :: vsigd, vun
        !real*8 :: cnu, a, trvsig, vmis

        !--- Loi de Coulomb (case 13)
        real*8, dimension(dime,dime) :: VVsig
        real*8 :: phi, psi, KP, KPS, RP

        !--- Loi de Coulomb "orientee" (case 14)
        real*8, dimension(dime) :: t1, t2, T, Tt
        real*8 :: Tt1, Tt2, sigma, tau

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

        !----- Criteres de rupture 
        vcrit1 = s1 - RT                      ! ... en ouverture pure
        vcrit = -1.d0

        !----- Calcul des criteres et des derivees
        if ((iedng==ie) .or. (ietat(ie) /=0 )) then
        
            !----- Rupture en ouverture
            if (vcrit1 >= 0.d0) then
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
            
    case(13)     !----- Loi de Tresca - Coulomb

        !----- Recuperation des parametres
        id = 6
        if (dime == 3) id = id-1
        RT = vprel(id)                        ! Contrainte limite
        !C  = 5. * vprel(id)                  ! Cohesion
        C  = vprel(id+1)                      ! Cohesion
        
        phi = pi*vprel(id+2)/180.             ! Angle de frottement
        psi = atan(RT/C)
       
        if (ietat(ie) /=0 ) RT = 0.0d0

        !----- Recuperation des contraintes principales 1 et 3                                               
        s1 = vsig(1)
        s3 = vsig(2)
                                                
        sigma = s1                               
        tau   = (s1-s3)/2.

        !----- Criteres de rupture 
        vcrit1  = sigma - RT                         ! ... en ouverture pure
        vcrit2a = abs(tau) - C                       ! ... en cisaillement maximal
        vcrit2b = abs(tau) + sigma * tan(phi)        ! ... Mohr-Coulomb
        vcrit3  = abs(tau) + sigma * tan(psi-pi/2.)
        A = C*(1-(tan(psi-pi/2.)/tan(phi)))
        vcrit4  = vcrit3 - A

        vcrit = -1.d0

        !----- Calcul des criteres et des derivees
        if ((iedng==ie) .or. (ietat(ie) /=0 )) then
           if ((vcrit3 < 0.d0) .or. (ietat(ie) < 0)) then
           
               vcrit2 = vcrit2a
               if (ietat(ie) /=0 ) vcrit2 = vcrit2b

               ! Critere en ouverture pure
               if ((vcrit1 > 0.d0) .or. (vcrit2 > 0.d0)) then
             
                 vcrit = vcrit1
                 if (ideriv) then                 
                    ietat(ie) = 1                    
                    if (dime == 2) then
                        vdfdsig = (/ 1.d0, 0.d0, 0.d0 /) 
                        vdgdsig = (/ 1.d0, 0.d0, 0.d0 /)
                        
                    elseif (dime == 3) then
                        vdfdsig = (/ 1.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0 /) 
                        vdgdsig = (/ 1.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0 /)
                    end if
                 end if                
               end if

           elseif ((vcrit3 > 0.d0) .and. (vcrit4 < 0.d0)) then
           
               vcrit2 = vcrit2a
               if (ietat(ie) /=0 ) vcrit2 = vcrit2b

               ! Critere en cisaillement
               if ((vcrit2 > 0.d0) .or. (vcrit1 > 0.d0)) then
                   
                  vcrit = abs(tau) + sigma * tan(phi)

                  if (ideriv) then
                    ietat(ie) = 3

                    if (dime == 2) then
                       ! dans le repere local...
                       vdfdsig = (/ tan(phi) + 0.5d0, -0.5d0, 0.d0 /)
                       if (vcrit2==vcrit2a) vdfdsig = (/ 0.5d0, -0.5d0, 0.d0 /)
                       vdgdsig = (/ tan(psi) + 0.5d0, -0.5d0, 0.d0/)

                    elseif (dime == 3) then
                       vdfdsig = (/ tan(phi) + 0.5d0, -0.5d0, 0.d0 , 0d0 , 0d0, 0d0 /)
                       if (vcrit2==vcrit2a) vdfdsig = (/ 0.5d0, -0.5d0 , 0.d0, 0d0 , 0d0, 0d0 /)
                       vdgdsig = (/ tan(psi) + 0.5d0, -0.5d0, 0.d0 ,  0d0 , 0d0, 0d0 /)    
                    end if
                  end if
               end if
           else

               vcrit2=vcrit2a
               ! Critere de Tresca
               if ((vcrit2 > 0.d0) .or. (vcrit1 > 0.d0)) then  
                  vcrit = abs(tau) - C

                  if (ideriv) then
                    ietat(ie) = 3
                    if (dime == 2) then
                       ! dans le repere local...
                       vdfdsig = (/      0.5d0      , -0.5d0, 0.d0 /)
                       vdgdsig = (/ tan(psi) + 0.5d0, -0.5d0, 0.d0/)
                       
                    elseif (dime == 3) then
                       vdfdsig = (/      0.5d0      , -0.5d0 , 0.d0 , 0d0 , 0d0 , 0d0 /)
                       vdgdsig = (/ tan(psi) + 0.5d0, -0.5d0 , 0.d0 , 0d0 , 0d0 , 0d0 /)    
                    end if
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

function fiss_changement_repere(V,n,ind) result (M)

!********************************************************!
!  Definition de la matrice de changement de repere :    !
!     du repere global au repere principal               !
!********************************************************!

        use variables, only : dime
        use math, only : norme
        implicit none

        real*8, dimension(dime,dime), intent(in) :: V
        integer, intent(in) :: n, ind

        real*8, dimension(n,n) :: M
        real*8 :: a, b, c, d, e, f, g, h, i

!********************************************************!
        a = V(1,1)
        b = V(2,1)

        d = V(1,2)
        e = V(2,2)
       
        if (dime == 3) then 
            c = V(3,1)
            f = V(3,2)
            g = V(1,3)
            h = V(2,3)
            i = V(3,3)
        end if

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
            end if

        !----- Calcul de la matrice M=Q : passage repere principal => repere global
        elseif (ind==2) then

            if (dime == 2) then
                M(1,:) = (/ (a**2), (d**2), (2.d0*a*d) /)
                M(2,:) = (/ (b**2), (e**2), (2.d0*b*e) /)
                M(3,:) = (/ (a*b) , (d*e), (b*d+a*e)  /)
            elseif (dime == 3) then
                M(1,:) = (/ (a**2), (d**2), (g**2), (2.d0*a*d), (2.d0*d*g), (2.d0*g*a) /)
                M(2,:) = (/ (b**2), (e**2), (h**2), (2.d0*b*e), (2.d0*e*h), (2.d0*h*b) /)
                M(3,:) = (/ (c**2), (f**2), (i**2), (2.d0*c*f), (2.d0*f*i), (2.d0*i*c) /)
                M(4,:) = (/ (a*b) , (d*e) , (g*h),  (b*d+a*e) , (e*g+d*h) , (g*b+a*h)  /)
                M(5,:) = (/ (b*c) , (e*f) , (h*i),  (c*e+b*f) , (h*f+e*i) , (h*c+b*i)  /)
                M(6,:) = (/ (a*c) , (d*f) , (i*g),  (c*d+a*f) , (g*f+d*i) , (g*c+a*i)  /)
            end if
            
        end if        

end function fiss_changement_repere
!********************************************************!

!********************************************************!

subroutine fiss_pilot(imetpilo,vsol,vduI,vduII,dlam,alpham,MOTm,elefiss)

!********************************************************!
!  Pilotage du calcul : gestion du facteur de chargement !
!********************************************************!

        use variables, only : dime, nelt, idmax, ktypel, nomtype, infele, &
                    & vprelg, kprop, kprop0, ietat, iedngm, fiss, pi,     &
                    & mrtrav, vrtrav, inorm, irloc, kloce
        use initialisation, only : init_vec, init_mat
        use lib_elem, only : elem_B, elem_kloce2, elem_hooke
        use math
        use utilitaire, only : princ
        implicit none

        ! Variables IN
        real*8, dimension(:), intent(in) :: vduII, vsol
        integer, intent(in) :: imetpilo

        ! Variables IN-OUT
        real*8, dimension(:), intent(in) :: vduI
        real*8, intent(in) :: dlam 
        real*8, intent(out) :: alpham
        character(len=8), intent(inout) :: MOTm
        integer, intent(inout) :: elefiss

        ! Variables locales 
        real*8, dimension(:), allocatable :: vdle, vdl0, vsi0, vdsi, vsi, &
                                        & vsi0loc, vdsiloc, alpg, vdlI, vdlII, ksigPG
        real*8, dimension(:,:), allocatable :: vn, vb, vh, ksig, P        
        real*8 :: vprel(idmax), alph(nelt), vbid(dime), vbid2(1,dime*dime), V(dime,dime)
        real*8 :: alpc, sign0, C, RT, phi, psi, sigma, tau, &
                  & vcrit1, vcrit2, vcrit3, detj, signe, coef, &
                  & s1, s2, s3, s01, s02, s03, &
                  & sigma0, dsigma, tau0, dtau
        character(len=5) :: typ
        character(len=8) :: MOT            
        character(len=8), dimension(nelt) :: rupt
        integer :: i, j, k, ie, ipg, ipgd, npg, id, iloi, ndle
        
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
            do ie = 1, nelt
            
               MOT = '  '
               typ = nomtype(ktypel(ie))
            
               !----- Pour les elements massifs vierges
               if ((typ=='MBT3'.or.typ=='MBQ4'.or.typ=='MBT6'.or.typ=='MTT4'.or. &
                  & typ=='MTP6'.or.typ=='MTH8') .and. (ietat(ie)==0)) then

                  !----- Proprietes elementaires
                  vprel = vprelg(kprop(ie),1:idmax)

                  !----- Recuperation des parametres de la loi
                  iloi = int(vprel(1))

                  if (iloi==1 .or. iloi==11 .or. iloi==201) goto 2000
                  
                  !-----  Recuperation des informations sur les elements
                  allocate(ksig(size(infele(ktypel(ie))%Q,1), size(infele(ktypel(ie))%Q,2)))
                  ksig = infele(ktypel(ie))%Q
                  npg = size(ksig,1)

                  ! Tenir compte le point Gauss au centre
                  allocate(ksigPG(size(infele(ktypel(ie))%Q,2)))
                  !ksigPG = sum(ksig,1) / size(ksig,1)
                  !npg = 1

                  !----- Initialisations relatives aux elements
                  call init_vec(vrtrav,npg*dime*dime)
                  allocate(alpg(npg))						
                  alpg = 1.d0
                        
                  do ipg = 1, 1
                    !- Detection de la premiere fissuration (ouverture)
                    !ipg = 1
                           
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
                    !call elem_B(vn,vb,detj,ksigPG,ie)    ! tenir compte le point Gauss au centre
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
                    if (iloi==1) then
                    
                    elseif (iloi==12) then

                       !----- Recuperation des parametres de la loi
                       id = 6
                       if (dime == 3) id = id-1
                                 
                       RT = 1.001 * vprel(id)                        ! Contrainte limite
                                                                  
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
                       end if
                       ! Fin du cas loi 12
                    else
                       stop 'FIDES_fiss_pilot : loi non encore programmee'
                    end if

                    deallocate(vsi0loc,vdsiloc)
                    
                  end do

                  !----- On conserve 
                  rupt(ie) = MOT
                  alph(ie) = minval(alpg)

                  !----- On conserve le repere local du point de Gauss le plus dangereux
                  ipgd = find_num(alpg,minval(alpg))
                  mrtrav(ie,:)=vrtrav((ipgd-1)*dime*dime+1:(ipgd*dime*dime))

                  deallocate(vrtrav)
                  deallocate(alpg,ksig,ksigPG)
                  
           2000   continue                                    
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


!********************************************************!

subroutine fiss_change_etat()

!********************************************************!
!                 Gestion du changement d'etat           !
!                  pour le calcul fissure                !
!********************************************************!

        use variables, only : nelt, vsol, ietat, detoscill, comptoscill, icompfiss, ouvermax, &
             & iouver1, iouver2, histetat1, histetat2, fiss, kloce, infele, ktypel, vprelg, nomtype, kprop
        use initialisation, only : init_vec
        use lib_elem, only : elem_kloce2
        implicit none

        real*8, dimension(:), allocatable :: vdle        
        real*8 :: ouv
        integer :: ie, oscill, ietat0, ietat1, ietat2, ndle, ipg, npg

!********************************************************!

        if (fiss==1) then
            do ie = 1, nelt                       ! Boucle sur les elements

                if (ietat(ie)/=0 .and. vprelg(kprop(ie),1)==icompfiss) then ! Element fissure

                    if (detoscill(ie)==6) then
                        ietat(ie) = -2
                        oscill=1
                        !print*, 'Oscillation de l''element',ie
                        detoscill(ie) = detoscill(ie)+1                         
                        comptoscill = comptoscill+1
                    end if

                    !----- Vecteur kloce de localisation pour assemblage
                    call elem_kloce2(ie,ndle)
                    allocate(vdle(ndle))

                    !----- Deplacement total
                    vdle = vsol(kloce(1:ndle))

                    !-----  Recuperation de l'ouverture de fissure
                    ouv = fiss_ouv(ie)

                    if ((ietat(ie)==1) .or. (ietat(ie)==2)) then
                        oscill = 0

                        iouver2(ie) = iouver1(ie)
                        iouver1(ie) = ouv
                     
                        ietat2 = histetat2(ie)
                        ietat1 = histetat1(ie)
                        ietat0 = ietat(ie)

                        if ((ietat0==ietat2) .and. (ietat0/=ietat1)) then
                            oscill=1
                            ietat(ie) = -2
                            !print*, 'Oscillation'
                            comptoscill=comptoscill+1
                        end if
                    
                        histetat2(ie) = histetat1(ie)
                        histetat1(ie) = ietat0

                        if ((ietat(ie)==1) .and. ouv > ouvermax) then
                           ietat(ie) = 2

                        elseif ((ietat(ie)==1) .and. (ouv <= 0.d0) .and. (oscill == 0)) then
                        !if ((ietat(ie)==1) .and. (ouv <= 0.d0) .and. (oscill == 0)) then
                           ietat(ie) = 2
                           detoscill(ie)=detoscill(ie)+1  
                        elseif ((ietat(ie)==2) .and. (ouv > 0.d0) .and. (oscill == 0)) then
                           ietat(ie) = 1
                           detoscill(ie)=detoscill(ie)+1  
                        end if
                    end if
                    deallocate(vdle)
                end if
            end do
        end if

end subroutine fiss_change_etat

!********************************************************!

function fiss_ouv(ie) result (ouv)

!********************************************************!
!          Calcul de l'ouverture moyenne d'une fissure   !
!********************************************************!

        use variables, only : dime, inorm, vcor, kconec, infele, ktypel, vsol, infnoe
        use math, only : norme                                
        
        implicit none
        integer, intent(in) :: ie

        real*8 :: XG, YG, ZG, AN, BN, CN, DN, DENO, dep1, dep2, ouv
        real*8, dimension(dime) :: nf, NN, PG

	    real*8, dimension(infele(ktypel(ie))%nnel) :: dist, depNorm
	    real*8, dimension(dime,infele(ktypel(ie))%nnel) :: depNds, vcore

        integer :: i, nnzp, nnzn, noel, ino, ndln
        integer, dimension(infele(ktypel(ie))%nnel,infele(ktypel(ie))%ndln) :: ddl

        logical, dimension(infele(ktypel(ie))%nnel) :: cotep
        
!*********************************************************
     
        if (dime==1) stop 'FIDES_fiss_ouv : cas dime = 1 non prevu !'      
        
        nf = inorm(ie,:)
        NN = nf/norme(nf)

        noel = infele(ktypel(ie))%nnel
        ndln = infele(ktypel(ie))%ndln

        vcore = vcor(:,kconec(ie,1:noel))
        
        do ino = 1, noel
            ddl(ino,:) = infnoe(kconec(ie,ino))%dln(1:ndln)
        end do

        depNds(1,:) = vsol(ddl(:,1))
        depNds(2,:) = vsol(ddl(:,2))
            
        if (dime==3) then
            depNds(3,:) = vsol(ddl(:,3))
        end if
                 
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

end function fiss_ouv

!********************************************************!

subroutine fiss_rupture(ie,ipg,vsig,vnle,wplael,vprel)

!********************************************************!
!           Fonction de detection de la rupture          !
!                pour le calcul fissure                  !
!********************************************************!

        use math, only : norme, find_vec
        use variables, only : nelt, calco, dime, ktypel, infele, fiss, ietat, inorm, &
                                & icompfiss, iedng
        use initialisation, only : init_vec, init_mat
        use utilitaire, only : princ
        use lib_elem, only : elem_B
        implicit none

        real*8, dimension(:,:), intent(inout) :: vnle, vsig
        integer, intent(in) :: ie, ipg
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
        real*8, dimension(dime,dime) :: vsigPGmat, V, VVsig
        real*8, dimension(dime) :: vsp, nf

        real*8 :: wplamax, detj, vspmax

        integer :: i, id, nnel, ncgr, nc1, iloi, NROT, itype

!********************************************************!

        itype = ktypel(ie)

        ! Recuperation des informations sur les elements
        call init_mat(ksin, size(infele(itype)%ksin,1), size(infele(itype)%ksin,2))
        ksin = infele(itype)%ksin
        nnel = infele(itype)%nnel
        nnel = infele(itype)%ncgr

        ! Recuperation des points d'integration
        call init_vec(vpg, size(infele(itype)%W,1))
        vpg = infele(itype)%W
        
        iloi=vprel(1)
          
        if ((fiss==1) .and. (iloi==icompfiss)) then
            if (iloi==12) then
                id=7
                if (dime==3) id=id-1
                wplamax=vprel(id)
            else
                stop 'Calcul impossible avec cette loi de comportement'
            end if

            !if ((wplael>wplamax) .and. (ietat(ie)==0)) then
            if (ie==iedng) then
            !if (ie==iedng .or. (ietat(ie)==0) ) then
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

!********************************************************!

subroutine fiss_modul(ie,ipg,vh,vprel)

!********************************************************!
!                Gestion de la raideur                   !
!********************************************************!

        use variables, only : dime, ietat, icompfiss, smear, depen
        implicit none

        real*8, dimension(:), intent(in) :: vprel        
        real*8, dimension(:,:), intent(inout) :: vh
        integer, intent(in) :: ie, ipg
        integer :: iloi        
        
!********************************************************!

        iloi = vprel(1)

        if (iloi==icompfiss) then

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

!********************************************************!

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
                stop 'Cas non prevu'
        end if

        ! Matrice Def = vh * N*inv(Dcr + N'*vh*N)*N' * vh;
        M0 = matmul(transpose(M),matmul(vh,M))
        M1 = INV(Dcr + M0)
        M2 = matmul(M,matmul(M1,transpose(M)))
        Def = matmul(vh,matmul(M2,vh))
        vh = vh - Def                       

end subroutine fiss_modul_smeared

!********************************************************!

subroutine fiss_loi(iloi,icomp,ie,ipg)

!========================================================!
!     Calcul des contraintes pour calcul fissure         !
!========================================================!

    use variables, only : ietat, icompfiss
    implicit none

    integer, intent(in) :: icomp, ie, ipg
    integer, intent(inout) :: iloi

!********************************************************!

        if (icomp==icompfiss) then
            if ((ietat(ie)==1) .or. (ietat(ie) == -1) .or. (ietat(ie) == -2)) then
            ! Element fissure => aiguillage vers calcul elastique avec raideur nulle        
                iloi=1
            elseif (ietat(ie)==2) then
                ! Element (fissure) referme => aiguillage vers calcul loi de Coulomb
                !iloi = 14
                 iloi = 1                ! Loi elastique
            else
                ! Element vierge ou en cours de fissuration => aiguillage vers calcul elasto-plastique
                iloi = 12
            end if
        else
               iloi = icomp        
        end if

end subroutine fiss_loi

!********************************************************!

subroutine fiss_init()

!********************************************************!
!  Initialisations liees au calcul fissure               !
!********************************************************!
    
    use variables, only : vprelg, nelt, dime, ietat, inorm, irloc, iouver1, iouver2, &
                          & histetat1, histetat2, icompfiss, interf_macro, &
                          & fiss, smear, depen, iprint, iedngm, ietatma, nomtype, infele, varg0, vargs, endo, ktypel
    use initialisation, only : init_mat, init_vec
    implicit none

    integer, dimension(size(nomtype,1)) :: npg
    integer, dimension(3) :: test_ind
    integer :: i, npgm, itype, ntype
    logical ismember
    character reply
    character(len=5) :: typ   

!********************************************************!

    fiss = 0
    smear = 0
    depen = 0
                
    icompfiss = 12   ! Numero de loi de comportement pour beton probabiliste
 
    ismember = .false.
    reply = ''            
  
    do i = 1, size(vprelg,1) 
        ! Cas particulier de la loi traduisant la fissuration du beton
        if (vprelg(i,1)==icompfiss) then
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
                   fiss = 1
                   smear = 1
                   if (iprint>0) print*, ' Ok, on va faire du SMEARED ! '
               case ('d')
                   fiss = 1
                   depen = 1
                   if (iprint>0) print*, ' Ok, on va faire de la DEPENALISATION ! '
            end select
        end do
    end if
 
    test_ind = (/ fiss, smear, depen /)
    if (all(test_ind /= 0)) stop 'Erreur, verifier les index fiss, depen !!'    
    call init_vec(ietat,nelt)

    ntype = maxval(ktypel)
    do itype = 1, ntype
       typ = nomtype(itype)
       if (typ=='MBT3'.or.typ=='MBQ4'.or.typ=='MBT6'.or. &     ! element 2D et ! element 3D
         & typ=='MTT4'.or.typ=='MTP6'.or.typ=='MTH8'.or.typ=='MBI4') then
         npg(itype) = size(infele(itype)%W,1)
       end if 
    end do

    if (fiss == 1) then        
        iedngm = 0
        call init_vec(iouver1,nelt)
        call init_vec(iouver2,nelt)
        call init_vec(histetat1,nelt)
        call init_vec(histetat2,nelt)
    end if

    if (fiss == 1 .or. interf_macro==1) then
        call init_mat(inorm,nelt,dime)
        call init_mat(irloc,nelt,dime*dime) ! JLT pour tenir compte de tout le repere
    end if

end subroutine fiss_init


!********************************************************!

subroutine fiss_stock()

!********************************************************!
!   Stockage des resultats à chaque pas de temps pour    !
!   le  modele  de  fissuration  avec  gestion  des      !
!   transitions etat(-1)--> etat(1), etat(-2)--> etat(2) !
!********************************************************!
    
    use variables, only :  fiss, ietat, nelt, kprop, vprelg, icompfiss
    implicit none
    integer :: ie

! *******************************************************!

    if (fiss==1) then
       do ie = 1, nelt
         if (ietat(ie)==-2 .and. vprelg(kprop(ie),1)==icompfiss) ietat(ie)=1
       end do
    end if

end subroutine fiss_stock

!********************************************************!

subroutine fiss_distal()

!---------------------------------------------------------------!
! Fonction de definition des proprietes mecaniques aleatoires   !
!   Donnees en entree :                                         !
!   - kprop : numerotation des groupes des elements             !
!   - vprelg : proprietes mecaniques (fonction du modele)       !
!                                                               !
!   Valeurs en sortie :                                         !
!   - kprop : actualise en fonction du numero de l'element      !
!             (1 element prob = 1 groupe)                       !
!   - vprelg : modifie pour prendre en compte la propriete      !
!      aleatoire consideree - Attention depend du modele choisi !
!---------------------------------------------------------------!
        
    use variables, only : nelt, dime, kconec, ktypel, nomtype, &
            & infele, vprelg, kprop, kprop0, interf, idmax, &
            !& prob, pi, fiss, icompfiss, Dg, fc, &
            & pi, fiss, icompfiss, Dg, fc, vitrav, &
            !--- gestion de l'alea (deb)
            !& alea,iloip, param, iparam
            & alea, young, resist, enerj, ietat
            !--- gestion de l'alea (fin)
    use lib_elem, only : elem_B
    use initialisation
    use proprietes_aleatoires
    use math

    implicit none
    
    real*8, dimension(:,:), allocatable :: vprelg0, ksig, vn, vb
    real*8, dimension(:), allocatable :: vpg, RR, EE, WW
    real*8, dimension(nelt) :: Vel
    real*8, dimension(1) :: Rt, C, E, W   
    real*8, dimension(5) :: loi1, loi2, loi3
    real*8 :: detj, epais, rp, &
              & poids, Ve, Vg, r, Vem, A, B, CD, b0, c0, iran(1) !, &
              !& alpha   ! alpha est utilise pour la determination des parametres du modele elasto-plastique. AM
    integer :: iebe, i, ie, ie0, nbg, ipg, npg, id, iloi, dPlay, igrp, igmp, ignp, ine, iefiss, ira
    logical :: resp, modp, nrjp, iprobVO, iprobFI, iprobEN, iprobMY

!********************************************************!

    !----- Preparation de la renumerotation de kprop et vprelg  !
    ! Principe : on conserve la numerotation pour tous les      !
    ! groupes dont l'ancien numero est inferieur a nga. Puis on !
    ! "decale" d'un cran dans la numerotation tous les groupes  !
    ! qui ont un ancien numero plus grand que nga. Enfin, tous  !
    ! les elements du groupe nga sont renumerote en nga+i.      !
    
    if ((fiss==1) .and. (alea)) then !---Pour le modèle de fissuration aleatoire
        rp = 1.e15
        ine = 0
        nbg = maxval(kprop)
        call init_mat(vprelg0,(size(kprop,1)+nbg),size(vprelg,2))
        vprelg0(1:nbg,:) = vprelg(1:nbg,:)

        Vel = 0.0
        call init_vec(RR,nelt)
        call init_vec(EE,nelt)
        call init_vec(WW,nelt)

        ie0 = 0
        iebe = 0
        iefiss = 0
        
        !----- Indices pour trace de distributions de proprietes aleatoires
        iprobVO = .false.
        iprobFI = .false.
        iprobEN = .false.
        iprobMY = .false.

        !-- Boucle sur les elements
        do ie = 1, nelt
        
            !-- detection d'un element probabiliste (module et/ou resistance et/ou energie)
            resp = .false.
            modp = .false.
            nrjp = .false.
            if (allocated(resist)) then
                if (count(resist%num==kprop0(ie))==1) resp=.true.
            end if
            if (allocated(young)) then
                if (count(young%num==kprop0(ie))==1) modp=.true.
            end if
            if (allocated(enerj)) then
                if (count(enerj%num==kprop0(ie))==1) nrjp=.true.
            end if
            
            !-- si l'element appartient a un groupe probabiliste (module et/ou resistance et/ou energie)
            if (resp .or. modp .or. nrjp) then
            
                !-- on verifie la compatibilite de la loi de comp avec sa formulation probabiliste
                iloi = vprelg0(kprop(ie),1)
                if (iloi==icompfiss .or. iloi==201 .or. iloi==202 .or. iloi==203) then   ! Loi du modele d'interface macro
            
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
                    if (dime == 2) epais = vprelg0(kprop(ie),idmax-2)

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
                        loi1(5) = resist(igrp)%loi      ! num

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
                            
                                if (r == 0.d0) then 
                                   print*, ie, 'volume nul'
                                   !stop
                                   r = 1.
                                end if

                                if (r < 0.d0) then 
                                   print*, ie, 'volume negative'
                                   !stop
                                   r = 1.
                                end if
                                
                                do while ((r< 0.004).or.(r>1000.))
                                   if (r<0.004) r = r*10.d0
                                   if (r>1000.) r = r/10.d0
                                end do

                                if ((r<.004d0).or.(r>1000.d0)) stop 'Sortie du domaine de validite du modele'
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
                             stop ' distal_fiss : Resistance negative'
                        end if
                        RR(iebe) = Rt(1)

                        !iran = randgp(1,(/ie0/))*10000

                        !if ((iefiss < 2000) .and. (iran(1) > 500) .and. (mod(int(iran(1)),2)==0)) then
                        !   ietat(ie) = 1
                        !   iefiss = iefiss + 1
                        !end if
                        
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
                        WW(iebe) = W(1)*Ve
                    end if   !fin if nrjp

                    !---- Stockage des valeurs
                    vprelg0(nbg + ie0,:) = vprelg(kprop(ie),:)
                    kprop(ie) = nbg + ie0

                    !----- Stockage resistance a la traction
                    id = 4
                    if (dime==3) id = id-1

                    if (resp) vprelg0(kprop(ie),id+2) = Rt(1)

                    !----- Stockage modules elastiques
                    if (modp) vprelg0(kprop(ie),id) = E(1)

                    !----- Stockage de l'energie
                    if (iloi==icompfiss) then
                       if (nrjp) vprelg0(kprop(ie),id+3) = W(1) * Ve   ! Energie de deformation plastique
                    end if

                end if
            end if
        end do

        !-- Sorties :
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
            !do i=1,size(resist%num)
                call alea_trac(RR,loi1,200,dPlay,'resistance')
            !end do
        end if
        
        ! Visualisation de la distribution de module elastique        
        if (iprobMY) then
            dPlay = 1
            call alea_trac(EE,loi3,200,dPlay,'module')
        end if
    
        ! Visualisation de la distribution de l'energie
        if (iprobEN) then  
            dPlay = 1
            call alea_trac(WW,loi2,200,dPlay,'energie')
        end if
        
        deallocate(vprelg0,EE,RR,WW)

        print*, 'Le rapport du volume le plus petit : ', rp
        print*, 'Elements prefissures : ',  iefiss
        
    end if

end subroutine fiss_distal
!-------------------------------------------------------!

end module fissuration

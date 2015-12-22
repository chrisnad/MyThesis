module acier

contains

!********************************************************!

subroutine acier_sig(vhep,vsig,vnle,vh,vb,vnl0,vprel,vdle,iloi,ie,ipg)

!********************************************************!
!                Calcul des contraintes                  !
!********************************************************!

        use variables, only : dime, inict, irloc, iedng
        use initialisation, only : init_vec, init_mat
        use utilitaire, only : princ
        use math, only : vecpro        
        implicit none

        !--- Variables IN
        real*8, dimension(:,:), intent(in) :: vh, vb
        real*8, dimension(:), intent(in) :: vnl0, vdle, vprel
        integer, intent(in) :: iloi, ie, ipg
        
        !--- Variables OUT
        real*8, dimension(:), intent(inout) :: vsig
        real*8, dimension(:), intent(out) :: vnle
        real*8, dimension(:,:), intent(out) :: vhep
        
        !--- Quantites globales
        real*8, dimension(size(vh,1),size(vh,2)) ::  P, Q
        real*8, dimension(dime,dime) :: V 
        real*8, dimension(dime) :: vsp  
        integer :: it, nc1
        character(len=5) :: calcu 

        !--- Quantites principales
        real*8, dimension(size(vh,1),size(vh,2)) :: vhloc
        real*8 :: vcrit, dlam
        real*8 :: f, plastic_rate, dqdep
                         
        ! Quantites principales
        real*8, dimension(size(vb,1)) :: vepspl, vsigm, vsigi
        real*8, dimension(size(vh,1)) :: vsigmloc, vsigiloc, vepsplloc
        real*8, dimension(size(vh,1)) :: vdfdsigloc, vdgdsigloc
        real*8, dimension(size(vh,1)) :: vdfdq, dedep
                                    
!********************************************************!

        !----- Recuperation des grandeurs de non-linearite
        vepspl = vnl0(1:size(vb,1))            ! deformations plastiques initiales
        nc1 = size(vb,1)
        
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
        P = acier_changement_repere(V,nc1,1)
        Q = acier_changement_repere(V,nc1,2)
        vsigi = 0.d0
                  
        ! Dans le repere principal...
        vhloc = matmul(P,matmul(vh,Q))
        vsigmloc = matmul(P,vsigm)
        vsigiloc = matmul(P,vsigi)
        vepsplloc = matmul(P,vepspl)
        plastic_rate = sqrt(dot_product(vepsplloc,vepsplloc))
        
        dedep = 0.
        if (plastic_rate /= 0.) dedep = vepsplloc / plastic_rate
        
        !----- Calcul des contraintes verifiant le critere         
        calcu = 'D1F'      
        call acier_crit(iloi,vprel,vsigmloc,ie,ipg,vcrit,calcu,plastic_rate,dqdep,vdfdsigloc,vdgdsigloc,vdfdq)
        calcu = 'D0F'

        !----- Calcul du multiplicateur plastique (algo semi-implicite)
        it = 0
        vsigiloc = vsigmloc
        dlam = 0.0d0
        f = 0.d0
        vhep = vhloc

        do while (vcrit > 1.d-20)
            f = dot_product(vdfdsigloc,matmul(vhloc,vdgdsigloc))! - dot_product(vdfdq,dedep)*dqdep
            dlam = dlam + vcrit / f
            vsigmloc = vsigiloc - dlam*matmul(vhloc,vdgdsigloc)            
            call acier_crit(iloi,vprel,vsigmloc,ie,ipg,vcrit,calcu,plastic_rate,dqdep)
            it = it + 1
            if (it > 50) exit            
            !vhep = vhloc - vecpro(matmul(vhloc,vdgdsigloc),matmul(vdfdsigloc,vhloc)) / f                       
        end do
                  
        vepsplloc = vepsplloc + dlam*vdgdsigloc

        !----- On retourne dans le repere global...
        vsigm = matmul(Q,vsigmloc)
        vepspl = matmul(Q,vepsplloc)
        vhep = matmul(Q,matmul(vhep,P))
                        
        !----- Sorties :
        vsig = vsigm
        vnle = vepspl

end subroutine acier_sig

!********************************************************!

subroutine acier_crit(iloi,vprel,vsig,ie,ipg,vcrit,calcu,plastic_rate,dqdep,vdfdsig,vdgdsig,vdfdq)

!*****************************************************************!
!     Calcul des criteres de rupture et des derivees              !
!         sorties :  vcrit, vdfdsig, vdgdsig                      !
!                                                                 ! 
! N.B. : ON EST DANS LE REPERE PRINCIPAL (sigma1, sigma2, sigma3) !
! sigma1 > sigma2 > sigma3                                        !
!*****************************************************************!

        use variables, only : dime, ktypel, nomtype, pi, &
                         &  ietat, limels, limrup, ecroui, iter

        use initialisation, only : init_vec, init_mat
        implicit none

        !--- Variables IN
        real*8, dimension(:), intent(in) :: vsig
        real*8, dimension(:), intent(in) :: vprel
        integer, intent(in) :: iloi, ie, ipg
        real*8, intent(in), optional :: plastic_rate        
        character(len=5), intent(in), optional :: calcu

        !--- Variables OUT
        real*8, intent(out) :: vcrit
        real*8, intent(out), optional :: dqdep
        real*8, dimension(:), intent(out), optional :: vdfdsig, vdgdsig, vdfdq
        
        !--- Variables locales procedure
        real*8 :: vcrit1
        logical :: ideriv
        character(len=5) :: typel

        !--- Loi de Von Mises (case 11)
        real*8 :: cnu, trvsig, vmis, vun(size(vsig)), vsigd(size(vsig)), an, maxsig, minsig, RT
        integer :: id
        
!********************************************************!

       !----- Test(s) de compatibilite des arguments d'entree
       ideriv = .false.

       if (calcu == 'D1F') then
            ideriv = .true.
            vdfdsig = 0.0d0
            vdgdsig = 0.0d0
            vdfdq = 0.d0
            dqdep = limels(ie)
       end if

       !----- Recuperation des parametres de la loi
       id = 6
       if (dime == 3) id = id-1
       
       if (plastic_rate==0.) then
           RT = limels(ie)       
       else
           if (ecroui(ie)==0) then
              RT = limels(ie)*(1. + 0.61*dexp(0.475*dlog10(plastic_rate**1.7)))
           else
              RT = limels(ie)*(1. + ecroui(ie)*plastic_rate)
           end if
       end if
       
       if (RT > limrup(ie)) RT = limrup(ie)
       
       vcrit = -1.d0

       maxsig = maxval(abs(vsig))
       minsig = minval(vsig)

       if (maxsig <= RT) then
          goto 10000
       end if

       !----- Calcul des derivees du critere
       if (maxsig == -minsig) then            ! compression  :  vsig(2) <-- max, signe -

          !if (abs(vsig(2)) > limrup(ie)) then
          !   ietat(ie) = 1
          !   goto 10000
          !end if
                 
          vcrit1 = abs(vsig(2)) - RT
                    
          if (vcrit1 >= 0.d0) then
             vcrit = vcrit1
             if (ideriv) then
                if (dime == 2) then
                   vdfdsig = (/ 0.d0, -1.d0, 0.d0 /)
                   vdgdsig = (/ 0.d0, -1.d0, 0.d0 /)
                   vdfdq   = (/ -1.d0, 0.d0, 0.d0 /)
                elseif (dime == 3) then
                   vdfdsig = (/ 0.d0, -1.d0, 0.d0, 0.d0, 0.d0, 0.d0 /) 
                   vdgdsig = (/ 0.d0, -1.d0, 0.d0, 0.d0, 0.d0, 0.d0 /)
                   vdfdq   = (/-1.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0 /)
                end if
             end if
          end if
       
       else              ! traction  :   vsig(1)  <-- max , signe +

          !if (abs(vsig(1)) > limrup(ie)) then
          !   ietat(ie) = 1
          !   goto 10000
          !end if

          vcrit1 = abs(vsig(1)) - RT
          if (vcrit1 >= 0.d0) then
             vcrit = vcrit1
             if (ideriv) then
                if (dime == 2) then
                   vdfdsig = (/ 1.d0, 0.d0, 0.d0 /)
                   vdgdsig = (/ 1.d0, 0.d0, 0.d0 /)
                   vdfdq   = (/-1.d0, 0.d0, 0.d0 /)
                elseif (dime == 3) then
                   vdfdsig = (/ 1.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0 /) 
                   vdgdsig = (/ 1.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0 /)
                   vdfdq   = (/-1.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0 /)
                end if
             end if
          end if
       end if

10000  continue

end subroutine acier_crit

!********************************************************!

function acier_changement_repere(V,n,ind) result (M)

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

end function acier_changement_repere
!********************************************************!

subroutine acier_rupture(ie,ipg,vsig,vnle,vprel)

!********************************************************!
!           Fonction de detection de la rupture          !
!********************************************************!

        use variables, only : ietat

        implicit none

        real*8, dimension(:,:), intent(inout) :: vnle, vsig
        integer, intent(in) :: ie, ipg
        real*8, dimension(:), intent(in) :: vprel 

        integer :: iloi

!********************************************************!
        
        if (vprel(1)==11 .and. ietat(ie)==1) then              
            vnle=0.0d0*vnle
            vsig=0.0d0*vsig
        end if
               
end subroutine acier_rupture

!********************************************************!

end module acier

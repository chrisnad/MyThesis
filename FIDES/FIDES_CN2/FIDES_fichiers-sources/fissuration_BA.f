!-----------------------------------------------------------------------------------
! MODULE: fissuration_BA
!
!> @author C Nader
!
!> @brief
!> Ensemble des routines relative au(x) modËle(s) macroscopique(s) Orthotrope
!> de fissuration.
!>
!-----------------------------------------------------------------------------------

        module fissuration_BA


        contains


!-----------------------------------------------------------------------------------
!> @author C Nader
!> @brief Routine de calcul des contraintes pour modËles d'endommagement orthotrope.
!
!> @details
!> ### DESCRIPTION:
!> Routine de calcul des contraintes pour modËles d'endommagement
!>
!>----------------------------------------------------------------------------------
            subroutine fiss_sig_endo_ortho(vsig,vnle,vin,vh,vb,vprel,vdle,iloi,ie)

                use variables  , only : dime, ietat, inorm, iedng, check, ipas

                implicit none

                !--- Variables IN
                real*8, dimension(:,:), intent(in) :: vb
                real*8, dimension(:,:), intent(inout) :: vh
                real*8, dimension(:), intent(in) :: vdle, vprel
                integer, intent(in) :: iloi, ie

                !--- Variables OUT
                real*8, dimension(:), intent(out) :: vnle
                
                !--- Variables INOUT
                real*8, dimension(:), intent(inout) :: vsig, vin

                !--- Quantites globales
                real*8, dimension(size(vh,1),size(vh,2)) :: C
                integer :: id, plag

                !--- Quantites principales
                real*8, dimension(size(vb,1)) :: veps, vepsloc, vepslim

                ! Variables pour le modele d'endommagement
                real*8 :: Dini, vars, D, epsilon1
                real*8 :: f, E1, Xa, Xb, Xc, Ya, Yb, Yc
                real*8 :: ouv, Dnul1, Dnul2, Dnul
                
                !--- Variables JLT
                real*8, dimension(size(vh,1),size(vh,2)) ::  P !, Q
                real*8, dimension(dime,dime) :: V
                real*8, dimension(size(vh,1)) :: vsigm, vsigmloc
                logical :: vcrit, knull
                integer :: nc1
                real*8 :: var0, var
                real*8 :: F_var  !varcrit, varc, F_var, RT, W, E, nu
                real*8 :: WB, WAB, sigpa, Ea, Wel, Efa, eps0, epsc, epsl
                integer :: hyp
                real*8 :: D0, D0ini, Dk, RBA

!**********************************************************************************!

                !----- Calcul de l'ouverture de fissure initiale (pour non-interpenetration)
                ouv = 0.D0 ; Dnul1 = 1.d0 ; Dnul2 = 1.d0 ; D = 0. ; Dk = 1.d0
                !if (ietat(ie) == 3) ouv = fiss_ouv(ie) ! Calcul de l'ouverture fait selon la normale au pas precedent
                !if (ouv < 0.D0) Dnul = 0.d0

                select case(iloi)
                    case(25)  !----- Loi d'endommagement orthotrope trilineaire

                        !---- Initialisations
                        Dini = vin(1) ! endommagement initial
                        plag = int(vin(2)) ! partie 1, 2 ou 3 du chargement
                        vars = vin(3) ! seuil courant

                        !----- Etat de deformation
                        veps  = matmul(vb,vdle)

                        !----- Parametres sepcifiques ‡ la loi
                        id = 1
                        if (dime == 3) STOP '3D NON IMPLIMENTE pour fiss_sig_endo_ortho'
                        
                        E1 = vprel(id+3)
                        
                        Xa = vprel(id+7)  !\
                        Ya = vprel(id+8)  ! \
                        Xb = vprel(id+9)  !  \ Parametres
                        Yb = vprel(id+10) !  / de la loi tri-lineaire
                        Xc = vprel(id+11) ! /
                        Yc = vprel(id+12) !/

                        !----- Dans le repere du renforcement principal...
                        vepsloc  = veps

                        !----- Recuperation de la contrainte principale et de la deformation correspondante
                        epsilon1 = vepsloc(1)
                        if(plag == 0) then
                            if(epsilon1 .gt. Xa) then
                            
                                plag = 1
                                f    = (Yb - Ya)/(Xb - Xa)*(epsilon1 - Xa) + Ya
                                D    = 1 - f/(epsilon1*E1)
                                vars = epsilon1
                            
                            endif
                        elseif(plag == 1) then
                            if(epsilon1 .gt. Xb) then
                            
                                plag = 2
                                ietat(ie) = 3
                                inorm(ie,1) = 1
                                f    = (Yc - Yb)/(Xc - Xb)*(epsilon1 - Xb) + Yb
                                D    = 1 - f/(epsilon1*E1)
                                vars = epsilon1
                                
                            else
                            
                                f = (Yb - Ya)/(Xb - Xa)*(epsilon1 - Xa) + Ya
                                D = 1 - f/(epsilon1*E1)
                                
                                if(D .gt. Dini) then

                                    vars = epsilon1
                                
                                else
                                    
                                    D = Dini
                                        
                                endif
                            endif
                        elseif(plag==2) then
                            if(epsilon1 .gt. Xc) then
                                
                                plag = 3
                                inorm(ie,1) = 1
                                Dini = 1 - Yc/(Xc*E1)
                                
                            else
                                
                                f = (Yc - Yb)/(Xc - Xb)*(epsilon1 - Xb) + Yb
                                D = 1 - f/(epsilon1*E1)
                                
                                if (D .gt. Dini) then

                                    vars = epsilon1
                                else
                                    
                                    D = Dini
                                        
                                endif
                            endif
                        endif
                        
!                        C(1,:) = (/(1.d0-Dnul*D)*vh(1,1),(1.d0-Dnul*D)*vh(1,2),(1.d0-Dnul*D)*vh(1,3)/)
!                        C(2,:) = (/(1.d0-Dnul*D)*vh(1,2),              vh(2,2),              vh(2,3)/)
!                        C(3,:) = (/(1.d0-Dnul*D)*vh(1,3),              vh(2,3),              vh(3,3)/)
                        
!                        vsig = matmul(C,veps)
                        vsig = (1.d0 - D) * matmul(vh,veps)
                        
                        vnle(1:size(veps,1)) = (-Dnul*D*matmul(vb,vdle))

                        vin(1) = D
                        vin(2) = plag
                        vin(3) = vars

                    case(27)  !----- Loi d'endommagement isotrope beton + armatures

                        !---- Initialisations
                        Dini = vin(1) ! endommagement initial
                        D0ini= vin(2) ! seuil initial d'endommagement
                        vars = vin(3) ! seuil courant d'endommagement

                        !----- Etat de contrainte et de deformation
                        veps  = matmul(vb,vdle)
                        vsigm = matmul(vh,veps)
                        vepslim = 0.d0

                        !----- Parametres specifiques a la loi
                        id = 1
                        if (dime == 3) STOP '3D NON IMPLEMENTE pour loi 27'

                        E1 = vprel(id+3)

                        Rba   = vprel(id+7)          !    \
                        eps0  = Rba/(E1)
                        Ea    = vprel(id+8)          !     | Parametres de la loi
                        epsl  = eps0 + vprel(id+9)   !    /

                        sigpa = (epsl - eps0) * Ea   ! debut de non-linearite du comportement des aciers
                        epsc = -Ea*eps0/(E1-Ea)

                        !----- Recuperation du repere principal de ferraillage
                        V(:,1) = (/1., 0./)
                        V(:,2) = (/0., 1./)

                        !----- Calcul des matrices de changement de repere (Global->ferraillage)
                        P = fissBA_changement_repere(V,nc1,1)
                        !Q = fissBA_changement_repere(V,nc1,2)

                        !----- Dans le repere du ferraillage...
                        !vsigmloc = matmul(P,vsigm)
                        vepsloc = matmul(P,veps)

                        !----- Recuperation de la deformation correspondante dans l'axe principal de ferraillage
                        !sigma1 = vsigmloc(1)
                        epsilon1 = vepsloc(1)
                        !if (epsilon1<=0.d0) then
                        !    Dnul1=0.d0
                        !    Dnul2=0.d0
                        !end if

                        !----- Test pour activation du calcul de l'endommagement
                        vcrit = .false.
                        knull = .false.
                        if (ie == iedng) then
                        !if ((epsilon1) >= eps0 .and. ietat(ie)==0) then
                            !----- Changement "vcrit"
                            vcrit = .true.
                            !inorm(ie,:) = V(:,1)
                            !----- Calcul du seuil de deformation initial
                            !if(dime == 2) then
                            !    coef = 1.d0 / (vh(1,1)*vh(2,2) - vh(1,2)*vh(2,1))
                            !    var0 = coef * (vh(2,2)*sigma1 - vh(1,2)*vsigmloc(2)) !contraintes planes uniquement
                            !else if(dime == 3) then
                                !var0 = (sigma1 - nu*(vsigmloc(2)+vsigmloc(3)))/E
                            !    stop 'loi 26 non encore implantee en 3D'
                            !end if
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
                                !Dnul1 = 1.d0
                                !Dnul2 = 1.d0
                            end if
                        end if

                        !---- Initialisation de l'endommagement
                        D0 = D0ini
                        D  = Dini

                        !----- Declanchement et evolution de l'endommagement
                        if(vcrit .eqv. .true.) then

                            !---- Actualisation du seuil de rupture
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

                            !---- Loi d'evolution d'endommagement
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

!                        if((var>varcrit) .or. (D == .9999d0)) then
!                            D = .9999d0
!                            ietat(ie) = 1
!                            inorm(ie,:) = V(:,1) ! la valeur est déjà stockee dans inorm (cf. plus haut)
!                        end if

                        ! Contraintes et Rigidité
                        vh(1,1) = (1.d0-Dnul1*D0)*(1.d0-Dnul2*D)*vh(1,1) * Dk
                        vh(1,2) = (1.d0-Dnul1*D0)*(1.d0-Dnul2*D)*vh(1,2) * Dk
                        vh(2,1) = (1.d0-Dnul1*D0)*(1.d0-Dnul2*D)*vh(2,1) * Dk
                        vh(2,2) = (1.d0-Dnul1*D0)*vh(2,2) * Dk
                        vh(3,3) = (1.d0-Dnul1*D0)*vh(3,3) * Dk
                        vsig = matmul(vh,veps-vepslim)
                        !vnle(1:size(veps,1)) = -Dnul*D*veps

!                        vh = vh + vhpetit

                        vin(1) = D
                        vin(2) = D0
                        vin(3) = vars


                    case(26)  !----- Loi d'endommagement isotrope beton + armatures

                        !---- Initialisations
                        Dini = vin(1) ! endommagement initial
                        var0 = vin(2) ! seuil initial d'endommagement
                        vars = vin(3) ! seuil courant d'endommagement

                        !----- Etat de contrainte et de deformation
                        veps  = matmul(vb,vdle)
                        vsigm = matmul(vh,veps)

                        !----- Parametres specifiques a la loi
                        id = 1
                        if (dime == 3) STOP '3D NON IMPLEMENTE pour loi 26'

                        E1 = vprel(id+3)

                        WB    = vprel(id+7)  ! Wab     \
                        sigpa = vprel(id+8)  ! sigma_pl | Parametres de la loi
                        Ea    = vprel(id+9)  ! Ea      /

                        Wel   = .5d0*sigpa*sigpa/Ea              ! Energie de deformation elastique des aciers
                        WAB   = WB + Wel
                        Efa   = Ea                               ! Module des acier (homogeneise)
                        eps0  = 2.d0*Efa*(WAB-Wel)/(E1*(E1-Efa))
                        eps0  = dsqrt(eps0)                      ! Seuil d'endommagement initial (en deformation)
                        epsc  = 2.d0*E1*(WAB-Wel)/(Efa*(E1-Efa))
                        epsc  = dsqrt(epsc)                      ! fin du palier (en deformation)
                        epsl  = sigpa/Ea                         ! debut de non-linearite du comportement des aciers

                        !----- Recuperation du repere principal de ferraillage
                        V(:,1) = (/1., 0./)
                        V(:,2) = (/0., 1./)

                        !----- Calcul des matrices de changement de repere (Global->ferraillage)
                        P = fissBA_changement_repere(V,nc1,1)
                        !Q = fissBA_changement_repere(V,nc1,2)

                        !----- Dans le repere du ferraillage...
                        !vsigmloc = matmul(P,vsigm)
                        vepsloc = matmul(P,veps)

                        !----- Recuperation de la deformation correspondante dans l'axe principal de ferraillage
                        !sigma1 = vsigmloc(1)
                        epsilon1 = vepsloc(1)
                        if (epsilon1<=0.d0) Dnul=0.d0

                        !----- Test pour activation du calcul de l'endommagement
                        vcrit = .false.
                        !if (ie == iedng) then
                        if ((Dnul*epsilon1) >= eps0) then
                            !----- Changement "vcrit"
                            vcrit = .true.
                            inorm(ie,:) = V(:,1)
                            !----- Calcul du seuil de deformation initial
                            !if(dime == 2) then
                            !    coef = 1.d0 / (vh(1,1)*vh(2,2) - vh(1,2)*vh(2,1))
                            !    var0 = coef * (vh(2,2)*sigma1 - vh(1,2)*vsigmloc(2)) !contraintes planes uniquement
                            !else if(dime == 3) then
                                !var0 = (sigma1 - nu*(vsigmloc(2)+vsigmloc(3)))/E
                            !    stop 'loi 26 non encore implantee en 3D'
                            !end if
                            var0 = eps0
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

                        !---- HYPOTHESE DE MODELE
                        hyp = 3

                        !----- Declanchement et evolution de l'endommagement
                        if(vcrit .eqv. .true.) then

                            !---- Variable de pilotage de l'endommagement
                            var = Dnul * epsilon1

                            if(var .ge. vars) then

                                !---- Actualisation du seuil d'endommagement
                                vars = var

                                !---- Loi d'evolution d'endommagement
                                if (hyp==1 .or. hyp==2) then
                                    !--- avec palier
                                    if (eps0 < var .and. var <= epsc) then
                                        F_var = 0.d0
                                    elseif (epsc < var .and. var < epsl) then
                                        F_var = 1.d0 - (Efa/E1)*var/var0
                                    elseif (epsl < var) then
                                        F_var = 1.d0 - (sigpa/(E1*var0))
                                    end if
                                elseif (hyp==3) then
                                    !--- sans palier
                                    if ( var < epsl) then
                                        F_var = 1.d0 - (Efa/E1)*var/var0
                                    elseif (epsl < var) then
                                        F_var = 1.d0 - (sigpa/(E1*var0))
                                    end if
                                end if

                                !---- Calcul de l'endommagement
                                D = 1.d0 - var0/var * (1.d0 - F_var)

                            else

                                D = Dini

                            end if

                        else

                            D = Dini

                        end if

                        if (D < 0.d0) D = 0.d0
                        D = max(Dini,D)
                        D = min(D,1.d0)

!                        if((var>varcrit) .or. (D == .9999d0)) then
!                            D = .9999d0
!                            ietat(ie) = 1
!                            inorm(ie,:) = V(:,1) ! la valeur est déjà stockee dans inorm (cf. plus haut)
!                        end if

                        ! Contraintes et Rigidité
                        if (hyp==1) then
                            vh(1,1:2) = (1.d0-Dnul*D)*vh(1,1:2)
                            vh(2,1) = (1.d0-Dnul*D)*vh(2,1)
                            vsig = matmul(vh,veps)
                        elseif (hyp==2) then
                            vh = (1.d0-Dnul*D)*vh
                            vsig = matmul(vh,veps)
                            !vnle(1:size(veps,1)) = -Dnul*D*veps
                        elseif (hyp==3) then
                            vh(1,1) = (1.d0-Dnul*D)*vh(1,1)
                            vh(1,2) = 0.d0
                            vh(2,1) = 0.d0
                            vh(3,3) = (1.d0-Dnul*D)*vh(3,3)
                            vsig = matmul(vh,veps)
                            !vnle(1:size(veps,1)) = -Dnul*D*veps
                        end if

                        !vh = (1.d0-Dnul*D)*vh
!                        vh = vh + vhpetit

                        vin(1) = D
                        vin(2) = var0
                        vin(3) = vars


                    case default

                        stop 'fiss_sig_endo_ortho : cas non implante'

                end select

            end subroutine fiss_sig_endo_ortho



!---------------------------------------------------------------------------
!> @author JL Tailhan
!> @brief Changement de repère (module fissBA)
!
!> @details
!> ### DESCRIPTION:
!> Definition de la matrice de changement de repere :
!> du repere global au repere du ferraillage
!>
!>--------------------------------------------------------------------------
function fissBA_changement_repere_old(V,n,ind) result (M)

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

end function fissBA_changement_repere_old



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

function fissBA_changement_repere(V,n,iopt,typind) result (M)

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

end function fissBA_changement_repere





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
subroutine fissBA_distal()


    use variables, only : nelt, dime, ktypel, &
            & infele, vprelg, kprop, kprop0, idmax, &
            !& prob, pi, fiss, icompfiss, Dg, fc, &
            & pi, fissBA, Dg, fc, vitrav, &
            !--- gestion de l'alea (deb)
            & alea, enerjba, enerjba2
            !--- gestion de l'alea (fin)
    use lib_elem, only : elem_B
    use initialisation
    use proprietes_aleatoires
    use math

    implicit none

    real*8, dimension(:,:), allocatable :: vprelg0, ksig, vn, vb
    real*8, dimension(:), allocatable   :: vpg, PPBA1, PPBA2
    integer, dimension(:), allocatable  :: IGPPBA1, IGPPBA2
    real*8, dimension(1) :: PBA1, PBA2
    real*8, dimension(5) :: loi1, loi2

    integer :: iebe, i, ie, ie0, nbg, ipg, npg, id, iloi, dPlay, igrp, igmp, ine
    logical :: nrjbap, nrjbap2, iprobPBA1, iprobPBA2

!********************************************************!

    !----- Preparation de la renumerotation de kprop et vprelg  !
    ! Principe : on conserve la numerotation pour tous les      !
    ! groupes dont l'ancien numero est inferieur a nga. Puis on !
    ! "decale" d'un cran dans la numerotation tous les groupes  !
    ! qui ont un ancien numero plus grand que nga. Enfin, tous  !
    ! les elements du groupe nga sont renumerote en nga+i.      !

    if ((fissBA==1) .and. (alea)) then !---Pour le modèle de fissuration aleatoire

        nbg = maxval(kprop)
        call init_mat(vprelg0,(size(kprop,1)+nbg),size(vprelg,2))
        vprelg0(1:nbg,:) = vprelg(1:nbg,:)

        call init_vec(PPBA1,nelt)
        call init_vec(IGPPBA1,nelt)

        call init_vec(PPBA2,nelt)
        call init_vec(IGPPBA2,nelt)

        ie0 = 0
        iebe = 0

        !----- Indices pour trace de distributions de proprietes aleatoires
        iprobPBA1  = .false.
        iprobPBA2 = .false.

        !-- Boucle sur les elements
        do ie = 1, nelt

            !-- detection d'un element probabiliste (module et/ou resistance et/ou energie)
            nrjbap = .false.
            if (allocated(enerjba)) then
                if (count(enerjba%num==kprop0(ie))==1) nrjbap=.true.
            end if
            
            nrjbap2 = .false.
            if (allocated(enerjba2)) then
                if (count(enerjba2%num==kprop0(ie))==1) nrjbap2=.true.
            end if

            !-- si l'element appartient a un groupe probabiliste (energie)
            if (nrjbap .or. nrjbap2) then

                !-- on verifie la compatibilite de la loi de comp avec sa formulation probabiliste
                iloi=int(vprelg0(kprop(ie),1))
                !if (iloi==icompfiss) then

                    ie0 = ie0 + 1
                    iebe = iebe + 1

                    !-- Distributions aleatoires de proprietes mecaniques

                    !-- Energie plastique (ce test a ete conserve au cas ou...)
                    if (nrjbap) then
                    
                        ! on recupere le numero de groupe probabiliste correspondant a l'element
                        allocate(vitrav(1))
                        vitrav= find_vec(enerjba%num,kprop0(ie))
                        igrp = vitrav(1)
                        deallocate(vitrav)

                        loi1(1) = 0.       ! b
                        loi1(2) = 0.       ! c
                        loi1(3) = 0.       ! moy
                        loi1(4) = 0.       ! ect
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

                        ! preparation pour trace
                        iprobPBA1 = .true.

                        PPBA1(iebe) = PBA1(1) ! On enleve le Ve pour la visualisation de W
                        IGPPBA1(iebe) = enerjba(igrp)%num

                    end if !fin if nrjp

                    !-- Energie plastique 2
                    if (nrjbap2) then
                    
                        ! on recupere le numero de groupe probabiliste correspondant a l'element
                        allocate(vitrav(1))
                        vitrav= find_vec(enerjba2%num,kprop0(ie))
                        igrp = vitrav(1)
                        deallocate(vitrav)

                        loi2(1) = 0.       ! b
                        loi2(2) = 0.       ! c
                        loi2(3) = 0.       ! moy
                        loi2(4) = 0.       ! ect
                        loi2(5) = enerjba2(igrp)%loi        ! num

                        if (loi2(5)==1) then
                            ! Loi normale
                            loi2(3) = enerjba2(igrp)%param(1)
                            loi2(4) = enerjba2(igrp)%param(2)
                        elseif (loi2(5)==2) then
                            ! Loi de Weibull
                            loi1(1) = enerjba(igrp)%param(1)
                            loi1(2) = enerjba(igrp)%param(2)
                        elseif (loi2(5)==3) then
                            ! Loi log-normale
                            loi2(3) = enerjba2(igrp)%param(1)
                            loi2(4) = enerjba2(igrp)%param(2)
                        elseif (loi2(5)==4) then
                            ! Loi log-normale
                            loi2(1) = enerjba2(igrp)%param(1)
                            loi2(2) = enerjba2(igrp)%param(2)
                            loi2(3) = enerjba2(igrp)%param(3)
                        end if

                        !====> Tirage aleatoire de l'energie
                        PBA2 = distr_alea(loi2,1,ie0)

                        ! preparation pour trace
                        iprobPBA2 = .true.
                        
                        PPBA2(iebe) = PBA2(1) ! On enleve le Ve pour la visualisation de W
                        IGPPBA2(iebe) = enerjba2(igrp)%num

                    end if !fin if nrjp

                    !---- Stockage des valeurs
                    vprelg0(nbg + ie0,:) = vprelg(kprop(ie),:)
                    kprop(ie) = nbg + ie0

                    id = 1
                    if (dime==3) id = id-1

                    !----- Stockage de l'energie BA
                    if (nrjbap .or. nrjbap2) then

                        if (iloi==26.or.iloi==27) then  ! Loi orthotrope endommageable
                             vprelg0(kprop(ie),id+7) = PBA1(1)
                        end if

                    end if

                    !----- Stockage de l'energie BA2
                    if (nrjbap2) then

                        if (iloi==26.or.iloi==27) then  ! Loi orthotrope endommageable
                             vprelg0(kprop(ie),id+8) = PBA2(1)
                        end if

                    end if
                    
                !end if
            end if
        end do


        !-- Sorties :
        deallocate(vprelg)
        call init_mat(vprelg,ie0+nbg,size(vprelg,2))
        vprelg(1:ie0+nbg,:) = vprelg0(1:ie0+nbg,:)

        !-- Traces graphiques eventuels

        ! Visualisation de la distribution de l'energie
        if (iprobPBA1) then
            dPlay = 1
            do i=1,size(enerjba%num)
                loi1(5)=enerjba(i)%loi
                call alea_trac(PPBA1,loi1,200,dPlay,'energie BA',IGPPBA1,enerjba(i)%num)
            end do
        end if

        ! Visualisation de la distribution de l'energie
        if (iprobPBA2) then
            dPlay = 1
            do i=1,size(enerjba2%num)
                loi2(5)=enerjba2(i)%loi
                call alea_trac(PPBA2,loi2,200,dPlay,'energie BA',IGPPBA2,enerjba2(i)%num)
            end do
        end if

        deallocate(vprelg0,PPBA1,IGPPBA1,PPBA2,IGPPBA2)
    end if

end subroutine fissBA_distal
!-------------------------------------------------------!







!---------------------------------------------------------------------------
!> @author C Nader
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
            subroutine fissBA_distal_old()

                use variables,             only : vprelg, fissBA, kprop, nbrgrp, listcompfissBA
                use proprietes_aleatoires, only : ECDF
                use utilitaire,            only : lecture_var
                use initialisation
                
                implicit none
                
                real*8, dimension(:,:), allocatable :: vprelg0
                real*8, dimension(:,:), allocatable :: matvar
                real*8, dimension(:), allocatable :: XX, YY, randx, qdiss
                real*8 :: randy, rand1, rand2, qcalc, pq
                integer :: i, j, k, kk, l, m, n, nbg, high, low, clock, nbrgrp0
                integer, dimension(:), allocatable :: rseed
                character(len = 1) :: kstring
                logical :: criteria

                nbg = maxval(kprop)
                call init_mat(vprelg0,(size(kprop,1)+nbg),size(vprelg,2))
                vprelg0(1:nbg,:) = vprelg(1:nbg,:)
                
                call system_clock(COUNT=clock)
                call random_seed(size = n)
                allocate(rseed(n), source = clock + 37 * [(i, i = 0,n-1)])
                call random_seed(PUT = rseed)
                
                !listcompfissBA = (/25, 26/)  ! JLT : ligne a enlever car deja definie dans le module variable
                kk = 0
                nbrgrp0 = nbrgrp
                do i = 1, nbrgrp
                    if (any(listcompfissBA==vprelg(i,1))) then
                        m = size(vprelg,1) + 1
                        kk = kk+1
                        do j = 1, size(kprop,1)
                            if (kprop(j)==i) then
                                nbrgrp  = nbrgrp + 1
                                kprop(j)= nbrgrp
                            end if
                        end do
                        kstring = achar(48 + kk)

                        call lecture_var('Variables', matvar, kstring)
                        call init_vec(randx, size(matvar,2))
                        allocate(qdiss(size(matvar,1)))
                        
!|||||||||||||||||||||||||| PARTIE UNIQUE AU CALCUL BA UNIDIRECTIONNEL, MODELE TRI-LINEAIR ||||||||||||||||||||||||||||||!
                        qdiss(1:size(qdiss,1)) = matvar(1:size(qdiss,1),1)*matvar(1:size(qdiss,1),2) *0.5 &              !
                                              + (matvar(1:size(qdiss,1),3)-matvar(1:size(qdiss,1),1))     &              !
                                              * (matvar(1:size(qdiss,1),4)+matvar(1:size(qdiss,1),2))*0.5 &              !  Calcul de l'energie dissipee (l'aire sous la courbe)
                                              + (matvar(1:size(qdiss,1),5)-matvar(1:size(qdiss,1),3))     &              !
                                              * (matvar(1:size(qdiss,1),6)+matvar(1:size(qdiss,1),4))*0.5                !
!|||||||||||||||||||||||||| PARTIE UNIQUE AU CALCUL BA UNIDIRECTIONNEL, MODELE TRI-LINEAIR ||||||||||||||||||||||||||||||!
                  
                        do n = m, nbrgrp + m - 1 - nbrgrp0
                            criteria = .false.
                            do while (criteria .eqv. .false.)
                                do l = 1, size(randx)
                                    call ecdf(matvar(:,l), xx, yy, k)
                                    call random_number(randy)
                                    
                                    if (k < 2) then
                                        randx(l) = xx(1)
                                        cycle
                                    end if

                                    do while ((randy >= maxval(yy)) .or. (randy <= minval(yy)))
                                        call random_number(randy)
                                    end do
                                    
                                    do j = 1, k
                                        if(yy(j) <= randy) low = j
                                    end do
                                    
                                    if (low == k) then
                                        randx(l) = xx(low)
                                    else
                                        high = low+1
                                        randx(l) = xx(low) + (randy - yy(low))*(xx(high)-xx(low))/(yy(high)-yy(low))
                                    end if
                                    
                                end do
                                
                                vprelg0(n, :) = vprelg(i,:)
                                
!|||||||||||||||||||||||||| PARTIE UNIQUE AU CALCUL BA UNIDIRECTIONNEL, MODELE TRI-LINEAIR ||||||||||||||||||||||||||||||!
                                                                                                                         !
                                vprelg0(n, 8) = randx(1)                                                                 !
                                vprelg0(n, 9) = vprelg0(n,4)*randx(1)                                                    !
                                vprelg0(n,10) = randx(3)                                                                 !
                                vprelg0(n,11) = randx(4)                                                                 !
                                vprelg0(n,12) = randx(5)                                                                 !
                                vprelg0(n,13) = randx(6)                                                                 !
                                                                                                                         !
                                qcalc = vprelg0(n, 8)*vprelg0(n, 9) *0.5 &                                               !
                                     + (vprelg0(n,10)-vprelg0(n, 8))     &                                               !
                                     * (vprelg0(n,11)+vprelg0(n, 9))*0.5 &                                               !
                                     + (vprelg0(n,12)-vprelg0(n,10))     &                                               !
                                     * (vprelg0(n,13)+vprelg0(n,11))*0.5                                                 !
                                                                                                                         !
                                call random_number(rand1)                                                                !
                                call random_number(rand2)                                                                !
                                call ecdf(qdiss, xx, yy, k)                                                              !
                                if ((qcalc >= maxval(xx)) .or. (qcalc <= minval(xx))) cycle                              !
                                                                                                                         !
                                do j = 1, k                                                                              !
                                    if(xx(j) <= qcalc) low = j                                                           !
                                end do                                                                                   !
                                high = low + 1                                                                           !
                                                                                                                         !
                                pq = yy(low) + (qcalc - xx(low))*(yy(high)-yy(low))/(xx(high)-xx(low))                   !
                                                                                                                         !
                                if (((pq <= rand2) .and. (pq >= rand1)) .or. ((pq >= rand2) .and. (pq <= rand1))) then   !
                                    criteria = .true.                                                                    !                                                                              
                                end if                                                                                   !
                                                                                                                         !
!|||||||||||||||||||||||||| PARTIE UNIQUE AU CALCUL BA UNIDIRECTIONNEL, MODELE TRI-LINEAIR ||||||||||||||||||||||||||||||!
                                
                            end do                            
                        end do
                        
                        deallocate(vprelg)
                        call init_mat(vprelg,nbrgrp + m - 1 - nbrgrp0,size(vprelg,2))
                        vprelg(1:size(vprelg,1),:) = vprelg0(1:size(vprelg,1),:)   
                    
                    end if
                end do
                
                deallocate(vprelg0, randx, qdiss, matvar)
                
            end subroutine fissBA_distal_old



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
        subroutine fissBA_pilot(imetpilo,vsol,vduI,vduII,dlam,alphaba,MOTba,elefiss,iedg)

        !********************************************************!
        !  Pilotage du calcul : gestion du facteur de chargement !
        !********************************************************!

            use variables, only : dime, nelt, idmax, ktypel, nomtype, infele, &
                        & vprelg, kprop, ietat, iedngba, fissBA, pi,     &
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
            !real*8, dimension(dime), intent(out) :: normfissba
            real*8, intent(in) :: dlam
            real*8, intent(out) :: alphaba
            character(len=8), intent(inout) :: MOTba
            integer, intent(inout) :: elefiss

            ! Variables locales
            real*8, dimension(:), allocatable :: vdle, vdl0, vsi0, vdsi, vsi, vep, &
                                                & vsi0loc, vdsiloc, veploc, alpg, vdlI, vdlII
            real*8, dimension(:,:), allocatable :: vn, vb, vh, ksig, Ps, Pe
            
            real*8, dimension(dime, dime) :: rotmatba
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
                    if (fissBA==1) then

                        !- Initialisations
                        alph = 1.d0
                        iedngba = 0

                        !- Boucle sur les elements
                        do ie = binf, bsup

                            MOT = '  '
                            alpt = 1.d0
                            alpc = 1.d0
                            typel = nomtype(ktypel(ie))

                            !----- Pour les elements BA vierges
                            !if ((typel == 'MBT3' .or. typel == 'MTT4' .or. typel == 'MBT6' .or. typel == 'MBQ4')  &
                            if ((typel == 'MBQ4' .or. typel == 'MBT3')  &
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
                                Ps = fissBA_changement_repere(V,size(vb,1),1,1) !--- changement de repere pour les contraintes
                                Pe = fissBA_changement_repere(V,size(vb,1),2,1) !--- changement de repere pour les deformations

                                ! Dans le repere principal...
                                allocate(vdsiloc(size(vh,1)))
                                allocate(vsi0loc(size(vh,1)))
                                allocate(veploc(size(vh,1)))
                                vdsiloc = matmul(Ps,vdsi)
                                vsi0loc = matmul(Ps,vsi0)
                                veploc  = matmul(Pe,vep)

                                deallocate(vdl0,vdlI,vdlII,vdle,vep,vh,vb,vn)
                                deallocate(vsi0,vdsi,vsi,Ps,Pe)
                                
                                !----- Cas de la loi 27
                                if (iloi==27) then

                                    !----- Recuperation des parametres de la loi
                                    id = 1
                                    if (dime == 3) stop 'FIDES_fissBA_pilot : non programmee en 3D'

                                    RT = 1.001d0 * vprel(id+7)     ! Contrainte limite en traction pure
                                    RC = 1.001d0 * fc              ! Contrainte limite en cisaillement

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

                                elseif (iloi==1) then
                                
                                elseif (iloi==25) then

                                else
                                    stop 'FIDES_fissBA_pilot : loi non encore programmee'
                                end if

                                deallocate(vsi0loc,vdsiloc,veploc)

                                !----- On conserve le repere local sur le point de gauss considere
                                vbid2=reshape(V,(/1,dime*dime/))
                                vrtrav((ipg-1)*dime*dime+1:ipg*dime*dime)=vbid2(1,:)

                                !----- On conserve
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

                        alphaba = minval(alph)

                        if (alphaba < 1.d0) then
                            !----- On stocke le numero de l'element le plus dangereux
                            iedngba = find_num(alph,alphaba)
                            elefiss = count(alph<=0)
                            MOTba = rupt(iedngba)
                        else
                            alphaba = 1.d0
                        end if
                    end if

                case(2)
                    !----- Pilotage sur l'increment de contrainte avec recalcul de "dlam"

                case default
                    stop 'FIDES_fissBA_pilot : cas non encore programme'
            end select

        end subroutine fissBA_pilot





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
        subroutine fissBA_loi(iloi,icomp,ie)


            use variables, only : ietat !, icompfiss

            implicit none

            integer, intent(in) :: icomp, ie
            integer, intent(inout) :: iloi

!********************************************************!
            iloi = icomp

        end subroutine fissBA_loi








!-------------------------------------------------------!
!*******************************************************!
!       Initialisations liees au calcul fissure BA      !
!*******************************************************!
!-------------------------------------------------------!
            subroutine fissBA_init()


                use variables,      only : vprelg, fissBA, nbrgrp, nelt, dime, ietat, inorm, listcompfissBA
                use initialisation, only : init_mat, init_vec
            
                implicit none
            
                integer :: i

                fissBA=0

                do i = 1, nbrgrp
                    ! Cas particulier de la loi traduisant la fissuration du beton
                    if (any(listcompfissBA==vprelg(i,1))) then
                        fissBA=1
                        if (.not.allocated(ietat)) call init_vec(ietat,nelt)
                        if (.not.allocated(inorm)) call init_mat(inorm,nelt,dime)
                        exit
                    endif
                end do

            end subroutine fissBA_init

!********************************************************!
!********************************************************!
!********************************************************!
        end module fissuration_BA

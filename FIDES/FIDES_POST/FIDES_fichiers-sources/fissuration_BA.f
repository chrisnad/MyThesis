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

                use variables  , only : dime, ietat, inorm

                implicit none

                !--- Variables IN
                real*8, dimension(:,:), intent(in) :: vb
                real*8, dimension(:,:), intent(in) :: vh
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
                real*8, dimension(size(vb,1)) :: veps, vepsloc

                ! Variables pour le modele d'endommagement
                real*8 :: Dini, vars, D, epsilon1
                real*8 :: f, E1, Xa, Xb, Xc, Ya, Yb, Yc
                real*8 :: ouv, Dnul

!**********************************************************************************!

                !----- Calcul de l'ouverture de fissure initiale (pour non-interpenetration)
                ouv = 0.D0 ; Dnul = 1.d0 ; D = 0. ; plag = 0
                !if (ietat(ie) == 3) ouv = fiss_ouv(ie) ! Calcul de l'ouverture fait selon la normale au pas precedent
                !if (ouv < 0.D0) Dnul = 0.d0

                select case(iloi)
                    case(25)  !----- Loi d'endommagement orthotrope trilineaire

                        !---- Initialisations
                        Dini = vin(1) ! endommagement initial
                        plag = vin(2) ! partie 1, 2 ou 3 du chargement
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
                                
                                if(D .gt. Dini) then

                                    vars = epsilon1
                                else
                                    
                                    D = Dini
                                        
                                endif
                            endif
                        endif
                        
                        C(1,:) = (/(1.d0-Dnul*D)*vh(1,1),(1.d0-Dnul*D)*vh(1,2),(1.d0-Dnul*D)*vh(1,3)/)
                        C(2,:) = (/(1.d0-Dnul*D)*vh(1,2),              vh(2,2),              vh(2,3)/)
                        C(3,:) = (/(1.d0-Dnul*D)*vh(1,3),              vh(2,3),              vh(3,3)/)
                        
                        vsig = matmul(C,veps)
                        
                        vnle(1:size(veps,1)) = (-Dnul*D*matmul(vb,vdle))

                        vin(1) = D
                        vin(2) = plag
                        vin(3) = vars

                    case default

                        stop 'fiss_sig_endo_ortho : cas non implante'

                end select

            end subroutine fiss_sig_endo_ortho















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
            subroutine fissBA_distal()

                use variables,             only : vprelg, fissBA, kprop, kprop0, nbrgrp, listcompfissBA
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
                
                listcompfissBA = (/25, 26, 27/)
                kk = 0
                nbrgrp0 = nbrgrp
                do i = 1, nbrgrp
                    if (any(listcompfissBA==vprelg(i,1))) then
                        m = size(vprelg,1) + 1
                        kk = kk+1
                        do j = 1, size(kprop,1)
                            if (kprop(j)==i) then
                                nbrgrp  = nbrgrp + 1
                                kprop(j)= nbg + nbrgrp - nbrgrp0
                            end if
                        end do
                        kstring = achar(48 + kk)

                        call lecture_var('Variables', matvar, kstring)
                        call init_vec(randx, size(matvar,2))
                        allocate(qdiss(size(matvar,1)))
                        
!|||||||||||||||||||||||||| PARTIE UNIQUE AU CALCUL BA UNIDIRECTIONNEL, MODELE TRI-LINEAIR ||||||||||||||||||||||||||||||!
                        qdiss(1:size(qdiss,1)) = matvar(1:size(qdiss,1),1)*matvar(1:size(qdiss,1),2) *0.5 &              !
                                              + (matvar(1:size(qdiss,1),3)-matvar(1:size(qdiss,1),1))     &              !
                                              * (matvar(1:size(qdiss,1),4)+matvar(1:size(qdiss,1),2))*0.5 &              !  Calcul de l'energie dissipee (l'air sous la courbe)
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
              
            end subroutine fissBA_distal



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

            iloi = 25

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

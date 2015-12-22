!------------------------------------------------------------------------------
! MODULE: maj_contrainte
!
!> @author JL Tailhan
!
!> @brief
!> Ensemble des routines relatives au calcul et stockage des contraintes
!> apres convergence.
!>
!------------------------------------------------------------------------------
module maj_contraintes

contains

!-------------------------------------------------------------------------!
!            Subroutine globale de mise a jour des contraintes            !
!     (boucle sur les elements et appel de la routine de calcul           !
!                   des contraintes elementaires)                         !
!  Entree :                                                               !
!                                                                         !
!  Sorties :                                                              !
!-------------------------------------------------------------------------!
    subroutine maj_contGlob()

        use variables, only : nelt, vsol, kloce, vprelg, idmax, kprop
        use lib_elem,  only : elem_kloce2
        use initialisation, only : init_vec

        implicit none

        real*8, dimension(:), allocatable :: vdle    ! deplacements elementaires
        real*8, dimension(idmax) :: vprel            ! proprietes elementaires

        integer :: ie, ndle

        !----- Allocation/initialisation des deplacements totaux elementaires
        call init_vec(vdle,ndle)

        !--------------------- Boucle sur les elements ---------------------
        do ie = 1, nelt
            call elem_kloce2(ie,ndle)

            !----- Recuperation des deplacements totaux elementaires
            vdle = vsol(kloce(1:ndle))

            !----- Proprietes elementaires
            vprel = vprelg(kprop(ie),1:idmax)

            !----- Calcul des contraintes elementaires
            call maj_contElem(ie,vprel,vdle)

        end do

        !------------ Desallocation(S) --------------
        deallocate(vdle)


    end subroutine maj_contGlob

!-------------------------------------------------------------------------!
!        Subroutine elementaire de mise a jour des contraintes            !
!                                                                         !
!  Entree :                                                               !
!                                                                         !
!  Sorties :                                                              !
!-------------------------------------------------------------------------!
    subroutine maj_contElem(ie,vprel,vdle)

        use variables, only : varint, infele, ktypel, nomtype, calco, inict, &
                            & vcont, vnoli, vrint, kcont, knoli, krint, &
                            ! Pour gestion des elements libres
                            & elemlibre
        use initialisation, only : init_mat, init_vec
        use lib_elem, only : elem_B
        use aiguillage, only : rupture
        use contraintes
        use math

        implicit none

        ! Variables IN
        real*8, dimension(:), intent(in) :: vprel, vdle
        integer, intent(in) :: ie

        ! Variables locales
        real*8, dimension(:,:), allocatable :: ksig, vn, vb, vhep
        real*8, dimension(:,:), allocatable :: vsig0, vnle0, vsig, vnle, vint
        real*8, dimension(:),   allocatable :: vsi0, vnl0, vsi, vnl, vin

        real*8  :: detj, pdet, wpla, wplael
        integer :: nnel, ndln, npg, ipg
        integer :: nc1, nc2, nc3, idim1, idim2, idim3, iloc1, iloc2, iloc3, idm1, idm2, idm3

        character(len=5) :: typel


        !------- Recuperation des informations sur les elements ------
        nnel = infele(ktypel(ie))%nnel     ! nombre de noeuds par element
        ndln = infele(ktypel(ie))%ndln     ! nombre de ddl par noeud
        nc1  = infele(ktypel(ie))%ncgr     ! nombre de composantes du gradient
        nc2  = nc1 + 1                     ! nombre de composantes non lineaires
        nc3  = varint                      ! nombre de variables internes

        !------- Recuperation des points d'integration ----------
        allocate(ksig(size(infele(ktypel(ie))%Q,1),size(infele(ktypel(ie))%Q,2)))
        ksig = infele(ktypel(ie))%Q
        npg  = size(ksig,1)

        idim1 = npg*nc1
        idim2 = npg*nc2
        idim3 = npg*nc3

        iloc1 = kcont(ie)
        iloc2 = knoli(ie)
        iloc3 = krint(ie)
        idm1  = idim1/nc1
        idm2  = idim2/nc2
        idm3  = idim3/nc3

        call init_mat(vsig0, nc1, idm1)
        call init_mat(vnle0, nc2, idm2)
        call init_mat(vint, nc3, idm3)

        vsig0 = reshape(vcont(iloc1 : (iloc1 + idim1 - 1)),(/ nc1, idm1/))
        vnle0 = reshape(vnoli(iloc2 : (iloc2 + idim2 - 1)), (/nc2, idm2/))
        vint  = reshape(vrint(iloc3 : (iloc3 + idim3 - 1)), (/nc3, idm3/))

        !------------ Si calcul des contraintes aux noeuds -----------
        if (calco == 'NOEUD') then
            ! Obsolete
        else
            call init_mat(vsig, nc1, npg); call init_mat(vnle, nc2, npg)
        end if

        typel = nomtype(ktypel(ie))

        !----------- Initialisation pour calcul de l'energie plastique -------
        wplael = 0.d0

        !--- Si calcul des contraintes aux points de Gauss pour integration et
        !      calcul du residu elementaire et de la matrice de masse -------
        do ipg = 1, npg

            !-----  Calcul des fonctions d'interpolation et des derivees ------
            call elem_B(vn,vb,detj,ksig(ipg,:),ie)

            !---- Recuperation (aux noeuds) ou calcul (aux PG) des contraintes ---
            if(calco == 'NOEUD') then
                ! Obsolete
            else

                call init_vec(vsi0, size(vsig0,1)); call init_vec(vnl0, size(vnle0,1));
                call init_vec(vsi,  size(vb,1))   ; call init_vec(vnl,  size(vb,1));
                call init_vec(vin,  size(vint,1));

                vsi0 = vsig0(:,ipg)
                vnl0 = vnle0(:,ipg)
                vin  = vint(:,ipg)                

                !----- Calcul des contraintes en fonction du modele
                if (inict) vsi = vsi0
                call elem_sig(vsi,vnl,vin,ie,vhep,vb,vnl0,vprel,vdle,ipg)

                vsig(1:size(vsi),ipg) = vsi
                vnle(1:size(vnl),ipg) = vnl
                vint(1:size(vin),ipg) = vin

                wpla = vnl0(nc2) + .5*dot_product(vsi0+vsi,vnl(1:nc1)-vnl0(1:nc1))
                vnle(nc2,ipg) = wpla

                !------------ Gestion de l'element libre  --------
                if ((allocated(elemlibre)).and.(elemlibre(ie)==1)) then
                    !print*,'annulation de la contribution de l''element detache',ie
                    vsig(1:size(vsi),ipg) = 0.d0
                    vnle(1:size(vnl),ipg) = 0.d0
                    vint(1:size(vin),ipg) = 0.d0
                    wpla = vnl0(nc2)
                    vnle(nc2,ipg) = 0.D0
                end if

                deallocate(vsi0,vnl0)
            end if

            !------------ Integration sur l'element de l'energie plastique --------
            wplael = wplael + (pdet*wpla)

            !---------------------------- Nettoyage -------------------------------
            deallocate(vsi,vnl,vin, vn,vb,vhep)

        end do

        !--------------- Detection de la rupture d'un element -----------------
        call rupture(ie,vsig,vnle,wplael,vprel)

        !----------------- Stockage dans les vecteurs globaux --------------------

        !  Stockage dans vcont des nouvelles contraintes pour PG de l'element ie
        vcont(iloc1 : (iloc1+ idim1 - 1)) = reshape(vsig,(/idim1/))

        !  Stockage dans vnoli des nouvelles defo plast pour PG de l'element ie
        vnoli(iloc2 : (iloc2+ idim2 - 1)) = reshape(vnle,(/idim2/))

        !  Stockage dans vrint des nouvelles variables internes pour PG de l'element ie
        vrint(iloc3 : (iloc3+ idim3 - 1)) = reshape(vint,(/idim3/))


        !---------------------------- desallocation --------------------------------
        deallocate(ksig,vsig0,vnle0,vsig,vnle,vint)


    end subroutine maj_contElem


end module maj_contraintes

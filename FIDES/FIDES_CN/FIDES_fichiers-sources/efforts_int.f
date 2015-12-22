!------------------------------------------------------------------------------
! MODULE: efforts_int
!
!> @author JL Tailhan
!
!> @brief
!> Ensemble des routines relatives au calcul et assemblage
!> d'un effort interieur connaissant l'etat de contrainte.
!> Peut être utilise pour un calcul de residu ou prendre en compte des
!> des contraintes initiales, ...
!>
!------------------------------------------------------------------------------
module efforts_int

contains

!-------------------------------------------------------------------------!
!            Subroutine d'assemblage des vecteurs elementaires            !
!         d'efforts interieurs (connaissant l'etat de contraintes)        !
!     (boucle sur les elements et appel de la routine de calcul           !
!                   des efforts int. elementaires)                        !
!  Entree :                                                               !
!                                                                         !
!  Sorties :                                                              !
!-------------------------------------------------------------------------!
    subroutine assem_Fint(vfin)

        use variables, only : nelt, kloce, vprelg, idmax, kprop
        use lib_elem,  only : elem_kloce2
        use initialisation, only : init_vec
        use math

        implicit none

        real*8, dimension(:), intent(inout) :: vfin  ! Vecteur effort global
        real*8, dimension(:), allocatable :: vfe     ! Vecteur effort elementaire
        real*8, dimension(idmax) :: vprel            ! Proprietes elementaires

        integer :: ie, ndle


        !----- Boucle sur les elements
        do ie = 1, nelt

            !----- Vecteur de localisation elementaire (kloce=variable globale)
            call elem_kloce2(ie,ndle)

            !----- Recuperation des proprietes materielles elementaires
            vprel = vprelg(kprop(ie),1:idmax)

            !----- Calcul des contraintes elementaires
            allocate(vfe(ndle))
            call integ_contElem(ie,vprel,vfe)

            !----- Assemblage de vre
            vfin(kloce(1:ndle)) = vfin(kloce(1:ndle)) + vfe

            !----- Desallocation
            deallocate(vfe)

        end do

    end subroutine assem_Fint

!-------------------------------------------------------------------------!
!        Subroutine elementaire de mise a jour des contraintes            !
!                                                                         !
!  Entree :                                                               !
!                                                                         !
!  Sorties :                                                              !
!-------------------------------------------------------------------------!
    subroutine integ_contElem(ie,vprel,vfe)

        use variables, only : dime, infele, vcont, kcont, ktypel,&
                            ! Pour gestion des elements libres
                            & elemlibre
        use initialisation, only : init_mat, init_vec
        use lib_elem, only : elem_B
        use math

        implicit none

        ! Variables IN
        real*8, dimension(:), intent(in) :: vprel
        integer, intent(in) :: ie

        ! Variables OUT
        real*8, dimension(:), intent(inout) :: vfe

        ! Variables locales
        real*8, dimension(:,:), allocatable :: ksig, vn, vb
        real*8, dimension(:,:), allocatable :: vsig
        real*8, dimension(:), allocatable   :: vpg

        real*8  :: epais, detj, poids, pdet
        integer :: id, nnel, ndln, npg, ipg
        integer :: nc1, idim1, iloc1


        !----- Recuperation de l'epaisseur (dime=2 et CP)
        id = size(vprel)
        epais = 1.d0
        if ((dime == 2).and.(vprel(2)==3)) epais = vprel(id-2)  ! contraintes planes (CP)

        !----- Recuperation des informations sur les elements
        nnel = infele(ktypel(ie))%nnel     ! nombre de noeuds par element
        ndln = infele(ktypel(ie))%ndln     ! nombre de ddl par noeud
        nc1  = infele(ktypel(ie))%ncgr     ! nombre de composantes du gradient

        !----- Recuperation des points d'integration
        allocate(vpg(size(infele(ktypel(ie))%W)))
        allocate(ksig(size(infele(ktypel(ie))%Q,1),size(infele(ktypel(ie))%Q,2)))
        vpg  = infele(ktypel(ie))%W
        ksig = infele(ktypel(ie))%Q
        npg  = size(ksig,1)

        idim1 = npg*nc1
        iloc1 = kcont(ie)


        !----- Recuperation des contraintes dans vecteur global

        call init_mat(vsig, nc1, npg)
        vsig  = reshape(vcont(iloc1 : (iloc1 + idim1 - 1)),(/ nc1, npg/))
        
        !----- Gestion de l'element libre
        if ((allocated(elemlibre)).and.(elemlibre(ie)==1)) then
            vsig = 0.d0
        end if

        !----- Initialisation du vecteur element d'effort int
        vfe = 0.d0

        !-----  boucle sur les points de Gauss pour integration du
        !       vecteur elementaire
        do ipg = 1, npg

            !-----  Calcul des fonctions d'interpolation et des derivees
            call elem_B(vn,vb,detj,ksig(ipg,:),ie)
            poids = vpg(ipg)
            pdet = detj*poids*epais
            
            !----- Calcul de l'effort interieur elementaire
            vfe = vfe + pdet*matmul(transpose(vb),vsig(1:size(vb,1),ipg))

            !----- Nettoyage
            deallocate(vn,vb)

        end do
        
        !----- desallocation
        deallocate(ksig,vsig,vpg)


    end subroutine integ_contElem


end module efforts_int

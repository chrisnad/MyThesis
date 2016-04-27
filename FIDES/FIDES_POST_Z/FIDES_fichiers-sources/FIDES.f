!----------------------------------------------------------------------------------------------------------
! PROGRAMME: FIDES
! Un programme elements finis pour la fissuration des structures
!
!> @author JL Tailhan
!> @author J. Goncalvez (version 1.0 - 2010)
!
!> @brief
!> Programme principal FIDES.
!
!> @details
!> FIDES (FIssuration DEs Structures) est un programme entierement dedie a la modelisation
!> du comportement fissurant des structures en beton arme
!----------------------------------------------------------------------------------------------------------
program FIDES

    !**************************************************************************************
    !*********************************** MODULES ******************************************
    !**************************************************************************************

    !--------------------------------------------------------------------------------------
    !--- Modules generaux

    use variables
    use utilitaire
    use initialisation
    use math
    use sparse_matrix
    use lib_elem, only : elem_info
    use element_interface
    use acier, only: acier_init
    use fissuration
    use fissuration_BA
    use post_traitement
    use aiguillage
    use efforts_int
    use assemblage
    use mumps
    use Maillage
    use formatCSRC
    use Assemblage_CSRC
    use maj_contraintes

    use gnufor2

    implicit none

    !----------------- MPI (car MUMPS est ici compile en version parallele ----------------
    !include 'mpif.h'

    !**************************************************************************************
    !********************************** VARIABLES *****************************************
    !**************************************************************************************

    !--------------------------------------------------------------------------------------
    !--- Algorithme
    integer :: icalc, icalc_iter, i, j, ic
    real*8 :: time0, time, time_it, time_kg, time_rs, time1, time2     ! Temps cpu
    integer, dimension(8) :: values
    character(len = 50) :: nomfichier, nomlist

    !--------------------------------------------------------------------------------------
    !--- vecteurs 1er et 2nd membres globaux et incrementaux
    !real*8, dimension(:), allocatable :: vdu, vdut, vduI, vduII
    real*8, dimension(:), allocatable :: vdut, vduI, vduII
    real*8, dimension(:), allocatable :: vfg ! sollicitation globale
    real*8, dimension(:), allocatable :: daval, aval0
    real*8, dimension(:), allocatable :: vimpII, vresI, vresII,vrco

    !--------------------------------------------------------------------------------------
    !--- Matrice vkg au format sparse
    type(sparse) :: vkg
    real*8 :: p_nzmax, valll = 1. ! valeur a inserer pour conditions aux limites methode directe

    !--------------------------------------------------------------------------------------
    !--- Conditions aux limites
    integer :: ncond
    integer, dimension(:), allocatable :: condi
    real*8, dimension(:), allocatable :: valCond, v1    

    !--------------------------------------------------------------------------------------
    !--- Convergence
    real*8 :: normvcd ! Norme de vcond
    real*8 :: normvfg ! Norme de vfg
    real*8 :: lam, dlam, dlam0, dlmin
    real*8 :: ndu, ndf, denou, denof, realise ! norme relative du sol increm
    integer :: minpos
    logical :: nconv, bpas

    !--------------------------------------------------------------------------------------
    !--- Pilotage indirect en force avec controle des deplacements (npilot=2)
    integer, dimension(:), allocatable :: kdimp
    real*8 :: dlnum, dlden, normimp, dlam1
    integer :: ndimp, k

    !--------------------------------------------------------------------------------------
    !--- Pilotage indirect (npilot=3)
    real*8, dimension(:), allocatable :: vimpo
    integer :: idir,inoimp

    !--------------------------------------------------------------------------------------
    !--- Pour le pilotage du calcul
    integer :: imetpilo, noimp1, noimp2
    real*8  :: xut, xdut, xdu1, xdu2, ouvdech, coefdech
    logical :: decharg, arret

    !--------------------------------------------------------------------------------------
    !--- Variables
    integer :: compt


    !**************************************************************************************
    !*********************************** CALCULS ******************************************
    !**************************************************************************************

    !--------------------------------------------------------------------------------------
    !--- Initialisation de MPI (car MUMPS est ici compile en version parallele
    !call mpi_init(ierr)

    !--------------------------------------------------------------------------------------
    !--- Quelques initialisations
    nrep     = .false.
    indirect = .false.
    ievtp    = .false.
    nconv    = .false.

    !--------------------------------------------------------------------------------------
    !--- gestion de l'alea
    alea = .false.
    redu = .false.

    !--------------------------------------------------------------------------------------
    !--- auto contraintes
    inict = .false.

    !--------------------------------------------------------------------------------------
    !--- pilotage indrect
    decharg = .false.
    arret   = .false.
    coefdech = 1.d0

    compt = 0

    do while (.not.nconv)

    compt = compt +1


        print*,' Numero ', compt

        !---------------------------------------------------------------------------!
        !                         Initialisation variables                          !
        !---------------------------------------------------------------------------!

        call date_and_time(VALUES=values)

        ! Lecture des informations necessaires
        call lecture_info(nomfichier, nomlist, p_nzmax)

        imetpilo = 1                 ! Pour le pilotage du calcul des elements d'interface

        ! Indicateur du temps cpu
        call CPU_TIME(time0)

        !---------------------------------------------------------------------------!
        !                               Affichage                                   !
        !---------------------------------------------------------------------------!

        if(iprint>0) then
            print*; print*; print*;
            print*, '                ****************************************'
            print*, '                *      PROGRAMME FIDES NON LINEAIRE    *'
            print*, '                ****************************************'
            print*; print*; print*
        endif

        !---------------------------------------------------------------------------!
        !                           lecture des donnees                             !
        !---------------------------------------------------------------------------!

        print*, 'Le fichier choisi est : ',  dirdata//nomfichier
        call lecture(dirdata//nomfichier)

        !---------------------------------------------------------------------------!
        !                       Initialisation des rigidites                        !
        !---------------------------------------------------------------------------!

        if (CHOIX==1) then
            print*,'FIDES avec la version CSR-MUMPS et CL directe'
            print*,' '
            call init_sparse(vkg, ndlt, ndlt, per = p_nzmax)
        endif

        if (CHOIX==2) then
            print*,'FIDES avec la version CSRC-MUMPS Cl-directe'
            print*,' '

            ta_max = ndlt*100 !taille max de ja dans creation du fichier csrc
            print*,'creation du fichier csrc...  '
            call Maillage_CSRC(nomfichier,ta_max)
        endif

        !---------------------------------------------------------------------------!
        !                   Reduction homothetique du maillage                      !
        !                   Reduction en consequence du chargement                  !
        !---------------------------------------------------------------------------!

        !call PROB_reducmail()

        !---------------------------------------------------------------------------!
        !                 Initialisation et dimensionnement des tableaux            !
        !---------------------------------------------------------------------------!

        call initva(vcont, kcont, 'GRAD')
        call initva(vnoli, knoli, 'NOLI')
        call initva(vrint, krint, 'VINT')

        !---------------------------------------------------------------------------!
        !                   Initialisation pour algo non-lineaire                   !
        !---------------------------------------------------------------------------!

        call init_vec(detoscill,nelt)   ! compte le nombre de fois ou chaque eleme change d'etat,
                                        ! (doit etre remis a zero a chaque pas de calcul)
        call init_vec(vsol,ndlt)
        call init_vec(vdu, ndlt)
        call init_vec(v1,ndlt)
        call init_vec(vfg,ndlt)
        vfg = vfcg
        call init_vec(convdu,niter)
        call init_vec(convdf,niter)

        !---------------------------------------------------------------------------!
        !         Initialisations liees au calcul des elements d'interface,         !
        !      au calcul beton fissurant, au calcul des aciers, au calcul BA        !
        !---------------------------------------------------------------------------!

        call interf_init
        if (interf==1) call init_vec(elemlibre,nelt) ! liste des elements libres
        call fiss_init
        call acier_init
        call fissBA_init

        !---------------------------------------------------------------------------!
        !                  Partie de la distribution aleatoire                      !
        !---------------------------------------------------------------------------!

        if (ipost==0) then
            if (.not.nrep) then
                call CPU_TIME(time1)
                call distal
                call CPU_TIME(time2)
                print*, 'Temps de distribution : ', time2 - time1

                !----- Sauvegarde de la distribution
                !  iopt = 2 : sauvegarde de la distribution de proprietes meca

                call ecriture_resu(nomfichier,2)

                !----- Initialisation des auto-contraintes du beton
                if (inict) call inicont
            endif
        endif

        !---------------------------------------------------------------------------!
        !                         Reprise d'un calcul                               !
        ! optrep = 1 : reprise complete                                             !
        ! optrep = 2 : reprise de la distribution de proprietes mecaniques seulemt  !
        !---------------------------------------------------------------------------!

        if (nrep) then
            call lecture_resu(nomfichrep,optrep)
        endif

        !---------------------------------------------------------------------------!
        !             Calcul du vecteur second membre dans les cas :                !
        !                    - d'une reprise des calculs                            !
        !                    - de contraintes initiales                             !
        !---------------------------------------------------------------------------!

        call init_vec(vrco, ndlt)
        !call assem_Fint(vrco)

        !---------------------------------------------------------------------------!
        !           Initialisation pour prise en compte de l'evolution              !
        !                des proprietes materiau (algo implicite)                   !
        !---------------------------------------------------------------------------!

        if (ievtp) then
            call init_mat(vprelg0,size(vprelg,1), size(vprelg,2))
            call init_mat(vprelg1,size(vprelg,1), size(vprelg,2))
            call init_vec(tps,npas)
            vprelg0 = vprelg
        endif

        !---------------------------------------------------------------------------!
        !                     Creation des vecteurs unitaires de pilotage           !
        !---------------------------------------------------------------------------!

        call init_vec(daval,npas)
        call init_vec(aval0,npas+1)
        daval(1) = 0.d0
        daval(2:npas) = aval(2:npas) - aval(1:(npas-1))
        aval0 =(/0.d0, aval/)

        call init_vec(vimpII,size(vcond))

        !--- Pilotage direct en forces imposees (direct ou indirect)
        if ((npilot == 0).or.(npilot == 4)) then
            normvfg = norme(vfg)
            vfg = vfg/normvfg
            vimpII = vcond
            daval = daval*normvfg
            aval0 = aval0*normvfg

        !--- Pilotage direct en deplacements imposes (direct ou indirect)
        else if ((npilot == 1).or.(npilot == 3)) then
            normvcd = norme(vcond)
            vcond = vcond/normvcd
            vimpII = vcond
            daval = daval*normvcd
            aval0 = aval0*normvcd

        !--- Pilotage indirect (en deplacement impose vfg == vcond ou en force imposee vfg ~= vcond)
        else if (npilot == 2) then
            vfg = vcond
            normvfg = norme(vfg)
            vfg = vfg/normvfg
            vimpII = vcond
            daval = daval*normvfg
            aval0 = aval0*normvfg

        endif
        if ((npilot==3).or.(npilot==4)) npas=600

        !--- On calcule l'increment non nul minimal (en valeur absolue)
        minpos = minloc(abs(daval), 1, mask=abs(daval).gt.0)
        dlmin = daval(minpos)

        !----------------------------------------------------------------------------!
        !                 Affichage de l'entete avant calcul sur ecran               !
        !----------------------------------------------------------------------------!

        if (ipost==0) then
            if (iprint>0) then
                print*
                print*, '    ============= MODULE NLIN MATERIAU, FIDES  ============= '
                print*
                print*, '----- Mode de pilotage : Npilot 0 = SOLLIC, 1 = CL -----'
                print'(a15,1x,i5,4x,a18,i5)', 'Nombre de pas =', npas, 'Mode de pilotage = ', npilot
                print*
            endif
        endif

        !----------------------------------------------------------------------------!
        !                            Initialisations                                 !
        !                       pour boucle sur pas de temps                         !
        !----------------------------------------------------------------------------!

        ! Allocation et initialisation
        call init_vec(vdut, ndlt); call init_vec(vduI,ndlt); call init_vec(vduII,ndlt)
        call init_vec(vresI, ndlt); call init_vec(vresII, ndlt)

        if (interf==1) then
            allocate(osc_interf(nelt))
            do i=1,nelt
                allocate(osc_interf(i)%el(size(infele(ktypel(i))%Q,1)))
            enddo
        endif

        !----------------------------------------------------------------------------!
        !                      Initialisation du parametre de chargement             !
        !----------------------------------------------------------------------------!

        lam = 0.d0

        ipas = 1
        icalc = 1
        icalc_iter = 1
        bpas = .true.

        !---------------------------------------------------------------------------!
        !                       Post-traitement  (ipost = 1)                        !
        !---------------------------------------------------------------------------!

        if (ipost==1) then
            call FIGID(nomfichier,nomlist)
            bpas = .false.
            nconv = .true.
        endif

        !======================================================================================
        !============= DEBUT DE LA BOUCLE PRINCIPALE SUR LES PAS DE CHARGEMENTS ===============
        !======================================================================================

        do while (bpas)

            nconv = .false.

            !----- Affichage du pas de temps et du numero de calcul
            print'(a12,1x,i5,5x,a8,1x,i5,4x,a1,i5,a1)',' Pas_numero: ', ipas,  &
            &       'Calcul: ', icalc_iter ,'(',icalc,')'

            !----- mise a 0 du compteur d'oscillations
            detoscill = 0
            if (interf==1) then
                do i=1,nelt
                    osc_interf(i)%el=0
                enddo
            endif
             
            !----- Initialisation de l'increment de deplacement
            vdu = 0.d0

            !==================================================================================
            !===================== DEBUT DE LA BOUCLE SUR LES ITERATIONS ======================
            !==================================================================================

            do iter = 1, niter

                if (iter==1)  call CPU_TIME(time_it)

                !----------------- Modification des proprietes mecaniques ----------------
                if (ievtp) then
                    if (iter==1) then
                        vprelg1 = vprelg
                    !else
                        !--- chute de resistance pour fluage des elements d'interface
                        call interf_resichute(ipas) !,icalc_iter,iter)
                    endif
                endif

                !-------------------------------------------------------------------------!
                !------------------- Assemblage : construire vkg et vres -----------------!
                !-------------------------------------------------------------------------!

                call CPU_TIME(time_kg)
                !vresI = vrco + vfg*lam
                vresI = vfg*lam
                vresII = vfg

                !---------- Arret de prise en compte des contraintes initiales -----------
                if (.not.((iter==1).and.(ipas==1).and.(icalc==1))) inict = .false.

                !-------------------- Assemblage des vecteurs globaux --------------------
                if (CHOIX==1) then
                    call assem(vkg,vresI) !, ipas)
                endif
                if (CHOIX==2) then
                    call Assem_CSRC(supervkg,vresI)
                endif

                call CPU_TIME(time)

                time_kg = time - time_kg
                !print*,'temps reel d assemblage' ,time_kg

                !-------------------- Prise en compte d'une precontrainte -----------------
                if (precont) vresI = vresI+vfpg

                call CPU_TIME(time1)

                !-------------------------------------------------------------------------!
                !------------------ Introduction conditions aux limites ------------------!
                !-------------------------------------------------------------------------!

                if (CHOIX==1) then
                    ndimp = count(vcond/=0.)
                    ncond = count(kcond)
                    !--- Pilotage indirect : pour ne tenir compte que des blocages
                    if (npilot==2) ncond = count(kcond)-ndimp

                    call init_vec(condi, ncond);  call init_vec(valCond, ncond)
                    call init_vec(kdimp, ndimp)

                    j = 1
                    k = 1
                    v1 = 1.d0
                    do ic = 1,ndlt

                        if (kcond(ic)) then

                            !--- Pilotage direct en deplacements imposes
                            if ((npilot==1).or.(npilot==3)) then

                            v1(ic) = 0.0d0

                            do i = 1, ndlt
                                vresII(i) = vresII(i) - get_sparse(vkg,i,ic)*vimpII(ic)
                            enddo

                            vresII(ic) = vimpII(ic) !JLT On remplacera vimpII par 0.d0
                            vresI(ic) = 0.0d0

                            call del_row(vkg, ic) ! On efface la ligne ic
                            call del_col(vkg, ic) ! On efface la colonne ic

                            ! On determine a quel endroit on doit modifier la valeur de vkg
                            condi(j) = ic
                            j = j + 1

                            !--- Pilotages en forces imposees (direct ou indirect)
                            else if ((npilot==0).or.(npilot==2).or.(npilot==4)) then

                                if(vcond(ic)==0.d0) then ! On ne travaille que sur les blocages

                                    v1(ic) = 0.0d0
                                    vresII(ic) = 0.d0
                                    vresI(ic) = 0.0d0

                                    call del_row(vkg, ic) ! On efface la ligne ic
                                    call del_col(vkg, ic) ! On efface la colonne ic

                                    ! On determine a quel endroit on doit modifier la valeur de vkg
                                    condi(j) = ic
                                    j = j + 1

                                else ! On comptabilise (eventuellement) les depl imposes non nuls
                                    kdimp(k) = ic
                                    k = k+1

                                endif

                            endif
                        endif

                    enddo

                    valCond = 1.0d0
                    call insert_values(vkg,condi,condi,valCond)

                    deallocate(condi,valCond)
                endif ! CHOIX==1

                ! Methode directe avec le format CSRC
                if (CHOIX==2) then

                    ndimp = count(vcond/=0.)
                    ncond = count(kcond)
                    !--- Pilotage indirect : pour ne tenir compte que des blocages
                    if (npilot==2) ncond = count(kcond)-ndimp

                    call init_vec(condi, ncond);
                    call init_vec(kdimp, ndimp)

                    j = 1
                    k = 1
                    v1 = 1.d0
                    do ic = 1,ndlt

                        if (kcond(ic)) then

                            !--- Pilotage direct en deplacements imposes
                            if ((npilot==1).or.(npilot==3)) then

                                v1(ic) = 0.0d0

                                do i = 1, ndlt
                                    vresII(i) = vresII(i) - CSRC_val(supervkg,i,ic)*vimpII(ic)
                                enddo

                                vresII(ic) = vimpII(ic) !JLT On remplacera vimpII par 0.d0
                                vresI(ic) = 0.0d0

                                call CSRC_lig_col_zero(supervkg,ic)
                                call CSRC_set(supervkg,ic,ic,valll)
                                ! On determine a quel endroit on doit modifier la valeur de vkg
                                condi(j) = ic
                                j = j + 1

                            !--- Pilotages en forces imposees (direct ou indirect)
                            else if ((npilot==0).or.(npilot==2).or.(npilot==4)) then

                                if(vcond(ic)==0.d0) then ! On ne travaille que sur les blocages

                                    v1(ic) = 0.0d0
                                    vresII(ic) = 0.d0
                                    vresI(ic) = 0.0d0
                                    call CSRC_lig_col_zero(supervkg,ic)
                                    call CSRC_set(supervkg,ic,ic,valll)
                                    ! On determine a quel endroit on doit modifier la valeur de vkg
                                    condi(j) = ic
                                    j = j + 1

                                else ! On comptabilise (eventuellement) les depl imposes non nuls
                                    kdimp(k) = ic
                                    k = k+1

                                endif

                            endif
                        endif

                    enddo

                    deallocate(condi)

                endif ! fin choix = 2
                call CPU_TIME(time2)
                !print*,'temps dinsertion des CL',time2-time1

                !-------------------------------------------------------------------------!
                !---------------- Critere de convergence sur le residu -------------------!
                !-----------------   et test de convergence globale   --------------------!
                !-------------------------------------------------------------------------!

                if (iter == 2) then
                    !-- preparation du test de convergence
                    if (dlam0==0.d0) dlam0 = dlmin
                    !-- Residu a la premiere iteration (on ajoute dlam0*vresII pour tenir compte
                    !   de l'increment de chargement en deplacements imposes au debut du pas de tps)
                    denof = norme((vresI+dlam0*vresII)*v1) + 1.D-20
                    ndf = 1.d0

                elseif (iter > 2) then
                    !-- Critere de convergence sur le residu
                    ndf = norme((vresI)*v1)/denof

                    !-- Test de convergence globale
                    if ((ndu < 1.d-06) .and. (ndf < 1.d-04)) nconv = .true.

                    if (iter>50) then
                        nconv = .true.
                    endif

                endif

                !-------------------------------------------------------------------------!
                !----------------- SI CONVERGENCE SORTIE DU CALCUL -----------------------!
                !-------------------------------------------------------------------------!

                !*************************************************************************!

                !-- Calcul du temps par iteration
                if (iter > 1) then
                    call CPU_TIME(time)
                    time_it = time - time_it

                    !-- Impression ecran des infos sur convergence
                    if (iprint > 0) then
                        !-- Il faut imprimer iter-1 car décalage d'une iteration
                        if (CHOIX==1) print'(a8,i4,a14,e13.6,a8,e12.6,a13,e8.2,a7,e8.2,a6,e8.2)', 'CSR_Iter =', iter-1, &
                        &' - norme ndu = ', ndu,'ndf = ', ndf,' | CPU: Kg = ', time_kg, ' - Rs = ', time_rs, ' (It) ',time_it

                        if (CHOIX==2) print'(a8,i4,a14,e13.6,a8,e12.6,a13,e8.2,a7,e8.2,a6,e8.2)', 'CRC_Iter =', iter-1, &
                        &' - norme ndu = ', ndu,'ndf = ', ndf,' | CPU: Kg = ', time_kg, ' - Rs = ', time_rs, ' (It) ',time_it
                       
                        if ((iter > 50) .and. (nconv)) print*,'Attention: convergence forcee !'
                    endif
                    call CPU_TIME(time_it)
                endif

                !-- Desallocations et sortie si convergence
                if (nconv) then
                    deallocate(kdimp)
                    exit
                endif

                !*************************************************************************!

                !-------------------------------------------------------------------------!
                !------------------------------- Resolution ------------------------------!
                !-------------------------------------------------------------------------!

                call CPU_TIME(time_rs)

                if (CHOIX==1) then
                    call Solve_MUMPS(vkg, vresI, vresII, vduI, vduII)
                endif

                if (CHOIX==2) then
                    !call  MUMPS_CSRC222(supervkg, vresI, vresII, vduI, vduII)
                    call  MUMPS_CSRC(supervkg, vresI, vresII, vduI, vduII)
                endif

                call CPU_TIME(time)
                time_rs = time - time_rs

                !----------------- Calcul de l'increment de facteur de charge -------------

                !--- Pilotage direct en forces imposees
                if (npilot==0) then
                    if (iter == 1) then
                        vdut = 0.d0
                        dlam  = aval0(ipas) - lam + daval(ipas)
                        if (icalc_iter ==1) dlam0 = daval(ipas)
                        if (.not.((iter==1).and.(ipas==1).and.(icalc==1))) &
                        & call pilot(imetpilo, vsol, vduI, vduII, dlam)
                    else
                        iedng=0
                        dlam = 0.0d0
                        !call pilot(imetpilo, vsol, vduI, vduII, dlam)
                    endif

                !--- Pilotage direct en deplacements imposes
                elseif (npilot==1) then
                    if (iter == 1) then
                        vdut = 0.d0
                        dlam  = aval0(ipas) - lam + daval(ipas)
                        if (icalc_iter ==1) dlam0 = daval(ipas)
                        if (.not.((iter==1).and.(ipas==1).and.(icalc==1))) &
                        & call pilot(imetpilo, vsol, vduI, vduII, dlam)
                    else
                        iedng=0
                        dlam = 0.0d0
                        !call pilot(imetpilo, vsol, vduI, vduII, dlam)
                    endif

                !--- Pilotage indirect en forces imposees par controle du deplacement
                elseif (npilot==2) then
                    if (iter == 1) then
                        vdut = 0.d0
                        normimp = norme(vduII(kdimp))
                        dlam1 = aval0(ipas) - lam * normimp + daval(ipas)
                        if (icalc_iter ==1) dlam0 = daval(ipas)
                        dlnum = dot_product(vimpII(kdimp),dlam1*vimpII(kdimp)-vdut(kdimp)-vduI(kdimp))
                        dlden = dot_product(vimpII(kdimp),vduII(kdimp))
                        dlam  = dlnum/dlden

                        if (.not.((iter==1).and.(ipas==1).and.(icalc==1))) &
                        & call pilot(imetpilo, vsol, vduI, vduII, dlam)
                        dlam1 = dot_product(vimpII(kdimp),vdut(kdimp)+vduI(kdimp)+dlam*vduII(kdimp))
                    else
                        iedng=0
                        dlnum = dot_product(vimpII(kdimp),dlam1*vimpII(kdimp)-vdut(kdimp)-vduI(kdimp))
                        dlden = dot_product(vimpII(kdimp),vduII(kdimp))
                        dlam  = dlnum/dlden
                    endif

                !--- Pilotage indirect particulier (fendage)
                elseif ((npilot==3 .and. ipas .ge. 2).or.(npilot==4 .and. ipas .ge. 2)) then
                    deallocate(kdimp)
                    inoimp = 2 ! Nombre de noeuds avec condition imposee (2 pour disque 2D)
                    if(dime == 3) inoimp = 4
                    idir = 1
                    call init_vec(vimpo, inoimp)
                    vimpo=0.d0
                    call init_vec(kdimp, inoimp)

                    noimp1 = 9177
                    noimp2 = 196
                    ouvdech = .08D0

                    if(dime==2) then
                        !pilotage sur deux points (cas 2D)
                        !kdimp = (/dime*(165-1)+idir,dime*(207-1)+idir/)
                        kdimp = (/dime*(noimp1-1)+idir,dime*(noimp2-1)+idir/)
                        vimpo=(/-1.d0,1.d0/)
                    else if(dime == 3) then
                        !pilotage sur 4 points (cas 3D)
                        kdimp = (/dime*(1579-1)+idir,dime*(1814-1)+idir,dime*(1877-1)+idir,dime*(2087-1)+idir/)
                        vimpo=(/-1.d0,1.d0,-1.d0,1.d0/)
                    endif

                    xdut = dot_product(vimpo,vdut(kdimp))
                    xdu1 = dot_product(vimpo,vduI(kdimp))
                    xdu2 = dot_product(vimpo,vduII(kdimp))
                    xut  = dot_product(vimpo,vsol(kdimp))

                    if (iter == 1) then
                        dlam1 = 0.001d0
                        if ((.not.decharg).and.(xut > ouvdech)) decharg = .true.
                        if (decharg) then
                            dlam1    = .005d0
                            coefdech = -1.d0
                        endif
                        vdut = 0.d0
                        !dlnum = dot_product(vimpo,dlam1*vimpo-vdut(kdimp)-vduI(kdimp))
                        !dlden = dot_product(vimpo,vduII(kdimp))
                        dlnum = abs(dlam1)
                        dlden = abs(xdu2)
                        dlam  = coefdech * dlnum/dlden
                        if (icalc_iter ==1) dlam0 = dlam
                        call pilot(imetpilo, vsol, vduI, vduII, dlam)
                        !dlam1 = dot_product(vimpo,vdut(kdimp)+vduI(kdimp)+dlam*vduII(kdimp))
                        !dlam1 = dlam1/dot_product(vimpo,vimpo)
                        dlam1 = dlam*dlden
                    else
                        iedng=0
                        !dlnum = dot_product(vimpo,dlam1*vimpo-vdut(kdimp)-vduI(kdimp))
                        !dlden = dot_product(vimpo,vduII(kdimp))
                        dlnum = (-1.d0)*xdut*xdu1
                        dlden = xdu2*xdu2
                        dlam  = dlnum/dlden
                    endif

                    deallocate(vimpo)

                endif

                deallocate(kdimp)

                !------- Petit ajout dans le cas du pilotage indirect et d'un dechargement
                if ((npilot==3).and.(decharg)) then
                    if (lam*(lam+dlam) <= 0.d0) then 
                        dlam = -lam
                        arret = .true.
                    endif
                    if (abs(lam) < dlam1*1.d-10) dlam = 0.d0
                endif

                !------- Calcul de l'increment de deplacement et du deplacement total -----
                vdu = vduI + dlam*vduII
                vdut = vdut + vdu
                vsol = vsol + vdu
                lam = lam + dlam

                if (iter==1) dlamsauv = dlam

                !----- Traitement des refermetures (non interpenetration) ------
                call change_etat

                !-------------------------------------------------------------------------!
                !---------------- Test de convergence sur le deplacement -----------------!
                !-------------------------------------------------------------------------!

                if (iter == 1) then
                    !denou = norme((vduI + dlmin*vduII)*v1) + 1.D-20
                    if (dlam0==0.d0) dlam0 = dlmin
                    denou = norme((vduI + dlam0*vduII)*v1) + 1.D-20
                    ndu = 1.d0
                else
                    ndu = norme(vdu*v1)/denou
                endif

                convdu(iter) = ndu
                convdf(iter) = ndf

            enddo

            !==================================================================================
            !======================= FIN DE LA BOUCLE SUR LES ITERATIONS ======================
            !==================================================================================

            if (.not.nconv) exit

            !----- Affichage du pourcentage de pas de chargement realise
            if (npilot==2) then
                realise = abs(lam * normimp)/(abs(aval0(ipas)+daval(ipas))+1.d-15)
                !--- La ligne suivante et pour le cas particulier de dechargement a zero:
                if ((abs(lam)<1.d-10).and.(abs(aval0(ipas)+daval(ipas))<1.d-10)) realise = 1.d0
                if ((ipas==1).and.(abs(lam*normimp)-abs(aval0(ipas)+daval(ipas))<1.d-10)) realise = 1.d0
            elseif ((npilot==3).or.(npilot==4)) then
                realise = 1.
            else
                realise = abs(lam) / (abs(aval0(ipas)+daval(ipas))+1.d-15)
                !--- La ligne suivante et pour le cas particulier de dechargement a zero:
                if ((abs(lam)<1.d-10).and.(abs(aval0(ipas)+daval(ipas))<1.d-10)) realise = 1.d0
                if ((ipas==1).and.(abs(lam)-abs(aval0(ipas)+daval(ipas))<1.d-10)) realise = 1.d0
            endif

            realise = realise * 100
            if (iprint>0) then
                print'(a35,i3,a3,f6.2,a2,/)',' (Pourcentage realise pour le pas ', &
                &          ipas,' : ',realise,'%)'
            endif

            icalc = icalc + 1
            icalc_iter = icalc_iter+1

            !----- stockage des resultats propres au modele de fissuration ------
            call stock()

            !----- Mise a jour et stockage des contraintes apres convergence ------
            !call maj_contGlob()

            if (npilot==2) then                        !--- Pilotage indirect
                !if ((abs(aval0(ipas)+daval(ipas)-lam * normimp )<1.D-10).and.(nelcritic==0)) then
                if ((realise>99.995d0).and.(nelcritic==0)) then
                    ipas = ipas+1
                    icalc_iter = 1
                    !------------- Sauvegarde des resultats a chaque pas de temps ----------!

                    if (isauvcalc==0) then
                        call evolglob()
                        call ecriture(nomfichier,ipas-1)
                    endif

                endif
            else if ((npilot==3).or.(npilot==4)) then   !--- Pilotage indirect
                if (nelcritic==0) then
                    ipas = ipas+1
                    icalc_iter = 1
                    !------------- Sauvegarde des resultats a chaque pas de temps ----------!

                    if (isauvcalc==0) then
                        call evolglob()
                        call ecriture(nomfichier,ipas-1)
                    endif

                    !------------- Sauvegarde du dernier pas de calcul pour reprise --------!
                    call ecriture_resu(nomfichier)
                endif
            else                                       !--- Pilotage direct (force et/ou deplacement)
                !if ((abs(aval0(ipas)+daval(ipas)-lam)<1.D-10)) then
                if ((realise>99.995d0).and.(nelcritic==0)) then
                    ipas = ipas+1
                    icalc_iter = 1
                    !------------- Sauvegarde des resultats a chaque pas de temps ----------!

                    if (isauvcalc==0) then
                        call evolglob()
                        call ecriture(nomfichier,ipas-1)
                    endif

                    !------------- Sauvegarde du dernier pas de calcul pour reprise --------!
                    call ecriture_resu(nomfichier)
!               else
!                     call evolglob()
!                     call ecriture(nomfichier,ipas-1)


                endif
            endif

            !----- Incrementation du compteur de pas de calcul
            if (ipas == npas + 1) bpas = .false.

            !---------------- Sauvegarde des resultats a chaque pas de temps --------------
            if (isauvcalc==1 .and. nelcritic==0) then
                call evolglob()
                call ecriture(nomfichier,icalc-1)
            endif

            !------- pilotage indirect (npilot==3): arret du calcul si decharge finie -------!
            if ((npilot==3).and.(arret)) exit

        enddo

        !======================================================================================
        !============== FIN DE LA BOUCLE PRINCIPALE SUR LES PAS DE CHARGEMENTS ================
        !======================================================================================

        if (ipost /= 1) then

            if(iprint>0) then

                print*;

                if(iprint>1) then
                print*, '---------------- Impression de la solution ddl valeur ----------------'
                do i = 1, nnt
                    print*, 'noeud',i,'depl',vsol(infnoe(i)%dln)
                enddo
                endif

                call CPU_TIME(time)
                print'(2x,a10,f20.6,1x,a8)', 'Temps cpu = ', time - time0, 'secondes'
                print '(2x,a,i2,a,i2,a,i2,a,i2,a,i4)', 'Le calcul a ete lance a ',values(5),'h', &
                & values(6),', le ', values(3),'/',values(2),'/',values(1)
                print*;
                !print*, 'nombre oscillations', comptoscill
                !print*;

            endif

            call ecriture_evgb(nomfichier)

        endif


        !======================================================================================
        !================================== DESALLOCATIONS ====================================
        !======================================================================================


        !--------------------------------------------------------------------------------------
        !--- Desallocation de infele
        do i = 1, ntypel
            deallocate(infele(i)%Q,infele(i)%W,infele(i)%ksin,infele(i)%face)
        enddo

         deallocate(infele)

        !--------------------------------------------------------------------------------------
        !--- Desallocation de infnoe
        do i = 1, nnt
            deallocate(infnoe(i)%dln)
        enddo
        deallocate(infnoe)

        !--------------------------------------------------------------------------------------
        !--- Autres desallocations
        if (CHOIX==1) call free_sparse(vkg)

        if (CHOIX==2) then
            call CSRC_free(supervkg)
            do i=1,nnt
                deallocate(listelemautour(i)%el)
            enddo
            deallocate(listelemautour)
        endif

        deallocate(vsol,vdu,vfg,daval,vfcg,aval,aval0)
        deallocate(convdu,convdf)
        deallocate(vcont,vnoli,vrint,kcond,vcond,kcont,knoli)
        deallocate(vimpII,vresI,vresII,vdut,vduI,vduII)
        deallocate(vrco)
        deallocate(vcor,ktypel,kconec,nomtype)
        deallocate(vprelg,kprop,kprop0)
        deallocate(knddl)


        if (ievtp) deallocate(vprelg0)

        if (alea) then
            if (allocated(young)) deallocate(young)
            if (allocated(resist)) deallocate(resist)
            if (allocated(enerj)) deallocate(enerj)
            if (allocated(enerjpf)) deallocate(enerjpf)
            if (allocated(nrjbeto)) deallocate(nrjbeto)
            if (allocated(nrjfibr)) deallocate(nrjfibr)
        endif

        if (precont) then
            if (allocated(vfpg)) deallocate(vfpg)
        endif

        deallocate(v1)

        deallocate(detoscill)
        if (interf==1) then
            do i=1,nelt
                deallocate(osc_interf(i)%el)
            enddo
            deallocate(osc_interf)
        endif

        if (interf==1) deallocate(ietatpg,histetatpg1,histetatpg2,endo,elemlibre)
        if (fiss==1) deallocate(ietat,iouver1,iouver2,histetat1,histetat2,inorm,irloc)
        if (acierl) deallocate(limels,limrup,ecroui)

        !nconv = .true.
    enddo

    !------------------- Fin de MPI (car MUMPS est ici compile en version parallele -----------------------------!
    ! call MPI_FINALIZE(ierr)

end program FIDES

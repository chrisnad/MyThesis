!------------------------------------------------------------------------------
! MODULE: utilitaires
!
!> @author JL Tailhan
!
!> @brief
!> Utilitaires généraux: lecture, sauvegarde, initialisations, ...
!>
!------------------------------------------------------------------------------
module utilitaire


contains

!====================================================================================
!-----------------------------------------------------------------------------------!
!    Lecture du nom du fichier de donnees et des informations utilisees pour la     !
!                                suite du calcul                                    !
!  Entrees :                                                                        !
!  - nomfichier : Nom du fichier de donnees                                         !
!  - p_nzmax : pourcentage maximal de valeurs non nulles                            !
!-----------------------------------------------------------------------------------!

! Lecture du fichier lecture_info.data contenant les informations suivantes :
! 1. Nom du fichier de donnees
! 2. Lieu de calcul des contraintes
! 3. Indicateur de calcul explicite
! 4. Lieu de calcul des contraintes
! 5. nzmax pour matrices sparses

   subroutine lecture_info(nomfichier,nomlist,p_nzmax)

       use variables, only : iprint, explicite, calco, ipost, CHOIX

       implicit none

       character(len=50), intent(out) :: nomfichier, nomlist ! Nom du fichier de donnes lu
       real*8, intent(out) :: p_nzmax ! pourcentage d'elements non nuls (max)

       integer :: num_exp   ! pour explicite
       integer :: num_calco ! pour calco

       open(unit = 1, file = '../lecture_info.data',form = 'formatted', status = 'old', action = 'read')

       !--------- Nom du fichier de donnees ----------!
       read(1,*)
       read(1,*) nomfichier

       !-------- Lieu de calcul des contraintes ------!
       read(1,*)
       read(1,*) num_calco
       if(num_calco == 1)then; calco = 'GAUSS'
       else; calco = 'NOEUD'
       end if

       !------------- Code pour l'affichage ----------!
       read(1,*)
       read(1,*) iprint

       !------ Indicateur de calcul explicite --------!
       read(1,*)
       read(1,*) num_exp
       if(num_exp == 0)then; explicite = .false.
       else; explicite = .true.
       end if

       !------------ Pourcentage de nnz --------------!
       read(1,*)
       read(1,*) p_nzmax

       !-------- Option du post-traitement -----------!
       ipost = 0
       read(1,*)
       read(1,*) ipost

       !--------   Nom du fichier *.list   -----------!
       if (ipost==1) then
            nomlist = nomfichier(1:len_trim(nomfichier)-5)//'.list'
       end if

       !-------- Option generale -----------!
       CHOIX = 1
       read(1,*)
       read(1,*) CHOIX

       !------------ Fin de lecture ------------------!
     close(1)

   end subroutine lecture_info



!====================================================================================
!-----------------------------------------------------------------------------------!
!                      Lecture du fichier de donnees                            !
!  Entree :                                                                         !
!  - nomfichier : nom du fichier de donnes                                          !
!-----------------------------------------------------------------------------------!

! Recuperation de l'ensemble des donnees necessaires a un calcul.
! Ces donnees sont lues sur un fichiers formate du type de ceux
! generes lors des calculs CESAR.

   subroutine lecture (nomfichier)

      use variables
      use initialisation
      use lib_elem
      use math, only : find_mat
      implicit none

      integer :: i, j
      character(len=*), intent(in) :: nomfichier
      integer :: err, test
      integer :: imp, imp1 ! Caractere d'impression
      character(len=50) :: CHAINE
      logical :: testmod, new

      ! Coordonnees
      real*8, dimension(:), allocatable :: coor
      ! Nombre de groupe d'element
      integer :: ngrpe
      ! pour connaitre le nombre de noeuds par element
      integer, dimension(:), allocatable :: nnoelem
      ! Liste des numeros de noeuds pour tous les elements
      integer, dimension(:), allocatable :: numno
      ! Nombre de noeuds par element
      integer, dimension(:), allocatable :: nbnoeud
      integer :: indice1, indice2, nbmax
      ! Longueur de tableau
      integer :: lmax, lmax2
      ! Types d'elements
      character(len=5), dimension(:), allocatable :: types
      character(len=5), dimension(10) :: unitypes ! vecteur types sans repetition

      integer :: id, ileas, inoli, iloi, inolimax
      ! Epaisseur (cas dim 2)
      real*8, dimension(:), allocatable :: epais

      ! Forces de volume
      real*8 :: fvx, fvy, fvz
      integer :: ig

      ! Chargement et conditions aux limites
      integer :: N, dbloq, dimpo, nbg
      integer, dimension(:), allocatable :: nbloq, nimpo
      real*8, dimension(:), allocatable :: vimpo
      integer :: M, igen

      ! Informations sur la probabilisation
      integer :: nbprob, nog
      real*8, dimension(2) :: parbet

      ! Nombre de monitoring points
      integer :: nbnoe

      open(unit=1, file=nomfichier, form='formatted', status='old', action='read')

!----------------------- Determination du type de calcul ----------------------------

      err = 0
      do while (err/=-1)
          read(1, *, iostat = err), CHAINE
          if ((CHAINE == 'MCNL').and.(err/=-1)) then
              testmod = .true.
              if(iprint>0) then
                  print*
                  print*, 'Calcul incremental en mecanique (non lineaire)'
              end if
           end if
      end do

      rewind(1)

!----------------------- Lecture du maillage ----------------------------------------

      if (iprint>0) then
          print*
          print*, 'Construction du maillage en cours ...'
          print*
      end if

      err = 0
      do while (err /= -1)

         !----- Recuperation du maillage --------

         !----------------- Travail sur les noeuds du maillage -------------------
         read(1, *, iostat = err), CHAINE

         if((CHAINE == 'COOR').and.(err/=-1)) then
             ! Dimension du probleme et nombre de noeuds
             read(1, *), imp, imp1, nnt, dime
             ! Lecture des coordonnees
             call init_vec(coor,dime*nnt); call init_mat(vcor,dime,nnt)
             read(1, *), coor
             vcor=reshape(coor,(/dime,nnt/))
             deallocate(coor)
         end if

         !--------------- Travail sur les elements du maillage ---------------------

         if((CHAINE == 'ELEM').and.(err/=-1)) then
             ! Nombre d'elements et nombre de groupe d'elements
             read(1, *), imp, imp1, nelt, ngrpe
             nbrgrp = ngrpe
             lmax= nelt + 1
             call init_vec(nnoelem,lmax)
             read(1, *), nnoelem
             ! Longueur du tableau des numeros de noeuds lus dans le fichier *.data
             lmax2=nnoelem(lmax)-1

             ! Calcul du nombre de noeuds par elements et du max
             call init_vec(nbnoeud,nelt)
             !do i = 1, lmax-1 ! modif JLT 05/05/2014
             do i = 1, nelt
                 nbnoeud(i) = nnoelem(i+1) - nnoelem(i)
             end do
             nbmax=maxval(nbnoeud)

             ! Liste des numeros de noeuds pour tous les elements
             call init_vec(numno,lmax2); call init_mat(kconec,nelt,nbmax)
             read(1, *), numno

             !do i = 1, lmax-1 ! modif JLT 05/05/2014
             do i = 1, nelt
                indice1 = nnoelem(i)
                indice2 = nnoelem(i+1) - 1
                kconec(i,1:nbnoeud(i))=numno(indice1:indice2)
             end do

             deallocate(numno,nnoelem,nbnoeud)
             ! Lecture des types d'elements
             call init_vec(types,nelt,5); call init_vec(kprop,nelt); call init_vec(ktypel,nelt)

             read(1, *), types
             read(1, *), kprop

             nbg = maxval(kprop)
             call init_vec(kprop0,size(kprop))
             kprop0 = kprop

             ! Fonction unique de matlab + ktypel
             !-----------------------------
             unitypes = ''
             ntypel = 1
             unitypes(1)=types(1)
             ktypel(1) = 1

             do i = 2, nelt
                 new = .true.
                 do j = 1, ntypel
                    if (types(i) == unitypes (j)) then
                        ktypel(i)= j
                        new = .false.
                        exit
                    end if
                 end do

                 if(ntypel < 10 .and. new) then
                    ntypel = ntypel + 1
                    unitypes(ntypel) = types(i)
                    ktypel(i) = ntypel

                 else if(ntypel >= 10) then
                    print*, 'Attention le nombre de types est superieur a 10'
                    stop
                 end if

             end do

             !------------------------------

             call init_vec(nomtype,ntypel,5)
             nomtype = unitypes(1 : ntypel)

             deallocate(types);

             !----- Recuperation des parametres materiaux et options de calcul par groupe ------

             idmax = 0
             inolimax = 0
             call init_mat(vprelg,ngrpe, 20) ! 20 : valeur maximale du nombre de parametres associes a la loi + loi
             call init_vec(epais,ngrpe)

             do i = 1, ngrpe
                 read(1, '(a20)', iostat = err), CHAINE
                 if(iprint>0) then
                     print*
                     print'(a7, 1x, i2, a9, a20)', 'Groupe ', i, ' - Nom : ', CHAINE
                 end if
                 read(1, *), iloi
                 backspace(1)
                 vprelg(i,1) = iloi ! Numero de la loi de comportement

                 id = 1
                 ileas = 4
                 if (dime == 3) ileas = ileas - 1
                 if (iloi == 1 .and. testmod) then
                    if(iprint>0) print*, '   * Modele mecanique : elasticite isotrope'
                    inoli = 0
                    id = id + ileas + inoli
                    ! Lecture de (CP/DP), rho, E, nu
                 else if (iloi == 2 .and. testmod) then
                    ileas = 6
                    if (dime == 3) ileas = 10
                    if(iprint>0) print*, '   * Modele mecanique : elasticite orthotrope'
                    inoli = 0
                    id = id + ileas + inoli
                    ! 2D : Lecture de (CP/DP), rho, E1, E2, G12, nu12
                    ! 3D : Lecture de          rho, E1, E2, E3, G12, G23, G31, nu12, nu23, nu31
                 else if (iloi == 5 .and. testmod) then
                    if(iprint>0) print*, '   * Modele mecanique : elasticite isotrope fragile'
                    inoli = 2
                    id = id + ileas + inoli
                    ! Lecture de (CP/DP), rho, E, nu, RT, RC
                 else if (iloi == 10 .and. testmod) then
                    if(iprint>0) print*, '   * Modele mecanique : comportement de type elasto-endommageable (Mazars)'
                    inoli = 6
                    id = id + ileas + inoli
                    ! Lecture de (CP/DP),rho,E,nu,ktr0,At,Bt,Ac,Bc Beta
                 else if (iloi == 11 .and. testmod) then
                    if(iprint>0) print*, '   * Modele mecanique : comportement de type Von Mises'
                    inoli = 3
                    id = id + ileas + inoli
                    acierl = .true.
                    ! Lecture de (CP/DP), rho, E, nu, sig0
                 else if (iloi == 12 .and. testmod) then
                    if(iprint>0) print*, '   * Modele mecanique : comportement fissurant'
                    inoli = 2
                    id = id + ileas + inoli
                    ! Lecture de (CP/DP), rho, E, nu, sig0, wplmax
                 else if (iloi == 13 .and. testmod) then
                    if(iprint>0) print*, '   * Modele mecanique : comportement Coulomb'
                    inoli = 4
                    id = id + ileas + inoli
                    ! Lecture de (CP/DP), rho, E, nu, Rt, C, phi, psi
                 else if (iloi == 15 .and. testmod) then ! Ajouté par Giuseppe Rastiello, Aug 12
                    if(iprint>0) print*, '   * Modele mecanique : comportement fissurant endommageable'
                    inoli = 2
                    id = id + ileas + inoli
                    ! Lecture de (CP/DP), rho, E, nu, sig0, Wpf_max
                 else if (iloi == 25 .and. testmod) then ! CNader
                    if(iprint>0) print*, '   * Modele mecanique : modele BA orthotrope'
                    ileas = 6
                    inoli = 6
                    id = id + ileas + inoli
                    ! Lecture de (CP/DP), rho, E1, E2, G12, nu12, Xa, Ya, Xb, Yb, Xc, Yc
                 else if (iloi == 26 .and. testmod) then ! JLT
                    if(iprint>0) print*, '   * Modele mecanique : modele macro BA orthotrope endommageable'
                    ileas = 6
                    inoli = 3
                    id = id + ileas + inoli
                    ! Lecture de (CP/DP), rho, E1, E2, G12, nu12, Wab, sigm, Ea
                 else if (iloi == 27 .and. testmod) then ! JLT
                    if(iprint>0) print*, '   * Modele mecanique : modele macro BA orthotrope endommageable'
                    ileas = 6
                    inoli = 3
                    id = id + ileas + inoli
                    ! Lecture de (CP/DP), rho, E1, E2, G12, nu12, RAB, Ea, epsl
                 else if (iloi == 100 .and. testmod) then
                    if(iprint>0) print*, '   * Modele mecanique : contact Mohr-Coulomb'
                       inoli = 4
                       id = id + ileas + inoli
                       ! Lecture de (CP/DP), rho, E, nu, Rt, C, phi, psi
                 else if (iloi == 101 .and. testmod) then
                    if(iprint>0) print*, '   * Modele mecanique : contact Mohr-Coulomb-Tresca'
                    inoli = 4
                    id = id + ileas + inoli
                    ! Lecture de (CP/DP), rho, E, nu, Rt, C, phi, psi
                 else if (iloi == 102 .and. testmod) then
                    if(iprint>0) print*, '   * Modele mecanique : contact endo'
                    inoli = 6
                    id = id + ileas + inoli
                    ! Lecture de (CP/DP), rho, E, nu, Rt, C, phi, psi, depcrin, depcrit
                 else if ((iloi == 105) .and. testmod) then
                    if(iprint>0) print*, '   * Modele mecanique : contact endo + frott'
                    inoli = 5
                    id = id + ileas + inoli
                    ! Lecture de (CP/DP), rho, E, nu, Rt, C, phi, psi, energie
                 else if ((iloi == 106) .and. testmod) then
                    if(iprint>0) print*, '   * Modele mecanique : contact BRF'
                    inoli = 9
                    id = id + ileas + inoli
                    ! Lecture de (CP/DP), rho, E, nu, Rt, C, phi, psi, energie, xpic, xc, aplha, alphar
                 else if ((iloi == 107) .and. testmod) then
                    if(iprint>0) print*, '   * Modele mecanique : contact endo + frott'
                    inoli = 5
                    id = id + ileas + inoli
                    ! Lecture de (CP/DP), rho, E, nu, Rt, C, phi, psi, beta (beta=coef reducteur rigidite)
                 else if (iloi==108 .and. testmod) then
                    if(iprint>0) print*, '   * Modele mecanique : contact BRF + Energie fissuration beton'
                    inoli = 10
                    id = id + ileas + inoli
                    ! Lecture de (CP/DP), rho, E, nu, Rt, C, phi, psi, energie1, energie2, xp, xc, aplha, alphar
                 else if (iloi==158 .and. testmod) then
                    if(iprint>0) print*, '   * Modele mecanique : contact BRF'
                    inoli = 9
                    id = id + ileas + inoli
                    ! Lecture de (CP/DP), rho, E, nu, Rt, C, phi, psi, energie1,energie2, xpic2, xc2, alpha2
                 else
                    stop 'Attention !!!! Modele non encore implante !!!'
                 end if

                 read(1,*), vprelg(i,1:id)

                 if (id > 16) then
                     stop 'Attention !!! Le nombre de paramètre associes a la loi est superieur a 16'
                 end if

                 if (inoli > inolimax) inolimax = inoli
                 ! Recuperation de l'epaisseur
                 if (dime == 2) then
                    if (vprelg(i,2) == 1) then
                        if(iprint>0) print*, '   * Calcul en deformations planes'
                        epais(i)=1.D0
                    else if (vprelg(i,2) == 3)then
                        if(iprint>0) print*, '   * Calcul en contraintes planes'
                        read(1,*), epais(i)    ! Lecture de l'epaisseur
                    else
                        stop 'Attention !!! CP/DP: Modele non encore implante !!!'
                    end if

                   id = id + 1

                 end if

                 if (id > idmax) idmax = id

             end do

             if (dime == 2) vprelg(:,idmax) = epais

             deallocate(epais)

         end if
      end do

      rewind(1)

      !----------------- Quelques initialisations compte tenu du maillage ----------------

      !----------------------------------------------------------------------------------!
      !             creation de la structure de donnees relative aux elements            !
      !----------------------------------------------------------------------------------!
      call init_infele(ntypel)

      do i = 1, ntypel
          call elem_info (infele(i), nomtype(i))
      end do

      !--- Calcul du nombre total de ddl (et des ddl par noeud - vecteur knddl) ----------
      call ddlnode(ndlt)

      !------------ Initialisation des vecteurs contenues des conditions limites ---------
      call init_vec(vcond,ndlt); call init_vec(kcond,ndlt); call init_vec(vfcg,ndlt)
      vcond = 0.; kcond = .false.; vfcg = 0.

      rewind(1)

      !------------------ Lecture des chargements et conditions aux limites --------------
      err=0
      do while (err /= -1)
          read(1, *, iostat = err), CHAINE

          !--------------- Recuperation des conditions aux limites type blocage ----------

          if ((CHAINE == 'NUL').and.(err/=-1)) then
              read(1, *), imp
              test=1
              do while (test/=0)
                  read(1,*), N
                  call init_vec(nbloq,N)
                  read(1,*), nbloq
                  read(1,*), dbloq

                  do i = 1, N
                      if (dbloq > size(infnoe(nbloq(i))%dln,1)) then
                          print*, 'Erreur : Impossible de definir la condition limite type blocage au noeud ',nbloq(i)
                          stop
                      end if

                      kcond(infnoe(nbloq(i))%dln(dbloq)) = .true.
                      vcond(infnoe(nbloq(i))%dln(dbloq)) = 0
                  end do

                  deallocate(nbloq)

                  read(1,*, iostat = err), test

              end do
          end if

          !--------------------- Recuperation des deplacements imposes -------------------
          if ((CHAINE == 'IMP').and.(err/=-1)) then
             npilot = 1
             read(1,*), imp
             test=1
             j = 1
             do while(test/=0)
                 read(1,*), N

                 if(N==0) then
                    test=0
                 else
                    call init_vec(nimpo,N); call init_vec(vimpo,N)
                    read(1,*), nimpo
                    read(1,*), dimpo
                    read(1,*), vimpo

                    do i = 1, N
                       if (dimpo > size(infnoe(nimpo(i))%dln,1)) then
                          print*, 'Erreur : Impossible de definir les deplacements imposes au noeud ',nimpo(i)
                          stop
                       end if

                       kcond(infnoe(nimpo(i))%dln(dimpo)) = .true.
                       vcond(infnoe(nimpo(i))%dln(dimpo)) = vimpo(i)
                    end do
                    dirimpo(j) = dimpo
                    j = j + 1
                    deallocate(nimpo,vimpo)
                    read(1,*, iostat = err), test
                 end if
              end do
          end if

          !------------------- Recuperation des forces ponctuelles imposes ---------------
          if ((CHAINE == 'SOL').and.(err/=-1)) then
               npilot = 0
               read(1,*), M
               if(M == 0) then ! Cf. la doc de CESAR
                  test = 1
                  j = 1
                  do while (test/=0)
                     read(1,*), igen
                     write(2, *), igen
                     if(igen == 0) then
                        test=0
                     elseif(igen == 1) then
                        print*, 'Erreur : Methode de definition du chargement impossible'
                        stop
                     else if(igen == 2) then ! Cf. la doc de CESAR
                        read(1,*), N
                        call init_vec(nimpo,N); call init_vec(vimpo,N)
                        read(1,*), nimpo
                        read(1,*), dimpo
                        read(1,*), vimpo

                        do i = 1, N
                           if (dimpo > size(infnoe(nimpo(i))%dln,1)) then
                            print*, 'Erreur : Impossible de definir les chargements imposes au noeud ',nimpo(i)
                             stop
                           end if
                           vfcg(infnoe(nimpo(i))%dln(dimpo)) = vimpo(i)
                        end do
                        dirimpo(j) = dimpo
                        j = j +1
                        deallocate(nimpo,vimpo)
                     end if
                  end do
               else
                  print*, 'Erreur : Methode de definition du chargement impossible'
                  stop
               end if

          end if

      end do

      rewind(1)

      !-------------------------------- Lecture du module --------------------------------
      !- Recuperation des donnees du calcul, liste des pas de temps, precision, etc ...

      err = 0
      do while (err/=-1)
          read(1, *, iostat = err), CHAINE

          if ((CHAINE == 'MCNL').and.(err/=-1)) then
              ! MODELE.CALC, MODELE.KTAN : Matlab
              read(1,*), imp, npas, niter, ndmax, imet
              call init_vec(aval,npas)
              read(1,*,iostat = err), aval
          end if

      end do

      rewind(1)

      !------------------------------ Verification de reprise ----------------------------

      err = 0
      do while (err/=-1)
          read(1, *, iostat = err), CHAINE
          if ((CHAINE == 'INIT').and.(err/=-1)) then
              nrep = .true.
              read(1,*), optrep
              if ((optrep > 2) .or. (optrep<1)) stop  'Type de reprise non programme'
              read(1, *, iostat = err), CHAINE
              nomfichrep=trim(CHAINE)
              if (iprint>0) then
                  print*
                  print*, 'Reprise de calcul '
              end if
           end if
      end do

      rewind(1)

      !---------------------------------- Pilotage indirect ------------------------------

      err = 0
      do while (err/=-1)
          read(1, *, iostat = err), CHAINE
          if ((CHAINE == 'INDI').and.(err/=-1)) then
              indirect = .true.
              !--- On peut prevoir ici plusieurs methodes (LINE SEARCH, ARC LENGTH, etc. ...
              !npilot=2
              npilot=3  !Controle du déplacement pour pilotage en deplacement indirect
              !npilot=4  !Controle du deplacement pour pilotage en force indirect
              !---
              if(iprint>0) then
                  print*
                  print*, 'Pilotage indirect en force imposee par controle du deplacement (npilot=',npilot,')'
              end if
           end if
      end do

      rewind(1)

      !------------------------------------- Precontrainte -------------------------------

      err = 0
      do while (err/=-1)
          read(1, *, iostat = err), CHAINE
          write(2, *), CHAINE
          if ((CHAINE == 'PRCT').and.(err/=-1)) then
              call init_vec(vfpg,ndlt)
              vfpg = 0.
              precont = .true.
              !---
              if(iprint>0) then
                  print*
                  print*, 'Prise en compte d''un effort de precontrainte'
              end if
              test = 1
              !j = 1
              do while (test/=0)
                 read(1,*), igen
                 write(2, *), igen
                 if(igen == 0) then
                    test=0
                 elseif(igen == 1) then
                    stop 'Methode de definition du chargement impossible'
                 else if(igen == 2) then ! Cf. la doc de CESAR
                    read(1,*), N
                    write(2, *), N
                    call init_vec(nimpo,N); call init_vec(vimpo,N)
                    read(1,*), nimpo
                    write(2, *), nimpo
                    read(1,*), dimpo
                    write(2, *), dimpo
                    read(1,*), vimpo
                    write(2, *), vimpo
                    do i = 1, N
                        if (dimpo > size(infnoe(nimpo(i))%dln,1)) then
                           print*, 'Erreur : Impossible de definir les chargements imposes au noeud ',nimpo(i)
                           stop
                        end if
                        vfpg(infnoe(nimpo(i))%dln(dimpo)) = vimpo(i)
                    end do
                    !dirimpo(j) = dimpo
                    !j = j +1
                    deallocate(nimpo,vimpo)
                 end if
              end do
           end if
      end do

      rewind(1)
      !-------------------------------- Evolution temporelle -----------------------------

      err = 0
      do while (err/=-1)
          read(1, *, iostat = err), CHAINE
          if ((CHAINE == 'EVTP').and.(err/=-1)) then
              ievtp=.true.
              call init_vec(dtps,npas)
              read(1,*,iostat = err), dtps
              if (size(dtps,1)/=size(aval,1)) then
                  print*,'Lecture d''une lise de pas de temps : '
                  print*,'  donnee incompatible : la lise de pas de temps doit avoir'
                  print*,'  la meme longueur que la liste de chargements'
                  stop
              end if
              !---
              if(iprint>0) then
                  print*
                  print*, 'Lecture d''une liste de pas de temps pour probleme evolutif'
              end if
           end if
      end do

      rewind(1)

      !----------------------------- Controle de l'alea (deb) ----------------------------
      alea = .false.

      err = 0
      do while (err/=-1)
          read(1, *, iostat = err), CHAINE

          if ((CHAINE == 'PROB').and.(err/=-1)) then
              alea=.true.

              read(1, *), parbet(1:2)
              fc = parbet(1)
              Dg = parbet(2)

              !--- lecture du nombre de groupes probabilistes
              read(1, *), nbprob

              !--- boucle sur le nombre de groupes
              do i = 1, nbprob

                  !--- lecture du type de groupe (massif, fissure ou interface)
                  read(1, *, iostat = err), CHAINE

                  if (CHAINE(1:4) == 'modu') then
                      !--- lecture des informations pour la definition d'un module d'Young
                      !     probabiliste :
                      !     nog : nbre de gpes (massifs) concernes
                      !     num loi ipa : numero groupe - numero loi - nbre de param pr loi
                      !     param(1) param(2) ... : parametres de la loi de proba
                      read(1, *), nog
                      allocate(young(nog))

                      do j=1,nog
                          call init_vec(vitrav,3)
                          read(1, *), vitrav(1:3)
                          young(j)%num = vitrav(1)
                          young(j)%loi = vitrav(2)
                          young(j)%ipa = vitrav(3)
                          deallocate(vitrav)
                          if (young(j)%ipa >5) stop 'lecture :: pas plus de 5 parametres par loi probabiliste !'

                          if (young(j)%ipa >0) then
                              call init_vec(vrtrav,young(j)%ipa)
                              read(1, *), vrtrav
                              young(j)%param(1:young(j)%ipa) = vrtrav
                              deallocate(vrtrav)
                          end if
                      end do
                  end if

                  if (CHAINE(1:4) == 'resi') then

                      !--- lecture des informations pour la definition d'une resistance
                      !     probabiliste :
                      !     nog : nbre de gpes (massifs) concernes
                      !     num loi ipa : numero groupe - numero loi - nbre de param pr loi
                      !     param(1) param(2) ... : parametres de la loi de proba
                      read(1, *), nog
                      allocate(resist(nog))
                      do j=1,nog
                          call init_vec(vitrav,3)
                          read(1, *), vitrav(1:3)
                          resist(j)%num = vitrav(1)
                          resist(j)%loi = vitrav(2)
                          resist(j)%ipa = vitrav(3)
                          deallocate(vitrav)
                          if (resist(j)%ipa >5) stop 'lecture :: pas plus de 5 parametres par loi probabiliste !'

                          if (resist(j)%ipa >0) then
                              call init_vec(vrtrav,resist(j)%ipa)
                              read(1, *), vrtrav
                              resist(j)%param(1:resist(j)%ipa) = vrtrav
                              deallocate(vrtrav)
                          end if

                      end do
                  end if

                  if (CHAINE(1:4) == 'ener') then
                      !--- lecture des informations pour la definition d'une energie
                      !     probabiliste :
                      !     nog : nbre de gpes (massifs) concernes
                      !     num loi ipa : numero groupe - numero loi - nbre de param pr loi
                      !     param(1) param(2) ... : parametres de la loi de proba
                      read(1, *), nog
                      allocate(enerj(nog))
                      do j=1,nog
                          call init_vec(vitrav,3)
                          read(1, *), vitrav(1:3)
                          enerj(j)%num = vitrav(1)
                          enerj(j)%loi = vitrav(2)
                          enerj(j)%ipa = vitrav(3)
                          deallocate(vitrav)
                          if (enerj(j)%ipa >5) stop 'lecture :: pas plus de 5 parametres par loi probabiliste !'

                          if (enerj(j)%ipa >0) then
                              call init_vec(vrtrav,enerj(j)%ipa)
                              read(1, *), vrtrav
                              enerj(j)%param(1:enerj(j)%ipa) = vrtrav
                              deallocate(vrtrav)
                          end if

                      end do
                  end if

                  if (CHAINE(1:4) == 'enpf') then
                      !--- lecture des informations pour la definition d'une energie
                      !     probabiliste post fissuration:
                      !     nog : nbre de gpes (massifs) concernes
                      !     num loi ipa : numero groupe - numero loi - nbre de param pr loi
                      !     param(1) param(2) ... : parametres de la loi de proba
                      read(1, *), nog
                      allocate(enerjpf(nog))
                      do j=1,nog
                          call init_vec(vitrav,3)
                          read(1, *), vitrav(1:3)
                          enerjpf(j)%num = vitrav(1)
                          enerjpf(j)%loi = vitrav(2)
                          enerjpf(j)%ipa = vitrav(3)
                          deallocate(vitrav)
                          if (enerjpf(j)%ipa >5) stop 'lecture :: pas plus de 5 parametres par loi probabiliste !'

                          if (enerjpf(j)%ipa >0) then
                              call init_vec(vrtrav,enerjpf(j)%ipa)
                              read(1, *), vrtrav
                              enerjpf(j)%param(1:enerjpf(j)%ipa) = vrtrav
                              deallocate(vrtrav)
                          end if

                      end do
                  end if


                  if (CHAINE(1:4) == 'pba1') then
                      !--- lecture des informations pour la definition d'une energie
                      !     probabiliste beton arme:
                      !     nog : nbre de gpes (massifs) concernes
                      !     num loi ipa : numero groupe - numero loi - nbre de param pr loi
                      !     param(1) param(2) ... : parametres de la loi de proba
                      read(1, *), nog
                      allocate(enerjba(nog))
                      do j=1,nog
                          call init_vec(vitrav,3)
                          read(1, *), vitrav(1:3)
                          enerjba(j)%num = vitrav(1)
                          enerjba(j)%loi = vitrav(2)
                          enerjba(j)%ipa = vitrav(3)
                          deallocate(vitrav)
                          if (enerjba(j)%ipa >5) stop 'lecture :: pas plus de 5 parametres par loi probabiliste !'

                          if (enerjba(j)%ipa >0) then
                              call init_vec(vrtrav,enerjba(j)%ipa)
                              read(1, *), vrtrav
                              enerjba(j)%param(1:enerjba(j)%ipa) = vrtrav
                              deallocate(vrtrav)
                          end if

                      end do
                  end if

                  if (CHAINE(1:4) == 'pba2') then
                      !--- lecture des informations pour la definition d'une energie
                      !     probabiliste beton arme:
                      !     nog : nbre de gpes (massifs) concernes
                      !     num loi ipa : numero groupe - numero loi - nbre de param pr loi
                      !     param(1) param(2) ... : parametres de la loi de proba
                      read(1, *), nog
                      allocate(enerjba2(nog))
                      do j=1,nog
                          call init_vec(vitrav,3)
                          read(1, *), vitrav(1:3)
                          enerjba2(j)%num = vitrav(1)
                          enerjba2(j)%loi = vitrav(2)
                          enerjba2(j)%ipa = vitrav(3)
                          deallocate(vitrav)
                          if (enerjba2(j)%ipa >5) stop 'lecture :: pas plus de 5 parametres par loi probabiliste !'

                          if (enerjba2(j)%ipa >0) then
                              call init_vec(vrtrav,enerjba2(j)%ipa)
                              read(1, *), vrtrav
                              enerjba2(j)%param(1:enerjba2(j)%ipa) = vrtrav
                              deallocate(vrtrav)
                          end if

                      end do
                  end if

                  if (CHAINE(1:4) == 'nrjb') then
                      !--- lecture des informations pour la definition d'une energie
                      !     de fissuration probabiliste du beton :
                      !     nog : nbre de gpes (interfaces) concernes
                      !     num loi ipa : numero groupe - numero loi - nbre de param pr loi
                      !     param(1) param(2) ... : parametres de la loi de proba
                      read(1, *), nog
                      allocate(nrjbeto(nog))
                      do j=1,nog
                          call init_vec(vitrav,3)
                          read(1, *), vitrav(1:3)
                          nrjbeto(j)%num = vitrav(1)
                          nrjbeto(j)%loi = vitrav(2)
                          nrjbeto(j)%ipa = vitrav(3)
                          deallocate(vitrav)
                          if (nrjbeto(j)%ipa >5) stop 'lecture :: pas plus de 5 parametres par loi probabiliste !'

                          if (nrjbeto(j)%ipa >0) then
                              call init_vec(vrtrav,nrjbeto(j)%ipa)
                              read(1, *), vrtrav
                              nrjbeto(j)%param(1:nrjbeto(j)%ipa) = vrtrav
                              deallocate(vrtrav)
                          end if
                      end do
                  end if

                  if (CHAINE(1:4) == 'nrjf') then
                      !--- lecture des informations pour la definition d'une energie
                      !     de pontage (fibre/matrice) probabiliste :
                      !     nog : nbre de gpes (interfaces) concernes
                      !     num loi ipa : numero groupe - numero loi - nbre de param pr loi
                      !     param(1) param(2) ... : parametres de la loi de proba
                      read(1, *), nog
                      allocate(nrjfibr(nog))
                      do j=1,nog
                          call init_vec(vitrav,3)
                          read(1, *), vitrav(1:3)
                          nrjfibr(j)%num = vitrav(1)
                          nrjfibr(j)%loi = vitrav(2)
                          nrjfibr(j)%ipa = vitrav(3)
                          deallocate(vitrav)
                          if (nrjfibr(j)%ipa >5) stop 'lecture :: pas plus de 5 parametres par loi probabiliste !'

                          if (nrjfibr(j)%ipa >0) then
                              call init_vec(vrtrav,nrjfibr(j)%ipa)
                              read(1, *), vrtrav
                              nrjfibr(j)%param(1:nrjfibr(j)%ipa) = vrtrav
                              deallocate(vrtrav)
                          end if
                      end do
                  end if
              end do
          end if
      end do
      rewind(1)

      !---------------------------- Controle de l'alea (Fin) -----------------------------

      !----------------------------- Calcul de reactions (deb) ---------------------------
      ! ATTENTION CECI EST VOUE A DISPARAITRE
      reacn=.false.
      err = 0
      do while (err/=-1)
          read(1, *, iostat = err), CHAINE

          if ((CHAINE == 'REAC').and.(err/=-1)) then
              reacn = .true.
              read(1,*), diri
              read(1,*), noe1, noe2
          end if
      end do

      rewind(1)
      ! ATTENTION CECI EST VOUE A DISPARAITRE
      !----------------------------- Calcul de reactions (Fin) ---------------------------

      !--------------------------- Reduction du maillage (deb) ---------------------------
      ! ATTENTION CECI EST VOUE A DISPARAITRE

      err = 0
      do while (err/=-1)
          read(1, *, iostat = err), CHAINE

          if ((CHAINE == 'REDU').and.(err/=-1)) then
              redu=.true.
              ! Lecture du rapport de reduction : irap
              read(1,*), rap
          end if

      end do

      rewind(1)
      ! ATTENTION CECI EST VOUE A DISPARAITRE
      !--------------------------- Reduction du maillage (Fin) ---------------------------

      !---------------------- Auto-contraintes initiales (deb) ---------------------------

      err = 0
      do while (err/=-1)
          read(1, *, iostat = err), CHAINE

          if ((CHAINE == 'SIG0').and.(err/=-1)) then
              inict=.true.
              read(1,*), nog
              allocate(gpinicont(nog))
              read(1,*), gpinicont
          end if
      end do

      rewind(1)

      !---------------------- Auto-contraintes initiales (Fin) ---------------------------

      !-------------------- Definition des forces de volume (deb) ------------------------

      fvx = 0.D0 ; fvy = 0.D0 ; fvz = 0.D0

      err = 0
      do while (err/=-1)
          read(1, *, iostat = err), CHAINE

          if ((CHAINE == 'FVOL').and.(err/=-1)) then
              ! Lecture du nombre de groupes
              read(1, *), nog
              ! lecture pour chaque groupe de: num groupe, fvx, fvy (, fvz)
              call init_vec(vrtrav,dime)
              do j=1,nog
                  read(1, *), ig, vrtrav(1:dime)
                  vprelg(ig, idmax + 1) = vrtrav(1) !fvx
                  vprelg(ig, idmax + 2) = vrtrav(2) !fvy
                  if (dime==3) vprelg(ig, idmax + 3) = vrtrav(3) !fvz
              end do
              deallocate(vrtrav)
          end if
      end do

      rewind(1)

      if(dime == 3) then
          idmax = idmax + 3
      end if
      if(dime == 2) then
          idmax = idmax + 2
      end if
      !-------------------- Definition des forces de volume (fin) ------------------------

      !--------------- Lecture du deplacemend aux noeuds (monitoring points) --------------------
      nbnoe = 0
      err = 0
      do while (err/=-1)
          read(1, *, iostat = err), CHAINE
          if ((CHAINE == 'DEPL').and.(err/=-1)) then
              read(1,*), nbnoe                  ! lecture du nombre de noeuds (monitoring points)
              call init_vec(vnonoedp,nbnoe)     ! initialisation du vecteur de noeuds
              read(1,*), vnonoedp               ! lecture des numeros de noeud
              read(1,*), idirdp                 ! lecture de la direction du deplacement
           end if
      end do

      rewind(1)
      !------------- Lecture du deplacemend aux noeuds (monitoring points) ------------------------


      !------------------- Quelques affichages ----------------------------------------

      if(iprint>0) then
           print*
           print*, 'Dimension du probleme : ', dime
           print*, 'Nombre de groupes     : ', ngrpe
           print*, 'Nombre d''elements     : ', nelt
           print*, 'Nombre de noeuds      : ', nnt
           print*
      end if

      close(1, iostat = err)
      if (err == -1) stop 'Erreur dans la fermeture du fichier de donnees'

      close(2, iostat = err)
      if (err == -1) stop 'Erreur dans la fermeture du fichier : verif_fichier.data'

      ! Fin de lecture du fichier de donnees

  end subroutine lecture
  

!====================================================================================
!-----------------------------------------------------------------------------------!
!    Lecture du fichier contenant les variables pour le modele fiss_BA              !
!                                                                                   !
!  Entrees :                                                                        !
!  - filemane : Nom du fichier de donnees                                           !
!  - kstring  : Numero du fichier 'variables.txt' ex: Variables1.txt                !
!  Sortie  :                                                                        !
!  - matval   : Matrice contenant tous les valeurs                                  !
!-----------------------------------------------------------------------------------!


        subroutine lecture_var(filename, matval, kstring)
            
            use initialisation , only : init_mat

            implicit none

            character(len=*), intent(in) :: filename
            character(len=1), intent(in) :: kstring
            real*8, dimension(:,:), allocatable, intent(out) :: matval
            
            integer :: i, j, m, n
            character(len=3):: back
            character(len=5):: front
            character(len=50):: f
                
                back = '../'
                front = kstring // '.txt'
                f = back // filename // front

            open(unit = 1, file = f, form = 'formatted', status = 'old', action = 'read')
            
            read(1,*) m, n
            
            call init_mat(matval, m, n)
            
            do i=1,m
                read(1,*) matval(i,:)
            end do
           
            close(1)

        end subroutine lecture_var
   


!====================================================================================
subroutine ecriture_distribution(RR,iefi,EE,iebe)

      use variables, only : dirresu

      implicit none

      real*8, dimension(:), allocatable, intent(in) :: RR, EE
      integer, intent(in) :: iefi, iebe

      open(unit = 8, file = dirresu//'distribution.data', form = 'formatted', status = 'replace')


      write(8,'(a10)') 'resistance'
      write(8,*) iefi
      write(8,*) RR(1:iefi)

      write(8,'(a6)') 'module'
      write(8,*) iebe
      write(8,*) EE(1:iebe)

      close(8)

end subroutine ecriture_distribution

!------------------------------------------------------------------------------------!
!         Ecriture des résultats dans un fichier a chaque pas de calcul              !
!  Entrees :                                                                         !
!  - ipas : pas de calcul                                                            !
!------------------------------------------------------------------------------------!

    subroutine ecriture(nomfichier,ipas)

       use variables, only : dirresu, vsol, vcont, vnoli, ndlt, nelt, ietatpg, &
                     & dpglobal, foglobal, interf, &
                     & ietat, inorm, fiss, fissBA, vnonoedp, dime, idirdp, &
                     !----- option (ievtp)
                     & ievtp, tps

       use initialisation, only : init_vec

       implicit none

       character(len=*), intent(in) :: nomfichier
       integer, intent(in) :: ipas
       integer :: n, m
       character(len=len_trim(nomfichier)) :: nometude, nomsauve
       !real*8 :: deltads, deltadr, deltadf, deltadff, deltadfr ! essai fendage
       !real*8 :: unode ! jlt
       !real*8, dimension(size(vnonoedp,1)) :: vnoedp      ! vecteur deplacement au noeud des monitoring points
       real*8, dimension(:), allocatable :: vnoedp
       integer :: inoe    ! compteur pour les monitoring points

       n = size(vcont)
       m = size(vnoli)

       nomsauve = nomfichier(1:len_trim(nomfichier)-5)//'.list'
       nometude = nomfichier(1:len_trim(nomfichier)-5)//'.reac'

       if(ipas == 1) then
          open(unit = 2, file = dirresu//nomsauve, form = 'formatted', status = 'replace')
          open(unit = 3, file = dirresu//nometude, form = 'formatted', status = 'replace')

          write(2,'(a7)') 'Donnees'
          write(2,*) nelt, ndlt, n , m
          write(3,'(a16)') '#Solution_globale'
       else
          open(unit = 2, file = dirresu//nomsauve, form = 'formatted', status = 'old', position = 'append')
          open(unit = 3, file = dirresu//nometude, form = 'formatted', status = 'old', position = 'append')
       end if


         !pour l'essai brezilien
!         idir = 1
!         if(dime == 2) then
!              deltads = -vsol(dime*(165-1)+idir)+vsol(dime*(207-1)+idir)
!              deltadf = deltads
!              deltadr = deltads
!              !deltadff = deltads
!              !deltadfr = deltads
!         else if (dime ==3) then
!              deltadf = - vsol(dime*(1579-1)+idir) + vsol(dime*(1814-1)+idir)
!              deltadr = - vsol(dime*(1877-1)+idir) + vsol(dime*(2087-1)+idir)
!              deltads = 0.5*(deltadf + deltadr)
!              !deltadff = deltads + (deltadf-deltadr)*2.5d0/3.d0
!              !deltadff = deltads - (deltadf-deltadr)*2.5d0/3.d0
!         end if

         ! poutre Casanova
         !idir = 2
         !unode = vsol(dime*(4781-1)+idir) !poutre_casa1
         !unode = vsol(dime*(2267-1)+idir) !poutre2_casa1
         !unode = vsol(dime*(1175-1)+idir) !poutre2_casa3
         !unode = vsol(dime*(9648-1)+idir) !poutre2_casa4
         !unode = vsol(dime*(4869-1)+idir) !poutre2_casa4_transfo
         !unode = vsol(dime*(4852-1)+idir) !poutre3_casa1
         !unode = vsol(dime*(4473-1)+idir)  !PBAF_prob_transfo

         if (allocated(vnonoedp)) then
             call init_vec(vnoedp,size(vnonoedp,1))

             ! ecriture du deplacement, au noeud défini (monitoring point)
             do inoe = 1, size(vnonoedp,1)
                vnoedp(inoe) = vsol(dime*(vnonoedp(inoe)-1)+idirdp)
             end do

             if (ievtp) then
                 write(3,*) dpglobal, foglobal, vnoedp, tps(ipas)
             else
                 !-------Enregistrement des deplacements, forces, et monitoring points dans .reac
                 write(3,*) dpglobal, foglobal, vnoedp !deltads, deltadf, deltadr!, deltadff, deltadfr
             end if

         else
             if (ievtp) then
                 write(3,*) dpglobal, foglobal, tps(ipas)
             else
                 write(3,*) dpglobal, foglobal
             end if
         end if

         !-------Enregistrement des deplacements, forces, et monitoring points dans .reac
         !write(3,*) dpglobal, foglobal, vnoedp !deltads, deltadf, deltadr!, deltadff, deltadfr

         !-------- Enregistrement de la solution globale----------!
         write(2,'(a3, 1x, i5)') 'pas', ipas
         write(2,'(a16)') 'Solution_globale'
         write(2,*) dpglobal
         write(2,*) foglobal

         !----- Sauvegarde des tableaux generaux
         write(2,'(a8)') 'Solution'
         write(2,*) vsol

         write(2,'(a11)') 'contraintes'
         write(2,*) vcont

         !write(2,'(a12)') 'deformations'
         !write(2,*) vnoli

         !----- Stockage specifique aux elements d'interface (à modifier dans le futur)
         if (interf==1) then
            write(2,'(a6)') 'etatpg'
            write(2,*) ietatpg
         end if

         !----- Stockage specifique au calcul beton fissurant (à modifier dans le futur)
         if ((fiss==1).or.(fissBA==1)) then
            write(2,'(a4)') 'etat'
            write(2,*) ietat

            write(2,'(a4)') 'norm'
            write(2,*) inorm
         end if

       close(unit = 2)
       close(unit = 3)

    end subroutine ecriture

!------------------------------------------------------------------------------------!
!                Ecriture du fichier maillage contenant les fissures                 !
!                 pour le traitement graphique par le logiciel GID                   !
!                              PTS - 03/04/2011                                      !
!------------------------------------------------------------------------------------!

    subroutine ecriture_fiss(nomfichier,nompas,Mcm)

       use variables, only : ietat, inorm, ietatpg, dirresu, vcor, ktypel, kprop, &
           & kconec, nomtype ,dime, nelt, nnt, fiss, fissBA
       use initialisation, only : init_vec, init_mat
       use math, only : find_vec

       implicit none

       character(len=*), intent(in) :: nomfichier, nompas
       real*8, intent(in) :: Mcm
       character(len=5) :: typel
       character(len=15) :: nomfis
       character(len=len_trim(nomfichier)+len_trim(nompas)+7) :: nommshgid1
       integer :: i, j, ino, ie, iek, couc, icouc, neft1, neft2, neft3,&
           & ief, ief1, ief2, ief3, inm, inf1, inf2, inf3, np
       integer :: inoeu1, inoeu2, inoeu3, ibef2, ibaf2, imtt4, imbq4, imbt3, iejq4, &
           & iejq6, iejt6, ibar2, ibef3

       real*8, dimension(nelt) :: ouv
       real*8 :: Mc, maxouv, coul
       integer, dimension(:,:), allocatable :: kconecf1, kconecf2, kconecf3
       real*8, dimension(:,:), allocatable :: vcorf1, vcorf2, vcorf3

       real*8, dimension(:,:), allocatable :: coor_fiss

       integer, dimension(:), allocatable :: left1, left2, left3, left33

       type couche
           integer, dimension(:), allocatable :: elem
           real*8, dimension(3) :: coul
       end type couche
       type(couche), dimension(20) :: lcouc

       ! --------------------------------------------------------------------!
       nommshgid1 = nomfichier(1:len_trim(nomfichier)-5)//'_pas'//trim(nompas)//'.gid.msh'
       open(unit = 10, file = dirresu//nommshgid1, form = 'formatted', status = 'replace')
       !
       ! Determiner les fissures
       !
       neft1 = count(ietat==1) ! nombre de fissures ouvertes

       if (neft1 == 0) then
          call init_vec(left1,1)
          call init_mat(vcorf1,3,dime)
          call init_mat(kconecf1,3,dime)
       elseif (neft1 >= 1) then
          call init_vec(left1,neft1)
          call init_mat(vcorf1,2*neft1*dime,dime)
          call init_mat(kconecf1,2*neft1,dime)
          left1 = find_vec(ietat,1)
       end if

       neft2 = count(ietat==2) ! nombre de fissures refermees

       if (neft2 == 0) then
          call init_vec(left2,1)
          call init_mat(vcorf2,3,dime)
          call init_mat(kconecf2,3,dime)
       elseif (neft2 >= 1) then
          call init_vec(left2,neft2)
          call init_mat(vcorf2,2*neft2*dime,dime)
          call init_mat(kconecf2,2*neft2,dime)
          left2 = find_vec(ietat,2)
       end if

       neft3 = count(ietat==3) ! nombre de fissures (en cours d'endommagement)

       if (neft3 == 0) then
          call init_vec(left3,1)
          call init_mat(vcorf3,3,dime)
          call init_mat(kconecf3,3,dime)
       elseif (neft3 >= 1) then
          call init_vec(left3,neft3)
          call init_mat(vcorf3,2*neft3*dime,dime)
          call init_mat(kconecf3,2*neft3,dime)
          left33 = find_vec(ietat,3)
          ! on ne retient de ces elements que ceux qui ont une normale non nulle
          left3 = pack(left33,mask=sqrt(sum(inorm(left33,:)*inorm(left33,:),2))>0)
          neft3 = count(left3/=0)
       end if

       inf1 = 0; inf2 = 0 ; inf3 = 0
       inoeu1 = 0 ; inoeu2 = 0 ; inoeu3 = 0

       ! Calcul de l'ouverture de toutes les fissures
       ouv = fiss_ouv()
       maxouv = maxval(ouv)
       Mc = maxval( (/Mcm, maxouv/) )

       do i = 1, size(lcouc,1)
!JLT Verfier la dimension du tableau lcouc(i)%elem
           !call init_vec(lcouc(i)%elem,3*neft1)
           call init_vec(lcouc(i)%elem,(neft1+neft3))
           lcouc(i)%coul = (/ 1.d0 , (i-1.d0)/(size(lcouc,1)-1) , 0.d0 /)
       end do

       if ((neft1 > 0) .and. (Mc>0)) then
          ! Calcul de la couleur des elements ayant etat==1
          do ief1 = 1, neft1
             ie = left1(ief1)

             if (ouv(ie) > 0.d0) then

                call fisscoor(ie,coor_fiss,np)

                ! Coordonnees des points
                do i = 1, np
                   vcorf1(inoeu1 + i,:) = coor_fiss(i,:)
                end do

                ! Table de connectivite
                do j = dime, np
                   inf1 = inf1 + 1
                   kconecf1(inf1,1:dime) = (/(inoeu1+j-dime+i,i=1,dime,1)/)

                   ! Couleur
                   coul = 1.d0 - ouv(ie)/Mc

                   ! Definition de la couleur des couches (20 couches)
                   if ((coul>=0.0d0) .and. (coul<=0.05d0)) lcouc(1)%elem(inf1) = inf1
                   if ((coul>0.05d0) .and. (coul<=0.10d0)) lcouc(2)%elem(inf1) = inf1
                   if ((coul>0.10d0) .and. (coul<=0.15d0)) lcouc(3)%elem(inf1) = inf1
                   if ((coul>0.15d0) .and. (coul<=0.20d0)) lcouc(4)%elem(inf1) = inf1
                   if ((coul>0.20d0) .and. (coul<=0.25d0)) lcouc(5)%elem(inf1) = inf1
                   if ((coul>0.25d0) .and. (coul<=0.30d0)) lcouc(6)%elem(inf1) = inf1
                   if ((coul>0.30d0) .and. (coul<=0.35d0)) lcouc(7)%elem(inf1) = inf1
                   if ((coul>0.35d0) .and. (coul<=0.40d0)) lcouc(8)%elem(inf1) = inf1
                   if ((coul>0.40d0) .and. (coul<=0.45d0)) lcouc(9)%elem(inf1) = inf1
                   if ((coul>0.45d0) .and. (coul<=0.50d0)) lcouc(10)%elem(inf1) = inf1
                   if ((coul>0.50d0) .and. (coul<=0.55d0)) lcouc(11)%elem(inf1) = inf1
                   if ((coul>0.55d0) .and. (coul<=0.60d0)) lcouc(12)%elem(inf1) = inf1
                   if ((coul>0.60d0) .and. (coul<=0.65d0)) lcouc(13)%elem(inf1) = inf1
                   if ((coul>0.65d0) .and. (coul<=0.70d0)) lcouc(14)%elem(inf1) = inf1
                   if ((coul>0.70d0) .and. (coul<=0.75d0)) lcouc(15)%elem(inf1) = inf1
                   if ((coul>0.75d0) .and. (coul<=0.80d0)) lcouc(16)%elem(inf1) = inf1
                   if ((coul>0.80d0) .and. (coul<=0.85d0)) lcouc(17)%elem(inf1) = inf1
                   if ((coul>0.85d0) .and. (coul<=0.90d0)) lcouc(18)%elem(inf1) = inf1
                   if ((coul>0.90d0) .and. (coul<=0.95d0)) lcouc(19)%elem(inf1) = inf1
                   if ((coul>0.95d0) .and. (coul<=1.00d0)) lcouc(20)%elem(inf1) = inf1
                end do

                inoeu1 = inoeu1 + np

                deallocate(coor_fiss)
             end if
          end do
       end if

       if (neft2 > 0) then
          do ief2 = 1, neft2
             ie = left2(ief2)

             call fisscoor(ie,coor_fiss,np)

             ! Coordonnees des points
             do i = 1, np
                 vcorf2(inoeu2 + i,:) = coor_fiss(i,:)
             end do

             ! Tables de connectivite
             do j = dime, np
                 inf2 = inf2 + 1
                 kconecf2(inf2,1:dime) = (/(inoeu2+j-dime+i,i=1,dime,1)/)
             end do

             inoeu2 = inoeu2 + np

             deallocate(coor_fiss)
          end do
       end if

       if ((neft3 > 0) .and. (Mc>0)) then
          do ief3 = 1, neft3
             ie = left3(ief3)
       
             call fisscoor(ie,coor_fiss,np)
       
             ! Coordonnees des points
             do i = 1, np
                vcorf3(inoeu3 + i,:) = coor_fiss(i,:)
             end do

             ! Tables de connectivite
             do j = dime, np
                 inf3 = inf3 + 1
                 kconecf3(inf3,1:dime) = (/(inoeu3+j-dime+i,i=1,dime,1)/)
                 
                   ! Couleur
                   coul = 1.d0 - ouv(ie)/Mc

                   ! Definition de la couleur des couches (20 couches)
                   if ((coul>=0.0d0) .and. (coul<=0.05d0)) lcouc(1)%elem(inf3) = inf3
                   if ((coul>0.05d0) .and. (coul<=0.10d0)) lcouc(2)%elem(inf3) = inf3
                   if ((coul>0.10d0) .and. (coul<=0.15d0)) lcouc(3)%elem(inf3) = inf3
                   if ((coul>0.15d0) .and. (coul<=0.20d0)) lcouc(4)%elem(inf3) = inf3
                   if ((coul>0.20d0) .and. (coul<=0.25d0)) lcouc(5)%elem(inf3) = inf3
                   if ((coul>0.25d0) .and. (coul<=0.30d0)) lcouc(6)%elem(inf3) = inf3
                   if ((coul>0.30d0) .and. (coul<=0.35d0)) lcouc(7)%elem(inf3) = inf3
                   if ((coul>0.35d0) .and. (coul<=0.40d0)) lcouc(8)%elem(inf3) = inf3
                   if ((coul>0.40d0) .and. (coul<=0.45d0)) lcouc(9)%elem(inf3) = inf3
                   if ((coul>0.45d0) .and. (coul<=0.50d0)) lcouc(10)%elem(inf3) = inf3
                   if ((coul>0.50d0) .and. (coul<=0.55d0)) lcouc(11)%elem(inf3) = inf3
                   if ((coul>0.55d0) .and. (coul<=0.60d0)) lcouc(12)%elem(inf3) = inf3
                   if ((coul>0.60d0) .and. (coul<=0.65d0)) lcouc(13)%elem(inf3) = inf3
                   if ((coul>0.65d0) .and. (coul<=0.70d0)) lcouc(14)%elem(inf3) = inf3
                   if ((coul>0.70d0) .and. (coul<=0.75d0)) lcouc(15)%elem(inf3) = inf3
                   if ((coul>0.75d0) .and. (coul<=0.80d0)) lcouc(16)%elem(inf3) = inf3
                   if ((coul>0.80d0) .and. (coul<=0.85d0)) lcouc(17)%elem(inf3) = inf3
                   if ((coul>0.85d0) .and. (coul<=0.90d0)) lcouc(18)%elem(inf3) = inf3
                   if ((coul>0.90d0) .and. (coul<=0.95d0)) lcouc(19)%elem(inf3) = inf3
                   if ((coul>0.95d0) .and. (coul<=1.00d0)) lcouc(20)%elem(inf3) = inf3

             end do

             inoeu3 = inoeu3 + np

             deallocate(coor_fiss)
          end do
       end if

       ! End determiner les fissures

       imtt4 = 0;  imbq4 = 0;  iejq4 = 0; iejq6 = 0; iejt6 = 0; imbt3 = 0;
       ibar2 = 0 ; ibef2 = 0; ibef3 = 0; ibaf2 = 0;
       do inm = 1, maxval(ktypel)
          if (nomtype(inm)=='MTT4') imtt4 = 1
          if (nomtype(inm)=='MBQ4') imbq4 = 1
          if (nomtype(inm)=='MBT3') imbt3 = 1
          if (nomtype(inm)=='EJQ4') iejq4 = 1
          if (nomtype(inm)=='EJQ6') iejq6 = 1
          if (nomtype(inm)=='EJT6') iejt6 = 1
          if (nomtype(inm)=='BAR2') ibar2 = 1
          if (nomtype(inm)=='BEF2') ibef2 = 1
          if (nomtype(inm)=='BAF2') ibaf2 = 1
          if (nomtype(inm)=='BEF3') ibef3 = 1
       end do

       ! Recuperation toutes les coordonnees des noeuds du maillage
       write(10,*) 'MESH "Layer 0 " dimension 3 ElemType Triangle Nnode 3'
       write(10,*) 'Coordinates'
       do ino = 1, nnt
          if (dime == 2) write(10,*) ino, vcor(1:2,ino), 0.d0
          if (dime == 3) write(10,*) ino, vcor(1:3,ino)
       end do
       write(10,*) 'end coordinates'
       write(10,*) 'Elements'
       write(10,*) 'end elements'
       ! fin de tracer des coordonnees

       ! Initialisation des couches et des elements
       couc = 0
       iek = 0

       if (dime==3) then
         ! Tracer les elements tetraedres
         if (imtt4 == 1) then
            couc = couc + 1
            write(10,*) 'MESH "Layer', couc ,'" dimension 3 ElemType Tetrahedra Nnode 4'
            write(10,*) '# Color 0.000000 0.000000 1.000000'

            write(10,*) 'Elements'
            do ie = 1, nelt
                typel = nomtype(ktypel(ie))
                if (typel == 'MTT4') then
                  if (kprop(ie)==1) then
                    iek = iek + 1
                    write(10,*) iek, kconec(ie,1:4)
                  end if
                end if
            end do
            write(10,*) 'end elements'

            if (maxval(kprop)>1) then
              couc = couc + 1
              write(10,*) 'MESH "Layer', couc ,'" dimension 3 ElemType Tetrahedra Nnode 4'
              write(10,*) '# Color 0.000000 0.000000 1.000000'

              write(10,*) 'Elements'
              do ie = 1, nelt
                 typel = nomtype(ktypel(ie))
                 if (typel == 'MTT4') then
                   if (kprop(ie)==2) then
                     iek = iek + 1
                     write(10,*) iek, kconec(ie,1:4)
                   end if
                 end if
              end do
              write(10,*) 'end elements'
            end if
     end if
     ! fin des elements tetraedres

     ! Tracer les elements d'interface 3D
     if (iejt6 == 1) then
         if (count(ietatpg(:,1)==1)>0) then
             couc = couc + 1
         write(10,*) 'MESH "Layer', couc ,'ouvert " dimension 3 ElemType Prism Nnode 6'
         write(10,*) '# Color 0.000000 0.5000000 0.000000'

         write(10,*) 'Elements'
         do ie = 1, nelt
           typel = nomtype(ktypel(ie))
           if (typel == 'EJT6') then
             if (ietatpg(ie,1)==1) then
               iek = iek + 1
               write(10,*) iek, kconec(ie,1:6)
             end if
           end if
         end do
         write(10,*) 'end elements'
         end if

             if (count(ietatpg(:,1)==2)>0) then
         couc = couc + 1
         write(10,*) 'MESH "Layer', couc ,'cisaille " dimension 3 ElemType Prism Nnode 6'
         write(10,*) '# Color 0.000000 0.500000 0.500000'

         write(10,*) 'Elements'
         do ie = 1, nelt
           typel = nomtype(ktypel(ie))
           if (typel == 'EJT6') then
             if (ietatpg(ie,1)==2) then
              iek = iek + 1
              write(10,*) iek, kconec(ie,1:6)
             end if
           end if
         end do
         write(10,*) 'end elements'
         end if

     end if
     ! fin des elements d'interface 3D

     ! Tracer les elements lineaires
     if (ibef3 == 1) then
         couc = couc + 1
         write(10,*) 'MESH "Layer', couc ,'" dimension 3 ElemType Linear Nnode 2'
         write(10,*) '# Color 0.000000 0.000000 1.000000'
         write(10,*) 'Elements'
         do ie = 1, nelt
           typel = nomtype(ktypel(ie))
           if (typel == 'BEF3') then
             iek = iek + 1
             write(10,*) iek, kconec(ie,1:2)
           end if
         end do
         write(10,*) 'end elements'
     end if

       elseif (dime==2) then
     ! Tracer les elements quadrilateres
     if (imbq4 == 1) then
         couc = couc + 1
         write(10,*) 'MESH "Layer', couc ,'" dimension 3 ElemType Quadrilateral Nnode 4'
             write(10,*) '# Color 0.000000 0.000000 1.000000'

             write(10,*) 'Elements'
             do ie = 1, nelt
               typel = nomtype(ktypel(ie))
               if (typel == 'MBQ4') then
                 iek = iek + 1
                 write(10,*) iek, kconec(ie,1:4)
               end if
             end do
             write(10,*) 'end elements'
     end if
     ! fin des elements quadrilateres

     ! Tracer les elements d'interface
     if (iejq4 == 1) then
        if (count(ietatpg(:,1)==1)>0) then
          couc = couc + 1
          write(10,*) 'MESH "Layer', couc ,'ouvert " dimension 3 ElemType Quadrilateral Nnode 4'
          write(10,*) '# Color 0.000000 0.5000000 0.000000'

          write(10,*) 'Elements'
          do ie = 1, nelt
             typel = nomtype(ktypel(ie))
           if (typel == 'EJQ4') then
             if (ietatpg(ie,1)==1) then
              iek = iek + 1
              write(10,*) iek, kconec(ie,1:4)
             end if
           end if
         end do
         write(10,*) 'end elements'
       end if

           if (count(ietatpg(:,1)==2)>0) then
         couc = couc + 1
         write(10,*) 'MESH "Layer', couc ,'cisaille " dimension 3 ElemType Quadrilateral Nnode 4'
         write(10,*) '# Color 0.000000 0.500000 0.500000'

         write(10,*) 'Elements'
         do ie = 1, nelt
           typel = nomtype(ktypel(ie))
           if (typel == 'EJQ4') then
             if (ietatpg(ie,1)==2) then
              iek = iek + 1
              write(10,*) iek, kconec(ie,1:4)
             end if
           end if
         end do
         write(10,*) 'end elements'
       end if
     end if

     if (iejq6 == 1) then
       if (count(ietatpg(:,1)==1)>0) then
         couc = couc + 1
         write(10,*) 'MESH "Layer', couc ,'ouvert " dimension 3 ElemType Quadrilateral Nnode 6'
         write(10,*) '# Color 0.000000 0.000000 1.000000'

         write(10,*) 'Elements'
         do ie = 1, nelt
           typel = nomtype(ktypel(ie))
           if (typel == 'EJQ6') then
             iek = iek + 1
             write(10,*) iek, kconec(ie,1:6)
           end if
         end do
         write(10,*) 'end elements'
       endif

       if (count(ietatpg(:,1)==2)>0) then
         couc = couc + 1
         write(10,*) 'MESH "Layer', couc ,'cisaille " dimension 3 ElemType Quadrilateral Nnode 6'
         write(10,*) '# Color 0.000000 0.500000 0.500000'

         write(10,*) 'Elements'
         do ie = 1, nelt
           typel = nomtype(ktypel(ie))
           if (typel == 'EJQ6') then
             if (ietatpg(ie,1)==2) then
               iek = iek + 1
               write(10,*) iek, kconec(ie,1:6)
             end if
           end if
         end do
         write(10,*) 'end elements'
       end if
     end if

     ! fin des elements d'interface

     ! Tracer les elements triangulaires
     if (imbt3 == 1) then
         couc = couc + 1
         write(10,*) 'MESH "Layer', couc ,'" dimension 3 ElemType Triangle Nnode 3'
         write(10,*) '# Color 0.000000 0.000000 1.000000'
         write(10,*) 'Elements'

         do ie = 1, nelt
           typel = nomtype(ktypel(ie))
           if (typel == 'MBT3') then
             iek = iek + 1
             write(10,*) iek, kconec(ie,1:3)
           end if
         end do
         write(10,*) 'end elements'
         end if
     ! fin des elements triangulaires

     ! Tracer les elements lineaires
     if ((ibar2 == 1) .or. (ibef2 == 1) .or. (ibaf2 == 1)) then
         couc = couc + 1
         write(10,*) 'MESH "Layer', couc ,'" dimension 3 ElemType Linear Nnode 2'
         write(10,*) '# Color 0.000000 0.000000 1.000000'
         write(10,*) 'Elements'
         do ie = 1, nelt
           typel = nomtype(ktypel(ie))
           if ((typel == 'BAR2') .or. (typel == 'BEF2') .or. (typel == 'BAF2')) then
             iek = iek + 1
             write(10,*) iek, kconec(ie,1:2)
           end if
         end do
         write(10,*) 'end elements'
     end if
     ! fin des elements lineaire
       end if
       ! end tracer le maillage

       !-----------------------------------------------------------------!

       ! Tracer des fissures
       if (((fiss==1).or.(fissBA==1)) .and. (neft1+neft2+neft3>0) .and. (dime>1)) then

      ! Recuperation toutes les coordonnees des noeuds des fissures
          write(10,*) 'Coordinates'

          ! Tracer des noeuds des elements fissures cisailles
          do ino = 1, inoeu3
             if (dime==2) write(10,*) nnt+ino, vcorf3(ino,1:dime), 0.d0
             if (dime==3) write(10,*) nnt+ino, vcorf3(ino,1:dime)
          end do

          ! Tracer des noeuds des elements fissures refermees
          do ino = 1, inoeu2
             if (dime==2) write(10,*) nnt+ino+inoeu3, vcorf2(ino,1:dime), 0.d0
             if (dime==3) write(10,*) nnt+ino+inoeu3, vcorf2(ino,1:dime)
          end do

          ! Tracer des noeuds des elements fissures
          do ino = 1, inoeu1
             if (dime==2) write(10,*) nnt+ino+inoeu3+inoeu2, vcorf1(ino,1:dime), 0.d0
             if (dime==3) write(10,*) nnt+ino+inoeu3+inoeu2, vcorf1(ino,1:dime)
          end do

          write(10,*) 'end coordinates'
      ! fin de tracer des coordonnees

      if (dime==2) nomfis = 'Linear'
      if (dime==3) nomfis = 'Triangle'

      ! Tracer les elements fissures 2D (Line) ou 3D (Triangle)
!      if (inf3 > 0) then
!        ! Elements endommages-fissures (etat==3)
!        write(10,*) 'MESH "Layer Fiss endommages fissures" dimension 3 ElemType ',nomfis,' Nnode ',dime,''
!        write(10,*) '# Color 1.00000 0.00000 1.000000 '
!        write(10,*) 'Elements'
!        do ief = 1, inf3
!          iek = iek + 1
!          write(10,*) iek, kconecf3(ief,1:dime)+nnt
!        end do
!        write(10,*) 'end elements'
!      end if

      if (inf3 > 0) then
          ! Elements endommages-fissures (etat==3)
          do icouc = 1,size(lcouc,1)
              if (count(lcouc(icouc)%elem>0)>0) then
                  write(10,*) 'MESH "Layer Endo-Fiss ', icouc ,' ouvert" dimension 3 ElemType ',nomfis,' Nnode ', dime, ''
                  write(10,*) '# Color ',lcouc(icouc)%coul , ''
                  write(10,*) 'Elements'
                  do ie = 1, inf3
                      if (lcouc(icouc)%elem(ie)/=0) then
                          iek = iek + 1
                          write(10,*) iek, kconecf3(lcouc(icouc)%elem(ie),1:dime)+nnt
                      end if
                  end do
                  write(10,*) 'end elements'
              end if
          end do
      end if

      if (inf2 > 0) then
        ! Elements fissures refermmes   (etat==2)
        write(10,*) 'MESH "Layer Fiss referme" dimension 3 ElemType ',nomfis,' Nnode ',dime,''
        write(10,*) '# Color 0.000000 1.000000 0.000000 '
        write(10,*) 'Elements'
        do ief = 1, inf2
          iek = iek + 1
          write(10,*) iek, kconecf2(ief,1:dime)+nnt+inoeu3
        end do
        write(10,*) 'end elements'
      end if

          ! Elements fissures (etat==1)
      if (inf1 > 0) then
       do icouc = 1,size(lcouc,1)
        if (count(lcouc(icouc)%elem>0)>0) then
          write(10,*) 'MESH "Layer Fiss ', icouc ,' ouvert" dimension 3 ElemType ',nomfis,' Nnode ', dime, ''

          write(10,*) '# Color ',lcouc(icouc)%coul , ''
          write(10,*) 'Elements'

          do ie = 1, inf1
            if (lcouc(icouc)%elem(ie)/=0) then
               iek = iek + 1
               write(10,*) iek, kconecf1(lcouc(icouc)%elem(ie),1:dime)+nnt+inoeu3+inoeu2
            end if
          end do
          write(10,*) 'end elements'
        end if
       end do
      end if
      ! fin des elements fissures 2D (Line) ou 3D (Triangle)

       end if
       ! fin de tracer des fissures

       deallocate(vcorf1,vcorf2, kconecf1, kconecf2, left1, left2)

       do i = 1, size(lcouc,1)
           deallocate(lcouc(i)%elem)
       end do

       close(10)

    end subroutine ecriture_fiss

!*********************************************************!
!  Fonction pour determiner les coordonnees des fissures  !
!*********************************************************!
    subroutine fisscoor(ie, coorfiss, NP)

       use variables, only : vcor, ktypel, nomtype, infele, kconec, dime, inorm
       use math, only : norme, intersect_mat, cross, det, inv
       use initialisation, only : init_vec, init_mat

       implicit none

       ! In
       integer, intent(in) :: ie

       ! Out
       real*8, dimension(:,:), allocatable, intent(out) :: coorfiss
       integer, intent(out) :: NP

       ! Commun
       integer :: i, id
       real*8, dimension(dime) :: nf, N
       real*8, dimension(:,:), allocatable :: vcore
       integer, dimension(:), allocatable :: vec
       character(len=5) :: typel

       !---- Cas 2D
       real*8, dimension(dime) :: xym
       real*8, dimension(dime,dime) :: xl
       real*8 :: aire, peri, lcar, lx, mx, ly, my

       !---- Cas 3D
       real*8, dimension(dime) :: PG

       integer, dimension(:,:), allocatable :: face, inter, ARRETE, INTFAC
       integer :: nface, nar, intf, NART, IL
       integer :: ns, j, k, l, face1, face2, NI, NJ

       integer, dimension(:), allocatable :: F1, F2, ar
       real*8, dimension(dime) :: PI, PJ, P1, P2, VI, V1, N1, V2, N2, V
       real*8 :: A1, B1, C1, D1, A2, B2, C2, D2, AN, BN, CN, DN, &
                & XG, YG, ZG, XI, YI, ZI, XJ, YJ, ZJ, DENO, KI, KJ, dM
       real*8, dimension(dime,dime) :: M
       real*8, dimension(:,:), allocatable :: PF

!-----------------------------------------------------------!

        !----- On recupere le nom de l'element
        typel = nomtype(ktypel(ie))
            

        !----- On recupere la normale (nf) à la fissure (N = normale normee)
        nf = inorm(ie,:)
        N = nf/norme(nf)

        !----- On recupere les coordonnées des noeuds de l'element
        call init_vec(vec, infele(ktypel(ie))%nnel)
        vec = kconec(ie,1:infele(ktypel(ie))%nnel)
        call init_mat(vcore,size(vec,1),dime)
        do i = 1, infele(ktypel(ie))%nnel
            vcore(i,:) = vcor(:,vec(i))
        end do

        !----- On calcule le centre de gravite de l'element
        do id = 1,dime
            PG(id) = sum(vcore(:,id))/size(vcore,1)
        end do

        if (dime==2) then !     <===================== 2D !!!

            xym = PG
            
            !----- cas d'un element triangulaire (les noeuds externes sont numerotes : 1, 2, 3)
            if (typel(1:4) == 'MBT3' .or. typel(1:4) == 'MBT6') then
                !----- Aire de l'element
                aire = 0.5 * (vcore(2,1)-vcore(1,1))*(vcore(3,2)-vcore(1,2)) &
                      & -(vcore(3,1)-vcore(1,1))*(vcore(2,2)-vcore(1,2))

                !----- Perimetre de l'element
                peri = sqrt((vcore(2,1)-vcore(1,1))**2 + (vcore(2,2)-vcore(1,2))**2) + &
                     & sqrt((vcore(3,1)-vcore(2,1))**2 + (vcore(3,2)-vcore(2,2))**2) + &
                     & sqrt((vcore(1,1)-vcore(3,1))**2 + (vcore(1,2)-vcore(3,2))**2)
            
            !----- cas d'un element quadrangulaire (les noeuds externes sont numerotes : 1, 2, 3, 4)
            elseif (typel(1:4) == 'MBQ4' .or. typel(1:4) == 'MBQ8') then

                !----- Aire de l'element
                aire = 0.5 * (vcore(2,1)-vcore(1,1))*(vcore(3,2)-vcore(1,2)) &
                      & -(vcore(3,1)-vcore(1,1))*(vcore(2,2)-vcore(1,2))

                !----- Perimetre de l'element
                peri = sqrt((vcore(2,1)-vcore(1,1))**2 + (vcore(2,2)-vcore(1,2))**2) + &
                     & sqrt((vcore(3,1)-vcore(2,1))**2 + (vcore(3,2)-vcore(2,2))**2) + &
                     & sqrt((vcore(4,1)-vcore(3,1))**2 + (vcore(4,2)-vcore(3,2))**2) + &
                     & sqrt((vcore(1,1)-vcore(4,1))**2 + (vcore(1,2)-vcore(4,2))**2)
            
            endif
            ! Longueur caracteritique
            lcar = 4*aire / peri

            lx = nf(1)
            mx = nf(2)
            ly = -mx
            my = lx

            ! Coordonnees des noeuds de la longueur caracteristiques
            xl(1,:) = (/ xym - lcar/2 * (/ly,my/) /)
            xl(2,:) = (/ xym + lcar/2 * (/ly,my/) /)

          call init_mat(coorfiss,size(xl,1),size(xl,2))

          NP = size(xl,1)
          coorfiss = xl

        elseif (dime==3) then !  <==================== 3D !!!
            ! Definition des liaisons de l'element
            call init_mat(face,size(infele(ktypel(ie))%face,1),size(infele(ktypel(ie))%face,2))
            face  = infele(ktypel(ie))%face
            nface = size(face,1)

            ! On cherche le nombre total d'arretes de l'element
            ns = 1
            do i = 1, nface-1
               ns = ns*i
            end do

            call init_mat(inter,ns,2)
            call init_mat(ARRETE,ns,2)
            call init_mat(INTFAC,ns,2)

            !--- Remplacer - nchoosek - MATLAB
            ! inter = nchoosek(1:1:nface,2)
            j = 0
            do k = 1, nface-1
              do l = k+1, nface
                j = j + 1
                inter(j,:) = (/k, l/)
              end do
            end do

            nar = 0

            ! Boucle sur toutes les faces trouvees
            do intf = 1, ns
                face1 = inter(intf,1)
                face2 = inter(intf,2)

                call init_vec(F1,size(face,2))
                call init_vec(F2,size(face,2))

                F1 = face(face1,:)
                F2 = face(face2,:)

                call init_vec(ar,size(face,2)-1)
                ar = intersect_mat(F1,F2)

                if (count(ar/=0)==2) then
                    nar = nar + 1
                    ARRETE(nar,:) = ar
                    INTFAC(nar,:) = (/face1,face2/)
                end if

                deallocate(F1, F2, ar)
            end do

            NART = nar
            NP = 0

            call init_mat(PF,nar,dime)

            ! Boucle sur toutes les arretes trouvees
            do IL = 1 , NART
                !
                NI = ARRETE(IL,1)  ! Numero du noeud I
                NJ = ARRETE(IL,2)  ! Numero du noeud J
                PI = vcore(NI,:)   !  Coordonnees du noeud I
                PJ = vcore(NJ,:)   !  Coordonnees du noeud J

                call init_vec(F1,size(face,2))
                call init_vec(F2,size(face,2))

                F1 = face(INTFAC(IL,1),:)
                F2 = face(INTFAC(IL,2),:)

                do id = 1, dime
                   P1(id) = sum(vcore(F1,id))/size(vcore(F1,1))  ! barycentre de la face 1
                   P2(id) = sum(vcore(F2,id))/size(vcore(F2,1))  ! barycentre de la face 2
                end do

                ! Vecteur directeur de l'arrete
                VI = PJ - PI
                ! Vecteur complementaire du plan face 1
                V1 = P1-PI
                N1 = cross(VI,V1)
                N1 = N1/norme(N1)
                ! Vecteur complementaire du plan face 2
                V2 = P2-PI
                N2 = cross(VI,V2)
                N2 = N2/norme(N2)
                !
                XG=PG(1); YG=PG(2); ZG=PG(3);
                XI=PI(1); YI=PI(2); ZI=PI(3);
                XJ=PJ(1); YJ=PJ(2); ZJ=PJ(3);
                !
                ! Plan Face 1 passant par PI
                !
                A1=N1(1)
                B1=N1(2)
                C1=N1(3)
                D1=-(A1*XI+B1*YI+C1*ZI)
                !
                ! Plan Face 2 passant par PI
                !
                A2=N2(1)
                B2=N2(2)
                C2=N2(3)
                D2=-(A2*XI+B2*YI+C2*ZI);
                !
                ! Plan de fissuration PN passant par PG
                !
                AN=N(1)
                BN=N(2)
                CN=N(3)
                DN=-(AN*XG + BN*YG + CN*ZG)

                ! On cherche le point intersection de ces trois plans
                M(1,:) = (/ A1 , B1 , C1 /)
                M(2,:) = (/ A2 , B2 , C2 /)
                M(3,:) = (/ AN , BN , CN /)

                dM = det(M)

                if (abs(dM) > 1.d-10) then
                    DENO = sqrt(AN**2 + BN**2 + CN**2)
                    ! Position des noeuds PI et PJ par rapport au plan (PN) :
                    ! si k1.k2<0 PI et PJ sont de part et d'autre du plan (PN)
                    KI = (AN*XI+BN*YI+CN*ZI+DN)/DENO
                    KJ = (AN*XJ+BN*YJ+CN*ZJ+DN)/DENO

                    if ((KI*KJ)<=0) then
                        V = (/ -D1,-D2,-DN /)
                        NP = NP + 1
                        PF(NP,:) = matmul(inv(M),V)
                    else
                       !print*,'intersection hors element'
                    end if
                else
                    !print*,'Pas de point d'intersection'
                end if

                deallocate(F1, F2)
            end do

            call init_mat(coorfiss,NP,size(PF,2))

            coorfiss = PF(1:NP,:)

            deallocate(inter, ARRETE, INTFAC)
        end if

        deallocate(vcore, vec)

    end subroutine fisscoor

!********************************************************!

    function fiss_ouv() result(ouver)

!********************************************************!
!          Calcul de l'ouverture moyenne d'une fissure   !
!********************************************************!

        use variables, only : dime, inorm, vcor, kconec, infele, &
            & ktypel, vsol, ietat, nelt, infnoe, fissBA
        use math, only : norme
        use initialisation, only : init_vec, init_mat

        implicit none

        real*8 :: XG, YG, ZG, AN, BN, CN, DN, DENO, dep1, dep2
        real*8, dimension(dime) :: nf, NN, PG
        real*8, dimension(:,:), allocatable :: depNds, vcore
        real*8, dimension(:), allocatable :: depNorm, dist

        integer :: ie, i, nnzp, nnzn, noel, ino
        integer, dimension(:), allocatable :: vec
        integer, dimension(:,:), allocatable :: ddl

        logical, dimension(:), allocatable :: cotep

        ! Out
        real*8, dimension(nelt) :: ouver

!*********************************************************

    ouver = 0.0d0

    do ie = 1, nelt
    
      if ((ietat(ie)==1) .or. (fissBA==1 .and. ietat(ie)==3 .and. norme(inorm(ie,:))/=0)) then

        nf = inorm(ie,:)
        NN = nf/norme(nf)

        noel = infele(ktypel(ie))%nnel
        call init_vec(dist,noel)
        call init_vec(cotep,noel)

        call init_vec(vec, noel)
        vec = kconec(ie,1:noel)

        call init_mat(vcore,dime,noel)
        do i = 1, noel
                vcore(:,i) = vcor(:,vec(i))
        end do

        if (dime==1) then
                stop 'FIDES_fiss_ouv : cas dime = 1 non prevu !'
        else
            call init_mat(ddl,noel,dime)
            do ino = 1, noel
            ddl(ino,:) = infnoe(vec(ino))%dln(1:dime)
        end do

                call init_mat(depNds, dime, dime*noel)
                depNds(1,:) = vsol(ddl(:,1))
                depNds(2,:) = vsol(ddl(:,2))
                if (dime==3) then
                    depNds(3,:) = vsol(ddl(:,3))
                end if
        end if

        call init_vec(depNorm, size(depNds,2))
        depNorm = matmul(NN,depNds)

        !----- On calcule les distances des sommets au plan de fissuration
        !----- pour detecter la position de ces points par rapport au plan de
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

        do i = 1, infele(ktypel(ie))%nnel
                if (cotep(i)) then
                        dep1 = dep1 + (depNorm(i))/nnzp
                else
                        dep2 = dep2 + (depNorm(i))/nnzn
                end if
        end do

        ouver(ie) = dep1 - dep2

        deallocate(cotep, dist, vec, vcore, depNds, depNorm, ddl)

      end if
    end do

    end function fiss_ouv
!********************************************************!


!------------------------------------------------------------------------------------!
!                  Ecriture du fichier de lecture de GNUPLOT                         !
!------------------------------------------------------------------------------------!

    subroutine ecriture_evgb(nomfichier)

      use variables, only : dirresu

      implicit none

      character(len=*), intent(in) :: nomfichier
      character(len=len_trim(nomfichier)) :: nometude

      nometude = nomfichier(1:len_trim(nomfichier)-5)//'.reac'

      open(unit = 4, file = 'evolglob.script', form = 'formatted', status = 'unknown')
      write(4,'(a)') '#evolglob.script'
      write(4,'(a)') 'set autoscale'
      write(4,'(a)') 'set pointsize 0.35'
      write(4,'(a)') 'set title "Force-Deplacement direction x"'
      write(4,'(a)') 'plot "'//dirresu//trim(nometude)//'" using 1:4 notitle wi linespoints pt 7'
      write(4,'(a)') 'pause -1'
      write(4,'(a)') 'set title "Force-Deplacement direction y"'
      write(4,'(a)') 'plot "'//dirresu//trim(nometude)//'" using 2:5 notitle wi linespoints pt 7'
      write(4,'(a)') 'pause -1'
      write(4,'(a)') 'set title "Force-Deplacement direction z"'
      write(4,'(a)') 'plot "'//dirresu//trim(nometude)//'" using 3:6 notitle wi linespoints pt 7'
      write(4,'(a)') 'pause -1'

      close(4)

    end subroutine ecriture_evgb

!***********************************************************************************!
!               Calcul de la force global et du deplacement global                  !
!***********************************************************************************!

subroutine evolglob()

    use variables
    use initialisation, only : init_mat, init_vec
    use lib_elem, only : elem_B, elem_kloce2

    implicit none

    integer :: ie, ipg, iloc, idir, i, j
    integer :: ndle, ndln, nc1, npg, nnel, idim1
    real*8 :: epais, detj, pdet, poids, vprel(idmax), vfg(ndlt)
    integer, dimension(:), allocatable :: kv
    real*8, dimension(:), allocatable :: vfint, vpg, vsg
    real*8, dimension(:,:), allocatable :: ksig, vn, vb, vsig
    character(len=5) :: typel

    dpglobal = 0.
    foglobal = 0.
    vfg = 0.

    do ie = 1, nelt
        !-----  Recuperation des informations sur les elements
        nnel = infele(ktypel(ie))%nnel
        ndln = infele(ktypel(ie))%ndln
        nc1  = infele(ktypel(ie))%ncgr

        !----- vecteur kloce de localisation pour assemblage
        call elem_kloce2( ie, ndle)

        !if ((count(kcond(kloce) .eqv. .true.) > 0) .or. (count(vfcg(kloce)/=0.) > 0)) then

        !----- initialisations vecteurs elementaires (calcules aux noeuds)
        call init_vec(vfint,ndle)

        !----- Caracteristiques geometriques des elements
        vprel = vprelg(kprop(ie),1:idmax)

        epais = 1.
        if ((dime == 2).and.(vprel(2)==3)) epais = vprel(idmax-2) ! En contraintes planes (CP)

        !----- Recuperation des contraintes et integration
        !----- Boucle sur les points de Gauss
        allocate(vpg(size(infele(ktypel(ie))%W)))
        allocate(ksig(size(infele(ktypel(ie))%Q,1),size(infele(ktypel(ie))%Q,2)))

        vpg = infele(ktypel(ie))%W
        ksig = infele(ktypel(ie))%Q
        npg = size(ksig,1)

        !----- localisation des donnees elementaires dans les vecteurs globaux -----
        if (calco == 'NOEUD') then
            idim1 = nnel*nc1
        elseif (calco == 'GAUSS') then
            idim1 = npg*nc1
        end if

        iloc = kcont(ie)

        call init_mat(vsig, nc1, npg)
        vsig  = reshape(vcont(iloc : (iloc+ idim1 - 1)),(/ nc1, npg/))

        do ipg = 1, npg

            poids = vpg(ipg)
            !----- Calcul des fonctions d'interpolation et des derivees

            call elem_B(vn,vb,detj,ksig(ipg,:),ie)
            pdet = detj*poids*epais

            if (calco == 'NOEUD') then
                !Cas non implante
            else
                !----- Recuperation des contraintes aux points de Gauss si calcul en ces points
                typel = nomtype(ktypel(ie))

                if ((typel=='EJQ4').or.(typel=='EJQ6').or.(typel=='EJT6')) then
                    if (dime == 2) then
                        !call init_vec(vsg, size(vsig,1)-1)
                        call init_vec(vsg, size(vsig,1))
                    elseif (dime == 3) then
                        !call init_vec(vsg, size(vsig,1)-3)
                        call init_vec(vsg, size(vsig,1))
                    endif
                else
                    call init_vec(vsg, size(vsig,1))
                endif

                vsg = vsig(1:size(vsg),ipg)
            end if

            !----- Integration des contraintes
            vfint = vfint + pdet*matmul(transpose(vb),vsg)

            deallocate(vn,vb,vsg)
        end do

        vfg(kloce(1:ndle)) = vfg(kloce(1:ndle)) + vfint

        deallocate(ksig, vfint, vpg, vsig)
        !end if

    end do

    if (reacn .eqv. .false.) then
        j = 1
        do idir = 1, size(dirimpo)
            if (idir == dirimpo(1) .or. idir == dirimpo(2) .or. idir == dirimpo(3)) then
                if ((npilot == 1).or.(npilot == 2).or.(npilot==3)) then
                    call init_vec(kv,min(count(vcond/=0),count(kcond .eqv. .true.)))
                    do i = 1, nnt
                        if ((kcond(infnoe(i)%dln(idir)) .eqv. .true.) .and. (vcond(infnoe(i)%dln(idir)) /=0 )) then
                            kv(j) = infnoe(i)%dln(idir)
                            j = j + 1
                        end if
                    end do
                else
                    call init_vec(kv,count(vfcg/=0.))
                    do i = 1,nnt
                        if (vfcg(infnoe(i)%dln(idir)) /=0. ) then
                            kv(j) = infnoe(i)%dln(idir)
                            j = j +1
                        end if
                    end do
                end if
                dpglobal(idir) = sum(vsol(kv))/count(kv /= 0)
                foglobal(idir) = sum(vfg(kv))
                deallocate(kv)
            end if
        end do
    else
        dpglobal(diri) = vsol(infnoe(noe1)%dln(diri))-vsol(infnoe(noe2)%dln(diri))
        foglobal(diri) = vfg(infnoe(noe1)%dln(diri))

        !dpglobal(diri) = vsol(infnoe(990)%dln(diri))-vsol(infnoe(351)%dln(diri))  ! maillage WST fi
        !foglobal(diri) = vfg(infnoe(7)%dln(diri))                                 ! maillage WST fi
        !dpglobal(diri) = vsol(infnoe(1002)%dln(diri))-vsol(infnoe(426)%dln(diri)) ! maillage WST mo
        !foglobal(diri) = vfg(infnoe(7)%dln(diri))                                 ! maillage WST mo
        !dpglobal(diri) = vsol(infnoe(1660)%dln(diri))-vsol(infnoe(33)%dln(diri))  ! maillage WST tgr
        !foglobal(diri) = vfg(infnoe(862)%dln(diri))                               ! maillage WST tgr
        !dpglobal(diri) = vsol(infnoe(284)%dln(diri))-vsol(infnoe(14)%dln(diri))   ! maillage Bruhwiller Gc/Rt
        !foglobal(diri) = vfg(infnoe(7)%dln(diri))                                 ! maillage Bruhwiller Gc:Rt
        !dpglobal(diri) = vsol(infnoe(296)%dln(diri))-vsol(infnoe(26)%dln(diri))   ! maillage Bruhwiller Gc2
        !foglobal(diri) = vfg(infnoe(13)%dln(diri))                                 ! maillage Bruhwiller Gc2
        !dpglobal(diri) = vsol(infnoe(284)%dln(diri))-vsol(infnoe(14)%dln(diri))   ! maillage Bruhwiller Gc3/Gc4
        !foglobal(diri) = vfg(infnoe(7)%dln(diri))                                 ! maillage Bruhwiller Gc3/Gc4
    end if

end subroutine evolglob

!--------------------------------------------------------------------!
!  SUBROUTINE DE SAUVEGARDE DU DERNIER PAS DE CALCUL REALISE         !
!  SAUVEGARDE EGALEMENT L'HETEROGENEITE DU MATERIAU DANS LE CAS      !
!  D'UN CALCUL PROBABILISTE                                          !
!--------------------------------------------------------------------!
    subroutine ecriture_resu(nomfichier,iopt)

       use variables, only : dirresu, kprop, vprelg, vsol, vcont, vnoli, &
                           & vrint, interf, ietatpg, histetatpg1, histetatpg2, &
                           & fiss, ietat, histetat1, histetat2, inorm

       implicit none

       character(len=*), intent(in) :: nomfichier
       integer, intent(in), optional :: iopt

       integer, parameter :: numfich=1
       character(len=len_trim(nomfichier)) :: nomsauve
       integer :: ierr, option

       !----- Option de lecture
       if (present(iopt)) then; option = iopt; else; option = 1; end if

       !----- Ouverture ou creation du fichier
       nomsauve = nomfichier(1:len_trim(nomfichier)-5)//'.resu'

       open(unit = numfich, file = dirresu//nomsauve, form = 'formatted', status = 'replace', iostat=ierr)

       if (ierr/=0) then
          open(unit = numfich, file = dirresu//nomsauve, form = 'formatted', status = 'new', iostat=ierr)
          if (ierr/=0) stop
       end if

       if (option>=1) then

          !----- Sauvegarde des tableaux de proprietes mecaniques
          write(numfich,*) size(kprop,1)
          write(numfich,*) kprop
          write(numfich,*) size(vprelg,1), size(vprelg,2)
          write(numfich,*) vprelg

       end if
       if (option<=1) then

          !----- Sauvegarde des tableaux generaux
          write(numfich,*) size(vsol,1)
          write(numfich,*) vsol

          write(numfich,*) size(vcont,1)
          write(numfich,*) vcont

          write(numfich,*) size(vnoli,1)
          write(numfich,*) vnoli

          write(numfich,*) size(vrint,1)
          write(numfich,*) vrint

          !----- Stockage specifique aux elements d'interface (à modifier dans le futur)
          if (interf==1) then
              write(numfich,*) size(ietatpg,1), size(ietatpg,2)
              write(numfich,*) ietatpg

              write(numfich,*) size(histetatpg1,1), size(histetatpg1,2)
              write(numfich,*) ietatpg

              write(numfich,*) size(histetatpg2,1), size(histetatpg2,2)
              write(numfich,*) ietatpg

          end if

          !----- Stockage specifique au modele de beton fissurant (à modifier dans le futur)
          if (fiss==1) then
              write(numfich,*) size(ietat,1)
              write(numfich,*) ietat

              write(numfich,*) size(histetat1,1)
              write(numfich,*) ietat

              write(numfich,*) size(histetat2,1)
              write(numfich,*) ietat

              write(numfich,*) size(inorm,1), size(inorm,2)
              write(numfich,*) inorm

          end if

       end if

       !----- Fermeture du fichier de sauvegarde
       close(numfich)

    end subroutine ecriture_resu

    subroutine lecture_resu(nomfichier,iopt)
!-------------------------------------------------------------------------!
!        Routine de lecture du fichier binaire de résultats au            !
!        dernier pas de calcul.                                           !
!-------------------------------------------------------------------------!
       use variables, only : dirresu, kprop, vprelg, vsol, vcont, vnoli, &
                           & vrint, interf, ietatpg, &
                           & fiss, ietat, inorm
       use initialisation

       implicit none

       character(len=*), intent(in) :: nomfichier
       integer, intent(in), optional :: iopt

       integer, parameter :: numfich=1
       character(len=len_trim(nomfichier)+5) :: nomsauve
       integer :: ierr, nc, nl, bid, option

       !----- Option de lecture
       if (present(iopt)) then; option = iopt; else; option = 1; end if

       !----- Ouverture ou creation des fichiers
       nomsauve = trim(nomfichier)//'.resu'
       open(unit = numfich, file = dirresu//nomsauve, form= 'formatted', status='old', action='read', iostat=ierr)

       if (ierr/=0) then
           print*,"lecture_alea : le fichier ",dirresu//nomsauve," n'existe pas !"
           print*,' '
           stop
       end if

       !----- LECTURE DE LA DISTRIBUTION ALEATOIRE DE PROPRIETES SI ELLE EXISTE
       if (option>=1) then
            !----- Liberation de l'espace memoire des anciens tableaux
            deallocate(kprop,vprelg)

            !----- Lecture de tailles et creation des tableaux de proprietes mecaniques
            read(numfich,*)  nc
            call init_vec(kprop,nc)
            read(numfich,*) kprop

            read(numfich,*)  nc, nl
            call init_mat(vprelg,nc,nl)
            read(numfich,*) vprelg

       end if
       !----- LECTURE COMPLETE DES DONNEES
       if (option<=1) then
           !----- Lecture des tableaux generaux
           read(numfich,*) bid
           read(numfich,*) vsol

           read(numfich,*) bid
           read(numfich,*) vcont

           read(numfich,*) bid
           read(numfich,*) vnoli

           read(numfich,*) bid
           read(numfich,*) vrint

           !----- Lecture specifique aux elements d'interface (à modifier dans le futur)
           if (interf==1) then
               read(numfich,*) bid, bid
               read(numfich,*) ietatpg

               read(numfich,*) bid, bid
               read(numfich,*) ietatpg

               read(numfich,*) bid, bid
               read(numfich,*) ietatpg

           end if

           !----- Lecture specifique au calcul beton fissurant (à modifier dans le futur)
           if (fiss==1) then
               read(numfich,*) bid
               read(numfich,*) ietat

               read(numfich,*) bid
               read(numfich,*) ietat

               read(numfich,*) bid
               read(numfich,*) ietat

               read(numfich,*) bid, bid
               read(numfich,*) inorm
           end if

       end if

       close(numfich)

    end subroutine lecture_resu


    subroutine inicont()
!---------------------------------------------------------------------
!  Subroutine d'initialisation des contraintes dans les elements
!  - utilisée pour distribuer des auto-contraintes initiales dans le
!    beton
!                    EN COURS DE DEVELOPPEMENT !!!
!---------------------------------------------------------------------
    use variables, only : nelt, nomtype, ktypel, infele, kcont, vcont, &
                        & dime, gpinicont, kprop0
    use initialisation
    use math

    implicit none

    !***** vecteurs locaux
    real*8, dimension(:), allocatable   :: vpg, vsi
    real*8, dimension(1) :: vs
    real*8, dimension(:,:), allocatable :: ksig, vsig

    !***** Variables locales
    integer :: ie, nc, ipg, npg, idim, iloc

    character(len=5) :: type

    !----- Boucle sur les elements
    do ie = 1, nelt
       type = nomtype(ktypel(ie))

       !----- Pour les elements massifs uniquement
       if ((type(1:2)=='MB').or.(type(1:2)=='MT')) then

           !----- Pour les elements du(des) groupe(s) concerne(s) par les auto-contraintes
           if (allocated(gpinicont) .and. (count(gpinicont==kprop0(ie))==1)) then

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

           !----- Initialisation des vecteurs contraintes locaux -----
           call init_mat(vsig, nc, npg)
           call init_vec(vsi,nc)
           !call init_vec(vs,1+ie)

           !----- Boucle sur les points de gauss
           do ipg = 1, npg

               !----- Cas des auto-contraintes du béton
               !vs = 2.d0 * (random(1+ie)-0.5)
               vs = 2.d0 * (randgp(1,(/ie/))-.5d0)

               !----- calcul des contraintes
               if (dime==2) then
                   vsi = vs(1) * (/1.d0,1.d0,0.d0/)
               else
                   vsi = vs(1) * (/1.d0,1.d0,1.d0,0.d0,0.d0,0.d0/)
               end if
               vsig(1:nc,npg) = vsi

           end do

           !----- Stockage des contraintes dans le vecteur global
           vcont(iloc:(iloc+idim-1)) = reshape(vsig,(/ nc*npg /))

           deallocate(ksig,vpg,vsig,vsi)

           end if
       else
           !stop 'Erreur subroutine inicont : type d''element non considere'
       end if
    end do

    end subroutine inicont

    subroutine princ(vmat0,vmat1,P1)
!---------------------------------------------------------------------
!  Subroutine de calcul et d'ordonnancement des contraintes
!  principales
!
!  En entree : mat0
!    vmat0 = matrice (sous forme vecteur) initiale dont on cherche les v.p.
!  En sortie : vmat1, P1
!    vmat1 = vecteur contenant les v.p. ordonnees (v(1)>v(2)>v(3))
!    P1 = matrice de passage rep ini => rep principal
!---------------------------------------------------------------------
    use variables, only : dime
    use math, only : Jacobi, find_vec, norme, sort

    implicit none

    !--- Variables IN
    real*8, dimension(:), intent(in) :: vmat0

    !--- Variables OUT
    real*8, dimension(dime,dime), intent(out) :: P1
    real*8, dimension(dime), intent(out)      :: vmat1

    !--- Variables locales
    real*8, dimension(dime,dime) :: MM, VP
    real*8, dimension(dime)      :: vmat
    integer, dimension(dime)     :: ilocp
    integer :: NROT

    if (dime == 1) then
        stop 'utilitaire - princ :: impossible en 1D'
    elseif (dime == 2) then
        MM(1,:) = (/ vmat0(1), vmat0(3) /)
        MM(2,:) = (/ vmat0(3), vmat0(2) /)
        call Jacobi(MM,dime,vmat,VP,NROT)
        ilocp = sort(vmat)
        vmat1 = (/ vmat(ilocp(2)), vmat(ilocp(1)) /)
        P1(:,1) = VP(:,ilocp(2))
        P1(:,2) = VP(:,ilocp(1))
    elseif (dime == 3) then
        MM(1,:) = (/ vmat0(1), vmat0(4), vmat0(6) /)
        MM(2,:) = (/ vmat0(4), vmat0(2), vmat0(5) /)
        MM(3,:) = (/ vmat0(6), vmat0(5), vmat0(3) /)
        call Jacobi(MM,dime,vmat,VP,NROT)
        ilocp = sort(vmat)
        vmat1 = (/ vmat(ilocp(3)), vmat(ilocp(2)), vmat(ilocp(1)) /)
        P1(:,1) = VP(:,ilocp(3))
        P1(:,2) = VP(:,ilocp(2))
        P1(:,3) = VP(:,ilocp(1))
    end if

    end subroutine princ


    subroutine ddlnode(ndlt)
!---------------------------------------------------------------------
!  La routine construit le vecteur kndll(nnt), variable globale,
!  qui pour chaque noeud du maillage donne son nombre de ddl.
!  La somme de tous les termes de kndll donne le nombre total
!  de ddl du probleme (ndlt).
!---------------------------------------------------------------------
    use variables, only : kconec, ktypel, nnt, nelt, infele, &
        & knddl, infnoe

    implicit none

    !--- Variable OUT
    integer, optional, intent(out) :: ndlt

    integer :: ie, noel, ino, inode, ndll, i, idll

    !--- Allocation et initialisation de kndll
    allocate(knddl(nnt)); knddl = 0

    !--- Boucle sur les elements
    do ie = 1, nelt
        noel = infele(ktypel(ie))%nnel
        ndll = infele(ktypel(ie))%ndln
        !--- Boucle sur les noeuds de d'element
        do ino = 1, noel
            inode = kconec(ie,ino)
            !--- kndll ne conserve que le plus grand nbre de ddl du noeud considere
            if (ndll > knddl(inode)) knddl(inode)=ndll
        end do
    end do

    ndlt = sum(knddl)

    ! Determination du numero des ddls par noeuds
    idll = 0
    allocate(infnoe(nnt))

    do ino = 1, nnt
       allocate(infnoe(ino)%dln(knddl(ino)))
       infnoe(ino)%dln = (/(i+idll, i=1, knddl(ino))/)
       idll = idll + knddl(ino)
    end do

    end subroutine ddlnode

end module utilitaire

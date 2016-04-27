!----------------------------------------------------------------------------------------------------------
! MODULE: Maillage
!
!> @author JL Tailhan
!> @author J Foulliaron  (version 1.0 - 2011)
!
!> @brief
!> Routines creant la structure de matrice creuse associee au maillage donne
!----------------------------------------------------------------------------------------------------------
module Maillage


contains


!------------------------------------------------------------------------------------------------------
!> @author JL Tailhan
!> @author J Foulliaron  (version 1.0 - 2011)
!
!> @brief
!> Routine creant la structure de matrice creuse associee au maillage donne     (josk1) \n\n
!
!> @details
!> #### DESCRIPTION:
!> Composee d'une boucle sur les noeuds du maillage (divisee en 5 parties) \n
!> puis d'une phase finale de creation de la matrice creuse et sortie dans un fichier texte si besoin
!------------------------------------------------------------------------------------------------------
subroutine  Maillage_CSRC(nomfichier,nzmax36)

    !**************************************************************************************
    !* DEFINITIONS
    !**************************************************************************************

    !--------------------------------------------------------------------------------------
    !--- Modules generaux
    use variables
    use formatCSRC
    use lib_elem

    implicit none

    !--------------------------------------------------------------------------------------
    !--- Variables in/out
    character(len = *), intent(in) :: nomfichier
    integer(kind=intkind), intent(in)::nzmax36  !-- intkind defini dans variables.f

    !--------------------------------------------------------------------------------------
    !--- Variables locales
    integer(kind=intkind) :: ta_kloce, ta_kloce_max, ind, A1, A2, A3, pos
    integer :: i, j, nddl, pos2, pos3, lol, ta_redu, ndle
    integer :: iel, jel, ino, jno, idl, ele
    integer :: ia(ndlt+1),ja(nzmax36), pcent(100), nbelemautour(nnt)
    integer, dimension(:), allocatable :: kloce_periph, kloce_interm

    character(len = len_trim(nomfichier)) :: nomfichierCSRC

    !**************************************************************************************
    !* CORPS DE PROCEDURE
    !**************************************************************************************

    !real*8 :: tdeb,tfin

    !--------------------------------------------------------------------------------------
    !--- Preparation ...

    !--------------------------------------------------------------------------------------
    !--- Creation du tableau permettant l'affichage du compteur en pourcentage
    do i=1,100
        pcent(i)=nnt/100*i
    enddo

    !--------------------------------------------------------------------------------------
    !--- Donnees initiales
    ia = 0
    ja = 0
    ia(1)=1
    pos=1;ind=1;ta_kloce=0;
    ta_kloce_max=10000

    allocate(kloce_periph(ta_kloce_max)) ; kloce_periph = 0
    allocate(kloce_interm(ta_kloce_max)) ; kloce_interm = 0
    allocate(listelemautour(nnt))

    !call CPU_TIME(tdeb)

    lol=1

    !--------------------------------------------------------------------------------------
    !--- On comptabilise le nombre d'elements auxquels appartiennent les noeuds (nbelemautour)
    nbelemautour=0
    do iel = 1,nelt
        do jno=1,infele(ktypel(iel))%nnel
            nbelemautour(kconec(iel,jno))=nbelemautour(kconec(iel,jno))+1
        enddo
    enddo

    !--------------------------------------------------------------------------------------
    !--- On cree la liste des elements connectes a chaque noeud (listelemautour)
    do ino = 1,nnt
        allocate(listelemautour(ino)%el(nbelemautour(ino)))
        listelemautour(ino)%el = 0
    enddo

    nbelemautour=0

    do iel=1,nelt
        do jno=1,infele(ktypel(iel))%nnel
            nbelemautour(kconec(iel,jno))=nbelemautour(kconec(iel,jno))+1
            listelemautour(kconec(iel,jno))%el(nbelemautour(kconec(iel,jno)))=iel
        enddo
    enddo

    if (ipost==0) then

        !----------------------------------------------------------------------------------
        !--- BOUCLE SUR LES NOEUDS
        do ino = 1, nnt

            if (ino==pcent(lol)) then
                !print*,'progression',lol,'% . .'
                lol=lol+1
            endif

            !------------------------------------------------------------------------------
            !--- 1) Calcul du nbre de ddls en relation avec ceux du noeud ino : ta_kloce
            !       (a partir des ddl des elements connectes au noeud)

            ta_kloce=0
            do jel = 1, nbelemautour(ino)
                ele = listelemautour(ino)%el(jel)
                ta_kloce = ta_kloce + infele(ktypel(ele))%nnel*infele(ktypel(ele))%ndln
            enddo
            if (ta_kloce>ta_kloce_max) stop 'pas assez despace pour kloce_periph: augmentez ta_kloce_max dans maillage.f ligne 42'

            !------------------------------------------------------------------------------
            !--- 2) Liste des ddls en relation avec ceux du noeud ino : kloce_periph

            pos2 = 1
            do jel = 1, nbelemautour(ino)
                call elem_kloce2(listelemautour(ino)%el(jel),ndle)
                kloce_periph(pos2:pos2+(ndle-1)) = kloce(1:ndle)
                pos2 = pos2 + ndle
            enddo

            !------------------------------------------------------------------------------
            !--- 3) Phase de suppression des doublons de kloce_periph
            kloce_periph(1:ta_kloce) = QuickSort(kloce_periph(1:ta_kloce))
            kloce_interm(1) = kloce_periph(1)

            pos3=1
            do i=2,ta_kloce
                if (kloce_periph(i) /= kloce_periph(i-1)) then
                    pos3=pos3+1
                    kloce_interm(pos3) = kloce_periph(i)  !on ne garde que les valeurs qui ne se repètent pas
                endif

            enddo
            ta_redu = pos3

            !------------------------------------------------------------------------------
            !--- 4)  Creation des valeurs de ja et ia associees au groupe de ddls associes au noeud ino
            !        grace a kloce_periph (indices pr remplissage du triangle inf de supervkg)

            nddl = knddl(ino)
            do idl = pos,pos+nddl-1  ! boucle sur les ddls du noeud ino
                                     ! idl = numero reel du ddl du noeud ino traite

                do j = 1,ta_redu   ! boucle sur les ddl peripheriques au ddl numero idl
                    if (idl > kloce_interm(j)) then
                        ja(ind) = kloce_interm(j)
                        ind = ind + 1
                        if (ind > nzmax36) stop 'pas assez d espace alloue pour ja,..aumgenter ta_max dans FIDES.f '
                    endif
                enddo
                ia(idl+1) = ind

            enddo

            pos = pos + nddl

        enddo
        !--- FIN DE LA BOUCLE SUR LES NOEUDS
        !----------------------------------------------------------------------------------

        deallocate(kloce_periph)
        deallocate(kloce_interm)

        !----------------------------------------------------------------------------------
        !--- 5 ) Phase finale de création de la matrice au format CSRC

        !----------------------------------------------------------------------------------
        !--- Ecriture de la matrice dans un fichier (obsolete, opt1 defini dans variables.f)
        if (opt1==1) then
            pos=0
            nomfichierCSRC = nomfichier(1:len_trim(nomfichier)-5)//'.csrc'
            print*,'Creation du fichier texte CSRC, fichier :',dircsrc//nomfichierCSRC
            open(unit = 2, file = dircsrc//nomfichierCSRC, form = 'formatted', status = 'replace')

            A1=(ndlt+1)/ta_texte
            pos=1
            A2=(ind-1)/ta_texte
            A3= nnt
            write(2,'(a7)') 'tabIA'
            write(2,*) ndlt+1
            do i=1,A1
                write(2,*) ia(pos:pos+ta_texte-1)
                pos=pos+ta_texte
            enddo
            write(2,*) ia(pos:(ndlt+1))

            pos=1
            write(2,'(a7)') 'tabJA'
            write(2,*) ind-1
            do i=1,A2
                write(2,*) ja(pos:pos+ta_texte-1)
                pos=pos+ta_texte
            enddo
            write(2,*) ja(pos:(ind-1))
            close(unit=2)

            nomfichierCSRC = nomfichier(1:len_trim(nomfichier)-5)//'.elem'
            print*,'creation du fichier texte CSRC, fichier :',dircsrc//nomfichierCSRC
            open(unit = 2, file = dircsrc//nomfichierCSRC, form = 'formatted', status = 'replace')

            write(2,'(a7)') 'kelem'
            do i=1,A3
                write(2,*)size(listelemautour(i)%el)
                write(2,*) listelemautour(i)%el
            enddo
            close(unit=2)
            print*,'Création fichier csrc et elem terminée'
        endif

        !----------------------------------------------------------------------------------
        !--- Creation de la structure tant attendue de la matrice creuse au format CSRC
        call CSRC_Init(supervkg,ndlt,ind-1)

        supervkg%ia=ia
        supervkg%ja=ja(1:ind-1)

    endif

end subroutine Maillage_CSRC


!------------------------------------------------------------------------------------------------------
!> @author JL Tailhan
!> @author J Foulliaron  (version 1.0 - 2011)
!
!> @brief
!> Fonction reordonnant les valeurs d'un vecteur d'entier dans l'ordre croissant
!
!> @details
!> #### DESCRIPTION:
!>  A faire
!------------------------------------------------------------------------------------------------------
recursive function QuickSort(InList) result(OutList)

    !**************************************************************************************
    !* DEFINITIONS
    !**************************************************************************************

    implicit none

    !--------------------------------------------------------------------------------------
    !--- Variables locales
    INTEGER,DIMENSION(:) :: InList
    INTEGER,DIMENSION(size(InList,1)) :: OutList
    INTEGER,DIMENSION(size(InList,1)) :: SupList, OrderedSupList, InfList, OrderedInfList
    INTEGER :: pivot
    INTEGER :: i, j, InfListSize, SupListSize

    !**************************************************************************************
    !* CORPS DE PROCEDURE
    !**************************************************************************************

    !--------------------------------------------------------------------------------------
    !--- S'il ne reste qu'un element dans la liste, on arrête la recursion
    if(size(InList,1) < 2) then
        OutList(1) = Inlist(1)
    else

        !----------------------------------------------------------------------------------
        !--- Le pivot sera le premier element de la liste
        pivot = InList(1)

        !----------------------------------------------------------------------------------
        !--- On trie la liste
        InfListSize = 0
        SupListSize = 0
        do i = 2, size(InList,1)
            if(InList(i) < Pivot) then
                InfListSize = InfListSize + 1
                InfList(InfListSize) = InList(i)
            elseif(InList(i) >= Pivot) then
                SupListSize = SupListSize + 1
                SupList(SupListSize) = InList(i)
            endif
        enddo

        !----------------------------------------------------------------------------------
        !--- On recompose la liste
        if(InfListSize < 1) then
            OrderedSupList = QuickSort(SupList(1:SupListSize))
            OutList = (/ Pivot, (OrderedSupList(j), j=1,SupListSize) /)
        elseif(SupListSize < 1) then
            OrderedInfList = QuickSort(InfList(1:InfListSize))
            OutList = (/ (OrderedInfList(j), j=1,InfListSize), Pivot /)
        else
            OrderedInfList = QuickSort(InfList(1:InfListSize))
            OrderedSupList = QuickSort(SupList(1:SupListSize))
            OutList = (/ (OrderedInfList(j), j=1,InfListSize), Pivot, (OrderedSupList(j), j=1,SupListSize) /)
        endif

    endif

end function QuickSort


!------------------------------------------------------------------------------------------------------
!> @author JL Tailhan
!> @author J Foulliaron  (version 1.0 - 2011)
!
!> @brief
!> Routine réordonnant les valeurs d'un vecteur d'entier dans l'ordre croissant
!> (non utilisee - prise d'internet)
!
!> @details
!> #### DESCRIPTION:
!>  A faire
!------------------------------------------------------------------------------------------------------
recursive function QuickSort2(InList) result(OutList)

    !**************************************************************************************
    !* DEFINITIONS
    !**************************************************************************************

    implicit none

    !--------------------------------------------------------------------------------------
    !--- Variables locales
    INTEGER*8,DIMENSION(:) :: InList
    INTEGER*8,DIMENSION(size(InList,1)) :: OutList
    INTEGER*8,DIMENSION(size(InList,1)) :: SupList, OrderedSupList, InfList, OrderedInfList
    INTEGER*8:: pivot
    INTEGER*8 :: i,j, InfListSize, SupListSize

    !--------------------------------------------------------------------------------------
    !--- S'il ne reste qu'un élément dans la liste, on arrête la récursion
    if(size(InList,1) < 2) then
        OutList(1) = Inlist(1)
    else

        !----------------------------------------------------------------------------------
        !--- Le pivot sera le premier élément de la liste
        pivot = InList(1)

        !----------------------------------------------------------------------------------
        !--- On trie la liste
        InfListSize = 0
        SupListSize = 0
        do i = 2, size(InList,1)
            if(InList(i) < Pivot) then
                InfListSize = InfListSize + 1
                InfList(InfListSize) = InList(i)
            elseif(InList(i) >= Pivot) then
                SupListSize = SupListSize + 1
                SupList(SupListSize) = InList(i)
            endif
        enddo

        !----------------------------------------------------------------------------------
        !--- On recompose la liste
        if(InfListSize < 1) then
            OrderedSupList = QuickSort2(SupList(1:SupListSize))
            OutList = (/ Pivot, (OrderedSupList(j), j=1,SupListSize) /)
        elseif(SupListSize < 1) then
            OrderedInfList = QuickSort2(InfList(1:InfListSize))
            OutList = (/ (OrderedInfList(j), j=1,InfListSize), Pivot /)
        else
            OrderedInfList = QuickSort2(InfList(1:InfListSize))
            OrderedSupList = QuickSort2(SupList(1:SupListSize))

            OutList = (/ (OrderedInfList(j), j=1,InfListSize), Pivot, (OrderedSupList(j), j=1,SupListSize) /)
        endif
    endif
end function QuickSort2

!------------------------------------------------------------------------------------------------------

end module Maillage

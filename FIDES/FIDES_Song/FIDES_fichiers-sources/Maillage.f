module Maillage
contains

!----------------------------------------------------------------------------------------------------------------------------------------------!
! Routine creeant la structure de matrice creuse associee au maillage donne     (josk1)                          !
!                                                                                                                                                                                  !
! Compose d'une boucle sur les noeuds du maillage (divisee en 5 parties)                                              !
!puis d'une phase finale de creation de la matrice creuse et sortie dans un fichier texte si besoin !
!----------------------------------------------------------------------------------------------------------------------------------------------!

subroutine  Maillage_CSRC(nomfichier,nzmax36)

!nzmax36 : ! nombre de position maximal a creer pour les valeurs non nulles (estimation)

use variables
use formatCSRC
use lib_elem
implicit none
 
 character(len = *), intent(in) :: nomfichier
 integer(kind=8), intent(in)::nzmax36

 integer(kind=8) :: ta_kloce, ta_kloce_max, ind, A1, A2, A3, pos
 integer :: i,j,k,y,a,b,z, nddl, pos2, toto, lol, repere, ta_redu, err1, ndle
 integer :: ia(ndlt+1), ja(nzmax36), pcent(100), nbelemautour(nnt)
 integer, dimension(:), allocatable :: kloce_periph, kloce_interm

 character(len = len_trim(nomfichier)) :: nomfichierCSRC

 real*8 :: tdeb,tfin

 ! Preparation . . 

 ! Creation du tableau permettant l'affichage du compteur en pourcentage
 do i=1,100
   pcent(i)=nnt/100*i
 end do
 
 ! Donnee initiales

 ia = 0
 ja = 0
 ia(1)=1      
 pos=1;ind=1;ta_kloce=0;
 ta_kloce_max=10000
 allocate(kloce_periph(ta_kloce_max))
 allocate(kloce_interm(ta_kloce_max))
 allocate(listelemautour(nnt))

 call CPU_TIME(tdeb)

nbelemautour=0

lol=1

do i=1,nelt
   do j=1,infele(ktypel(i))%nnel
     nbelemautour(kconec(i,j))=nbelemautour(kconec(i,j))+1
  end do
end do

do i = 1,nnt
   allocate(listelemautour(i)%el(nbelemautour(i)))
end do

nbelemautour=0

do i=1,nelt
   do j=1,infele(ktypel(i))%nnel    
      nbelemautour(kconec(i,j))=nbelemautour(kconec(i,j))+1
      listelemautour(kconec(i,j))%el(nbelemautour(kconec(i,j)))=i
   end do
end do

!============BOUCLE SUR LES NOEUDS======!

if (ipost==0) then

do i = 1, nnt

      if (i==pcent(lol)) then
        ! print*,'progression',lol,'% . .'
         lol=lol+1
      end if

      nddl=knddl(i)
      ta_kloce=0

!  1) Creation de listelem associee au noeud i
      do j = 1, nbelemautour(i)
        toto = listelemautour(i)%el(j)
        ta_kloce=ta_kloce+infele(ktypel(toto))%nnel*infele(ktypel(toto))%ndln
      end do

      if (ta_kloce>ta_kloce_max) stop 'pas assez despace pour kloce_periph: augmentez ta_kloce_max dans maillage.f ligne 42'

!  2) creation de kloce_periph( ddls en relation avec ceux du noeud i) 

       pos2 = 1
      
       do z = 1, nbelemautour(i)
            call elem_kloce2(listelemautour(i)%el(z),ndle)
            kloce_periph(pos2:pos2+(ndle-1)) = kloce(1:ndle)
            pos2 = pos2 + ndle     
       end do

!3) Phase de suppression des doublons de kloce_periph

        kloce_periph(1:ta_kloce) = QuickSort(kloce_periph(1:ta_kloce))
        repere=1

        kloce_interm(1)=kloce_periph(1)


        do y=2,ta_kloce
  
           if(kloce_periph(y)/=kloce_periph(y-1))then
            
             repere=repere+1
             kloce_interm(repere) = kloce_periph(y)  !on ne garde que les valeurs qui ne se repètent pas 
          
           end if

        end do
        ta_redu = repere

!4)  Creation des valeur de ja et ia associees au groupe de ddls associes au noeud i , grace a kloce_periph 

        do a=pos,pos+nddl-1    !boucle sur les ddls du noeud i
                               !a= numero reel du ddl du noeud i traite 
                
                do b =1,ta_redu ! boucle sur les ddl peripheriques au ddl numero a 
                     
                      if( a>kloce_interm(b))then
                        
                           ja(ind)=kloce_interm(b)
                           ind=ind+1
	
                           if (ind>nzmax36) stop 'pas assez despace alloue pour ja,..aumgenter ta_max dans FIDES.f '

                      end if          
                   
                 end do       
                 ia(a+1)=ind  
         end do

       pos=pos+nddl

   end do 
   !============== FIN DE LA BOUCLE SUR LES NOEUDS===========!


   deallocate(kloce_periph)
   deallocate(kloce_interm)
 
! 5 ) Phase finale de création de la matrice au format CSRC

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
   end do
   write(2,*) ia(pos:(ndlt+1))

   pos=1
   write(2,'(a7)') 'tabJA'
   write(2,*) ind-1
   do i=1,A2
       write(2,*) ja(pos:pos+ta_texte-1)
       pos=pos+ta_texte
   end do
   write(2,*) ja(pos:(ind-1))   
   close(unit=2)
  
   nomfichierCSRC = nomfichier(1:len_trim(nomfichier)-5)//'.elem'
   print*,'creation du fichier texte CSRC, fichier :',dircsrc//nomfichierCSRC 
   open(unit = 2, file = dircsrc//nomfichierCSRC, form = 'formatted', status = 'replace')

    write(2,'(a7)') 'kelem'
    do i=1,A3
       write(2,*)size(listelemautour(i)%el)
       write(2,*) listelemautour(i)%el
    end do
   close(unit=2)
   print*,'Création fichier csrc et elem terminée' 
end if   
 
   ! Creation de la structure tant attendue de la matrice creuse au format CSRC

   call CSRC_Init(supervkg,ndlt,ind-1)

   supervkg%ia=ia
   supervkg%ja=ja(1:ind-1)
   
end if

end subroutine Maillage_CSRC

!========================================================================================================!

!--------------------------------------------------------------------------------------------------------!
! Routine reordonnant les valeurs d'un vecteur d'entier dans l'ordre croissant                           !      
!--------------------------------------------------------------------------------------------------------!

recursive function QuickSort(InList) result(OutList)
    INTEGER,DIMENSION(:) :: InList
    INTEGER,DIMENSION(size(InList,1)) :: OutList
    INTEGER,DIMENSION(size(InList,1)) :: SupList, OrderedSupList, InfList, OrderedInfList
    INTEGER:: pivot
    INTEGER :: i,j, InfListSize, SupListSize
 
    ! S'il ne reste qu'un element dans la liste, on arrête la recursion
    if(size(InList,1) < 2) then
       OutList(1) = Inlist(1)
    else
       ! Le pivot sera le premier element de la liste
       pivot = InList(1)
 
       ! On trie la liste
       InfListSize = 0
       SupListSize = 0
       do i = 2, size(InList,1)
          if(InList(i) < Pivot) then
             InfListSize = InfListSize + 1
             InfList(InfListSize) = InList(i)
          elseif(InList(i) >= Pivot) then
             SupListSize = SupListSize + 1
             SupList(SupListSize) = InList(i)
          end if
       enddo
 
       ! On recompose la liste
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
    end if
end function QuickSort
!=================================================================================================!

!----------------------------------------------------------------------------------------------------------------------------------!
! Routine réordonnant les valeurs d'un vecteur d'entier dans l'ordre croissant (utilisée dans Maillage_CSRC : (prise d'internet)   !      
!----------------------------------------------------------------------------------------------------------------------------------!

recursive function QuickSort2(InList) result(OutList)
    INTEGER*8,DIMENSION(:) :: InList
    INTEGER*8,DIMENSION(size(InList,1)) :: OutList
    INTEGER*8,DIMENSION(size(InList,1)) :: SupList, OrderedSupList, InfList, OrderedInfList
    INTEGER*8:: pivot
    INTEGER*8 :: i,j, InfListSize, SupListSize
 
    ! S'il ne reste qu'un élément dans la liste, on arrête la récursion
    if(size(InList,1) < 2) then
       OutList(1) = Inlist(1)
    else
       ! Le pivot sera le premier élément de la liste
       pivot = InList(1)
 
       ! On trie la liste
       InfListSize = 0
       SupListSize = 0
       do i = 2, size(InList,1)
          if(InList(i) < Pivot) then
             InfListSize = InfListSize + 1
             InfList(InfListSize) = InList(i)
          elseif(InList(i) >= Pivot) then
             SupListSize = SupListSize + 1
             SupList(SupListSize) = InList(i)
          end if
       enddo
 
       ! On recompose la liste
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
    end if
end function QuickSort2


end module Maillage

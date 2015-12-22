module sparse_matrix

!*********************************************************************************!
!            Module permettant la gestion des matrices creuses                    !
!*********************************************************************************!

!---------------------------------------------------------------------------------!
!                        Type sparse : Stockage au format CSR                     !
!---------------------------------------------------------------------------------!

   type sparse

    ! Matrice A stokee au fromat CSR
    real*8,  dimension(:), allocatable :: a
    integer, dimension(:), allocatable :: ja, ia

    ! Nombre de valeurs non nulles maximales
    integer :: nzmax

    ! nombre de ligne de A
    integer :: n

    ! nombre de colonne de A
    integer :: nc

   end type sparse

contains

!====================== Foncions et procedures generales =========================!

!---------------------------------------------------------------------------------!
!       Allocation et initialisation d'une matrice au format CSR                  !
!   Entrees :                                                                     !
!   - nbl : nombre de ligne de la matrice                                         !
!   - nbc : nombre de colonne de la matrice                                       !
!   - nzmax : nombre de valeurs non nulles maximales                              !
!   - per : pourcentage de nnz (compris entre 0 et 1)                             !
!                                                                                 !
!   Sortie :                                                                      !
!   - mat : matrice au format CSR                                                 !
!---------------------------------------------------------------------------------!

   subroutine init_sparse(mat, nbl, nbc, nzmax, per)

      use variables, only : nnt, dime

      implicit none

      type (sparse), intent(out) :: mat
      integer, intent(in) :: nbl, nbc
      real*8, intent(in), optional :: per
      integer, intent(in), optional :: nzmax

      real*8 :: per1
      integer :: nzmax1

      if (present(per) .and. per<1. .and. per>0.) then; per1 = per; else; per1 = 1.; end if
      if (dime==2) then
         per1 = 7.0*nnt**(-0.85)   ! PTS  : deduire a partir des tests, jusqu'a 100000 ddl
      else
         per1 = 5.0*nnt**(-0.8)
      end if

      if (present(nzmax)) then; nzmax1 = nzmax; else; nzmax1 = nint(nbl*per1*nbc); end if
      ! nzmax1 < 2.e9 : environ 1.5e6 - 2.e6 ddl

      allocate(mat%a(nzmax1)); allocate(mat%ja(nzmax1)); allocate(mat%ia(nbl+1))
      mat%a = 0.
      mat%ja = 0
      mat%ia = 1
      mat%nzmax = nzmax1
      mat%n = nbl
      mat%nc = nbc

   end subroutine init_sparse

!---------------------------------------------------------------------------------!
!       Fonction qui retourne la position de A(i,j) dans les tableaux ja at a     !
!                                                                                 !
!  Entrees :                                                                      !
!  - A : matrice au format CSR                                                    !
!  - i : indice de ligne                                                          !
!  - j : indice de colonne                                                        !
!---------------------------------------------------------------------------------!

   function get_indice(A,i,j) result(indice)

       implicit none

       ! retourne la valeur 0, si A(i,j) n'est pas stocke sinon sa position

       type(sparse), intent(in) :: A
     integer, intent(in) :: i
     integer, intent(in) :: j
     integer :: indice

     integer :: indicei,indicej, nbcl
    ! nbcl : nb de nnz pour la ligne i
    logical :: zero

    ! verfification des indices
    if(i > A%n) then
        print*, 'get_indice : l''indice de ligne est trop grand'
        stop
    end if

    if(j > A%nc) then
        print*, 'get_sparse : l''indice de colonne est trop grand'
        stop
    end if

    nbcl = A%ia(i+1)-A%ia(i)
    indicei = A%ia(i)
    zero = .true.

    do indicej = 1, nbcl
       if(A%ja(indicei + indicej - 1) == j) then
          zero = .false.
          indice = indicei + indicej - 1
          exit ! sortie de la boucle
       end if
    end do

    if (zero) indice = 0

   end function get_indice



!---------------------------------------------------------------------------------!
!                Fonction qui retourne la valeur A(i,j)           !
!  Entrees :                                                                      !
!  - A : matrice au format CSR                                                    !
!  - i : indice de ligne                                                          !
!  - j : indice de colonne                                                        !
!---------------------------------------------------------------------------------!


   function get_sparse(A,i,j) result(valeur)

      implicit none

    type(sparse), intent(in) :: A
    integer, intent(in) :: i
    integer, intent(in) :: j
    real*8 :: valeur

    integer :: ind

    ! verfification des indices
    if(i > A%n) then
        print*, 'get_sparse : l''indice de ligne est trop grand'
        stop
    end if

    if(j > A%nc) then
        print*, 'get_sparse : l''indice de colonne est trop grand'
        stop
    end if

    ind = get_indice(A,i,j)

    if(ind/=0) then
        valeur = A%a(ind)
    else
        valeur = 0.
    end if

   end function

!-----------------------------------------------------------------------------------!
!                Procedure qui modifie la valeur A(i,j)                 !
!   Entrees :                                                                       !
!  - A : matrice au format CSR                                                      !
!  - i : indice de ligne                                                            !
!  - j : indice de colonne                                                          !
!  - valeur : nouvelle valeur de A(i,j)                                             !
!                                                                                   !
!  Sortie :                                                                         !
!  - A : matrice apres modification                                                 !
!-----------------------------------------------------------------------------------!

   subroutine set_sparse(A,i,j,valeur)

      implicit none

      type(sparse), intent(inout) :: A
    integer, intent(in) :: i
    integer, intent(in) :: j
    real*8, intent(in) :: valeur

    integer :: tailleA, ind, indicei,indicej, nbcl
    ! nbcl : nb de nnz pour la ligne i
    logical :: zero

    ! verfification des indices
    if(i > A%n) then
        print*, 'set_sparse : l''indice de ligne est trop grand'
        stop
    end if

    if(j > A%nc) then
        print*, 'set_sparse : l''indice de colonne est trop grand'
        stop
    end if

    nbcl = A%ia(i+1)-A%ia(i)
    indicei = i+1
    tailleA = A%ia(A%n+1) - 1
    zero = .true.

    ! On determine si A(i,j)/=0 et si oui on determine sa place dans a

    indicej = get_indice(A,i,j)

    if(indicej/=0) zero = .false.  ! La valeur etait stockee

    ! Premier cas : valeur = 0
    if(valeur == 0.) then

         ! Si A(i,j)/=0, On retire A(i,j) des tableaux d'indices et de a
         if(.not.zero .and. tailleA/=0) then   ! A(i,j) etait non nul

            do while (indicej <= tailleA)
                 A%ja(indicej) = A%ja(indicej+1)
                 A%a(indicej) = A%a(indicej + 1)
                 indicej = indicej + 1
            end do

            do while (indicei <= a%n+1)
                 A%ia(indicei) = A%ia(indicei) - 1
                 indicei = indicei + 1
            end do

         end if

         ! Si A(i,j) etait deja nul on ne fait rien

     end if

     ! Deuxieme cas : valeur /=0
     if(valeur /=0.) then

         ! Si A(i,j) = 0., On ajoute A(i,j) aux tableaux d'indices et a a

         if(tailleA >= A%nzmax) then
             print*, 'set_sparse : La place memoire allouee pour la matrice sparse A n ''est pas suffisante'
             stop
         end if

         if(zero) then   ! A(i,j) etait nul

            ! On determine a quel endroit on se situe dans la matrice

            indicej = A%ia(i)

            do while(A%ja(indicej ) < j .and. indicej < A%ia(i+1))
                indicej = indicej + 1
            end do


            ! On introduit le terme A(i,j)

            ind = tailleA

          do while (ind >= indicej)

             A%ja(ind+1) = A%ja(ind)
             A%a(ind+1) = A%a(ind)
             ind = ind - 1
          end do

          A%ja(indicej) = j
          A%a(indicej) = valeur

            do while (indicei <= a%n+1 )
               A%ia(indicei) = A%ia(indicei) + 1
               indicei = indicei + 1
            end do

         ! Si A(i,j) etait deja nul on remplace la valeur dans a
         else
            A%a(indicej) = valeur

         end if

     end if


   end subroutine


!------------------------------------------------------------------------------------!
!         Procedure qui met a zero toutes les valeurs d'une ligne            !
!  Entrees :                                                                         !
!  - A : matrice au format CSR                                                       !
!  - i : indice de la ligne                                                          !
!                                                                                    !
!  Sortie :                                                                          !
!  - A : matrice apres modification                                                  !
!------------------------------------------------------------------------------------!

   subroutine del_row(A, i)

       implicit none

       type(sparse), intent(inout) :: A
       integer, intent(in) :: i

       integer :: nb  ! Nombre d'element non nul de la ligne
       integer :: j, tailleA

       nb = A%ia(i+1)-A%ia(i)
       tailleA = A%ia(A%n+1) - 1

       do j = A%ia(i), tailleA - nb
           A%ja(j) = A%ja(j+nb)
           A%a(j) = A%a(j+nb)
       end do

       do j = i+1, A%n + 1
           A%ia(j) = A%ia(j) - nb
       end do

       A%a(tailleA-nb+1 : A%nzmax) = 0.
       A%ja(tailleA-nb+1 : A%nzmax) = 0


   end subroutine del_row

!------------------------------------------------------------------------------------!
!         Procedure qui met a zero toutes les valeurs d'une colonne          !
!  Entrees :                                                                         !
!  - A : matrice au format CSR                                                       !
!  - j : indice de la colonne                                                        !
!                                                                                    !
!  Sortie :                                                                          !
!  - A : matrice apres modification                                                  !
!------------------------------------------------------------------------------------!

   subroutine del_col(A, j)

       implicit none

       type(sparse), intent(inout) :: A
       integer, intent(in) :: j

       integer, dimension(A%n+1) :: indice
       integer :: i, nb, tailleA, k, ind
       integer, dimension(A%n) :: modifc

       tailleA = A%ia(A%n+1) - 1
       nb = 0

       ! On determine la place des colonnes
       do i = 1, A%n

           ind = get_indice(A,i,j)
           if( ind /=0) then
                indice(nb+1) = ind - nb
                  nb = nb+1

           end if

           modifc(i) = nb

       end do

       indice(nb+1) = tailleA + 1

       ! on decalle les indices
       do i = 1, nb
           do k = indice(i), indice(i+1)-1
              A%ja(k) = A%ja(k+i)
              A%a(k) = A%a(k+i)
           end do
       end do

       A%ia(2: A%n+1) = A%ia(2: A%n+1) - modifc

       A%a(tailleA-nb+1 : A%nzmax) = 0.
       A%ja(tailleA-nb+1 : A%nzmax) = 0


   end subroutine del_col

!----------------------------------------------------------------------------------------------!
!   Procedure qui permet l'elimination d'eventuels zeros stockes dans la matrice               !
!              (elle permet de retablir un "reel" stockage au format CSR)                      !
!  Entree :                                                                                    !
!  - A : matrice avant modification                                                            !
!                                                                                              !
!  Sortie :                                                                                    !
!  - A : matrice au format CSR obtenue apres modification                                      !
!----------------------------------------------------------------------------------------------!
   subroutine del_zeros(A)

    implicit none

    type(sparse), intent(inout) :: A

    integer, dimension(A%ia(A%n+1)-1) :: zeros
      integer :: i, j, nb, tailleA, k
      integer, dimension(A%n) :: modifc

       tailleA = A%ia(A%n+1) - 1
       nb = 0

     !  mem = A%ia

       ! On determine la place des colonnes
       do i = 1, A%n

           do j = A%ia(i), A%ia(i+1)-1

              if( A%a(j) == 0.) then
            zeros(nb+1) = j - nb
                nb = nb+1

              end if

       end do

           modifc(i) = nb

       end do

       zeros(nb+1) = tailleA + 1

       ! on decalle les indices
       do i = 1, nb
           do k = zeros(i), zeros(i+1)-1
              A%ja(k) = A%ja(k+i)
              A%a(k) = A%a(k+i)
           end do
       end do

       A%ia(2: A%n+1) = A%ia(2: A%n+1) - modifc(1:A%n)

       A%a(tailleA-nb+1 : A%nzmax) = 0.
       A%ja(tailleA-nb+1 : A%nzmax) = 0

   end subroutine del_zeros

!----------------------------------------------------------------------------------------!
!               Procedure qui permet de au minimum reduire la place memoire              !
!                          reservee a la matrice (nzmax)                                 !
!  Entree :                                                                              !
!  - mat : matrice stockee au format CSR                                                 !
!                                                                                        !
!  Sortie :                                                                              !
!  - mat : matrice apres reduction de nzmax                                              !
!----------------------------------------------------------------------------------------!

   subroutine Sparse_reduc_taille(mat)

       implicit none

       type(sparse), intent(inout) :: mat
       real*8, dimension(:), allocatable :: aux1
       integer, dimension(:), allocatable :: aux2
       integer :: taille

       taille = mat%ia(mat%n+1) - 1
       taille = int(taille*1.25) ! Ceci est à prendre avec des pincettes ... (JLT+Josk1)

       ! Reduction de la taille du tableau de valeurs
       allocate(aux1(taille));
       aux1 = mat%a(1:taille);
       deallocate(mat%a);
       allocate(mat%a(taille));
       mat%a = aux1
       deallocate(aux1)

       ! Reduction de la taille du tableau des indices des colonnes
       allocate(aux2(taille));
       aux2 = mat%ja(1:taille);
       deallocate(mat%ja);
       allocate(mat%ja(taille));
       mat%ja = aux2
       deallocate(aux2)

       mat%nzmax = taille

   end subroutine Sparse_reduc_taille


!----------------------------------------------------------------------------------------!
!         Procedure qui libere la place memoire reservee a une matrice stockee   !
!                                         au format CSR                                  !
!  Entree :                                                                              !
!  - mat : matrice stockee au format CSR pour laquelle la place memoire doit etre libere !
!----------------------------------------------------------------------------------------!

   subroutine free_sparse(mat)

       implicit none

       type(sparse), intent(inout) :: mat
       deallocate(mat%a,mat%ia,mat%ja)

   end subroutine free_sparse



!============================ Fonctions liées à l'assemblage =============================!

!-----------------------------------------------------------------------------------------!
!         Procedure qui permet d'assembler plus rapidement la matrice globale             !
!                          a partir des matrices elementaires                             !
!  Entrees :                                                                              !
!  - matg : matrice globale a assembler                                                   !
!  - matel : matrice elementaire a inserer                                                !
!  - indi : indices de ligne ou doivent etre inserer les elements de matel                !
!  - indj : indices de colonne ou doivent etre inserer les elements de matel              !
!                                                                                         !
!  Sortie :                                                                               !
!  - matg : matrice globale apres insertion de la matrice elementaire                     !
!-----------------------------------------------------------------------------------------!

   subroutine assem_sparse(matg, matel, indi, indj)

    use math, only : sort
    implicit none

    type(sparse), intent(inout) :: matg
    real*8, dimension(:,:), intent(in) :: matel ! matrice elementaire
    integer, dimension(:), intent(in) :: indi, indj ! tableau d'indice

    integer, dimension(size(indi)) :: pi
    integer, dimension(size(indj)) :: pj

    ! On reordonne les indices
    pi = sort(indi)
    pj = sort(indj)

    ! On ajoute val a matg
    call add_sparse(matg, indi(pi), indj(pj), matel(pi,pj))


   end subroutine assem_sparse

!------------------------------------------------------------------------------------------!
!              Procedure qui permet d'ajouter des valeurs plus rapidement              !
!                          dans une matrice stockee au format CSR                          !
!  Entrees :                                                                               !
!  - A : matrice stockee au format CSR                                                     !
!  - indi : indices de ligne de A ou les valeurs sont ajoutees                             !
!  - indj : indices de colonne de A ou les valeurs sont ajoutees                           !
!  - val : matrice contenant les valeurs a ajouter                                         !
!                                                                                          !
!  Sortie :                                                                                !
!  - A : matrice apres ajout des valeurs                                                   !
!------------------------------------------------------------------------------------------!
   subroutine add_sparse (A, indi, indj, val)

      implicit none

      type(sparse), intent(inout) :: A
      integer, dimension(:), intent(in) :: indi, indj
      real*8, dimension(:,:), intent(in) :: val

      integer :: n,m
      integer, dimension(size(val)+1) :: newP    ! Position du coefficient a inserer
      integer, dimension(size(val)) :: CnewP ! Numero de la colonne
      real*8, dimension(size(val)) :: Vnew
      integer, dimension(size(indi)) :: ModifC  ! Nombre de colonnes ajoutees
      integer :: indice, i, j, indicej, ind
      real*8 :: valeurA

      n = size(indi)
      m = size(indj)
      indice = 0

      ! On determine parmis les positions dans K qui doivent stocker les valeurs de val celles qui
      ! doivent etre crees

      do i = 1, n

         indicej = A%ia(indi(i))

         do j = 1, m

            valeurA = get_sparse(A, indi(i), indj(j))
            ind = get_indice(A,indi(i),indj(j))

            if(ind/= 0 ) then
              ! La valeur etait stockee
               if( val(i, j)/=0) then
              A%a(ind) = valeurA + val(i,j)
                !  call set_sparse(A, indi(i), indj(j), valeurA)
               end if

            else

               if( val(i, j)/=0) then
                ! La valeur n'etait pas stockee et val /=0
                do while(A%ja(indicej ) < indj(j) .and. indicej < A%ia(indi(i)+1))
                    indicej = indicej + 1
                end do
                indice = indice + 1
                newP(indice) = indicej
                CnewP(indice) = indj(j)
                Vnew(indice) = val(i,j)
               end if

            end if

         end do

         ModifC(i)=indice

      end do

      if(A%ia(A%n+1)+indice -1 > A%nzmax) then
          print*, 'add_sparse : La place memoire allouee pour la matrice sparse A n ''est pas suffisante'
          stop
      end if

      newP(indice+1) = A%ia(A%n+1)

      ! On actualise A%ia en fonction du nb d'element qui va etre ajoute a chaque ligne
      do i = 1, n - 1
         A%ia(indi(i)+1:indi(i+1)) = A%ia(indi(i)+1:indi(i+1)) + ModifC(i)
      end do

      A%ia(indi(n)+1: A%n+1) = A%ia(indi(n)+1:A%n+1) + ModifC(n)

      do i = indice, 1, -1

         do j = newP(i+1)-1,newP(i),-1
             A%a(j+i) = A%a(j)
             A%ja(j+i) = A%ja(j)
         end do

      end do


      do i = 1, indice
        ! La valeur n'etait pas stockee

         A%a(newP(i)+i-1) = Vnew(i)
         A%ja(newP(i)+i-1) = CnewP(i)


      end do

   end subroutine add_sparse

!------------------------------------------------------------------------------------------!
!              Procedure qui permet d'inserer plusieurs valeurs                            !
!                     dans une matrice stockee au format CSR                               !
!  Entrees :                                                                               !
!  - A : matrice stockee au format CSR                                                     !
!  - indi : indices de ligne de A ou les valeurs sont inserees                             !
!  - indj : indices de colonne de A ou les valeurs sont inserees                           !
!  - val : tableau contenant les valeurs a inserees                                        !
!                                                                                          !
!  Sortie :                                                                                !
!  - A : matrice apres insertion des valeurs                                               !
!------------------------------------------------------------------------------------------!
  subroutine insert_values(A, indi, indj, val)

  ! indi(i) : indice de ligne de la ieme valeur a inserer
  ! indj(i) : indice de colonne de la ieme valeur a inserer
  ! val(i)  : ieme valeur a inserer

     implicit none

     type(sparse), intent(inout) :: A
     integer, dimension(:) :: indi ,indj ! Indice de ligne et de colonne
     real*8, dimension(:) :: val ! Tableau de valeurs

      integer :: n
      integer, dimension(size(val)+1) :: newP    ! Position du coefficient a inserer
      integer, dimension(size(val)) :: CnewP ! Numero de la colonne
      real*8, dimension(size(val)) :: Vnew
      integer, dimension(size(indi)) :: ModifC  ! Nombre de colonnes ajoutees
      integer :: indice, i, j, indicej, ind

      if(size(indi)/=size(indj)) then
         print*, 'insert_values :les tableaux d''indices doivent avoir la meme taille'
         stop
      end if

      n = size(indi)

      indice = 0

      ! On determine parmi les positions de A celles qui doivent etre crees

      do i = 1, n

         indicej = A%ia(indi(i))

         ind = get_indice(A,indi(i),indj(i))

         ! La valeur etait stockee
         if(ind /= 0 ) then
               A%a(ind) = val(i)

         ! La valeur n'etait pas stockee
         else

               if( val(i)/=0) then
                ! La valeur n'etait pas stockee et val /=0
                do while(A%ja(indicej ) < indj(i) .and. indicej < A%ia(indi(i)+1))
                    indicej = indicej + 1
                end do
                indice = indice + 1
                newP(indice) = indicej
                CnewP(indice) = indj(i)
                Vnew(indice) = val(i)
               end if

            end if


         ModifC(i)=indice

      end do

      if(A%ia(A%n+1)+indice > A%nzmax) then
          print*, 'insert_values : La place memoire allouee pour la matrice sparse A n ''est pas suffisante'
          stop
      end if

      newP(indice+1) = A%ia(A%n+1)

      ! On actualise A%ia en fonction du nb d'element qui va etre ajoute a chaque ligne
      do i = 1, n - 1
         A%ia(indi(i)+1:indi(i+1)) = A%ia(indi(i)+1:indi(i+1)) + ModifC(i)
      end do

      A%ia(indi(n)+1: A%n+1) = A%ia(indi(n)+1:A%n+1) + ModifC(n)

      do i = indice, 1, -1

         do j = newP(i+1)-1,newP(i),-1
             A%a(j+i) = A%a(j)
             A%ja(j+i) = A%ja(j)
         end do

      end do


      do i = 1, indice
        ! La valeur n'etait pas stockee

         A%a(newP(i)+i-1) = Vnew(i)
         A%ja(newP(i)+i-1) = CnewP(i)

      end do

   end subroutine insert_values

!-----------------------------------------------------------------------------------------------!
end module sparse_matrix

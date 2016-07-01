!----------------------------------------------------------------------------------------------------------
! MODULE: math
!
!> @author JL Tailhan
!> @author J. Goncalvez (version 1.0 - 2010)
!
!> @brief
!> Ensemble de fonctions et procédures mathématiques générales
!----------------------------------------------------------------------------------------------------------
module math

    implicit none

    interface find_vec
        module procedure find_int_vec, find_real_vec
    end interface
    interface find_num
        module procedure find_int_num, find_real_num
    end interface
    interface sort
        module procedure sort_int, sort_real
    end interface
    interface min_vec
        module procedure min_int_vec, min_rea_vec
    end interface

contains


!------------------------------------------------------------------------------------------------------
!> @author JL Tailhan
!> @author J. Goncalvez (version 1.0 - 2010)
!
!> @brief
!> Routine de calcul du produit membre a membre de deux vecteurs
!
!> @details
!> #### DESCRIPTION:
!>  A faire
!------------------------------------------------------------------------------------------------------
function vecpro(vec1,vec2)

    !**************************************************************************************
    !* DEFINITIONS
    !**************************************************************************************

    implicit none

    !--------------------------------------------------------------------------------------
    !--- Variables in/out
    real*8, dimension(:),intent(in) :: vec1,vec2

    !--------------------------------------------------------------------------------------
    !--- Variables locales
    real*8, dimension(size(vec1),size(vec2)) :: vecpro

    !**************************************************************************************
    !* CORPS DE PROCEDURE
    !**************************************************************************************

    if ((size(vec1) /= size(vec2)) .or. (size(vec1) > 3)) then
        print*,' Dimension de deux sont incompatibles'
        stop

    else
        if (size(vec1)==2) then
            vecpro(1,1) = vec1(1)*vec2(1)
            vecpro(1,2) = vec1(1)*vec2(2)
            vecpro(2,1) = vec1(2)*vec2(1)
            vecpro(2,2) = vec1(2)*vec2(2)
        elseif (size(vec1)==3) then
            vecpro(1,1) = vec1(1)*vec2(1)
            vecpro(1,2) = vec1(1)*vec2(2)
            vecpro(1,3) = vec1(1)*vec2(3)
            vecpro(2,1) = vec1(2)*vec2(1)
            vecpro(2,2) = vec1(2)*vec2(2)
            vecpro(2,3) = vec1(2)*vec2(3)
            vecpro(3,1) = vec1(3)*vec2(1)
            vecpro(3,2) = vec1(3)*vec2(2)
            vecpro(3,3) = vec1(3)*vec2(3)
        endif
    endif

end function vecpro


!------------------------------------------------------------------------------------------------------
!> @author JL Tailhan
!> @author J. Goncalvez (version 1.0 - 2010)
!
!> @brief
!> Routine de calcul du produit vectoriel de deux vecteurs
!
!> @details
!> #### DESCRIPTION:
!>  A faire
!------------------------------------------------------------------------------------------------------
function cross(vec1,vec2)

    !**************************************************************************************
    !* DEFINITIONS
    !**************************************************************************************

    implicit none

    !--------------------------------------------------------------------------------------
    !--- Variables in/out
    real*8, dimension(:),intent(in) :: vec1,vec2

    !--------------------------------------------------------------------------------------
    !--- Variables locales
    real*8, dimension(size(vec1)) :: cross

    !**************************************************************************************
    !* CORPS DE PROCEDURE
    !**************************************************************************************

    if ((size(vec1) /= size(vec2)) .or. (size(vec1)/= 3)) then
        stop ' Dimension de deux sont incompatibles'
    else
        cross(1) = vec1(2)*vec2(3) - vec1(3)*vec2(2)
        cross(2) = vec1(3)*vec2(1) - vec1(1)*vec2(3)
        cross(3) = vec1(1)*vec2(2) - vec1(2)*vec2(1)
    endif

end function cross

!------------------------------------------------------------------------------!
!                    Pour generer un nombre aleatoire ( 0 ->1 )                !
!------------------------------------------------------------------------------!

!   function random(N)

!        implicit none

!        integer, intent(in) :: N
!        real*8 :: random(N), rand
!        integer :: i, j, values(8), N0 = 1, N1 = 10000

!        call date_and_time(VALUES=values)
!        i = rand(values(1)+values(2)+values(3)+values(4)+values(5)+values(6)+values(7)+values(8))

!        do j = 1, N
!             random(j) = (real(rand(0) * (N1 + 1 - N0)) + N0)/N1
!        enddo

!   end function random


!------------------------------------------------------------------------------------------------------
!> @author JL Tailhan
!> @author J. Goncalvez (version 1.0 - 2010)
!
!> @brief
!> Generateur de nombre aleatoire a gde periode (part I)
!
!> @details
!> #### DESCRIPTION:
!>  A faire
!------------------------------------------------------------------------------------------------------
function randgp0(x,b)

    !**************************************************************************************
    !* DEFINITIONS
    !**************************************************************************************

    implicit none

    !--------------------------------------------------------------------------------------
    !--- Variables in/out
    integer, intent(in) :: b
    real*8, dimension(:), intent(in) :: x

    !--------------------------------------------------------------------------------------
    !--- Variables locales
    real*8, dimension(size(x)+1) :: randgp0

    integer :: n
    real*8 :: c, delta

    !**************************************************************************************
    !* CORPS DE PROCEDURE
    !**************************************************************************************

    n=size(x)
    if (n<24) stop 'math - randgp0 : mauvais dimentionnement de x'

    c=0.d0
    randgp0(1:23)=x(2:24)
    delta=x(15)-x(1)-c
    if (delta>=0) then
        randgp0(24)=delta
        randgp0(25)=0.d0
    else
        randgp0(24)=delta+b
        randgp0(25)=1.d0
    endif

end function randgp0


!------------------------------------------------------------------------------------------------------
!> @author JL Tailhan
!> @author J. Goncalvez (version 1.0 - 2010)
!
!> @brief
!> Generateur de nombre aleatoire a gde periode (part II)
!
!> @details
!> #### DESCRIPTION:
!>  A faire
!------------------------------------------------------------------------------------------------------
function randgp(n,seed)

    !**************************************************************************************
    !* DEFINITIONS
    !**************************************************************************************

    implicit none

    !--------------------------------------------------------------------------------------
    !--- Variables in/out
    integer, intent(in)    :: n
    integer, dimension(:),intent(in), optional :: seed

    !--------------------------------------------------------------------------------------
    !--- Variables locales
    integer :: b, i, values(8)
    integer :: germe(12) !, germe(12) pour gfortran 4.4 et later
    real*8 :: randgp(n), x(24), y(25)

    !**************************************************************************************
    !* CORPS DE PROCEDURE
    !**************************************************************************************

    !--------------------------------------------------------------------------------------
    !--- Initialisation du generateur
    germe = 1
    if (present(seed)) then
        if (size(seed)<1) stop 'math.f - randgp :: erreur de dimensionnement de seed'

        call date_and_time(VALUES=values)
        germe(1) = values(2)*values(8)+seed(1)
        germe(2) = values(3)*values(7)*seed(1)
        germe(3) = values(7)*values(8)/seed(1)
        germe(4) = germe(1)+germe(2)
        germe(5) = germe(2)+germe(3)
        germe(6) = germe(3)+germe(1)
        germe(7) = values(7)
        germe(8) = values(8)
        germe(7) = germe(4)+germe(5)
        germe(8) = germe(5)+germe(6)
        germe(9) = germe(6)+germe(4)
        germe(10) = germe(7)+germe(8)*seed(1)
        germe(11) = values(7)
        germe(12) = values(8)
    endif

    call random_seed(put=germe)

    !--------------------------------------------------------------------------------------
    !--- Activation du generateur
    b=2**24
    call random_number(x)
    x=(b-1)*x
    randgp = 0.d0
    do i=1,n
        randgp(i)=x(1)/real(b)
        y=randgp0(x,b)
        x=y(1:24)
    enddo

end function randgp


!------------------------------------------------------------------------------------------------------
!> @author JL Tailhan
!> @author J. Goncalvez (version 1.0 - 2010)
!
!> @brief
!> Determination de la position (si elle existe) de n dans un vecteur - ENTIERS
!
!> @details
!> #### DESCRIPTION:
!>  A faire
!------------------------------------------------------------------------------------------------------
function find_int_vec(a,n)

    !**************************************************************************************
    !* DEFINITIONS
    !**************************************************************************************

    implicit none

    !--------------------------------------------------------------------------------------
    !--- Variables in/out
    integer, dimension(:), intent(in) :: a
    integer, intent(in) :: n

    !--------------------------------------------------------------------------------------
    !--- Variables locales
    integer, dimension(count(a==n)) ::find_int_vec
    integer :: i,k

    !**************************************************************************************
    !* CORPS DE PROCEDURE
    !**************************************************************************************

    find_int_vec = 0

    k = 1
    do i =1, size(a)
        if(a(i)==n) then
            find_int_vec(k)=i
            k=k+1
        endif
    enddo

end function find_int_vec


!------------------------------------------------------------------------------------------------------
!> @author JL Tailhan
!> @author J. Goncalvez (version 1.0 - 2010)
!
!> @brief
!> Determination de la position (si elle existe) de n dans un vecteur - REELS
!
!> @details
!> #### DESCRIPTION:
!>  A faire
!------------------------------------------------------------------------------------------------------
function find_real_vec(a,n)

    !**************************************************************************************
    !* DEFINITIONS
    !**************************************************************************************

    implicit none

    !--------------------------------------------------------------------------------------
    !--- Variables in/out
    real*8, dimension(:), intent(in) :: a
    real*8, intent(in) :: n

    !--------------------------------------------------------------------------------------
    !--- Variables locales
    integer, dimension(count(a==n)) ::find_real_vec
    integer :: i,k

    !**************************************************************************************
    !* CORPS DE PROCEDURE
    !**************************************************************************************

    find_real_vec = 0

    k=1
    do i = 1, size(a)
        if(a(i)==n) then
            find_real_vec(k)=i
            k=k+1
        endif
    enddo

end function find_real_vec


!------------------------------------------------------------------------------------------------------
!> @author JL Tailhan
!> @author J. Goncalvez (version 1.0 - 2010)
!> @author Anna Maria   (version 1.0 - 2011)
!
!> @brief
!> Determination des indices des composantes non nulles d'un vecteur du type int
!
!> @details
!> #### DESCRIPTION:
!>  A faire
!------------------------------------------------------------------------------------------------------
function find_int_non_vec(a,n)

    !**************************************************************************************
    !* DEFINITIONS
    !**************************************************************************************

    implicit none

    !--------------------------------------------------------------------------------------
    !--- Variables in/out
    integer, dimension(:), intent(in) :: a
    integer, intent(in) :: n

    !--------------------------------------------------------------------------------------
    !--- Variables locales
    integer, dimension(count(a/=n)) ::find_int_non_vec
    integer :: i,k

    !**************************************************************************************
    !* CORPS DE PROCEDURE
    !**************************************************************************************

    find_int_non_vec = 0

    k=1
    do i =1, size(a)
        if(a(i)==n) then
             find_int_non_vec(k)=i
             k=k+1
        endif
    enddo

end function find_int_non_vec


!------------------------------------------------------------------------------------------------------
!> @author JL Tailhan
!> @author J. Goncalvez (version 1.0 - 2010)
!
!> @brief
!> Determination des indices des composantes non nulles d'un vecteur du type booleen
!
!> @details
!> #### DESCRIPTION:
!>  A faire
!------------------------------------------------------------------------------------------------------
function find_log_vec(a)

    !**************************************************************************************
    !* DEFINITIONS
    !**************************************************************************************

    implicit none

    !--------------------------------------------------------------------------------------
    !--- Variables in/out
    logical, dimension(:), intent(in) :: a

    !--------------------------------------------------------------------------------------
    !--- Variables locales
    integer, dimension(count(a .eqv. .true.)) :: find_log_vec
    integer :: i,k

    !**************************************************************************************
    !* CORPS DE PROCEDURE
    !**************************************************************************************

    find_log_vec = 0
    k=1
    do i = 1, size(a)
        if (a(i)) then
           find_log_vec(k)=i
            k=k+1
        endif
    enddo

end function find_log_vec


!------------------------------------------------------------------------------------------------------
!> @author JL Tailhan
!> @author J. Goncalvez (version 1.0 - 2010)
!
!> @brief
!> Determination des composantes non nulles d'un vecteur du type integer
!
!> @details
!> #### DESCRIPTION:
!>  A faire
!------------------------------------------------------------------------------------------------------
function find_int_num(a,n)

    !**************************************************************************************
    !* DEFINITIONS
    !**************************************************************************************

    implicit none

    !--------------------------------------------------------------------------------------
    !--- Variables in/out
    integer, dimension(:), intent(in) :: a
    integer, intent(in) :: n

    !--------------------------------------------------------------------------------------
    !--- Variables locales
    integer :: find_int_num
    integer :: i

    !**************************************************************************************
    !* CORPS DE PROCEDURE
    !**************************************************************************************

    find_int_num = 0

    do i = 1, size(a)
        if (a(i)==n) find_int_num=i
    enddo

   end function find_int_num


!------------------------------------------------------------------------------------------------------
!> @author JL Tailhan
!> @author J. Goncalvez (version 1.0 - 2010)
!
!> @brief
!> Determination des composantes non nulles d'un vecteur du type reel
!
!> @details
!> #### DESCRIPTION:
!>  A faire
!------------------------------------------------------------------------------------------------------
function find_real_num(a,n)

    !**************************************************************************************
    !* DEFINITIONS
    !**************************************************************************************

    implicit none

    !--------------------------------------------------------------------------------------
    !--- Variables in/out
    real*8, dimension(:), intent(in) :: a
    real*8, intent(in) :: n

    !--------------------------------------------------------------------------------------
    !--- Variables locales
    integer :: find_real_num
    integer :: i

    !**************************************************************************************
    !* CORPS DE PROCEDURE
    !**************************************************************************************

    find_real_num = 0

    do i = 1, size(a)
       if(a(i)==n) find_real_num=i
    enddo

end function find_real_num


!------------------------------------------------------------------------------------------------------
!> @author JL Tailhan
!> @author J. Goncalvez (version 1.0 - 2010)
!
!> @brief
!> Fonction unique (avec reordonnancement)
!
!> @details
!> #### DESCRIPTION:
!>  A faire
!------------------------------------------------------------------------------------------------------
subroutine unique(v,k,u)

    !**************************************************************************************
    !* DEFINITIONS
    !**************************************************************************************

    implicit none

    !--------------------------------------------------------------------------------------
    !--- Variables in/out
    integer, dimension(:), intent(in) :: v
    integer, dimension(size(v)), intent(out) :: u
    integer, intent(out) :: k ! Taille de unique

    !--------------------------------------------------------------------------------------
    !--- Variables locales
    integer, dimension(size(v)) :: tmp1
    integer :: i, n

    !**************************************************************************************
    !* CORPS DE PROCEDURE
    !**************************************************************************************

    n = size(v)

    tmp1 = v(sort(v))
    u(1) = tmp1(1)
    k = 1

    do i = 1, n-1
        if (tmp1(i+1) /= tmp1(i)) then
            k = k + 1
            u(k) = tmp1(i+1)
        endif
    enddo


end subroutine unique


!------------------------------------------------------------------------------------------------------
!> @author JL Tailhan
!> @author J. Goncalvez (version 1.0 - 2010)
!
!> @brief
!> Fonction unique_r (QUE POUR LES VECTEUR REEL ORDONNEES)
!
!> @details
!> #### DESCRIPTION:
!>  A faire
!------------------------------------------------------------------------------------------------------
subroutine unique_r(v,k,u)

    !**************************************************************************************
    !* DEFINITIONS
    !**************************************************************************************

    implicit none

    !--------------------------------------------------------------------------------------
    !--- Variables in/out
    real*8 , dimension(:), intent(in) :: v
    real*8 , dimension(:), allocatable, intent(out) :: u
    integer, intent(out) :: k ! Taille de unique_r

    !--------------------------------------------------------------------------------------
    !--- Variables locales
    real*8, dimension(size(v)) :: tmp1
    integer :: i, n

    !**************************************************************************************
    !* CORPS DE PROCEDURE
    !**************************************************************************************

    n = size(v)

    tmp1(1) = v(1)
    k = 1

    do i = 1, n-1
        if (v(i+1) /= v(i)) then
            k = k + 1
            tmp1(k) = v(i+1)
        endif
    enddo
      
    allocate(u(k))
      
    u = tmp1(1:k)

end subroutine unique_r


!------------------------------------------------------------------------------------------------------
!> @author JL Tailhan
!> @author J. Goncalvez (version 1.0 - 2010)
!
!> @brief
!> Fonction de reordonnancement des valeurs d'un vecteur d'entier
!
!> @details
!> #### DESCRIPTION:
!>  A faire
!------------------------------------------------------------------------------------------------------
function sort_int(v) result(p)

    !**************************************************************************************
    !* DEFINITIONS
    !**************************************************************************************

    implicit none

    !--------------------------------------------------------------------------------------
    !--- Variables in/out
    integer, dimension(:), intent(in) :: v
    integer, dimension(size(v)) :: p

    !--------------------------------------------------------------------------------------
    !--- Variables locales
    integer :: i, j, n, minv, posj, aux

    !**************************************************************************************
    !* CORPS DE PROCEDURE
    !**************************************************************************************

    n = size(v)

    p = (/(i, i = 1, n )/)

    do i = 1, n
        minv = v(p(i))
        posj = i
        do j = i + 1, n
            if(v(p(j)) < minv) then
            minv = v(p(j))
            posj = j
          endif
        enddo

        ! Inversion des deux valeurs
        aux = p(i)
        p(i) = p(posj)
        p(posj) = aux

    enddo

end function sort_int


!------------------------------------------------------------------------------------------------------
!> @author JL Tailhan
!> @author J. Goncalvez (version 1.0 - 2010)
!
!> @brief
!> Fonction de reordonnancement des valeurs d'un vecteur de reels
!
!> @details
!> #### DESCRIPTION:
!>  A faire
!------------------------------------------------------------------------------------------------------
function sort_real(v) result(p)

    !**************************************************************************************
    !* DEFINITIONS
    !**************************************************************************************

    implicit none

    !--------------------------------------------------------------------------------------
    !--- Variables in/out
    real*8, dimension(:), intent(in) :: v
    integer, dimension(size(v)) :: p

    !--------------------------------------------------------------------------------------
    !--- Variables locales
    real*8 :: minv
    integer :: i, j, n, posj, aux

    !**************************************************************************************
    !* CORPS DE PROCEDURE
    !**************************************************************************************

    n = size(v)

    p = (/(i, i = 1, n )/)

    do i = 1, n
        minv = v(p(i))
        posj = i
        do j = i + 1, n
            if(v(p(j)) < minv) then
                minv = v(p(j))
                posj = j
            endif
        enddo

        ! Inversion des deux valeurs
        aux = p(i)
        p(i) = p(posj)
        p(posj) = aux

    enddo

end function sort_real


!------------------------------------------------------------------------------------------------------
!> @author JL Tailhan
!> @author J. Goncalvez (version 1.0 - 2010)
!
!> @brief
!> Sorts an array RA, in ascending order by the Heapsort method
!
!> @details
!> #### DESCRIPTION:
!>  INPUTS:
!>     RA      1D array to be sorted
!>  OUTPUT:
!>     RA    1D array sorted in ascending order
!>  NOTE: The Heapsort method is a N Log2 N routine, and can be used for very large arrays.
!------------------------------------------------------------------------------------------------------
subroutine HPSORT(RA)

    !**************************************************************************************
    !* DEFINITIONS
    !**************************************************************************************

    implicit none

    !--------------------------------------------------------------------------------------
    !--- Variables in/out
    real*8, dimension(:), intent(inout) :: RA

    !--------------------------------------------------------------------------------------
    !--- Variables locales
    real*8  :: RRA
    integer :: i,j,l,n,ir

    !**************************************************************************************
    !* CORPS DE PROCEDURE
    !**************************************************************************************

    N = SIZE(RA)
    L = N/2+1
    IR= N
    !The index L will be decremented from its initial value during the
    !"hiring" (heap creation) phase. Once it reaches 1, the index IR 
    !will be decremented from its initial value down to 1 during the
    !"retirement-and-promotion" (heap selection) phase.
    10  continue
        if(L > 1)then
            L=L-1
            RRA=RA(L)
        else
            RRA=RA(IR)
            RA(IR)=RA(1)
            IR=IR-1
            if(IR.eq.1)then
                RA(1)=RRA
                return
            endif
        endif
        I=L
        J=L+L
    20  if(J.le.IR)then
            if(J < IR)then
                if(RA(J) < RA(J+1))  J=J+1
            endif
            if(RRA < RA(J))then
                RA(I)=RA(J)
                I=J; J=J+J
            else
                J=IR+1
            endif
            goto 20
        endif
        RA(I)=RRA
        goto 10

end subroutine HPSORT


!================================== Manipulation sur les matrices ====================================!


!------------------------------------------------------------------------------------------------------
!> @author JL Tailhan
!> @author J. Goncalvez (version 1.0 - 2010)
!
!> @brief
!> Determine les positions dans une matrice d'une valeur n donnee
!
!> @details
!> #### DESCRIPTION:
!>  A faire
!------------------------------------------------------------------------------------------------------
subroutine find_mat(A,n,indi)

    !**************************************************************************************
    !* DEFINITIONS
    !**************************************************************************************

    implicit none

    !--------------------------------------------------------------------------------------
    !--- Variables in/out
    integer, dimension(:,:), intent(in) :: A
    integer, intent(in) :: n
    integer, dimension(count(A==n)), intent(inout) :: indi

    !--------------------------------------------------------------------------------------
    !--- Variables locales
    integer :: i,j,k

    !**************************************************************************************
    !* CORPS DE PROCEDURE
    !**************************************************************************************

    indi = 0

    k=1
    do i = 1, size(A,1)
        do j = 1,size(A,2)
            if (A(i,j)==n) then
                indi(k) = i
                k=k+1
                exit
            endif
        enddo
    enddo

end subroutine find_mat


!================================== Comparaison de vecteurs et matrices ==============================!


!------------------------------------------------------------------------------------------------------
!> @author JL Tailhan
!> @author J. Goncalvez (version 1.0 - 2010)
!
!> @brief
!> Determine les valeurs communes de deux vecteurs d'entiers
!
!> @details
!> #### DESCRIPTION:
!>  A faire
!------------------------------------------------------------------------------------------------------
function intersect_mat(vec1, vec2)

    !**************************************************************************************
    !* DEFINITIONS
    !**************************************************************************************

    implicit none

    !--------------------------------------------------------------------------------------
    !--- Variables in/out
    integer, dimension(:), intent(in) :: vec1, vec2

    !--------------------------------------------------------------------------------------
    !--- Variables locales
    integer, dimension(min(size(vec1),size(vec2))):: intersect_mat
    integer :: i,j,k

    !**************************************************************************************
    !* CORPS DE PROCEDURE
    !**************************************************************************************

    k = 0
    intersect_mat = 0

    do i = 1, size(vec1)
        do j = 1, size(vec2)
            if (vec1(i)==vec2(j)) then
                k = k + 1
                intersect_mat(k) = vec1(i)
            endif
        enddo
    enddo

end function intersect_mat


!------------------------------------------------------------------------------------------------------
!> @author JL Tailhan
!> @author J. Goncalvez (version 1.0 - 2010)
!
!> @brief
!> Determine la valeur commune de deux vecteurs d'entiers
!
!> @details
!> #### DESCRIPTION:
!>  A faire
!------------------------------------------------------------------------------------------------------
function intersect_num(vec1, vec2)

    !**************************************************************************************
    !* DEFINITIONS
    !**************************************************************************************

    implicit none

    !--------------------------------------------------------------------------------------
    !--- Variables in/out
    integer, dimension(:), intent(in) :: vec1, vec2

    !--------------------------------------------------------------------------------------
    !--- Variables locales
    integer :: intersect_num
    integer :: i,j

    !**************************************************************************************
    !* CORPS DE PROCEDURE
    !**************************************************************************************

    intersect_num = 0

    do i = 1, size(vec1)
        do j = 1, size(vec2)
            if(vec1(i)==vec2(j)) then
                intersect_num = vec1(i)
            endif
        enddo
    enddo

end function intersect_num


!------------------------------------------------------------------------------------------------------
!> @author JL Tailhan
!> @author J. Goncalvez (version 1.0 - 2010)
!
!> @brief
!> Determine les lignes communes dans un tableau pour des valeurs d'entiers
!
!> @details
!> #### DESCRIPTION:
!>  A faire
!------------------------------------------------------------------------------------------------------
function intersection(A,vec)

    !**************************************************************************************
    !* DEFINITIONS
    !**************************************************************************************

    implicit none

    !--------------------------------------------------------------------------------------
    !--- Variables in/out
    integer, dimension(:,:), intent(in) :: A
    integer, dimension(:), intent(in) :: vec

    !--------------------------------------------------------------------------------------
    !--- Variables locales
    integer, dimension(2) :: intersection
    integer :: i, j, k, l, m, n, t

    !**************************************************************************************
    !* CORPS DE PROCEDURE
    !**************************************************************************************

    t = 1
    intersection = 0

    do i = 1, size(A,1)
        k = 0 ; l = 0 ; m = 0 ; n = 0
        do j = 1, size(A,2)
            if (A(i,j)==vec(1)) k = 1
            if (A(i,j)==vec(2)) l = 1

            if (size(vec)==3) then
                if (A(i,j)==vec(3)) m = 1
            endif

            if (size(vec)==4) then
                if (A(i,j)==vec(3)) m = 1
                if (A(i,j)==vec(4)) n = 1
            endif
        enddo

        if ((k==1 .and. l==1 .and. size(vec)==2).or.(k==1 .and. l==1 .and. m==1 .and. size(vec)==3) &
        &    .or. (k==1 .and. l==1 .and. m==1 .and. n==1 .and. size(vec)==4)) then
            intersection(t) = i
            t = t + 1
        endif
    enddo

end function intersection


!------------------------------------------------------------------------------------------------------
!> @author JL Tailhan
!> @author J. Goncalvez (version 1.0 - 2010)
!
!> @brief
!> Determine les valeurs differentes de deux vecteurs d'entiers
!
!> @details
!> #### DESCRIPTION:
!>  A faire
!------------------------------------------------------------------------------------------------------
function setxor_mat(vec1, vec2)

    !**************************************************************************************
    !* DEFINITIONS
    !**************************************************************************************

    implicit none

    !--------------------------------------------------------------------------------------
    !--- Variables in/out
    integer, dimension(:), intent(in) :: vec1, vec2

    !--------------------------------------------------------------------------------------
    !--- Variables locales
    integer, dimension(size(vec1)+size(vec2)):: setxor_mat
    integer :: i,j,k

    !**************************************************************************************
    !* CORPS DE PROCEDURE
    !**************************************************************************************

    k = 0
    setxor_mat = 0

    do i = 1, size(vec1)
        do j = 1, size(vec2)
            if(vec1(i)/=vec2(j)) then
                k = k + 1
                setxor_mat(k) = vec1(i)
            endif
        enddo
    enddo

end function setxor_mat


!------------------------------------------------------------------------------------------------------
!> @author JL Tailhan
!> @author J. Goncalvez (version 1.0 - 2010)
!
!> @brief
!> Determine la difference de deux vecteurs d'entiers
!
!> @details
!> #### DESCRIPTION:
!>  A faire
!------------------------------------------------------------------------------------------------------
function setxor_num(vec1, vec2)

    !**************************************************************************************
    !* DEFINITIONS
    !**************************************************************************************

    implicit none

    !--------------------------------------------------------------------------------------
    !--- Variables in/out
    integer, dimension(:), intent(in) :: vec1, vec2

    !--------------------------------------------------------------------------------------
    !--- Variables locales
    integer :: setxor_num
    integer :: i,j

    !**************************************************************************************
    !* CORPS DE PROCEDURE
    !**************************************************************************************

    setxor_num = 0

    do i = 1, size(vec1)
        do j = 1, size(vec2)
            if ((vec1(i) /=0) .and. (vec2(j) /=0) .and. (vec1(i)/=vec2(j))) then
                setxor_num = vec1(i)
            endif
        enddo
    enddo

end function setxor_num


!========================= Resolution de systeme lineaire ==================!

!*******************************************************
!*    LU decomposition routines used by test_lu.f90    *
!*                                                     *
!*                 F90 version by J-P Moreau, Paris    *
!* --------------------------------------------------- *
!* Reference:                                          *
!*                                                     *
!* "Numerical Recipes By W.H. Press, B. P. Flannery,   *
!*  S.A. Teukolsky and W.T. Vetterling, Cambridge      *
!*  University Press, 1986" [BIBLI 08].                *
!*                                                     *
!*******************************************************

!--------------------------------------------------------------------------------

!  ***************************************************************
!  * Given an N x N matrix A, this routine replaces it by the LU *
!  * decomposition of a rowwise permutation of itself. A and N   *
!  * are input. INDX is an output vector which records the row   *
!  * permutation effected by the partial pivoting; D is output   *
!  * as -1 or 1, depending on whether the number of row inter-   *
!  * changes was even or odd, respectively. This routine is used *
!  * in combination with LUBKSB to solve linear equations or to  *
!  * invert a matrix. Return code is 1, if matrix is singular.   *
!  ***************************************************************
 Subroutine LUDCMP(A,N,INDX,D,CODE)

 real*8, intent(inout) :: A(N,N)
 integer, intent(inout) :: CODE, D, INDX(N)
 integer, intent(in) :: N

 integer, parameter :: NMAX=300
 real*8, parameter  :: TINY=1.5D-16

 real*8 :: AMAX,DUM, SUM
 real*8 :: VV(N)
 integer :: I, J, K, IMAX
 D=1; CODE=0

 ! Ajout DEB JLT 2015-10-14
 IMAX = 0
 ! Ajout FIN JLT 2015-10-14

 DO I=1,N
   AMAX=0.d0
   DO J=1,N
     IF (DABS(A(I,J)).GT.AMAX) AMAX=DABS(A(I,J))
   enddo ! j loop

   IF(AMAX.LT.TINY) THEN
     CODE = 1
     RETURN
   endif

   VV(I) = 1.d0 / AMAX

 enddo ! i loop

 DO J=1,N
   DO I=1,J-1
     SUM = A(I,J)
     DO K=1,I-1
       SUM = SUM - A(I,K)*A(K,J)
     enddo ! k loop
     A(I,J) = SUM
   enddo ! i loop
   AMAX = 0.d0
   DO I=J,N
     SUM = A(I,J)
     DO K=1,J-1
       SUM = SUM - A(I,K)*A(K,J)
     enddo ! k loop
     A(I,J) = SUM
     DUM = VV(I)*DABS(SUM)
     IF(DUM.GE.AMAX) THEN
       IMAX = I
       AMAX = DUM
     endif
   enddo ! i loop

   IF(J.NE.IMAX) THEN
     DO K=1,N
       DUM = A(IMAX,K)
       A(IMAX,K) = A(J,K)
       A(J,K) = DUM
     enddo ! k loop
     D = -D
     VV(IMAX) = VV(J)
   endif

   INDX(J) = IMAX
   IF(DABS(A(J,J)) < TINY) A(J,J) = TINY

   IF(J.NE.N) THEN
     DUM = 1.d0 / A(J,J)
     DO I=J+1,N
       A(I,J) = A(I,J)*DUM
     enddo ! i loop
   endif
 enddo ! j loop

 RETURN
 END subroutine LUDCMP

!-----------------------------------------------------------------------------

!  ******************************************************************
!  * Solves the set of N linear equations A . X = B.  Here A is     *
!  * input, not as the matrix A but rather as its LU decomposition, *
!  * determined by the routine LUDCMP. INDX is input as the permuta-*
!  * tion vector returned by LUDCMP. B is input as the right-hand   *
!  * side vector B, and returns with the solution vector X. A, N and*
!  * INDX are not modified by this routine and can be used for suc- *
!  * cessive calls with different right-hand sides. This routine is *
!  * also efficient for plain matrix inversion.                     *
!  ******************************************************************

 Subroutine LUBKSB(A,N,INDX,B)
 integer :: N
 REAL*8, intent(in) :: A(N,N)
 real*8,intent(inout) :: B(N)
 INTEGER, intent(in) :: INDX(N)
 integer :: I,J,II, LL
 real*8 :: sum

 II = 0

 DO I=1,N
   LL = INDX(I)
   SUM = B(LL)
   B(LL) = B(I)
   IF(II.NE.0) THEN
     DO J=II,I-1
       SUM = SUM - A(I,J)*B(J)
     enddo ! j loop
   ELSE IF(SUM.NE.0.d0) THEN
     II = I
   endif
   B(I) = SUM
 enddo ! i loop

 DO I=N,1,-1
   SUM = B(I)
   IF(I < N) THEN
     DO J=I+1,N
       SUM = SUM - A(I,J)*B(J)
     enddo ! j loop
   endif
   B(I) = SUM / A(I,I)
 enddo ! i loop

 RETURN
 END subroutine LUBKSB


!------------------------------------------------------------------------------------------------------
!> @author JL Tailhan
!> @author J. Goncalvez (version 1.0 - 2010)
!
!> @brief
!> Calcul du determinant d'une matrice A de dimension NxN
!
!> @details
!> #### DESCRIPTION:
!>  A faire
!------------------------------------------------------------------------------------------------------
function DET(A)

 real*8, dimension(:,:), intent(in) :: A
 real*8 :: DET
 integer :: N
 integer ID, INDX(size(A,1)), icode
 real*8, dimension(size(A,1),size(A,2)):: A1
 integer :: J

 if(size(A,1)/=size(A,2)) then
        print*, 'det : la matrice doit-être carre'
        stop
 endif

 A1 = A

 N=size(A,1)
 call LUDCMP(A1,N,INDX,ID,ICODE)

 DET = DFLOAT(ID)

 do J=1, N
   DET=DET*A1(J,J)
 enddo

end function DET


!------------------------------------------------------------------------------------------------------
!> @author JL Tailhan
!> @author J. Goncalvez (version 1.0 - 2010)
!
!> @brief
!> Inverse de matrice
!
!> @details
!> #### DESCRIPTION:
!>  A faire
!------------------------------------------------------------------------------------------------------
function INV(A)

 real*8, dimension(:,:), intent(in) :: A
 real*8, dimension(size(A,1),size(A,1)) :: INV
 integer :: M, N
 integer :: INFO
 real*8, dimension(size(A,1)):: IPIV
 real*8, dimension(size(A,1)) :: WORK

  N=size(A,1); M=size(A,2);

  if(N/=M) then
          print*, 'Inv : la matrice doit etre carre'
          stop
  endif

  INV = A

  call DGETRF(N,N,INV,N,IPIV,INFO)

  call DGETRI(N,INV,N,IPIV,WORK,N,INFO)

end function INV


!------------------------------------------------------------------------------------------------------
!> @author JL Tailhan
!> @author J. Goncalvez (version 1.0 - 2010)
!
!> @brief
!> Resolution Ax =B (format dense) (A et B non modifie par l'appel) 
!
!> @details
!> #### DESCRIPTION:
!>  A faire
!------------------------------------------------------------------------------------------------------
function Solve_dns(A,B) result (B1)

 real*8, dimension(:,:), intent(in) :: A
 real*8, dimension(:), intent(in) :: B
 integer :: N
 integer :: icode
 real*8, dimension(size(A,1),size(A,2)):: A1
 real*8, dimension(size(A,1)):: B1
 real*8, dimension(size(A,1)):: IPIV

 A1 = A; B1 = B

 N = size(A,1)

 call DGESV(N, 1, A1, N, IPIV, B1, N, icode)

 end function Solve_dns


!------------------------------------------------------------------------------------------------------
!> @author JL Tailhan
!> @author J. Goncalvez (version 1.0 - 2010)
!
!> @brief
!> Calcul des valeurs propres et des vecteurs propres d'une matrice reelle symetrique
!
!> @details
!> #### DESCRIPTION:
!>  A faire
!------------------------------------------------------------------------------------------------------
!=============================================================================
!---------- Calcul des valeurs propres et des vecteurs propres  --------------
!----------        d'une matrice reelle symetrique             ---------------
!---------- ATTENTION : valable uniquement pour les matrices SYMETRIQUES -----
!*************************************************************
!* This subroutine computes all eigenvalues and eigenvectors *
!* of a real SYMMETRIC square matrix A(N,N). On output, ele- *
!* ments of A above the diagonal are destroyed. D(N) returns *
!* the eigenvalues of matrix A. V(N,N) contains, on output,  *
!* the eigenvectors of A by columns. THe normalization to    *
!* unity is made by main program before printing results.    *
!* NROT returns the number of Jacobi matrix rotations which  *
!* were required.                                            *
!* --------------------------------------------------------- *
!* Ref.:"NUMERICAL RECIPES, Cambridge University Press, 1986,*
!*       chap. 11, pages 346-348" [BIBLI 08].                *
!*************************************************************
subroutine Jacobi(A,N,D,V,NROT)
implicit none

integer, intent(in) :: N
real*8, dimension(N,N), intent(inout) :: A
real*8, dimension(N,N), intent(inout) :: V
real*8, dimension(N), intent(inout) :: D
integer, intent(out) :: NROT

real*8, dimension(:) ,allocatable :: B, Z
real*8  c,g,h,s,sm,t,tau,theta,tresh
integer :: ip,iq,i,j

allocate(B(100)) ; B = 0.d0
allocate(Z(100)) ; Z = 0.d0

  do ip=1, N    !initialize V to identity matrix
    do iq=1, N
      V(ip,iq)=0.d0
    enddo
      V(ip,ip)=1.d0
  enddo
  do ip=1, N
    B(ip)=A(ip,ip)
    D(ip)=B(ip)
    Z(ip)=0.d0
  enddo
  NROT=0
  do i=1, 50
    sm=0.d0
    do ip=1, N-1     !sum off-diagonal elements
      do iq=ip+1, N
        sm=sm+DABS(A(ip,iq))
      enddo
    enddo
    if(sm==0.d0) return  !normal return
    if(i.lt.4) then
      tresh=0.2d0*sm**2
    else
      tresh=0.d0
    endif
    do ip=1, N-1
      do iq=ip+1, N
        g=100.d0*DABS(A(ip,iq))
! after 4 sweeps, skip the rotation if the off-diagonal element is small
        if((i.gt.4).and.(DABS(D(ip))+g.eq.DABS(D(ip))) &
                .and.(DABS(D(iq))+g.eq.DABS(D(iq)))) then
                  A(ip,iq)=0.d0
        else if(DABS(A(ip,iq)).gt.tresh) then
          h=D(iq)-D(ip)
          if(DABS(h)+g.eq.DABS(h)) then
            t=A(ip,iq)/h
          else
            theta=0.5d0*h/A(ip,iq)
            t=1.d0/(DABS(theta)+DSQRT(1.d0+theta**2))
            if(theta.lt.0.d0) t=-t
          endif
          c=1.d0/DSQRT(1.d0+t**2)
          s=t*c
          tau=s/(1.d0+c)
          h=t*A(ip,iq)
          Z(ip)=Z(ip)-h
          Z(iq)=Z(iq)+h
          D(ip)=D(ip)-h
          D(iq)=D(iq)+h
          A(ip,iq)=0.d0
          do j=1, ip-1
            g=A(j,ip)
            h=A(j,iq)
            A(j,ip)=g-s*(h+g*tau)
            A(j,iq)=h+s*(g-h*tau)
          enddo
          do j=ip+1, iq-1
            g=A(ip,j)
            h=A(j,iq)
            A(ip,j)=g-s*(h+g*tau)
            A(j,iq)=h+s*(g-h*tau)
          enddo
          do j=iq+1, N
            g=A(ip,j)
            h=A(iq,j)
            A(ip,j)=g-s*(h+g*tau)
            A(iq,j)=h+s*(g-h*tau)
          enddo
          do j=1, N
            g=V(j,ip)
            h=V(j,iq)
            V(j,ip)=g-s*(h+g*tau)
            V(j,iq)=h+s*(g-h*tau)
          enddo
          NROT=NROT+1
        endif !if ((i.gt.4)...
      enddo !main iq loop
    enddo !main ip loop
    do ip=1, N
      B(ip)=B(ip)+Z(ip)
      D(ip)=B(ip)
      Z(ip)=0.d0
    enddo
  enddo !main i loop
  !pause ' 50 iterations !'
  deallocate (B,Z)
  return
END subroutine Jacobi


!------------------------------------------------------------------------------------------------------
!> @author JL Tailhan
!> @author J. Goncalvez (version 1.0 - 2010)
!
!> @brief
!> Calcul de la norme euclidienne d'un vecteur
!
!> @details
!> #### DESCRIPTION:
!>  A faire
!------------------------------------------------------------------------------------------------------
function norme(v)

    !**************************************************************************************
    !* DEFINITIONS
    !**************************************************************************************

    implicit none

    !--------------------------------------------------------------------------------------
    !--- Variables in/out
    real*8, dimension(:), intent(in) :: v

    !--------------------------------------------------------------------------------------
    !--- Variables locales
    real*8 :: norme

    !**************************************************************************************
    !* CORPS DE PROCEDURE
    !**************************************************************************************

    norme = sqrt(dot_product(v,v))

end function norme


!------------------------------------------------------------------------------------------------------
!> @author JL Tailhan
!> @author J. Goncalvez (version 1.0 - 2010)
!
!> @brief
!> Calcul du maximum d'un vecteur d'entier
!
!> @details
!> #### DESCRIPTION:
!>  A faire
!------------------------------------------------------------------------------------------------------
function max_vec(v)

    !**************************************************************************************
    !* DEFINITIONS
    !**************************************************************************************

    implicit none

    !--------------------------------------------------------------------------------------
    !--- Variables in/out
    integer, dimension(:), intent(in) :: v

    !--------------------------------------------------------------------------------------
    !--- Variables locales
    integer :: max_vec
    integer :: i

    !**************************************************************************************
    !* CORPS DE PROCEDURE
    !**************************************************************************************

    max_vec = v(1)

    do i = 2, size(v)
        max_vec = max(v(i),max_vec)
    enddo

end function max_vec


!------------------------------------------------------------------------------------------------------
!> @author JL Tailhan
!> @author J. Goncalvez (version 1.0 - 2010)
!
!> @brief
!> Calcul du minimum d'un vecteur d'entier
!
!> @details
!> #### DESCRIPTION:
!>  A faire
!------------------------------------------------------------------------------------------------------
function min_int_vec(v)

    !**************************************************************************************
    !* DEFINITIONS
    !**************************************************************************************

    implicit none

    !--------------------------------------------------------------------------------------
    !--- Variables in/out
    integer, dimension(:), intent(in) :: v

    !--------------------------------------------------------------------------------------
    !--- Variables locales
    integer :: min_int_vec
    integer :: i

    !**************************************************************************************
    !* CORPS DE PROCEDURE
    !**************************************************************************************

    min_int_vec = v(1)

    do i = 2, size(v)
        min_int_vec = min(v(i),min_int_vec)
    enddo

end function min_int_vec


!------------------------------------------------------------------------------------------------------
!> @author JL Tailhan
!> @author J. Goncalvez (version 1.0 - 2010)
!
!> @brief
!> Calcul du minimum d'un vecteur de reels
!
!> @details
!> #### DESCRIPTION:
!>  A faire
!------------------------------------------------------------------------------------------------------
function min_rea_vec(v)

    !**************************************************************************************
    !* DEFINITIONS
    !**************************************************************************************

    implicit none

    !--------------------------------------------------------------------------------------
    !--- Variables in/out
    real*8, dimension(:), intent(in) :: v

    !--------------------------------------------------------------------------------------
    !--- Variables locales
    real*8 :: min_rea_vec
    integer :: i

    !**************************************************************************************
    !* CORPS DE PROCEDURE
    !**************************************************************************************

    min_rea_vec = v(1)

    do i = 2, size(v)
        min_rea_vec = min(v(i),min_rea_vec)
    enddo

end function min_rea_vec


!------------------------------------------------------------------------------------------------------
!> @author JL Tailhan
!> @author J. Goncalvez (version 1.0 - 2010)
!
!> @brief
!> Calcul du minimum minimum en valeur absolue et non nul d'un vecteur de reels
!
!> @details
!> #### DESCRIPTION:
!>  A faire
!------------------------------------------------------------------------------------------------------
function min_abs_rea_vec(v)

    !**************************************************************************************
    !* DEFINITIONS
    !**************************************************************************************

    implicit none

    !--------------------------------------------------------------------------------------
    !--- Variables in/out
    real*8, dimension(:), intent(in) :: v

    !--------------------------------------------------------------------------------------
    !--- Variables locales
    real*8 :: min_abs_rea_vec
    logical, dimension(size(v)) :: t
    integer :: i,k,pos

    !**************************************************************************************
    !* CORPS DE PROCEDURE
    !**************************************************************************************

    t = (v/=0.d0)
    min_abs_rea_vec = abs(v(1))

    k = 0
    pos = 0
    do i = 1, size(v)
        if (t(i)) then
            k = k+1
            if (k==1) then
                min_abs_rea_vec = abs(v(i))
            else
                min_abs_rea_vec = min(abs(v(i)),min_abs_rea_vec)
                pos=i
            endif
        endif
    enddo
    min_abs_rea_vec=v(pos)

end function min_abs_rea_vec
    
!---------------------------------------------------------------------------
!> @author CN
!> @brief Matrice de rotation d'un vecteur autour d'un axe qq (modules fiss et fissBA)
!
!> @details
!> ### DESCRIPTION:
!> En 2D : Matrice de rotation usuelle (angle = theta)
!> En 3D : Matrice de rotation autour d'un axe (a, b, c)
!>
!>--------------------------------------------------------------------------
function rotmat(theta, dime, axe) result (M)

        implicit none

        real*8, dimension(dime), optional, intent(in) :: axe
        real*8, intent(in) :: theta
        integer, intent(in) :: dime

        real*8, dimension(dime,dime) :: M
        real*8 :: a, b, c

!********************************************************!
        M = 0.d0
        if (dime == 2) then
            M = reshape((/ cos(theta), sin(theta), &
                           -sin(theta), cos(theta) /), shape(M))
        else
            a = axe(1)
            b = axe(2)
            c = axe(3)
            M = reshape((/ a**2+(1-a**2)*cos(theta), a*b*(1-cos(theta))+c*sin(theta), &
                           a*c*(1-cos(theta))-b*sin(theta), a*b*(1-cos(theta))-c*sin(theta), &
                           b**2+(1-b**2)*cos(theta), b*c*(1-cos(theta))+a*sin(theta), &
                           a*c*(1-cos(theta))+b*sin(theta), b*c*(1-cos(theta))-a*sin(theta), & 
                           c**2+(1-c**2)*cos(theta) /), shape(M))
        end if


end function rotmat

!-------------------------------------------------------------------------------!
!         Pour filtrer un n. reel                                                !
!-------------------------------------------------------------------------------!

!    subroutine filtre_reel(x,d)

!        implicit none

!        real*8, intent(inout) :: x
!        integer, intent(in), optional :: d

!        real*8 :: pippo !, m
!        integer :: y, d1, dd = 5

!        y = floor(log10(x))

!        print*, 'y', y

!        !m = x*10.**(-y)

!        !print*, m

!        d1 = dd - y
!        if(present(d)) d1 = d - y

!        print*, d1

!        x = floor((10.**d1)*x) / 10.**(d1)

!        print*, '10.**(-d1)', 10.**(-d1)
!        print*, 'floor((10.**d1)*x) *1', floor((10.**d1)*x) * 1.

!        !x = m*(10.**y)
!        print*, 'x', x

!    end subroutine filtre_reel

!------------------------------------------------------------------------------------------------------

end module math

!----------------------------------------------------------------------------------------------------------
! MODULE: formatCSRC
!
!> @author JL Tailhan
!> @author J  Foulliaron (version 1.0 - 2011)
!
!> @brief
!> Module permettant l'utilisation du format CSRC pour le stockage
!> des matrices creuses.
!----------------------------------------------------------------------------------------------------------
module formatCSRC

!------------------------------------------------------------------------------------------------------
!> @author JL Tailhan
!> @author J  Foulliaron (version 1.0 - 2011)
!
!> @brief
!> Declaration du type matrice creuse carree au format CSRC
!------------------------------------------------------------------------------------------------------
    type matrice_CSRC

        ! Tableaux des valeurs non nulles

        real*8,  dimension(:), allocatable :: Da ! valeurs diagonales
        real*8,  dimension(:), allocatable :: Ua ! valeurs dans le triangle superieur
        real*8,  dimension(:), allocatable :: La ! valeurs dans le triangle inferieur
        integer, dimension(:), allocatable :: ja, ia ! tableaux concernants les indices de reperage des valeurs

        integer :: nvalmax  ! Nombre de valeurs non nulles maximales dans un bloc triangulaire
        integer :: nl       ! nombre de lignes de la matrice

    end type matrice_CSRC

contains


!------------------------------------------------------------------------------------------------------
!> @author JL Tailhan
!> @author J  Foulliaron (version 1.0 - 2011)
!
!> @brief
!> Routine permettant d'initialiser juste en memoire une matrice creuse au format CSRC
!
!> @details
!> #### DESCRIPTION:
!>  A faire
!------------------------------------------------------------------------------------------------------
subroutine CSRC_init(mat,nlig,nzmax)

    !**************************************************************************************
    !* DEFINITIONS
    !**************************************************************************************

    implicit none

    !--------------------------------------------------------------------------------------
    !--- Variables in/out
    type (matrice_CSRC), intent(inout) :: mat ! matrice creuse à declarer en memoire
    integer, intent(in) :: nlig  ! nombre de lignes de la matrice à declarer
    integer, intent(in) :: nzmax ! nombre de valeurs non nulles dans un des triangles inf ou sup.

    !**************************************************************************************
    !* CORPS DE PROCEDURE
    !**************************************************************************************

    mat%nl = nlig
    mat%nvalmax = nzmax

    allocate(mat%Da(nlig))  ; allocate(mat%Ua(nzmax))  ;  allocate (mat%La(nzmax))
    allocate(mat%ja(nzmax)) ; allocate(mat%ia(nlig+1))

    mat%Da = 0.d0 ; mat%Ua = 0.d0 ; mat%La = 0.d0 ; mat%ia = 1 ; mat%ja = 0

end subroutine CSRC_init


!------------------------------------------------------------------------------------------------------
!> @author JL Tailhan
!> @author J  Foulliaron (version 1.0 - 2011)
!
!> @brief
!> Procedure qui libere la place memoire reservee à une matrice stockee au format CSRC
!
!> @details
!> #### DESCRIPTION:
!>  A faire
!------------------------------------------------------------------------------------------------------
subroutine CSRC_free(mat)

    !**************************************************************************************
    !* DEFINITIONS
    !**************************************************************************************

    implicit none

    !--------------------------------------------------------------------------------------
    !--- Variables in/out
    type(matrice_CSRC), intent(inout) :: mat ! Matrice creuse dont on veut liberer la memoire

    !**************************************************************************************
    !* CORPS DE PROCEDURE
    !**************************************************************************************

    deallocate(mat%Da)
    deallocate(mat%ja)
    deallocate(mat%Ua)
    deallocate(mat%La)
    deallocate(mat%ia)

end subroutine CSRC_free


!------------------------------------------------------------------------------------------------------
!> @author JL Tailhan
!> @author J  Foulliaron (version 1.0 - 2011)
!
!> @brief
!> Procedure reduisant au minimum, la taille memoire de la matrice creuse au format CSRC
!
!> @warning
!> En theorie inutile !
!
!> @details
!> #### DESCRIPTION:
!>  A faire
!------------------------------------------------------------------------------------------------------
subroutine CSRC_reduc_taille(mat)

    !**************************************************************************************
    !* DEFINITIONS
    !**************************************************************************************

    implicit none

    !--------------------------------------------------------------------------------------
    !--- Variables in/out
    type(matrice_CSRC), intent(inout) :: mat ! matrice dont on veut reduire la taille

    !--------------------------------------------------------------------------------------
    !--- Variables locales
    real*8, dimension(:), allocatable :: aux1 ! tab de passage pour la copie des valeurs de Ua et La
    integer, dimension(:), allocatable :: aux2! tab de passage pour la copie des valeurs de ja
    integer :: taille

    !**************************************************************************************
    !* CORPS DE PROCEDURE
    !**************************************************************************************

    taille = mat%ia(mat%nl+1) - 1 ! Donne le nombre de valeurs non nulles dans un des triangles

    !--------------------------------------------------------------------------------------
    !--- Reduction de la taille des tableaux de valeurs mat%La et mat%Ua
    allocate(aux1(taille));
    aux1 = mat%La(1:taille);
    deallocate(mat%La);
    allocate(mat%La(taille));
    mat%La = aux1
    aux1 = mat%Ua(1:taille);
    deallocate(mat%Ua);
    allocate(mat%Ua(taille));
    mat%Ua=aux1
    deallocate(aux1)

    !--------------------------------------------------------------------------------------
    !--- Reduction de la taille du tableau des indices des colonnes ja
    allocate(aux2(taille));
    aux2 = mat%ja(1:taille);
    deallocate(mat%ja);
    allocate(mat%ja(taille));
    mat%ja = aux2
    deallocate(aux2)

    mat%nvalmax = taille

end subroutine CSRC_reduc_taille


!------------------------------------------------------------------------------------------------------
!> @author JL Tailhan
!> @author J  Foulliaron (version 1.0 - 2011)
!
!> @brief
!> Routine effectuant le produit d'une matrice au format CSRC avec un vecteur
!
!> @details
!> #### DESCRIPTION:
!>  A faire
!------------------------------------------------------------------------------------------------------
subroutine ProdMatVect(A,u,p)

    !**************************************************************************************
    !* DEFINITIONS
    !**************************************************************************************

    implicit none

    !--------------------------------------------------------------------------------------
    !--- Variables in/out
    type(matrice_CSRC), intent(in) :: A   ! matrice creuse au format CSRC
    real*8, dimension(:), intent(in) :: u   ! vecteur d'entree
    real*8, dimension(:), intent(inout) :: p ! vecteur resultat

    !--------------------------------------------------------------------------------------
    !--- Variables in/out
    integer :: i, k
    real*8 :: t

    !**************************************************************************************
    !* CORPS DE PROCEDURE
    !**************************************************************************************

    p = 0.d0
    do i = 1, A%nl
        t = A%Da(i)*u(i)     ! calcul de la contribution diagonale

        do k = A%ia(i),A%ia(i+1) - 1
            t = t + A%La(k)*u(A%ja(k)) ! contribution de la portion de la i eme ligne du triangle inferieur de A sur P(i)
            p(A%ja(k)) = p(A%ja(k))+A%Ua(k)*u(i) !  contribution de la colonne associee sur les autres termes de p
        enddo
        p(i) = t
    enddo

end subroutine ProdMatVect


!------------------------------------------------------------------------------------------------------
!> @author JL Tailhan
!> @author J  Foulliaron (version 1.0 - 2011)
!
!> @brief
!> Routine effectuant le produit d'une matrice symetrique au format CSRC avec un vecteur
!
!> @details
!> #### DESCRIPTION:
!>  A faire
!------------------------------------------------------------------------------------------------------
subroutine ProdMatSymVect(A,u,p)

    !**************************************************************************************
    !* DEFINITIONS
    !**************************************************************************************

    implicit none

    !--------------------------------------------------------------------------------------
    !--- Variables in/out
    type(matrice_CSRC), intent(in) :: A      ! matrice creuse au format CSRC
    real*8, dimension(:), intent(in) :: u    ! vecteur d'entree
    real*8, dimension(:), intent(inout) :: p ! vecteur resultat

    !--------------------------------------------------------------------------------------
    !--- Variables locales
    integer :: i, k
    real*8 :: t

    !**************************************************************************************
    !* CORPS DE PROCEDURE
    !**************************************************************************************

    p = 0.d0

    do i = 1, A%nl
        t = A%Da(i)*u(i) ! calcul de la contribution diagonale

        do k = A%ia(i),A%ia(i+1) - 1
            t=t+A%La(k)*u(A%ja(k))
            p(A%ja(k))=p(A%ja(k))+A%La(k)*u(i)
        enddo
        p(i) = t
    enddo

end subroutine ProdMatSymVect


!------------------------------------------------------------------------------------------------------
!> @author JL Tailhan
!> @author J  Foulliaron (version 1.0 - 2011)
!
!> @brief
!> Routine permettant l'assemblage rapide d'une matrice elementaire dense
!> dans une matrice au format CSRC
!
!> @details
!> #### DESCRIPTION:
!>  A faire
!------------------------------------------------------------------------------------------------------
subroutine Assem_rapide22(matg,matelem,indice)

    !**************************************************************************************
    !* DEFINITIONS
    !**************************************************************************************

    implicit none

    !--------------------------------------------------------------------------------------
    !--- Variables in/out
    type(matrice_CSRC), intent(inout) :: matg  ! matrice creuse au format CSRC
    real*8, dimension(:,:), intent(in):: matelem ! matrice elementaire
    integer, dimension(:), intent(in) :: indice ! indices ou l'on doit inserer les valeurs de la matrice

    !--------------------------------------------------------------------------------------
    !--- Variables locales
    integer :: i, j, k, tailleind

    !**************************************************************************************
    !* CORPS DE PROCEDURE
    !**************************************************************************************

    !--------------------------------------------------------------------------------------
    !--- test, on verifie que les indices dinsertion sont corrects
    !if (any(indice > matg%nl)) STOP 'formatCSRC: mauvais indicage pour insertion dans la matrice globale'

    do i=1, size(indice)
        if(indice(i)>matg%nl) STOP 'formatCSRC: mauvais indicage pour insertion dans la matrice globale'
    enddo

    tailleind = size(indice)

    !--------------------------------------------------------------------------------------
    !--- On traite d'abord les valeurs diagonales de matelem
    do i = 1, tailleind
        matg%Da(indice(i)) = matg%Da(indice(i)) + matelem(i,i)
    enddo

    !--------------------------------------------------------------------------------------
    !--- On passe aux valeurs non diagonales de matelem
    do i = 2, tailleind ! On parcourt sur les valeurs du triangle inferieur de matelem
        do j = 1, i-1
            !------------------------------------------------------------------------------
            !--- Cas ou l'on insere la valeur dans le triangle inferieur de matg
            if (indice(j)<indice(i)) then
                do k = matg%ia(indice(i)),matg%ia(indice(i)+1)-1
                    if (matg%ja(k)==indice(j)) then
                        if (matelem(i,j)/=0.d0) matg%La(k)=matg%La(k)+matelem(i,j)
                        if (matelem(j,i)/=0.d0) matg%Ua(k)=matg%Ua(k)+matelem(j,i)!tant qu a faire on insere aussi la valeur symetrique de matelem
                        exit
                    endif
                enddo
            endif

            !------------------------------------------------------------------------------
            !--- Cas ou l'on insere la valeur dans le triangle superieur de matg
            if (indice(i)<indice(j)) then
                do k = matg%ia(indice(j)),matg%ia(indice(j)+1)-1
                    if( matg%ja(k)==indice(i)) then
                        if (matelem(i,j)/=0.d0) matg%Ua(k)=matg%Ua(k)+matelem(i,j)
                        if (matelem(j,i)/=0.d0) matg%La(k)=matg%La(k)+matelem(j,i) !tant qu a faire on insère aussi la valeur symetrique de matelem
                        exit
                    endif
                enddo
            endif
        enddo
    enddo

end subroutine Assem_rapide22


!------------------------------------------------------------------------------------------------------
!> @author JL Tailhan
!> @author J  Foulliaron (version 1.0 - 2011)
!
!> @brief
!> Teste si une matrice est bien symetrique
!
!> @details
!> #### DESCRIPTION:
!>  A faire
!------------------------------------------------------------------------------------------------------
subroutine testSYM_CSRC(A,sym)

    !**************************************************************************************
    !* DEFINITIONS
    !**************************************************************************************

    implicit none

    !--------------------------------------------------------------------------------------
    !--- Variables in/out
    logical,intent(inout):: SYM
    type (matrice_CSRC),intent(in)::A

    !--------------------------------------------------------------------------------------
    !--- Variables locales
    integer :: i

    !**************************************************************************************
    !* CORPS DE PROCEDURE
    !**************************************************************************************

    SYM = .true.

    do i = 1, A%ia(A%nl+1)-1
        if (A%La(i)/=A%Ua(i)) then
            SYM = .false.
            exit
        endif
    enddo

    if (SYM) then
        print*, 'la matrice est symetrique'
    else
        print*,'la matrice nest pas symetrique'
    endif

end subroutine testSYM_CSRC


!------------------------------------------------------------------------------------------------------
!> @author JL Tailhan
!> @author J  Foulliaron (version 1.0 - 2011)
!
!> @brief
!> Routine accedant à l'element i,j d'une matrice CSRC
!
!> @details
!> #### DESCRIPTION:
!>  A faire
!------------------------------------------------------------------------------------------------------
function CSRC_val(A,i,j) result(valeur)

    !**************************************************************************************
    !* DEFINITIONS
    !**************************************************************************************

    implicit none

    !--------------------------------------------------------------------------------------
    !--- Variables in/out
    type(matrice_CSRC), intent(in) :: A
    integer, intent(in) :: i, j

    !--------------------------------------------------------------------------------------
    !--- Variables locales
    real*8 :: valeur
    integer :: k

    !**************************************************************************************
    !* CORPS DE PROCEDURE
    !**************************************************************************************

    valeur = 0.d0

    if (i>A%nl.or.i<0.or.i>A%nl.or.i<0.or.i>A%nl) &
        & stop 'CSRC_val : erreur, les valeurs sont incorrectes, hors champs'

    if (i==j) then
        valeur = A%Da(i)
    elseif (i>j) then
        do k = A%ia(i),A%ia(i+1)-1
            if (j==A%ja(k)) then
                valeur = A%La(k)
                exit
            endif
        enddo
    else
        do k = A%ia(j), A%ia(j+1)-1
            if (i==A%ja(k)) then
                valeur = A%Ua(k)
                exit
            endif
        enddo
    endif

end function CSRC_val


!------------------------------------------------------------------------------------------------------
!> @author JL Tailhan
!> @author J  Foulliaron (version 1.0 - 2011)
!
!> @brief
!> Routine modifiant l'element i,j d'une matrice CSRC
!
!> @details
!> #### DESCRIPTION:
!>  (si la position existe deja sinon ne fais rien)
!------------------------------------------------------------------------------------------------------
subroutine CSRC_set(A,i,j,valeur)

    !**************************************************************************************
    !* DEFINITIONS
    !**************************************************************************************

    implicit none

    !--------------------------------------------------------------------------------------
    !--- Variables in/out
    type(matrice_CSRC), intent(inout) :: A
    integer, intent(in) :: i, j
    real*8, intent(in) :: valeur ! nouvelle valeur

    !--------------------------------------------------------------------------------------
    !--- Variables locales
    integer :: k

    !**************************************************************************************
    !* CORPS DE PROCEDURE
    !**************************************************************************************

    if (i>A%nl.or.i<0.or.i>A%nl.or.i<0.or.i>A%nl) &
        &  stop 'CSRC_set :erreur, les valeurs sont incorrectes, hors champs'

    if (i==j) then
        A%Da(i) = valeur
    elseif (i>j) then
        do k = A%ia(i),A%ia(i+1)-1
            if (j==A%ja(k)) then
                A%La(k) = valeur
                exit
            endif
        enddo
    else
        do k = A%ia(j), A%ia(j+1)-1
            if (i==A%ja(k)) then
                A%Ua(k) = valeur
                exit
            endif
        enddo
    endif

end subroutine CSRC_set


!------------------------------------------------------------------------------------------------------
!> @author JL Tailhan
!> @author J  Foulliaron (version 1.0 - 2011)
!
!> @brief
!> Routine mettant a 0. la ligne et la colonne "ic" de la matrice CSRC
!
!> @details
!> #### DESCRIPTION:
!>  A faire
!------------------------------------------------------------------------------------------------------
subroutine CSRC_lig_col_zero(A,ic)

    !**************************************************************************************
    !* DEFINITIONS
    !**************************************************************************************

    implicit none

    !--------------------------------------------------------------------------------------
    !--- Variables in/out
    type(matrice_CSRC), intent(inout) :: A
    integer, intent(in) :: ic

    !--------------------------------------------------------------------------------------
    !--- Variables locales
    integer :: k, m

    !**************************************************************************************
    !* CORPS DE PROCEDURE
    !**************************************************************************************

    A%Da(ic) = 0.d0

    !--------------------------------------------------------------------------------------
    !--- 1ere partie des valeurs mises à 0 (axes -x,+y)
    do k = A%ia(ic), A%ia(ic+1)-1
        A%La(k) = 0.d0
        A%Ua(k) = 0.d0
    enddo

    !--------------------------------------------------------------------------------------
    !--- 2eme partie des valeurs mises a 0 (axes -y,+x)
    do k = ic+1,A%nl
        do m=A%ia(k),A%ia(k+1)-1
            if(A%ja(m)==ic) then
                A%La(m)=0.d0
                A%Ua(m)=0.d0
                exit
            endif
        enddo
    enddo

end subroutine CSRC_lig_col_zero

!------------------------------------------------------------------------------------------------------

end module formatCSRC

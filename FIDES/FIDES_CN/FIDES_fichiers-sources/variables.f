module variables

!**************************************************************!
!       Declaration des variables globales et des stuctures    !
!**************************************************************!
    use formatCSRC

    implicit none
    save

!------------------------------------ Configuration -----------------------------!

    character(len = 14) :: dirdata='../FIDES_data/'
    character(len = 19) :: dirresu='../FIDES_resultats/'
    character(len = 14) :: dircsrc='../FIDES_CSRC/'

    integer, parameter :: intkind = 4    ! (4=32 bits, 8=64 bits) pour fixer le type d'entiers

!------------------------------------- Parametres -------------------------------!

    integer (kind=intkind) :: nnt    ! Nombre total de noeuds
    integer (kind=intkind) :: nelt   ! Nombre total d'elements
    integer (kind=intkind) :: ndlt   ! Nombre total de ddl (donne par knddl)

    integer :: ntypel ! Nombre de type d'element
    integer :: dime   ! Dimension du pb
    integer :: varint = 10 ! Nombre de variables internes maximales pour un modele
    integer, dimension(3) :: dirimpo = 0  ! Direction du chargement impose (en depl ou en force)

!------------------------------------- Constantes -------------------------------!

    real*8, parameter :: pi = 4.*atan(1.d0)

!--------------------------- Vecteurs et matrices de travail --------------------!
    real*8, dimension(:), allocatable   :: vrtrav
    real*8, dimension(:,:), allocatable :: mrtrav
    integer (kind=intkind), dimension(:), allocatable   :: vitrav
    integer (kind=intkind), dimension(:,:), allocatable :: mitrav
    integer (kind=intkind), dimension(100) :: kloce
    logical :: interne ! sert a choisir si les Cl sont internes ou externes (pas implentée pour l'instant)

!------------- Variables pour test avec nouveaux outils de FIDES-----------------!

    type (matrice_CSRC) :: supervkg  ! matrice principale au format CSRC
    integer :: CHOIX
    integer :: opt1 = 0
    integer (kind=intkind) :: ta_max
    integer :: ta_texte = 16 ! nombre de valeurs par ligne lors de la création des fichiers .csrc externes .

!----------------------------- Info relatives aux elements ----------------------!

! Definition d'un type representant le type des elements (ex MBT3, MBQ4,...)
    type forme

        integer :: nnel  ! Nombre de noeuds par element
        integer :: ndln  ! Nombre de ddl par noeud
        integer :: order ! Ordre d'integration element normal
        integer :: ncgr  ! Nombre de composantes du gradient

        real*8, dimension(:,:), allocatable :: Q     ! Point de Gauss
        real*8, dimension(:), allocatable :: W       ! Point de Gauss
        real*8, dimension(:,:), allocatable :: ksin  ! Coordonnees reduites des noeuds
        integer (kind=intkind), dimension(:,:), allocatable :: face ! Numerotation des points par face

    end type forme

    type(forme), dimension(:), allocatable :: infele


    type ddl_noeud
        integer (kind=intkind), dimension(:), allocatable :: dln    ! Numerotation des ddl par noeuds
    end type ddl_noeud

    type(ddl_noeud), dimension(:), allocatable :: infnoe

    type kelem
        integer (kind=intkind), dimension(:), allocatable :: el   ! Liste des elements voisins d'un element
        integer :: cote
    end type kelem

    type(kelem), dimension(:), allocatable :: kconelem

    type knoe
        integer(kind=intkind), dimension(:), allocatable :: no   ! Liste des elements voisins d'un noeud
    end type knoe

    type(knoe), dimension(:), allocatable :: kconoeud


!-------------------------------- Maillage et allocation -------------------------!

    real*8, dimension(:,:), allocatable :: vcor    ! vcor(dime, nnt) : coordonnees
    integer (kind=intkind), dimension(:), allocatable  :: ktypel   ! ktypel(nelt) : numero du type d'element
    integer, dimension(:), allocatable  :: knddl    ! knddl(nnt) : nombre de ddl par noeud, donne egalement ndlt.
    integer (kind=intkind), dimension(:,:), allocatable :: kconec ! kconec(nelt, nnel): connectivites (sparse)
    character(len=5), dimension(:), allocatable :: nomtype ! nomtype(ntypel)
    type(kelem), dimension(:), allocatable :: listelemautour ! tableau de la liste des elements autour de chaque noeuds
    integer (kind=intkind), dimension(:), allocatable  :: elemlibre ! tableau de la liste des elements libres

!-------------------------------------- Proprietes -------------------------------!

    integer :: idmax, nbrgrp
    real*8, dimension(:,:), allocatable :: vprelg  ! vprel(nprop, idmax) : Proprietes elementaires
    real*8, dimension(:,:), allocatable :: vprelg0 ! vprel(nprop, idmax) : Proprietes elementaires pour algo implicite
    real*8, dimension(:,:), allocatable :: vprelg1 ! vprel(nprop, idmax) : Proprietes elementaires pour algo implicite
    integer (kind=intkind), dimension(:), allocatable :: kprop, kprop0   ! kprop(nelt) : Numero du type de prop elementaire

!------------------------------------- Chargement et CL --------------------------!

    logical, dimension(:), allocatable :: kcond ! kcond(ndlt) : indice condition limite imposee
    real*8, dimension(:), allocatable :: vcond  ! vcond(ndlt) : valeur codition limite imposee (sparse)
    real*8, dimension(:), allocatable :: vfcg   ! vfcg(ndlt) : sollicitation concentree lecture (sparse)
    real*8, dimension(:), allocatable :: aval   ! aval(npas)
    integer :: npilot

!----------------------------------------- Solution ------------------------------!

    real*8, dimension(:), allocatable :: vsol   ! vsol(ndlt) : solution globale
    real*8, dimension(:), allocatable :: vdu    ! vdu(ndlt)  : increment de solution globale dans processus iteratif
    real*8, dimension(:), allocatable :: vcont  ! contraintes  : sxx, syy, sxy (pt de Gauss)
    real*8, dimension(:), allocatable :: vnoli  ! deformation plastique : epxx, epyy, epxy et energ. post-pic : Wp (pt de Gauss)
    real*8, dimension(:), allocatable :: vrint  ! autres variables internes des modeles (pt de Gauss)
    integer (kind=intkind), dimension(:), allocatable :: kcont ! Pour situer les elems dans vcont
    integer (kind=intkind), dimension(:), allocatable :: knoli ! Pour situer les elems dans vnoli
    integer (kind=intkind), dimension(:), allocatable :: krint ! Pour situer les elems dans vrint

    type resultat
       real*8, dimension(:), allocatable :: vsolu
       real*8, dimension(:), allocatable :: vcontr
       real*8, dimension(:), allocatable :: vnolin
       integer, dimension(:,:), allocatable :: etatpg
       integer, dimension(:), allocatable :: etat
       real*8, dimension(:,:), allocatable :: inorme
       real*8, dimension(:), allocatable :: dpglobal, foglobal
    end type resultat

    type(resultat), dimension(:), allocatable :: resu

!----------------------------------- Structure de resulats -----------------------!

    real*8, dimension(3) :: dpglobal,foglobal
    logical :: reacn
    integer :: diri
    integer (kind=intkind) :: noe1, noe2

!------------------------------------- Algorithme --------------------------------!

    real*8, dimension(:), allocatable :: convdu,convdf   ! vecteurs de convergence sur u et f
    integer :: niter, npas, ipas, iter
    integer :: imet
    real*8 :: ndmax
    logical :: explicite

!----------------------------- Algorithme (pb evolutif) --------------------------!

    real*8, dimension(:), allocatable :: dtps,tps   ! (d)tps(npas) : listes de pas de temps physiques pour pb evolutif (implicite)
    real*8 :: Dtm1, Dtm2, Dtrest, alp_Dt
    logical :: ievtp

!--------------------------------------- Reprise ---------------------------------!

    logical :: nrep
    integer :: optrep
    character(len=50) :: nomfichrep

!----------------------------------- Pilotage indirect ---------------------------!

    logical :: indirect

!--------------------------------- Contraintes initiales -------------------------!

    logical :: contini

!--------------------- Definition d'un effort de precontrainte -------------------!

    real*8, dimension(:), allocatable :: vfpg   ! vfpg(ndlt) : sollicitation concentree
    logical :: precont=.false.

!------------------------------ Lieu de calcul des contraintes -------------------!

    character(len=5) :: calco

!---------------------------------------- Impression -----------------------------!

    integer :: iprint

!----------- Sauvegarde des resultats (a chaque calcul : isauvcal=1)--------------!

    integer :: isauvcalc=0
    real*8  :: dlamsauv

!------------------------------------Post-traitement -----------------------------!

    integer :: ipost

!--------------- Informations sur l'elements le plus dangereux -------------------!

    integer (kind=intkind) :: iedng
    
!------------------ Informations sur les elements critiques ----------------------!

    integer (kind=intkind) :: iecritic

!----------- Definition d'un facies de fissuration comme graphe ------------------!

    type facies2D
       integer, dimension(:), allocatable   :: vertex   ! Extremite d'une branche de fissure
       real*8, dimension(:,:), allocatable  :: vxyz     ! Coordonnées des vertex
       integer, dimension(:), allocatable   :: vdeg     ! Degre du vertex (nbre d'elements connectes à l'extremite)
       integer, dimension(:,:), allocatable :: edge     ! Branche de fissure
    end type facies2D

    type(facies2D), dimension(:), allocatable :: fissures2D

!-------------------- Informations sur les elements d'interface ------------------!

    integer, dimension(:,:), allocatable  :: ietatpg, histetatpg1, histetatpg2
    real*8, dimension(:,:), allocatable  :: endo
    integer (kind=intkind) :: iedngi
    integer :: interf

!-------------------- Informations pour le modele beton fissurant ----------------!

    real*8, dimension(:,:), allocatable :: inorm, irloc
    real*8, dimension(:), allocatable :: iouver1, iouver2
    integer, dimension(:), allocatable  :: ietat, histetat1, histetat2
    integer (kind=intkind), dimension(3) :: listcompfiss=(/5, 12, 15/)
    integer (kind=intkind) :: iedngm
    integer :: fiss, smear, depen !, icompfiss

    real*8 :: ouvermax = 20.      !  verification du deplacement des noeuds

!------------------------- Informations pour l'acier ------------------------------!

    real*8, dimension(:), allocatable :: limels, limrup, ecroui
    logical :: acierl = .false.

!------------------ Informations pour le modele beton arme fissurant --------------!
    integer :: fissBA !, smear, depen, icompfiss
    integer (kind=intkind), dimension(3) :: listcompfissBA=(/25, 26, 27/)
    integer (kind=intkind) :: iedngba
    integer, dimension(1500) :: check = 0
    logical :: dist_triline

!-------------------- Informations sur les groupes probabilistes ------------------!

    type probab
        integer :: num       ! numero du groupe
        integer :: loi       ! numero de loi proba
        integer :: ipa       ! nbre de parametre de la loi proba
        real*8, dimension(6) :: param    ! liste des parametres
    end type probab

    type(probab), dimension(:), allocatable :: young,resist,enerj,enerjpf,nrjbeto,nrjfibr
    type(probab), dimension(:), allocatable :: enerjba, enerjba2
    logical :: alea

!--------------------------- Auto contraintes initiales --------------------------!

    integer, dimension(:), allocatable :: gpinicont
    logical :: inict

!-------------------- Informations sur la probabilisation des comp. béton  -------!

    real*8  :: Dg,fc
    logical :: redu
    real*8  :: rap

!------- Pour detecter les oscillations dans le cas des elements d interface

    integer, dimension(:), allocatable :: detoscill, detoscillp, deboscill
    integer :: comptoscill = 0

    type(kelem),dimension(:),allocatable::osc_interf

!------ Vecteur des numeros de noeuds pour deplacement (monitoring points) -------!
    integer, dimension(:), allocatable :: vnonoedp
    integer :: idirdp


    logical :: affiche = .false.


end module variables

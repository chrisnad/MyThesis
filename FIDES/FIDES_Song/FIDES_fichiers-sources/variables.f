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
   character(len =14 ) :: dircsrc='../FIDES_CSRC/' 
   
!------------------------------------- Parametres -------------------------------!

   integer :: nnt    ! Nombre de nds total
   integer :: nelt   ! Nombre d'elemnts total
   integer :: ndlt   ! Nombre de ddl total (donne par knddl)
   integer :: ntypel ! Nombre de type d'element
   integer :: dime   ! Dimension du pb
   integer, dimension(3) :: dirimpo = 0  ! Direction du chargement impose (en depl ou en force)

!------------------------------------- Constantes -------------------------------!
   real(kind=8), parameter :: pi = 4.*atan(1.d0)

!--------------------------- Vecteurs et matrices de travail --------------------!
   real*8, dimension(:), allocatable   :: vrtrav
   real*8, dimension(:,:), allocatable :: mrtrav
   integer, dimension(:), allocatable   :: vitrav
   integer, dimension(:,:), allocatable :: mitrav
   integer, dimension(100) :: kloce
   logical :: interne ! sert a choisir si les Cl sont internes ou externes (pas implentée pour l'instant)

!------------- Variables pour test avec nouveaux outils de FIDES-----------------!

   type (matrice_CSRC) :: supervkg  ! matrice principale au format CSRC
   integer :: CHOIX 
   integer :: opt1 = 0 
   integer(kind = 8) :: ta_max  
   integer :: ta_texte = 16 ! nombre de valeurs par ligne lors de la création des fichiers .csrc externes .
!----------------------------- Info relatives aux elements ----------------------!

! Definition d'un type representant le type des elements (ex MBT3, MBQ4,...)
   type forme 

       integer :: nnel  ! Nombre de noeuds 
       integer :: ndln  ! Nombre de ddl par noeud
       integer :: order ! Ordre d'integration element normal
       integer :: split ! Ordre d'integration element split
       integer :: ncgr  ! Nombre de composantes du gradient

       real*8, dimension(:,:), allocatable :: Q     ! Point de Gauss 
       real*8, dimension(:), allocatable :: W       ! Point de Gauss 
       real*8, dimension(:,:), allocatable :: ksin  ! Coordonnees reduites des noueds
       integer, dimension(:,:), allocatable :: face    ! Numerotation des points par face 

   end type forme

   type(forme), dimension(:), allocatable :: infele


   type ddl_noeud 
       integer, dimension(:), allocatable :: dln    ! Numerotation des ddl par noeuds
   end type ddl_noeud 

   type(ddl_noeud), dimension(:), allocatable :: infnoe

   type kelem
       integer, dimension(:), allocatable :: el   ! List des elements voisins d'un element
       integer :: cote   
   end type kelem

   type(kelem), dimension(:), allocatable :: kconelem

   type knoe
       integer, dimension(:), allocatable :: no   ! List des elements voisins d'un noeud
   end type knoe

   type(knoe), dimension(:), allocatable :: kconoeud


!-------------------------------- Maillage et allocation -------------------------!

   real*8, dimension(:,:), allocatable :: vcor    ! vcor(dime, nnt) : coordonnees
   integer, dimension(:), allocatable :: ktypel   ! ktypel(nelt) : numero du type d'element
   integer, dimension(:), allocatable :: knddl    ! knddl(nnt) : nombre de ddl par noeud, donne egalement ndlt.
   integer, dimension(:,:), allocatable :: kconec ! kconec(nelt, nnel): connectivites (sparse)
   character(len=5), dimension(:), allocatable :: nomtype ! nomtype(ntypel)
   type(kelem), dimension(:), allocatable :: listelemautour ! tableau de la liste des elements autour de chaque noeuds
!-------------------------------------- Proprietes -------------------------------!
   
   integer :: idmax 
   real*8, dimension(:,:), allocatable :: vprelg  ! vprel(nprop, idmax) : Proprietes elementaires
   real*8, dimension(:,:), allocatable :: vprelg0 ! vprel(nprop, idmax) : Proprietes elementaires pour algo implicite
   real*8, dimension(:,:), allocatable :: vprelg1 ! vprel(nprop, idmax) : Proprietes elementaires pour algo implicite
   integer, dimension(:), allocatable :: kprop, kprop0   ! kprop(nelt) : Numero du type de prop elementaire
   
!------------------------------------- Chargement et CL --------------------------!

   logical, dimension(:), allocatable :: kcond ! kcond(ndlt) : indice condition limite imposee
   real*8, dimension(:), allocatable :: vcond  ! vcond(ndlt) : valeur codition limite imposee (sparse)
   real*8, dimension(:), allocatable :: vfcg   ! vfcg(ndlt) : sollicitation concentree lecture (sparse)   
   real*8, dimension(:), allocatable :: aval   ! aval(npas)
   integer :: npilot

!----------------------------------------- Solution ------------------------------!
   real*8, dimension(:), allocatable :: vfgl   ! pour calcul anticipe dans evolglob
   real*8, dimension(:), allocatable :: vsol   ! vsol(ndlt) : solution globale 
   real*8, dimension(:), allocatable :: vcont  ! contraintes  : sxx, syy, sxy (pt de Gauss)
   real*8, dimension(:), allocatable :: vnoli  ! deformation plastique : epxx, epyy, epxy et energ. post-pic : Wp (pt de Gauss)
   real*8, dimension(:), allocatable :: vdu0  
   real*8, dimension(:), allocatable :: vdle0
   integer, dimension(:), allocatable :: kcont ! Pour situer les elems dans vcont
   integer, dimension(:), allocatable :: knoli ! Pour situer les elems dans vnoli
   
   type resultat
       real*8, dimension(:), allocatable :: vsolu
       real*8, dimension(:), allocatable :: vcontr
       real*8, dimension(:), allocatable :: vnolin
       integer, dimension(:,:), allocatable :: etatpg
       integer, dimension(:,:), allocatable :: etatma
       integer, dimension(:), allocatable :: etat
       real*8, dimension(:,:), allocatable :: inorme
       real*8, dimension(:), allocatable :: dpglobal, foglobal
   end type resultat
      
   type(resultat), dimension(:), allocatable :: resu
   
!----------------------------------- Structure de resulats -----------------------!
   real*8, dimension(3) :: dpglobal,foglobal
   logical :: reacn
   integer :: diri, noe1, noe2   
   
!------------------------------------- Algorithme --------------------------------!

   integer :: niter, npas, ipas, iter
   integer :: imet
   real*8 :: ndmax
   logical :: explicite
     
!----------------------------- Algorithme (pb evolutif) --------------------------!

   real*8, dimension(:), allocatable :: dtps,tps   ! (d)tps(npas) : listes de pas de temps physiques pour pb evolutif (implicite)
   logical :: ievtp
   
!--------------------------------------- Reprise ---------------------------------!

   logical :: nrep
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

!------------------------------------Post-traitement -----------------------------!

   integer :: ipost
    
!------------- Informations sur les elements le plus dangereux ------------------!
   integer :: iedng  

!------------- Informations sur les elements macro-interface --------------------!
   integer :: interf_macro, iedngma
   real*8 :: Emax = 150000.
   character(len=3), dimension(:), allocatable :: mode
   integer, dimension(:,:), allocatable  :: ietatma, histetatma1, histetatma2
   logical :: imacro = .false.

!-------------------- Informations sur les elements d'interface ------------------!
   integer, dimension(:,:), allocatable  :: ietatpg, histetatpg1, histetatpg2
   real*8, dimension(:,:), allocatable  :: endo, varg0, vargs
   integer :: iedngi, ieendo, interf
	
!-------------------- Informations pour le modele beton fissurant ----------------!
   real*8, dimension(:,:), allocatable :: inorm, irloc
   real*8, dimension(:), allocatable :: iouver1, iouver2
   integer, dimension(:), allocatable  :: ietat, histetat1, histetat2
   integer :: fiss, smear, depen, icompfiss, iedngm
   
   integer, dimension(:), allocatable :: detoscill
   integer :: comptoscill = 0
   real*8 :: ouvermax = 0.02      !  verification du deplacement des noeuds

!-------------------- Informations pour l'acier ----------------!
   real*8, dimension(:), allocatable :: limels, limrup, ecroui
   
!-------------------- Informations sur les groupes probabilistes ------------------!

   type probab
       integer :: num       ! numero du groupe
       integer :: loi       ! numero de loi proba
       integer :: ipa       ! nbre de parametre de la loi proba
       real*8, dimension(5) :: param    ! liste des parametres
   end type probab

   type(probab), dimension(:), allocatable :: young, resist, interfa, enerj
   logical :: alea

!--------------------------- Auto contraintes initiales --------------------------!

   integer, dimension(:), allocatable :: gpinicont
   logical :: inict

!-------------------- Informations sur la probabilisation des comp. béton  -------!
   real*8 :: Dg,fc
   logical :: redu
   real*8  :: rap

end module variables

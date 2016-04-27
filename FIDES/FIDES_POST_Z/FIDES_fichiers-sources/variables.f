!----------------------------------------------------------------------------------------------------------
! MODULE: variables
!
!> @author JL Tailhan
!> @author J. Goncalvez (version 1.0 - 2010)
!
!> @brief
!> Declaration des variables globales et des stuctures
!----------------------------------------------------------------------------------------------------------
module variables


    use formatCSRC

    implicit none
    save

    !--------------------------------------------------------------------------------------
    !--- Configuration
    character(len = 14) :: dirdata='../FIDES_data/'                 !< Repertoire du fichier de donnees *.data
    character(len = 19) :: dirresu='../FIDES_resultats/'            !< Repertoire des fichiers de resultats *.resu, *.list
    character(len = 14) :: dircsrc='../FIDES_CSRC/'                 !< Repertoire des fichiers des matrices CSR (obsolete)

    integer, parameter :: intkind = 4                               !< Type d'entiers (4=32 bits, 8=64 bits)

    !--------------------------------------------------------------------------------------
    !--- Parametres
    integer (kind=intkind) :: nnt                                   !< Nombre total de noeuds
    integer (kind=intkind) :: nelt                                  !< Nombre total d'elements
    integer (kind=intkind) :: ndlt                                  !< Nombre total de ddl (donne par knddl)

    integer :: ntypel                                               !< Nombre de type d'element
    integer :: dime                                                 !< Dimension du pb
    integer :: varint = 10                                          !< Nombre de variables internes maximales pour un modele
    integer, dimension(3) :: dirimpo = 0                            !< Direction du chargement impose (en depl ou en force)

    !--------------------------------------------------------------------------------------
    !--- Constantes
    real*8, parameter :: pi = 4.*atan(1.d0)                         !< Constante PI = 3.141926...

    !--------------------------------------------------------------------------------------
    !--- Vecteurs et matrices de travail
    real*8, dimension(:), allocatable   :: vrtrav                   !< Vecteurs de travail (reels)
    real*8, dimension(:,:), allocatable :: mrtrav                   !< Matrice de travail (reels)
    integer (kind=intkind), dimension(:), allocatable   :: vitrav   !< Vecteur de travail (entiers)
    integer (kind=intkind), dimension(:,:), allocatable :: mitrav   !< Matrice de travail (entiers)
    integer (kind=intkind), dimension(100) :: kloce                 !< Vecteur de localisation des ddls elementaires
    logical :: interne                      !< Booleen de choix pour Cl sont internes ou externes (pas implentée pour l'instant)

    !--------------------------------------------------------------------------------------
    !--- Variables pour test avec nouveaux outils de FIDES
    type (matrice_CSRC) :: supervkg                                 !< Matrice principale au format CSRC
    integer :: CHOIX                                                !< Entier pour choix 
    integer :: opt1 = 0                                             !< Indice d'option
    integer (kind=intkind) :: ta_max                                !< Taille maximale de matrice
    integer :: ta_texte = 16                !< Nombre de valeurs par ligne lors de la création des fichiers .csrc externes

    !--------------------------------------------------------------------------------------
    !--- Info relatives aux elements
    !    Definition d'un type representant le type des elements (ex MBT3, MBQ4,...)
    type forme                                                      !< Type forme representant le type d'element
        integer :: nnel                                             !< Nombre de noeuds par element
        integer :: ndln                                             !< Nombre de ddl par noeud
        integer :: order                                            !< Ordre d'integration element normal
        integer :: ncgr                                             !< Nombre de composantes du gradient

        real*8, dimension(:,:), allocatable :: Q                    !< Point de Gauss elementaires
        real*8, dimension(:), allocatable :: W                      !< Poids de Gauss elementaire
        real*8, dimension(:,:), allocatable :: ksin                 !< Coordonnees reduites des noeuds de l'element
        integer (kind=intkind), dimension(:,:), allocatable :: face !< Numerotation des points par face de l'element
    end type forme

    type(forme), dimension(:), allocatable :: infele                !< Informations sur les elements


    type ddl_noeud                                                  !< Structure de numerotation des ddls par noeud
        integer (kind=intkind), dimension(:), allocatable :: dln    !< Numeros ddl par noeud
    end type ddl_noeud

    type(ddl_noeud), dimension(:), allocatable :: infnoe            !< Informations sur les elements

    type kelem                                                      !< Structure de definition des voisins d'un element
        integer (kind=intkind), dimension(:), allocatable :: el     !< Liste des elements voisins d'un element
        integer :: cote                                             !< Liste des cotes de l'element
    end type kelem

    type(kelem), dimension(:), allocatable :: kconelem

    type knoe                                                       !< Structure de definition des elements voisins
        integer(kind=intkind), dimension(:), allocatable :: no      !< Liste des elements voisins d'un noeud
    end type knoe

    type(knoe), dimension(:), allocatable :: kconoeud


    !--------------------------------------------------------------------------------------
    !--- Maillage et allocation
    real*8, dimension(:,:), allocatable :: vcor                     !< vcor(dime, nnt) : vecteur des coordonnees nodales
    integer (kind=intkind), dimension(:), allocatable  :: ktypel    !< ktypel(nelt) : numero du type d'element
    integer, dimension(:), allocatable  :: knddl                    !< knddl(nnt) : nombre de ddl par noeud, donne egalement ndlt.
    integer (kind=intkind), dimension(:,:), allocatable :: kconec   !< kconec(nelt, nnel): connectivites (sparse)
    character(len=5), dimension(:), allocatable :: nomtype          !< nomtype(ntypel)
    type(kelem), dimension(:), allocatable :: listelemautour        !< tableau de la liste des elements autour de chaque noeuds
    integer (kind=intkind), dimension(:), allocatable  :: elemlibre !< tableau de la liste des elements libres

    !--------------------------------------------------------------------------------------
    !--- Proprietes
    integer :: idmax, nbrgrp
    real*8, dimension(:,:), allocatable :: vprelg   !< vprel(nprop, idmax) : Proprietes elementaires
    real*8, dimension(:,:), allocatable :: vprelg0  !< vprel(nprop, idmax) : Proprietes elementaires pour algo implicite
    real*8, dimension(:,:), allocatable :: vprelg1  !< vprel(nprop, idmax) : Proprietes elementaires pour algo implicite
    integer (kind=intkind), dimension(:), allocatable :: kprop, kprop0   !< kprop(nelt) : Numero du type de prop elementaire

    !--------------------------------------------------------------------------------------
    !--- Chargement et CL
    logical, dimension(:), allocatable :: kcond     !< kcond(ndlt) : indice condition limite imposee
    real*8, dimension(:), allocatable :: vcond      !< vcond(ndlt) : valeur codition limite imposee (sparse)
    real*8, dimension(:), allocatable :: vfcg       !< vfcg(ndlt) : sollicitation concentree lecture (sparse)
    real*8, dimension(:), allocatable :: aval       !< aval(npas)
    integer :: npilot

    !--------------------------------------------------------------------------------------
    !--- Solution
    real*8, dimension(:), allocatable :: vsol       !< vsol(ndlt) : solution globale
    real*8, dimension(:), allocatable :: vdu        !< vdu(ndlt)  : increment de solution globale dans processus iteratif
    real*8, dimension(:), allocatable :: vcont      !< contraintes  : sxx, syy, sxy (pt de Gauss)
    real*8, dimension(:), allocatable :: vnoli  !< deformation plastique : epxx, epyy, epxy et energ. post-pic : Wp (pt de Gauss)
    real*8, dimension(:), allocatable :: vrint      !< autres variables internes des modeles (pt de Gauss)
    integer (kind=intkind), dimension(:), allocatable :: kcont !< Pour situer les elems dans vcont
    integer (kind=intkind), dimension(:), allocatable :: knoli !< Pour situer les elems dans vnoli
    integer (kind=intkind), dimension(:), allocatable :: krint !< Pour situer les elems dans vrint

    type resultat                                       !< Structure de resultats
        real*8, dimension(:), allocatable :: vsolu      !< Vecteur solution
        real*8, dimension(:), allocatable :: vcontr     !< Vecteur des contraintes
        real*8, dimension(:), allocatable :: vnolin     !< Vecteur des quantites non-lineaire
        integer, dimension(:,:), allocatable :: etatpg  !< Vecteur des etats de fissuration des elements de contact
        integer, dimension(:), allocatable :: etat      !< Vecteur des etats de fissuration des elements massifs
        real*8, dimension(:,:), allocatable :: inorme   !< Vecteur d'orientation de la normale a la fissure
        real*8, dimension(:), allocatable :: dpglobal   !< Deplacements globaux
        real*8, dimension(:), allocatable :: foglobal   !< Efforts globaux
    end type resultat

    type(resultat), dimension(:), allocatable :: resu

    !--------------------------------------------------------------------------------------
    !--- Structure de resulats
    real*8, dimension(3) :: dpglobal
    real*8, dimension(3) :: foglobal
    logical :: reacn
    integer :: diri
    integer (kind=intkind) :: noe1, noe2

    !--------------------------------------------------------------------------------------
    !--- Algorithme
    real*8, dimension(:), allocatable :: convdu         !< vecteurs de convergence sur u et f
    real*8, dimension(:), allocatable :: convdf         !< vecteurs de convergence sur u et f
    integer :: niter                                    !< Nombre d'iterations
    integer :: npas                                     !< Nombre de pas de calcul
    integer :: ipas                                     !< Numero du pas de calcul courant
    integer :: iter                                     !< Numero de l'iteration courante
    integer :: imet                                     !< Numero de la methode de resolution (tangente, initiale, etc. ...)
    real*8  :: ndmax                                    !< Valeur du critere de convergence (non utilise)
    logical :: explicite                                !< Type de resolution (explicite, implicite)

    !--------------------------------------------------------------------------------------
    !--- Algorithme (pb evolutif)
    real*8, dimension(:), allocatable :: dtps,tps   !< (d)tps(npas) : listes de pas de temps physiques pour pb evolutif (implicite)
    real*8 :: Dtm1
    real*8 :: Dtm2
    real*8 :: Dtrest
    real*8 :: alp_Dt
    logical :: ievtp

    !--------------------------------------------------------------------------------------
    !--- Reprise
    logical :: nrep                                     !< Indice de reprise d'un calcul
    integer :: optrep                                   !< Option de reprise
    character(len=50) :: nomfichrep                     !< Nom du fichier de reprise

    !--------------------------------------------------------------------------------------
    !--- Pilotage indirect
    logical :: indirect                                 !< Booleen pour pilotage indirect

    !--------------------------------------------------------------------------------------
    !--- Contraintes initiales
    logical :: contini                                  !< Boolee pour calcul de contraintes initiales

    !--------------------------------------------------------------------------------------
    !--- Definition d'un effort de precontrainte
    real*8, dimension(:), allocatable :: vfpg           !< vfpg(ndlt) : sollicitation concentree
    logical :: precont=.false.                          !< Booleen pour prise en compte d'un effort de precontrainte

    !--------------------------------------------------------------------------------------
    !--- Lieu de calcul des contraintes
    character(len=5) :: calco                           !< Indice de lieu de calcul des contraintes (pts d'integr. ou noeuds)

    !--------------------------------------------------------------------------------------
    !--- Impression
    integer :: iprint                                   !< Indice d'impression

    !--------------------------------------------------------------------------------------
    !--- Sauvegarde des resultats (a chaque calcul : isauvcal=1)
    integer :: isauvcalc=0                              !< Indice de sauvegarde des resultats (a chaque calcul : isauvcal=1)
    real*8  :: dlamsauv                                 !< valeur de dlam a la sauvegarde

    !--------------------------------------------------------------------------------------
    !--- Post-traitement
    integer :: ipost                                    !< Indice de post-traitement

    !--------------------------------------------------------------------------------------
    !--- Informations sur l'elements le plus dangereux
    integer (kind=intkind) :: iedng                     !< Numero de l'element le plus dangereux
    
    !--------------------------------------------------------------------------------------
    !--- Informations sur les elements critiques
    integer (kind=intkind) :: nelcritic                 !< Nombre d'elements critiques (en cours de fissuration)

    !--------------------------------------------------------------------------------------
    !--- Definition d'un facies de fissuration comme graphe
    type facies2D                                       !< Structure de definition d'un facies de fissuration
       integer, dimension(:), allocatable   :: vertex   !< Extremite d'une branche de fissure
       real*8, dimension(:,:), allocatable  :: vxyz     !< Coordonnées des vertex
       integer, dimension(:), allocatable   :: vdeg     !< Degre du vertex (nbre d'elements connectes à l'extremite)
       integer, dimension(:,:), allocatable :: edge     !< Branche de fissure
    end type facies2D

    type(facies2D), dimension(:), allocatable :: fissures2D

    !--------------------------------------------------------------------------------------
    !--- Informations sur les elements d'interface
    integer, dimension(:,:), allocatable  :: ietatpg        !< Vecteur des etats de fissuration des elements
    integer, dimension(:,:), allocatable  :: histetatpg1    !< Vecteur des etats de fissuration des elements (pas n-1)
    integer, dimension(:,:), allocatable  :: histetatpg2    !< Vecteur des etats de fissuration des elements (pas n-2)
    real*8, dimension(:,:), allocatable   :: endo           !< Vecteur de stockage de l'endommagement (ele interface)
    integer (kind=intkind) :: iedngi                        !< Numero de l'element d'interface le plus dangereux
    integer :: interf                                       !< Switch d'utilisation des elements d'interface

    !--------------------------------------------------------------------------------------
    !--- Informations pour le modele beton fissurant
    real*8, dimension(:,:), allocatable :: inorm            !< Vecteur normal a la fissure
    real*8, dimension(:,:), allocatable :: irloc            !< Repere normal a la fissure
    real*8, dimension(:), allocatable   :: iouver1          !< ???
    real*8, dimension(:), allocatable   :: iouver2          !< ???
    integer, dimension(:), allocatable  :: ietat            !< Vecteur des etats de fissuration des elements
    integer, dimension(:), allocatable  :: histetat1        !< Vecteur des etats de fissuration des elements (pas n-1)
    integer, dimension(:), allocatable  :: histetat2        !< Vecteur des etats de fissuration des elements (pas n-2)
    integer (kind=intkind), dimension(2):: listcompfiss=(/12, 15/)  !< Liste des lois beton fissurant
    integer (kind=intkind) :: iedngm                        !< Numero de l'element massif le plus dangereux
    integer :: fiss                                         !< Switch pour utilisation du modele de beton fissurant
    integer :: smear                                        !< Switch pour approche de type smeared (obsolete)
    integer :: depen                                        !< Switch pour approche par depenalisation (obsolete)
    !integer :: icompfiss

    real*8 :: ouvermax = 20.      !  verification du deplacement des noeuds

    !--------------------------------------------------------------------------------------
    !--- Informations pour l'acier
    real*8, dimension(:), allocatable :: limels             !< ???
    real*8, dimension(:), allocatable :: limrup             !< ???
    real*8, dimension(:), allocatable :: ecroui             !< ???
    logical :: acierl = .false.

    !--------------------------------------------------------------------------------------
    !--- Informations pour le modele beton arme fissurant
    integer :: fissBA                                       !< Switch pour utilisation du modele de BA fissurant
    integer (kind=intkind), dimension(3) :: listcompfissBA=(/25, 26, 27/)   !< Liste des lois BA fissurant
    integer (kind=intkind) :: iedngba                       !< Numero de l'element BA le plus dangereux

    !--------------------------------------------------------------------------------------
    !--- Informations sur les groupes probabilistes
    type probab                                             !< Structure d'information sur les grpes probabilistes
        integer :: num                                      !< numero du groupe
        integer :: loi                                      !< numero de loi proba
        integer :: ipa                                      !< nbre de parametre de la loi proba
        real*8, dimension(5) :: param                       !< liste des parametres
    end type probab

    type(probab), dimension(:), allocatable :: young,resist,enerj,enerjpf,nrjbeto,nrjfibr
    type(probab), dimension(:), allocatable :: enerjba, enerjba2
    logical :: alea

    !--------------------------------------------------------------------------------------
    !--- Auto contraintes initiales
    integer, dimension(:), allocatable :: gpinicont         !< ??? (groupe a autocontraintes initiales ?)
    logical :: inict                                        !< Switch pour activation autocontraintes initiales

    !--------------------------------------------------------------------------------------
    !--- Informations sur la probabilisation des comp. béton
    real*8  :: Dg                                           !< Diametre du plus gros granulat
    real*8  :: fc                                           !< Resistance a la compression
    logical :: redu                                         !< ???
    real*8  :: rap                                          !< ??? (rapport ve/vg ?)

    !--------------------------------------------------------------------------------------
    !--- Pour detecter les oscillations dans le cas des elements d interface
    integer, dimension(:), allocatable :: detoscill         !< Tableau de comptage des oscillations
    integer :: comptoscill = 0                              !< Compteur des oscillations
    type(kelem),dimension(:),allocatable::osc_interf        !< Elements d'interface en oscillation

    !--------------------------------------------------------------------------------------
    !--- Vecteur des numeros de noeuds pour deplacement (monitoring points)
    integer, dimension(:), allocatable :: vnonoedp          !< ???
    integer :: idirdp                                       !< Indice de direction de deplacement (???)
    
!------ Numero du fichier a post-traité
    integer :: nfich

end module variables

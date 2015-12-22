module interface_macro

contains

!********************************************************!

subroutine interf_macro_endo(vsig,vnle,vh,vb,vnl0,vprel,vdle,iloi,ie,ipg)

!********************************************************!
!     Calcul des contraintes pour le calcul fissure      !
!********************************************************!

        use variables, only : dime, inict, irloc, iedng, ietatma, inorm, ipas, endo, mode, varg0, vargs, ktypel, nomtype
        use initialisation, only : init_vec, init_mat
        use utilitaire, only : princ
        use fissuration, only : fiss_changement_repere
        implicit none

        !--- Variables IN
        real*8, dimension(:,:), intent(in) :: vb
        real*8, dimension(:,:), intent(inout) :: vh
        real*8, dimension(:), intent(in) :: vnl0, vdle, vprel
        integer, intent(in) :: iloi, ie, ipg
        
        !--- Variables OUT
        real*8, dimension(:), intent(inout) :: vsig
        real*8, dimension(:), intent(out) :: vnle

 !****************************************************************************
        
        !--- Quantites globales
        real*8, dimension(size(vh,1),size(vh,2)) ::  P, Q
        real*8, dimension(dime,dime) :: V 
        real*8, dimension(dime) :: vsp 
        integer :: nc1, nc2, id

        !--- Quantites principales
        real*8, dimension(size(vh,1),size(vh,2)) :: vhloc
        logical :: vcrit
                
        ! Quantites principales
        real*8, dimension(size(vb,1)) :: vepspl, vsigm, vsigi, veps
        real*8, dimension(size(vh,1)) :: vsigmloc, vsigiloc, vepsloc
        real*8, dimension(size(vh,1)) :: vdfdsigloc, vdgdsigloc

        ! Variables por le modele d'endommagement 
        real*8 :: Dini, eps_0, eps_s, eps_l, D, sigma1, epsilon1, vcrit1, varcrit, eps_c, F_var, RT, C, W, E, nu, &
                 & eps_n, eps_t, tau, sigma, vcrit_ouv, vcrit_cis, D_t, eps_t_cri
                            
!********************************************************!
 
        ! Endommagement initial 
        Dini = endo(ie,ipg)
        eps_0 = varg0(ie,ipg) ! seuil initial
        eps_s = vargs(ie,ipg) ! seuil courant
        eps_l = 1000. ! deformation limite

        !----- Recuperation des parametres pour le critere
        id = 6
        if (dime == 3) id = id-1
        RT = vprel(id)    ! Resistance
        C  = vprel(id+1)  ! Cohesion
        eps_t_cri = vprel(id+2)  ! Deformation tangente critique
        
        if (C < 0.01) W = 0.d0

        E = vprel(id-2)
        nu = vprel(id-1)

        D = Dini

        ! Etat de contrainte et deformation
        vsigm = matmul((1.d0-D)*vh,(matmul(vb,vdle)))
        veps = (matmul(vb,vdle))

        V=reshape(irloc(ie,:),(/dime,dime/))
        !if (ipg==1) print*, ie, nomtype(ktypel(ie)), irloc(ie,:)
        !----- Recuperation de la matrice de changement de repere global --> local principal
        nc1 = size(vb,1)
        P = fiss_changement_repere(V,nc1,1)
        Q = fiss_changement_repere(V,nc1,2)
        vsigi = 0.d0
                  
        ! Dans le repere principal...
        vhloc = matmul(P,matmul(vh,Q))
        vsigmloc = matmul(P,vsigm)
        vsigiloc = matmul(P,vsigi)
        vepsloc = matmul(P,veps)

        !---- Mode rupture en ouverture pure ------
        sigma = vsigmloc(1)       ! Contrainte normale perpendiculaire avec l'acier
        eps_n = vepsloc(1)        ! Deformation dans la direction perpendiculaire avec l'acier
        vcrit_ouv = sigma - 1.0001*RT

        !----- Mode rupture en cisaillement -----
        tau = vsigmloc(2)         ! Contrainte tangentielle dans la direction de l'acier
        eps_t = abs(vepsloc(2))   ! Deformation suivante la direction de l'acier - Positive
        vcrit_cis = abs(tau) - 1.0001*C

        ! Choix entre les deux modes de rupture
        ! Verifier le critere au point central (pour la premiere fois de rupture)
        ! Fixer le mode de rupture 

        vcrit = .false.

        !if (ietat(ie) == 0) then
        !if (ie==iedng .and. ietatma(ie,ipg)==0) then
        if (ietatma(ie,ipg)==0) then
        
           if (vcrit_ouv > 0.d0) then
              !if (ie ==iedng) then
                ietatma(ie,ipg) = 1
                D = 1.d0

                if (ipg==1) then
                    ietatma(ie,:) = 1
                    endo(ie,:) = 1.d0
                end if
              !end if
              
           elseif (vcrit_cis > 0.d0) then
           !if (ie==iedng) then                
           
              vcrit = .true.              
              ietatma(ie,ipg) = 4    !-Changement etat de l'element

              eps_0 = eps_t
              eps_s = eps_0       ! Seuil de deformation pour tous les points Gauss.
              
              if (ipg==1) then
                 ietatma(ie,:) = 4   ! Point central                 
              end if
              
            !end if
           end if
           
        elseif (ietatma(ie,ipg)==1) then
           D = 1.d0

        elseif (ietatma(ie,ipg)==4) then
           vcrit = .true.
           
           if (eps_0==0.d0) eps_0 = varg0(ie,1)
           if (eps_s==0.d0) eps_s = vargs(ie,1)
           
        end if

        !----- Declanchement et evolution de l'endommagement
        if (vcrit) then

             !---- Valeur critique de la variable d'endommagement
             if (eps_t_cri < eps_0) eps_t_cri = eps_0

             if (eps_t >= eps_s) then

                  !---- Actualisation du seuil d'endommagement
                  eps_s = eps_t

                  !---- Loi d'evolution d'endommagement
                  if (eps_0 == eps_t_cri) then
                      D_t = 1.d0          ! Elastique fragile
                  else
                      D_t = (eps_t - eps_0)/(eps_c - eps_0 + 1.d-20)
                      if (D_t > 1.d0) D_t = 1.d0
                      if (D_t < 0.d0) D_t = 0.d0
                  end if

             end if

             !---- Calcul de l'endommagement (updated)
             D = 1.d0 - eps_0/eps_t * (1.d0 - D_t)

             if (D < 0.d0) D = 0.d0
             D = max(Dini,D)
             D = min(D,.9999d0)

             if ((eps_t > eps_t_cri) .or. (D .ge. 0.9999d0)) then
                  !if(eps_t > eps_c) pause
                  D=.9999d0
                  ietatma(ie,ipg) = 3
                  
                  if (ipg==1) then
                     ietatma(ie,:) = 3
                     endo(ie,:) = D
                  end if
             end if

        end if

        D = floor(1.d4 * D) * 1.d-4

        ! Actualisation de la contrainte et rigidite 
        vsig = matmul((1.d0-D)*vh,(matmul(vb,vdle)))
        vh = (1.d0-D)*vh

        endo(ie,ipg) = D
        varg0(ie,ipg) = eps_0
        vargs(ie,ipg) = eps_s
        
end subroutine interf_macro_endo


!----------------------------------------------------------------

subroutine interf_macro_pilot(imetpilo,vsol,vduI,vduII,dlam,alpham,MOTm,elefiss)

!********************************************************!
!  Pilotage du calcul : gestion du facteur de chargement !
!********************************************************!

        use variables, only : dime, nelt, idmax, ktypel, nomtype, infele, &
                    & vprelg, kprop, kprop0, ietat, iedngma, interf_macro, pi,     &
                    & mrtrav, vrtrav, inorm, irloc, kloce, ietatma
        use initialisation, only : init_vec, init_mat
        use lib_elem, only : elem_B, elem_kloce2, elem_hooke
        use fissuration, only : fiss_changement_repere
        use math
        use utilitaire, only : princ
        implicit none

        ! Variables IN
        real*8, dimension(:), intent(in) :: vduII, vsol
        integer, intent(in) :: imetpilo

        ! Variables IN-OUT
        real*8, dimension(:), intent(in) :: vduI
        real*8, intent(in) :: dlam 
        real*8, intent(out) :: alpham
        character(len=8), intent(inout) :: MOTm
        integer, intent(inout) :: elefiss

        ! Variables locales 
        real*8, dimension(:), allocatable :: vdle, vdl0, vsi0, vdsi, vsi, &
                                        & vsi0loc, vdsiloc, vdlI, vdlII, ksigPG
        real*8, dimension(:,:), allocatable :: vn, vb, vh, P        
        real*8 :: vprel(idmax), alph(nelt), vbid(dime), vbid2(1,dime*dime), V(dime,dime)
        real*8 :: alpc, sign0, C, RT, phi, psi, sigma, tau, &
                  & vcrit1, vcrit2, vcrit3, detj, signe, coef, alpg(1)
        character(len=5) :: typ
        character(len=8) :: MOT, rupt(nelt)

        integer :: i, j, k, ie, ipg, ipgd, npg, id, iloi, ndle
        
!----- Switch en fonction de la methode de pilotage

select case(imetpilo)
  
    case(1)
        !----- Pilotage sur element le plus dangereux par recalcul de vduI

        !- Pour le cas du modele d'interface macro
        if (interf_macro==1) then

            !- Initialisations
            alph = 1.d0
            iedngma = 0

            !- Boucle sur les elements
            do ie = 1, nelt
            
               MOT = '  '
               typ = nomtype(ktypel(ie))
            
               !----- Pour les elements d'interface macro vierges
               if ((typ=='MBI4') .and. (all(ietatma(ie,:)==0))) then

                  !----- Proprietes elementaires
                  vprel = vprelg(kprop(ie),1:idmax)

                  !----- Recuperation des parametres de la loi
                  iloi = int(vprel(1))

                  if (iloi/=201) goto 10000

                  alpg = 1.d0

                  do ipg = 1, 1   ! Verifier les contraintes au centre de l'element

                    !----- Vecteur kloce de localisation dans vecteur global
                    call elem_kloce2( ie, ndle)
                    allocate(vdl0(ndle));  allocate(vdle(ndle))
                    allocate(vdlI(ndle));  allocate(vdlII(ndle))

                    !----- Deplacements             
                    vdl0  = vsol(kloce(1:ndle))
                    vdlI  = vduI(kloce(1:ndle))
                    vdlII = vduII(kloce(1:ndle))
                    vdle  = vdlI + dlam*vdlII

                    !----- Matrice d'elasticite vh, en fonction du type d'element
                    call elem_hooke(vh,nomtype(ktypel(ie)),vprel)

                    !----- Calcul des contraintes au centre de gravite de l'element
                    call elem_B(vn,vb,detj,infele(ktypel(ie))%Q(ipg,:),ie)

                    ! dans le repere global
                    allocate(vsi0(size(vb,1)))
                    allocate(vdsi(size(vb,1)))
                    vsi0  = matmul(vh,matmul(vb,vdl0))
                    vdsi  = matmul(vh,matmul(vb,vdle))

                    allocate(vsi(size(vsi0,1)))
                    vsi = vsi0+vdsi
                    
                    ! Recuperation du repere de l'acier
                    V=reshape(irloc(ie,:),(/dime,dime/))

                    !----- Recuperation de la matrice de changement de repere global --> local principal (repere de l'acier)
                    call init_mat(P,size(vh,1),size(vh,2))
                    P = fiss_changement_repere(V,size(vb,1),1)
                              
                    ! Dans le repere principal...
                    allocate(vdsiloc(size(vh,1)))
                    allocate(vsi0loc(size(vh,1)))
                    vdsiloc = matmul(P,vdsi)
                    vsi0loc = matmul(P,vsi0)
                                      
                    deallocate(vdl0,vdlI,vdlII,vdle,vh,vb,vn)
                    deallocate(vsi0,vdsi,vsi,P)

                    if (iloi==201) then

                       !----- Recuperation des parametres de la loi
                       id = 6
                       if (dime == 3) id = id-1
                                 
                       RT = 1.001 * vprel(id)
                       C  = 1.001 * vprel(id+1)

                       !----- Recuperation de la contrainte normale et tangente
                       sigma = vsi0loc(1) + vdsiloc(1)
                       tau   = vsi0loc(2) + vdsiloc(2)

                       !----- Critere de rupture
                       vcrit1 = sigma - RT
                       vcrit2 = abs(tau) - C

                       if (vcrit1 > 0.d0) then            !----- Critere de rupture en ouverture pure
                           !- Critere en traction pure
                           MOT = 'ouvert'
                           alpg(ipg) = (RT - vsi0loc(1)) / vdsiloc(1)
                           !print*, sigma, RT, tau, C, ie

                       elseif (vcrit2 > 0.d0) then        !----- Critere de rupture en cisaillement
                           !- Critere en cisaillement
                           MOT = 'cisaille'
                           signe = tau/abs(tau)
                           alpg(ipg) = (C - signe*vsi0loc(2)) / (signe*vdsiloc(2))
                           !print*, sigma, RT, tau, C
                       end if

                    else
                       stop 'FIDES_fiss_pilot : loi non encore programmee'
                    end if

                    deallocate(vsi0loc,vdsiloc)
                    
                  end do

                  !----- On conserve 
                  rupt(ie) = MOT
                  alph(ie) = minval(alpg)

           10000   continue                                    
               end if 
            end do   
            ! Fin de la boucle sur les elements

            alpham = minval(alph)
         
            if (alpham < 1.d0) then
              !----- On stocke le numero de l'element le plus dangereux
              iedngma = find_num(alph,alpham)               
              elefiss = count(alph<=0)
              MOTm = rupt(iedngma)
           else
              alpham = 1.d0
           end if
        end if
        
    case(2)
        !----- Pilotage sur l'increment de contrainte avec recalcul de "dlam"
        
    case default
        stop 'FIDES_interf_macro_pilot : cas non encore programme'
    end select

end subroutine interf_macro_pilot


!********************************************************!


subroutine interf_macro_init()

!********************************************************!
!  Initialisations liees a la gestion des etats de       !
!               elements d'interface                     !
!********************************************************!
    
    use variables, only : dime, nelt, infele, iedngma,ietatma, interf_macro, &
           & kconec, ktypel, kprop, vcor, inorm, irloc, Emax, ntypel, vprelg, imacro, nomtype, &
           & varg0, vargs, endo, interf

    use initialisation
    use math
    implicit none

    integer :: id, ie, i, j, k, l, m , n, ncote, icote, ino, nbg, nbel, noel, npgm, itype, ntype
    integer, dimension(:,:), allocatable :: cote 
    integer, dimension(:), allocatable :: noec
    real*8, dimension(dime) :: vt, vn, cengr, cengr_r
    real*8, dimension(:,:), allocatable :: vcorn
    character(len=5) :: typel

!********************************************************!

    if (interf==0 .and. interf_macro==1)then

       ntype = maxval(ktypel)
       
       do itype = 1, ntype
         if (nomtype(itype)=='MBI4') then
           npgm = size(infele(itype)%W,1)
         end if 
       end do
       
       iedngma = 0
       call init_mat(endo,nelt,npgm)
       call init_mat(varg0,nelt,npgm)
       call init_mat(vargs,nelt,npgm)
       call init_mat(ietatma,nelt,npgm)
       endo = 0.d0
       varg0 = 0.d0
       vargs = 0.d0
    
       id = 4
       if (dime==3) id = id - 1
       nbg = maxval(kprop)

       ! Chercher les elements voisins d'un element         
       do ie = 1, nelt
         
         if (vprelg(kprop(ie),id) < 0.1) then                         ! Element poutre
            allocate(cote(size(infele(ktypel(ie))%face,1),size(infele(ktypel(ie))%face,2)))
            cote = infele(ktypel(ie))%face
            ncote = size(cote,1)

            allocate(noec(size(cote,2)))
            
            nbel = 0
       
            do icote = 1, ncote
              do i = 1, nelt
                noec = 0
                ino  = 0
                if (i/=ie .and. nomtype(ktypel(i))=='MBI4') then
                  k = 0 ; l = 0 ; m = 0 ; n = 0
                  do j = 1, size(kconec,2)
                    if (kconec(i,j)==kconec(ie,cote(icote,1))) then
                       k = 1
                       ino = ino + 1
                       noec(ino) = kconec(i,j)
                    end if
                     
                    if (kconec(i,j)==kconec(ie,cote(icote,2))) then
                       l = 1
                       ino = ino + 1
                       noec(ino) = kconec(i,j)
                    end if
                  
                    if (size(cote,2)==3) then
                       if (kconec(i,j)==kconec(ie,cote(icote,3))) m = 1
                    end if
             
                    if (size(cote,2)==4) then
                       if (kconec(i,j)==kconec(ie,cote(icote,3))) m = 1
                       if (kconec(i,j)==kconec(ie,cote(icote,4))) n = 1
                    end if
                  end do
         
                  if ((k==1.and.l==1.and.size(cote,2)==2).or.(k==1.and.l==1.and.m==1.and.size(cote,2)==3) &
                      & .or. (k==1 .and. l==1 .and. m==1 .and. n==1 .and. size(cote,2)==4)) then
                      
                      nbel = nbel + 1
                      
                      ! Changer les proprietes pour cet element
                      if (imacro .eqv. .false.) then  ! Si le groupe d'interface n'est pas encore declare
                         kprop(i)  = nbg+1            ! Loi interface macro, pour ces elements, on ne calcule plus la normale
                         ktypel(i) = ntypel           ! Changer le type d'element
                      end if

                      ! Definir la normale de cet element
                      ! (la vecteur normale doit etre orientee vers l'element de beton)

                      !         *-----------*-------*-----------*
                      !         |           |       |           |
                      !         |  Beton    |       |  Beton    |
                      !         |      <----| Acier |---->      |
                      !         |       n   |       |  n        |
                      !         |           |       |           |
                      !         *-----------*-------*-----------*

                      ! Calculer le centre de gravite de l'element
                      allocate(vcorn(dime,infele(ktypel(i))%nnel))
                      noel = infele(ktypel(i))%nnel

                      vcorn = vcor(:,kconec(i,1:noel))
                      cengr = sum(vcorn,2) / size(vcorn,2)                   ! centre d'element de massif
                      cengr_r = sum(vcor(:,noec),2) / size(vcor(:,noec),2)   ! centre d'element d'acier

                      vn = (cengr - cengr_r) / norme(cengr - cengr_r)
                      vt = (/ vn(2), -vn(1)/)

                      inorm(i,:) = vn
                      irloc(i,:) = (/vn, vt/)
                      !irloc(i,:) = (/vt, vn/)

                      deallocate(vcorn)
                      if (nbel == 2) exit
                  end if
                end if
              end do
           end do

           deallocate(cote,noec)

         
         
         ! Construire la liste d'element entourne l'acier
         elseif (vprelg(kprop(ie),id) >= Emax) then                   ! Element massif

            allocate(cote(size(infele(ktypel(ie))%face,1),size(infele(ktypel(ie))%face,2)))
            cote = infele(ktypel(ie))%face
            ncote = size(cote,1)

            allocate(noec(size(cote,2)))

            do icote = 1, ncote
              do i = 1, nelt
                noec = 0
                ino  = 0
                if (i/=ie .and. nomtype(ktypel(i))=='MBI4') then
                  k = 0 ; l = 0 ; m = 0 ; n = 0
                  do j = 1, size(kconec,2)
                    if (kconec(i,j)==kconec(ie,cote(icote,1))) then
                       k = 1
                       ino = ino + 1
                       noec(ino) = kconec(i,j)
                    end if
                     
                    if (kconec(i,j)==kconec(ie,cote(icote,2))) then
                       l = 1
                       ino = ino + 1
                       noec(ino) = kconec(i,j)
                    end if
                  
                    if (size(cote,2)==3) then
                       if (kconec(i,j)==kconec(ie,cote(icote,3))) m = 1
                    end if
             
                    if (size(cote,2)==4) then
                       if (kconec(i,j)==kconec(ie,cote(icote,3))) m = 1
                       if (kconec(i,j)==kconec(ie,cote(icote,4))) n = 1
                    end if
                  end do
         
                  if ((k==1.and.l==1.and.size(cote,2)==2).or.(k==1.and.l==1.and.m==1.and.size(cote,2)==3) &
                      & .or. (k==1 .and. l==1 .and. m==1 .and. n==1 .and. size(cote,2)==4)) then

                      ! Changer les proprietes pour cet element
                      if (imacro .eqv. .false.) then  ! Si le groupe d'interface n'est pas encore declare
                         kprop(i)  = nbg+1            ! Loi interface macro, pour ces elements, on ne calcule plus la normale
                         ktypel(i) = ntypel           ! Changer le type d'element
                      end if


                      ! Definir la normale de cet element
                      ! (la vecteur normale doit etre orientee vers l'element de beton)

                      !         *-----------*-------*-----------*
                      !         |           |       |           |
                      !         |  Beton    |       |  Beton    |
                      !         |      <----| Acier |---->      |
                      !         |       n   |       |  n        |
                      !         |           |       |           |
                      !         *-----------*-------*-----------*

                      ! Calculer le centre de gravite de l'element
                      allocate(vcorn(dime,infele(ktypel(i))%nnel))
                      noel = infele(ktypel(i))%nnel

                      vcorn = vcor(:,kconec(i,1:noel))
                      cengr = sum(vcorn,2) / size(vcorn,2)                   ! centre d'element de massif
                      cengr_r = sum(vcor(:,noec),2) / size(vcor(:,noec),2)   ! centre d'element d'acier

                      vn = (cengr - cengr_r) / norme(cengr - cengr_r)
                      vt = (/ vn(2), -vn(1)/)

                      inorm(i,:) = vn
                      irloc(i,:) = (/vn, vt/)
                      !irloc(i,:) = (/vt, vn/)

                      deallocate(vcorn)
                      
                      exit
                  end if
                end if
              end do
           end do

           deallocate(cote,noec)
         end if
      end do

    end if

end subroutine interf_macro_init


end module interface_macro

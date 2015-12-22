module depouil_lance

contains

subroutine depouill(nomfichier)

   use variables
   use post_traitement
   use initialisation
   implicit none

   character(len = *), intent(in) :: nomfichier
   character(len = len_trim(nomfichier)) :: nomsauve
   character(len = len_trim(nomfichier)+2) :: nomsauve1
   character(len = 100) :: nomlist
   character(len = len_trim(nomfichier)+3) :: nompoui
   
   real*8, dimension(:), allocatable :: C0
   real*8 :: ouvmin00, dp(3), fo(3), criouvmin(3)
 
   integer :: N, ncal, ical, i, dirimp, opt, numfich=1
      
   !----- Entree des donnees
   criouvmin = (/100, 200, 300/)
     
   N = 4
   allocate(C0(N))
   C0 = (/5, 5, 5, 5/)

   !----- Data -------------------
   dirimp = 2
   ncal = 1          
   opt = 2   ! 1 ---- tirant,  2 ----- poutre
      
   !----- Fin de partie des donnees
   allocate(resu(npas))
   do i = 1, npas
      call init_vec(resu(i)%vsolu,ndlt)
      call init_vec(resu(i)%vcontr,size(vcont))
      call init_vec(resu(i)%vnolin,size(vnoli))
      call init_mat(resu(i)%etatpg,nelt,size(ietatpg,2))
      call init_vec(resu(i)%etat,nelt)
      call init_mat(resu(i)%inorme,nelt,dime)
      call init_vec(resu(i)%dpglobal,3)
      call init_vec(resu(i)%foglobal,3)              
   end do  
           
   do ical = 1, ncal
        
      print*,' '
      print*,'Calcul : ', ical
        
      if (ncal==1) then
         nomsauve = nomfichier(1:len_trim(nomfichier)-5)//'.list'
         nomlist = dirresu//nomsauve
      end if
        
      if (ncal>1) then
         call filelist(nomfichier,nomsauve1,ical)
         nomlist = dirresu//nomfichier(1:len_trim(nomfichier)-5)//nomsauve1
      end if
         
      print*, nomlist
                     
      call lecture_list(nomlist)
               
      do ipas = 1, npas        
          dp(ipas) = resu(ipas)%dpglobal(dirimp)  ;   fo(ipas) = resu(ipas)%foglobal(dirimp)
      end do
            
      call depl_moy(opt)

      nompoui = nomfichier(1:len_trim(nomfichier)-5)//'.depouil'       
      open(unit = numfich, file = dirresu//nompoui, form = 'formatted', status = 'replace')
      
      if (opt==1) then
         write(numfich,'(10x, a3, 10x, a6, 10x, a7, 10x, a7, 10x, a6, 10x, a7, 10x, a7)') &
                  & 'Pas', 'Nb inf', 'Ouv inf', 'Esp inf', 'Nb sup', 'Ouv sup', 'Esp sup'
      else
         write(numfich,'(10x, a3, 10x, a6, 10x, a6, 10x, a03)') 'Pas', 'Nombre', 'Ouvmoy', 'Esp'
      end if
      
      do i = 1, size(criouvmin)
          ouvmin00 = criouvmin(i) * 1.e-6
          print*, 'Ouvmin = ', criouvmin(i)
          write(numfich,'(a6,3x,f8.2 )') 'Ouvmin', criouvmin(i)
          if (fiss==1) then
             if (dime==2) then
             
                if (opt==1) then
                    call detect_fissures_macro_tir(ouvmin00,2)
                else
                    call detect_fissures_macro_poutre(ouvmin00,2)
                end if
                 
                 
             elseif (dime==3) then
                 !call detect_fissures_macro_3d(fissure,ouvmin00)
             end if
          else
             !call detect_fissures_new(fissure,ouvmin00)
          end if      
      end do

     close(numfich)

   end do
     
   ! Depouillement
   print*, ' '
   print*, 'Fichier sauvegarde", Depouillement en cours ... '

   !call depouill
   
   do i = 1, npas
      deallocate(resu(i)%vsolu,resu(i)%vcontr,resu(i)%vnolin,resu(i)%etatpg,resu(i)%etat)
      deallocate(resu(i)%inorme,resu(i)%dpglobal,resu(i)%foglobal)
   end do           
   deallocate(resu)
   
end subroutine depouill


subroutine filelist(nomfichier,nomsauve,ical)

   implicit none
   integer, intent(in) :: ical
   character(len = *), intent(in) :: nomfichier
   character(len = len_trim(nomfichier)+2) :: nomsauve
    
   if (ical==1) nomsauve = nomfichier(1:len_trim(nomfichier)-5)//'-1.list'
   if (ical==2) nomsauve = nomfichier(1:len_trim(nomfichier)-5)//'-2.list'
   if (ical==3) nomsauve = nomfichier(1:len_trim(nomfichier)-5)//'-3.list'
   if (ical==4) nomsauve = nomfichier(1:len_trim(nomfichier)-5)//'-4.list'
   if (ical==5) nomsauve = nomfichier(1:len_trim(nomfichier)-5)//'-5.list'
   if (ical==6) nomsauve = nomfichier(1:len_trim(nomfichier)-5)//'-6.list'
   if (ical==7) nomsauve = nomfichier(1:len_trim(nomfichier)-5)//'-7.list'
   if (ical==8) nomsauve = nomfichier(1:len_trim(nomfichier)-5)//'-8.list'
   if (ical==9) nomsauve = nomfichier(1:len_trim(nomfichier)-5)//'-9.list'

end subroutine filelist

! ---------------------------------------------------

subroutine depl_moy(opt)

  use variables, only : dime, npas

  implicit none
  integer, intent(in) :: opt
  real*8, dimension(npas,dime) :: dep1s, dep2s, dep1i, dep2i
  real*8, dimension(npas) :: dep_sup, dep_inf, fle_moy, dep_moy
  real*8 :: x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4
  
  !-------------------------------------------------------------------------%
  !   Calculer le deplacement relatif entre les deux points sur le beton    %
  !-------------------------------------------------------------------------%

  
  if (dime==2) then

    if (opt==1) then
       !x1 = 0.1d0     ;    y1 = 0.d0
       !x2 = 1.6d0     ;    y2 = 0.d0
       !x3 = 0.1d0     ;    y3 = 0.1d0
       !x4 = 1.6d0     ;    y4 = 0.1d0
    
       x1 = 0.1d0     ;    y1 = 0.04825d0
       x2 = 1.6d0     ;    y2 = 0.04825d0
       x3 = 0.1d0     ;    y3 = -0.05175d0
       x4 = 1.6d0     ;    y4 = -0.05175d0
       
       call trace_depl(dep1s,x1,y1,0.d0,0.d0)
       call trace_depl(dep2s,x2,y2,0.d0,0.d0)
       call trace_depl(dep1i,x3,y3,0.d0,0.d0)
       call trace_depl(dep2i,x4,y4,0.d0,0.d0)
    
    elseif (opt==2) then
    
       x1 = 1.65d0    ;    y1 =  0.16d0
       x2 = 1.65d0    ;    y2 =  0.00d0
   
       call trace_depl(dep1s,x1,y1,0.d0,0.d0)
       call trace_depl(dep2s,x2,y2,0.d0,0.d0)       
    end if
    
  else
    stop ' Pas encore implante'
    
    z1 = 0.d0      ;    z2 = 0.d0
  
    call trace_depl(dep1s,x1,y1,z1,0.d0)
    call trace_depl(dep2s,x2,y2,z2,0.d0)    
    call trace_depl(dep1i,x3,y3,z3,0.d0)
    call trace_depl(dep2i,x4,y4,z4,0.d0)
  end if

  ! Tirant
  if (opt==1) then
    dep_sup = abs(dep2s(:,1)-dep1s(:,1))
    dep_inf = abs(dep2i(:,1)-dep1i(:,1))
    dep_moy = (dep_sup + dep_inf) / 2.
  end if
  
  ! Poutre
  if (opt==2) fle_moy = (dep1s(:,2) + dep2s(:,2))/2.

end subroutine depl_moy

! ---------------------------------------------------

subroutine trace_depl(dep,x,y,z,t)

! Trace des deplacements en un point quelconque du maillage
!       [dep]=trace_depl(resu,x,y,z)
!
! Entrees :
!   x, y, z : coordonnees du point considere (1D, 2D ou 3D)
!   itrac   : option de trace graphique (si=1)
!
! Sorties :
!   dep(:,1) : deplacements selon x pour tous les pas de calcul
!   dep(:,2) : deplacements selon y pour tous les pas de calcul
!

   use variables, only : vcor, kconec, ktypel, nomtype, nelt, dime, kloce, npas, infele, vsol, resu
   use math
   use lib_elem, only : elem_kloce2, elem_interp
   implicit none
   
   real*8, dimension(:,:), allocatable :: N, dNdxi, xyz, xni
   real*8, dimension(:), allocatable :: vdle
   real*8, dimension(dime) :: vpt, ksi, dksi, R
   real*8, dimension(dime) :: x0
   real*8, dimension(dime,1) :: x1
   real*8, dimension(dime,dime) :: dR

   real*8, dimension(npas,dime), intent(inout) :: dep
   real*8, intent(in) :: x, y, z, t
   real*8 :: n1, no
   
   integer :: ie, ele, ndle, ipas, iter, iconv, noel, inh
   character(len=4) :: typel
   
!-----------------------------------------
   if (dime==1) then
      vpt = (/x/)
   elseif (dime==2) then
      vpt = (/x , y /)
   else
      vpt = (/x, y, z /)
   end if
   
   inh = 0
   
   do ie = 1, nelt
      noel = infele(ktypel(ie))%nnel
      allocate(xyz(noel,dime))
      xyz = transpose(vcor(:,kconec(ie,1:noel)))
      
      if ((xyz(1,1)>vpt(1) .and. xyz(2,1)>vpt(1) .and. xyz(3,1)>vpt(1))) goto 100
      if ((xyz(1,1)<vpt(1) .and. xyz(2,1)<vpt(1) .and. xyz(3,1)<vpt(1))) goto 100
    
      if ((xyz(1,2)>vpt(2) .and. xyz(2,2)>vpt(2) .and. xyz(3,2)>vpt(2))) goto 100
      if ((xyz(1,2)<vpt(2) .and. xyz(2,2)<vpt(2) .and. xyz(3,2)<vpt(2))) goto 100
    
      if (dime==3) then
         if((xyz(1,3)>vpt(3) .and. xyz(2,3)>vpt(3) .and. xyz(3,3)>vpt(3))) goto 100
         if((xyz(1,3)<vpt(3) .and. xyz(2,3)<vpt(3) .and. xyz(3,3)<vpt(3))) goto 100
      end if

      call inhull(vpt,xyz(1,:),xyz(2,:),xyz(3,:),inh)

      if (inh==1) then
         ele=ie
         exit
      end if
      
100   deallocate(xyz)
   end do

   typel = nomtype(ktypel(ele))

   ! Calcul des coordonnees reduites du point dans l'element
   ! le calcul propose ici est non lineaire.
   ksi = sum(xyz,1)/noel
   allocate(xni(dime,noel))
   xni = transpose(xyz)
   x0 =  vpt

   if (inh==1) deallocate(xyz)

   iconv = 1
   iter = 1 

   do while (iconv==1)
     call elem_interp(N, dNdxi, typel, ksi)    
     x1 = matmul(xni,N)    
     R =  x1(:,1) - x0
     
     if (iter==1) then
        n1   = norme(R)       
        dR   = matmul(xni,dNdxi)    
        dksi =-matmul(inv(dR),R)
        ksi  = ksi + dksi
     else
        no = norme(R)/n1
        if (no < 1.e-04) then
           iconv = 0
        else
           dR   = matmul(xni,dNdxi)    
           dksi =-matmul(inv(dR),R)
           ksi  = ksi + dksi
        end if
     end if
     iter = iter + 1
     
     deallocate(N, dNdxi)
   end do

   call elem_kloce2(ele,ndle)
   allocate(vdle(ndle))
        
   do ipas = 1, npas
       ! Recherche des deplacements des noeuds de l'element
       vsol = resu(ipas)%vsolu
       vdle(1:ndle) = vsol(kloce(1:ndle))
       
       ! Interpolation des deplacements nodaux au point considere      
       call elem_interp(N, dNdxi, typel, ksi, dime)
       dep(ipas,:) = matmul((/vdle/),N)
       deallocate(N, dNdxi)
   end do
   
   deallocate(vdle)

end subroutine trace_depl

! -------------------------------------------------------------

subroutine inhull(P,A,B,C,inh)	

  implicit none

  real*8, dimension(2) :: P,A,B,C,v0,v1,v2
  real*8 :: dot00, dot01, dot02, dot11, dot12, invDenom, u, v
  integer, intent(out) :: inh
  
  inh = 0
  ! Compute vectors        
  v0 = C - A
  v1 = B - A
  v2 = P - A

  ! Compute dot products
  dot00 = dot_product(v0, v0)
  dot01 = dot_product(v0, v1)
  dot02 = dot_product(v0, v2)
  dot11 = dot_product(v1, v1)
  dot12 = dot_product(v1, v2)

  ! Compute barycentric coordinates
  invDenom = 1.d0 / (dot00 * dot11 - dot01 * dot01)
  u = (dot11 * dot02 - dot01 * dot12) * invDenom
  v = (dot00 * dot12 - dot01 * dot02) * invDenom

  ! Check if point is in triangle
  if ((u >= 0) .and. (v >= 0) .and. (u + v <= 1)) inh = 1

end subroutine inhull


! -----------------------------------------------------------------------------------

subroutine detect_fissures_macro_tir(ouvmin0,dir)

   use variables, only : dime, vcor, kconec, infele, ktypel, nomtype, nelt, &
                        & resu, npas, nnt, vsol, ietat, inorm
   use utilitaire, only : fiss_ouv
   use initialisation
   implicit none
   
   integer, intent(in) :: dir
   real*8, intent(in) :: ouvmin0
   real*8 :: ouv(nelt), elesup(10000), eleinf(10000), mcor(4)
   integer :: k, m, n, i, j, ie, idir, noel, ino, ipas, nbfiss, numfich=1

   integer :: sup_totnb, inf_totnb, totnb
   real*8 :: sup_totouv, inf_totouv, sup_moyouv, inf_moyouv, sup_espmoy, inf_espmoy, totouv
   
   !-----------------------------------------------------
   i = 0
   j = 0

   ! Pour les tirants
   mcor = (/ maxval(vcor(1,:)) , maxval(vcor(2,:)) , minval(vcor(1,:)) , minval(vcor(2,:))  /)

   do ie = 1, nelt
       if (dime==2) then
           if (nomtype(ktypel(ie))=='MBT3') then
               noel = infele(ktypel(ie))%nnel    
               if (count(vcor(dir,kconec(ie,1:noel))==mcor(dir))==2) then
                   i = i + 1;     elesup(i) = ie
               end if
               if (count(vcor(dir,kconec(ie,1:noel))==mcor(dir+2))==2) then
                   j = j + 1;     eleinf(j) = ie
               end if
           end if
       else
           stop 'depouill, Cas non implante'
       end if
   end do

   idir = 1
   if (dir==1) idir = 2

   ! ----- Determiner le nombre de fissures et l'ouverture de ces fissures
   do ipas = 1, npas
       !---- Initialisation
      sup_totnb = 0
      inf_totnb = 0
       
      sup_totouv = 0.d0
      inf_totouv = 0.d0
    
      sup_moyouv = 0.d0
      inf_moyouv = 0.d0

      sup_espmoy = 1.5d0          
      inf_espmoy = 1.5d0
                                
      !---- Detecter les fissures de tous les elements
      vsol = resu(ipas)%vsolu
      ietat = resu(ipas)%etat
      inorm = resu(ipas)%inorme
      ouv = fiss_ouv()
      
      !----- Face inferieure -------------------------------------------------------
      nbfiss = count(ouv(eleinf) >= ouvmin0)
      
      if (nbfiss > 0) then

        inf_totnb = nbfiss
                      
        m = 0
        do k = 1, j
            if (ouv(eleinf(k)) >= ouvmin0) then
                m = m + 1
                inf_totouv = inf_totouv + ouv(eleinf(k))
                if (m==nbfiss) exit
            end if
        end do
        
        inf_espmoy = 1.5d0/(nbfiss+1)        
        inf_moyouv = inf_totouv / inf_totnb
      end if

      !----- Face superieure -------------------------------------------------------
      nbfiss = count(ouv(elesup) >= ouvmin0)
      
      if (nbfiss > 0) then

        sup_totnb = nbfiss
                      
        n = 0
        do k = 1, i
            if (ouv(elesup(k)) >= ouvmin0) then
                n = n + 1
                inf_totouv = inf_totouv + ouv(elesup(k))
                if (n==nbfiss) exit
            end if
        end do
        
        inf_espmoy = 1.5d0/(nbfiss+1)        
        sup_moyouv = sup_totouv / sup_totnb       
      end if
      
      !---- Les valeurs totales   
      totnb  = sup_totnb + inf_totnb
      totouv = sup_totouv + inf_totouv
      
      write(numfich,*) ipas, inf_totnb, sup_totnb, inf_moyouv, sup_moyouv, inf_espmoy, sup_espmoy
      
   end do

end subroutine detect_fissures_macro_tir


!----------------------------------------------------------------------

! -----------------------------------------------------------------------------------

subroutine detect_fissures_macro_poutre(ouvmin0,dir)

   use variables, only : dime, vcor, kconec, infele, ktypel, nomtype, nelt, &
                  & resu, npas, nnt, vsol, ietat, inorm, dirresu
   use utilitaire, only : fiss_ouv
   use initialisation
   implicit none
   
   integer, intent(in) :: dir
   real*8, intent(in) :: ouvmin0
   
   real*8 :: ouv(nelt), mcor(4), totouv, moyouv, espmoy
   integer :: i, j, m, ie, idir, noel, ino, ipas, nbfiss, eleinf(10000), numfich=1

   !-----------------------------------------------------
   i = 0

   ! correc poutre rond
   do ino = 1, nnt
       if (vcor(2,ino) > -1.e-4 .and. vcor(2,ino) < 1.e-4) vcor(2,ino) = 0.0d0
   end do

   ! Pour les poutres
   mcor = (/ maxval(vcor(1,:)) , 0.16d0 , minval(vcor(1,:)) , 0.d0  /)   
   eleinf = 0
   
   do ie = 1, nelt
       if (dime==2) then
           if (nomtype(ktypel(ie))=='MBT3') then
               noel = infele(ktypel(ie))%nnel    
               if (count(vcor(dir,kconec(ie,1:noel))==mcor(dir+2))==2) then
                   i = i + 1;     eleinf(i) = ie
               end if
           end if
       else
           stop 'depouill, Cas non implante'
       end if
   end do

   idir = 1
   if (dir==1) idir = 2

   ! ----- Determiner le nombre de fissures et l'ouverture de ces fissures
   do ipas = 1, npas
      !---- Initialisation
      nbfiss = 0
      totouv = 0.d0
      moyouv = 0.d0
      espmoy = 3.d0

      !---- Detecter les fissures de tous les elements
      vsol = resu(ipas)%vsolu
      ietat = resu(ipas)%etat
      inorm = resu(ipas)%inorme
      ouv = fiss_ouv()

      !----- Face inferieure -------------------------------------------------------
      nbfiss = count(ouv(eleinf) >= ouvmin0)

      if (nbfiss > 0) then

        m = 0
        do j = 1, i
            if (ouv(eleinf(j)) >= ouvmin0) then
                m = m + 1
                totouv = totouv + ouv(eleinf(j))
                
                if (m==nbfiss) exit
            end if
        end do    
               
        moyouv = totouv / nbfiss
        espmoy = 3.d0/(nbfiss+1)
        print*, moyouv
      end if

      write(numfich,*) ipas, nbfiss, moyouv, espmoy

   end do



end subroutine detect_fissures_macro_poutre

!----------------------------------------------------------------

end module depouil_lance

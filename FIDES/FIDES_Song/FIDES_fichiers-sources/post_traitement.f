module post_traitement

!******************************************************************************
!*  Programme simple de post-traitement des resultats *.list crees par FIDES  *
!******************************************************************************

contains

subroutine FIGID(nomfichier,nomlist)

   use variables
   use utilitaire
   use initialisation

   implicit none
   
! -----------------------------------------------------------------------------------!
        logical :: ismember
        character :: reply
        integer :: i, pas, pas1, pas2, dpas, opt, ino
        character(len = 50), intent(in) :: nomfichier, nomlist
        character(len = 50) ::  nompas, nom_pas
        real*8 :: Mcm, maxouv, ampl = 1.
        real*8, dimension(nelt) :: ouv, Vel
        real*8 :: vcor0(dime,nnt)

! -----------------------------------------------------------------------------------!
         ! Entete
           print*,'  '
           print*,'   ****************************************************************************'
           print*,'   *                                                                          *'
           print*,'   * Programme simple de post-traitement des resultats *.list crees par FIDES *'
           print*,'   *                                                                          *'
           print*,'   ****************************************************************************'
           print*,'  '

         ! Calcul du volume de tous elements      
           Vel = volume_elem()

           call ecriture_mesh(nomfichier) 

           !call det_kconec()

         ! Lecture du fichier *.list

           allocate(resu(npas))

           vcor0 = vcor

           do i = 1, npas
              call init_vec(resu(i)%vsolu,ndlt)
              call init_vec(resu(i)%vcontr,size(vcont))
              call init_vec(resu(i)%vnolin,size(vnoli))
              call init_mat(resu(i)%etatpg,nelt,size(ietatpg,2))
              if (interf_macro==1) call init_mat(resu(i)%etatma,nelt,size(ietatma,2))
              call init_vec(resu(i)%etat,nelt)
              call init_mat(resu(i)%inorme,nelt,dime)           ! correction de Josquin. 09/06/2011
           end do     
           
           print *,'   Le fichier *.list choisi est : ', dirresu//nomlist
           print *,'  '
           print *,'   Lecture du fichier *.list en cours ......... '
           print *,'  '
           call lecture_list(nomlist)
           
           ! Chercher l'ouverture maximale
           maxouv = 0.d0
           Mcm = 0.d0
           print*,' Traitement ... '        
           
           !interf_macro=1
           
           do i = 1, npas
           
             print*,'Pas numero ', i
             
             vsol = resu(i)%vsolu
             vcont= resu(i)%vcontr
             vnoli = resu(i)%vnolin

             if (fiss==1 .or. interf_macro==1) then
                ietat = resu(i)%etat
                if (interf_macro==1) ietatma = resu(i)%etatma
                inorm = resu(i)%inorme
         
                ouv = fiss_ouv()
                maxouv = maxval(ouv)                
                Mcm = maxval( (/ Mcm,maxouv/) )
                
                if (i==npas) print*, '   Ouverture maximale : ', Mcm
             end if
             
             call ecriture_res(nomfichier,i,Vel)
               
           end do

           print*,'  '

           !if (fiss==0) goto 600
           
           write (nom_pas,'(i3)') npas
           ismember = .true.
           
       100 do while (ismember)
       
              write(*,'(a,$)') '    Entrez le pas pour le post-traitement (Ctrl + C pour exit) [ max = '//trim(nom_pas)//' ]: '
              read(*,*) pas              
              opt = 1              
              write(*,'(a,$)') '    Option de traitement (default opt = 1 ) : '
              read(*,*) opt  
                       
              if (pas < 1) goto 200
              
              write (nompas,*) pas
              
              nompas = adjustl(nompas)
              
              ! Recuperation des informations du pas concerne
              vsol = resu(pas)%vsolu
              vcont = resu(pas)%vcontr
              vnoli = resu(pas)%vnolin
             
              if (fiss==1 .or. interf_macro==1) then
                 ietat = resu(pas)%etat
                 if (interf==1) ietatpg = resu(pas)%etatpg
                 if (interf_macro==1) ietatma = resu(pas)%etatma
                 inorm = resu(pas)%inorme
              end if
                         
              if (fiss==0 .and. interf==1) then
                 ietatpg = resu(pas)%etatpg
                 call ecriture_interf(nomfichier,nompas)
                 call ecriture_res_interf(nomfichier,nompas,pas,Vel)
              end if
              
              if (fiss==1 .or. interf_macro==1) then
                 if (opt ==1) then
                    ampl = 0.
                    write(*,'(a,$)') '    Valeur d''amplification (default = 1 ) : '
                    read(*,*) ampl
                    
                    do ino = 1, nnt
                       if (dime == 2) vcor(1:2,ino) = vcor(1:2,ino) + ampl*vsol(infnoe(ino)%dln(1:2))
                       if (dime == 3) vcor(1:3,ino) = vcor(1:3,ino) + ampl*vsol(infnoe(ino)%dln(1:3))
                    end do
                    
                    call ecriture_fiss(nomfichier,nompas,Mcm)
                    vcor = vcor0
                 else
                    call ecriture_fiss2(nomfichier,nompas,Mcm) 
                 end if
              end if

              print*,'  '
              print*, '   **********************************************************'
              print*, '==>* Ecriture du fichier *.gid.msh du pas ',pas, ' est finie '
              print*, '   **********************************************************'
              print*,'  '
         300  write(*,'(a,$)') '    Voulez-vous continuer (1-continuer, autre-quitter, Ctrl+C pour exit) :  '                      
              read(*,*) reply  
              print*,'   '        
            
              if (reply /= '1') then
                 ismember = .false.
                 print*,'   '
                 print*,'==> exit program'
                 print*,'   '
              elseif (reply == '1') then
                 ismember = .true.
              end if
              
              go to 100
              
         200  if (pas==0) then
               print*,'  '
               print*,'   Creer une serie de fichier  '
               write(*,'(a,$)') '    Entrez le premier pas de temps [ min = 1 ]:  '; read(*,*) pas1
               write(*,'(a,$)') '    Entrez le deuxieme pas de temps [ max = '//trim(nom_pas)//' ]  '; 
               read(*,*) pas2                
               print*,'  '
               
               do pas = pas1, pas2

                    write (nompas,*) pas
                    nompas = adjustl(nompas) 
                                 
                    ! Recuperation des informations du pas concerne
                    vsol = resu(pas)%vsolu
                    vcont = resu(pas)%vcontr
             
                    if (fiss==1 .or. interf_macro==1) then
                       ietat = resu(pas)%etat
                       inorm = resu(pas)%inorme
                    end if
                    
                    if (fiss==0 .and. interf==1) then
                       ietatpg = resu(pas)%etatpg
                       call ecriture_interf(nomfichier,nompas)
                       call ecriture_res_interf(nomfichier,nompas,pas,Vel)
                    end if

                  if (fiss==1 .or. interf_macro==1) then
                    if (opt ==1) then
                       call ecriture_fiss(nomfichier,nompas,Mcm) 
                    else
                       call ecriture_fiss2(nomfichier,nompas,Mcm) 
                    end if
                  end if
               end do

              elseif (pas < 0) then
                 print*, '   **********************************************'
                 print*, '   ** ', pas, ' ce n''est pas un pas correct     '
                 print*, '   **********************************************'  
              end if

              go to 300

           end do

     600   print*,

end subroutine FIGID
!--------------------------------------------------------------------------!    

subroutine lecture_list(nomlist)
!-------------------------------------------------------------------------!
!        Routine de lecture du fichier binaire de résultats au            !
!        dernier pas de calcul.                                           !
!-------------------------------------------------------------------------!
       use variables
       use initialisation

       character(len = *), intent(in) :: nomlist
       integer, parameter :: numfich = 1
       integer :: ierr, bid, err, iforma
       real*8 :: bid0
       character(len=50) :: CHAINE
       
       iforma = 0                 ! Formatted et unformatted

       if (iforma==1) then
           open(unit = numfich, file = dirresu//nomlist, form='formatted', status='old', action='read', iostat=ierr)
       else
           open(unit = numfich, file = dirresu//nomlist, form='unformatted', status='old', action='read', iostat=ierr)
       end if

       if (ierr/=0) then
           print*,"lecture_list : le fichier ",dirresu//nomlist," n'existe pas !"
           stop
       end if

       ipas = 0
       err = 0

       do while (err /= -1)      

         if (iforma==1) then
           !----------------- Travail sur les noeuds du maillage -------------------
           read(numfich, *, iostat = err), CHAINE
         
           if((CHAINE(1:3) == 'pas').and.(err/=-1)) then           
             ipas = ipas + 1
                     
             print*,'Lecture du pas',ipas

             read(numfich,*) CHAINE
             read(numfich,*) bid0, bid0, bid0
             read(numfich,*) bid0, bid0, bid0
           end if
         
           if (CHAINE=='Solution' .and. err/=-1)      read(numfich,*) resu(ipas)%vsolu         
           if (CHAINE=='contraintes' .and. err/=-1)   read(numfich,*) resu(ipas)%vcontr         
           !if (CHAINE=='deformations' .and. err/=-1) read(numfich,*) resu(ipas)%vnolin
         
           !----- Lecture specifique aux elements d'interface
           if (CHAINE=='etatpg' .and. err/=-1)        read(numfich,*) resu(ipas)%etatpg
         
           !----- Lecture specifique au calcul beton fissurant
           if (CHAINE=='etat' .and. err/=-1)          read(numfich,*) resu(ipas)%etat               
           if (CHAINE=='norm' .and. err/=-1)          read(numfich,*) resu(ipas)%inorme

         else
           !----------------- Travail sur les noeuds du maillage -------------------
           read(numfich, iostat = err), CHAINE

           if ((CHAINE(1:3) == 'pas').and.(err/=-1)) then           
             ipas = ipas + 1
             read(numfich) bid

             print*,'Lecture du pas',ipas

             read(numfich) bid0, bid0, bid0
             read(numfich) bid0, bid0, bid0

             read(numfich) resu(ipas)%vsolu
             read(numfich) resu(ipas)%vcontr         
             !read(numfich) resu(ipas)%vnolin
         
             !----- Lecture specifique aux elements d'interface
             if (interf==1) read(numfich) resu(ipas)%etatpg
             if (interf_macro==1) read(numfich) resu(ipas)%etatma

             !----- Lecture specifique au calcul beton fissurant
             if (fiss==1) read(numfich) resu(ipas)%etat               
             if (fiss==1 .or. interf_macro==1) read(numfich) resu(ipas)%inorme
           end if
         end if

      end do

 1000 close(numfich)

      npas = ipas
     
end subroutine lecture_list

!------------------------------------------------------------------------------------!
!   Ecriture du fichier maillage pour le traitement graphique par le logiciel GID    !
!------------------------------------------------------------------------------------!

    subroutine ecriture_mesh(nomfichier)

       use variables, only : dirresu, vcor, ktypel, kconec, nomtype ,dime, nelt, nnt, kprop, kprop0

       implicit none

       character(len = *), intent(in) :: nomfichier
       character(len = 4) :: typel
       character(len = len_trim(nomfichier)+3) :: nommshgid
       integer :: ino, ie, iek, i

       nommshgid = nomfichier(1:len_trim(nomfichier)-5)//'.gid.msh'
       open(unit = 11, file = dirresu//nommshgid, form = 'formatted', status = 'replace')
       
       write(11,*) 'MESH dimension 3 ElemType Quadrilateral Nnode 4'
       write(11,*) 'Coordinates'

       do ino = 1, nnt
          if (dime == 2) write(11,*) ino, vcor(1:2,ino), 0.d0
          if (dime == 3) write(11,*) ino, vcor(1:3,ino)
       end do

       write(11,*) 'end coordinates'
       
       iek = 0
       
       write(11,*) 'Elements'
       write(11,*) 'end elements' 
       
       do i=1, maxval(kprop)
       write(11,*) 'MESH dimension 3 ElemType Triangle Nnode 3'
       write(11,*) 'Elements'
       do ie = 1, nelt
         typel = nomtype(ktypel(ie))
         if (typel == 'MBT3') then         
            if (kprop(ie)==i) then
              iek = iek + 1
              write(11,*) iek, kconec(ie,1:3)
            end if
         end if
       end do
       write(11,*) 'end elements'
       end do
       
       do i=1, maxval(kprop)
       write(11,*) 'MESH dimension 3 ElemType Quadrilateral Nnode 4'
       write(11,*) 'Elements'  
       do ie = 1, nelt
         typel = nomtype(ktypel(ie))
         if (typel == 'MBI4' .or. typel == 'MBQ4') then
           if (kprop(ie)==i) then
            iek = iek + 1
            write(11,*) iek, kconec(ie,1:4)
           end if
         end if
       end do       
       write(11,*) 'end elements' 
       end do
       
       write(11,*) 'MESH dimension 3 ElemType Tetrahedra Nnode 4'  
       write(11,*) 'Elements'
       do ie = 1, nelt
         typel = nomtype(ktypel(ie))        
         if (typel == 'MTT4') then
            if (kprop0(ie)==1) then
              iek = iek + 1
              write(11,*) iek, kconec(ie,1:4)
            end if
         end if
       end do
       write(11,*) 'end elements'
       
       write(11,*) 'MESH dimension 3 ElemType Tetrahedra Nnode 4'  
       write(11,*) 'Elements'
       do ie = 1, nelt
         typel = nomtype(ktypel(ie))        
         if (typel == 'MTT4') then
            if (kprop0(ie)==2) then
              iek = iek + 1
              write(11,*) iek, kconec(ie,1:4)
            end if
         end if
       end do
       write(11,*) 'end elements'

       write(11,*) 'MESH dimension 3 ElemType Tetrahedra Nnode 4'  
       write(11,*) 'Elements'
       do ie = 1, nelt
         typel = nomtype(ktypel(ie))        
         if (typel == 'MTT4') then
            if (kprop0(ie)==3) then
              iek = iek + 1
              write(11,*) iek, kconec(ie,1:4)
            end if
         end if
       end do
       write(11,*) 'end elements'
       

       do i=1, maxval(kprop)
       write(11,*) 'MESH dimension 3 ElemType Triangle Nnode 6'
       write(11,*) 'Elements'
       do ie = 1, nelt
         typel = nomtype(ktypel(ie))
         if (typel == 'MBT6') then
            if (kprop(ie)==i) then
                iek = iek + 1
                write(11,*) iek, kconec(ie,1:6)
            end if
         end if
       end do
       write(11,*) 'end elements'
       end do        

       write(11,*) 'MESH dimension 3 ElemType Quadrilateral Nnode 4'
       write(11,*) 'Elements'
       do ie = 1, nelt
         typel = nomtype(ktypel(ie))
         if (typel == 'EJQ4') then      
            iek = iek + 1
            write(11,*) iek, kconec(ie,1:4)
         end if
       end do
       write(11,*) 'end elements' 
       
       write(11,*) 'MESH dimension 3 ElemType Prism Nnode 6'
       write(11,*) 'Elements'
       do ie = 1, nelt
         typel = nomtype(ktypel(ie))
         if (typel == 'EJQ6') then         
            iek = iek + 1
            write(11,*) iek, kconec(ie,1:6)
         end if
       end do
       write(11,*) 'end elements'  

       write(11,*) 'MESH dimension 3 ElemType Hexahedra Nnode 8'
       write(11,*) 'Elements'
       do ie = 1, nelt
         typel = nomtype(ktypel(ie))
         if (typel == 'MTH8') then         
            iek = iek + 1
            write(11,*) iek, kconec(ie,1:8)
         end if
       end do
       write(11,*) 'end elements'  

       write(11,*) 'MESH dimension 3 ElemType Hexahedra Nnode 8'
       write(11,*) 'Elements'
       do ie = 1, nelt
         typel = nomtype(ktypel(ie))
         if (typel == 'EJQ8') then         
            iek = iek + 1
            write(11,*) iek, kconec(ie,1:8)
         end if
       end do
       write(11,*) 'end elements'  
              
       write(11,*) 'MESH dimension 3 ElemType Linear Nnode 2'
       write(11,*) 'Elements'
       do ie = 1, nelt
         typel = nomtype(ktypel(ie))
         if ((typel == 'BAR2').or.(typel == 'BEF2').or.(typel == 'BEF3').or.(typel == 'BAF2')) then
            iek = iek + 1
            write(11,*) iek, kconec(ie,1:2)
         end if
       end do
       write(11,*) 'end elements'

       close(11)

    end subroutine ecriture_mesh
   
   
!------------------------------------------------------------------------------------!
!   Ecriture du fichier maillage pour le traitement graphique par le logiciel GID    !
!------------------------------------------------------------------------------------!

    subroutine ecriture_interf(nomfichier, nompas)

       use variables, only : dirresu, vcor, ktypel, kconec, nomtype ,dime, nelt, nnt, kprop, ietatpg

       implicit none

       character(len = *), intent(in) :: nomfichier, nompas
       character(len = 4) :: typel
       character(len = len_trim(nomfichier)+len_trim(nompas)+7) :: nommshgid
       integer :: ino, ie, iek,inm

       integer :: ibef2, ibaf2, imtt4, imbq4, imbt3, iejq4, &
           & iejq6, iejt6, ibar2, ibef3, imth8, imtp6, iejt8

       nommshgid = nomfichier(1:len_trim(nomfichier)-5)//'_pas'//trim(nompas)//'.gid.msh'
       open(unit = 11, file = dirresu//nommshgid, form = 'formatted', status = 'replace')

       write(11,*) 'MESH dimension 3 ElemType Quadrilateral Nnode 4'
       write(11,*) 'Coordinates'


       imtt4 = 0;  imbq4 = 0;  iejq4 = 0; iejq6 = 0; iejt6 = 0; imbt3 = 0;
       ibar2 = 0 ; ibef2 = 0; ibef3 = 0; ibaf2 = 0; imth8 = 0;  imtp6 = 0; iejt8 = 0
       do inm = 1, maxval(ktypel)
          if (nomtype(inm)=='MTT4') imtt4 = 1
          if (nomtype(inm)=='MBQ4') imbq4 = 1
          if (nomtype(inm)=='MBT3') imbt3 = 1          
          if (nomtype(inm)=='EJQ4') iejq4 = 1
          if (nomtype(inm)=='EJQ6') iejq6 = 1
          if (nomtype(inm)=='EJT6') iejt6 = 1
          if (nomtype(inm)=='BAR2') ibar2 = 1
          if (nomtype(inm)=='BEF2') ibef2 = 1
          if (nomtype(inm)=='BAF2') ibaf2 = 1
          if (nomtype(inm)=='BEF3') ibef3 = 1
          if (nomtype(inm)=='MTH8') imth8 = 1
          if (nomtype(inm)=='MTP6') imtp6 = 1
          if (nomtype(inm)=='EJT8') iejt8 = 1
       end do

       do ino = 1, nnt
          if (dime == 2) write(11,*) ino, vcor(1:2,ino), 0.d0
          if (dime == 3) write(11,*) ino, vcor(1:3,ino)
       end do

       write(11,*) 'end coordinates'
       
       iek = 0
       
       write(11,*) 'Elements'  
       do ie = 1, nelt
         typel = nomtype(ktypel(ie))        
         if (typel == 'MBQ4') then
            iek = iek + 1
            write(11,*) iek, kconec(ie,1:4)
         end if
       end do       
       write(11,*) 'end elements' 
       
       write(11,*) 'MESH dimension 3 ElemType Tetrahedra Nnode 4'  
       write(11,*) 'Elements'
       do ie = 1, nelt
         typel = nomtype(ktypel(ie))        
         if (typel == 'MTT4') then
            if (kprop(ie)==1) then
              iek = iek + 1
              write(11,*) iek, kconec(ie,1:4)
            end if
         end if
       end do
       write(11,*) 'end elements'
       
       write(11,*) 'MESH dimension 3 ElemType Tetrahedra Nnode 4'  
       write(11,*) 'Elements'
       do ie = 1, nelt
         typel = nomtype(ktypel(ie))        
         if (typel == 'MTT4') then
            if (kprop(ie)==2) then
              iek = iek + 1
              write(11,*) iek, kconec(ie,1:4)
            end if
         end if
       end do
       write(11,*) 'end elements'

       write(11,*) 'MESH dimension 3 ElemType Tetrahedra Nnode 4'  
       write(11,*) 'Elements'
       do ie = 1, nelt
         typel = nomtype(ktypel(ie))        
         if (typel == 'MTT4') then
            if (kprop(ie)==3) then
              iek = iek + 1
              write(11,*) iek, kconec(ie,1:4)
            end if
         end if
       end do
       write(11,*) 'end elements'
       
       write(11,*) 'MESH dimension 3 ElemType Triangle Nnode 3'
       write(11,*) 'Elements'
       do ie = 1, nelt
         typel = nomtype(ktypel(ie))
         if (typel == 'MBT3') then         
            if (kprop(ie)==1) then
              iek = iek + 1
              write(11,*) iek, kconec(ie,1:3)
            end if
         end if
       end do
       write(11,*) 'end elements'

      write(11,*) 'MESH dimension 3 ElemType Triangle Nnode 3'
       write(11,*) 'Elements'
       do ie = 1, nelt
         typel = nomtype(ktypel(ie))
         if (typel == 'MBT3') then         
            if (kprop(ie)==2) then
            iek = iek + 1
            write(11,*) iek, kconec(ie,1:3)
            end if
         end if
       end do
       write(11,*) 'end elements'

       write(11,*) 'MESH dimension 3 ElemType Triangle Nnode 3'
       write(11,*) 'Elements'
       do ie = 1, nelt
         typel = nomtype(ktypel(ie))
         if (typel == 'MBT3') then         
            if (kprop(ie)==3) then
              iek = iek + 1
              write(11,*) iek, kconec(ie,1:3)
            end if
         end if
       end do
       write(11,*) 'end elements'

       write(11,*) 'MESH dimension 3 ElemType Triangle Nnode 6'
       write(11,*) 'Elements'
       do ie = 1, nelt
         typel = nomtype(ktypel(ie))
         if (typel == 'MBT6') then         
                iek = iek + 1
                write(11,*) iek, kconec(ie,1:6)
         end if
       end do
       write(11,*) 'end elements'

       write(11,*) 'MESH dimension 3 ElemType Triangle Nnode 6'
       write(11,*) 'Elements'
       do ie = 1, nelt
         typel = nomtype(ktypel(ie))
         if (typel == 'MBT6') then         
            if (kprop(ie)==2) then
            iek = iek + 1
            write(11,*) iek, kconec(ie,1:6)
            end if
         end if
       end do
       write(11,*) 'end elements'
       
       write(11,*) 'MESH dimension 3 ElemType Triangle Nnode 6'
       write(11,*) 'Elements'
       do ie = 1, nelt
         typel = nomtype(ktypel(ie))
         if (typel == 'MBT6') then         
            if (kprop(ie)==3) then
            iek = iek + 1
            write(11,*) iek, kconec(ie,1:6)
            end if
         end if
       end do
       write(11,*) 'end elements'          


	 ! Tracer les elements d'interface
	 if (iejq4 == 1) then
	   if (count(ietatpg(:,1)==1)>0) then
	     write(11,*) 'MESH "Layer ouvert" dimension 3 ElemType Quadrilateral Nnode 4'
	     write(11,*) '# Color 0.000000 0.5000000 0.000000'
			   
	     write(11,*) 'Elements'
	     do ie = 1, nelt
	 	   typel = nomtype(ktypel(ie))        
	 	   if ((typel == 'EJQ4').and.(ietatpg(ie,1)==1)) then
	 		  iek = iek + 1
	 		  write(11,*) iek, kconec(ie,1:4)
	 	   end if
	     end do
	     write(11,*) 'end elements'
	   end if
	     
    	   if (count(ietatpg(:,1)==2)>0) then
	     write(11,*) 'MESH "Layer cisaille" dimension 3 ElemType Quadrilateral Nnode 4'
	     write(11,*) '# Color 0.000000 0.500000 0.500000'

	     write(11,*) 'Elements'
	     do ie = 1, nelt
	 	   typel = nomtype(ktypel(ie))        
	 	   if ((typel == 'EJQ4').and.(ietatpg(ie,1)==2)) then
	 		  iek = iek + 1
	 		  write(11,*) iek, kconec(ie,1:4)
	 	   end if
	     end do
	     write(11,*) 'end elements'
	   end if

        if (count(ietatpg(:,1)==3)>0) then
          write(11,*) 'MESH "Layer endomage" dimension 3 ElemType Quadrilateral Nnode 4'
	     write(11,*) '# Color 0.500000 0.500000 0.000000'
	     
          write(11,*) 'Elements'
          do ie = 1, nelt
            typel = nomtype(ktypel(ie))
            if ((typel == 'EJQ4') .and. (ietatpg(ie,1)==3)) then         
              iek = iek + 1
              write(11,*) iek, kconec(ie,1:4)
            end if
          end do
          write(11,*) 'end elements'
        end if
       
        if (count(ietatpg(:,1)==0)>0) then
          write(11,*) 'MESH "Layer vierge" dimension 3 ElemType Quadrilateral Nnode 4'
	     write(11,*) '# Color 0.500000 0.000000 0.500000'
          write(11,*) 'Elements'
          do ie = 1, nelt
            typel = nomtype(ktypel(ie))
            if ((typel == 'EJQ4') .and. (ietatpg(ie,1)==0)) then
              iek = iek + 1
              write(11,*) iek, kconec(ie,1:4)
            end if
          end do
          write(11,*) 'end elements'
        end if

	 end if

	 if (iejq6 == 1) then
       if (count(ietatpg(:,1)==1)>0) then
	     write(11,*) 'MESH "Layer ouvert" dimension 3 ElemType Quadrilateral Nnode 6'
	     write(11,*) '# Color 0.000000 0.000000 1.000000'
			   
	     write(11,*) 'Elements'
	     do ie = 1, nelt
	 	   typel = nomtype(ktypel(ie))        
	 	   if ((typel == 'EJQ6') .and. (ietatpg(ie,1)==1)) then
	 		   iek = iek + 1
	 		   write(11,*) iek, kconec(ie,1:6)
	 	   end if
	     end do
	     write(11,*) 'end elements'
       endif
           
       if (count(ietatpg(:,1)==2)>0) then
	     write(11,*) 'MESH "Layer cisaille" dimension 3 ElemType Quadrilateral Nnode 6'
	     write(11,*) '# Color 0.000000 0.500000 0.500000'

	     write(11,*) 'Elements'
	     do ie = 1, nelt
	 	   typel = nomtype(ktypel(ie))        
	 	   if ((typel == 'EJQ6') .and. (ietatpg(ie,1)==2)) then
	 		   iek = iek + 1
	 		   write(11,*) iek, kconec(ie,1:6)
	 	   end if
	     end do
	     write(11,*) 'end elements'
	   end if

        if (count(ietatpg(:,1)==3)>0) then
          write(11,*) 'MESH "Layer endomage" dimension 3 ElemType Quadrilateral Nnode 6'
	     write(11,*) '# Color 0.500000 0.500000 0.000000'
          write(11,*) 'Elements'
          do ie = 1, nelt
            typel = nomtype(ktypel(ie))
            if ((typel == 'EJQ6') .and. (ietatpg(ie,1)==3)) then
              iek = iek + 1
              write(11,*) iek, kconec(ie,1:6)
            end if
          end do
          write(11,*) 'end elements'
        end if

        if (count(ietatpg(:,1)==0)>0) then
          write(11,*) 'MESH "Layer vierge" dimension 3 ElemType Quadrilateral Nnode 6'
	     write(11,*) '# Color 0.500000 0.000000 0.500000'
          write(11,*) 'Elements'
          do ie = 1, nelt
            typel = nomtype(ktypel(ie))
            if ((typel == 'EJQ6') .and. (ietatpg(ie,1)==0)) then
              iek = iek + 1
              write(11,*) iek, kconec(ie,1:6)
            end if
          end do
          write(11,*) 'end elements'
        end if

      end if

	 ! fin des elements d'interface

       write(11,*) 'MESH dimension 3 ElemType Hexahedra Nnode 8'
       write(11,*) 'Elements'
       do ie = 1, nelt
         typel = nomtype(ktypel(ie))
         if (typel == 'MTH8') then
            iek = iek + 1
            write(11,*) iek, kconec(ie,1:8)
         end if
       end do
       write(11,*) 'end elements'  

       write(11,*) 'MESH dimension 3 ElemType Hexahedra Nnode 8'
       write(11,*) 'Elements'
       do ie = 1, nelt
         typel = nomtype(ktypel(ie))
         if (typel == 'EJQ8') then         
            iek = iek + 1
            write(11,*) iek, kconec(ie,1:8)
         end if
       end do
       write(11,*) 'end elements'  
              
       write(11,*) 'MESH dimension 3 ElemType Linear Nnode 2'
       write(11,*) 'Elements'
       do ie = 1, nelt
         typel = nomtype(ktypel(ie))
         if ((typel == 'BAR2').or.(typel == 'BEF2').or.(typel == 'BEF3').or.(typel == 'BAF2')) then
            iek = iek + 1
            write(11,*) iek, kconec(ie,1:2)
         end if
       end do
       write(11,*) 'end elements'

       close(11)

    end subroutine ecriture_interf
    
!------------------------------------------------------------------------------------!
!         Ecriture des résultats dans un fichier a chaque pas de calcul              !
!  Entrees :                                                                         !
!  - ipas : pas de calcul                                                            !      
!------------------------------------------------------------------------------------!

    subroutine ecriture_res(nomfichier,ipas,Vel)
    
       use variables, only : dirresu, vsol, vcont, vnoli, ndlt, nelt, infnoe, knddl, &
             & calco, fiss, interf, dime, nnt, ktypel, infele, nomtype, kprop, kconec, kcont, listelemautour
       use math
       use lib_elem, only : elem_extrap
       use initialisation
       
       implicit none

       character(len = *), intent(in) :: nomfichier
       integer, intent(in) :: ipas
       integer :: i, n, m, iopt, maxnoe, inm, ic
       integer :: ierr    
       character(len = len_trim(nomfichier)) :: nometude, nomsauve
       character(len = len_trim(nomfichier)+3) :: nomresgid
       
       real*8, dimension(nelt), intent(in) :: Vel
       
       integer :: ino, iel, ie, nc1, npg, idim1, iloc1, noel, ip, t, j, nip(nnt), liel(100)
       real*8, dimension(:,:), allocatable ::  ksig, ksin, vsig, vcont_no, vsigt
       real*8 :: Vet
       real*8, dimension(:), allocatable :: vn
                 
       nomresgid = nomfichier(1:len_trim(nomfichier)-5)//'.gid.res'

       if(ipas == 1) then        
          open(unit = 9, file = dirresu//nomresgid, form = 'formatted', status = 'replace')       
       else
          open(unit = 9, file = dirresu//nomresgid, form = 'formatted', status = 'old', position = 'append')                        
       end if
       
       maxnoe = maxval(knddl)    
       
       !-------- Enregistrement des deplacement ---------!
       write(9,'(a13,1x,a3,1x,i4,1x,a12)') 'Deplacement  ',' 1  ',ipas, '  2   1   1 '
       write(9,'(a13)') 'Deplacement x'
       write(9,'(a13)') 'Deplacement y'
                                         
       if (dime==2) then
           if (maxnoe==3) write(9,'(a8)') 'Rotation xy'
          
       elseif (dime==3) then              
          write(9,'(a13)') 'Deplacement z'
       end if       

       do i = 1, nnt       
          if (dime==2) then  
            if ((knddl(i)==2).and.(maxnoe==3)) then 
                write(9,*)  i, vsol(infnoe(i)%dln), 0.d0                
            else
                write(9,*)  i, vsol(infnoe(i)%dln)                
            end if          

          elseif (dime == 3) then 
              write(9,*) i, vsol(infnoe(i)%dln)      
          end if
       end do       
      
       !-------- Enregistrement des contraintes ---------!       
       if (dime==2) call init_mat(vcont_no,nnt,3)
       if (dime==3) call init_mat(vcont_no,nnt,6)

    liel = 0

    if ((fiss==1) .or. (interf==0 .and. fiss==0)) then
       do ino = 1, nnt

          liel(1:size(listelemautour(ino)%el)) = listelemautour(ino)%el
          t = size(listelemautour(ino)%el)
          if (dime==2) call init_mat(vsigt,size(listelemautour(ino)%el),3)
          if (dime==3) call init_mat(vsigt,size(listelemautour(ino)%el),6)

          Vet = 0.d0
          
          !----- Boucle sur les elements
          do iel = 1, t
            
            ie = liel(iel)
            
            !----- Recuperation des points d'integration
            call init_mat(ksig, size(infele(ktypel(ie))%Q,1), size(infele(ktypel(ie))%Q,2))            
            ksig = infele(ktypel(ie))%Q
            nc1  = infele(ktypel(ie))%ncgr     ! nombre de composantes du gradient
            
            npg = size(ksig,1)

            if (calco == 'GAUSS') then
                idim1 = npg*nc1
            else
                stop ' Calcul des vecteurs contrainte au point de Gauss '    
            end if
            
            iloc1 = kcont(ie)
            
            !----- Recuperation des vecteurs elementaires dans vcont
            call init_mat(vsig, nc1, npg)
            vsig  = reshape(vcont(iloc1 : (iloc1 + idim1 - 1)),(/ nc1, npg /))     
            
            vsigt(iel,:) = Vel(ie) * vsig(1:nc1,1)
            
            Vet = Vet + Vel(ie)

            deallocate(vsig, ksig)        
            
          end do

          ! Calcul de la moyenne des composantes extrapolees aux noeuds
          do ic = 1, nc1
              vcont_no(ino,ic) = sum(vsigt(:,ic))/Vet
          end do
             
          deallocate(vsigt)  
       end do 
       
    else

       do ie = 1, nelt
            
         if ((nomtype(ktypel(ie))/='EJQ4').and.(nomtype(ktypel(ie))/='EJQ6').and.(nomtype(ktypel(ie))/='EJT6')) then
          call init_mat(ksig, size(infele(ktypel(ie))%Q,1), size(infele(ktypel(ie))%Q,2))
          call init_mat(ksin, size(infele(ktypel(ie))%ksin,1), size(infele(ktypel(ie))%ksin,2))

          noel = infele(ktypel(ie))%nnel
          ksin = infele(ktypel(ie))%ksin
          ksig = infele(ktypel(ie))%Q
          nc1  = infele(ktypel(ie))%ncgr     ! nombre de composantes du gradient
            
          npg = size(ksig,1)
              
          if (calco == 'GAUSS') then
               idim1 = npg*nc1
          else
               stop ' Calcul des vecteurs contrainte au point de Gauss '
          end if   
                            
          !----- Recuperation des vecteurs elementaires dans vcont
          call init_mat(vsig, nc1, npg)
          iloc1 = kcont(ie)       

          vsig  = reshape(vcont(iloc1 : (iloc1 + idim1 - 1)),(/ nc1, npg /))
          
          do ino = 1, noel
          
             select case (nomtype(ktypel(ie)))
                case('MBT3','MTT4')
                    call elem_extrap(vn,1,ksin(ino,:))
                case('MBT6')
                    call elem_extrap(vn,3,ksin(ino,:))  
                case('MBQ4')
                    call elem_extrap(vn,4,ksin(ino,:))
                case('MTP6')
                    call elem_extrap(vn,6,ksin(ino,:))
                case('MTH8')
                    call elem_extrap(vn,8,ksin(ino,:))
	        case default
    	            stop 'Post_traitement : cas non encore implante'                  
             end select
                       
             vcont_no(kconec(ie,ino),:) = vcont_no(kconec(ie,ino),:) + matmul(vsig,vn)

             deallocate(vn)
            
          end do
            
          deallocate(vsig, ksig, ksin)
         end if                                       
       end do

       ! Calcul des occurences du noeud ip dans kconec
       do ip = 1, nnt

           ! Calcul de la moyenne des composantes extrapolees aux noeuds
           vcont_no(ip,:) = vcont_no(ip,:)/size(listelemautour(ip)%el)
       end do
                  
     end if
   
       !-------- Enregistrement des contraintes ---------!
       write(9,'(a13,1x,a3,1x,i4,1x,i1,a9)') 'Contraintes  ',' 1  ',ipas, dime,' 1   1 '
       write(9,'(a13)') 'Contraintes x'
       write(9,'(a13)') 'Contraintes y'
                                         
       if (dime==2) then
          write(9,'(a14)') 'Contraintes xy'
          
       elseif (dime==3) then              
          write(9,'(a13)') 'Contraintes z'
          write(9,'(a14)') 'Contraintes xy'
          write(9,'(a14)') 'Contraintes yz'
          write(9,'(a14)') 'Contraintes xz'
       end if       
             
       do i = 1, nnt
          if (dime==2) then             
              write(9,*)  i, vcont_no(i,:)
    
          elseif (dime == 3) then
              write(9,*) i, vcont_no(i,:)
          end if
       end do 
       
       deallocate(vcont_no)
           
       close(unit = 9)
       
    end subroutine ecriture_res

!------------------------------------------------------------------------------------!
!         Ecriture des résultats dans un fichier a chaque pas de calcul              !
!  Entrees :                                                                         !
!  - ipas : pas de calcul                                                            !      
!------------------------------------------------------------------------------------!

    subroutine ecriture_res_interf(nomfichier,nompas,ipas,Vel)
    
       use variables, only : dirresu, vsol, vcont, vnoli, ndlt, nelt, infnoe, knddl, &
             & calco, fiss, interf, dime, nnt, ktypel, infele, nomtype, kprop, kconec, kcont
       use math
       use lib_elem, only : elem_extrap
       use initialisation
       
       implicit none

       character(len = *), intent(in) :: nomfichier, nompas
       integer, intent(in) :: ipas
       integer :: i, j, t, n, m, iopt, maxnoe, inm, ic
       integer :: ierr
       character(len = len_trim(nomfichier)) :: nometude, nomsauve
       character(len = len_trim(nomfichier)+len_trim(nompas)+7) :: nomresgid


       real*8, dimension(nelt), intent(in) :: Vel
       integer :: nip(nnt), liel(100)
       
       integer :: ino, iel, ie, nc1, npg, idim1, iloc1, noel, ip
       real*8, dimension(:,:), allocatable ::  ksig, ksin, vsig, vcont_no, vsigt
       real*8 :: Vet, time2, time1
       real*8, dimension(:), allocatable :: vn
                 
       nomresgid = nomfichier(1:len_trim(nomfichier)-5)//'_pas'//trim(nompas)//'.gid.res'
       open(unit = 9, file = dirresu//nomresgid, form = 'formatted', status = 'replace')

       maxnoe = maxval(knddl)    
       
       !-------- Enregistrement des deplacement ---------!
       write(9,'(a13,1x,a3,1x,i4,1x,a12)') 'Deplacement  ',' 1  ',ipas, '  2   1   1 '
       write(9,'(a13)') 'Deplacement x'
       write(9,'(a13)') 'Deplacement y'
                                         
       if (dime==2) then
           if (maxnoe==3) write(9,'(a8)') 'Rotation xy'
          
       elseif (dime==3) then              
          write(9,'(a13)') 'Deplacement z'
       end if       

       do i = 1, nnt       
          if (dime==2) then  
            if ((knddl(i)==2).and.(maxnoe==3)) then 
                write(9,*)  i, vsol(infnoe(i)%dln), 0.d0                
            else
                write(9,*)  i, vsol(infnoe(i)%dln)                
            end if          

          elseif (dime == 3) then 
              write(9,*) i, vsol(infnoe(i)%dln)      
          end if
       end do       
      
       !-------- Enregistrement des contraintes ---------!
       if (dime==2) call init_mat(vcont_no,nnt,3)
       if (dime==3) call init_mat(vcont_no,nnt,6)

       liel = 0

       do ie = 1, nelt
            
         if ((nomtype(ktypel(ie))/='EJQ4').and.(nomtype(ktypel(ie))/='EJQ6').and.(nomtype(ktypel(ie))/='EJT6')) then
          call init_mat(ksig, size(infele(ktypel(ie))%Q,1), size(infele(ktypel(ie))%Q,2))
          call init_mat(ksin, size(infele(ktypel(ie))%ksin,1), size(infele(ktypel(ie))%ksin,2))

          noel = infele(ktypel(ie))%nnel
          ksin = infele(ktypel(ie))%ksin
          ksig = infele(ktypel(ie))%Q
          nc1  = infele(ktypel(ie))%ncgr     ! nombre de composantes du gradient
            
          npg = size(ksig,1)
              
          if (calco == 'GAUSS') then
               idim1 = npg*nc1
          else
               stop ' Calcul des vecteurs contrainte au point de Gauss '
          end if   
                            
          !----- Recuperation des vecteurs elementaires dans vcont
          call init_mat(vsig, nc1, npg)
          iloc1 = kcont(ie)       

          vsig  = reshape(vcont(iloc1 : (iloc1 + idim1 - 1)),(/ nc1, npg /))
          
          do ino = 1, noel
          
             select case (nomtype(ktypel(ie)))
                case('MBT3','MTT4')
                    call elem_extrap(vn,1,ksin(ino,:))
                case('MBT6')
                    call elem_extrap(vn,3,ksin(ino,:))  
                case('MBQ4')
                    call elem_extrap(vn,4,ksin(ino,:))
                case('MTP6')
                    call elem_extrap(vn,6,ksin(ino,:))
                case('MTH8')
                    call elem_extrap(vn,8,ksin(ino,:))
	        case default
    	            stop 'Post_traitement : cas non encore implante'                  
             end select
                       
             vcont_no(kconec(ie,ino),:) = vcont_no(kconec(ie,ino),:) + matmul(vsig,vn)

             deallocate(vn)
            
          end do
            
          deallocate(vsig, ksig, ksin)
         end if                                       
       end do

       ! Calcul des occurences du noeud ip dans kconec
       do ip = 1, nnt
           nip(ip) = count(kconec==ip)

           ! Calcul de la moyenne des composantes extrapolees aux noeuds
           vcont_no(ip,:) = vcont_no(ip,:)/nip(ip)
       end do

       !-------- Enregistrement des contraintes ---------!
       write(9,'(a13,1x,a3,1x,i4,1x,i1,a9)') 'Contraintes  ',' 1  ',ipas, dime,' 1   1 '
       write(9,'(a13)') 'Contraintes x'
       write(9,'(a13)') 'Contraintes y'
                                         
       if (dime==2) then
          write(9,'(a14)') 'Contraintes xy'
          
       elseif (dime==3) then              
          write(9,'(a13)') 'Contraintes z'
          write(9,'(a14)') 'Contraintes xy'
          write(9,'(a14)') 'Contraintes yz'
          write(9,'(a14)') 'Contraintes xz'
       end if       
             
       do i = 1, nnt
          if (dime==2) then             
              write(9,*)  i, vcont_no(i,:)
    
          elseif (dime == 3) then
              write(9,*) i, vcont_no(i,:)
          end if
       end do 
       
       deallocate(vcont_no)
           
       close(unit = 9)
       
    end subroutine ecriture_res_interf

!-----------------------------------------------------------------------------!
  
    function volume_elem() result(Vel)
    
       use variables, only : nelt, interf, fiss ,dime, ktypel, idmax, infele, kprop, vprelg
       use lib_elem, only : elem_B
       use initialisation
    
       implicit none

       real*8, dimension(nelt) :: Vel
       real*8, dimension(:,:), allocatable :: vn, vb, ksig
       real*8, dimension(:), allocatable ::  vpg

       real*8 :: detj, epais, poids, Ve
       integer :: i, ie, nbg, ipg, npg, id
          
       !--------------------------------------------------!
          
       do ie = 1, nelt   
                       
          !-- Calcul du volume de chaque element
          call init_vec(vpg, size(infele(ktypel(ie))%W))
          call init_mat(ksig, size(infele(ktypel(ie))%Q,1), size(infele(ktypel(ie))%Q,2))
          vpg = infele(ktypel(ie))%W
          ksig = infele(ktypel(ie))%Q
          npg = size(ksig,1)
                
          epais = 1.d0

          ! recuperation de l'epaisseur si 2D
          if (dime == 2) epais = vprelg(kprop(ie),idmax-2)
          
          Ve = 0
          do ipg = 1, npg
             poids = vpg(ipg)
             ! calcul des fonctions d'interpolation et des derivees 
             call elem_B(vn,vb,detj,ksig(ipg,:),ie)
             deallocate(vb,vn)
             Ve = Ve + detj*poids*epais
          end do

          Vel(ie) = Ve
                         
          deallocate(vpg, ksig)
         
       end do 
       
    end function volume_elem

end module post_traitement
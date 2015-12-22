module proprietes_aleatoires

!********************************************************!
!     Proprietes mecaniques aleatoire                 ! 
!********************************************************!

contains

!********************************************************!

subroutine alea_trac(eta,loi,interv,dPlay,mot)
!------------------------------------------------------------!
!       Fonction de trace d'une serie de valeurs tirees      !
!   aleatoirement selon une densite de probabilites donnee.  !
!------------------------------------------------------------!
    
    use variables, only : nelt
    use initialisation, only : init_vec
    use gnufor2

    integer, intent(in)  :: dPlay, interv
    integer :: ie, i, nnz, j, nEch, nEch_prec
    real*8 :: moy, ect, tmp, mx, mn, dm
    real*8, dimension(5) :: loi
    real*8, dimension(:), intent(in) :: eta
    real*8, dimension(:), allocatable :: etannz 
    real*8, dimension(interv+1) :: rgt, cumul, dist
    character(len=15) :: nom
    character(len=*), optional :: mot
    
    nnz = count(eta /= 0.)
    call init_vec(etannz,nnz)
    i = 1

    do ie = 1, nelt
       if (eta(ie) /= 0.) then
           etannz(i) = eta(ie)
           i = i + 1
       end if
    end do   

    if (dPlay /= 0) then
       if (loi(5) == 1) then
          nom = 'gaussienne'
       elseif (loi(5) == 2) then
          nom = 'de Weibull'
       elseif (loi(5) == 3) then
          nom = 'log-normale'
       else
          nom = 'inconnue'
       end if

       moy = (sum(etannz)/size(etannz))

       tmp = 0.
       do i = 1,size(etannz)
        tmp = tmp + (etannz(i) - moy)**2 
       end do

       if (size(etannz)==1) then
        ect = sqrt(tmp)
       else
               ect = sqrt(tmp/(size(etannz)-1))
       end if

       mx = maxval(etannz)
       mn = minval(etannz)
       dm = (mx - mn)/interv

       do i = 0, interv
           rgt(i+1) = mn + dm*i
       end do

       nEch_prec = 0
       nEch = 0

       do i = 1,(interv+1)
         nEch = count((etannz <= rgt(i)) .eqv. .true.)
         cumul(i) = nEch
         dist(i) = nEch - nEch_prec
         nEch_prec = nEch
       end do
 
       cumul = cumul/nEch
       dist = dist/nEch

       if (dPlay>= 1) then
           print*,
           print*,'===> ',mot,' :'
           print'(a11,e12.5,1x,a19,e12.5)','Moyenne : ', moy ,' -   Ecart-type : ',  ect
           print'(a13,a13,a16,i9,a3)','Repartition ',nom,' ( population : ',size(etannz),' ) '
           print*, 

           if (dPlay == 2) then
               !call  hist0(rgt,dist,color='#779921',pause=2.,persist='no')
               call  hist0(rgt,cumul,color='#779921',pause=4.,persist='no', moy=moy, ect=ect)
           end if
       end if

    end if

    deallocate(etannz)

  end subroutine alea_trac

!------------------------------------------------------------!

subroutine distrib_volumes(Ve,Dg,interv0,dPlay)
!------------------------------------------------------------!
!     Trace graphique de la distribution des volumes des     !
!       elements pour analyse du probleme traite         !
!------------------------------------------------------------!
    use variables, only : nelt,pi
    use initialisation, only : init_vec
    use gnufor2

    implicit none

    integer, intent(in) :: dPlay, interv0
    real*8, intent(in) :: Dg
    real*8 :: moy, ect, tmp, xpeti,  mx, mn, dm, rap, dia
    integer :: ie, i, nnz, j, nEch, nEch_prec, interv
    real*8, dimension(:), intent(in) :: Ve
    real*8, dimension(:), allocatable :: Vennz
    real*8, dimension(interv0+1) :: rgt, dist

    nnz = count(Ve /= 0.)

    call init_vec(Vennz,nnz)
    i = 1

    do ie = 1, nelt
       if (Ve(ie) /= 0.) then
           Vennz(i) = Ve(ie)
       i = i + 1
       end if
    end do   


    if (dPlay /= 0) then
        moy = (sum(Vennz)/size(Vennz))
          tmp = 0.
          do i = 1,size(Vennz)
         tmp = tmp + (Vennz(i) - moy)**2 
          end do

    if (size(Vennz)==1) then
        ect = sqrt(tmp)
    else
        ect = sqrt(tmp/(size(Vennz)-1))
    end if

    xpeti = 1.d-10
    interv = interv0
    mx = maxval(Vennz)
        mn = minval(Vennz)
        if (abs((mx-mn)/(mn+xpeti)) < xpeti) then
       mn = 0.
       interv = 1
    end if
        dm = (mx - mn)/interv

           do i = 0, interv
           rgt(i+1) = mn + dm*i
           end do

           nEch_prec = 0
           nEch = 0
        dist = 0

           do i = 1,(interv+1)
      nEch = count((Vennz <= rgt(i)) .eqv. .true.)
          dist(i) = nEch - nEch_prec
          nEch_prec = nEch
           end do

    rap = (sum(Vennz)/size(Vennz))/((4./3.)*pi*(Dg/2.)**3)

        print*,'===> Volumes elementaires :'
        print'(a18,e12.5,1x,a19,e12.5)','Valeur moyenne : ', moy ,' -   Ecart type : ', ect
        print*, 'Diametre du grain (m)     :  ', Dg
        print*, 'Rapport moyen R = Ve / Vg :  ', rap
        print*,

        if (dPlay == 2) then
            call hist0(rgt,dist,color='#779921',pause=4.,persist='no', moy=moy, ect=ect,dia=Dg, rap=rap)
        end if
    end if
    
    deallocate(Vennz) 
   
  end subroutine distrib_volumes

  !-------------------------------------------------------!

  subroutine rossi_coef(m,s,r,fc,mot,youn)
  !-------------------------------------------------------
  ! Calcul des parametres m et s de la loi normale 
  ! en fonction de r = V/Vg et de fc
  ! valable pour :
  ! * la resistance : mot = 'RESI'
  ! * et le module  : mot = 'MODU'
  !-------------------------------------------------------

  implicit none

  real*8, intent(out) :: m, s
  real*8, intent(in)  :: r, fc
  real*8, intent(in), optional :: youn
  character*4 :: mot
  
  real*8 :: A, B, AA, BB, C, D

  !----- Verification des donnees d'entree
  !if ((r<1.d0).or.(r>100000.d0)) stop 'prop_alea.f - rossi_coef :: donnees hors domaine de validite du rapport V/Vg'
  if ((fc<25.d0).or.(fc>130.d0)) stop 'prop_alea.f - rossi_coef :: donnees hors domaine de validite de la resistance fc'
  
  if (mot=='RESI') then
      !----- Definition des parametres (RESISTANCE)
      A = 6.5d0
      B = 0.25d0 - 3.6d-03 * fc + 1.3d-05 * (fc**2.)
      AA = 0.35d0
      BB = 4.5d-02 + 4.5d-03 * fc - 1.8d-5 * (fc**2.)

      !----- valeur moyenne et ecart type
      m = A*exp(-B*log(r))
      s = AA*exp(-BB*log(r))
      s = s*m
  end if
  if (mot=='MODU') then
      if(present(youn)) then
          m = youn
      else
          m = 11000.d0 * (fc **(1.d0/3.d0))
      end if
      !----- Definition des parametres (MODULE)
      AA = .15d0
      BB = 0.16d0 + 2.7d-03 * fc- 3.4d-06 * (fc**2.)
      
      !----- valeur moyenne et ecart type
      s = AA*exp(-BB*log(r))
      s = s*m
  end if
 
  end subroutine rossi_coef
  !-------------------------------------------------------!

  subroutine wbl_coef(b,c,r,fc)
  !-------------------------------------------------------
  ! Calcul des parametres b et c de la loi de Weibull 
  ! en fonction de r = V/Vg et de fc
  !-------------------------------------------------------

  implicit none

  real*8, intent(out) :: b, c
  real*8, intent(in) :: r, fc  
  
  integer :: id
  real*8 :: X, AA, BB, CC, log_r, a, a1, a2, b1, b2, bref, &
             & C_fc, coef, cref, d, D_fc, e, f, g

  !-------------------------------------------------------

  id=0

  if (id==0) then
         X=log(r)
    
         AA = 2.8711d+04;  BB = -1.2265d-05; CC = -2.8710d+04
        !AA = 2.5198d+04; BB = -1.3430d-05; CC = -2.5197d+04 ! B 50
        !AA = 2.7200d+04; BB = -1.3241d-05; CC = -2.7199d+04 ! B 25
        !AA = 2.3809d+04; BB = -1.4557d-05; CC = -2.3808d+04 ! B 25 1er fittage
        c = exp(AA*exp(-BB*X)+CC)
    
         AA = 2.5687d+00;  BB = 1.155d-01;  CC = -9.4046d-02 
        !AA = 2.2272d+00; BB = 1.2007d-01; CC = 1.7918d-01  ! B 50
        !AA = 2.8602d+00; BB = 1.1022d-01; CC = -7.1673d-01 ! B 25
        !AA = 3.6736d+00; BB = 8.2271d-02; CC = -1.3251d+00 ! B 25 1er fittage
        b = exp(AA*exp(-BB*X)+CC)

  else
    ! CI DESSOUS : EN ATTENTE DE VALIDATION

    if (r>10000.) stop 'FIDES_wblcoef : attention ! Rapport de volumes hors du domaine de validite (>10000) ! Changer le maillage !'
    if (r<0.01) stop 'FIDES_wblcoef : attention ! Rapport de volumes hors du domaine de validite (<0.01) ! Changer le maillage !'

    log_r = log(r)/log(10.)

    ! Parametres pour l'expression de b(r,fc)
    !a = 7.909
    a = 5.
    b = 1.074
    c = 3.136
    !d = 0.4
    d = 0.01
    a1 = 0.0084
    a2 = 0.4942

    C_fc = 1.-(a1*fc+a2) 
    coef = 1.-C_fc*(1-exp(-d*r))
    bref = a*exp(-b*log_r)+c
    b = bref*coef

    ! Parametres pour l'expression de c(r,fc)
    e = 46.39
    f = 6.26
    g = 3.625
    b1 = 0.0076
    b2 = 0.5398

    D_fc = b1*fc+b2
    !cref = e*exp(-((log_r-f)/g)**2)
    cref = 45.*bref**(-1.2)

    coef = D_fc
    c = cref*coef
  end if
    
  end subroutine wbl_coef

  !-------------------------------------------------------!

  function distr_alea(loi,nEchant,ie0)
  !-------------------------------------------------------
  !  TU_alea(loi,nEchant,ie0);
  !
  !  Fonction de distribution aleatoire (Gaussienne)
  !
  !  Parametres :
  !     loi     : tableau contenant des informations sur la loi de distribution
  !     nEchant : nombre d'echantillons generes
  !
  !  Remarque :
  !
  !-------------------------------------------------------
        use math, only : randgp

        implicit none
        real*8, dimension(5), intent(in) :: loi

        integer, intent(in) :: nEchant,ie0
        real*8, dimension(nEchant) :: distr_alea
        real*8 :: et, s, m, s0, m0, v, b, c
        integer :: nEch, interv, i
        logical :: itest

        real*8 :: somme
        !real*8, dimension(12+ie0) :: alea
        real*8, dimension(12) :: alea

        !-------------------------------------------------------
        ! Traitement des donnees d'entree

        itest = .false.

        if (loi(5) == 1) then
            ! Tirage aleatoire de variables a distrib gaussienne
            m = loi(3) ! Valeur moyenne
            s = loi(4) ! Ecart type

            nEch = 0
            itest = .true.

            do while (itest)
                somme = 0.

                do i =1,12
                    alea = randgp(12,(/ie0/))
                    somme = somme + (alea(i)-0.5)
                end do

                et = somme*s+m

                !if (et <= 0) et = abs(et)
                if (et>=0) then

                    nEch = nEch+1
                    distr_alea(nEch) = et
                    if (nEch == nEchant)  itest = .false.
                    
                end if

            end do

        elseif (loi(5) == 2) then
            ! Tirage aleatoire de variables a distrib de Weibull
            b = loi(1)
            c = loi(2)
            nEch = 0
            itest = .true.

            do while (itest)
                alea = randgp(1,(/ie0/))
                et = b*(-log(1-alea(1)))**(1./c)

                if (et >= 0) then
                    nEch = nEch+1
                    distr_alea(nEch) = et
                    if (nEch == nEchant)  itest = .false.
                end if
            end do

        elseif (loi(5) == 3) then
            ! Tirage aleatoire de variables a distrib log-Normale
            m0 = loi(3)
            s0 = loi(4)

            v = log(1. + (s0**2.)/(m0**2.))
            s = sqrt(v)
            m = log(m0)-1./2.*v

            nEch = 0
            itest = .true.

            do while (itest)
                somme = 0.

                do i =1,12            
                    alea = randgp(12,(/ie0/))
                    somme = somme + (alea(i)-0.5)
                end do

                et = somme*s+m
                nEch = nEch + 1
                distr_alea(nEch) = et
                if (nEch == nEchant)  itest = .false.
            end do

            distr_alea = exp(distr_alea)
        end if

  end function distr_alea

!-------------------------------------------------------!

end module proprietes_aleatoires


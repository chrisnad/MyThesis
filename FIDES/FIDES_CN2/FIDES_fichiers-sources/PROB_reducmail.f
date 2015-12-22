subroutine PROB_reducmail()

! Script de modification du maillage

    use variables, only : dime, vcor, vprelg, aval, pi, idmax, Dg, &
            & redu,rap
    real*8  :: V, Vg,VsVg

        if (redu) then

        ! Modification des coordonnées
        vcor=rap*vcor

        ! Modification de l'epaisseur
            if (dime==2) then
            vprelg(:,idmax-2)=vprelg(:,idmax-2)*rap
            end if

            ! fin modification de l'epaisseur
            aval = aval * rap


        V=rap**3*1.d0
        Vg=pi*(Dg**3.)/6.
        VsVg=V/Vg

        print*
        print*,'Reduction automatique de maillage :'
        print'(a30,E12.5)','  - rapport de reduction    :',rap
        print'(a30,E12.5)','  - rapport de volumes V/Vg :',VsVg
        print*

    else
        rap = 0.
    end if

end subroutine PROB_reducmail


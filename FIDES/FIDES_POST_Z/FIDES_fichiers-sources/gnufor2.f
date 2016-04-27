! Author: Alexey Kuznetsov
! Modified: 28/12/2008
! this Fortran90 module contains a collection of subroutines for plotting data,
! including 2D, 3D plots, surfaces, polar coordinates, histograms
! it is a modification of the GNUFOR interface written by John Burkardt:
! http://orion.math.iastate.edu/burkardt/g_src/gnufor/gnufor.html
!***********************************************************************************
    module gnufor2
    implicit none
!***********************************************************************************
! these are default parameters which control linewidth, colors and terminal
!***********************************************************************************
    character(len=3), parameter :: default_linewidth='1'
    character(len=100), parameter   :: default_color1='blue'
    character(len=100), parameter   :: default_color2='dark-green'
    character(len=100), parameter   :: default_color3='orange-red'
    character(len=100), parameter   :: default_color4='dark-salmon'
    character(len=100), parameter   :: default_terminal='wxt'
    character(len=100), parameter   :: default_palette='CMY'

contains

!***********************************************************************************

function my_date_and_time() result(f_result)

!***********************************************************************************
! this function creates a string with current date and time
! it is a default method to name output files
!***********************************************************************************

    implicit none
        character(len=8)  :: date
        character(len=10) :: time
        character(len=33) :: f_result

    call date_and_time(date,time)
    f_result= 'date_'//date(7:8)//'-'//date(5:6)//'-'//date(1:4)//'_time_'//time(1:2)//':'//time(3:4)//':'//time(5:10)

end function my_date_and_time

!***********************************************************************************

function output_terminal(terminal) result(f_result)

    implicit none
    character(len=*),intent(in) :: terminal
    integer, parameter      :: Nc=35
    character(len=Nc)       :: f_result

    select case(terminal)
        case('ps')
            f_result='postscript landscape color'
        case default
            f_result=terminal
    end select

end function output_terminal

!***********************************************************************************

subroutine hist0(x,y,pause,color,persist,moy,ect,dia,rap)

!***********************************************************************************
! this subroutine plots the histogram of data contained in array x, using n bins
!***********************************************************************************

    implicit none
    real(kind=8), intent(in)    :: x(:),y(:) !the data to plot
    real(kind=8), intent(in), optional  :: dia, rap , moy, ect
    real(kind=4), optional      :: pause
    character(len=*),optional   :: color, persist
    integer             :: i, ierror, ios, file_unit, n
    character(len=100)      :: data_file_name, command_file_name, my_color, &
                & my_pause, my_persist, xrange1, xrange2, yrange, moyen, ecart, diam, rapp

!***********************************************************************************
    write (xrange1,'(e15.7)') minval(x)*0.995
    write (xrange2,'(e15.7)') maxval(x)*1.005
    write (yrange,'(e15.7)') maxval(y)*1.1

    write (moyen,'(e15.7)') moy
    write (ecart,'(e15.7)') ect
    write (diam,'(e15.7)') dia
    write (rapp,'(e15.7)') rap

!***********************************************************************************
    data_file_name='data_file.txt'
    command_file_name='command_file.txt'

!***********************************************************************************
    ierror=0
    call get_unit(file_unit)
    if (file_unit==0) then
        ierror=1
        print *,'write_vector_date - fatal error! Could not get a free FORTRAN unit.'
        stop
    endif
    open (unit=file_unit, file=data_file_name, status='replace', iostat=ios)
    if (ios/=0) then
        ierror=2
        print *,'write_vector_data - fatal error! Could not open the terminal data file.'
        stop
    endif
!***********************************************************************************
! here we write the date to the data_file - the gnuplot will read this data later
!***********************************************************************************
    n = size(x)
    do i = 1,n
        write (file_unit,'(2E15.7)') x(i), y(i)
    enddo
!***********************************************************************************
    close (unit=file_unit)
!***********************************************************************************
    ierror = 0
    call get_unit(file_unit)
    if (file_unit==0) then
        ierror=1
        print *,'write_vector_date - fatal error! Could not get a free FORTRAN unit.'
        stop
    endif
    open (unit=file_unit, file=command_file_name, status='replace', iostat=ios)
    if (ios/=0) then
        ierror=2
        print *,'write_vector_data - fatal error! Could not open the terminal command file.'
        stop
    endif
!***********************************************************************************
! here we write the commands to the commands file which gnuplot will execute
!***********************************************************************************
    my_persist='persist'
    if (present(persist).and.(persist=='no')) my_persist=' '

    if (present(moy)) then
       write ( file_unit, '(a)' ) 'set title " Moyenne : ' // trim(moyen) // '      Ecart-type ' &
            & //trim(ecart) // '"'
    endif
    if (present(dia)) then
           write ( file_unit, '(a)' ) 'set key top ti "Diametre :  '// trim(diam) // &
            & '   Rapport Ve/Vg :  ' //trim(rapp) // '"'
    else
       write ( file_unit, '(a)' ) 'set nokey'
    endif
!***********************************************************************************

    write ( file_unit, '(a)' ) 'set yrange [0.0:'// trim(yrange) //']'
    write ( file_unit, '(a)' ) 'set xrange ['// trim(xrange1) // ':'// trim(xrange2) //']'

    write ( file_unit, '(a)' ) 'set style data histograms'
    write ( file_unit, '(a)' ) 'set style fill solid border -1'
!***********************************************************************************
    if (present(color)) then
        my_color='"'//color//'"'
    else
        my_color='"'//trim(default_color1)//'"'
    endif
    write ( file_unit, '(a)' ) 'plot "' // trim ( data_file_name ) // &
    '" using 1:2 notitle with boxes linecolor rgb' // trim(my_color)
!***********************************************************************************
    if (present(pause)) then
        if (pause<0.0) then
            write ( file_unit, '(a)' ) 'pause -1 "press RETURN to continue"'
        else
            write ( my_pause,'(e9.3)') pause
            write ( file_unit, '(a)' ) 'pause ' // trim(my_pause)
        endif
    else
        write ( file_unit, '(a)' ) 'pause 0'
    endif
!***********************************************************************************
    write ( file_unit, '(a)' ) 'q'
    close ( unit = file_unit )
!***********************************************************************************
    call run_gnuplot (command_file_name)

end subroutine hist0

!***********************************************************************************

subroutine hist(x,n,pause,color,terminal,filename,persist,input)

!***********************************************************************************
! this subroutine plots the histogram of data contained in array x, using n bins
!***********************************************************************************

    implicit none
    real(kind=8), intent(in)    :: x(:) !the data to plot
    integer, intent(in)     :: n !the number of intervals
    real(kind=4), optional      :: pause
    character(len=*),optional   :: color, terminal, filename, persist, input
    integer             :: i, ierror, ios, file_unit, nx
    character(len=100)      :: data_file_name, command_file_name, yrange, xrange1, xrange2, my_color, &
                & xtic_start, dxtic, xtic_end, my_pause, my_persist, nbelem
    real(kind=8)            :: xmin, xmax, dx

!***********************************************************************************
! prepare the data
    nx=size(x)
    xmin=minval(x)
    xmax=maxval(x)
    dx=(xmax-xmin)/n

!***********************************************************************************
    write (nbelem,'(i5)') n     ;   write (dxtic,'(i5)') n-1
    write (xrange1,'(i5)') 0    ;   write (xrange2,'(i5)') n+1
    write (xtic_start,'(i5)') 1 ;   write (xtic_end,'(i5)') n

    write (yrange,'(e15.7)') maxval(x)*1.05

!***********************************************************************************
    if (present(input)) then
        data_file_name='data_file_'//input//'.txt'
        command_file_name='command_file_'//input//'.txt'
    else
        data_file_name='data_file.txt'
        command_file_name='command_file.txt'
    endif
!***********************************************************************************
    ierror=0
    call get_unit(file_unit)
    if (file_unit==0) then
        ierror=1
        print *,'write_vector_date - fatal error! Could not get a free FORTRAN unit.'
        stop
    endif
    open (unit=file_unit, file=data_file_name, status='replace', iostat=ios)
    if (ios/=0) then
        ierror=2
        print *,'write_vector_data - fatal error! Could not open the terminal data file.'
        stop
    endif
!***********************************************************************************
! here we write the date to the data_file - the gnuplot will read this data later
!***********************************************************************************
    do i=1,n
        write (file_unit,'(i5,E15.7)') i, x(i)
    enddo
!***********************************************************************************
    close (unit=file_unit)
!***********************************************************************************
    ierror = 0
    call get_unit(file_unit)
    if (file_unit==0) then
        ierror=1
        print *,'write_vector_date - fatal error! Could not get a free FORTRAN unit.'
        stop
    endif
    open (unit=file_unit, file=command_file_name, status='replace', iostat=ios)
    if (ios/=0) then
        ierror=2
        print *,'write_vector_data - fatal error! Could not open the terminal command file.'
        stop
    endif
!***********************************************************************************
! here we write the commands to the commands file which gnuplot will execute
!***********************************************************************************
    my_persist='persist'
    if (present(persist).and.(persist=='no')) my_persist=' '
    if (present(terminal)) then
        write ( file_unit, '(a)' ) 'set terminal '// trim(output_terminal(terminal))
    if (present(filename)) then
        write ( file_unit, '(a)' ) 'set output "'// trim(filename) //'"'
    else
        write ( file_unit, '(a)' ) 'set output "'//my_date_and_time()//'"'
    endif
    else
        write ( file_unit, '(a)' ) 'set terminal ' // trim(default_terminal) // ' ' &
            & //trim(my_persist) // ' title  "Gnuplot"'
    endif
!***********************************************************************************
    write ( file_unit, '(a)' ) 'set nokey'
    write ( file_unit, '(a)' ) 'set xlabel "' //  trim(nbelem)  // ' elements "'
    write ( file_unit, '(a)' ) 'set yrange [0.0:'// trim(yrange) //']'
    write ( file_unit, '(a)' ) 'set xrange ['// trim(xrange1) // ':'// trim(xrange2) //']'
    write ( file_unit, '(a)' ) 'set noxtics'
    write ( file_unit, '(a)' ) 'set style data histograms'
    write ( file_unit, '(a)' ) 'set style fill solid border -1'
!***********************************************************************************
    if (present(color)) then
        my_color='"'//color//'"'
    else
        my_color='"'//trim(default_color1)//'"'
    endif
    write ( file_unit, '(a)' ) 'plot "' // trim ( data_file_name ) // &
    '" using 1:2 with boxes linecolor rgb' // trim(my_color)
!***********************************************************************************
    if (present(pause)) then
        if (pause<0.0) then
            write ( file_unit, '(a)' ) 'pause -1 "press RETURN to continue"'
        else
            write ( my_pause,'(e9.3)') pause
            write ( file_unit, '(a)' ) 'pause ' // trim(my_pause)
        endif
    else
        write ( file_unit, '(a)' ) 'pause 0'
    endif
!***********************************************************************************
    write ( file_unit, '(a)' ) 'q'
    close ( unit = file_unit )
!***********************************************************************************
    call run_gnuplot (command_file_name)

end subroutine hist


!***********************************************************************************


subroutine plot_1(x1,y1,style,pause,color1,terminal,filename,polar,persist,input,linewidth)
!***********************************************************************************
! this subroutine plots a two-dimensional graph
!***********************************************************************************
    implicit none
    real(kind=8), intent(in)    :: x1(:), y1(:)
    real(kind=4), optional      :: pause,linewidth
    character(len=*),optional   :: style, color1, terminal, filename, polar, persist, input
    integer             :: i, ierror, ios, file_unit, Nx1
    character(len=100)      :: data_file_name, command_file_name, my_linewidth
    integer, parameter      :: Nc=20
    character(len=Nc)       :: my_line_type1, my_color1, my_range, my_pause, my_persist
!***********************************************************************************
    if (present(input)) then
        data_file_name='data_file_'//input//'.txt'
        command_file_name='command_file_'//input//'.txt'
    else
        data_file_name='data_file.txt'
        command_file_name='command_file.txt'
    endif
!***********************************************************************************
    Nx1=size(x1)
    if ((size(x1).ne.size(y1))) then
        print *,'subroutine plot ERROR: size(x) is not equal to size(y)'
        stop
    endif
    if (present(style).and.(len(style).ne.3)) then
        print *,'subroutine plot ERROR: argument "style" has wrong number of characters'
        stop
    endif
!***********************************************************************************
    ierror=0
    call get_unit(file_unit)
    if (file_unit==0) then
        ierror=1
        print *,'write_vector_data - fatal error! Could not get a free FORTRAN unit.'
        stop
    endif
    open (unit=file_unit, file=data_file_name, status='replace', iostat=ios)
    if (ios/=0) then
        ierror=2
        print *,'write_vector_data - fatal error! Could not open the terminal data file.'
        stop
    endif
!***********************************************************************************
! here we write the date to the data_file - the gnuplot will read this data later
!***********************************************************************************
    do i=1,Nx1
        write (file_unit,'(2E15.7)') x1(i), y1(i)
    enddo
!***********************************************************************************
    close (unit=file_unit)
!***********************************************************************************
    ierror = 0
    call get_unit(file_unit)
    if (file_unit==0) then
        ierror=1
        print *,'write_vector_data - fatal error! Could not get a free FORTRAN unit.'
        stop
    endif
    open (unit=file_unit, file=command_file_name, status='replace', iostat=ios)
    if (ios/=0) then
        ierror=2
        print *,'write_vector_data - fatal error! Could not open the terminal command file.'
        stop
    endif
!***********************************************************************************
! here we write the commands to the commands file which gnuplot will execute
!***********************************************************************************
    my_line_type1='lines'
    if (present(style)) then
    if ((style(3:3)=='-')) then
        my_line_type1='linespoints'
    else
        my_line_type1='points'
    endif
    endif
    if (present(linewidth)) then
        write ( my_linewidth,'(e9.3)') linewidth
    else
        my_linewidth=trim(default_linewidth)
    endif
    if (present(color1)) then
        my_color1='"'//trim(color1)//'"'
    else
        my_color1='"'//trim(default_color1)//'"'
    endif
!***********************************************************************************
    my_persist='persist '
    if (present(persist).and.(persist=='no')) my_persist=' '
    if (present(terminal)) then
        write ( file_unit, '(a)' ) 'set terminal '// trim(output_terminal(terminal))
    if (present(filename)) then
        write ( file_unit, '(a)' ) 'set output "'// trim(filename) //'"'
    else
        write ( file_unit, '(a)' ) 'set output "'//my_date_and_time()//'"'
    endif
    else
        write ( file_unit, '(a)' ) 'set terminal ' // trim(default_terminal) // ' ' &
            & //trim(my_persist) //' title  "Gnuplot"'
    endif
!***********************************************************************************
    write ( file_unit, '(a)' ) 'unset key'
    if (present(polar).and.(polar=='yes')) then
        write (my_range,'(e15.7)') maxval(abs(y1))
        write ( file_unit, '(a)' ) 'set xrange [-'//trim(my_range)//':'//trim(my_range)//']'
        write ( file_unit, '(a)' ) 'set yrange [-'//trim(my_range)//':'//trim(my_range)//']'
        write ( file_unit, '(a)' ) 'set size square'
        write ( file_unit, '(a)' ) 'set polar'
        write ( file_unit, '(a)' ) 'set grid polar'
    else
        write ( file_unit, '(a)' ) 'set grid'
    endif
!***********************************************************************************
    if (present(style)) then
        write ( file_unit, '(a,i2,a)' ) 'plot "' // trim (data_file_name) &
        &//'" using 1:2 with ' // trim(my_line_type1) // ' pointtype ' // &
        & style(1:2) // ' linecolor rgb ' // trim(my_color1) // ' linewidth '// trim(my_linewidth)
    else
        write ( file_unit, '(a,i2,a)' ) 'plot "' // trim (data_file_name) &
        & //'" using 1:2 with ' // trim(my_line_type1)  // ' linecolor rgb '&
        & // trim(my_color1) // ' linewidth '// trim(my_linewidth)
    endif
!***********************************************************************************
    if (present(pause)) then
        if (pause<0.0) then
            write ( file_unit, '(a)' ) 'pause -1 "press RETURN to continue"'
        else
            write ( my_pause,'(e9.3)') pause
            write ( file_unit, '(a)' ) 'pause ' // trim(my_pause)
        endif
    else
        write ( file_unit, '(a)' ) 'pause 0'
    endif
!***********************************************************************************
    write ( file_unit, '(a)' ) 'q'
    close ( unit = file_unit )
!***********************************************************************************
    call run_gnuplot (command_file_name)
!***********************************************************************************
end subroutine plot_1
!***********************************************************************************



!***********************************************************************************

subroutine run_gnuplot(command_file_name)

    implicit none
    character (len = 100) command
    character (len = *) command_file_name
    integer status
    integer system
!***********************************************************************************
!  Issue a command to the system that will startup GNUPLOT, using
!  the file we just wrote as input.
!***********************************************************************************
    write (command, *) 'gnuplot ' // trim (command_file_name)
    status=system(trim(command))
    if (status.ne.0) then
        print *,'RUN_GNUPLOT - Fatal error!'
        stop
    endif
    return

end subroutine run_gnuplot

!***********************************************************************************

subroutine get_unit(iunit)

    implicit none
    integer i
    integer ios
    integer iunit
    logical lopen

    iunit=0
    do i=1,99
        if (i/= 5 .and. i/=6) then
            inquire (unit=i, opened=lopen, iostat=ios)
            if (ios==0) then
                if (.not.lopen) then
                    iunit=i
                    return
                endif
            endif

        endif
    enddo
    return
end subroutine get_unit
!***********************************************************************************
end module gnufor2

program main
    use mod_rodas3_standAlone, only: &
         read_input &
        ,chem_rodas3_dyndt
    implicit none

    integer :: processor
    real :: time
    integer :: iUnit, err
    logical :: ex

    namelist /TARGET/ time,processor

    inquire(file='target.nml', exist=ex)
    if (.not. ex) then
        print *, "Error: file target.nml don't exist!"
        print *, "Please, verify!"
        stop
    end if
    open (newunit=iUnit, FILE='target.nml', STATUS='OLD')
    read (iunit, iostat=err, NML=TARGET)
    if (err /= 0) then
        print *, "Error: reading target.nml!"
        print *, "Please, check syntax!"
        stop   
    else
        print *,"Time      = ",time
        print *,"Processor = ",processor
    end if
    close(iUnit)

    call read_input(time,processor)
    call chem_rodas3_dyndt(time, processor)

end program main
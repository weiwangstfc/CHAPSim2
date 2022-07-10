!> statistics.f90
!> \brief A file containing modules  to enable on the fly
!> statistics generation for chapsim 2. 
!> The modules are:
!>  * stat_types: module containing low-level routines
!>                for processing statistics
!>  * chapsim_stats: high-level module containing 
!>                   routines which the user should use to
!>                   generate statistics
!>  * statistics: Interface for use with the rest of chapsim 2
!_____________________________________________________________

module stat_types
    use geometry_mod
    use domain_decomposition_mod
    use mpi_mod
    use udf_type_mod
    use vars_df_mod
    use decomp_2d
    use parameters_constant_mod
    use iso_fortran_env
    use MPI
    implicit none

    public 

    integer, parameter, public :: ISTAT_ACCUM_ASSYM = 1
    integer, parameter, public :: ISTAT_ACCUM_LOCAL = 2

    integer, parameter, public :: ISTAT_AVG_XZ = 3
    integer, parameter, public :: ISTAT_AVG_Z = 4
    integer, parameter, public :: ISTAT_AVG_N = 5

    type t_stats_base
        character(len=10)                    :: name
        integer                              :: idomain
        integer                              :: accumulation
        integer                              :: xz_avg
        real(wp)                             :: dt_calc
        real(wp)                             :: dt_write
        real(wp)                             :: start_time
        integer                              :: size
        integer                              :: nstatis

    end type t_stats_base

    type, extends(t_stats_base) :: t_stats_xz
        procedure(stat_calc_xz), pointer, nopass   :: stat_calc
        real(wp), dimension(:,:), allocatable :: data_array 
    end type t_stats_xz

    type, extends(t_stats_base) :: t_stats_z
        procedure(stat_calc_z), pointer, nopass   :: stat_calc
        real(wp), dimension(:,:,:), allocatable :: data_array         
    end type t_stats_z

    type, extends(t_stats_base) :: t_stats_N
        procedure(stat_calc_N), pointer, nopass   :: stat_calc
        real(wp), dimension(:,:,:,:), allocatable :: data_array         
    end type t_stats_N

    abstract interface
        subroutine stat_calc_xz(fl,dm,size, array)
            import :: t_flow, t_domain, wp
            type(t_flow) :: fl
            type(t_domain) :: dm
            integer :: size
            real(wp), dimension(size, &
                                dm%dccc%ysz(2)) :: array
        end subroutine
    end interface

    abstract interface
        subroutine stat_calc_z(fl,dm,size, array)
            import :: t_flow, t_domain, wp
            type(t_flow) :: fl
            type(t_domain) :: dm
            integer :: size
            real(wp), dimension(size, &
                                dm%dccc%zsz(1), &
                                dm%dccc%zsz(2)) :: array
        end subroutine
    end interface

    abstract interface
        subroutine stat_calc_N(fl,dm,size, array)
            import :: t_flow, t_domain, wp
            type(t_flow) :: fl
            type(t_domain) :: dm
            integer :: size
            real(wp), dimension(size, &
                                dm%dccc%xsz(1), &
                                dm%dccc%xsz(2), &
                                dm%dccc%xsz(3)) :: array
        end subroutine
    end interface

    interface allocate_data_array
        module procedure allocate_data_array_N, allocate_data_array_xz, allocate_data_array_z
    end interface

    interface average_stats
        module procedure average_stats_N, average_stats_xz, average_stats_z
    end interface


    interface write_stat
        module procedure write_stat_N, write_stat_xz, write_stat_z
    end interface
    
    contains
        subroutine set_stat_info(stat, stat_name, size, idm, &
            accumulation_type, xz_avg, &
            dt_calc, dt_write, start_time)

            class(t_stats_base), intent(INOUT) :: stat
            character(len=*),                   intent(in) :: stat_name
            integer,                            intent(in) :: size
            integer,                            intent(in) :: idm
            integer,                            intent(in) :: accumulation_type
            integer,                            intent(in) :: xz_avg
            real(wp),                            intent(in) :: dt_calc
            real(wp),                            intent(in) :: dt_write
            real(wp),                           intent(in) :: start_time

            stat%name = stat_name
            stat%idomain = idm
            stat%xz_avg = xz_avg
            stat%accumulation = accumulation_type
            stat%dt_calc = dt_calc
            stat%dt_write = dt_write
            stat%start_time = start_time
            stat%size = size
            stat%nstatis = 0

        end subroutine

        subroutine allocate_data_array_xz(stat)
            type(t_stats_xz) :: stat
            integer :: s1, s2

            s1 = stat%size
            s2 = domain(stat%idomain)%nc(2)
            
            allocate(stat%data_array(s1,s2))

        end subroutine allocate_data_array_xz
        subroutine allocate_data_array_z(stat)
            type(t_stats_z) :: stat
            integer :: s1, s2, s3

            s1 = stat%size
            s2 = domain(stat%idomain)%dccc%zsz(1)
            s3 = domain(stat%idomain)%dccc%zsz(2)

            allocate(stat%data_array(s1,s2,s3))

        end subroutine allocate_data_array_z
        subroutine allocate_data_array_N(stat)
            type(t_stats_N) :: stat
            integer :: s1, s2, s3, s4

            s1 = stat%size
            s2 = domain(stat%idomain)%dccc%xsz(1)
            s3 = domain(stat%idomain)%dccc%xsz(2)
            s4 = domain(stat%idomain)%dccc%xsz(3)

            allocate(stat%data_array(s1,s2,s3,s4))

        end subroutine allocate_data_array_N


        subroutine average_stats_xz(stat)
            type(t_stats_xz) :: stat

            real(wp), allocatable :: avg(:,:)
            
            integer :: s2_y

            if (.not. check_time_calc(stat)) return

            s2_y = domain(stat%idomain)%dccc%ysz(2)

            allocate(avg(stat%size, s2_y))

            call stat%stat_calc(flow(stat%idomain),&
                                domain(stat%idomain),&
                                stat%size, avg)
        

            call accumulate_xz(stat, avg)                                 

            deallocate(avg)

        end subroutine average_stats_xz

        logical function check_time_calc(stat) 
            class(t_stats_base) :: stat

            real(wp) :: dt, time

            dt = domain(stat%idomain)%dt
            time = flow(stat%idomain)%time

            if (dmod(time, stat%dt_calc) .lt. dt .and.&
                         time.gt.stat%start_time) then
                check_time_calc = .true.
            else
                check_time_calc = .false.
            endif

        end function

        logical function check_time_write(stat) 
            class(t_stats_base) :: stat

            real(wp) :: dt, time

            dt = domain(stat%idomain)%dt
            time = flow(stat%idomain)%time

            if (dmod(time, stat%dt_write) .lt. dt) then
                check_time_write = .true.
            else
                check_time_write = .false.
            endif
        end function

        subroutine accumulate_xz(stat, avg)
            type(t_stats_xz) :: stat
            real(wp) :: avg(stat%size, domain(stat%idomain)%dccc%ysz(2))
            real(wp) :: coe1, coe2
            stat%nstatis = stat%nstatis + 1
            select case (stat%accumulation)
                case(ISTAT_ACCUM_ASSYM)
                    coe1 = dble(stat%nstatis - 1)/ dble(stat%nstatis)
                    coe2 = 1.0_wp/dble(stat%nstatis)

                    stat%data_array(:,:) = coe1*stat%data_array(:,:) &
                                                + coe2*avg

                case default
                    stat%data_array(:,:) = avg
            end select

        end subroutine

        subroutine average_stats_z(stat)
            type(t_stats_z) :: stat

            real(wp), allocatable :: avg_local(:,:,:)
            
            integer :: s1_z, s2_z

            if (.not. check_time_calc(stat)) return

            s1_z = domain(stat%idomain)%dccc%zsz(1)
            s2_z = domain(stat%idomain)%dccc%zsz(2)
            
            allocate(avg_local(stat%size,s1_z,s2_z))

            avg_local(:,:,:) = ZERO

            call stat%stat_calc(flow(stat%idomain),&
                                domain(stat%idomain),&
                                stat%size, avg_local)


            call accumulate_z(stat,avg_local)

            deallocate(avg_local)

        end subroutine
        
        subroutine accumulate_z(stat, avg)
            type(t_stats_z), intent(inout) :: stat
            real(wp), intent(in) :: avg(stat%size,&
                             domain(stat%idomain)%dccc%zsz(1),&
                             domain(stat%idomain)%dccc%zsz(2))
            real(wp) :: coe1, coe2

            stat%nstatis = stat%nstatis + 1
            select case (stat%accumulation)
                case(ISTAT_ACCUM_ASSYM)

                    coe1 = dble(stat%nstatis - 1)/ dble(stat%nstatis)
                    coe2 = 1.0_wp/dble(stat%nstatis)

                    stat%data_array(:,:,:) = coe1*stat%data_array(:,:,:) &
                                                + coe2*avg(:,:,:)

                case (ISTAT_ACCUM_LOCAL)
                    stat%data_array(:,:,:) = avg
                    
                case default
                    call Print_error_msg("Somethink has gone wrong")
                
            end select

        end subroutine

        subroutine average_stats_N(stat)
            type(t_stats_N) :: stat

            real(wp), allocatable :: array(:,:,:,:)
            integer :: s1_x, s2_x, s3_x
            
            if (.not. check_time_calc(stat)) return

            s1_x = domain(stat%idomain)%dccc%xsz(1)
            s2_x = domain(stat%idomain)%dccc%xsz(2)
            s3_x = domain(stat%idomain)%dccc%xsz(3)
            
            allocate(array(stat%size, s1_x, s2_x, s3_x))
            call stat%stat_calc(flow(stat%idomain),&
                                domain(stat%idomain),&
                                stat%size, array)
            
            call accumulate_N(stat,array)

            deallocate(array)

        end subroutine

        subroutine accumulate_N(stat, avg)
            type(t_stats_N) :: stat
            real(wp) :: avg(stat%size,&
                             domain(stat%idomain)%dccc%xsz(1),&
                             domain(stat%idomain)%dccc%xsz(2),&
                             domain(stat%idomain)%dccc%xsz(3))
            real(wp) :: coe1, coe2


            stat%nstatis = stat%nstatis + 1
            select case (stat%accumulation)
                case(ISTAT_ACCUM_ASSYM)
                    coe1 = dble(stat%nstatis - 1)/ dble(stat%nstatis)
                    coe2 = 1.0_wp/dble(stat%nstatis)

                    stat%data_array(:,:,:,:) = coe1*stat%data_array(:,:,:,:) &
                                                + coe2*avg(:,:,:,:)

                case default
                    stat%data_array(:,:,:,:) = avg
            end select

        end subroutine

        subroutine write_stat_xz(stat)
            type(t_stats_xz) :: stat

            character(len=256) :: file_name

            integer :: int_info(3), ashape(2)
            real(wp) :: rl_info(1)

            integer :: unit = 13, iostat

            if (.not.check_time_write(stat)) return 

            call create_file_name(stat,file_name)

            if (nrank.eq.0) then
                int_info(:) = [stat%accumulation,&
                               stat%nstatis,&
                               stat%xz_avg]

                ashape(:) = shape(stat%data_array)

                rl_info(:) = [flow(stat%idomain)%ren]

                open(unit,file=file_name, status='REPLACE',&
                            form='UNFORMATTED',action='WRITE',&
                            iostat=iostat, access='STREAM')

                if (iostat.ne.0) &
                    call Print_error_msg("There has been a problem writing file "//trim(file_name))

                write(unit) int_info
                write(unit) ashape
                write(unit) rl_info

                write(unit) stat%data_array

                close(unit)
            endif
        end subroutine

        subroutine write_stat_z(stat)
            type(t_stats_z) :: stat
            character(len=256) :: file_name

            integer :: int_info(3), ashape(3)
            real(wp) :: rl_info(1)

            integer :: sizes_array(3), subsizes(3), starts(3)
            integer :: unit = 13, iostat, write_id
            integer(MPI_OFFSET_KIND) :: offset
            integer :: mpi_type

            if (.not.check_time_write(stat)) return 

            call create_file_name(stat,file_name)

            if (nrank.eq.0) then
                int_info(:) = [stat%accumulation,&
                               stat%nstatis,&
                               stat%xz_avg]

                ashape(:) = [stat%size,&
                            domain(stat%idomain)%nc(1),&
                            domain(stat%idomain)%nc(2)]

                rl_info(:) = [flow(stat%idomain)%ren]

                open(unit,file=file_name, status='REPLACE',&
                            form='UNFORMATTED',action='WRITE',&
                            iostat=iostat, access='STREAM')

                if (iostat.ne.0) &
                    call Print_error_msg("There has been a problem writing file "//trim(file_name))

                write(unit) int_info
                write(unit) ashape
                write(unit) rl_info

                close(unit)

                offset = kind(int_info)*size(int_info) + &
                     kind(ashape)*size(ashape) + &
                     kind(rl_info)*size(rl_info)
            endif

            call MPI_Bcast(offset,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)

            sizes_array(1) = stat%size
            sizes_array(2) = domain(stat%idomain)%nc(1)
            sizes_array(3) = domain(stat%idomain)%nc(2)

            subsizes(1) = stat%size
            subsizes(2) = domain(stat%idomain)%dccc%zsz(1)
            subsizes(3) = domain(stat%idomain)%dccc%zsz(2)

            ! write(*,*) domain(stat%idomain)%nc(1)
            starts(1) = 0
            starts(2) = domain(stat%idomain)%dccc%zst(1) - 1
            starts(3) = domain(stat%idomain)%dccc%zst(2) - 1

            call  MPI_Type_create_subarray(3, sizes_array, subsizes, starts,&
                                        MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION,&
                                        mpi_type, IERROR)

            CALL MPI_TYPE_COMMIT(mpi_type, IERROR)

            call MPI_File_open(DECOMP_2D_COMM_CART_Z,&
                                file_name,MPI_MODE_WRONLY,&
                                 MPI_INFO_NULL,&
                                 write_id, ierror)

            call MPI_File_set_view(write_id, offset ,&
                                 MPI_DOUBLE_PRECISION,&
                                 mpi_type,"NATIVE",&
                                 MPI_INFO_NULL,&
                                 ierror)
            
            call MPI_File_write_all(write_id,&
                                    stat%data_array,&
                                    size(stat%data_array),&
                                    MPI_DOUBLE_PRECISION,&
                                    MPI_STATUS_IGNORE,&
                                    IERROR)

            call MPI_File_close(write_id, IERROR)

        end subroutine

        subroutine write_stat_N(stat)
            type(t_stats_N) :: stat
            
            character(len=256) :: file_name

            integer :: int_info(3), ashape(4)
            real(wp) :: rl_info(1)

            integer :: sizes_array(4), subsizes(4), starts(4)
            integer :: unit = 13, iostat, write_id
            integer(MPI_OFFSET_KIND) :: offset
            integer :: mpi_type

            if (.not.check_time_write(stat)) return 

            call create_file_name(stat,file_name)

            if (nrank.eq.0) then
                int_info(:) = [stat%accumulation,&
                               stat%nstatis,&
                               stat%xz_avg]

                ashape(:) = [stat%size,&
                            domain(stat%idomain)%nc(1),&
                            domain(stat%idomain)%nc(2),&
                            domain(stat%idomain)%nc(3)]

                rl_info(:) = [flow(stat%idomain)%ren]

                open(unit,file=file_name, status='REPLACE',&
                            form='UNFORMATTED',action='WRITE',&
                            iostat=iostat, access='STREAM')

                if (iostat.ne.0) &
                    call Print_error_msg("There has been a problem writing file "//trim(file_name))

                write(unit) int_info
                write(unit) ashape
                write(unit) rl_info

                close(unit)

                offset = kind(int_info)*size(int_info) + &
                     kind(ashape)*size(ashape) + &
                     kind(rl_info)*size(rl_info)
            endif

            call MPI_Bcast(offset,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)

            sizes_array(1) = stat%size
            sizes_array(2) = domain(stat%idomain)%nc(1)
            sizes_array(3) = domain(stat%idomain)%nc(2)
            sizes_array(4) = domain(stat%idomain)%nc(3)

            subsizes(1) = stat%size
            subsizes(2) = domain(stat%idomain)%dccc%xsz(1)
            subsizes(3) = domain(stat%idomain)%dccc%xsz(2)
            subsizes(4) = domain(stat%idomain)%dccc%xsz(3)

            ! write(*,*) domain(stat%idomain)%nc(1)
            starts(1) = 0
            starts(2) = domain(stat%idomain)%dccc%xst(1) - 1
            starts(3) = domain(stat%idomain)%dccc%xst(2) - 1
            starts(4) = domain(stat%idomain)%dccc%xst(3) - 1

            call  MPI_Type_create_subarray(4, sizes_array, subsizes, starts,&
                                        MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION,&
                                        mpi_type, IERROR)

            CALL MPI_TYPE_COMMIT(mpi_type, IERROR)

            call MPI_File_open(DECOMP_2D_COMM_CART_X,&
                                file_name,MPI_MODE_WRONLY,&
                                MPI_INFO_NULL,&
                                write_id, ierror)

            call MPI_File_set_view(write_id, offset ,&
                                MPI_DOUBLE_PRECISION,&
                                mpi_type,"NATIVE",&
                                MPI_INFO_NULL,&
                                ierror)
           
           call MPI_File_write_all(write_id,&
                                   stat%data_array,&
                                   size(stat%data_array),&
                                   MPI_DOUBLE_PRECISION,&
                                   MPI_STATUS_IGNORE,&
                                   IERROR)

           call MPI_File_close(write_id, IERROR)                                
        end subroutine

        subroutine restart_avg(stat,time)
            class(t_stats_base), intent(inout) :: stat
            real(wp), intent(in) :: time

            character(len=256) :: file_name
            logical :: exists

            call create_file_name(stat,file_name,time=time)

            inquire(file=file_name,exist=exists)

            if (.not.exists) &
                call Print_error_msg("File "//trim(file_name)&
                                    //" does not exist")
            
            select type (stat)
                type is (t_stats_xz)
                    call read_stat_xz(stat,file_name)

                type is (t_stats_z)
                    call read_stat_z(stat,file_name)

                type is (t_stats_N)
                    call read_stat_N(stat,file_name)
            end select
        end subroutine

        subroutine read_stat_xz(stat, file_name)
            type(t_stats_xz) :: stat
            character(len=*) :: file_name

            integer :: int_info(3), ashape(2)
            real(wp) :: rl_info(1)

            integer :: unit = 13, iostat, a_size

            if (nrank .eq. 0) then
                open(unit,file=file_name, status='old',&
                            form='UNFORMATTED',action='read',&
                            iostat=iostat, access='STREAM')
                
                read(unit) int_info
                read(unit) ashape
                read(unit) rl_info

                if (int_info(1).ne.stat%accumulation) &
                        call Print_error_msg("Different accumultion types on average "&
                                            // trim(stat%name))

                stat%nstatis  = int_info(2)

                if (.not.all(ashape.eq.shape(stat%data_array))) &
                        call Print_error_msg("Different shape in file "//trim(file_name)//" on average "&
                                            // trim(stat%name))

                read(unit) stat%data_array
            endif

            call MPI_Bcast(stat%nstatis,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
            
            a_size = size(stat%data_array)
            call MPI_Bcast(stat%data_array,a_size,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierror)

        end subroutine

        subroutine read_stat_z(stat, file_name)
            type(t_stats_z) :: stat
            character(len=*) :: file_name

            integer :: int_info(3), ashape(3), ashape_file(3)
            real(wp) :: rl_info(1)
            integer :: sizes_array(3), subsizes(3), starts(3)
            integer :: unit = 13, iostat, read_id
            integer(MPI_OFFSET_KIND) :: offset
            integer :: mpi_type


            if (nrank .eq. 0) then
                open(unit,file=file_name, status='old',&
                            form='UNFORMATTED',action='read',&
                            iostat=iostat, access='STREAM')
                
                read(unit) int_info
                read(unit) ashape_file
                read(unit) rl_info

                if (int_info(1).ne.stat%accumulation) &
                        call Print_error_msg("Different accumultion types on average "&
                                            // trim(stat%name))

                stat%nstatis  = int_info(2)

                ashape(:) = [stat%size,&
                            domain(stat%idomain)%nc(1),&
                            domain(stat%idomain)%nc(2)]

                if (.not.all(ashape.eq.ashape_file)) &
                        call Print_error_msg("Different shape in file "//trim(file_name)//" on average "&
                                            // trim(stat%name))

                close(unit)
                offset = kind(int_info)*size(int_info) + &
                     kind(ashape)*size(ashape) + &
                     kind(rl_info)*size(rl_info)

            endif

            call MPI_Bcast(stat%nstatis,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
            call MPI_Bcast(offset,1,MPI_INTEGER8,0,MPI_COMM_WORLD,ierror)

            sizes_array(1) = stat%size
            sizes_array(2) = domain(stat%idomain)%nc(1)
            sizes_array(3) = domain(stat%idomain)%nc(2)

            subsizes(1) = stat%size
            subsizes(2) = domain(stat%idomain)%dccc%zsz(1)
            subsizes(3) = domain(stat%idomain)%dccc%zsz(2)

            ! write(*,*) domain(stat%idomain)%nc(1)
            starts(1) = 0
            starts(2) = domain(stat%idomain)%dccc%zst(1) - 1
            starts(3) = domain(stat%idomain)%dccc%zst(2) - 1

            call  MPI_Type_create_subarray(3, sizes_array, subsizes, starts,&
                                        MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION,&
                                        mpi_type, IERROR)

            CALL MPI_TYPE_COMMIT(mpi_type, IERROR)

            call MPI_File_open(DECOMP_2D_COMM_CART_Z,&
                                file_name,MPI_MODE_RDONLY,&
                                 MPI_INFO_NULL,&
                                 read_id, ierror)

            call MPI_File_set_view(read_id, offset ,&
                                 MPI_DOUBLE_PRECISION,&
                                 mpi_type,"NATIVE",&
                                 MPI_INFO_NULL,&
                                 ierror)
            
            call MPI_File_read_all(read_id,&
                                 stat%data_array,&
                                 size(stat%data_array),&
                                 MPI_DOUBLE_PRECISION,&
                                 MPI_STATUS_IGNORE,&
                                 IERROR)

            call MPI_File_close(read_id, IERROR)

        end subroutine

        subroutine read_stat_N(stat, file_name)
            type(t_stats_N), intent(inout) :: stat
            character(len=*), intent(in) :: file_name

            integer :: int_info(3), ashape(4), ashape_file(4)
            real(wp) :: rl_info(1)
            integer :: sizes_array(4), subsizes(4), starts(4)
            integer :: unit = 13, iostat, read_id
            integer(MPI_OFFSET_KIND) :: offset
            integer :: mpi_type

            if (nrank .eq. 0) then
                open(unit,file=file_name, status='old',&
                            form='UNFORMATTED',action='read',&
                            iostat=iostat, access='STREAM')
                
                read(unit) int_info
                read(unit) ashape_file
                read(unit) rl_info

                if (int_info(1).ne.stat%accumulation) &
                        call Print_error_msg("Different accumultion types on average "&
                                            // trim(stat%name))

                stat%nstatis  = int_info(2)

                ashape(:) = [stat%size,&
                            domain(stat%idomain)%nc(1),&
                            domain(stat%idomain)%nc(2),&
                            domain(stat%idomain)%nc(3)]
                            
                if (.not.all(ashape.eq.ashape_file)) &
                        call Print_error_msg("Different shape in file "//trim(file_name)//" on average "&
                                            // trim(stat%name))


                close(unit)
                offset = kind(int_info)*size(int_info) + &
                     kind(ashape)*size(ashape) + &
                     kind(rl_info)*size(rl_info)

            endif

            call MPI_Bcast(offset,1,MPI_INTEGER8,0,MPI_COMM_WORLD,ierror)
            call MPI_Bcast(stat%nstatis,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)

            write(*,*) nrank, offset, MPI_OFFSET_KIND
            sizes_array(1) = stat%size
            sizes_array(2) = domain(stat%idomain)%nc(1)
            sizes_array(3) = domain(stat%idomain)%nc(2)
            sizes_array(4) = domain(stat%idomain)%nc(3)

            subsizes(1) = stat%size
            subsizes(2) = domain(stat%idomain)%dccc%xsz(1)
            subsizes(3) = domain(stat%idomain)%dccc%xsz(2)
            subsizes(4) = domain(stat%idomain)%dccc%xsz(3)

            starts(1) = 0
            starts(2) = domain(stat%idomain)%dccc%xst(1) - 1
            starts(3) = domain(stat%idomain)%dccc%xst(2) - 1
            starts(4) = domain(stat%idomain)%dccc%xst(3) - 1

            call  MPI_Type_create_subarray(4, sizes_array, subsizes, starts,&
                                        MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION,&
                                        mpi_type, IERROR)

            CALL MPI_TYPE_COMMIT(mpi_type, IERROR)

            call MPI_File_open(DECOMP_2D_COMM_CART_X,&
                                file_name,MPI_MODE_RDONLY,&
                                MPI_INFO_NULL,&
                                read_id, ierror)

            call MPI_File_set_view(read_id, offset ,&
                                MPI_DOUBLE_PRECISION,&
                                mpi_type,"NATIVE",&
                                MPI_INFO_NULL,&
                                ierror)
           
           call MPI_File_read_all(read_id,&
                                   stat%data_array,&
                                   size(stat%data_array),&
                                   MPI_DOUBLE_PRECISION,&
                                   MPI_STATUS_IGNORE,&
                                   IERROR)

            call MPI_File_close(read_id, IERROR)                                

        end subroutine

        subroutine create_file_name(stat,file_name,time)
            class(t_stats_base) :: stat
            character(len=256) :: file_name
            real(wp), optional :: time
            character(len=4) :: ndm
            character(len=15) :: t_stamp
            character(len=2) :: type

            write(ndm,'(I0)') stat%idomain

            if (present(time)) then
                write(t_stamp,'(1ES15.9)') time
            else
                write(t_stamp,'(1ES15.9)') flow(stat%idomain)%time
            endif
            select case(stat%xz_avg)
                case (ISTAT_AVG_XZ)
                    type = "XZ"
                case (ISTAT_AVG_Z)
                    type = "Z"
                case (ISTAT_AVG_N)
                    type = "N"
            end select

            file_name = "DNS_DOMAIN"//trim(ndm)//"_"//trim(type)//"_T"&
                    //trim(t_stamp)//"_"//trim(stat%name)//".D"


        end subroutine
end module

!> # Module chapsim_stats
!> \brief A high-level module for postprocessing statistics for CHAPSim2.0
!> Brief documentation is given for the public API

module chapsim_stats
    use stat_types
    use iso_fortran_env
    use operations

    implicit none

    private

    integer :: nstat_xz = 0
    integer :: nstat_z = 0
    integer :: nstat_N = 0

    type(t_stats_xz), allocatable :: chapsim_stats_xz(:)
    type(t_stats_z), allocatable  :: chapsim_stats_z(:) 
    type(t_stats_N), allocatable  :: chapsim_stats_N(:) 
    
    public :: chapsim_stats_init
    public :: chapsim_stats_finalize
    public :: chapsim_register_stat_xz
    public :: chapsim_register_stat_z
    public :: chapsim_register_stat_N
    public :: chapsim_stats_calculate
    public :: chapsim_avg_xz
    public :: chapsim_avg_z
    public :: chapsim_stats_write

    public :: ISTAT_ACCUM_ASSYM
    public :: ISTAT_ACCUM_LOCAL


    contains

        !> \brief Initialises chapsim_stats
        !-----------------------------------------------------------------------
        !  Arguments
        !_______________________________________________________________________
        !  mode         name    role
        !-----------------------------------------------------------------------
        !> \param[in]   nstats  Number of elements the statistic type arrays are
        !>                      allocated to. If exceeded the code with crash   
        !_______________________________________________________________________
        subroutine chapsim_stats_init()
            
            call allocate_stats()

        end subroutine
        
        !---------------------------------------------------------------
        !> \brief Finalises and deallocates arrays in chapsim_stats
        !________________________________________________________________
        subroutine chapsim_stats_finalize()
            call deallocate_stats
        end subroutine

        subroutine allocate_stats
         
            allocate(chapsim_stats_xz(10))
            allocate(chapsim_stats_z(10))
            allocate(chapsim_stats_N(10))


        end subroutine

        subroutine deallocate_stats
            integer :: i

            do i = 1, nstat_xz
                deallocate(chapsim_stats_xz(i)%data_array)
            enddo   

            do i = 1, nstat_z
                deallocate(chapsim_stats_z(i)%data_array)
            enddo

            do i = 1, nstat_N
                deallocate(chapsim_stats_N(i)%data_array)
            enddo

            deallocate(chapsim_stats_xz)
            deallocate(chapsim_stats_z)
            deallocate(chapsim_stats_N)
        end subroutine

        !----------------------------------------------------------------------
        !> \brief Registers information about statistics to be processed later.
        !----------------------------------------------------------------------
        ! Arguments
        !______________________________________________________________________
        ! mode          name                role
        !----------------------------------------------------------------------
        !> \param[in]   name               string containing the name of the 
        !>                                 variable
        !> \param[in]   idm                domain index for the statistics
        !> \param[in]   file_prefix        prefix of the file the stats will 
        !>                                 be saved to
        !> \param[in]   accumulation_type  Flag indicating assymptotic or local 
        !>                                 average
        !> \param[in]   xz_avg             Flag indicating the type of spatial 
        !>                                 average used
        !> \param[in]   stats_calc         Subroutine calculating the centered 
        !>                                 time local value of the statistic 
        !>                                 without space averaging
        !> \param[in]   dt_calc            Time interval between statistic 
        !>                                 calculations
        !> \param[in]   size               Size of the first index of the 
        !>                                 statistics
        !> \param[in]   dt_save            Time interval between file writes
        !______________________________________________________________________

        subroutine chapsim_register_stat_xz(name, idm, nsize,&
                                             accumulation_type, stats_calc, &
                                             dt_calc,dt_write,start_time,&
                                             restart,restart_time)

            character(len=*),                  intent(in) :: name
            integer,                           intent(in) :: nsize
            integer,                           intent(in) :: idm
            integer,                           intent(in) :: accumulation_type
            procedure(stat_calc_xz) :: stats_calc
            real(wp),                           intent(in) :: dt_calc
            real(wp),                           intent(in) :: dt_write
            real(wp), optional, intent(in) :: start_time
            logical, optional, intent(in) :: restart
            real(wp), optional, intent(in) :: restart_time

            logical :: restart_local
            real(wp) :: start

            if (check_list_size(chapsim_stats_xz,nstat_xz+1)) then
                call reallocate_stats_xz(nstat_xz+10)
            endif
            nstat_xz = nstat_xz + 1

            if (present(start_time)) then
                start = start_time
            else
                start = -1.0_wp
            endif

            call set_stat_info(chapsim_stats_xz(nstat_xz),&
                                        name, nsize, idm, &
                                        accumulation_type, ISTAT_AVG_XZ, &
                                        dt_calc,dt_write,start)

            call allocate_data_array(chapsim_stats_xz(nstat_xz))  
            
            if (present(restart)) then
                restart_local = restart
                if (.not.present(restart_time)) &
                    call chapsim_stats_error("if restart is present, a restart time must be present")
            else
                restart_local = .false.
            endif

            if (restart_local.and.restart_time.gt.start_time) &
                call restart_avg(chapsim_stats_xz(nstat_xz),restart_time)
            
            chapsim_stats_xz(nstat_xz)%stat_calc => stats_calc

        end subroutine

        subroutine chapsim_register_stat_z(name, idm, nsize,&
                                            accumulation_type, stats_calc, &
                                            dt_calc,dt_write,start_time,&
                                             restart,restart_time)

            character(len=*),                  intent(in) :: name
            integer,                           intent(in) :: nsize
            integer,                           intent(in) :: idm
            integer,                           intent(in) :: accumulation_type
            procedure(stat_calc_z) :: stats_calc
            real(wp),                           intent(in) :: dt_calc
            real(wp),                           intent(in) :: dt_write
            real(wp), optional, intent(in) :: start_time
            logical, optional, intent(in) :: restart
            real(wp), optional, intent(in) :: restart_time

            logical :: restart_local
            real(wp) :: start

            if (check_list_size(chapsim_stats_z,nstat_z+1)) then
                call reallocate_stats_z(nstat_z+10)
            endif
            nstat_z = nstat_z + 1

            if (present(start_time)) then
                start = start_time
            else
                start = -1.0_wp
            endif

            call set_stat_info(chapsim_stats_z(nstat_z),&
                    name, nsize, idm, &
                    accumulation_type, ISTAT_AVG_Z, &
                    dt_calc, dt_write,start)

            call allocate_data_array(chapsim_stats_z(nstat_z))   
            
            if (present(restart)) then
                restart_local = restart
                if (.not.present(restart_time)) &
                    call chapsim_stats_error("if restart is present, a restart time must be present")
            else
                restart_local = .false.
            endif

            if (restart_local.and.restart_time.gt.start_time) &
                call restart_avg(chapsim_stats_z(nstat_z),restart_time)

            chapsim_stats_z(nstat_z)%stat_calc => stats_calc

        end subroutine

        subroutine chapsim_register_stat_N(name, idm, nsize,&
                                            accumulation_type, stats_calc, &
                                            dt_calc,dt_write,start_time,&
                                             restart,restart_time)

            character(len=*),                  intent(in) :: name
            integer,                           intent(in) :: nsize
            integer,                           intent(in) :: idm
            integer,                           intent(in) :: accumulation_type
            procedure(stat_calc_N) :: stats_calc
            real(wp),                           intent(in) :: dt_calc
            real(wp),                           intent(in) :: dt_write
            real(wp), optional, intent(in) :: start_time
            logical, optional, intent(in) :: restart
            real(wp), optional, intent(in) :: restart_time

            logical :: restart_local
            real(wp) :: start

            if (check_list_size(chapsim_stats_N,nstat_N+1)) then
                call reallocate_stats_N(nstat_N+10)                        
            endif
            nstat_N = nstat_N + 1

            if (present(start_time)) then
                start = start_time
            else
                start = -1.0_wp
            endif

            call set_stat_info(chapsim_stats_N(nstat_N),&
                                name, nsize, idm, &
                                accumulation_type, ISTAT_AVG_N, &
                                dt_calc, dt_write,start)

            call allocate_data_array(chapsim_stats_N(nstat_N))

            if (present(restart)) then
                restart_local = restart
                if (.not.present(restart_time)) &
                    call chapsim_stats_error("if restart is present, a restart time must be present")
            else
                restart_local = .false.
            endif

            if (restart_local.and.restart_time.gt.start_time) &
                call restart_avg(chapsim_stats_N(nstat_N),restart_time)

            chapsim_stats_N(nstat_N)%stat_calc => stats_calc

        end subroutine

        subroutine chapsim_avg_xz(dm,data_array,avg_array)
            type(t_domain), intent(in) :: dm
            real(wp), intent(in) :: data_array(:,:,:)
            real(wp), intent(out) :: avg_array(:)

            real(wp), allocatable :: y_pencil_array(:,:,:)
            real(wp), allocatable :: avg_local(:)

            integer :: s1_y, s2_y, s3_y
            integer :: j, k, l
            integer :: reduce_size, comm

            real(wp) :: xz_div

            s1_y = dm%dccc%ysz(1)
            s2_y = dm%dccc%ysz(2)
            s3_y = dm%dccc%ysz(3)

            allocate(y_pencil_array(s1_y, s2_y, s3_y)) 
            allocate(avg_local(s2_y))

            avg_local(:) = ZERO

            call transpose_x_to_y(data_array,&
                                  y_pencil_array,&
                                  dm%dccc)

            xz_div = dble(s1_y*s3_y*nproc)

            do j = 1, s1_y
                do k = 1, s2_y
                    do l = 1, s3_y
                        avg_local(k) = avg_local(k) &
                                        + y_pencil_array(j,k,l)/xz_div
                    enddo
                enddo
            enddo

            reduce_size = s2_y

            call MPI_Allreduce(avg_local, avg_array,&
                                 reduce_size, &
                                 MPI_DOUBLE_PRECISION, &
                                 MPI_SUM,&
                                 DECOMP_2D_COMM_CART_Y,&
                                 ierror)

            deallocate(avg_local, y_pencil_array)
        end subroutine

        subroutine chapsim_avg_z(dm, data_array, avg_array)
            type(t_domain), intent(in) :: dm
            real(wp), intent(in)  :: data_array(:,:,:)
            real(wp), intent(out) :: avg_array(:,:)

            real(wp), allocatable :: z_pencil_array(:,:,:)
            real(wp), allocatable :: y_pencil_array(:,:,:)
            
            integer :: s1_y, s2_y, s3_y
            integer :: s1_z, s2_z, s3_z
            integer :: j, k, l
            real(wp) :: z_div    

            s1_z = dm%dccc%zsz(1)
            s2_z = dm%dccc%zsz(2)
            s3_z = dm%dccc%zsz(3)

            s1_y = dm%dccc%ysz(1)
            s2_y = dm%dccc%ysz(2)
            s3_y = dm%dccc%ysz(3)

            allocate(z_pencil_array(s1_z, s2_z, s3_z)) 
            allocate(y_pencil_array(s1_y, s2_y, s3_y)) 

            avg_array(:,:) = ZERO

            z_div = dble(s3_z)

            call transpose_x_to_y(data_array,&
                                 y_pencil_array,&
                                 dm%dccc)

            call transpose_y_to_z(y_pencil_array,&
                                  z_pencil_array,&
                                  dm%dccc)

            do j = 1, s1_z
                do k = 1, s2_z
                    do l = 1, s3_z
                        avg_array(j,k) = avg_array(j,k) &
                                            + z_pencil_array(j,k,l)/z_div
                    enddo
                enddo
            enddo    
            
            deallocate(y_pencil_array, z_pencil_array)
        end subroutine

        logical function check_list_size(stat,array_size)
            class(t_stats_base), dimension(:) :: stat
            integer :: array_size

            if (size(stat).lt.array_size) then
                check_list_size = .true.
            else
                check_list_size = .false.
            endif
        end function

        subroutine reallocate_stats_xz(array_size)
            integer, intent(in) :: array_size
            integer :: old_size

            type(t_stats_xz), allocatable :: stat_tmp(:)

            if (nstat_xz.gt.0) then
                allocate(stat_tmp(nstat_xz))

                stat_tmp(:) = chapsim_stats_xz(1:nstat_xz)
            endif

            deallocate(chapsim_stats_xz)
            allocate(chapsim_stats_xz(array_size))

            if (nstat_xz.gt.0) then
                chapsim_stats_xz(1:nstat_xz) = stat_tmp(:)
                deallocate(stat_tmp)
            endif

        end subroutine

        subroutine reallocate_stats_z(array_size)
            integer, intent(in) :: array_size
            integer :: old_size

            type(t_stats_z), allocatable :: stat_tmp(:)

            if (nstat_z.gt.0) then
                allocate(stat_tmp(nstat_z))

                stat_tmp(:) = chapsim_stats_z(1:nstat_z)
            endif

            deallocate(chapsim_stats_z)
            allocate(chapsim_stats_z(array_size))

            if (nstat_z.gt.0) then
                chapsim_stats_z(1:nstat_z) = stat_tmp(:)
                deallocate(stat_tmp)
            endif
            
        end subroutine

        subroutine reallocate_stats_N(array_size)
            integer, intent(in) :: array_size

            type(t_stats_N), allocatable :: stat_tmp(:)

            if (nstat_N.gt.0) then
                allocate(stat_tmp(nstat_N))

                stat_tmp(:) = chapsim_stats_N(1:nstat_N)
            endif

            deallocate(chapsim_stats_N)
            allocate(chapsim_stats_N(array_size))

            if (nstat_N.gt.0) then
                chapsim_stats_N(1:nstat_N) = stat_tmp(:)
                deallocate(stat_tmp)
            endif
            
        end subroutine

        !-------------------------------------------------------------
        !> \brief Calulcates the statistics if at the appropriate time
        !_____________________________________________________________
        subroutine chapsim_stats_calculate
            integer :: i

            do i = 1, nstat_xz
                call average_stats(chapsim_stats_xz(i))
            enddo

            do i = 1, nstat_z
                call average_stats(chapsim_stats_z(i))
            enddo

            do i = 1, nstat_N
                call average_stats(chapsim_stats_N(i))
            enddo


        end subroutine

        !------------------------------------------------------------
        !> \brief calls write subroutines on all statistics, 
        !> if time matches dt_write, statistic will write to file 
        !_____________________________________________________________
        subroutine chapsim_stats_write
            integer :: i

            do i = 1, nstat_xz
                call write_stat(chapsim_stats_xz(i))
            enddo

            do i = 1, nstat_z
                call write_stat(chapsim_stats_z(i))
            enddo

            do i = 1, nstat_N
                call write_stat(chapsim_stats_N(i))
            enddo
        end subroutine

        subroutine chapsim_stats_error(msg)
            character(len=*) :: msg

            call Print_error_msg("**CHAPSim statistics reported error: "&
                                    //trim(adjustl(msg)))
        end subroutine
end module chapsim_stats

module statistics
    use chapsim_stats
    use parameters_constant_mod
    use decomp_2d
    use mpi_mod
    use udf_type_mod
    use vars_df_mod
    use operations

    implicit none

    private

    public :: chapsim_user_stats
    public :: chapsim_stats_calculate
    public :: chapsim_stats_write
    public :: chapsim_stats_finalize

    contains
        !> \brief subroutine called from within chapsim
        !> Additional statistics should be added in here
        subroutine chapsim_user_stats

            real(wp) :: dt_save
            real(wp) :: dt_calc
            real(wp) :: start_time
            logical  :: restart
            real(wp) :: restart_time 

            dt_save = 2.0_wp
            dt_calc = 0.01_wp
            start_time = 0.0_wp
            restart = .false.
            restart_time = 0.0_wp

            call chapsim_stats_init()

            ! Here are some example statistics, subroutines defined below
            ! These are for x-z periodic flows

            ! [u, v, w, p]
            call chapsim_register_stat_xz('u',&
                                          1,&
                                          4,&
                                          ISTAT_ACCUM_ASSYM,&
                                          u_mean_calc_xz,&
                                          dt_calc,&
                                          dt_save,&
                                          restart=restart,&
                                          restart_time=restart_time,&
                                          start_time=start_time)
            ! [uu, vv, ww, uv, uw, vw]
            call chapsim_register_stat_xz('uu',&
                                          1,&
                                          6,&
                                          ISTAT_ACCUM_ASSYM,&
                                          uu_mean_calc_xz,&
                                          dt_calc,&
                                          dt_save,&
                                          restart=restart,&
                                          restart_time=restart_time,&
                                          start_time=start_time)


        end subroutine

        !-------------------------------------------------------------------
        ! place subroutines below here
        !-------------------------------------------------------------------            

        ! example subroutine for creating a vector [u, v, w, p]^T
        subroutine u_mean_calc_xz(fl,dm,size, array)
            type(t_flow) :: fl
            type(t_domain) :: dm
            integer :: size
            real(wp), dimension(size, &
                                dm%dccc%ysz(2)) :: array

            real(wp), allocatable :: x_pencil_ccc(:,:,:)
            real(wp), allocatable :: y_pencil_cpc(:,:,:)
            real(wp), allocatable :: y_pencil_ccc(:,:,:)
            real(wp), allocatable :: y_pencil_ccp(:,:,:)
            real(wp), allocatable :: z_pencil_ccp(:,:,:)
            real(wp), allocatable :: z_pencil_ccc(:,:,:)

            real(wp), allocatable :: array_tmp1(:)

            allocate(x_pencil_ccc(dm%dccc%xsz(1),dm%dccc%xsz(2), dm%dccc%xsz(3)))
            allocate(y_pencil_cpc(dm%dcpc%ysz(1),dm%dcpc%ysz(2), dm%dcpc%ysz(3)))
            allocate(y_pencil_ccc(dm%dccc%ysz(1),dm%dccc%ysz(2), dm%dccc%ysz(3)))
            allocate(y_pencil_ccp(dm%dccp%ysz(1),dm%dccp%ysz(2), dm%dccp%ysz(3)))
            allocate(z_pencil_ccp(dm%dccp%zsz(1),dm%dccp%zsz(2), dm%dccp%zsz(3)))
            allocate(z_pencil_ccc(dm%dccc%zsz(1),dm%dccc%zsz(2), dm%dccc%zsz(3)))

            allocate(array_tmp1(dm%dccc%ysz(2)))

            if (size.ne.4) call Print_error_msg("Size of array incorrect")
            
            ! u_mean 
            call Get_x_midp_P2C_3D(fl%qx,x_pencil_ccc,dm,dm%ibcx(:,1))

            call chapsim_avg_xz(dm,x_pencil_ccc,array_tmp1)

            array(1,:) = array_tmp1

            ! v mean
            call transpose_x_to_y(fl%qy, y_pencil_cpc,dm%dcpc)
            call Get_y_midp_P2C_3D(y_pencil_cpc,y_pencil_ccc,dm,dm%ibcy(:,2))
            call transpose_y_to_x(y_pencil_ccc, x_pencil_ccc, dm%dccc)

            call chapsim_avg_xz(dm,x_pencil_ccc,array_tmp1)

            array(2,:) = array_tmp1

            ! w mean

            call transpose_x_to_y(fl%qz, y_pencil_ccp,dm%dccp)
            call transpose_y_to_z(y_pencil_ccp, z_pencil_ccp,dm%dccp)

            call Get_z_midp_P2C_3D(z_pencil_ccp,z_pencil_ccc,dm,dm%ibcz(:,3))

            call transpose_z_to_y(z_pencil_ccc, y_pencil_ccc,dm%dccc)
            call transpose_y_to_x(y_pencil_ccc, x_pencil_ccc,dm%dccc)

            call chapsim_avg_xz(dm,x_pencil_ccc,array_tmp1)
            array(3,:) = array_tmp1

            ! p mean

            call chapsim_avg_xz(dm,fl%pres,array_tmp1)
            array(4,:) = array_tmp1

            deallocate(x_pencil_ccc)
            deallocate(y_pencil_cpc)
            deallocate(y_pencil_ccc)
            deallocate(y_pencil_ccp)
            deallocate(z_pencil_ccp)
            deallocate(z_pencil_ccc)
            deallocate(array_tmp1)

        end subroutine u_mean_calc_xz

        ! example subroutine for creating a array [uu, vv, ww, uv, uw, vw]
        subroutine uu_mean_calc_xz(fl,dm,size, array)
            type(t_flow) :: fl
            type(t_domain) :: dm
            integer :: size
            real(wp), dimension(size, &
                                dm%dccc%ysz(2)) :: array
            
            real(wp), allocatable :: u_array(:,:,:)
            real(wp), allocatable :: v_array(:,:,:)
            real(wp), allocatable :: w_array(:,:,:)

            real(wp), allocatable ::uu_work(:,:,:)

            real(wp), allocatable :: y_pencil_cpc(:,:,:)
            real(wp), allocatable :: y_pencil_ccc(:,:,:)
            real(wp), allocatable :: y_pencil_ccp(:,:,:)
            real(wp), allocatable :: z_pencil_ccp(:,:,:)
            real(wp), allocatable :: z_pencil_ccc(:,:,:)

            real(wp), allocatable :: array_tmp1(:)

            allocate(u_array(dm%dccc%xsz(1),dm%dccc%xsz(2), dm%dccc%xsz(3)))
            allocate(v_array(dm%dccc%xsz(1),dm%dccc%xsz(2), dm%dccc%xsz(3)))
            allocate(w_array(dm%dccc%xsz(1),dm%dccc%xsz(2), dm%dccc%xsz(3)))
            allocate(uu_work(dm%dccc%xsz(1),dm%dccc%xsz(2), dm%dccc%xsz(3)))

            allocate(y_pencil_cpc(dm%dcpc%ysz(1),dm%dcpc%ysz(2), dm%dcpc%ysz(3)))
            allocate(y_pencil_ccc(dm%dccc%ysz(1),dm%dccc%ysz(2), dm%dccc%ysz(3)))
            allocate(y_pencil_ccp(dm%dccp%ysz(1),dm%dccp%ysz(2), dm%dccp%ysz(3)))
            allocate(z_pencil_ccp(dm%dccp%zsz(1),dm%dccp%zsz(2), dm%dccp%zsz(3)))
            allocate(z_pencil_ccc(dm%dccc%zsz(1),dm%dccc%zsz(2), dm%dccc%zsz(3)))

            allocate(array_tmp1(dm%dccc%ysz(2)))

            if (size.ne.6) call Print_error_msg("Size of array incorrect")

            ! u
            call Get_x_midp_P2C_3D(fl%qx,u_array,dm,dm%ibcx(:,1))

            ! v
            call transpose_x_to_y(fl%qy, y_pencil_cpc,dm%dcpc)
            call Get_y_midp_P2C_3D(y_pencil_cpc,y_pencil_ccc,dm,dm%ibcy(:,2))
            call transpose_y_to_x(y_pencil_ccc, v_array, dm%dccc)

            ! w

            call transpose_x_to_y(fl%qz, y_pencil_ccp,dm%dccp)
            call transpose_y_to_z(y_pencil_ccp, z_pencil_ccp,dm%dccp)

            call Get_z_midp_P2C_3D(z_pencil_ccp,z_pencil_ccc,dm,dm%ibcz(:,3))

            call transpose_z_to_y(z_pencil_ccc, y_pencil_ccc,dm%dccc)
            call transpose_y_to_x(y_pencil_ccc, w_array,dm%dccc)

            ! uu
            uu_work(:,:,:) = u_array(:,:,:)*u_array(:,:,:)
            call chapsim_avg_xz(dm,uu_work,array_tmp1)
            array(1,:) = array_tmp1(:)
            
            ! vv
            uu_work(:,:,:) = v_array(:,:,:)*v_array(:,:,:)
            call chapsim_avg_xz(dm,uu_work,array_tmp1)
            array(2,:) = array_tmp1(:)

            ! ww
            uu_work(:,:,:) = w_array(:,:,:)*w_array(:,:,:)
            call chapsim_avg_xz(dm,uu_work,array_tmp1)
            array(3,:) = array_tmp1(:)

            ! uv
            uu_work(:,:,:) = u_array(:,:,:)*v_array(:,:,:)
            call chapsim_avg_xz(dm,uu_work,array_tmp1)
            array(4,:) = array_tmp1(:)

            ! uw
            uu_work(:,:,:) = u_array(:,:,:)*w_array(:,:,:)
            call chapsim_avg_xz(dm,uu_work,array_tmp1)
            array(5,:) = array_tmp1(:)

            ! vw
            uu_work(:,:,:) = v_array(:,:,:)*w_array(:,:,:)
            call chapsim_avg_xz(dm,uu_work,array_tmp1)
            array(6,:) = array_tmp1(:)

            deallocate(u_array)
            deallocate(v_array)
            deallocate(w_array)
            deallocate(uu_work)

            deallocate(y_pencil_cpc)
            deallocate(y_pencil_ccc)
            deallocate(y_pencil_ccp)
            deallocate(z_pencil_ccp)
            deallocate(z_pencil_ccc)

            deallocate(array_tmp1)

        end subroutine
end module statistics
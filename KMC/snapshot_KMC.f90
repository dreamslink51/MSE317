      program snapshot_KMC
      use kinetic_monte_carlo
      use overlapFINAL


      implicit none
      real(8), allocatable :: atcoord(:,:), atcoord2(:,:), epots(:)
      real(8), allocatable :: com(:,:), megacom(:,:)
      real(8), allocatable :: mass(:), coeff(:,:), jacob(:,:)
      real(8), allocatable :: karray(:,:), mol(:,:)
      integer, allocatable :: neighblistpre(:,:), neighblist(:,:)
      integer, allocatable :: connectlist(:,:), neworder(:)
      real*8, allocatable :: distarraypre(:,:,:), distarray(:,:,:)
      real(8), dimension(3,3) :: vectmatrix
      real(8), dimension(3) :: abc, angles, fake, Ein
      integer :: nat, nmol, ntot, i, length=3, elposold, elposnew 
      integer :: elposinit, rom, ierr, jmax, incr=1, ioerr, j, k, termin
      integer :: degen, ntot2, snapcount, scrap, aoinum, ioread, snapmax
      integer, allocatable :: atomlist(:)
      real(8) :: pi, distance(3)=0, tmig=0, v_drift(3), v_drift_ave(3)=0 
      real(8) :: finish, start, timing=0, timbit, size
      character(len=100) :: thesnapshot, geominput, junk, Epotfile, junk2
      character(len=100) :: coeffile, aoifile, ordfile, inputfile
      character(len=1) :: elemchar
      logical :: isepot, isreord, isouter
      integer :: snapfreq, nKMC, snapnum
      real*8  :: tmax, lambdai, omega, cutoff, T, epss, epsop, rVdW
      real*8, allocatable :: tempcoord(:,:)

!     Separating an input section to generalise the programme

      !geominput='PCBM-108_MD_NVT_short.xyz'
      !ordfile='neworder'
      !Epotfile='epottots'
      !coeffile='coeffpcbm'
      !aoifile='pcbm.aoi'

!     Lattice parameters
      !abc=(/39.87, 44.79, 58.44/) 
      !angles=(/90.0, 180-106.9, 90.0/)

      !nat=88 ! number of atoms per molecule
      !aoinum=60 ! number of relevant atoms per molecule
      !nmol=108  ! number of molecules per snapshot
      !degen=2   ! degeneracy

      !print*,'Please type in the name of the input file'
 
      !read(*,*) inputfile
 
      inputfile='input.inp'

      epss = 0.0D0
      epsop = 0.0D0
      rVdW = 0.0D0

      open(unit=21, file=inputfile, status='old', action='read', iostat=ioerr)
      if (ioerr .ne. 0) then
        print*, 'Could not open file', inputfile
        stop
      else
         read(21,*, iostat=ioread) geominput
         read(21,*, iostat=ioread) snapmax
         read(21,*, iostat=ioread) isepot
         if(isepot .eqv. .true.) read(21,*, iostat=ioread) Epotfile
         read(21,*, iostat=ioread) coeffile
         read(21,*, iostat=ioread) aoifile
         read(21,*, iostat=ioread) abc(1), abc(2), abc(3)
         read(21,*, iostat=ioread) angles(1), angles(2), angles(3)
         read(21,*, iostat=ioread) degen
         read(21,*, iostat=ioread) isreord
         if(isreord .eqv. .true.) read(21,*, iostat=ioread) ordfile
         read(21,*, iostat=ioread) cutoff
         read(21,*, iostat=ioread) lambdai
         read(21,*, iostat=ioread) omega
         read(21,*, iostat=ioread) Ein(1), Ein(2), Ein(3)
         read(21,*, iostat=ioread) T
         read(21,*, iostat=ioread) tmax
         read(21,*, iostat=ioread) nKMC
         read(21,*, iostat=ioread) snapfreq
         read(21,*, iostat=ioread) isouter
         if(isouter .eqv. .true.) read(21,*, iostat=ioread) epss, epsop, rVdW
      end if
      close(21)

      open(unit=33, file=geominput, status='old', action='read',iostat=ioerr)

      if (ioerr .ne. 0) then

          print*, 'Could not open file:', geominput
          stop

      else
          read(33,*) ntot

          close(33) 

      end if

!     Finding out number of atoms in single molecule based on coefficient file 

      open(unit=11, file=coeffile, status='old', action='read', iostat=ioerr)
      if (ioerr .ne. 0) then
        print*, 'Could not open file', coeffile
        stop
      else
        nat = 0
        ioread = 0
        do while (ioread .eq. 0)
          read(11,*, iostat=ioread) junk
          nat = nat + 1
        end do
        nat = nat - 1
        close(11)
      end if

      nmol=ntot/nat

!     Finding out number of atom of interests based on aoi file

      open(unit=11, file=aoifile, status='old', action='read', iostat=ioerr)
      if (ioerr .ne. 0) then
        print*, 'Could not open file', aoifile
        stop
      else
        aoinum = 0
        ioread = 0
        do while (ioread .eq. 0)
          read(11,*, iostat=ioread) junk
          aoinum = aoinum + 1
        end do
        aoinum = aoinum - 1
        close(11)
      end if
 
      open(unit=11, file=aoifile, status='old', action='read', iostat=ioerr)
      if (ioerr .ne. 0) then
        print*, 'Could not open file', aoifile
        stop
      else
        aoinum = 0
        ioread = 0
        do while (ioread .eq. 0)
          read(11,*, iostat=ioread) junk
          aoinum = aoinum + 1
        end do
        aoinum = aoinum - 1
        close(11)
      end if

      
      allocate(mass(nat))

      mass=12
      pi=4.d0*datan(1.d0)
      angles(:)=angles(:)*pi/180

!     Geometry trajectory file (xyz format)

      open(unit=33, file=geominput, status='old', action='read',iostat=ioerr)

      if (ioerr .ne. 0) then

          print*, 'Could not open file:', geominput
          stop

      end if

      if(isepot .eqv. .true.) then

!         Potential energies for each snapshot

          open(unit=34, file=Epotfile, status='old', action='read',iostat=ioerr)


          if (ioerr .ne. 0) then
              print*, 'Could not open file:', Epotfile
              stop
          end if

      end if

      call lattice_vectors(abc, angles, vectmatrix)

      allocate(com(3,nmol),stat=ierr)
      if (ierr .ne. 0) then
          write(*,*) 'Memory problem at com'
          stop
      end if

      allocate(connectlist(4,aoinum),stat=ierr)
      if (ierr .ne. 0) then
          write(*,*) 'Memory problem at com'
          stop
      end if

      allocate(atomlist(aoinum),stat=ierr)
      if (ierr .ne. 0) then
          write(*,*) 'Memory problem at atcoord'
          stop
      end if
            
      allocate(atcoord(4,ntot),stat=ierr)
      if (ierr .ne. 0) then
          write(*,*) 'Memory problem at atcoord'
          stop
      end if

      allocate(tempcoord(4,nat),stat=ierr)
      if (ierr .ne. 0) then
          write(*,*) 'Memory problem at atcoord'
          stop
      end if


      if(isreord .eqv. .true.) then

          allocate(neworder(ntot),stat=ierr)
          if (ierr .ne. 0) then
              write(*,*) 'Memory problem at atcoord'
              stop
          end if

      end if

      allocate(megacom(3,nmol*27), stat=ierr)
      if (ierr .ne. 0) then
          write(*,*) 'Memory problem at megacom'
          stop
      end if

      allocate(coeff(nat,degen),stat=ierr)
      if (ierr .ne. 0) then
          write(*,*) 'Memory problem at coeff'
          stop
      end if

      allocate(epots(nmol),stat=ierr)
      if (ierr .ne. 0) then
          write(*,*) 'Memory problem at epots'
          stop
      end if

      epots=0.0D0

!     Read in coefficients 
      open(unit=12, file=coeffile, status='old', action='read',iostat=ioerr)
      if (ioerr .ne. 0) then
          print*, 'Could not open file coeffile', coeffile
          stop
      else
        
          do i=1,nat
              read(12,*) (coeff(i,k), k=1,degen)
          end do
      end if
      close(12)

!     Read in AOI 
      open(unit=12, file=aoifile, status='old', action='read',iostat=ioerr)
      if (ioerr .ne. 0) then

          print*, 'Could not open file coeffile', aoifile
          stop
      else
          do i=1,aoinum
              read(12,*) atomlist(i)
          end do
      end if
      close(12)


!     This is only needed if the atoms in the molecules need reordering
!     for orbital coefficient analysis reasons.

      if(isreord .eqv. .true.) then

          open(unit=12, file=ordfile, status='old', action='read',iostat=ioerr)
          if (ioerr .ne. 0) then
              print*, 'Could not open file', ordfile
              stop
          else

              do i=1,nat
                 read(12,*) neworder(i)
              end do
          end if
          close(12)

!         Creating new order for every molecule in the snapshot

          do i=1,nmol-1

              neworder((nat*i)+1:nat*(i+1))=neworder(1:nat)+i*nat

          end do

      end if

      snapnum=0

!     Output files for certain data elelments  

      open(unit=43, file='rates.dat', action='write')
      open(unit=44, file='tcheck.dat', action='write')
      !open(unit=45, file='distance1.dat', action='write')
      !open(unit=46, file='distance2.dat', action='write')
      !open(unit=47, file='distance3.dat', action='write')
      !open(unit=48, file='distance4.dat', action='write')
      !open(unit=49, file='distance5.dat', action='write')
      !open(unit=50, file='distance6.dat', action='write')
      !open(unit=51, file='distance7.dat', action='write')
      !open(unit=52, file='distance8.dat', action='write')
      !open(unit=53, file='distance9.dat', action='write')
      open(unit=54, file='distance.dat', action='write')

!     Looping all over the trajectory reading it line by line.

      do snapcount=1,snapmax

          read(33,*) ntot2 

          if(ntot2 .ne. ntot) then
              write(*,*) 'Atom number mismatch check input files'
              stop
          end if

          read(33,*) junk !comment line in xyz file

          do j=1,ntot2
        
               if(isreord .eqv. .true.) then
   
                   i=neworder(j)
              
               else

                   i=j

               end if

               read(33,*) elemchar, atcoord(2,i), atcoord(3,i), atcoord(4,i)
  
               if     (elemchar .eq. 'C') then
                       atcoord(1,i) = 6
               elseif (elemchar .eq. 'N') then
                       atcoord(1,i) = 7
               elseif (elemchar .eq. 'O') then
                       atcoord(1,i) = 8
               elseif (elemchar .eq. 'H') then
                       atcoord(1,i) = 1
               else
               print*, 'Element not yet implemented'
              stop
              end if

          end do
 
          !do i=1,size(coeff,1)

          !    print*, (atcoord(k, i), k=1,4) , coeff(i,1)

          !end do

!         Calculating connectivity list for the first molecule    

          if (snapcount .eq. 1) then

             call connect_list(atcoord(2:4, 1:nat)*1.889725989D0, atomlist, connectlist)

             tempcoord(1,:)=atcoord(1,1:nat)
             tempcoord(2:4,:)=atcoord(2:4,1:nat)*1.889725989D0

!         Oh glorious first normalisation, only useful if the structures are idealised


             do i=1,degen

                call calc_sab(tempcoord, connectlist, coeff(:,i))

             end do

          end if

          !do i=1,size(coeff,1)
 
          !    print*, (atcoord(k, i), k=1,4) , coeff(i,1)

          !end do

          mass=atcoord(1,1:nat)

          if(isepot .eqv. .true.) then

              !print*, 'vomit'
              !print*, Epotfile

              read(34,*) (epots(k), k=1,nmol)

          end if

!         Only run KMC simulation on certain snapshots
        
          if(mod(snapcount,snapfreq) .eq. 1) then


              snapnum=snapnum+1
  
               write(43,*) 'Begin snapshot', snapcount 
              ! write(44,*) 'Begin snapshot', snapcount
        
              do i=0,8

!                  write(45+i,*) 'Begin snapshot', snapcount

              end do

               call full_KMC(atcoord, nat, mass, nmol, vectmatrix, coeff,&
&                            epots, connectlist, atomlist, lambdai, &
&                             omega, tmax, nKMC, cutoff, Ein, T, isouter,&
&                             epss, epsop, rVdW)


          end if

      end do

      call cpu_time(finish)

      timing=finish-start

       open(unit=29, file='environment.dat', action='write')

!     Print external data for post-processing

       write(29,*) T, Ein(1), Ein(2), Ein(3), nmol, snapnum


       close(29)


      write(*,*) 'Here is the timing', timing

      deallocate(com, atcoord, mass, megacom, epots)

      if(isepot .eqv. .true.) then

         close(34)

      end if

      close(33)
      close(43)
      close(44)
      close(54)

      end program
      
     

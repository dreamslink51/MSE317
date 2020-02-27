      module kinetic_monte_carlo

      implicit none
      
      contains

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                                !
      !               SUBROUTINES FOR TRAJECTORY FILES                 !
      !                                                                !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine read_snapshots(filename, atcoords)

      real(8), intent(out) :: atcoords(:,:)
      character(len=100) :: filename
      integer :: j, ioerr, nrows
      integer :: termin

      nrows=size(atcoords, 2)

      open(unit=22,file=filename, action='read', iostat=ioerr)

      if(ioerr .ne. 0) then
       write(*,*) 'No snapshot error error'
       stop
      else
         do j=1,nrows
            read(22,*) termin, atcoords(1,j), atcoords(2,j),atcoords(3,j)
         end do
      end if
      close(22)
      return

      end subroutine


      subroutine read_traj_file(nat,nmol,atcoords)

      ! reads in the trajectory file sent to me by Jacob needs a  bit of
      ! polishing

      ! nat : number of atoms per molecule
      ! nmol : number of molecules per unit cell
      ! atcoords: 3 x N array of atomic coordinates in the entire unit cell      


      implicit none

      real(8), intent(out) :: atcoords(:,:)
      integer :: i, j, nrows, ioerr
      integer, intent(in) :: nat, nmol
      open(unit=22,file='C60-108_sample.mdcrd',action='read', iostat=ioerr)
!      open(unit=22,file='mdcrd1',action='read', iostat=ioerr)


      if(ioerr .ne. 0) then
       write(*,*) 'No mdcrd file error error'
       stop
      else
         nrows=floor(dble(nat*nmol)/10)
!         write(*,*) 'This is nrows', nrows
         do i=1,nrows
            j=((i-1)*10)+1
            read(22,*) atcoords(1,j), atcoords(2,j),atcoords(3,j),&
&                    atcoords(1,j+1), atcoords(2,j+1),atcoords(3,j+1),&
&                    atcoords(1,j+2), atcoords(2,j+2),atcoords(3,j+2),&
&                    atcoords(1,j+3)
            read(22,*) atcoords(2,j+3),atcoords(3,j+3),&
&                    atcoords(1,j+4),&
&                    atcoords(2,j+4), atcoords(3,j+4),atcoords(1,j+5),&
&                    atcoords(2,j+5), atcoords(3,j+5),atcoords(1,j+6),&
&                    atcoords(2,j+6)
            read(22,*) atcoords(3,j+6),atcoords(1,j+7),&
&                    atcoords(2,j+7),&
&                    atcoords(3,j+7), atcoords(1,j+8),atcoords(2,j+8),&
&                    atcoords(3,j+8), atcoords(1,j+9),atcoords(2,j+9),&
&                    atcoords(3,j+9)
         end do
         
      end if
      close(22)
      return
      end subroutine

      subroutine moleculechopper(fullatcoord, reduced, nat, scrap)

!     In case one does not need the entirity of the molecule for the analysis (e.g. PCBM)
!     this subroutin chops off the first 'scrap' atom of each molecule

      implicit none

      real*8, intent(in) :: fullatcoord(:,:)
      real*8, intent(out) :: reduced(:,:)
      integer, intent(in) :: scrap, nat
      integer :: i, k
      real*8 :: size

      k=1

      do i=1,size(fullatcoord,2)

         !write(*,*) mod(i,nat)

         if (mod(i,nat) > scrap .or. mod(i,nat) .eq. 0) then
       
            reduced(:,k) = fullatcoord(:,i)
            k = k+1
           
         end if
            
      end do


      end subroutine
      

      subroutine centre_of_mass_old(xyz, natpermol, mass, com)

      !calculate and returns the centre of masses of each molecule in the slab
      !completely general nmol can be different from nmol earlier (e.g. PCBM and
      !other such molecules)

      !xyz : atomic coordinates for the entire unit cell
      !natpermol : numberof atoms per molecules
      !mass : a vector od the mass of the atoms
      !com  : output vector for centre of mass

      implicit none

      real(8), dimension(:,:), intent(in) :: xyz
      real(8), intent(in) :: mass(:)
      real(8), intent(out) :: com(:,:)
      real(8), allocatable :: current(:,:)
 

      integer, intent(in) :: natpermol
      integer :: i, nat, nmol
      integer :: stat_err

      
      nat = size(xyz,2)
      nmol=nat/natpermol

      allocate(current(3,natpermol))

      do i=1,nmol
        current=xyz(:,((i-1)*natpermol+1):(i*natpermol))
        call DGEMV('n', 3, natpermol, 1.0d0, current, 3, mass, 1, 0.0d0, com(:,i), 1)

        com(:,i) = com(:,i) / sum(mass)
      end do
      deallocate(current)

      !print*, sum(mass)

      return
      end subroutine

      subroutine centre_of_mass(xyz, atomlist, natpermol, mass, com)

!     Calculates and returns the centre of masses of each molecule in the slab.
!     Takes atomlist and only uses the relevant atoms to calculate 'centre of mass' 

      !xyz : atomic coordinates for the entire unit cell
      !natpermol : numberof atoms per molecules
      !mass : a vector od the mass of the atoms
      !com  : output vector for centre of mass

      implicit none

      real(8), dimension(:,:), intent(in) :: xyz
      real(8), intent(in) :: mass(:)
      integer, intent(in) :: atomlist(:)
      real(8), intent(out) :: com(:,:)
      real(8), allocatable :: current(:,:)
      real(8), allocatable :: relmol(:,:)
      real(8), allocatable :: relmass(:)


      integer, intent(in) :: natpermol
      integer :: i, j, nat, nmol
      integer :: stat_err


      nat = size(xyz,2)
      nmol=nat/natpermol

      allocate(relmol(3,size(atomlist)))
      allocate(relmass(size(atomlist)))
      allocate(current(3,natpermol))

      do i=1,nmol

         current=xyz(:,((i-1)*natpermol+1):(i*natpermol))

         do j=1,size(atomlist)
       
            relmol(:,j) = current(:,atomlist(j))
 
            relmass(j)  = mass(atomlist(j))       
   
         end do

         call DGEMV('n', 3, size(atomlist), 1.0d0, relmol, 3, relmass, 1, 0.0d0, com(:,i), 1)

         com(:,i) = com(:,i) / sum(relmass)
      end do

      deallocate(current)
      return
      end subroutine



      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                                           !
      !                           SUBROUTINES FOR PBC                             !
      !                                                                           !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine lattice_vectors(abc, angles, vectmatrix)
      ! Creates Cartesian vectors from the length and angularity of lattice  vectors. 
      ! The conversion is a,b,c, alpha beta gamma -> a(ax,ay,az), b(bx,by,bz), c(cx,cy,cz)

      implicit none
      
      real(8), intent(in) :: abc(:), angles(:)
      real(8), intent(out) :: vectmatrix(:,:)
      real(8) :: cosa, cosb, cosg, sing

      cosa = cos(angles(1))
      cosb = cos(angles(2))
      cosg = cos(angles(3))
      sing = sin(angles(3))
      ! lattice vector a
      vectmatrix(1,1) = abc(1) !ax
      vectmatrix(2,1) = 0      !ay 
      vectmatrix(3,1) = 0      !az 
      ! lattice vector b
      vectmatrix(1,2) = abc(2)*cosg !bx
      vectmatrix(2,2) = abc(2)*sing !by 
      vectmatrix(3,2) = 0           !bz 
      ! lattice vector c
      vectmatrix(1,3) = abc(3)*cosb !cx
      vectmatrix(2,3) = abc(3)*(cosa - cosb * cosg)/sing !cy 
      vectmatrix(3,3) = abc(3)*(sqrt(1 - cosa**2 - cosb**2 - cosg**2 &
&                       + 2*(cosa*cosb*cosg)))/sing      !cz  
      
      return
      end subroutine
      


      subroutine com_pbc(com, megacom, vectmatrix)

      !This soubroutine creates the 26 periodic images of the initial cell 
      !to look for all possible neighbours
     
      implicit none

      integer :: array(3), nmol, i, k, l, m, theindex, ierr
      real*8, intent(in) :: vectmatrix(:,:), com(:,:)
      real*8, intent(out) :: megacom(:,:)
   
      nmol = size(com, 2)

      array=(/0, 1, -1/)

!     megacom is the middle cell and the neighbouring 26 the centre of mass points
!     are filled in with the centre of masses form the middle cell using translation
!     along the lattice vectors

      megacom = 0

      do l =1,3
      do k =1,3
      do m =1,3

      do i=1,nmol

      theindex = ((9*(l-1)+3*(k-1)+(m-1))*nmol+i)

      megacom(:, theindex) = &
&     com(:,i)+array(m)*vectmatrix(:, 1)+array(k)*vectmatrix(:,2)&
&     +array(l)*vectmatrix(:,3)

      end do
      end do
      end do
      end do

 !      do i=1,size(megacom,2)

  !       write(*,*) megacom(1,i), megacom(2,i), megacom(3,i)

   !   end do


      return

      end subroutine

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                                            !
      !     THESE SUBROUTINES ARE HERE TO CREATE THE STATIC NEIGHBOUR              !
      !                   LIST AND THE CORRESPONDING RATES                         !
      !                                                                            !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      subroutine static_neighblist(com, megacom, neighblist, jmax, distarray, cutoff)

      implicit none

      real*8, intent(in) :: com(:,:), megacom(:,:)
      integer, intent(out) :: neighblist(:,:)
      real*8, intent(out) :: distarray(:,:,:)
      real*8 :: distvec(3), posvec(3), DNRM2
      integer :: i, nmol, j, k
      integer, intent(out) :: jmax
      real*8, allocatable :: deltaEad(:)
      real(8) :: cutoff

      nmol = size(com,2)


      do k=1,nmol

            posvec(:) = com(:,k)

            j=0

            do i=1,nmol*27

            distvec(:) = megacom(:,i) - posvec(:)

            !write(*,*) 'this is j',j


!     Identifying nearest neighbours anything that is closer than 'cutoff' angs is
!     stored so as its index

            if (DNRM2(3,distvec,1) < cutoff .and. k .ne. i) then
                         j = j+1
             neighblist(j+1,k) = i
             distarray(:,j,k)= distvec(:)

            end if

           end do

!     The fisrst lin contains the number of the neighbours for each molecule (will be 
!     important in amorphous systems) 

            neighblist(1,k)=j

            if (j > jmax) jmax=j

      end do

      return

      end subroutine

!     Old version of ratelist creator using radial p_pi directions

      subroutine ratelist_creator(distarray, neighblist, karray, atcoord,& 
&     natpermol, vectmatrix, coeff, megacom, epots)
     
      use overlapFINAL
 
      real(8), intent(in) :: distarray(:,:,:), atcoord(:,:) 
      real(8), intent(in) :: coeff(:,:), epots(:)
      real(8), allocatable :: coeffa(:), coeffb(:), deltaEad(:)
      integer, intent(in) :: natpermol, neighblist(:,:)
      real(8), intent(out) :: karray(:,:)
      real(8), intent(in) :: vectmatrix(:,:), megacom(:,:)
      integer :: n, i, z, r, k, j, ioerr
      integer :: ierr, blom, degen
      integer :: theindex, theindex2
      real(8) :: neighbourcom(3)
      real(8), allocatable :: sabmat(:,:)
      real(8), allocatable :: neighbourmol(:,:), current(:,:)
      real(8) :: sab, hab, rate, check, deltaE, Ebarrier

      !do i=1,size(distarray,2)

      !    write(*,*) distarray(1,i,2), distarray(2,i,2), distarray(3,i,2)

      !end do

      n=size(neighblist, 2)     

      degen=size(coeff,2)
 
      allocate(neighbourmol(3,natpermol), stat=ierr)

      if (ierr .ne. 0) then
          write(*,*) 'neighbourmol is baaaad'
          stop
      end if     

      allocate(current(3,natpermol), stat=ierr)
      if (ierr .ne. 0) then
          write(*,*) 'current is baaaad'
          stop
      end if
      
      allocate(coeffa(natpermol), stat=ierr)
      if (ierr .ne. 0) then
          write(*,*) 'memory problem with coeffa'
          stop
      end if

      allocate(coeffb(natpermol), stat=ierr)
      if (ierr .ne. 0) then
          write(*,*) 'current is baaaad'
          stop
      end if

      allocate(sabmat(degen,degen), stat=ierr)
      if (ierr .ne. 0) then
          write(*,*) 'memory problem with sabmat'
          stop
      end if

      allocate(deltaEad((size(neighblist, 1)-1)), stat=ierr)

      if (ierr .ne. 0) then
          write(*,*) 'neighbourmol is baaaad'
          stop
      end if


      do z= 1,n 
 
!     Looking up the donor molecule on which we are performing the rate caclulation

         do r=1,natpermol

            theindex = ((z-1)*natpermol+r)
            current(:,r)=atcoord(:,theindex) !*1.889725989

         end do
 
         !if(z .eq. 1) then

         !    do r=1,natpermol

         !    write(*,*) current(1,r), current(2,r), current(3,r)

         !    end do

         !end if

      
         current=current*1.889725989D0

         do i=2,neighblist(1,z)+1


!        Performing translational tranformation on the full molecule to match its real position 
!        (aka creating acceptor coordinates) 

            do r=1,natpermol

               if(mod(neighblist(i,z), n) .eq. 0) then

                   blom=n

               else 

                   blom = mod(neighblist(i,z), n)

               end if
           
               theindex2 = ((blom-1)*natpermol+r)

               neighbourmol(:,r)=atcoord(:,theindex2)&
&              +(megacom(:, neighblist(i,z)) - megacom(:,blom))
     
            end do

            !if(z .eq. 1 .and. (i-1) .eq. 1) then

            !   do r=1,natpermol

            !      write(*,*) current(1,r), current(2,r), current(3,r)

            !   end do

            ! end if



!     Calculating the overlap for the dimers

        !if(z .eq. 1 .and. (i-1) .eq. 1) then

            !do r=1,natpermol

            !write(*,*) neighbourmol(1,r), neighbourmol(2,r), neighbourmol(3,r)

            !end do

        !end if

        neighbourmol=neighbourmol*1.889725989D0

        deltaE=(epots(blom)-epots(z))*0.0016D0 ! now in Hartree

        do k=1,degen
           
           coeffa=coeff(:,k)           

           !if(k .eq. 1) write(*,*) coeffa

           do j=1,degen

              coeffb=coeff(:,j)
             
              !if(j .eq. 1) write(*,*) coeffb

              !call calc_sab(current, coeffa, sab, neighbourmol, coeffb)

              sabmat(j,k)=sab

           end do
        
        end do

        sab=0

        do k=1,degen

           do j=1,degen

              sab = sab+sabmat(j,k)**2

           end do

        end do
       
        !if(z .eq. 1 .and. i .eq. 2) then

        !do k=1,3

        !        write(*,*) (sabmat(k,j), j=1,3)

        !end do

        !end if


        !sab=abs(sabmat(2,3))

        hab=sqrt(sab)*(-1819.0D0)/27211.396D0/dble(degen)

        !sab=sabmat(1,1)

        
        !hab = sabmat(1,1)*(-1819.0D0)/27211.396D0


        neighbourcom=distarray(:,i-1,z)

        !if(z .eq. 36) then

        !write(*,*) hab*27211.396, neighbourcom(1), neighbourcom(2), neighbourcom(3)

        !end if



!     Caclualting the hopping rate
        !call hopping_rate(neighbourcom, hab, rate, check, deltaE, Ebarrier)
  
        !if(z .eq. 99) then

        !write(43,*) rate

        !end if

        deltaEad(i-1)=Ebarrier

        karray(i-1,z)=rate

       end do

       !write(43,*) (deltaEad(k), k=1,(size(neighblist,1)-1))
      end do

       !write(43,*) 'blab' !'end of snapshot'

      end subroutine


      subroutine ratelist_creator2(distarray, neighblist, karray, atcoord,&
&                                  natpermol, vectmatrix, coeff, megacom, &
&                                  epots, connectlist, lambdai, omega,    &
&                                  Ein, T, outer, epss, epsop, rVdW)

!     New version of rate list creator uses atomlists and the new version
!     of the overlap calculator.


      use overlapFINAL

      real(8), intent(in) :: distarray(:,:,:), atcoord(:,:), Ein(:)
      integer, intent(in) :: connectlist(:,:)
      real(8), intent(in) :: coeff(:,:), epots(:)
      real(8), allocatable :: coeffa(:,:), coeffb(:,:), deltaEad(:)
      integer, intent(in) :: natpermol, neighblist(:,:)
      real(8), intent(out) :: karray(:,:)
      real(8), intent(in) :: vectmatrix(:,:), megacom(:,:)
      integer :: n, i, z, r, k, j, ioerr
      integer :: ierr, blom, degen
      integer :: theindex, theindex2
      real(8) :: neighbourcom(3)
      real(8), allocatable :: sabmat(:,:)
      real(8), allocatable :: neighbourmol(:,:), current(:,:)
      real(8) :: sab, hab, rate, check, deltaE, Ebarrier
      real(8) :: lambdai, omega, T, epss, epsop, rVdW
      logical, intent(in) :: outer

      n=size(neighblist, 2)

      degen=size(coeff,2)

      allocate(neighbourmol(4,natpermol), stat=ierr)

      if (ierr .ne. 0) then
          write(*,*) 'neighbourmol is baaaad'
          stop
      end if

      allocate(current(4,natpermol), stat=ierr)
      if (ierr .ne. 0) then
          write(*,*) 'current is baaaad'
          stop
      end if

      allocate(coeffa(natpermol,degen), stat=ierr)
      if (ierr .ne. 0) then
          write(*,*) 'memory problem with coeffa'
          stop
      end if

      allocate(coeffb(natpermol, degen), stat=ierr)
      if (ierr .ne. 0) then
          write(*,*) 'current is baaaad'
          stop
      end if

      allocate(sabmat(degen,degen), stat=ierr)
      if (ierr .ne. 0) then
          write(*,*) 'memory problem with sabmat'
          stop
      end if

      allocate(deltaEad((size(neighblist, 1)-1)), stat=ierr)

      if (ierr .ne. 0) then
          write(*,*) 'neighbourmol is baaaad'
          stop
      end if


      do z= 1,n


!     Looking up the donor molecule on which we are performing the rate caclulation

         do r=1,natpermol

            theindex = ((z-1)*natpermol+r)
            current(:,r)=atcoord(:,theindex) 

         end do

         current(2:4,:)=current(2:4,:)*1.889725989D0


         coeffa=coeff

         ! Normalising non-idealised structure coeffs

         do k=1,degen

            call calc_sab(current, connectlist, coeffa(:,k))

         end do


         do i=2,neighblist(1,z)+1


!        Performing translational tranformation on the full molecule to match its real position 
!        (aka creating acceptor coordinates) 

            do r=1,natpermol

               if(mod(neighblist(i,z), n) .eq. 0) then

                   blom=n

                else

                   blom = mod(neighblist(i,z), n)

                end if

                theindex2 = ((blom-1)*natpermol+r)

                neighbourmol(:,r)=atcoord(:,theindex2)
                neighbourmol(2:4,r)= neighbourmol(2:4,r)&
&               +(megacom(:, neighblist(i,z)) - megacom(:,blom))

             end do

!            Calculating the overlap for the dimers

             neighbourmol(2:4,:)=neighbourmol(2:4,:)*1.889725989D0

             deltaE=(epots(blom)-epots(z))*0.0016D0 ! now in Hartree

             coeffb=coeff
          
             do k=1,degen

                 call calc_sab(neighbourmol, connectlist, coeffb(:,k))

             end do


             do k=1,degen


                do j=1,degen

                   call calc_sab(current, connectlist, coeffa(:,k), sab, &
&                           neighbourmol, connectlist,coeffb(:,j))

                   sabmat(j,k)=sab

                end do
             end do

             sab=0

!            Calculating RMS Hab a'la Newton

             do k=1,degen

                do j=1,degen

                   sab = sab+sabmat(j,k)**2

                end do

             end do


             hab=sqrt(sab)*(1819.0D0)/27211.396132D0/dble(degen)

!            Enable following lines for non-degenerate calculations

             !sab=sabmat(1,1)

             !hab = sabmat(1,1)*(-1819.0D0)/27211.396D0


             neighbourcom=distarray(:,i-1,z)


!            Caclualting the hopping rate

             call hopping_rate(neighbourcom, hab, rate, check, deltaE, & 
&                  Ebarrier, lambdai, omega, Ein, T, outer, epss, epsop,&
&                  rVdW)

             write(43,*) rate, hab*27211.3961320D0, neighbourcom,deltaE*27211.3961320D0, Ebarrier, z, neighblist(i,z)

             deltaEad(i-1)=Ebarrier

             karray(i-1,z)=rate

          end do

      end do

      end subroutine



      subroutine hopping_rate(neighbourcom, hab, rate, check, deltaEin, Ebarrier,&
                              lambdain, omega, Ein, T, outer, epss, epsop, rVdW)


      real(8), intent(in) :: neighbourcom(:),  deltaEin, Ein(:), T
      integer :: i, ierr
      real(8), intent(in) :: hab, omega, lambdain, epss, epsop, rVdW
      real(8), intent(out) :: rate, check, Ebarrier
      real(8) :: lambdai, kb, qe, E(3), omega_nau, DDOT, DNRM2, pi, omega_n
      real(8) :: lambda, twopigamma, PLZ, kappa_el, Ebar
      real(8) :: neighbcomloc(3), deltaEddagger, Delta, deltaE
      logical, intent(in) :: outer  

!     Hopping rate calculation algorithm using the formula from PCCP 2012
    

      lambdai=lambdain/27211.396132D0 !from meV to Hartree

      kb = 8.617343D-5/27.211396132D0

      qe = -1.0D0

      !T = 300.0D0

      E=Ein

      E = (E/27.211396132D0)/(1.889725989D8)

      deltaE=deltaEin !*3.0D0

      pi=4.D0*datan(1.D0)

      neighbcomloc = neighbourcom(:)

      lambda=lambdai

      if (outer .eqv. .TRUE.) then

      lambda=lambdai+((1.0D0/epsop)-(1.0D0/epss)) *(1.0D0/(rVdW)&
&            -1.0D0/(DNRM2(3,neighbcomloc(:),1))) * 1.439964521677473e+04      &
&            /27211.396132D0

      end if

      neighbcomloc = neighbourcom(:)*1.889725989D0

      omega_nau=omega/27211.396132D0 !61.5
      omega_n=omega/(1000.D0*6.5821189D-16)

      twopigamma = (pi**(1.5D0))*(hab**2)/(omega_nau*dsqrt(lambda*kb*T))

      PLZ=1-dexp(-twopigamma)

      kappa_el=(2.D0*PLZ)/(PLZ+1)

!     With potential energy .neq. 0

       deltaEddagger=(lambda+deltaE+qe*DDOT(3, neighbcomloc(:),1, E, 1))**2/(4.0D0*lambda)

!     With potential energy = 0

      !deltaEddagger=(lambda+qe*DDOT(3, neighbcomloc(:),1, E, 1))**2/(4.0D0*lambda)

!     With potential energy .neq. 0

      Delta=dabs(hab)+(lambda+deltaE+(qe*DDOT(3, neighbcomloc(:),1, E, 1)))/2.0D0 &
&     -dsqrt(hab**2 +(lambda+deltaE+qe*DDOT(3, neighbcomloc(:),1, E, 1))**2.D0/4.0D0)

!     With potential energy = 0

      !Delta=abs(hab)+(lambda+(qe*DDOT(3, neighbcomloc(:),1, E, 1)))/2.0D0 &
      !&   -sqrt(hab**2 +(lambda+qe*DDOT(3, neighbcomloc(:),1, E, 1))**2.D0/4.0D0)


!      Ebarrier = deltaE*27211.396132D0
      Ebarrier = (deltaEddagger-Delta)*27211.396132D0

!     Setting the energy barrier 0 where it would be negative

      Ebar=((deltaEddagger-Delta)+abs(deltaEddagger-Delta))/2.0D0

!      Ebar=(deltaEddagger-Delta)


      rate=omega_n/(2.D0*pi)*kappa_el*dexp(-1.0D0*(Ebar)/(kb*T))

      return

      end subroutine

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                                          !
      !           THESE SUBROUTINES ARE FOR MONTE CARLO PROPAGATION              !
      !                                                                          !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine hopper(natpermol, elposold, elposnew, &
&     tmig, distance, neighblist, distarray, karray)


      ! Propagates the electron in time and space registering and updating 
      ! the migration time, the displacement of the electron, and the position 
      
      implicit none
 
      integer,intent(in) :: neighblist(:,:)
      real*8,intent(in) :: distarray(:,:,:), karray(:,:)
      integer, intent(in) :: elposold, natpermol
      integer, intent(out) :: elposnew
      real(8), allocatable :: neighbourcom(:,:), neighbours(:,:)
      real(8), allocatable :: rates(:)
      integer :: i, nmol,  r, j
      integer, allocatable :: ind(:)
      integer :: ierr, neigh, colno 
      real(8) :: distvec(3), DNRM2, posvec(3), tjump, tmig, distance(3)
      logical :: indic = .TRUE.
      
      nmol=size(karray,2)

!     Creating periodic images in all the directions (first nearest neighbour
!     approach can be extended further)

      j = neighblist(1,elposold)

      allocate(ind(j), stat=ierr)
      if (ierr .ne. 0) then
         write(*,*) 'ind is baaaad'
         stop
      end if


      allocate(rates(j), stat=ierr)
      if (ierr .ne. 0) then
         write(*,*) 'rates is baaaad'
         stop
      end if

      rates = karray(1:j,elposold)

      ind = neighblist(2:j+1,elposold)

!     One Monte Carlo jump according to Voter

      call the_monte_carlo(rates, ind, colno, elposnew, tjump)

!     Collecting and summing up migration time and ditance

      distvec = distarray(:,colno,elposold)

      distance(:)=distance(:)+distvec(:)
       
      tmig =tmig + tjump

!     Applying periodic boundary condition: folding back the image cell to the
!      unit cell

      if (elposnew > nmol) then

      elposnew = mod(elposnew,nmol)

            if (elposnew .eq. 0) then

                 elposnew=nmol

            end if

      end if


      deallocate(ind)
      deallocate(rates)
      
      return

      end subroutine


      subroutine the_monte_carlo(rates, ind, colno, chosen, tjump)

      real*8, intent(in) :: rates(:)
      real*8, intent(out) :: tjump
      real*8 :: r1, r2, sum, incid 
      integer :: n, i
      integer, intent(out) :: chosen, colno
      integer, intent(in) :: ind(:)
      

      call random_number(r1)

!      r1 = rand() !call random_number(r1)

!     Since the random number is generated on the 0 =< r < 1 interval in order to
!     avoid log(0) I redfine the random numbers on the 0 < r =< 1 interval

      r1 = 1.0d0 - r1 

!     First random number defines how long the charge sat on a site

      tjump= (-1.0d0/sum(rates))*dlog(r1) 

      call random_number(r2)

!      r2 = rand() !call random_number(r2)

      r2 = 1.0d0 - r2 

!     Second random number decides the next site

      n = size(rates)

      i=0

      chosen = 0

      do while (chosen .eq. 0)

        i=i+1
        !incid = (sum(rates(1:i)))/(sum(rates))
        !write(*,*) "Boundaries are ", incid
        if (r2 .le. (sum(rates(1:i)))/(sum(rates))) then
            
             colno=i 
             chosen = ind(i)

        end if

      end do
    
      return

      end subroutine

      subroutine simple_gen(com,elposinit)

!     Generating initial position of the charge in the unit cell

      real*8, intent(in) :: com(:,:)
      integer, intent(out) :: elposinit
      real*8 :: r
      integer :: n


      n=size(com,2)

!      r=rand() !random_number(r)

      call random_number(r)
    
      elposinit = ceiling(n * r)

      end subroutine


!     Leftover from FCC calculation

      subroutine habs_fcc_test(neighbourcom, j, habs)

      implicit none

      real*8, intent(in) :: neighbourcom (:,:)
      real*8, intent(out) :: habs(:) 
      integer, intent(in) :: j
      integer :: i
      real*8 :: pothabs(6)

      !write(*,*) 'Small!!!', epsilon(0.5d0)
     
      do i=1,j

         if(abs(neighbourcom(2,i)) < 100.0*epsilon(0.5d0) .and. &
&           abs(neighbourcom(3,i)) < 100.0*epsilon(0.5d0))  then
        
            habs(i)=4.92/27211.396132

         else if(abs(neighbourcom(1,i)) < 100.0*epsilon(0.5d0) .and. &
&                abs(neighbourcom(3,i)) < 100.0*epsilon(0.5d0))  then


            habs(i)=2.044/27211.396132

         else if(abs(neighbourcom(1,i)) < 100.0*epsilon(0.5d0) .and. &
&                abs(neighbourcom(2,i)) < 100.0*epsilon(0.5d0))  then

            habs(i)=1.478/27211.396132

         else if(abs(neighbourcom(1,i)) > 100.0*epsilon(0.5d0) .and. &
&                abs(neighbourcom(2,i)) > 100.0*epsilon(0.5d0))  then

            habs(i)=0.597/27211.396132

         else if(abs(neighbourcom(1,i)) > 100.0*epsilon(0.5d0) .and. &
&                abs(neighbourcom(3,i)) > 100.0*epsilon(0.5d0))  then

            habs(i)=13.227/27211.396132

         else if(abs(neighbourcom(2,i)) > 100.0*epsilon(0.5d0) .and. &
&                abs(neighbourcom(3,i)) > 100.0*epsilon(0.5d0))  then

            habs(i)=1.089/27211.396132         

         else 
         
         write(*,*) 'Uwaga, uwaga what what what?'
         
         end if

      !habs = habs*(1/27211.396132)

      end do

      return

      end subroutine

      subroutine build_lattice(n,unitcell)

      integer, intent(in) :: n
      real*8, intent(out) :: unitcell(:,:)
      real*8 :: a=14.3d0
      integer :: i, edge=3, j, elemunic


      unitcell(:,1) = (/0.0d0, 0.0d0, 0.0d0/)
      unitcell(:,2) = (/0.0d0, 0.5d0*a, 0.5d0*a/)
      unitcell(:,3) = (/0.5d0*a, 0.0d0, 0.5d0*a/)
      unitcell(:,4) = (/0.5d0*a, 0.5d0*a, 0.0d0/)

      

      do i=1,(edge-1)
      
      elemunic=4

      do j=1,elemunic
      
      unitcell(:,(elemunic*i)+j) = unitcell(:,j)+(/i*a, 0.0d0, 0.0d0/)
            
      end do
      end do

      do i=1,(edge-1)

      do j=1,elemunic*edge

      unitcell(:,(elemunic*edge*i)+j) = unitcell(:,j)+(/0.0d0, i*a, 0.0d0/)

      end do
      end do
      
      do i=1,(edge-1)

      do j=1,elemunic*edge*edge

      unitcell(:,(elemunic*edge*edge*i)+j) = unitcell(:,j)+(/0.0d0, 0.0d0, i*a/)
      
      !write(*,*)
      end do
      end do

      return
      end subroutine

      subroutine init_random_seed()
      implicit none
      integer, allocatable :: seed(:)
      integer :: i, n, un, istat, dt(8), pid, t(2), s
      integer(8) :: count, tms
          
      call random_seed(size = n)
      allocate(seed(n))
      ! First try if the OS provides a random number generator
      open(un, file="/dev/urandom", access="stream", &
           form="unformatted", action="read", status="old", iostat=istat)
      if (istat == 0) then
          read(un) seed
          close(un)
          else
      ! Fallback to XOR:ing the current time and pid. The PID is
      ! useful in case one launches multiple instances of the same
      ! program in parallel.
          call system_clock(count)
          if (count /= 0) then
          t = transfer(count, t)
          else
          call date_and_time(values=dt)
          tms = (dt(1) - 1970) * 365_8 * 24 * 60 * 60 * 1000 &
&                 + dt(2) * 31_8 * 24 * 60 * 60 * 1000 &
&                 + dt(3) * 24 * 60 * 60 * 60 * 1000 &
&                       + dt(5) * 60 * 60 * 1000 &
&                       + dt(6) * 60 * 1000 + dt(7) * 1000 &
&                       + dt(8)
                  t = transfer(tms, t)
               end if
               s = ieor(t(1), t(2))
               pid = getpid() + 1099279 ! Add a prime
               s = ieor(s, pid)
               if (n >= 3) then
                  seed(1) = t(1) + 36269
                  seed(2) = t(2) + 72551
                  seed(3) = pid
                  if (n > 3) then
                     seed(4:) = s + 37 * (/ (i, i = 0, n - 4) /)
                  end if
               else
                  seed = s + 37 * (/ (i, i = 0, n - 1 ) /)
               end if
            end if
            call random_seed(put=seed)
       end subroutine init_random_seed

!     LEftover from FCC calculation

      subroutine build_crystal(mol, crystal)

      real*8, intent(in) :: mol(:,:)
      real*8 :: unitcell(3,4), com(3,1)
      real*8, allocatable :: molloc(:,:), mass(:)
      real*8, intent(out) :: crystal(:,:)
      real*8 :: a=14.30d0
      integer :: i, edge=3, j, elemunic, k, natpermol, ierr

      unitcell(:,2) = (/0.0d0, 0.5d0*a, 0.5d0*a/)
      unitcell(:,3) = (/0.5d0*a, 0.0d0, 0.5d0*a/)
      unitcell(:,4) = (/0.5d0*a, 0.5d0*a, 0.0d0/)

      natpermol=size(mol,2)

      allocate(mass(natpermol),stat=ierr)
      if (ierr .ne. 0) then
          write(*,*) 'Memory issue'
          stop
      end if

      mass=12

      allocate(molloc(3,natpermol),stat=ierr)
      if (ierr .ne. 0) then
          write(*,*) 'Memory issue'
          stop
      end if

      call centre_of_mass_old(mol, natpermol, mass, com)
    
      do i=1,natpermol

         molloc(:,i)=mol(:,i)-com(:,1)
     
      end do

      do i=0,(edge-1)

      elemunic=4

      do j=1,elemunic

      do k=1,natpermol

      crystal(:, (i*elemunic+(j-1))*natpermol+k) = molloc(:,k)+unitcell(:,j)+(/i*a, 0.0d0, 0.0d0/)

      end do
      end do
      end do

      do i=1,(edge-1)

         do j=1,elemunic*edge

            do k=1,natpermol

               crystal(:,((elemunic*edge*i)+(j-1))*natpermol+k) = & 
&              crystal(:,(k+(j-1)*natpermol))+(/0.0d0, i*a, 0.0d0/)

            end do
         end do
      end do

      do i=1,(edge-1)

         do j=1,elemunic*edge*edge

            do k=1,natpermol

               crystal(:,((elemunic*edge*edge*i)+(j-1))*natpermol+k) =& 
&             crystal(:,(k+(j-1)*natpermol))+(/0.0d0, 0.0d0, i*a/)

            end do
         end do
      end do

      return
      end subroutine

      subroutine full_KMC(atcoord, nat, mass, nmol, vectmatrix,&
&                         coeff, epots, connectlist, atomlist, &
&                         lambdai, omega, tmax, nKMC, cutoff, Ein, &
&                         T, outer, epss, epsop, rVdW)

      real(8), dimension(:,:), intent(in) :: atcoord
      integer,                 intent(in) :: nat, nmol
      real(8), dimension(:),   intent(in) :: mass, Ein
      real(8), dimension(:,:), intent(in) :: vectmatrix
      real(8), dimension(:),   intent(in) :: epots
      real(8), dimension(:,:), intent(in) :: coeff
      integer, dimension(:,:), intent(in) :: connectlist
      integer, dimension(:),   intent(in) :: atomlist

      real(8), allocatable, dimension(:,:) :: com, megacom
      integer, allocatable, dimension(:,:) :: neighblistpre, neighblist
      real(8), allocatable, dimension(:,:,:) :: distarraypre, distarray 
      real(8), allocatable, dimension(:,:) :: karray
   
      integer :: ierr, jmax, elposinit, elposnew, elposold, rom
      real(8) :: tmig, distance(3), start, finish, timing 
      integer :: nKMC, dump
      real(8) :: lambdai, omega, tmax, cutoff, T
      real(8) :: epsop, epss, rVdW  
      logical, intent(in) :: outer


      allocate(com(3,nmol),stat=ierr)
      if (ierr .ne. 0) then
         write(*,*) 'Memory problem at com'
      stop
      end if

      allocate(megacom(3,nmol*27), stat=ierr)
      if (ierr .ne. 0) then
         write(*,*) 'Memory problem at megacom'
      stop
      end if


!     Centre of mass calculation and periodic image building

      call centre_of_mass(atcoord(2:4,:), atomlist, nat, mass, com)

      call com_pbc(com, megacom, vectmatrix)

      allocate(neighblistpre(100,nmol), stat=ierr)
      if (ierr .ne. 0) then
          write(*,*) 'neighblistpre is baaaad'
          stop
      end if

      neighblistpre = 0

      allocate(distarraypre(3,100,nmol), stat=ierr)
      if (ierr .ne. 0) then

          write(*,*) 'distarraypre is baaaad'
          stop

      end if

      distarraypre = 0

      jmax=0

!     Calculating neighbourlist for given snapshot

      call static_neighblist(com,megacom, neighblistpre, &
&                            jmax, distarraypre, cutoff)

      allocate(neighblist(jmax+1,nmol), stat=ierr)
      if (ierr .ne. 0) then
          write(*,*) 'neighblist is baaaad'
          stop
      end if

      neighblist = 0

      allocate(distarray(3,jmax,nmol), stat=ierr)
      if (ierr .ne. 0) then
          write(*,*) 'distarray is baaaad'
          stop
      end if

      distarray = 0

      allocate(karray(jmax,nmol), stat=ierr)
      if (ierr .ne. 0) then
          write(*,*) 'karray is baaaad'
          stop
      end if

      karray = 0

      neighblist(:,:)= neighblistpre(1:(jmax+1),:)
      distarray(:,:,:)=distarraypre(:,1:jmax,:)

      deallocate(neighblistpre, distarraypre)


!     Using static molecular neighbour list calculating the static ratelist 

      call ratelist_creator2(distarray, neighblist, karray, atcoord,&
&                   nat, vectmatrix, coeff, megacom, epots, connectlist,&
&                   lambdai, omega, Ein, T, outer, epss, epsop, rVdW)

      call init_random_seed()

      call cpu_time(start)

!    Do nKMC trajectories using different initial positions and the
!    static rate list

      do rom=1,nKMC

         call simple_gen(com, elposinit)

         elposold=elposinit

         distance = 0

         tmig = 0

!        Running simulation until system time reaches maximum: tmax

         dump=1

         do while (abs(tmig) < 9.0D0*tmax)


              call hopper(nat, elposold, elposnew, tmig, distance, &
&                        neighblist, distarray, karray)

              elposold = elposnew

              !if (tmig > dump*tmax) then

              !       write(44+dump,*)  distance
    
              !      dump=dump+1

              !end if

        end do

!       Dumping simulation time and final position into output files 

        write(44,*),  tmig
        write(54,*),  distance

      end do

!     FYI CPU time about calculation on STD(OUT)

      call cpu_time(finish)

      timing=finish-start

      write(*,*) 'Here is the timing', timing

      return     

      end subroutine full_KMC


      end module kinetic_monte_carlo

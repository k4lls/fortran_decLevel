program declevel

implicit none

integer                     :: nfft                  ! data window length
integer                     :: ndfc                  ! decimation filter length
integer                     :: mdl                   ! maximum number decimation levels
integer                     :: nhfreq                ! # harmonic frequencies per data window 
integer                     :: noverlap              ! data window overlap, 75% for 4PiProlate & Kaiser tapers  
integer                     :: noffset               ! data window offset
real                        :: tsfrq                 ! time series length
integer                     :: npnt                  ! Number of points
integer                     :: nfc                   ! fourier coefficient estimates
real                     :: fileSize     
integer                     :: i , ndl , nhfrq , nw , sec
integer,dimension(:,:),allocatable    :: SUMMARY1
real,dimension(:,:),allocatable       :: SUMMARY2
integer                               :: s,m,h,d


nfft=64                     ! data window length
ndfc=9                      ! decimation filter length
mdl=25                      ! maximum number decimation levels
nhfreq=3                    ! # harmonic frequencies per data window 
noverlap=3*nfft/4           ! data window overlap, 75% for 4PiProlate & Kaiser tapers  
noffset=nfft-noverlap       ! data window offset

allocate(SUMMARY1(mdl,3))
allocate(SUMMARY2(mdl,2))

! time series length: 256 hertz AD * 10 hours * min/hour * sec/min
write (*,*), '---------------------------------------------------------------------'
write (*,*), 'Select S/R (hertz) [ex:256,1024,4096]'
read *, tsfrq               ! time series sampling frequency (hertz=1/sec)
write (*,*),''
write (*,*), '---------------------------------------------------------------------'
write (*,*), 'Select Timeseries duration (sec) [ex: 1hour=3600s]'
read *, sec


d=sec/86400
h=mod(sec,86400)/3600
m=mod(sec,3600)/60
s=mod(sec,60)


fileSize=((tsfrq*4+32)*sec+512)/(1000**2)
npnt  = tsfrq*sec

   nfc = 0;        ! initialize # fourier coefficient estimates
   do i = 1,mdl
     if (npnt<nfft) then
       exit
     end if

     ndl = i
     nw = 1 + ((npnt-nfft)/noffset) ! # data windows
     nfc = nfc + nw*nhfrq
     SUMMARY1(i,1) = ndl             ! Decimation level
     SUMMARY1(i,2) = nw              ! one estimate per window at each frequency
     SUMMARY1(i,3) = npnt
     SUMMARY2(i,1) = 6*tsfrq/32      ! using same harmonics from 32 & 64 pnt windows
     SUMMARY2(i,2) = 8*tsfrq/32
     npnt = ((npnt-ndfc)/2) + 1
     tsfrq = tsfrq/2;
     if (SUMMARY2(i,1)<0.0007) then
     ndl=ndl-1
       exit
     end if
   end do
     
   ! DISPLAY
   100 format(f10.4,f10.4,12x,i8,12x,i8)
   101 format(A,A,A)
   99 format(A,i4,A,i2.2,A,i2.2,A,i2.2,A)

   write (*,*),''
   write  (*,*),'---------------------------------------------------------------------'
   write(*,fmt=101), '|     Frequencies     | ' , ' Number of estimation(s)  | ' , 'Decimation level  | '
   write (*,*),'---------------------------------------------------------------------'
   write (*,*),''
   do i=1,ndl
   write(*,fmt=100) SUMMARY2(i,2), SUMMARY2(i,1), SUMMARY1(i,2), SUMMARY1(i,1)
   end do
   write (*,*), '---------------------------------------------------------------------'
   write (*,'(A,f10.2,A)'), 'Z3D file size : ' , fileSize , ' MB'
   write (*,fmt=99), 'Duration      : ',d,' day(s) ',h,'h',m,'m',s,'s'
   write (*,*), '---------------------------------------------------------------------'
   
deallocate(SUMMARY1)
deallocate(SUMMARY2)

end program declevel
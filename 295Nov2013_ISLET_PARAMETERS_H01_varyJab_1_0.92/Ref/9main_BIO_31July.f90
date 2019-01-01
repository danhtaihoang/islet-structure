!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!     HOANG Danh Tai - Asia Pacific Center for Theoretical Physics, 
!     Hogil Kim Memorial Building #501 POSTECH,
!     San 31, Hyoja-dong, Namgum, Pohang, Gyeongbuk 790-784, Korea.
!     Personal site: http://hoangdanhtai.com
!-----------------------------------------------------------------------------------------!
!!!   Estimate parameter Jaa, Jab, Jbb from data
!!!   16.6: Bo sung T, de chon gia tri E be nhat
!!!!  17.6: Check xem neu dung ic la nearest neighbor thi co dung khong???
!!!!  30.07: Sua lai doc experiment data tu thu muc
!!    31.07.2013: Sua lai r(i) de tim nearest neighbor
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%! 
      PROGRAM main_bio
      IMPLICIT NONE

      CHARACTER (LEN=3)  :: SAt
      CHARACTER (LEN=150):: name_data,CONFIG_INI
      CHARACTER (LEN=50) :: name
      CHARACTER (LEN=15) :: tmp
      CHARACTER (256)    :: Ligne23
      REAL    (KIND=8),PARAMETER :: nul=0.

      INTEGER (KIND=8):: i,j,n_cell,n_max,i_n,ic,S,i_Jbb,n_Jbb,i_Jab,n_Jab,number
      INTEGER (KIND=8):: i_delN,i_loop,i_loop1,i_loop2,n_equi1,n_equi2,n_average,i_times,n_times
      INTEGER (KIND=8):: nn_total,nn_max,na,nb,n_T,i_T,ia,ib,j_ic,n_min,nn_min
     
      REAL    (KIND=8):: r02,r2,n_cell_p,delT,T,Tmin,Tmax,rdn_config
      
      REAL    (KIND=8):: energy,E1_old,E2_old,E1_new,E2_new,Naa,Nbb,Nab,H1a,H1b,H2a,H2b
      REAL    (KIND=8):: Jaa,Jbb,Jab,del_Jbb,Jbbmax,Jbbmin,del_Jab,Jabmax,Jabmin
      REAL    (KIND=8):: Sa_mtp,Sb_mtp,rdn_ic,rdn_mtp,Naa_av,Nbb_av,Nab_av
      REAL    (KIND=8):: E_av,Nabc_total,delNabc,Naa0,Nbb0,Nab0,Nabc0_total
     
      INTEGER (KIND=8),DIMENSION(:),ALLOCATABLE :: nn
      INTEGER (KIND=8),DIMENSION(:,:),ALLOCATABLE :: name_in
      REAL (KIND=8),DIMENSION(:),ALLOCATABLE :: Sa,Sb,x,y,z,r
      REAL (KIND=8),DIMENSION(:,:),ALLOCATABLE :: delNabc_tab                    
!!!=======================================================================================
!!!=======================================================================================
      CALL system('rm -r config_ini_3D')
      CALL system('mkdir config_ini_3D')
      CALL system('rm -r config_3D')
      CALL system('mkdir config_3D')
      CALL system('rm *.dat*')
      
      CALL ini_rdm_number()
      CALL read_input_parameter()
      CALL open_data()
          
      ALLOCATE(x(n_cell),y(n_cell),z(n_cell),r(n_cell))
      ALLOCATE(Sa(n_cell),Sb(n_cell))
      ALLOCATE(nn(n_cell))
      ALLOCATE(delNabc_tab(n_Jbb,n_Jab))
      
      CALL read_data()
      
      IF (n_T==1) THEN
            delT=0.
      ELSE
            delT=(Tmax-Tmin)/real(n_T-1)
      END IF

      IF (n_Jbb==1) THEN
            del_Jbb=0.
      ELSE
            del_Jbb=(Jbbmax-Jbbmin)/real(n_Jbb-1)
      END IF
  
      IF (n_Jab==1) THEN
            del_Jab=0.
      ELSE
            del_Jab=(Jabmax-Jabmin)/real(n_Jab-1)
      END IF

      n_cell_p=real(n_cell)+1.

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! ====== MAIN PROGRAM === MAIN PROGRAM === MAIN PROGRAM === MAIN PROGRAM ======
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
      
      CALL write_config_ini_3D()
      CALL number_neighbor_ini()

      OPEN(unit=21,file='average_thermal.dat')
      OPEN(unit=22,file='del_Nabc2_min.dat')
      OPEN(unit=23,file='energy.dat')        
      
      i_delN=0

      DO i_T=1,n_T
      WRITE(*,*)'i_T = ', i_T                       
      T=Tmin+delT*real(i_T-1)

      DO i_Jbb=1,n_Jbb
            WRITE(*,*)'i_Jbb=',i_Jbb
            Jbb=Jbbmin+del_Jbb*real(i_Jbb-1)
                  
            DO i_Jab=1,n_Jab

            i_delN=i_delN+1
            !WRITE(*,*)'i_delN=',i_delN
            Jab=Jabmin+del_Jab*real(i_Jab-1)

            IF (CONFIG_INI=='NO') THEN
                  CALL load_config_ini()
            END IF
            
            CALL average_thermal()
       
            !CALL value_thermal()
                  
            !WRITE(*,*)'Naa=',Naa,'Nbb=',Nbb,'Nab=',Nab
            
            CALL write_config_3D()

            delNabc_tab(i_Jbb,i_Jab)=delNabc

            END DO
      END DO
      END DO
      
      WRITE(22,'(2I5,5X,F16.3)')MINLOC(delNabc_tab),MINVAL(delNabc_tab)
               
      WRITE(*,*)Jbb,MINLOC(delNabc_tab),MINVAL(delNabc_tab)

      CLOSE(21)
      CLOSE(22)
      CLOSE(23)
            
      CONTAINS

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! SUBROUTINE init_rdm_number()
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE ini_rdm_number()
      IMPLICIT NONE

      INTEGER (KIND=8),DIMENSION(8) :: time
      INTEGER (KIND=8),DIMENSION(50) :: seed

      CALL DATE_AND_TIME(values=time)     ! Get the current time
      seed(1) = time(4)*(360000*time(5) + 6000*time(6) + 100*time(7) + time(8))
      CALL RANDOM_SEED(PUT=seed)

      END SUBROUTINE ini_rdm_number 
      
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! OPEN the parameter from file "parameter.in"
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE read_input_parameter()
      IMPLICIT NONE
      CHARACTER (LEN=150) :: tamp
      OPEN(11,file='1parameter.in')
 
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A30,(A10))')   tamp, CONFIG_INI
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A30,(A10))')   tamp, name_data
      READ(11, '(A30,(F12.6))') tamp, r02
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A30,(F12.6))') tamp, Jaa
      READ(11, '(A50)')         tamp
      READ(11, '(A30,(F12.6))') tamp, Jbbmin
      READ(11, '(A30,(F12.6))') tamp, Jbbmax
      READ(11, '(A30,(I5))')    tamp, n_Jbb
      READ(11, '(A50)')         tamp
      READ(11, '(A30,(F12.6))') tamp, Jabmin
      READ(11, '(A30,(F12.6))') tamp, Jabmax
      READ(11, '(A30,(I5))')    tamp, n_Jab
      READ(11, '(A50)')         tamp
      READ(11, '(A30,(F12.6))') tamp, Tmin
      READ(11, '(A30,(F12.6))') tamp, Tmax
      READ(11, '(A30,(I5))')    tamp, n_T
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A30,(I12))')   tamp,n_equi1
      READ(11, '(A30,(I12))')   tamp,n_equi2
      READ(11, '(A30,(I12))')   tamp,n_average
      READ(11, '(A30,(I12))')   tamp,n_times
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp

      CLOSE(11) 

      END SUBROUTINE read_input_parameter
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! Change config as random from data
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE load_config_ini()
      IMPLICIT NONE

      Sa(:)=0. ; Sb(:)=0.

      ia=0 ; ib=0
      
      DO WHILE (ia<na)
            i=0

            DO WHILE(i==0)
                  CALL random_number(rdn_config)
                  i=int(rdn_config*real(n_cell_p))
            ENDDO

            IF (int(Sa(i))==0) THEN                      
                  Sa(i)=1.
                  ia=ia+1
            END IF

      END DO

      DO i=1,n_cell

            IF (int(Sa(i))==0) THEN   
                  Sb(i)=1.
                  ib=ib+1
            END IF

      ENDDO

      !WRITE(*,*)'ia:',ia,'ib:',ib

      END SUBROUTINE load_config_ini
   
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! WRITE initial position configuration in 3D
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE write_config_ini_3D()
      IMPLICIT NONE

      OPEN(unit=13,file='config_ini_3D/config_ini_3D.pdb')
      
      DO i=1,n_cell
           
            IF (int(Sb(i))==1) THEN 
                  SAt='Cu'
                  WRITE(13,'("ATOM",9X,A3,14X,3F8.3,1X,"1.00",1X,F8.3)') &
                        SAt,x(i),y(i),z(i),nul
            ELSE
                              
            IF (int(Sa(i))==1) THEN 
                  SAt='Au'
                  WRITE(13,'("ATOM",9X,A3,14X,3F8.3,1X,"1.00",1X,F8.3)') &
                        SAt,x(i),y(i),z(i),nul
            ELSE
                  SAt='H'
                  WRITE(13,'("ATOM",9X,A3,14X,3F8.3,1X,"1.00",1X,F8.3)') &
                        SAt,x(i),y(i),z(i),nul

            END IF
            END IF

      END DO

      CLOSE(13)

      END SUBROUTINE write_config_ini_3D
      
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! WRITE initial position configuration in 3D
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE write_config_3D()
      IMPLICIT NONE

      number=10000000+i_delN
      
      WRITE(tmp,'(I8)') number

      name='_BIO_'//TRIM(tmp)
      
      OPEN(unit=14,file='config_3D/config_3D'//trim(name)//'.pdb')

      DO i=1,n_cell
           
            IF (int(Sb(i))==1) THEN 
                  SAt='Cu'
                  WRITE(14,'("ATOM",9X,A3,14X,3F8.3,1X,"1.00",1X,F8.3)') &
                        SAt,x(i),y(i),z(i),nul
            ELSE
                              
            IF (int(Sa(i))==1) THEN 
                  SAt='Au'
                  WRITE(14,'("ATOM",9X,A3,14X,3F8.3,1X,"1.00",1X,F8.3)') &
                        SAt,x(i),y(i),z(i),nul
            ELSE
                  SAt='H'
                  WRITE(14,'("ATOM",9X,A3,14X,3F8.3,1X,"1.00",1X,F8.3)') &
                        SAt,x(i),y(i),z(i),nul

            END IF
            END IF

      END DO

      CLOSE(14)

      END SUBROUTINE write_config_3D
 !!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!  SUBROUTINE Open data and Read_data() Read experimental data
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE open_data()
      IMPLICIT NONE
      
      !!----------- For TEST ------------------
      IF (name_data=='GS1') THEN     
            OPEN(unit=12,file='GS1-40-60.txt')
            n_cell=100
      END IF

      IF (name_data=='GS3') THEN     
            OPEN(unit=12,file='GS3-50-50.txt')
            n_cell=100
      END IF

      !!! ---------- For Human --------------
      IF (name_data=='H01') THEN     
            OPEN(unit=12,file='H01.txt')
            n_cell=681
      END IF   

      IF (name_data=='H02') THEN     
            OPEN(unit=12,file='H02.txt')
            n_cell=1178
      END IF  
      
      IF (name_data=='H03') THEN     
            OPEN(unit=12,file='H03.txt')
            n_cell=340
      END IF
      
      IF (name_data=='H04') THEN     
            OPEN(unit=12,file='H04.txt')
            n_cell=2472
      END IF
      
      IF (name_data=='H06') THEN     
            OPEN(unit=12,file='H06.txt')
            n_cell=916
      END IF
      
      !!! ---------- For MOUSE --------------
      IF (name_data=='M01') THEN     
            OPEN(unit=12,file='DATA_MOUSE/M01.txt')
            n_cell=893
      END IF 

      IF (name_data=='M02') THEN     
            OPEN(unit=12,file='DATA_MOUSE/M02.txt')
            n_cell=2598
      END IF 
      
      IF (name_data=='M03') THEN     
            OPEN(unit=12,file='DATA_MOUSE/M03.txt')
            n_cell=1870
      END IF
      
      IF (name_data=='M04') THEN     
            OPEN(unit=12,file='DATA_MOUSE/M04.txt')
            n_cell=1680
      END IF
      
      IF (name_data=='M05') THEN     
            OPEN(unit=12,file='DATA_MOUSE/M05.txt')
            n_cell=518
      END IF
      
      IF (name_data=='M06') THEN     
            OPEN(unit=12,file='DATA_MOUSE/M06.txt')
            n_cell=1895
      END IF 

      IF (name_data=='M07') THEN     
            OPEN(unit=12,file='DATA_MOUSE/M07.txt')
            n_cell=571
      END IF 

      IF (name_data=='M08') THEN     
            OPEN(unit=12,file='DATA_MOUSE/M08.txt')
            n_cell=1192
      END IF

      IF (name_data=='M09') THEN     
            OPEN(unit=12,file='DATA_MOUSE/M09.txt')
            n_cell=1153
      END IF 

      IF (name_data=='M10') THEN     
            OPEN(unit=12,file='DATA_MOUSE/M10.txt')
            n_cell=1063
      END IF 

      IF (name_data=='M11') THEN     
            OPEN(unit=12,file='DATA_MOUSE/M11.txt')
            n_cell=2143
      END IF

      IF (name_data=='M12') THEN     
            OPEN(unit=12,file='DATA_MOUSE/M12.txt')
            n_cell=1576
      END IF

      IF (name_data=='M13') THEN     
            OPEN(unit=12,file='DATA_MOUSE/M13.txt')
            n_cell=2471
      END IF

      IF (name_data=='M14') THEN     
            OPEN(unit=12,file='DATA_MOUSE/M14.txt')
            n_cell=1113
      END IF

      IF (name_data=='M15') THEN     
            OPEN(unit=12,file='DATA_MOUSE/M15.txt')
            n_cell=1230
      END IF

      IF (name_data=='M16') THEN     
            OPEN(unit=12,file='DATA_MOUSE/M16.txt')
            n_cell=3294
      END IF

      IF (name_data=='M17') THEN     
            OPEN(unit=12,file='DATA_MOUSE/M17.txt')
            n_cell=1645
      END IF

      IF (name_data=='M18') THEN     
            OPEN(unit=12,file='DATA_MOUSE/M18.txt')
            n_cell=4159
      END IF

      IF (name_data=='M19') THEN     
            OPEN(unit=12,file='DATA_MOUSE/M19.txt')
            n_cell=2248
      END IF

      IF (name_data=='M20') THEN     
            OPEN(unit=12,file='DATA_MOUSE/M20.txt')
            n_cell=1182
      END IF

      IF (name_data=='M21') THEN     
            OPEN(unit=12,file='DATA_MOUSE/M21.txt')
            n_cell=4010
      END IF

      IF (name_data=='M22') THEN     
            OPEN(unit=12,file='DATA_MOUSE/M22.txt')
            n_cell=3884
      END IF

      IF (name_data=='M23') THEN     
            OPEN(unit=12,file='DATA_MOUSE/M23.txt')
            n_cell=3956
      END IF

      IF (name_data=='M24') THEN     
            OPEN(unit=12,file='DATA_MOUSE/M24.txt')
            n_cell=2097
      END IF

      IF (name_data=='M25') THEN     
            OPEN(unit=12,file='DATA_MOUSE/M25.txt')
            n_cell=3209
      END IF

      IF (name_data=='M26') THEN     
            OPEN(unit=12,file='DATA_MOUSE/M26.txt')
            n_cell=9314
      END IF

      IF (name_data=='M27') THEN     
            OPEN(unit=12,file='DATA_MOUSE/M27.txt')
            n_cell=1736
      END IF

      IF (name_data=='M28') THEN     
            OPEN(unit=12,file='DATA_MOUSE/M28.txt')
            n_cell=1430
      END IF

      IF (name_data=='M29') THEN     
            OPEN(unit=12,file='DATA_MOUSE/M29.txt')
            n_cell=2224
      END IF

      IF (name_data=='M30') THEN     
            OPEN(unit=12,file='DATA_MOUSE/M30.txt')
            n_cell=1631
      END IF

      END SUBROUTINE open_data    
!!!=======================================================================
      SUBROUTINE read_data()
      IMPLICIT NONE
  
      Sa(:)=0. ; Sb(:)=0.
      na=0 ; nb=0
      
      DO i=1,n_cell
            READ(12,*)S,x(i),y(i),z(i)
            
            IF (S==11) THEN
                  Sa(i)=0. ; Sb(i)=1.
                  
                  ELSE
                  
                  IF (S==12) THEN
                  Sa(i)=1. ; Sb(i)=0.
                  
                  END IF
                  
            END IF
          na=na+int(Sa(i)) ; nb=nb+int(Sb(i))
                  
      END DO
      
      CLOSE(12)     
    
      END SUBROUTINE read_data           
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! Number nearest neighbor from experimental data
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE number_neighbor_ini()
      IMPLICIT NONE
 
      OPEN(unit=15,file='nn_i_cell.dat')
      OPEN(unit=16,file='Nab0_ini.dat')
  
!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      nn_total=0 ;  nn(:)=0 ; r(:)=0.
        
      DO i=1,n_cell

            !DO WHILE ((nn(i)<12).and.(r(i)<(r02+50.)))
            
                  DO i_loop=1,10
                  r(i)=r02+i_loop
            
                        DO j=1,n_cell
                       
                        r2=(x(i)-x(j))**2.+(y(i)-y(j))**2.+(z(i)-z(j))**2.
                       
                        IF (((0.<r2).and.(r2<r(i))).and.(nn(i)<12)) THEN
                              nn(i)=nn(i)+1

                        END IF

                        END DO
                        
                  END DO      
            
            !END DO
         
            nn_total=nn_total+nn(i)   
            
            WRITE(15,*)i,nn(i),r(i)

      END DO
  
  
      
!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@       
      !!!-------------------------------------------------
      !!! maximum of number nearest neighbor
      nn_max=MAXVAL(nn) ; nn_min=MINVAL(nn)
      n_max=0 ; n_min=0
      DO i=1,n_cell
      
            IF (nn(i)==nn_max) THEN                  
                  n_max=n_max+1            
            END IF
            
            IF (nn(i)==nn_min) THEN                  
                  n_min=n_min+1            
            END IF
      
      END DO   
      
      !!!--------------------------------------------------------
      !!! Assign name for nearest neighbor of cell i
      
      ALLOCATE(name_in(n_cell,0:n_cell))
      
      name_in(:,:)=0
      
      DO i=1,n_cell
           i_n=0
           
            DO j=1,n_cell
           
                  r2=(x(i)-x(j))**2.+(y(i)-y(j))**2.+(z(i)-z(j))**2.
                 
                  IF ((0.<r2).and.(r2<r02)) THEN
                        i_n=i_n+1
                        name_in(i,i_n)=j         

                  END IF
                    
            END DO

      END DO    
      

      !!!-------------------------------------------------------
      !!! Calculate initial number interaction Naa0, Nbb0, Nab0
      Naa0=0. ; Nbb0=0. ; Nab0=0.
      DO i=1,n_cell
            CALL value_H1()
            
            Naa0=Naa0+Sa(i)*H1a ; Nbb0=Nbb0+Sb(i)*H1b ; Nab0=Nab0+Sa(i)*H1b+Sb(i)*H1a
      
      END DO 
      Naa0=Naa0/2. ; Nbb0=Nbb0/2. ; Nab0=Nab0/2.
      Nabc0_total=Naa0+Nbb0+Nab0          
      
      WRITE(16,'(A6,I8,6X,A6,I8,6X,A6,I8)')'na:',na,'nb:',nb,'ntotal:',na+nb
      WRITE(*,'(A6,I8,6X,A6,I8,6X,A6,I8)')'na:',na,'nb:',nb,'ntotal:',na+nb

      
      WRITE(16,'(A6,F7.1,6X,A6,F7.1,6X,A6,F7.1,6X,A13,F7.1)')&
               'Naa0:',Naa0,'Nbb0:',Nbb0,'Nab0:',Nab0,'N_total:',Nabc0_total
      WRITE(*,'(A6,F7.1,6X,A6,F7.1,6X,A6,F7.1,6X,A13,F7.1)')&
               'Naa0:',Naa0,'Nbb0:',Nbb0,'Nab0:',Nab0,'N_total:',Nabc0_total
      

      WRITE(16,'(A6,F7.2,6X,A6,F7.3,6X,A7,I8,6X,A6,I8,6X,A7,I2,6X,A6,I2)')&
      'r02:',r02,'nn_av:',real(nn_total)/real(n_cell),'nn_max:',nn_max,'n_max:',n_max,&
      'nn_min:',nn_min,'n_min:',n_min      
      WRITE(*,'(A6,F7.2,6X,A6,F7.3,6X,A7,I8,6X,A6,I8,6X,A7,I2,6X,A6,I2)')&
      'r02:',r02,'nn_av:',real(nn_total)/real(n_cell),'nn_max:',nn_max,'n_max:',n_max,&
      'nn_min:',nn_min,'n_min:',n_min
      

      CLOSE(15)
      CLOSE(16)
      
      END SUBROUTINE number_neighbor_ini

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!   SUBROUTINE value_H()
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE value_H1()
      IMPLICIT NONE
      
      !!! value H1 
      H1a=0.  ; H1b=0.
      DO j=1,nn(i)
            H1a=H1a+Sa(name_in(i,j))
            H1b=H1b+Sb(name_in(i,j))
      END DO

      END SUBROUTINE value_H1
!!!================================================
      SUBROUTINE value_H2()
      IMPLICIT NONE

      !!! value H2
      H2a=0.  ; H2b=0.
      DO j=1,nn(ic)
            H2a=H2a+Sa(name_in(ic,j))
            H2b=H2b+Sb(name_in(ic,j))
      END DO

      END SUBROUTINE value_H2

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!   SUBROUTINE equi_lattice1()
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE equi_lattice1()
      IMPLICIT NONE

      DO i_loop=1,n_equi1
            CALL equi_lattice()
            !CALL value_thermal()
      END DO

      END SUBROUTINE equi_lattice1

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!   SUBROUTINE equi_lattice()
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE equi_lattice()
      IMPLICIT NONE
          
      E1_old=0. ; E2_old=0. ; E1_new=0. ; E2_new=0.

      DO i=1,n_cell
      
           !!! su dung ic bat ky
            !ic=0      
            !DO WHILE (int(Sa(ic)+Sb(ic))==0)
           
            !CALL random_number(rdn_ic)
            !ic=int(rdn_ic*n_cell_p)
                      
            !END DO
            
            !!! Su dung ic la neighbor cua i:
            IF (nn(i)>0) THEN
            
            ic=0 
            DO WHILE (int(Sa(ic)+Sb(ic))==0)             
                  CALL random_number(rdn_ic)
                  j_ic=int(rdn_ic*real(nn(i)+1))
                  ic=name_in(i,j_ic)
            END DO
            
            CALL value_H1()
                
            E1_old=-Sa(i)*(Jaa*H1a+Jab*H1b)-Sb(i)*(Jab*H1a+Jbb*H1b)   

!!========================================================================================        
            IF (int(Sa(i))/=int(Sa(ic))) THEN           

                 CALL value_H2()
                 E2_old=-Sa(ic)*(Jaa*H2a+Jab*H2b)-Sb(ic)*(Jab*H2a+Jbb*H2b)

                 Sa_mtp=Sa(i)  ;  Sb_mtp=Sb(i)
                 Sa(i)=Sa(ic)  ;  Sb(i)=Sb(ic)
                 Sa(ic)=Sa_mtp ;  Sb(ic)=Sb_mtp

                 CALL value_H1()
                 CALL value_H2()            
                 E1_new=-Sa(i)*(Jaa*H1a+Jab*H1b)-Sb(i)*(Jab*H1a+Jbb*H1b) 
                 E2_new=-Sa(ic)*(Jaa*H2a+Jab*H2b)-Sb(ic)*(Jab*H2a+Jbb*H2b)

                 CALL random_number(rdn_mtp)

                 IF (exp((E1_old+E2_old-E1_new-E2_new)/T) > rdn_mtp) THEN
                        E1_old=E1_new

                 ELSE         
                        Sa(ic)=Sa(i) ;  Sb(ic)=Sb(i)
                        Sa(i)=Sa_mtp ;  Sb(i)=Sb_mtp

                 END IF

           END IF           
!!========================================================================================                    
      END IF
      
      END DO     

                     
      END SUBROUTINE equi_lattice   

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!   SUBROUTINE value_thermal()
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE value_thermal()
      IMPLICIT NONE
          
      energy=0. ; E1_old=0. ; E2_old=0. ; E1_new=0. ; E2_new=0.
      Naa=0. ; Nbb=0. ; Nab=0.

      DO i=1,n_cell
      
            !!! su dung ic bat ky
            !ic=0 
            !DO WHILE (int(Sa(ic)+Sb(ic))==0)
           
            !CALL random_number(rdn_ic)
            !ic=int(rdn_ic*n_cell_p)
                      
            !END DO
            !!!-----------------
            
            IF (nn(i)>0) THEN
            
            !!! Su dung ic la neighbor cua i:
            ic=0 
            DO WHILE (int(Sa(ic)+Sb(ic))==0)
             
            CALL random_number(rdn_ic)
            j_ic=int(rdn_ic*real(nn(i)+1))            
            ic=name_in(i,j_ic)

            END DO 


            !WRITE(*,*)'HAHHAHHHAHHH j_ic=',j_ic
            CALL value_H1()
                
            E1_old=-Sa(i)*(Jaa*H1a+Jab*H1b)-Sb(i)*(Jab*H1a+Jbb*H1b)   

!!========================================================================================        
            IF (int(Sa(i))/=int(Sa(ic))) THEN           

                 CALL value_H2()
                 E2_old=-Sa(ic)*(Jaa*H2a+Jab*H2b)-Sb(ic)*(Jab*H2a+Jbb*H2b)

                 Sa_mtp=Sa(i)  ;  Sb_mtp=Sb(i)
                 Sa(i)=Sa(ic)  ;  Sb(i)=Sb(ic)
                 Sa(ic)=Sa_mtp ;  Sb(ic)=Sb_mtp

                 CALL value_H1()
                 CALL value_H2()            
                 E1_new=-Sa(i)*(Jaa*H1a+Jab*H1b)-Sb(i)*(Jab*H1a+Jbb*H1b) 
                 E2_new=-Sa(ic)*(Jaa*H2a+Jab*H2b)-Sb(ic)*(Jab*H2a+Jbb*H2b)

                 CALL random_number(rdn_mtp)

                 IF (exp((E1_old+E2_old-E1_new-E2_new)/T) > rdn_mtp) THEN
                        E1_old=E1_new

                 ELSE         
                        Sa(ic)=Sa(i) ;  Sb(ic)=Sb(i)
                        Sa(i)=Sa_mtp ;  Sb(i)=Sb_mtp

                 END IF

           END IF
!!========================================================================================                    
           
           energy=energy+E1_old

            !!! Tinh Naa, Nbb, Nab : -----------------------------------
            CALL value_H1()
            Naa=Naa+Sa(i)*H1a ; Nbb=Nbb+Sb(i)*H1b ; Nab=Nab+Sa(i)*H1b+Sb(i)*H1a
            !!!-----------------------------------------------------

            END IF
            
      END DO     

      energy=energy/2./n_cell

      Naa=Naa/2. ; Nbb=Nbb/2. ; Nab=Nab/2.
      
      END SUBROUTINE value_thermal

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!   SUBROUTINE average_thermal()
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE average_thermal()

      IMPLICIT NONE

      E_av=0.; Naa_av=0.; Nbb_av=0. ; Nab_av=0.

      DO i_times=1,n_times

            DO i_loop1=1,n_equi2
                  CALL equi_lattice()
            END DO
            
            DO i_loop2=1,n_average
                  CALL value_thermal()        
                 
                  E_av=E_av+energy
                  Naa_av=Naa_av+Naa ; Nbb_av=Nbb_av+Nbb ; Nab_av=Nab_av+Nab
                                      
            END DO
      
      END DO

      E_av=E_av/real(n_times*n_average)
      
      Naa_av=Naa_av/real(n_times*n_average)
      Nbb_av=Nbb_av/real(n_times*n_average)
      Nab_av=Nab_av/real(n_times*n_average)   
      
      Nabc_total=(Naa_av+Nbb_av+Nab_av)
      
      delNabc=(Naa0-Naa_av)**2.+(Nbb0-Nbb_av)**2.+(Nab0-Nab_av)**2.
            
      !WRITE(*,*)'Naa-av=',Naa_av,Nbb_av,Nab_av,Nabc_total
      
      WRITE(21,'(I3,2X,F8.4,2X,F8.4,2X,F8.4,2X,4F12.4,1X,F16.3)') &
                i_delN,Jaa,Jbb,Jab,Naa_av,Nbb_av,Nab_av,Nabc_total,delNabc
      
      WRITE(Ligne23,*) T,E_av
      WRITE(23,'(a)') trim(Ligne23)

     
      END SUBROUTINE average_thermal

      END PROGRAM main_bio
      

      

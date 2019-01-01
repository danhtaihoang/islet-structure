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
!!!   06.08.2013: Xac dinh so neighbor bang mat, bo sung i0
!!!   09.08.2013: Xac dinh nearest neighbor bang PP sap xep khoang cach tu be den lon
!!!   25.08.2013: Su dung ic de so sanh la random   
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%! 
      PROGRAM main_bio
      IMPLICIT NONE

      CHARACTER (LEN=3)  :: SAt
      CHARACTER (LEN=150):: name_data,CONFIG_INI,FIND_NN
      CHARACTER (LEN=50) :: name
      CHARACTER (LEN=15) :: tmp
      CHARACTER (256)    :: Ligne23,Ligne24
      REAL    (KIND=8),PARAMETER :: nul=0.

      INTEGER (KIND=8):: i,j,n_cell,n_max,i_n,ic,S,i_Jbb,n_Jbb,i_Jab,n_Jab,number,i0,i0_n,k,i_n2
      INTEGER (KIND=8):: i_delN,i_loop,i_loop1,i_loop2,n_equi1,n_equi2,n_average,i_times,n_times
      INTEGER (KIND=8):: nn_total,nn_max,na,nb,n_T,i_T,ia,ib,n_min,nn_min,tmp2,j_ic,name_tmp,j_i
     
      REAL    (KIND=8):: r02,r2,n_cell_p,delT,T,Tmin,Tmax,rdn_config,r02max
      
      REAL    (KIND=8):: energy,E1_old,E2_old,E1_new,E2_new,Naa,Nbb,Nab,H1a,H1b,H2a,H2b
      REAL    (KIND=8):: Jaa,Jbb,Jab,del_Jbb,Jbbmax,Jbbmin,del_Jab,Jabmax,Jabmin
      REAL    (KIND=8):: rdn_ic,rdn_mtp,Naa_av,Nbb_av,Nab_av
      REAL    (KIND=8):: E_av,Nabc_total,delNabc,Naa0,Nbb0,Nab0,Nabc0_total
      REAL    (KIND=8):: pi,alpha,alpha0,OA2,OB2,AB2,tmp1,na_p,real_nb

      INTEGER (KIND=8),DIMENSION(:),ALLOCATABLE :: nn,spin1,spin2,nn1,name_a,name_b
      INTEGER (KIND=8),DIMENSION(:,:),ALLOCATABLE :: name_in,a
      REAL (KIND=8),DIMENSION(:),ALLOCATABLE :: Sa,Sb,x,y,z,Sa0,Sb0
      REAL (KIND=8),DIMENSION(:,:),ALLOCATABLE :: delNabc_tab,r                    
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
          
      ALLOCATE(x(n_cell),y(n_cell),z(n_cell))
      ALLOCATE(Sa(n_cell),Sb(n_cell))
      ALLOCATE(nn(n_cell))
      ALLOCATE(nn1(0:n_cell))
      ALLOCATE(delNabc_tab(n_Jbb,n_Jab))
      ALLOCATE(spin1(n_cell),spin2(n_cell))

      ALLOCATE(r(n_cell,0:n_cell))
      ALLOCATE(a(n_cell,0:n_cell))
      ALLOCATE(Sa0(n_cell),Sb0(n_cell))

      CALL read_data()
      CALL read_name_ab()     
   

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
      na_p=real(na)+1.
      real_nb=real(nb)
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! ====== MAIN PROGRAM === MAIN PROGRAM === MAIN PROGRAM === MAIN PROGRAM ======
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
      
      CALL write_config_ini_3D()
      CALL number_neighbor_ini()
      CALL nearest_neighbor()
      CALL number_interaction_ini()
      
      !!! Doi voi i-consider (i0)
      IF (FIND_NN=='YES') THEN
            CALL neighbor_i0()
            WRITE(*,*)'Finish find NN of i0'
      END IF

      OPEN(unit=21,file='average_thermal.dat')
      OPEN(unit=22,file='del_Nabc2_min.dat')
      OPEN(unit=23,file='energy.dat')
      OPEN(unit=24,file='E_times.dat')             
      
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

            CALL load_config_ini()
            
            CALL average_thermal()
       
            !CALL value_thermal()
                  
            !WRITE(*,*)'Naa=',Naa,'Nbb=',Nbb,'Nab=',Nab
            
            !CALL write_config_3D()

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
      
      DO j_i=1,nb
            i=name_b(j_i)
      
            !!! su dung ic la a-cell bat ky
            j_ic=0 
            DO WHILE (j_ic==0)           
                  CALL random_number(rdn_ic)
                  j_ic=int(rdn_ic*na_p)
                  ic=name_a(j_ic)
            END DO

            CALL value_H1()
            CALL value_H2()
    
            E1_old=-Jab*H1a-Jbb*H1b
            E2_old=-Jaa*H2a-Jab*H2b
                            
                 Sa(i)=1.  ;  Sb(i)=0.
                 Sa(ic)=0. ;  Sb(ic)=1.

                 CALL value_H1()
                 CALL value_H2()            
                 E1_new=-Jaa*H1a-Jab*H1b
                 E2_new=-Jab*H2a-Jbb*H2b

                 CALL random_number(rdn_mtp)

                 IF (exp((E1_old+E2_old-E1_new-E2_new)/T) > rdn_mtp) THEN
                        name_tmp=name_a(j_ic)
                        name_a(j_ic)=name_b(j_i)
                        name_b(j_i)=name_tmp
                 ELSE         
                        Sa(i)=0.  ;  Sb(i)=1.
                        Sa(ic)=1. ;  Sb(ic)=0.

                 END IF                  
            
      END DO            
      END SUBROUTINE equi_lattice   

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!   SUBROUTINE value_thermal()
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE value_thermal()
      IMPLICIT NONE
          
      energy=0. ;  Naa=0. ; Nbb=0. ; Nab=0.

      DO j_i=1,nb
            i=name_b(j_i)
      
            !!! su dung ic la a-cell bat ky
            j_ic=0 
            DO WHILE (j_ic==0)           
                  CALL random_number(rdn_ic)
                  j_ic=int(rdn_ic*na_p)
                  ic=name_a(j_ic)
            END DO

            CALL value_H1()
            CALL value_H2()
    
            E1_old=-Jab*H1a-Jbb*H1b
            E2_old=-Jaa*H2a-Jab*H2b
                            
                 Sa(i)=1.  ;  Sb(i)=0.
                 Sa(ic)=0. ;  Sb(ic)=1.

                 CALL value_H1()
                 CALL value_H2()            
                 E1_new=-Jaa*H1a-Jab*H1b
                 E2_new=-Jab*H2a-Jbb*H2b

                 CALL random_number(rdn_mtp)

                 IF (exp((E1_old+E2_old-E1_new-E2_new)/T) > rdn_mtp) THEN
                        E1_old=E1_new
                        name_tmp=name_a(j_ic)
                        name_a(j_ic)=name_b(j_i)
                        name_b(j_i)=name_tmp
                 ELSE         
                        Sa(i)=0.  ;  Sb(i)=1.
                        Sa(ic)=1. ;  Sb(ic)=0.

                 END IF                  
         
           energy=energy+E1_old
            
      END DO
      energy=energy/2./real_nb

      !!! Tinh Naa, Nbb, Nab : -----------------------------------
      DO i=1,n_cell
            CALL value_H1()
            Naa=Naa+Sa(i)*H1a ; Nbb=Nbb+Sb(i)*H1b ; Nab=Nab+Sa(i)*H1b+Sb(i)*H1a
      !!!-----------------------------------------------------
      END DO   

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

      WRITE(Ligne24,*) i_times,E_av/real(i_times*n_average)
      WRITE(24,'(a)') trim(Ligne24)
      
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
      READ(11, '(A30,(A10))')   tamp, FIND_NN
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A30,(A10))')   tamp, name_data
      READ(11, '(A30,(F12.6))') tamp, r02
      READ(11, '(A30,(F12.6))') tamp, r02max
      READ(11, '(A30,(I5))')    tamp, i0
      READ(11, '(A30,(F12.6))') tamp, alpha0
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

 !!!==========================================================
!!! Hinh dang ban dau nhu data goc
      IF (CONFIG_INI=='YES') THEN

      Sa(:)=0. ; Sb(:)=0.
      
      DO i=1,n_cell
            Sa(i)=Sa0(i)
            Sb(i)=Sb0(i)
      END DO

      END IF
      
!!!========================================================== 
!!! Hinh dang ban dau bat ky
      IF (CONFIG_INI=='NO') THEN

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

      END IF
!!!==========================================================  

      END SUBROUTINE load_config_ini
   
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! WRITE initial position configuration in 3D
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE write_config_ini_3D()
      IMPLICIT NONE

      OPEN(unit=13,file='config_ini_3D/config_ini_3D.pdb')
      
      DO i=1,n_cell
           
            IF (i==i0) THEN 
                  SAt='H'
                  WRITE(13,'("ATOM",9X,A3,14X,3F8.3,1X,"1.00",1X,F8.3)') &
                        SAt,x(i),y(i),z(i),nul
            ELSE          
           
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
            OPEN(unit=12,file='M01.txt')
            n_cell=893
      END IF 

      IF (name_data=='M02') THEN     
            OPEN(unit=12,file='M02.txt')
            n_cell=2598
      END IF 
      
      IF (name_data=='M03') THEN     
            OPEN(unit=12,file='M03.txt')
            n_cell=1870
      END IF
      
      IF (name_data=='M04') THEN     
            OPEN(unit=12,file='M04.txt')
            n_cell=1680
      END IF
      
      IF (name_data=='M05') THEN     
            OPEN(unit=12,file='M05.txt')
            n_cell=518
      END IF
      
      IF (name_data=='M06') THEN     
            OPEN(unit=12,file='M06.txt')
            n_cell=1895
      END IF 

      IF (name_data=='M07') THEN     
            OPEN(unit=12,file='M07.txt')
            n_cell=571
      END IF 

      IF (name_data=='M08') THEN     
            OPEN(unit=12,file='M08.txt')
            n_cell=1192
      END IF

      IF (name_data=='M09') THEN     
            OPEN(unit=12,file='M09.txt')
            n_cell=1153
      END IF 

      IF (name_data=='M10') THEN     
            OPEN(unit=12,file='M10.txt')
            n_cell=1063
      END IF 

      IF (name_data=='M11') THEN     
            OPEN(unit=12,file='M11.txt')
            n_cell=2143
      END IF

      IF (name_data=='M12') THEN     
            OPEN(unit=12,file='M12.txt')
            n_cell=1576
      END IF

      IF (name_data=='M13') THEN     
            OPEN(unit=12,file='M13.txt')
            n_cell=2471
      END IF

      IF (name_data=='M14') THEN     
            OPEN(unit=12,file='M14.txt')
            n_cell=1113
      END IF

      IF (name_data=='M15') THEN     
            OPEN(unit=12,file='M15.txt')
            n_cell=1230
      END IF

      IF (name_data=='M16') THEN     
            OPEN(unit=12,file='M16.txt')
            n_cell=3294
      END IF

      IF (name_data=='M17') THEN     
            OPEN(unit=12,file='M17.txt')
            n_cell=1645
      END IF

      IF (name_data=='M18') THEN     
            OPEN(unit=12,file='M18.txt')
            n_cell=4159
      END IF

      IF (name_data=='M19') THEN     
            OPEN(unit=12,file='M19.txt')
            n_cell=2248
      END IF

      IF (name_data=='M20') THEN     
            OPEN(unit=12,file='M20.txt')
            n_cell=1182
      END IF

      IF (name_data=='M21') THEN     
            OPEN(unit=12,file='M21.txt')
            n_cell=4010
      END IF

      IF (name_data=='M22') THEN     
            OPEN(unit=12,file='M22.txt')
            n_cell=3884
      END IF

      IF (name_data=='M23') THEN     
            OPEN(unit=12,file='M23.txt')
            n_cell=3956
      END IF

      IF (name_data=='M24') THEN     
            OPEN(unit=12,file='M24.txt')
            n_cell=2097
      END IF

      IF (name_data=='M25') THEN     
            OPEN(unit=12,file='M25.txt')
            n_cell=3209
      END IF

      IF (name_data=='M26') THEN     
            OPEN(unit=12,file='M26.txt')
            n_cell=9314
      END IF

      IF (name_data=='M27') THEN     
            OPEN(unit=12,file='M27.txt')
            n_cell=1736
      END IF

      IF (name_data=='M28') THEN     
            OPEN(unit=12,file='M28.txt')
            n_cell=1430
      END IF

      IF (name_data=='M29') THEN     
            OPEN(unit=12,file='M29.txt')
            n_cell=2224
      END IF

      IF (name_data=='M30') THEN     
            OPEN(unit=12,file='M30.txt')
            n_cell=1631
      END IF

      END SUBROUTINE open_data    
!!!=======================================================================
      SUBROUTINE read_data()
      IMPLICIT NONE
  
      Sa(:)=0. ; Sb(:)=0.
      na=0 ; nb=0
      Sa0(:)=0. ; Sb0(:)=0.
      DO i=1,n_cell
            READ(12,*)S,x(i),y(i),z(i)
            
            IF (S==11) THEN
                  Sa(i)=0. ; Sb(i)=1.
                  Sa0(i)=0. ; Sb0(i)=1.    
                  
                  ELSE
                  
                  IF (S==12) THEN
                  Sa(i)=1. ; Sb(i)=0.
                  Sa0(i)=1. ; Sb0(i)=0.
                  END IF
                  
            END IF
          na=na+int(Sa(i)) ; nb=nb+int(Sb(i))        
              
      END DO
      
      CLOSE(12)       

      END SUBROUTINE read_data

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!  SUBROUTINE read_name_ab()
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE read_name_ab()
      IMPLICIT NONE
  
      ALLOCATE(name_a(na))
      ALLOCATE(name_b(nb))

      ia=0 ; ib=0
      name_a(:)=0 ; name_b(:)=0

      DO i=1,n_cell
            IF (int(Sa0(i))==1) THEN
                  ia=ia+1
                  name_a(ia)=i 
            ELSE
                  IF (int(Sb0(i))==1) THEN
                  ib=ib+1
                  name_b(ib)=i 
                  END IF
            END IF
      END DO     

      END SUBROUTINE read_name_ab

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! Number nearest neighbor from experimental data
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE number_neighbor_ini()
      IMPLICIT NONE
 
      OPEN(unit=15,file='nn_i_cell.dat') 
      
      !!!--------------------------------------------------------------
      !!! 
      nn_total=0 ;  nn(:)=0
      
      DO i=1,n_cell
            DO j=1,n_cell
           
            r2=(x(i)-x(j))**2.+(y(i)-y(j))**2.+(z(i)-z(j))**2.
           
            IF ((0.1<r2).and.(r2<r02)) THEN
                  nn(i)=nn(i)+1
            END IF

            END DO
                                                !!! nn(i): number neigbors of i cell
            nn_total=nn_total+nn(i)             !!! number neigbors of system
            
            WRITE(15,*)i,nn(i)

      END DO

      !!!--------------------------------------------------------
      !!! Assign name for nearest neighbor of cell i
      
      ALLOCATE(name_in(n_cell,0:n_cell))
      
      name_in(:,:)=0
      
      DO i=1,n_cell
           i_n=0
           
            DO j=1,n_cell
           
                  r2=(x(i)-x(j))**2.+(y(i)-y(j))**2.+(z(i)-z(j))**2.
                 
                  IF ((0.1<r2).and.(r2<r02)) THEN
                        i_n=i_n+1
                        name_in(i,i_n)=j         

                  END IF
                    
            END DO

      END DO

      CLOSE(15)

      
      END SUBROUTINE number_neighbor_ini

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! Number interaction ini
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE number_interaction_ini()
      IMPLICIT NONE

      OPEN(unit=16,file='Nab0_ini.dat')    

      nn_total=0
 
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

      WRITE(16,'(A6,F7.2,6X,A6,F7.3,6X,A7,I8,6X,A6,I8,6X,A7,I2,6X,A6,I2)')&
      'r02:',r02,'nn_av:',2.*real(Nabc0_total)/real(n_cell),'nn_max:',nn_max,'n_max:',n_max,&
      'nn_min:',nn_min,'n_min:',n_min      
      WRITE(*,'(A6,F7.2,6X,A6,F7.3,6X,A7,I8,6X,A6,I8,6X,A7,I2,6X,A6,I2)')&
      'r02:',r02,'nn_av:',2.*real(Nabc0_total)/real(n_cell),'nn_max:',nn_max,'n_max:',n_max,&
      'nn_min:',nn_min,'n_min:',n_min

      WRITE(16,'(A6,F7.1,6X,A6,F7.1,6X,A6,F7.1,6X,A13,F7.1)')&
               'Naa0:',Naa0,'Nbb0:',Nbb0,'Nab0:',Nab0,'N_total:',Nabc0_total
      WRITE(*,'(A6,F7.1,6X,A6,F7.1,6X,A6,F7.1,6X,A13,F7.1)')&
               'Naa0:',Naa0,'Nbb0:',Nbb0,'Nab0:',Nab0,'N_total:',Nabc0_total
      
      CLOSE(16)

      WRITE(*,'(A6,F7.1,6X,A6,F7.1,6X,A6,F7.1)')&
                 'Naa-eq=',Nabc0_total*(real(na)/real(na+nb))**2.,&
                'Nbb-eq=',Nabc0_total*(real(nb)/real(na+nb))**2.,&
                'Nab-eq=',Nabc0_total*2.*real(na*nb)/real(na+nb)**2.
               
     ! WRITE(*,*)'NN i-star:', nn(i0)
      
      
      END SUBROUTINE number_interaction_ini


!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! 09.08.2013: nearest neighbor
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE nearest_neighbor()
      IMPLICIT NONE

      OPEN(unit=31,file='ri_in.dat')
      OPEN(unit=32,file='config_final_nn.pdb')
      OPEN(unit=33,file='nn_i_final.dat')

      !!!---------------------------------------------------------
      !!! Tim khoang cach r(i,i_n)
      DO i=1,n_cell
            DO i_n=1,nn(i) 
            r(i,i_n)=(x(i)-x(name_in(i,i_n)))**2.+(y(i)-y(name_in(i,i_n)))**2.+(z(i)-z(name_in(i,i_n)))**2.
            END DO
      END DO

      !!!---------------------------------------------------------
      !!! Sap xep theo thu tu tang dan khoang cach r(i,i_n)
      DO i=1,n_cell
            DO i_n=1,nn(i)-1

            k=i_n

            DO j=i_n+1,nn(i)
                  IF (r(i,j)<r(i,k)) THEN
                        k=j
                  END IF
            END DO

            IF (k/=i_n) THEN
                  tmp1=r(i,i_n)
                  r(i,i_n)=r(i,k)
                  r(i,k)=tmp1

                  tmp2=name_in(i,i_n)
                  name_in(i,i_n)=name_in(i,k)
                  name_in(i,k)=tmp2

            END IF

            END DO
      END DO

      !! Test for i0:
      DO i_n=1,nn(i0)
            WRITE(31,*)i_n,r(i0,i_n),name_in(i0,i_n)
            !WRITE(*,*)i_n,r(i0,i_n),name_in(i0,i_n)
      END DO

      !!!---------------------------------------------------------
      !!! Kiem tra goc de tim cac nearest neighbor
      pi=acos(-1.)

      a(:,:)=1
      DO i=1,n_cell
      DO i_n=1,nn(i)-1

      IF (a(i,i_n)==1) THEN
      DO i_n2=i_n+1,nn(i)

      OA2=(x(i)-x(name_in(i,i_n)))**2.+(y(i)-y(name_in(i,i_n)))**2.+(z(i)-z(name_in(i,i_n)))**2.
      OB2=(x(i)-x(name_in(i,i_n2)))**2.+(y(i)-y(name_in(i,i_n2)))**2.+(z(i)-z(name_in(i,i_n2)))**2.
      AB2=(x(name_in(i,i_n))-x(name_in(i,i_n2)))**2.+(y(name_in(i,i_n))-y(name_in(i,i_n2)))**2.&
          +(z(name_in(i,i_n))-z(name_in(i,i_n2)))**2.

      alpha=acos((OA2+OB2-AB2)/(2.*sqrt(OA2)*sqrt(OB2)))*180./pi

      !IF ((i==112).and.(i_n==1).and.(i_n2==3)) THEN
       !     WRITE(*,*)'alpha=',alpha
     ! END IF
      IF (alpha<alpha0) THEN
            a(i,i_n2)=0
      END IF

      END DO
      END IF

      END DO
      END DO

      !!!---------------------------------------------------------
      !!! Sap xep lai (theo thu tu khoang cach tang dan) sau khi loai bo
      DO i=1,n_cell
            k=0
            DO i_n=1,nn(i)

            IF (a(i,i_n)==1) THEN
                  k=k+1
                  r(i,k)=r(i,i_n)
                  name_in(i,k)=name_in(i,i_n)
            END IF
            END DO
            
            !! Update nn(i):
            nn(i)=k            

      END DO
      
      !DO i_n=1,nn(i0)
            !WRITE(*,*)'Final:',i_n,r(i0,i_n),name_in(i0,i_n)
      !END DO

      !!!Write for i0

      
           SAt='H'
           WRITE(32,'("ATOM",9X,A3,14X,3F8.3,1X,"1.00",1X,F8.3)') &
           SAt,x(i0),y(i0),z(i0),nul


          DO i_n=1,nn(i0)
           SAt='Au'
                  WRITE(32,'("ATOM",9X,A3,14X,3F8.3,1X,"1.00",1X,F8.3)') &
                        SAt,x(name_in(i0,i_n)),y(name_in(i0,i_n)),z(name_in(i0,i_n)),nul

          END DO
                            
      !!! Number neighbors
      nn1(:)=0

      DO j=1,MAXVAL(nn)
            DO i=1,n_cell
                  IF (nn(i)==j) THEN
                        nn1(j)=nn1(j)+1            
                  END IF
            END DO
      WRITE(33,*)j,nn1(j)
      WRITE(*,*)j,nn1(j)

      END DO

      CLOSE(31)
      CLOSE(32)
      CLOSE(33)

      END SUBROUTINE nearest_neighbor

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! neighbor of i0
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE neighbor_i0()
      IMPLICIT NONE
 
      spin1(:)=0 ; spin2(:)=0 ; i0_n=0
       
      OPEN(unit=17,file='config_ini_i0.pdb')
      OPEN(unit=18,file='ri0.dat')
           
      DO i=1,n_cell            
            r2=(x(i0)-x(i))**2.+(y(i0)-y(i))**2.+(z(i0)-z(i))**2.
           
            IF ((0.1<r2).and.(r2<r02)) THEN
                  spin1(i)=1
                  i0_n=i0_n+1                
            ELSE                   
            IF ((0.1<r2).and.(r2<r02max)) THEN
                  spin2(i)=1                                                                    
            END IF
            END IF

      END DO
      
      WRITE(*,*)'Number neigbor of i0 cell',i0_n

      WRITE(18,*)i0,r02,i0_n

!!!--------------------------------------------------
!!!Write config neigbor of i0      
      DO i=1,n_cell           
            IF (i==i0) THEN 
                  SAt='H'
                  WRITE(17,'("ATOM",9X,A3,14X,3F8.3,1X,"1.00",1X,F8.3)') &
                        SAt,x(i),y(i),z(i),nul
            END IF            
      END DO
      
      DO i=1,n_cell
            IF (spin1(i)==1) THEN
                  WRITE(*,*)i

                  SAt='Au'
                  WRITE(17,'("ATOM",9X,A3,14X,3F8.3,1X,"1.00",1X,F8.3)') &
                        SAt,x(i),y(i),z(i),nul
            END IF
                              
            IF (spin2(i)==1) THEN 
                  SAt='Cu'
                  WRITE(17,'("ATOM",9X,A3,14X,3F8.3,1X,"1.00",1X,F8.3)') &
                        SAt,x(i),y(i),z(i),nul       
            END IF                        
            
      END DO
                              
      CLOSE(17)
      CLOSE(18)      

      
      END SUBROUTINE neighbor_i0
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END PROGRAM main_bio
      

      

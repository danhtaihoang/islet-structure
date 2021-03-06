      PROGRAM min_delNabc
      IMPLICIT NONE

      INTEGER (KIND=8),PARAMETER::n_line=2601
      INTEGER (KIND=8):: i
      REAL    (KIND=8):: Naa_av,Nbb_av,Nab_av,Jbb_av,Jab_av
      REAL    (KIND=8):: compt,delNabc_av,delJbb,delJab

      REAL (KIND=8),DIMENSION(n_line) :: Jbb,Jab,FNaa,FNbb,FNab,delNabc,cost

!!!=======================================================================================
!!!======================================================================================= 
      OPEN(unit=12,file='average_thermal.dat')
      OPEN(unit=13,file='average_thermal_min.dat')
      OPEN(unit=14,file='average_thermal_plot.dat')

      compt=0. ; Jbb_av=0.; Jab_av=0. ; delNabc_av=0.
      
      DO i=1,n_line
            READ(12,*)Jbb(i),Jab(i),Naa_av,Nbb_av,Nab_av,FNaa(i),FNbb(i),FNab(i),delNabc(i)
            
            cost(i)=delNabc(i)/(FNaa(i)+FNbb(i)+FNab(i))/2.
            
            !IF (cost(i)<1.) THEN
            IF (delNabc(i)<1.) THEN
                  !WRITE(13,'(F8.4,3X,F8.4,4X,F16.4,3X,F8.4)') &
                  !Jbb(i),Jab(i),delNabc(i),cost(i)
                  !WRITE(*,'(F8.4,3X,F8.4,4X,F16.4,3X,F8.4)') &
                  !Jbb(i),Jab(i),delNabc(i),cost(i)
                  
                  WRITE(13,'(F8.4,3X,F8.4,4X,F16.4)') &
                  Jbb(i),Jab(i),delNabc(i)
                  WRITE(*,'(F8.4,3X,F8.4,4X,F16.4)') &
                  Jbb(i),Jab(i),delNabc(i)                 
                                   
                  compt=compt+1.
                  Jbb_av=Jbb_av+Jbb(i)
                  Jab_av=Jab_av+Jab(i)
                  delNabc_av=delNabc_av+delNabc(i)
                              
            END IF
            
            
            WRITE(14,'(F8.4,3X,F8.4,4X,F16.4)') &
                  Jbb(i),Jab(i),delNabc(i)

      END DO
      
      Jbb_av=Jbb_av/compt
      Jab_av=Jab_av/compt
      delNabc_av=delNabc_av/compt
      
     
!!!=======================================================================================
!!! Tinh delJbb, delJab

      delJbb=0. ; delJab=0.            
      DO i=1,n_line                 
            !IF (cost(i)<1.) THEN
            IF (delNabc(i)<1.) THEN                                  
                  delJbb=delJbb+(Jbb(i)-Jbb_av)**2.
                  delJab=delJab+(Jab(i)-Jab_av)**2.                             
            END IF

      END DO
      
      delJbb=sqrt(delJbb/compt)
      delJab=sqrt(delJab/compt)

      WRITE(13,*)'Jbb_av:',Jbb_av,'Jab_av:',Jab_av,'delJbb:',delJbb,'delJab:',delJab,'delNabc_av:',delNabc_av
      WRITE(*,*)'Jbb_av:',Jbb_av,'Jab_av:',Jab_av,'delJbb:',delJbb,'delJab:',delJab,'delNabc_av:',delNabc_av

      CLOSE(12)
      CLOSE(13)
      CLOSE(14)       
            

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END PROGRAM min_delNabc
      

      

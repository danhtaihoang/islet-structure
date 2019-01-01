      PROGRAM min_cost
      IMPLICIT NONE

      INTEGER (KIND=8),PARAMETER::n_line=441
      REAL (KIND=8),PARAMETER::cost0=2.
      INTEGER (KIND=8):: i
      REAL    (KIND=8):: Naa_av,Nbb_av,Nab_av
      REAL    (KIND=8):: Z,Jbb_av,Jab_av,Jbb2_av,Jab2_av,delJbb,delJab

      REAL (KIND=8),DIMENSION(n_line) :: Jbb,Jab,FNaa,FNbb,FNab,cost

!!!=======================================================================================
!!!======================================================================================= 
      OPEN(unit=12,file='average_thermal.dat')
      OPEN(unit=13,file='average_thermal_min.dat')
      OPEN(unit=14,file='average_thermal_plot.dat')
      OPEN(unit=15,file='result_final.dat')
!!!=======================================================================================
!!! Partition function: 
      Z=0.     
      DO i=1,n_line
            READ(12,*)Jbb(i),Jab(i),Naa_av,Nbb_av,Nab_av,FNaa(i),FNbb(i),FNab(i),cost(i)
                              
            Z=Z+exp(-cost(i))
                                                     
            IF (cost(i)<cost0) THEN                                
                  WRITE(13,'(F8.4,3X,F8.4,4X,F16.4)') Jbb(i),Jab(i),cost(i)
                  WRITE(*,'(F8.4,3X,F8.4,4X,F16.4)') Jbb(i),Jab(i),cost(i)                                       
            END IF

            WRITE(14,'(F8.4,3X,F8.4,4X,F16.4)') Jbb(i),Jab(i),cost(i)

      END DO
 
!!!=======================================================================================
!!! Expectation value by using normal distribution 
      Jbb_av=0. ;  Jab_av=0. ; Jbb2_av=0. ;  Jab2_av=0.
      DO i=1,n_line                                                
                  Jbb_av=Jbb_av+Jbb(i)*exp(-cost(i))
                  Jab_av=Jab_av+Jab(i)*exp(-cost(i))
                  
                  Jbb2_av=Jbb2_av+Jbb(i)*Jbb(i)*exp(-cost(i))
                  Jab2_av=Jab2_av+Jab(i)*Jab(i)*exp(-cost(i))
                                                          
      END DO 
 
      Jbb_av=Jbb_av/Z ; Jab_av=Jab_av/Z
      Jbb2_av=Jbb2_av/Z ; Jab2_av=Jab2_av/Z
      
      delJbb=sqrt(Jbb2_av-Jbb_av**2.)
      delJab=sqrt(Jab2_av-Jab_av**2.)
      
      WRITE(13,*)'Jbb_av:',Jbb_av,'Jab_av:',Jab_av,'delJbb:',delJbb,'delJab:',delJab
      WRITE(*,*)'Jbb_av:',Jbb_av,'delJbb:',delJbb,'Jab_av:',Jab_av,'delJab:',delJab
      
      WRITE(15,'(4F9.3,4X,4F9.3)')Jbb_av,delJbb,Jbb_av-delJbb,Jbb_av+delJbb,Jab_av,delJab,Jab_av-delJab,Jab_av+delJab
      WRITE(*,'(4F9.3,4X,4F9.3)')Jbb_av,delJbb,Jbb_av-delJbb,Jbb_av+delJbb,Jab_av,delJab,Jab_av-delJab,Jab_av+delJab

      CLOSE(12)
      CLOSE(13)
      CLOSE(14)
      CLOSE(15)       
            

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END PROGRAM min_cost
      

      

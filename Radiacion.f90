!     Last change:  EF   11 Mar 2005   10:58 am

MODULE adapt
   
   USE typePrecision
   USE typeMalla
   USE typeArreglo3D
   USE Propiedades
   USE typeBorde
   
   !VARIABLES QUE SE USAN PARA LA RADIACI�N
   REAL(nP), parameter :: PI  = 3.1415926
   REAL(nP), parameter :: Sig = 5.667E-08
   REAL(nP), parameter :: Gra = 0.01116
   REAL(nP), parameter :: Pr  = 0.7
   REAL(nP), parameter :: Amu = 0.0000185
   
   REAL(nP), parameter :: PVa = 496.637
   REAL(nP), parameter :: PVb = -3.4998
   REAL(nP), parameter :: PVc = 0.006187
   
   !Variables usadas en el m�dulo de radiaci�n
   INTEGER           :: NumDat, ContDia
   REAL(nP)          :: TPR,EmiPE,EmiPW,EmiPN,EmiPS,EmiTe,VAir
   REAL(nP)          :: HR,NnPiso
   REAL(nP), POINTER :: TI(:),TA(:),IrrPE(:),IrrTe(:),IrrPW(:),IrrPN(:),IrrPS(:),HRe(:),VAire(:)
   REAL(nP)          :: ICal,INub,CAir,FacN
   REAL(nP)          :: TAmin,TAmax,TiMin,TiMax,TAprom,TiProm,TiCom,DifTA
   REAL(nP)          :: Lat,Decl,Dif,IrrMaxPE,TiMaxPE,IrrMaxPW,TiMaxPW,IrrMaxPN,TiMaxPN,IrrMaxPS,TiMaxPS,IrrMaxTe
   INTEGER           :: ICBTE,ICBK1,ICBN1,ICBI1,ICBL1
   REAL(nP)          :: RefPE,RefTe,RefPW,RefPN,RefPS
   REAL(nP)          :: LatR,DecR,Calc1,Calc2,CalAPS,APS,TiSol,TiSS
   REAL(nP)          :: TiMaxSPE,TiMaxSPW,TiMaxSPN,TiMaxSPS
   REAL(nP)          :: Hinf!,VolT
   
   REAL(nP)          :: TiTotal,Tamb,TiLM,TiSM,Dint
   REAL(nP)          :: IrrGPE,IrrGTe,IrrGPW,IrrGPN,IrrGPS,QPE,QTe,QPW,QPN,QPS
   REAL(nP)          :: Pg,PVair,Tcielo,HSC,HTC
   REAL(nP)          :: R1,R2
   REAL(nP)          :: EMI,EmiPu,EmiV,RefPu,RefV
   
   
   !VENTILACI�N
   INTEGER           :: x1, x2, z1, z2
   LOGICAL           :: Vent
   
   !Aire
   REAL(nP)          :: Area1, Area2, V1, V2
   LOGICAL           :: Aire
   CONTAINS
   
   SUBROUTINE begin(u,v,w,t,TIMER)
      
      USE LeerData
      
      TYPE(arreglo3D) :: u,v,w,t
      INTEGER         :: i,j,k,a,b,c
      REAL(nP)        :: TIMER
      
      INTENT(IN OUT)  :: u,v,w,t,TIMER
      
      ALLOCATE(Tipo%value(ni,nj,nk))
      
      !Tipo = 0  ----->  Aire
      !Tipo = 1  ----->  Pared exterior
      !Tipo = 2  ----->  Ventana
      !Tipo = 3  ----->  Piso
      !Tipo = 4  ----->  Techo
      !Tipo = 5  ----->  Puerta
      !Tipo = 8  ----->  Pared interior
      !Tipo = 9  ----->  Entrepiso
      
      !Condiciones de contorno de velocidad
      KBC_I1 = "INFLOW"
      KBC_L1 = "INFLOW"
      KBC_J1 = "INFLOW"
      KBC_M1 = "INFLOW"
      KBC_K1 = "INFLOW"
      KBC_N1 = "INFLOW"
      ContDia = 0
      
      OPEN (UNIT = 11, FILE = nomRad, STATUS ='OLD')
      
      READ (11,*)
      
      !Lectura de las propiedades nodo por nodo, desde el archivo de radiaci�n
      DO a=1,objmalla%l1
         DO b=1,objmalla%m1
            DO c=1,objmalla%n1
               
               READ (11,*) i,j,k,Tipo%value(i,j,k),gamT%value(i,j,k),Cp%value(i,j,k),rho%value(i,j,k)
               
               IF (Tipo%value(i,j,k).gt.0) THEN
                  visc%value(i,j,k) = 1.0D30
                  rho%value(i,j,k) = rho%value(i,j,k)*Cp%value(i,j,k) * 1000
               ELSE
                  visc%value(i,j,k) = 1.85D-05
                  GamT%value(i,j,k) = GamT%value(i,j,k)/(1000*Cp%value(i,j,k))
               END IF
               
            END DO
         END DO
      END DO
      
      !Calcula el numero de nodos en el piso
      NnPiso=0
      DO b=1,nj
         IF (Tipo%value(Imon,b,Kmon).eq.3) NnPiso = NnPiso + 1
      END DO
      
      READ (11,*)
      !Emisividad de las paredes, puertas y ventanas
      READ (11,*) EmiPE,EmiTe,EmiPW,EmiPN,EmiPS,EmiV,EmiPu
      
      READ (11,*)
      !Condiciones de borde
      READ (11,*) ICBTE,ICBK1,ICBN1,ICBI1,ICBL1
      
      READ (11,*)
      !Reflectividad de las paredes, puertas y ventanas
      READ (11,*) RefPE,RefTe,RefPW,RefPN,RefPS,RefV,RefPu
      
      READ (11,*)
      !Forma de calculo, indicador de nubosidad, y cambios de aire
      READ (11,*) ICal,INub,CAir
      
      IF (ICal.eq.0) THEN
         
         READ (11,*)
         READ (11,*) NumDat
         
         !Lectura de datos de irradiancia,temperatura ambiente, HUMEDAD RELATIVA Y
         !S1 VERTICAL Y HORIZONTAL
         
         ALLOCATE(TI(NumDat) ,TA(NumDat) ,IrrPE(NumDat) ,IrrTe(NumDat) ,IrrPW(NumDat) , &
         IrrPN(NumDat) ,IrrPS(NumDat) ,HRe(NumDat) ,VAire(NumDat))
         
         READ (11,*)
         
         !Lectura de los datos horarios de temperaturas e irradiancia
         DO a=1,NumDat
            
            READ (11,*) TI(a),TA(a),IrrPE(a),IrrTe(a),IrrPW(a),IrrPN(a),IrrPS(a),HRe(a),VAire(a)
            
         END DO
         
      END IF
      
      FacN=(1.0-(0.056*INub))       !Factor de nubosidad
      
      !C.+.+.+.+.+  CALCULOS PREVIOS DE LA TEMPERATURA AMBIENTE  .+.+.+.+.+
      IF (ICal.eq.1) THEN
         
         READ (11,*)
         !Lectura de temperatura maximas y minimas
         READ (11,*) TAmin,TAmax,TiMin,TiMax
         
         TAprom=(TAmin+TAmax)/2.
         TiProm=(TiMin+TiMax)/2.
         TiCom=1440.-(TiMax-TiMin)
         DifTA=TAmax-TAprom
         
         READ (11,*)
         !Lectura de irradiancias maximas y minimas, latitud, declinanacion, diferencaia, velocidad
         !del aire, humedad relativa
         READ (11,*) Lat,Decl,Dif,IrrMaxPE,TiMaxPE,IrrMaxPW,TiMaxPW,IrrMaxPN,TiMaxPN,IrrMaxPS,TiMaxPS,IrrMaxTe,Vair,HR
         
         LatR=Lat*PI/180     !Transforma de grados a radianes
         DecR=Decl*PI/180
         !c.....  calculo del angulo puesta del sol..........
         Calc1=-(TAN(LatR)*TAN(DecR))
         Calc2=SQRT(1./(Calc1)**2.-1.)
         CalAPS=ATAN(Calc2)
         
         IF (CalAPS.LT.0.) CalAPS=CalAPS+PI
         
         !Angulo de la Puesta del Sol en grados
         APS=CalAPS/(PI/180)
         !c....................................................
         TiSol=2.*(APS/15.)*60.
         TiSS=TiSol/2.
         TiMaxSPE=TiMaxPE+Dif
         TiMaxSPW=TiMaxPW+Dif
         TiMaxSPN=TiMaxPN+Dif
         TiMaxSPS=TiMaxPS+Dif
      END IF
      
      !C/////// FIN CALCULOS PREVIOS DE LA IRRADIANCIA //////////////
      
      !C.+.+.+.+.+  CALCULO PREVIO  COEFICIENTE CONVECTIVO .+.+.+.+.+
      Hinf=2.839+3.833*VAir
      
      CLOSE(11)
      
      !Calculo del volumen de aire dentro de la vivienda
      !VolT = 0.0D0
      !DO a=1,ni
      !DO b=1,nj
      !DO c=1,nk
      
      !   IF (tipo%value(a,b,c).eq.0) then
      !      VolT = VolT + objmalla%vol(a,b,c)
      !   END IF
      
      !END DO
      !END DO
      !END DO
      
      !Temperatura inicial para todos los nodos
      t%value(:,:,:) = 303.15
      t%value(:,1,:) = 303.35
      t%value(:,2,:) = 303.35
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !         Vent = .false.
      !         IF (Vent) THEN
      !         DO i = 1, objmalla%l1
      !            IF ((objmalla%x(16) - 0.5).GE.objmalla%x(i).AND.(objmalla%x(16) - 0.5).LT.objmalla%x(i+1)) THEN
      !               x1 = i
      !            END IF
      !            IF ((objmalla%x(16) + 0.5).GT.objmalla%x(i).AND.(objmalla%x(16) + 0.5).LE.objmalla%x(i+1)) THEN
      !               x2 = i + 1
      !            END IF
      !         END DO
      !
      !         DO k = 1, objmalla%n1
      !            IF ((objmalla%z(14) - 0.5).GE.objmalla%z(k).AND.(objmalla%z(14) - 0.5).LT.objmalla%z(k+1)) THEN
      !               z1 = k
      !            END IF
      !            IF ((objmalla%z(14) + 0.5).GT.objmalla%z(k).AND.(objmalla%z(14) + 0.5).LE.objmalla%z(k+1)) THEN
      !               z2 = k + 1
      !            END IF
      !         END DO
      !
      !         DO i = x1, x2
      !         DO k = z1, z2
      !            tipo%value(i,17,k) = 6
      !         END DO
      !         END DO
      !
      !         END IF
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      !         Aire = .false.
      !         IF (Aire) THEN
      !
      !         Area1 = 0.0
      !         Area2 = 0.0
      !         DO i = 16,21
      !         DO j = 16,17
      !         DO k = 18,19
      !            tipo%value(i,j,k) = 7
      !            Area1 = Area1 + objmalla%AreaZ(i,j)
      !         END DO
      !         END DO
      !         END DO
      !
      !         DO i = 16,21
      !         DO j = 14,15
      !         DO k = 18,19
      !            tipo%value(i,j,k) = 14
      !            Area2 = Area2 + objmalla%AreaZ(i,j)
      !         END DO
      !         END DO
      !         END DO
      !          V1 = 5.0
      !          V2 = V1*Area1/Area2
      !          !write (*,*) V2, Area1, Area2
      !         END IF
      !
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111
      
      !Lectura de datos previos obtenidos de corridas anteriores
      !Se leen los valores de u,v,w, presion y temperatura, desde un
      !archivo de resultados obtenido en una corrida previa
      IF (DatPrev) then
         
         OPEN (UNIT=15,FILE= nomPrev)
         
         READ (15,*)
         READ (15,*)
         READ (15,*)
         READ (15,*)
         READ (15,*)
         READ (15,*)
         READ (15,*)
         READ (15,*)
         READ (15,*)
         READ (15,*)
         READ (15,*)
         READ (15,*)
         READ (15,*)
         
         DO c=1,objmalla%n1
            DO b=1,objmalla%m1
               DO a=1,objmalla%l1
                  
                  READ (15,*) PRU,PRU,PRU,u%value(a,b,c),v%value(a,b,c),w%value(a,b,c),temp1%value(a,b,c),t%value(a,b,c),PRU
                  
               END DO
            END DO
         END DO
         
         READ (15,*)
         READ (15,*) TIMER
         CLOSE (15)
         
      END IF
      
      DEALLOCATE (Cp%value)
      
   END SUBROUTINE begin
   
   !SUBROUTINE bound(f,nf,iter)
   
   !   USE typePrecision
   !   USE typeMalla
   !   USE typeArreglo3D
   !   USE UVWT
   !   USE typeBorde
   
   !   INTEGER          :: iter
   
   !   TYPE(arreglo4D)  :: f(nfmax)
   
   !   INTENT (IN OUT)  :: f
   
   !END SUBROUTINE bound
   
   SUBROUTINE phi(nf,timer,u,v,w,t)
      
      USE Coeficientes
      
      INTEGER  :: i,j,k,kk,nf,cont
      REAL(nP) :: TIMER
      
      TYPE(arreglo3D) :: u,v,w,t
      
      INTENT (IN)          :: nf,timer
      INTENT (IN OUT)      :: u,v,w,t
      
      IF (nf.EQ.5) THEN
         gamX.value = gamT.value
      END IF
      
      IF (nf.EQ.0) THEN
         gamX%value(:,:,:) = visc%value(:,:,:)
      END IF
      
      IF(NF.EQ.5)THEN
         
         !Mantiene el valor del tiempo legal entre 0 y 1440 minutos (1 dia)
         TiLM=timer/60. - 1440*INT(TIMER/86400)
         
         IF (ICal.EQ.1)THEN
            !C.+.+.+.+.+  CALCULO DE LA TEMPERATURA AMBIENTE  .+.+.+.+.+
            TiTotal=timer/60.
            Tamb=TEMAMB(TiTotal,PI,TAprom,DifTA,TiProm,TiCom,TiMax,TiMin)
            
            !C/////// FIN CALCULO DE TEMPERATURA AMBIENTE //////////////
            TiSM=TiLM+Dif
            !C.+.+.+.+.+  CALCULO DE LA IRRADIANCIA TECHO Y PAREDES  +.+
            IF(TiLM.GE.450.0.AND.TiLM.LE.1170.0)THEN
               IrrGTe=IRRATE(IrrMaxTe,PI,TiSM,TiSS,TiSol)
               IrrGPE=IRRAPA(TiSol,TiMaxSPE,TiSS,TiLM,TiMaxPE,IrrMaxPE,PI,TiSM)
               IrrGPW=IRRAPA(TiSol,TiMaxSPW,TiSS,TiLM,TiMaxPW,IrrMaxPW,PI,TiSM)
               IrrGPN=IRRAPA(TiSol,TiMaxSPN,TiSS,TiLM,TiMaxPN,IrrMaxPN,PI,TiSM)
               IrrGPS=IRRAPA(TiSol,TiMaxSPS,TiSS,TiLM,TiMaxPS,IrrMaxPS,PI,TiSM)
            ELSE
               IrrGTe=0.
               IrrGPE=0.
               IrrGPW=0.
               IrrGPN=0.
               IrrGPS=0.
            END IF
            !C/////// FIN CALCULO DE IRRADIANCIA TECHO Y PAREDES   /////
         END IF
         
         IF (ICal.EQ.0)THEN
            
            IF(TiLM.EQ.0.) THEN
               
               KK=1
               !C.....Interpolaciones de datos horarios expresados cada media hora.
            END IF
            
            Dint=(TiLM-TI(KK))/(TI(KK+1)-TI(KK))
            IF(TiLM.GE.TI(KK+1)) KK=KK+1
            Tamb=TA(KK)+(TA(KK+1)-TA(KK))*Dint
            HR=HRe(KK)
            Hinf=2.839+3.833*VAire(KK)
            IrrGTe=IrrTe(KK)+(IrrTe(KK+1)-IrrTe(KK))*Dint
            IrrGPE=IrrPE(KK)+(IrrPE(KK+1)-IrrPE(KK))*Dint
            IrrGPW=IrrPW(KK)+(IrrPW(KK+1)-IrrPW(KK))*Dint
            IrrGPN=IrrPN(KK)+(IrrPN(KK+1)-IrrPN(KK))*Dint
            IrrGPS=IrrPS(KK)+(IrrPS(KK+1)-IrrPS(KK))*Dint
         END IF
         
         !C.+.+.+.+.+  CALCULO DE LA TEMPERATURA DEL CIELO  +.+
         Pg=PVa+PVb*Tamb+PVc*Tamb**2.
         PVair=HR*Pg
         TPR=(-PVb+(PVb**2.-4.0*PVc*(PVa-PVair))**0.5)/(2.0*PVc)
         Tcielo=Tamb*(0.8+((TPR-273.)/250.))**0.25
         !C////// FIN CALCULO DE LA TEMPERATURA DEL CIELO   /////
         
         !C.+.+.+.+.+  CONDICIONES DE CONTORNO  +.+
         !c.....techo-piso
         QTe = IrrGTe*(1-RefTe)
         DO k=2,objMalla%n2
            DO i=2,objMalla%l2
               KBC.J1(i,k)=1
               KBC.M1(i,k)=1
               !c..... condici�n de flujo conocido y calculado"2" o  adiabatico"1"
               IF(ICBTE.NE.0)THEN
                  KBC.M1(i,k)=2
                  IF(ICBTE.EQ.2)THEN
                     HSC=EmiTe*Sig*(3.0*t%value(i,objMalla%m1,k)**4.+Tcielo**4.)*FacN &
                     + Hinf*Tamb
                     HTC=4.0*EmiTe*Sig*t%value(i,objMalla%m1,k)**3.*FacN+Hinf
                     IF(HTC.LT.0.) WRITE(*,*)'CUIDADO VALOR DE HTC NEGATIVO'
                     FLX.PM1(I,K)=-HTC
                     FLX.CM1(I,K)=QTe+HSC
                  END IF
               END IF
            END DO
         END DO
         
         !c.....paredes externas
         k = 1
         DO i=2,objMalla%l2
            DO j=2,objMalla%m2
               IF(ICBK1.NE.0)THEN
                  KBC.K1(i,j)=2
                  IF(ICBK1.EQ.2)THEN
                     IF (j.GT.NnPiso)THEN
                        IF (tipo%value(i,j,k).eq.1.or.tipo%value(i,j,k).eq.4) then
                           EMI = EmiPN
                           QPN = IrrGPN*(1-RefPN)
                        END IF
                        IF (tipo%value(i,j,k).eq.2) then
                           EMI = EmiV
                           QPN = IrrGPN*(1-RefV)
                        END IF
                        IF (tipo%value(i,j,k).eq.5) then
                           EMI = EmiPu
                           QPN = IrrGPN*(1-RefPu)
                        END IF
                        R1=-4.*EMI*Sig*t%value(i,j,k)**3./(1.+EMI)
                        R2=EMI*Sig*(3.*t%value(i,j,k)**4.+Tcielo**4.)/(1.+EMI)
                        FLX.PK1(I,J)=-Hinf+R1*FacN
                        FLX.CK1(I,J)=QPN+Hinf*Tamb+R2*FacN
                     END IF
                  END IF
               END IF
            END DO
         END DO
         
         k = objmalla%n1
         DO i=2,objMalla%l2
            DO j=2,objMalla%m2
               IF(ICBN1.NE.0)THEN
                  KBC.N1(i,j)=2
                  IF(ICBN1.EQ.2)THEN
                     IF (j.GT.NnPiso)THEN
                        IF (tipo%value(i,j,k).eq.1.or.tipo%value(i,j,k).eq.4) then
                           EMI = EmiPS
                           QPS = IrrGPS*(1-RefPS)
                        END IF
                        IF (tipo%value(i,j,k).eq.2) then
                           EMI = EmiV
                           QPS = IrrGPS*(1-RefV)
                        END IF
                        IF (tipo%value(i,j,k).eq.5) then
                           EMI = EmiPu
                           QPS = IrrGPS*(1-RefPu)
                        END IF
                        R1=-4.*EMI*Sig*t%value(i,j,k)**3./(1.+EMI)
                        R2=EMI*Sig*(3.*t%value(i,j,k)**4.+Tcielo**4.)/(1.+EMI)
                        FLX.PN1(I,J)=-Hinf+R1*FacN
                        FLX.CN1(I,J)=QPS+Hinf*Tamb+R2*FacN
                     END IF
                  END IF
               END IF
            END DO
         END DO
         
         i = 1
         DO k=2,objMalla%n2
            DO j=2,objMalla%m2
               IF(ICBI1.NE.0)THEN
                  KBC.I1(j,k)=2
                  IF(ICBI1.EQ.2)THEN
                     IF (j.GT.NnPiso)THEN
                        IF (tipo%value(i,j,k).eq.1.or.tipo%value(i,j,k).eq.4) then
                           EMI = EmiPE
                           QPE = IrrGPE*(1-RefPE)
                        END IF
                        IF (tipo%value(i,j,k).eq.2) then
                           EMI = EmiV
                           QPE = IrrGPE*(1-RefV)
                        END IF
                        IF (tipo%value(i,j,k).eq.5) then
                           EMI = EmiPu
                           QPE = IrrGPE*(1-RefPu)
                        END IF
                        R1=-4.*EMI*Sig*t%value(i,j,k)**3./(1.+EMI)
                        R2=EMI*Sig*(3.*t%value(i,j,k)**4.+Tcielo**4.)/(1.+EMI)
                        FLX.PI1(j,k)=-Hinf+R1*FacN
                        FLX.CI1(j,k)=QPE+Hinf*Tamb+R2*FacN
                     END IF
                  END IF
               END IF
            END DO
         END DO
         
         i = objmalla%l1
         DO k=2,objMalla%n2
            DO j=2,objMalla%m2
               IF(ICBL1.NE.0)THEN
                  KBC.L1(j,k)=2
                  IF(ICBL1.EQ.2)THEN
                     IF (j.GT.NnPiso)THEN
                        IF (tipo%value(i,j,k).eq.1.or.tipo%value(i,j,k).eq.4) then
                           EMI = EmiPW
                           QPW = IrrGPW*(1-RefPW)
                        END IF
                        IF (tipo%value(i,j,k).eq.2) then
                           EMI = EmiV
                           QPW = IrrGPW*(1-RefV)
                        END IF
                        IF (tipo%value(i,j,k).eq.5) then
                           EMI = EmiPu
                           QPW = IrrGPW*(1-RefPu)
                        END IF
                        R1=-4.*EMI*Sig*t%value(i,j,k)**3./(1.+EMI)
                        R2=EMI*Sig*(3.*t%value(i,j,k)**4.+Tcielo**4.)/(1.+EMI)
                        FLX.PL1(J,K)=-Hinf+R1*FacN
                        FLX.CL1(J,K)=QPW+Hinf*Tamb+R2*FacN
                     END IF
                  END IF
               END IF
            END DO
         END DO
         
         !c.....velocidad en paredes internas y techos internos
         
         DO i=1,objMalla%l1
            DO j=1,objMalla%m1
               DO k=1,objMalla%n1
                  
                  IF (Tipo%value(i,j,k).eq.8.or.Tipo%value(i,j,k).eq.9) THEN
                     
                     u%value(i,j,k)  = 0.0D0
                     v%value(i,j,k)  = 0.0D0
                     w%value(i,j,k)  = 0.0D0
                     
                  END IF
                  
               END DO
            END DO
         END DO
         
         
         !c.....velocidad en paredes externas, piso y techo externo
         
         DO i=1,objMalla%l1
            DO j=1,objMalla%m1
               DO k=1,objMalla%n1
                  
                  IF (Tipo%value(i,j,k).eq.1.or.Tipo%value(i,j,k).eq.3.or.Tipo%value(i,j,k).eq.4) THEN
                     
                     u%value(i,j,k)  = 0.0D0
                     v%value(i,j,k)  = 0.0D0
                     w%value(i,j,k)  = 0.0D0
                     
                  END IF
                  
               END DO
            END DO
         END DO
         
         
         !C.....Cambios de aire............................................
         IF (CAir.NE.0.)THEN
            DO K=2,objMalla%n2
               DO J=2,objMalla%m2
                  DO I=2,objMalla%l2
                     IF (j.GT.NnPiso)THEN
                        con%value(i,j,k)=(Tamb-t%value(i,j,k))*0.00033*CAir
                     END IF
                  END DO
               END DO
            END DO
         END IF
         
      END IF
      !C......FINALIZA LAS CONDICIONES PARA EL CALCULO DE LA TEMPERATURA
      
      !c.....cambio de velocidad por flotaci�n.......
      IF(NF.EQ.2)THEN
         DO k=2,objMalla%n2
            DO j=2,objMalla%m2
               DO i=2,objMalla%l2
                  IF(tipo%value(i,j,k).eq.0) con%value(i,j,k)=(t%value(i,j,k)-273.)*Gra
                  !IF(tipo%value(i,j,k).eq.6) THEN
                  !   ap%value(i,j,k) = -1.0D30
                  !   con%value(i,j,k) = 1.0D30*(-2.0)
                  !   gamX%value(i,j,k) = 1.0D20
                  !END IF
               END DO
            END DO
         END DO
      END IF
      
      !C/////// FIN CONDICIONES DE CONTORNO  /////
      
      IF (aire) THEN
         
         IF(NF.EQ.3)THEN
            DO k=2,objMalla%n2
               DO j=2,objMalla%m2
                  DO i=2,objMalla%l2
                     IF(tipo%value(i,j,k).eq.7) THEN
                        ap%value(i,j,k) = -1.0D30
                        con%value(i,j,k) = 1.0D30*(-V1)
                        gamX%value(i,j,k) = 1.0D20
                     END IF
                     
                     IF(tipo%value(i,j,k).eq.14) THEN
                        ap%value(i,j,k) = -1.0D30
                        con%value(i,j,k) = 1.0D30*(V2)
                        gamX%value(i,j,k) = 1.0D20
                     END IF
                  END DO
               END DO
            END DO
         END IF
         
         IF(NF.EQ.5)THEN
            DO k=2,objMalla%n2
               DO j=2,objMalla%m2
                  DO i=2,objMalla%l2
                     IF(tipo%value(i,j,k).eq.7) THEN
                        ap%value(i,j,k) = -1.0D30
                        con%value(i,j,k) = 1.0D30*(288)
                        gamX%value(i,j,k) = 1.10D26 !1.15
                     END IF
                  END DO
               END DO
            END DO
         END IF
         
         
      END IF
      
      
   END SUBROUTINE phi
   !  Sc = con
   !  Sp = ap
   ! ------------------------
   
   
   !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
   !C.....FUNCION DE CALCULO DE TEMPERATURA AMBIENTE..................
   FUNCTION TEMAMB (TiTotal,PI,TAprom,DifTA,TiProm,TiCom,TiMax,TiMin)
      USE typePrecision
      REAL(nP) TiTotal,PI,TAprom,DifTA,TiProm,TiCom,TiMax,TiMin
      
      !C-----ESTA RUTINA DEVUELVE EL VALOR DE TiTotal AL PRIMER CICLO ------
      IF (TiTotal.GT.1440.)THEN
         IAITY=(TiTotal-360.)/1440.
         TiTotal=TiTotal-IAITY*1440.
      ENDIF
      !C-----------------------------------------------------------------
      IF (TiTotal.GE.0.AND.TiTotal.LT.TiMin)THEN
         AIT1=TiTotal+1440.
         Tamb=TAprom+DifTA*SIN(PI*((1440.-(TiMax-TiMin))/2.-(AIT1-TiMax))/TiCom)
      END IF
      IF (TiTotal.GE.TiMin.AND.TiTotal.LE.TiMax)THEN
         Tamb=TAprom+DifTA*SIN(PI*(TiTotal-TiProm)/(TiMax-TiMin))
      END IF
      IF (TiTotal.GT.TiMax)THEN
         AIT1=TiTotal
         Tamb=TAprom+DifTA*SIN(PI*((1440.-(TiMax-TiMin))/2.-(AIT1-TiMax))/TiCom)
      END IF
      TEMAMB=Tamb
   END FUNCTION TEMAMB
   !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
   
   !C.....FUNCION DE CALCULO DE IRRADIANCIA TECHO..................
   FUNCTION IRRATE (IrrMaxTe,PI,TiSM,TiSS,TiSol)
      USE typePrecision
      REAL(nP) IrrMaxTe,PI,TiSM,TiSS,TiSol
      IRRATE=IrrMaxTe*SIN(PI*(TiSM-TiSS)/TiSol)
      IF(IRRATE.LT.0.)IRRATE=0.
   END FUNCTION IRRATE
   !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
   
   !C.....FUNCION DE CALCULO DE IRRADIANCIA PAREDES................
   FUNCTION IRRAPA (TiSol,TMAXS,TiSS,TiLM,TMAX,QMAXV,PI,TiSM)
      USE typePrecision
      REAL(nP) TiSol,TMAXS,TiSS,TiLM,TMAX,QMAXV,PI,TiSM
      AKA=TiSol/(2.*(TMAXS-TiSS))
      IF(TiLM.GT.TMAX)THEN
         IRRAPA=QMAXV*SIN((PI*0.5)*(TiSol-(TiSM-TiSS))/(TiSol-(TMAXS-TiSS)))
      ELSE
         IRRAPA=QMAXV*SIN(PI*AKA*(TiSM-TiSS)/TiSol)
      END IF
      IF(IRRAPA.LT.0.)IRRAPA=0.
   END FUNCTION IRRAPA
   
END MODULE adapt

!     Last change:  EF    3 Feb 2005   11:43 am

MODULE UVWT
   !Calcula los coeficientes de la ecuacion discretizada
   !para las velocidades y la temperatura.
   
   USE typePrecision
   USE typeArreglo3D
   USE typeMalla
   
   TYPE(arreglo3D) :: u,v,w,t
   
   CONTAINS
   
   SUBROUTINE setUVWT()
      
      !Dimensiona din�micamente la velocidad y la temperatura (en memoria)
      
      ALLOCATE(u%value(ni,nj,nk))
      ALLOCATE(v%value(ni,nj,nk))
      ALLOCATE(w%value(ni,nj,nk))
      ALLOCATE(t%value(ni,nj,nk))
      
      u%title = "Velocidad U"
      v%title = "Velocidad V"
      w%title = "Velocidad W"
      t%title = "Temperatura"
      
      !Inicializa las velocidades y la temperatura
      u%value(:,:,:)  = 0.0D0
      v%value(:,:,:)  = 0.0D0
      w%value(:,:,:)  = 0.0D0
      T%value(:,:,:)  = 0.0D0
      
   END SUBROUTINE setUVWT
   
   SUBROUTINE initVelFace()
      
      USE Propiedades
      USE presion
      
      REAL(nP) :: uface,vface,wface,rhoface
      
      DO i = 1,objMalla%l2
         DO j = 1,objMalla%m2
            DO k = 1,objMalla%n2
               
               !Calcula el flujo masico solo en sectores donde haya aire
               IF (visc%value(i,j,k).LT.1.0D25.AND.visc%value(i+1,j,k).LT.1.0D25) THEN
                  
                  ! Direcci�n X
                  iDir = 1
                  uface    = InterLin(u%value(i,j,k),u%value(i+1,j,k),i,iDir)
                  rhoFace  = InterLin(rho%value(i,j,k),rho%value(i+1,j,k),i,iDir)
                  F1%value(i,j,k) = rhoface*uface*objMalla%areaX(j,k)
               else
                  
                  F1%value(i,j,k) = 0.00D0
                  
               END IF
               
               !Calcula el flujo masico solo en sectores donde haya aire
               IF (visc%value(i,j,k).LT.1.0D25.AND.visc%value(i,j+1,k).LT.1.0D25) THEN
                  
                  !Direcci�n Y
                  iDir = 2
                  vface    = InterLin(v%value(i,j,k),v%value(i,j+1,k),j,iDir)
                  rhoFace  = InterLin(rho%value(i,j,k),rho%value(i,j+1,k),j,iDir)
                  F2%value(i,j,k) = rhoface*vface*objMalla%areaY(i,k)
               else
                  
                  F2%value(i,j,k) = 0.00D0
                  
               END IF
               
               !Calcula el flujo masico solo en sectores donde haya aire
               IF (visc%value(i,j,k).LT.1.0D25.AND.visc%value(i,j,k+1).LT.1.0D25) THEN
                  
                  !Direcci�n Z
                  iDir = 3
                  wface    = InterLin(w%value(i,j,k),w%value(i,j,k+1),k,iDir)
                  rhoFace  = InterLin(rho%value(i,j,k),rho%value(i,j,k+1),k,iDir)
                  F3%value(i,j,k) = rhoface*wface*objMalla%areaZ(i,j)
               else
                  
                  F3%value(i,j,k) = 0.00D0
                  
               END IF
               
            END DO
         END DO
      END DO
      
      
   END SUBROUTINE initVelFace
   
   SUBROUTINE calUVW(timer,DT)
      
      USE Coeficientes
      USE typeBorde
      USE solvers
      USE Propiedades
      USE adapt
      USE presion
      
      REAL(nP)        :: DT,timer
      
      INTENT(IN)      ::  timer,DT
      
      ! Dimensiona Dinamicamente el arreglo correspondiente a el coeficiente de difusion gamX
      ALLOCATE(gamX%value(ni,nj,nk))
      
      !Dimensiona Dinamicamente los arreglos correspondientes a los flujo de la variable f(:,:,:,1)
      ALLOCATE(FLX%CI1(nj,nk),FLX%PI1(nj,nk))
      ALLOCATE(FLX%CL1(nj,nk),FLX%PL1(nj,nk))
      
      ALLOCATE(FLX%CJ1(ni,nk),FLX%PJ1(ni,nk))
      ALLOCATE(FLX%CM1(ni,nk),FLX%PM1(ni,nk))
      
      ALLOCATE(FLX%CK1(ni,nj),FLX%PK1(ni,nj))
      ALLOCATE(FLX%CN1(ni,nj),FLX%PN1(ni,nj))
      
      !Inicializa los arreglos correspondientes los flujos.
      FLX%CI1 = 0.0D0
      FLX%PI1 = 0.0D0
      FLX%CL1 = 0.0D0
      FLX%PL1 = 0.0D0
      
      FLX%CJ1 = 0.0D0
      FLX%PJ1 = 0.0D0
      FLX%CM1 = 0.0D0
      FLX%PM1 = 0.0D0
      
      FLX%CK1 = 0.0D0
      FLX%PK1 = 0.0D0
      FLX%CN1 = 0.0D0
      FLX%PN1 = 0.0D0
      
      CALL initcoef()
      nf = 0
      CALL phi(nf,TIMER,u,v,w,t)
      CALL getCoefDominio()
      
      ! C�lculo de los coeficientes para la ecuaci�n de cantidad de movimiento
      ! Cantidad de variables a resolver ==> velocidades U, V, y W
      
      ! Construcci�n del lazo para todas las componentes de velocidad
      DO nf = 1,3
         SELECT CASE (nf)
         CASE (1) ! Velocidad U
            
            ! Inicializa todos los coeficientes de los bordes a cero
            CALL initCoefBordes()
            ! Condiciones de contorno que varian por iteraci�n
            !CALL bound(u,nf,iter)
            ! Calculo de los Coeficientes de Difusion
            CALL phi(nf,TIMER,u,v,w,t)
            ! Discretizacion de la Ecuacion de Cantidad de Movimiento
            CALL getCoefVol(u,DT)       !Coeficientes volumetricos
            CALL getCoefBordesUVW(u,nf) !Coeficientes en los bordes
            ! Aplica el esquema de Solucion (Alto orden)
            CALL scheme(u,nf)
            ! Resuelve el sistema de Ecuaciones
            CALL solveEq(u)
            ! Extrapola la Velocidad U hacia los bordes y las esquinas
            !CALL getBondaryValuesAndFluxs(u,nf)
            CALL setEdgesValues(u)
            
         CASE (2) ! Velocidad V
            
            ! Inicializa todos los coeficientes de los bordes a cero
            CALL initCoefBordes()
            ! Condiciones de contorno que varian por iteraci�n
            !CALL bound(v,nf,iter)
            ! Calculo de los Coeficientes de Difusion
            CALL phi(nf,TIMER,u,v,w,t)
            ! Discretizacion de la Ecuacion de Cantidad de Movimiento
            !CALL getCoefDominio()
            CALL getCoefVol(v,DT)        !Coeficientes volumetricos
            CALL getCoefBordesUVW(v,nf)  !Coeficientes en los borde
            ! Aplica el esquema de Solucion (Alto orden)
            CALL scheme(v,nf)
            ! Resuelve el sistema de Ecuaciones
            CALL solveEq(v)
            ! Extrapola la Velocidad V hacia los bordes y las esquinas
            !CALL getBondaryValuesAndFluxs(v,nf)
            CALL setEdgesValues(v)
            
         CASE (3) ! Velocidad W
            
            ! Inicializa todos los coeficientes a cero
            CALL initCoefBordes()
            ! Condiciones de contorno que varian por iteraci�n
            !CALL bound(w,nf,iter)
            ! Calculo de los Coeficientes de Difusion
            CALL phi(nf,TIMER,u,v,w,t)
            ! Discretizacion de la Ecuacion de Cantidad de Movimiento
            !CALL getCoefDominio()
            CALL getCoefVol(w,DT)        !Coeficientes volumetricos
            CALL getCoefBordesUVW(w,nf)  !Coeficientes en los borde
            ! Aplica el esquema de Solucion (Alto orden)
            CALL scheme(w,nf)
            ! Resuelve el sistema de Ecuaciones
            CALL solveEq(w)
            ! Extrapola la Velocidad W hacia los bordes y las esquinas
            !CALL getBondaryValuesAndFluxs(w,nf)
            CALL setEdgesValues(w)
            
         END SELECT
      END DO
      
      IF (lcalp) THEN
         CALL calP(u,v,w)
      ELSE
         CALL initVelFace()
      END IF
      
      call outflowBC()
      
      !Libera la memoria donde se encontraba almacenada la conductividad.
      DEALLOCATE(gamX%value)
      
      !Libera la memoria en la que se encontraban almecenados los arreglos de la variable f(:,:,:,1)
      DEALLOCATE(FLX%CI1,FLX%PI1)
      DEALLOCATE(FLX%CL1,FLX%PL1)
      
      DEALLOCATE(FLX%CJ1,FLX%PJ1)
      DEALLOCATE(FLX%CM1,FLX%PM1)
      
      DEALLOCATE(FLX%CK1,FLX%PK1)
      DEALLOCATE(FLX%CN1,FLX%PN1)
      
   END SUBROUTINE calUVW
   
   SUBROUTINE calT(DT,timer)
      
      USE Coeficientes
      USE typeBorde
      USE solvers
      USE Propiedades
      USE adapt
      USE presion
      
      INTEGER  :: nf
      
      REAL(nP) :: DT, timer
      
      INTENT(IN) :: DT,timer
      
      ! Dimensiona Dinamicamente el arreglo correspondiente a el coeficiente de difusion gamX
      ALLOCATE(gamX%value(ni,nj,nk))
      
      ! Dimensiona Dinamicamente los arreglos correspondientes a los flujo de la variable f(:,:,:,1)
      ! Dimensiona Dinamicamente los arreglos correspondientes a la condicion de borde
      ALLOCATE(FLX%CI1(nj,nk),FLX%PI1(nj,nk))
      ALLOCATE(FLX%CL1(nj,nk),FLX%PL1(nj,nk))
      
      ALLOCATE(FLX%CJ1(ni,nk),FLX%PJ1(ni,nk))
      ALLOCATE(FLX%CM1(ni,nk),FLX%PM1(ni,nk))
      
      ALLOCATE(FLX%CK1(ni,nj),FLX%PK1(ni,nj))
      ALLOCATE(FLX%CN1(ni,nj),FLX%PN1(ni,nj))
      
      ALLOCATE(kbc%I1(nj,nk),kbc%J1(ni,nk))
      ALLOCATE(kbc%K1(ni,nj),kbc%L1(nj,nk))
      ALLOCATE(kbc%M1(ni,nk),kbc%N1(ni,nj))
      
      !Inicializa por defecto la condicion de borde como conocida
      kbc%I1=1
      kbc%L1=1
      
      kbc%J1=1
      kbc%M1=1
      
      kbc%K1=1
      kbc%N1=1
      
      !Inicializa los arrglos de flujos correspondiente a la variable f(:,:,:,1)
      FLX%CI1 = 0.0
      FLX%PI1 = 0.0
      FLX%CL1 = 0.0
      FLX%PL1 = 0.0
      
      FLX%CJ1 = 0.0
      FLX%PJ1 = 0.0
      FLX%CM1 = 0.0
      FLX%PM1 = 0.0
      
      FLX%CK1 = 0.0
      FLX%PK1 = 0.0
      FLX%CN1 = 0.0
      FLX%PN1 = 0.0
      
      ! C�lculo de la Temperatura (ecuaci�n de energ�a)
      nf = 5
      ! Inicializa todos los coeficientes a cero
      CALL initCoef()
      ! Condiciones de contorno que varian por iteraci�n
      !CALL bound(T,nf,iter)
      ! Calculo de los Coeficientes de Difusion
      CALL phi(nf,timer,u,v,w,t)
      ! Discretizacion de la Ecuacion
      CALL getCoefDominio()
      CALL getCoefVol(T,DT)  !Coeficientes volumetricos
      CALL getCoefBounds(T)     !Coeficientes en los borde
      ! Aplica el esquema de Solucion (Alto orden)
      CALL scheme(T,nf)
      ! Resuelve el sistema de Ecuaciones
      CALL solveEq(T)
      ! Extrapola la temperatura hacia los bordes y las esquinas
      CALL getBondaryValuesAndFluxs(T,nf)
      CALL setEdgesValues(T)
      
      DEALLOCATE(gamX%value)
      
      DEALLOCATE(FLX%CI1,FLX%PI1)
      DEALLOCATE(FLX%CL1,FLX%PL1)
      
      DEALLOCATE(FLX%CJ1,FLX%PJ1)
      DEALLOCATE(FLX%CM1,FLX%PM1)
      
      DEALLOCATE(FLX%CK1,FLX%PK1)
      DEALLOCATE(FLX%CN1,FLX%PN1)
      
      DEALLOCATE(kbc%I1,kbc%J1)
      DEALLOCATE(kbc%K1,kbc%L1)
      DEALLOCATE(kbc%M1,kbc%N1)
      
   END SUBROUTINE calT
   
   SUBROUTINE getCoefVol(f,DT)
      
      USE Coeficientes
      USE Propiedades
      
      REAL(nP) :: DT,APT
      
      TYPE(arreglo3D)   :: f
      
      INTENT (IN)       :: DT
      
      DO i = 2,objMalla%l2
         DO j = 2,objMalla%m2
            DO k = 2,objMalla%n2
               
               ! C�lculo del t�rmino 1/ap
               APT = rho%value(i,j,k)/DT
               
               ! C�lculo del t�rmino fuente
               con%value(i,j,k) = (con%value(i,j,k) + APT*f%value(i,j,k))*objMalla%vol(i,j,k)
               
               ! C�lculo del t�rmino ap
               ap%value(i,j,k)  = (APT-ap%value(i,j,k))*objMalla%vol(i,j,k)
               
            END DO
         END DO
      END DO
      
   END SUBROUTINE getCoefVol
   
   SUBROUTINE getCoefDominio()
      
      USE Coeficientes
      USE Propiedades
      USE presion
      
      REAL(nP)     :: diff,gam,flow
      
      DO i = 2,objMalla%l2
         DO j = 2,objMalla%m2
            DO k = 2,objMalla%n2
               
               ! C�lculo de los coeficentes de dominio de los nodos vecinos en la direcci�n X,
               ! excluy�ndo los bordes
               IF (i.LE.objMalla%l2-1) THEN
                  flow = F1%value(i,j,k)
                  gam  = InterArmonica(GamX%value(i,j,k),GamX%value(i+1,j,k),objMalla%dx(i),objMalla%dx(i+1))
                  diff = objMalla%areaX(j,k)*gam
                  aim%value(i+1,j,k) = diff + MAX(flow,0.0D0)
                  aip%value(i,j,k) = aim%value(i+1,j,k) - flow
               ENDIF
               
               ! C�lculo de los coeficientes de dominio de los nodos vecinos en la direcci�n Y
               ! excluy�ndo los bordes.
               IF (j.LE.objMalla%m2-1) THEN
                  flow = F2%value(i,j,k)
                  gam  = InterArmonica(GamX%value(i,j,k),GamX%value(i,j+1,k),objMalla%dy(j),objMalla%dy(j+1))
                  diff = objMalla%areaY(i,k)*gam
                  ajm%value(i,j+1,k) = diff + MAX(flow,0.0D0)
                  ajp%value(i,j,k) = ajm%value(i,j+1,k) - flow
               ENDIF
               
               ! C�lculo de los coeficientes de dominio de los nodos vecinos en la direcci�n Z
               ! excluy�ndo los bordes.
               IF (k.LE.objMalla%n2-1) THEN
                  flow = F3%value(i,j,k)
                  gam  = InterArmonica(GamX%value(i,j,k),GamX%value(i,j,k+1),objMalla%dz(k),objMalla%dz(k+1))
                  diff = objMalla%areaZ(i,j)*gam
                  akm%value(i,j,k+1) = diff + MAX(flow,0.0D0)
                  akp%value(i,j,k) = akm%value(i,j,k+1) - flow
               ENDIF
               
            END DO
         END DO
      END DO
      
   END SUBROUTINE getCoefDominio
   
   SUBROUTINE getCoefBounds(f)
      
      USE typeBorde
      USE Coeficientes
      USE Propiedades
      USE presion
      
      REAL(nP)     :: diff,temp,flow
      
      TYPE(arreglo3D)      :: f
      
      ! Direccion X
      DO j = 2,objMalla%m2
         DO k = 2,objMalla%n2
            
            !C�lculo de coeficientes de los nodos vecinos en en borde izquierdo
            ! Borde I = 1
            Diff = 2.0D0 * GamX%value(2,j,k)/(objMalla%dx(2))+SMALL
            flow = F1%value(1,j,k)
            aim%value(2,j,k) = Diff + MAX(flow,0.0D0)
            aip%value(1,j,k) = aim%value(2,j,k) - flow
            aim%value(2,j,k) = aim%value(2,j,k) * objMalla%areaX(j,k)
            aim%value(1,j,k) = 0.0D0
            
            !
            IF (kbc%I1(j,k).EQ.1) THEN
               con%value(2,j,k) = con%value(2,j,k)+aim%value(2,j,k)*f%value(1,j,k)
            ELSE
               ap%value(1,j,k)  = aip%value(1,j,k) - FLX%PI1(j,k)
               con%value(1,j,k) = FLX%CI1(j,k)
               temp             = aim%value(2,j,k)/ap%value(1,j,k)
               ap%value(2,j,k)  = ap%value(2,j,k)-aip%value(1,j,k)*temp
               aip%value (2,j,k)= aip%value(2,j,k)-aim%value(1,j,k)*temp
               con%value(2,j,k) = con%value(2,j,k)+con%value(1,j,k)*temp
            ENDIF
            ap%value(2,j,k)  = ap%value(2,j,k)+aim%value(2,j,k)
            aim%value(2,j,k) = 0.0D0
            
            ! Borde I = L1
            Diff = 2.0D0 * GamX%value(objMalla%l2,j,k)/(objMalla%dx(objMalla%l2))+SMALL
            flow = F1%value(objMalla%l2,j,k)
            aip%value(objMalla%l2,j,k) = Diff + MAX(-flow,0.0D0)
            aim%value(objMalla%l1,j,k) = aip%value(objMalla%l2,j,k) + flow
            aip%value(objMalla%l2,j,k) = aip%value(objMalla%l2,j,k)*objMalla%areaX(j,k)
            aip%value(objMalla%l1,j,k)  = 0.0D0
            
            IF (kbc%L1(j,k).EQ.1) THEN
               con%value(objMalla%l2,j,k) = con%value(objMalla%l2,j,k)+aip%value(objMalla%l2,j,k)*f%value(objMalla%l1,j,k)
            ELSE
               ap%value(objMalla%l1,j,k)  = aim%value(objMalla%l1,j,k)-FLX%PL1(j,k)
               con%value(objMalla%l1,j,k) = FLX%CL1(j,k)
               temp                       = aip%value(objMalla%l2,j,k)/ap%value(objMalla%l1,j,k)
               ap%value(objMalla%l2,j,k)  = ap%value(objMalla%l2,j,k)-aim%value(objMalla%l1,j,k)*temp
               aim%value(objMalla%l2,j,k) = aim%value(objMalla%l2,j,k)-aip%value(objMalla%l1,j,k)*temp
               con%value(objMalla%l2,j,k) = con%value(objMalla%l2,j,k)+con%value(objMalla%l1,j,k)*temp
            ENDIF
            ap%value(objMalla%l2,j,k)  = ap%value(objMalla%l2,j,k)+aip%value(objMalla%l2,j,k)
            aip%value(objMalla%l2,j,k) = 0.0D0
         END DO
      END DO
      
      ! Direccion Y
      DO i = 2,objMalla%l2
         DO k = 2,objMalla%n2
            
            ! Borde J = 1
            Diff = 2.0D0 * GamX%value(i,2,k)/(objMalla%dy(2))+SMALL
            flow = F2%value(i,1,k)
            ajm%value(i,2,k) = Diff + MAX(flow,0.0D0)
            ajp%value(i,1,k) = ajm%value(i,2,k) - flow
            ajm%value(i,2,k) = ajm%value(i,2,k)*objMalla%areaY(i,k)
            ajm%value(i,1,k) = 0.0D0
            
            IF (kbc%J1(i,k).EQ.1) THEN
               con%value(i,2,k) = con%value(i,2,k)+ajm%value(i,2,k)*f%value(i,1,k)
            ELSE
               ap%value(i,1,k)  = ajp%value(i,1,k) - FLX%PJ1(i,k)
               con%value(i,1,k) = FLX%CJ1(i,k)
               temp             = ajm%value (i,2,k)/ap%value(i,1,k)
               ap%value(i,2,k)  = ap%value(i,2,k)-ajp%value (i,1,k)*temp
               ajp%value(i,2,k) = ajp%value (i,2,k)-ajm%value (i,1,k)*temp
               con%value(i,2,k) = con%value(i,2,k)+con%value(i,1,k)*temp
            ENDIF
            ap%value(i,2,k)  = ap%value(i,2,k)+ajm%value(i,2,k)
            ajm%value(i,2,k) = 0.0D0
            
            ! Borde J = M1
            Diff = 2.0D0 * GamX%value(i,objMalla%m2,k)/(objMalla%dy(objMalla%m2))+small
            flow = F2%value(i,objMalla%m2,k)
            ajp%value(i,objMalla%m2,k) = Diff + MAX(-flow,0.0D0)
            ajm%value(i,objMalla%m1,k) = ajp%value(i,objMalla%m2,k) + flow
            ajp%value(i,objMalla%m2,k) = ajp%value(i,objMalla%m2,k)*objMalla%areaY(i,k)
            ajp%value(i,objMalla%m1,k) = 0.0D0
            
            IF(kbc%M1(i,k).EQ.1) THEN
               con%value(i,objMalla%m2,k) = con%value(i,objMalla%m2,k)+ajp%value(i,objMalla%m2,k)*f%value(i,objMalla%m1,k)
            ELSE
               ap%value(i,objMalla%m1,k)  = ajm%value(i,objMalla%m1,k)-FLX%PM1(i,k)
               con%value(i,objMalla%m1,k) = FLX%CM1(i,k)
               temp                       = ajp%value(i,objMalla%m2,k)/ap%value(i,objMalla%m1,k)
               ap%value(i,objMalla%m2,k)  = ap%value(i,objMalla%m2,k)-ajm%value(i,objMalla%m1,k)*temp
               ajm%value(i,objMalla%m2,k) = ajm%value(i,objMalla%m2,k)-ajp%value(i,objMalla%m1,k)*temp
               con%value(i,objMalla%m2,k) = con%value(i,objMalla%m2,k)+con%value(i,objMalla%m1,k)*temp
            ENDIF
            ap%value(i,objMalla%m2,k)  = ap%value(i,objMalla%m2,k)+ajp%value(i,objMalla%m2,k)
            ajp%value(i,objMalla%m2,k) = 0.0D0
            
         END DO
      END DO
      
      ! Direccion Z
      DO i = 2,objMalla%l2
         DO j = 2,objMalla%m2
            
            ! Borde Z = 1
            Diff = 2.0D0 * GamX%value(i,j,2)/(objMalla%dz(2))+SMALL
            flow = F3%value(i,j,1)
            akm%value(i,j,2) = Diff + MAX(flow,0.0D0)
            akp%value(i,j,1) = akm%value(i,j,2) - flow
            akm%value(i,j,2) = akm%value(i,j,2)*objMalla%areaZ(i,j)
            akm%value(i,j,1) = 0.0D0
            
            IF (kbc%K1(i,j).EQ.1) THEN
               con%value(i,j,2) = con%value(i,j,2)+akm%value(i,j,2)*f%value(i,j,1)
            ELSE
               ap%value(i,j,1)  = akp%value(i,j,1) - FLX%PK1(i,j)
               con%value(i,j,1) = FLX%CK1(i,j)
               temp             = akm%value (i,j,2)/ap%value(i,j,1)
               ap%value(i,j,2)  = ap%value(i,j,2)-akp%value (i,j,1)*temp
               akp%value(i,j,2) = akp%value (i,j,2)-akm%value (i,j,1)*temp
               con%value(i,j,2) = con%value(i,j,2)+con%value(i,j,1)*temp
            ENDIF
            ap%value(i,j,2)  = ap%value(i,j,2)+akm%value(i,j,2)
            akm%value(i,j,2) = 0.0D0
            
            ! Borde Z = N1
            Diff = 2.0D0 * GamX%value(i,j,objMalla%n2)/(objMalla%dz(objMalla%n2))+small
            flow = F3%value(i,j,objMalla%n2)
            akp%value(i,j,objMalla%n2) = Diff + MAX(-flow,0.0D0)
            akm%value(i,j,objMalla%n1) = akp%value(i,j,objMalla%n2) + flow
            akp%value(i,j,objMalla%n2) = akp%value(i,j,objMalla%n2)*objMalla%areaZ(i,j)
            akp%value(i,j,objMalla%n1) = 0.0D0
            
            IF(kbc%N1(i,j).EQ.1) THEN
               con%value(i,j,objMalla%n2)  = con%value(i,j,objMalla%n2)+akp%value(i,j,objMalla%n2)*f%value(i,j,objMalla%n1)
            ELSE
               ap%value(i,j,objMalla%n1)   = akm%value(i,j,objMalla%n1)-FLX%PN1(i,j)
               con%value(i,j,objMalla%n1)  = FLX%CN1(i,j)
               temp                        = akp%value(i,j,objMalla%n2)/ap%value(i,j,objMalla%n1)
               ap%value(i,j,objMalla%n2)   = ap%value(i,j,objMalla%n2)-akm%value(i,j,objMalla%n1)*temp
               akm%value(i,j,objMalla%n2)  = akm%value(i,j,objMalla%n2)-akp%value(i,j,objMalla%n1)*temp
               con%value(i,j,objMalla%n2)  = con%value(i,j,objMalla%n2)+con%value(i,j,objMalla%n1)*temp
            ENDIF
            ap%value(i,j,objMalla%n2)  = ap%value(i,j,objMalla%n2)+akp%value(i,j,objMalla%n2)
            akp%value(i,j,objMalla%n2) = 0.0D0
            
         END DO
      END DO
      
   END SUBROUTINE getCoefBounds
   
   SUBROUTINE getCoefBordesUVW(f,nf)
      
      USE typeBorde
      USE Coeficientes
      USE Propiedades
      USE presion
      
      REAL(nP)     :: diff,temp,flow
      INTEGER      :: nf
      INTENT (IN)  :: nf
      
      TYPE(arreglo3D)      :: f
      
      DO nDir = 1,3
         
         SELECT CASE (nDir)
            
         CASE (1)    ! Direccion X
            
            DO j = 2,objMalla%m2
               DO k = 2,objMalla%n2
                  
                  ! Borde I = 1
                  Diff = 2.0D0 * GamX%value(2,j,k)/(objMalla%dx(2))+SMALL
                  flow = F1%value(1,j,k)
                  aim%value(2,j,k) = Diff + MAX(flow,0.0D0)
                  aip%value(1,j,k) = aim%value(2,j,k) - flow
                  aim%value(2,j,k) = aim%value(2,j,k) * objMalla%areaX(j,k)
                  aim%value(1,j,k) = 0.0D0
                  
                  IF (get_BC_I1(j,k,nDir,nf).eq.1) THEN
                     con%value(2,j,k) = con%value(2,j,k)+aim%value(2,j,k)*f%value(1,j,k)
                  ELSE
                     ap%value(1,j,k)  = aip%value(1,j,k) - FLX%PI1(j,k)
                     con%value(1,j,k) = FLX%CI1(j,k)
                     temp             = aim%value(2,j,k)/ap%value(1,j,k)
                     ap%value(2,j,k)  = ap%value(2,j,k)-aip%value(1,j,k)*temp
                     aip%value (2,j,k)= aip%value(2,j,k)-aim%value(1,j,k)*temp
                     con%value(2,j,k) = con%value(2,j,k)+con%value(1,j,k)*temp
                  ENDIF
                  ap%value(2,j,k)  = ap%value(2,j,k)+aim%value(2,j,k)
                  aim%value(2,j,k) = 0.0D0
                  
                  ! Borde I = L1
                  Diff = 2.0D0 * GamX%value(objMalla%l2,j,k)/(objMalla%dx(objMalla%l2))+SMALL
                  flow = F1%value(objMalla%l2,j,k)
                  aip%value(objMalla%l2,j,k) = Diff + MAX(-flow,0.0D0)
                  aim%value(objMalla%l1,j,k) = aip%value(objMalla%l2,j,k) + flow
                  aip%value(objMalla%l2,j,k) = aip%value(objMalla%l2,j,k)*objMalla%areaX(j,k)
                  aip%value(objMalla%l1,j,k)  = 0.0D0
                  
                  IF (get_BC_L1(j,k,nDir,nf).eq.1) THEN
                     con%value(objMalla%l2,j,k) = con%value(objMalla%l2,j,k)+aip%value(objMalla%l2,j,k)*f%value(objMalla%l1,j,k)
                  ELSE
                     ap%value(objMalla%l1,j,k)  = aim%value(objMalla%l1,j,k)-FLX%PL1(j,k)
                     con%value(objMalla%l1,j,k) = FLX%CL1(j,k)
                     temp                       = aip%value(objMalla%l2,j,k)/ap%value(objMalla%l1,j,k)
                     ap%value(objMalla%l2,j,k)  = ap%value(objMalla%l2,j,k)-aim%value(objMalla%l1,j,k)*temp
                     aim%value(objMalla%l2,j,k) = aim%value(objMalla%l2,j,k)-aip%value(objMalla%l1,j,k)*temp
                     con%value(objMalla%l2,j,k) = con%value(objMalla%l2,j,k)+con%value(objMalla%l1,j,k)*temp
                  ENDIF
                  ap%value(objMalla%l2,j,k)  = ap%value(objMalla%l2,j,k)+aip%value(objMalla%l2,j,k)
                  aip%value(objMalla%l2,j,k) = 0.0D0
               END DO
            END DO
            
         CASE (2)   ! Direccion Y
            
            DO i = 2,objMalla%l2
               DO k = 2,objMalla%n2
                  
                  ! Borde J = 1
                  Diff = 2.0D0 * GamX%value(i,2,k)/(objMalla%dy(2))+SMALL
                  flow = F2%value(i,1,k)
                  ajm%value(i,2,k) = Diff + MAX(flow,0.0D0)
                  ajp%value(i,1,k) = ajm%value(i,2,k) - flow
                  ajm%value(i,2,k) = ajm%value(i,2,k)*objMalla%areaY(i,k)
                  ajm%value(i,1,k) = 0.0D0
                  
                  IF (get_BC_J1(i,k,nDir,nf).eq.1) THEN
                     con%value(i,2,k) = con%value(i,2,k)+ajm%value(i,2,k)*f%value(i,1,k)
                  ELSE
                     ap%value(i,1,k)  = ajp%value(i,1,k) - FLX%PJ1(i,k)
                     con%value(i,1,k) = FLX%CJ1(i,k)
                     temp             = ajm%value (i,2,k)/ap%value(i,1,k)
                     ap%value(i,2,k)  = ap%value(i,2,k)-ajp%value (i,1,k)*temp
                     ajp%value(i,2,k) = ajp%value (i,2,k)-ajm%value (i,1,k)*temp
                     con%value(i,2,k) = con%value(i,2,k)+con%value(i,1,k)*temp
                  ENDIF
                  ap%value(i,2,k)  = ap%value(i,2,k)+ajm%value(i,2,k)
                  ajm%value(i,2,k) = 0.0D0
                  
                  ! Borde J = M1
                  Diff = 2.0D0 * GamX%value(i,objMalla%m2,k)/(objMalla%dy(objMalla%m2))+small
                  flow = F2%value(i,objMalla%m2,k)
                  ajp%value(i,objMalla%m2,k) = Diff + MAX(-flow,0.0D0)
                  ajm%value(i,objMalla%m1,k) = ajp%value(i,objMalla%m2,k) + flow
                  ajp%value(i,objMalla%m2,k) = ajp%value(i,objMalla%m2,k)*objMalla%areaY(i,k)
                  ajp%value(i,objMalla%m1,k) = 0.0D0
                  
                  IF(get_BC_M1(i,k,nDir,nf).eq.1) THEN
                     con%value(i,objMalla%m2,k) = con%value(i,objMalla%m2,k)+ajp%value(i,objMalla%m2,k)*f%value(i,objMalla%m1,k)
                  ELSE
                     ap%value(i,objMalla%m1,k)  = ajm%value(i,objMalla%m1,k)-FLX%PM1(i,k)
                     con%value(i,objMalla%m1,k) = FLX%CM1(i,k)
                     temp                       = ajp%value(i,objMalla%m2,k)/ap%value(i,objMalla%m1,k)
                     ap%value(i,objMalla%m2,k)  = ap%value(i,objMalla%m2,k)-ajm%value(i,objMalla%m1,k)*temp
                     ajm%value(i,objMalla%m2,k) = ajm%value(i,objMalla%m2,k)-ajp%value(i,objMalla%m1,k)*temp
                     con%value(i,objMalla%m2,k) = con%value(i,objMalla%m2,k)+con%value(i,objMalla%m1,k)*temp
                  ENDIF
                  ap%value(i,objMalla%m2,k)  = ap%value(i,objMalla%m2,k)+ajp%value(i,objMalla%m2,k)
                  ajp%value(i,objMalla%m2,k) = 0.0D0
                  
               END DO
            END DO
            
         CASE (3)  ! Direccion Z
            
            DO i = 2,objMalla%l2
               DO j = 2,objMalla%m2
                  
                  ! Borde Z = 1
                  Diff = 2.0D0 * GamX%value(i,j,2)/(objMalla%dz(2))+SMALL
                  flow = F3%value(i,j,1)
                  akm%value(i,j,2) = Diff + MAX(flow,0.0D0)
                  akp%value(i,j,1) = akm%value(i,j,2) - flow
                  akm%value(i,j,2) = akm%value(i,j,2)*objMalla%areaZ(i,j)
                  akm%value(i,j,1) = 0.0D0
                  
                  IF (get_BC_K1(i,j,nDir,nf).eq.1) THEN
                     con%value(i,j,2) = con%value(i,j,2)+akm%value(i,j,2)*f%value(i,j,1)
                  ELSE
                     ap%value(i,j,1)  = akp%value(i,j,1) - FLX%PK1(i,j)
                     con%value(i,j,1) = FLX%CK1(i,j)
                     temp             = akm%value (i,j,2)/ap%value(i,j,1)
                     ap%value(i,j,2)  = ap%value(i,j,2)-akp%value (i,j,1)*temp
                     akp%value(i,j,2) = akp%value (i,j,2)-akm%value (i,j,1)*temp
                     con%value(i,j,2) = con%value(i,j,2)+con%value(i,j,1)*temp
                  ENDIF
                  ap%value(i,j,2)  = ap%value(i,j,2)+akm%value(i,j,2)
                  akm%value(i,j,2) = 0.0D0
                  
                  ! Borde K = N1
                  Diff = 2.0D0 * GamX%value(i,j,objMalla%n2)/(objMalla%dz(objMalla%n2))+small
                  flow = F3%value(i,j,objMalla%n2)
                  akp%value(i,j,objMalla%n2) = Diff + MAX(-flow,0.0D0)
                  akm%value(i,j,objMalla%n1) = akp%value(i,j,objMalla%n2) + flow
                  akp%value(i,j,objMalla%n2) = akp%value(i,j,objMalla%n2)*objMalla%areaZ(i,j)
                  akp%value(i,j,objMalla%n1) = 0.0
                  
                  IF(get_BC_N1(i,j,nDir,nf).eq.1) THEN
                     con%value(i,j,objMalla%n2)  = con%value(i,j,objMalla%n2)+akp%value(i,j,objMalla%n2)*f%value(i,j,objMalla%n1)
                  ELSE
                     ap%value(i,j,objMalla%n1)   = akm%value(i,j,objMalla%n1)-FLX%PN1(i,j)
                     con%value(i,j,objMalla%n1)  = FLX%CN1(i,j)
                     temp                        = akp%value(i,j,objMalla%n2)/ap%value(i,j,objMalla%n1)
                     ap%value(i,j,objMalla%n2)   = ap%value(i,j,objMalla%n2)-akm%value(i,j,objMalla%n1)*temp
                     akm%value(i,j,objMalla%n2)  = akm%value(i,j,objMalla%n2)-akp%value(i,j,objMalla%n1)*temp
                     con%value(i,j,objMalla%n2)  = con%value(i,j,objMalla%n2)+con%value(i,j,objMalla%n1)*temp
                  ENDIF
                  ap%value(i,j,objMalla%n2)  = ap%value(i,j,objMalla%n2)+akp%value(i,j,objMalla%n2)
                  akp%value(i,j,objMalla%n2) = 0.0D0
                  
               END DO
            END DO
            
         END SELECT
         
      END DO
      
   END SUBROUTINE getCoefBordesUVW
   
   
   SUBROUTINE scheme(f,nf)
      
      USE Coeficientes
      USE LeerData
      USE presion
      
      INTEGER           :: nf
      TYPE(arreglo3D)   :: f
      
      REAL(nP)          :: rlx, anb, ainr
      
      INTENT(IN)        :: nf
      
      ! Subrelajaci�n
      rlx = (1.0D0 - relax%value(nf))/relax%value(nf)
      
      DO i=2,objMalla%l2
         DO j=2,objMalla%m2
            DO k=2,objMalla%n2
               
               anb  = aip%value(i,j,k) + aim%value(i,j,k) + ajp%value(i,j,k) + &
               ajm%value(i,j,k) + akp%value(i,j,k) + akm%value(i,j,k)
               
               ainr = anb * rlx
               
               IF (lcalP) then
                  IF (NF.EQ.1) dup%value(i,j,k)  = 1.0D0/(ap%value(i,j,k)+anb)
                  IF (NF.EQ.2) dvp%value(i,j,k)  = 1.0D0/(ap%value(i,j,k)+anb)
                  IF (NF.EQ.3) dwp%value(i,j,k)  = 1.0D0/(ap%value(i,j,k)+anb)
                  
                  IF (nf.eq.1) con%value(i,j,k) = con%value(i,j,k) - objMalla%vol(i,j,k)*dpx%value(i,j,k)
                  IF (nf.eq.2) con%value(i,j,k) = con%value(i,j,k) - objMalla%vol(i,j,k)*dpy%value(i,j,k)
                  IF (nf.eq.3) con%value(i,j,k) = con%value(i,j,k) - objMalla%vol(i,j,k)*dpz%value(i,j,k)
               END IF
               
               ap%value(i,j,k)   = ap%value(i,j,k)  + anb + ainr
               
               con%value(i,j,k)  = con%value(i,j,k) + ainr * f%value(i,j,k)
               
            END DO
         END DO
      END DO
      
   END SUBROUTINE scheme
   
   SUBROUTINE getBondaryValuesAndFluxs(f,nf)
      
      USE Coeficientes
      USE typeBorde
      
      INTEGER           :: nf
      TYPE(arreglo3D)   :: f
      
      REAL(nP) :: temp
      
      INTENT(IN)     :: nf
      
      ! Bordes en X
      DO j = 2,objMalla%m2
         DO k = 2,objMalla%n2
            
            i = 1
            temp = aim%value(i,j,k) * (f%value(i+2,j,k) - f%value(i+1,j,k))
            IF (kbc%I1(j,k).EQ.2) THEN
               f%value(i,j,k) = (aip%value(i,j,k) * f%value(i+1,j,k) - temp + con%value(i,j,k)) / ap%value(i,j,k)
               !else
               !IF (kbc_I1.EQ."OUTFLOW") T%value(i,j,k) = T%value(i+1,j,k)
            END IF
            !           flux%I1(j,k,nf) = aip%value(i,j,k) * (f%value(i,j,k)-f%value(i+1,j,k)) + temp
            !            IF (tipo%value(i,j,k).EQ.-1) t%value(i,j,k) = t%value(i+1,j,k)
            
            i = objMalla%l1
            temp = aip%value(i,j,k) * (f%value(i-2,j,k) - f%value(i-1,j,k))
            IF (kbc%L1(j,k).EQ.2) THEN
               f%value(i,j,k) = (aim%value(i,j,k) * f%value(i-1,j,k) - temp + con%value(i,j,k)) / ap%value(i,j,k)
               !else
               !IF (kbc_L1.eq."OUTFLOW") T%value(i,j,k) = T%value(i-1,j,k)
            END IF
            !           flux%L1(j,k,nf) = aim%value(i,j,k) * (f%value(i,j,k)-f%value(i-1,j,k)) + temp
            !            IF (tipo%value(i,j,k).EQ.-1) t%value(i,j,k) = t%value(i-1,j,k)
            
         END DO
      END DO
      
      ! Bordes en Y
      DO i = 2,objMalla%l2
         DO k = 2,objMalla%n2
            
            j = 1
            temp = ajm%value(i,j,k) * (f%value(i,j+2,k) - f%value(i,j+1,k))
            IF (kbc%J1(i,k).EQ.2) THEN
               f%value(i,j,k) = (ajp%value(i,j,k) * f%value(i,j+1,k) - temp + con%value(i,j,k)) / ap%value(i,j,k)
               !else
               !IF (kbc_J1.EQ."OUTFLOW") T%value(i,j,k) = T%value(i,j+1,k)
            END IF
            !           flux%J1(i,k,nf) = ajp%value(i,j,k) * (f%value(i,j,k)-f%value(i,j+1,k)) + temp
            !            IF (tipo%value(i,j,k).EQ.-1) t%value(i,j,k) = t%value(i,j+1,k)
            
            j = objMalla%m1
            temp = ajp%value(i,j,k) * (f%value(i,j-2,k) - f%value(i,j-1,k))
            IF (kbc%M1(i,k).EQ.2) THEN
               f%value(i,j,k) = (ajm%value(i,j,k) * f%value(i,j-1,k) - temp + con%value(i,j,k)) / ap%value(i,j,k)
               !else
               !IF (kbc_M1.EQ."OUTFLOW") T%value(i,j,k) = T%value(i,j-1,k)
            END IF
            !           flux%M1(i,k,nf) = ajm%value(i,j,k) * (f%value(i,j,k)-f%value(i,j-1,k)) + temp
            !            IF (tipo%value(i,j,k).EQ.-1) t%value(i,j,k) = t%value(i,j-1,k)
            
         END DO
      END DO
      
      ! Bordes en Z
      DO i = 2,objMalla%l2
         DO j = 2,objMalla%m2
            
            k = 1
            temp = akm%value(i,j,k) * (f%value(i,j,k+2) - f%value(i,j,k+1))
            IF (kbc%K1(i,j).EQ.2) THEN
               f%value(i,j,k) = (akp%value(i,j,k) * f%value(i,j,k+1) - temp + con%value(i,j,k)) / ap%value(i,j,k)
               !else
               !IF (kbc_K1.eq."OUTFLOW") T%value(i,j,k) = T%value(i,j,k+1)
            END IF
            !            IF (tipo%value(i,j,k).EQ.-1) t%value(i,j,k) = t%value(i,j,k+1)
            
            !           flux%K1(i,j,nf) = akp%value(i,j,k) * (f%value(i,j,k)-f%value(i,j,k+1)) + temp
            
            k = objMalla%n1
            temp = akp%value(i,j,k) * (f%value(i,j,k-2) - f%value(i,j,k-1))
            IF (kbc%N1(i,j).EQ.2) THEN
               f%value(i,j,k) = (akm%value(i,j,k) * f%value(i,j,k-1) - temp + con%value(i,j,k)) / ap%value(i,j,k)
               !ELSE
               !IF (kbc_N1.eq."OUTFLOW") T%value(i,j,k) = T%value(i,j,k-1)
            END IF
            !            IF (tipo%value(i,j,k).EQ.-1) t%value(i,j,k) = t%value(i,j,k-1)
            !           flux%N1(i,j,nf) = akm%value(i,j,k) * (f%value(i,j,k)-f%value(i,j,k-1)) + temp
            
         END DO
      END DO
      
   END SUBROUTINE getBondaryValuesAndFluxs
   
   SUBROUTINE outflowBC()
      
      USE typePrecision
      USE typeMalla
      USE propiedades
      USE typeBorde
      USE presion
      
      ! Determina el flujo masa  que entra  y sale del  dominio
      ! Corrige el valor de la velocidad en el contorno OUTFLOW
      ! a fin de satisfacer el balance global de masa
      
      REAL(nP) :: flowIn, flowOut, factor
      
      factor   = 0.0D0
      flowIn   = 0.0D0
      flowOut  = 0.0D0
      
      ! Plano YZ
      DO j = 2, objMalla%m2
         DO k = 2, objMalla%n2
            i = 2 	       ! Cara E
            IF (bc_I1(j,k).eq.INFLOW) then
               flowIn   =  flowIn + objMalla%areaX(j,k) * rho%value(i,j,k) * u%value(i-1,j,k)
            else IF (bc_I1(j,k).eq.OUTFLOW) then
               flowOut  = flowOut + objMalla%areaX(j,k) * rho%value(i,j,k) * u%value(i,j,k)
            END IF
            i = objMalla%l1  ! Cara W
            IF (bc_L1(j,k).eq.INFLOW) then
               flowIn   =  flowIn + objMalla%areaX(j,k) * rho%value(i-1,j,k) * u%value(i,j,k)
            else IF (bc_L1(j,k).eq.OUTFLOW) then
               flowOut  = flowOut + objMalla%areaX(j,k) * rho%value(i-1,j,k) * u%value(i-1,j,k)
            END IF
         END DO
      END DO
      
      !	   ! Plano XZ
      DO i = 2, objMalla%l2
         DO k = 2, objMalla%n2
            j = 2             ! Cara S
            IF (bc_J1(i,k).eq.INFLOW) then
               flowIn   =  flowIn + objMalla%areaY(i,k) * rho%value(i,j,k) * v%value(i,j-1,k)
            else IF (bc_J1(i,k).eq.OUTFLOW) then
               flowOut  = flowOut + objMalla%areaY(i,k) * rho%value(i,j,k) * v%value(i,j,k)
            END IF
            j = objMalla%m1   ! Cara N
            IF (bc_M1(i,k).eq.INFLOW) then
               flowIn   =  flowIn + objMalla%areaY(i,k) * rho%value(i,j-1,k) * v%value(i,j,k)
            else IF (bc_M1(i,k).eq.OUTFLOW) then
               flowOut  = flowOut + objMalla%areaY(i,k) * rho%value(i,j-1,k) * v%value(i,j-1,k)
            END IF
         END DO
      END DO
      
      ! Plano XY
      DO i = 2, objMalla%l2
         DO j = 2, objMalla%m2
            k = 2             ! Cara B
            IF (bc_K1(i,j).eq.INFLOW) then
               flowIn   =  flowIn + objMalla%areaZ(i,j) * rho%value(i,j,k) * w%value(i,j,k-1)
            else IF (bc_K1(i,j).eq.OUTFLOW) then
               flowOut  = flowOut + objMalla%areaZ(i,j) * rho%value(i,j,k) * w%value(i,j,k)
            END IF
            k = objMalla%n1   ! Cara T
            IF (bc_N1(i,j).eq.INFLOW) then
               flowIn   =  flowIn + objMalla%areaZ(i,j) * rho%value(i,j,k-1) * w%value(i,j,k)
            else IF (bc_N1(i,j).eq.OUTFLOW) then
               flowOut  = flowOut + objMalla%areaZ(i,j) * rho%value(i,j,k-1) * w%value(i,j,k-1)
            END IF
         END DO
      END DO
      
      ! Corrige el valor de la velocidad y el flujo de masa a la salida de cada contorno.
      ! La correcci�n es realizada si la condici�n de contorno es igual a OUTPUT
      factor = flowIn / (flowOut+SMALL)
      
      ! Plano YZ
      DO j = 2, objMalla%m2
         DO k = 2, objMalla%n2
            i = 1
            IF (bc_I1(j,k).eq.OUTFLOW) then
               u%value(i,j,k)    = u%value(i+1,j,k) * factor
               v%value(i,j,k)    = v%value(i+1,j,k)
               w%value(i,j,k)    = w%value(i+1,j,k)
               F1%value(i,j,k)   = u%value(i,j,k) * objMalla%areaX(j,k) * rho%value(i+1,j,k)
            END IF
            
            i = objMalla%l1
            IF (bc_L1(j,k).eq.OUTFLOW) then
               u%value(i,j,k)    = u%value(i-1,j,k) * factor
               v%value(i,j,k)    = v%value(i-1,j,k)
               w%value(i,j,k)    = w%value(i-1,j,k)
               F1%value(i-1,j,k) = u%value(i,j,k) * objMalla%areaX(j,k) * rho%value(i-1,j,k)
            END IF
         END DO
      END DO
      
      !Plano XZ
      DO i = 2, objMalla%l2
         DO k = 2, objMalla%n2
            j = 1
            IF (bc_J1(i,k).eq.OUTFLOW) then
               u%value(i,j,k)    = u%value(i,j+1,k)
               v%value(i,j,k)    = v%value(i,j+1,k) * factor
               w%value(i,j,k)    = w%value(i,j+1,k)
               F2%value(i,j,k)   = v%value(i,j,k) * objMalla%areaY(i,k) * rho%value(i,j+1,k)
            END IF
            
            j = objMalla%m1
            IF (bc_M1(i,k).eq.OUTFLOW) then
               u%value(i,j,k)    = u%value(i,j-1,k)
               v%value(i,j,k)    = v%value(i,j-1,k) * factor
               w%value(i,j,k)    = w%value(i,j-1,k)
               F2%value(i,j-1,k) = v%value(i,j,k) * objMalla%areaY(i,k) * rho%value(i,j-1,k)
            END IF
            
         END DO
      END DO
      
      ! Plano XY
      DO i = 2, objMalla%l2
         DO j = 2, objMalla%m2
            k = 1
            IF (bc_K1(i,j).eq.OUTFLOW) then
               u%value(i,j,k)    = u%value(i,j,k+1)
               v%value(i,j,k)    = v%value(i,j,k+1)
               w%value(i,j,k)    = w%value(i,j,k+1) * factor
               F3%value(i,j,k)   = w%value(i,j,k) * objMalla%areaZ(i,j) * rho%value(i,j,k+1)
            END IF
            
            k = objMalla%n1
            IF (bc_N1(i,j).eq.OUTFLOW) then
               u%value(i,j,k) 	= u%value(i,j,k-1)
               v%value(i,j,k) 	= v%value(i,j,k-1)
               w%value(i,j,k) 	= w%value(i,j,k-1) * factor
               F3%value(i,j,k-1) = w%value(i,j,k) * objMalla%areaZ(i,j) * rho%value(i,j,k-1)
            END IF
         END DO
      END DO
      
   END SUBROUTINE outflowBC
   
   
END MODULE UVWT

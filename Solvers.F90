!     Last change:  EF   10 Mar 2005   11:10 am

MODULE solvers
   !Este modulo crea la matriz para la solucion de la
   !ecuacion discretizada, luego llama al metodo de
   !solucion utilizado para la resolucion de las variables
   
   USE typePrecision
   USE typeMalla
   USE LeerData
   
   TYPE coordLoc
   INTEGER       :: i,j,k
   INTEGER       :: nodesAct
   INTEGER       :: nodeVec(6)
   END TYPE coordLoc
   
   TYPE nodeLoc
   INTEGER  :: node
   END TYPE nodeLoc
   
   TYPE(coordLoc),   POINTER  :: getCoord(:)
   TYPE(nodeLoc),    POINTER  :: getNode(:,:,:)
   INTEGER, ALLOCATABLE       :: JA(:),IA(:)
   REAL(nP), POINTER          :: AA(:)
   
   INTEGER :: numNodes,nodesActTotal
   
   CONTAINS
   
   SUBROUTINE setupSolver()
      
      nodesActTotal = 0
      numNodes      = 0
      
      !Construye la estructura getCoord & getNode
      ALLOCATE (getCoord((ni-2)*(nj-2)*(nk-2)))
      ALLOCATE (getNode((ni-1),(nj-1),(nk-1)))
      CALL makeArr4Solve()
      
      !Construye los arreglos JA & IA
      ALLOCATE (JA(nodesActTotal))
      ALLOCATE (IA(numNodes + 1))
      ALLOCATE (AA(nodesActTotal))
      CALL getArr4Solve ()
      
   END SUBROUTINE setupSolver
   
   SUBROUTINE  makeArr4Solve()
      
      INTEGER              :: i,j,k,l2,m2,n2,node,nodesAct,nVecino
      
      numNodes = (ni-2)*(nj-2)*(nk-2)
      
      l2 = ni - 1
      m2 = nj - 1
      n2 = nk - 1
      
      ! Llena la estructura getCoord & getNode
      node = 0
      
      DO k = 2, n2
         DO j = 2, m2
            DO i = 2, l2
               
               node = node + 1
               
               getCoord(node)%i = i
               getCoord(node)%j = j
               getCoord(node)%k = k
               
               getNode(i,j,k)%node = node
               
            END DO
         END DO
      END DO
      
      nodesActTotal = 0
      
      DO node = 1, numNodes
         
         nodesAct = 0
         
         i = getCoord(node)%i
         j = getCoord(node)%j
         k = getCoord(node)%k
         
         DO nVecino = 1, 6
            
            ! nVecino = 1 -- Cara E
            ! nVecino = 2 -- Cara W
            ! nVecino = 3 -- Cara N
            ! nVecino = 4 -- Cara S
            ! nVecino = 5 -- Cara T
            ! nVecino = 6 -- Cara B
            
            SELECT CASE (nVecino)
               
            CASE (1)
               
               IF (i.EQ.l2) THEN
                  getCoord(node)%nodeVec(nVecino) = 0
               ELSE
                  getCoord(node)%nodeVec(nVecino) = getNode(i+1,j,k)%node
                  nodesAct = nodesAct + 1
               END IF
               
            CASE (2)
               
               IF (i.EQ.2) THEN
                  getCoord(node)%nodeVec(nVecino) = 0
               ELSE
                  getCoord(node)%nodeVec(nVecino) = getNode(i-1,j,k)%node
                  nodesAct = nodesAct + 1
               END IF
               
            CASE (3)
               
               IF (j.EQ.m2) THEN
                  getCoord(node)%nodeVec(nVecino) = 0
               ELSE
                  getCoord(node)%nodeVec(nVecino) = getNode(i,j+1,k)%node
                  nodesAct = nodesAct + 1
               END IF
               
            CASE (4)
               
               IF (j.EQ.2) THEN
                  getCoord(node)%nodeVec(nVecino) = 0
               ELSE
                  getCoord(node)%nodeVec(nVecino) = getNode(i,j-1,k)%node
                  nodesAct = nodesAct + 1
               END IF
               
            CASE (5)
               
               IF (k.EQ.n2) THEN
                  getCoord(node)%nodeVec(nVecino) = 0
               ELSE
                  getCoord(node)%nodeVec(nVecino) = getNode(i,j,k+1)%node
                  nodesAct = nodesAct + 1
               END IF
               
            CASE (6)
               
               IF (k.EQ.2) THEN
                  getCoord(node)%nodeVec(nVecino) = 0
               ELSE
                  getCoord(node)%nodeVec(nVecino) = getNode(i,j,k-1)%node
                  nodesAct = nodesAct + 1
               END IF
               
            END SELECT
            
            getCoord(node)%nodesAct = nodesAct + 1
            
         END DO
         
         nodesActTotal = nodesActTotal + getCoord(node)%nodesAct
         
      END DO
      
   END SUBROUTINE makeArr4Solve
   
   SUBROUTINE getArr4Solve ()
      
      INTEGER     :: n, nFil, nVecino, nPos
      
      ! Establece el tama�o de los arreglos JA & IA
      n = 0
      
      DO nFil = 1, numNodes
         
         n = n + 1
         
         JA(n) = nFil
         
         DO nVecino = 1, 6
            
            IF (getCoord(nFil)%nodeVec(nVecino).NE.0) THEN
               
               n = n + 1
               
               nPos = getCoord(nFil)%nodeVec(nVecino)
               
               ! nVecino = 1 -- Cara E
               ! nVecino = 2 -- Cara W
               ! nVecino = 3 -- Cara N
               ! nVecino = 4 -- Cara S
               ! nVecino = 5 -- Cara T
               ! nVecino = 6 -- Cara B
               
               JA(n) = nPos
               
            END IF
            
         END DO
         
      END DO
      
      ! No se encuentra optimizado, pero es entendible
      IA(1) = 1
      DO nFil = 1, numNodes+1
         IF (nFil.GE.2) THEN
            IA(nFil) = IA(nFil-1) + getCoord(nFil-1)%nodesAct
         END IF
      END DO
      
   END SUBROUTINE getArr4Solve
   
   SUBROUTINE solveEq(f)
      
      USE typeArreglo3D
      
      INTEGER 	  :: icon,ntimes
      
      TYPE(arreglo3D)  :: f
      
      ! indAlgSol = 1 --> TDMA+GS (Default)
      ! indAlgSol = 2 --> SSORSI
      ! indAlgSol = 3 --> AMG
      
      SELECT CASE (indAlgSol)
         
      CASE default
         n = 1
         ntimes = 10        !numero de barridos realizados
         DO WHILE (n.NE.ntimes)
            CALL blockCorrectionScheme(f,icon)
            CALL tdma3D(f)
            
            IF (icon.EQ.0) THEN
               n = n + 1
            ELSE
               n = ntimes
            END IF
         END DO
         
      CASE (2,3,4)
         
         CALL useItpack2(f)
         
      END SELECT
      
   END SUBROUTINE solveEq
   
   SUBROUTINE blockCorrectionScheme(f,icon)
      
      USE typeArreglo3D
      USE Coeficientes
      
      INTEGER         :: icon
      TYPE(arreglo3D) :: f
      
      REAL(nP) :: denom
      
      REAL(nP), ALLOCATABLE :: ptx(:),qtx(:)
      REAL(nP), ALLOCATABLE :: pty(:),qty(:)
      REAL(nP), ALLOCATABLE :: ptz(:),qtz(:)
      
      REAL(nP) :: bl,blp,blm,blc
      REAL(nP) :: res,term,rt(8)
      REAL(nP) :: BIG, tmp, SMALL1
      
      INTENT (IN OUT)  :: f,icon
      
      ALLOCATE(ptx(ni),qtx(ni))
      ALLOCATE(pty(nj),qty(nj))
      ALLOCATE(ptz(nk),qtz(nk))
      
      SMALL1= 1.0E-20
      BIG   = 1.0E10
      term  = 1.0E-8
      tmp   = 1.0E-10
      
      icon  = 1
      
      ! Direccion X
      
      ptx(1) = 0.0
      qtx(1) = 0.0
      
      DO i = 2,objMalla%l2
         
         bl  = SMALL1
         blp = 0.0
         blm = 0.0
         blc = 0.0
         
         DO k = 2,objMalla%n2
            DO j = 2,objMalla%m2
               
               IF (ap%value(i,j,k).LT.BIG) THEN
                  
                  bl = bl + ap%value(i,j,k)
                  
                  IF (ap%value(i,j+1,k).LT.BIG) bl  = bl  - ajp%value(i,j,k)
                  IF (ap%value(i,j-1,k).LT.BIG) bl  = bl  - ajm%value(i,j,k)
                  IF (ap%value(i,j,k+1).LT.BIG) bl  = bl  - akp%value(i,j,k)
                  IF (ap%value(i,j,k-1).LT.BIG) bl  = bl  - akm%value(i,j,k)
                  IF (ap%value(i+1,j,k).LT.BIG) blp = blp + aip%value(i,j,k)
                  IF (ap%value(i-1,j,k).LT.BIG) blm = blm + aim%value(i,j,k)
                  
                  ! Vector de Residuales
                  rt(1) = aip%value(i,j,k)*f%value(i+1,j,k)
                  rt(2) = aim%value(i,j,k)*f%value(i-1,j,k)
                  rt(3) = ajp%value(i,j,k)*f%value(i,j+1,k)
                  rt(4) = ajm%value(i,j,k)*f%value(i,j-1,k)
                  rt(5) = akp%value(i,j,k)*f%value(i,j,k+1)
                  rt(6) = akm%value(i,j,k)*f%value(i,j,k-1)
                  rt(7) = -ap%value(i,j,k)*f%value(i,j,k)
                  rt(8) = con%value(i,j,k)
                  
                  ! Residual
                  res  = SUM(rt)
                  
                  DO n = 1,8
                     term = MAX(term,ABS(rt(n)))
                  END DO
                  
                  IF (ABS(res/term).GT.tmp) THEN
                     icon = 0
                  END IF
                  
                  blc  = blc + res
                  
               END IF
               
            END DO
         END DO
         
         denom = bl - ptx(i-1) * blm
         
         IF (ABS(denom/bl).LT.SMALL1) THEN
            denom = BIG
         END IF
         
         ptx(i) = blp / denom
         qtx(i) = (blc+blm*qtx(i-1))/denom
         
      END DO
      
      bl = 0.0
      
      DO i = objMalla%l2,2,-1
         
         bl = bl * ptx(i) + qtx(i)
         
         DO k = 2,objMalla%n2
            DO j = 2,objMalla%m2
               IF (ap%value(i,j,k).LT.BIG) THEN
                  f%value(i,j,k) = f%value(i,j,k) + bl
               ENDIF
            END DO
         END DO
         
      END DO
      
      !WRITE(*,*) "Block Correction Scheme Performed in X"
      
      ! Direccion Y
      
      pty(1) = 0.0
      qty(1) = 0.0
      
      DO j = 2,objMalla%m2
         
         bl  = SMALL1
         blp = 0.0
         blm = 0.0
         blc = 0.0
         
         DO k = 2,objMalla%n2
            DO i = 2,objMalla%l2
               
               IF (ap%value(i,j,k).LT.BIG) THEN
                  
                  bl = bl + ap%value(i,j,k)
                  
                  IF (ap%value(i+1,j,k).LT.BIG) bl  = bl  - aip%value(i,j,k)
                  IF (ap%value(i-1,j,k).LT.BIG) bl  = bl  - aim%value(i,j,k)
                  IF (ap%value(i,j,k+1).LT.BIG) bl  = bl  - akp%value(i,j,k)
                  IF (ap%value(i,j,k-1).LT.BIG) bl  = bl  - akm%value(i,j,k)
                  IF (ap%value(i,j+1,k).LT.BIG) blp = blp + ajp%value(i,j,k)
                  IF (ap%value(i,j-1,k).LT.BIG) blm = blm + ajm%value(i,j,k)
                  
                  ! Vector de Residuales
                  rt(1) = aip%value(i,j,k)*f%value(i+1,j,k)
                  rt(2) = aim%value(i,j,k)*f%value(i-1,j,k)
                  rt(3) = ajp%value(i,j,k)*f%value(i,j+1,k)
                  rt(4) = ajm%value(i,j,k)*f%value(i,j-1,k)
                  rt(5) = akp%value(i,j,k)*f%value(i,j,k+1)
                  rt(6) = akm%value(i,j,k)*f%value(i,j,k-1)
                  rt(7) = -ap%value(i,j,k)*f%value(i,j,k)
                  rt(8) = con%value(i,j,k)
                  
                  ! Residual
                  res  = SUM(rt)
                  
                  DO n = 1,8
                     term = MAX(term,ABS(rt(n)))
                  END DO
                  
                  IF (ABS(res/term).GT.tmp) THEN
                     icon = 0
                  END IF
                  
                  blc  = blc + res
                  
               END IF
               
            END DO
         END DO
         
         denom = bl - pty(j-1) * blm
         
         IF (ABS(denom/bl).LT.SMALL1) THEN
            denom = BIG
         END IF
         
         pty(j) = blp / denom
         qty(j) = (blc+blm*qty(j-1))/denom
         
      END DO
      
      bl = 0.0
      
      DO j = objMalla%m2,2,-1
         
         bl = bl * pty(j) + qty(j)
         
         DO k = 2,objMalla%n2
            DO i = 2,objMalla%l2
               IF (ap%value(i,j,k).LT.BIG) THEN
                  f%value(i,j,k) = f%value(i,j,k) + bl
               ENDIF
            END DO
         END DO
         
      END DO
      
      !WRITE(*,*) "Block Correction Scheme Performed in Y"
      
      ! Direccion Z
      
      ptz(1) = 0.0
      qtz(1) = 0.0
      
      DO k = 2,objMalla%n2
         
         bl  = SMALL1
         blp = 0.0
         blm = 0.0
         blc = 0.0
         
         DO j = 2,objMalla%m2
            DO i = 2,objMalla%l2
               
               IF (ap%value(i,j,k).LT.BIG) THEN
                  
                  bl = bl + ap%value(i,j,k)
                  
                  IF (ap%value(i+1,j,k).LT.BIG) bl  = bl  - aip%value(i,j,k)
                  IF (ap%value(i-1,j,k).LT.BIG) bl  = bl  - aim%value(i,j,k)
                  IF (ap%value(i,j+1,k).LT.BIG) bl  = bl  - ajp%value(i,j,k)
                  IF (ap%value(i,j-1,k).LT.BIG) bl  = bl  - ajm%value(i,j,k)
                  IF (ap%value(i,j,k+1).LT.BIG) blp = blp + akp%value(i,j,k)
                  IF (ap%value(i,j,k-1).LT.BIG) blm = blm + akm%value(i,j,k)
                  
                  ! Vector de Residuales
                  rt(1) = aip%value(i,j,k)*f%value(i+1,j,k)
                  rt(2) = aim%value(i,j,k)*f%value(i-1,j,k)
                  rt(3) = ajp%value(i,j,k)*f%value(i,j+1,k)
                  rt(4) = ajm%value(i,j,k)*f%value(i,j-1,k)
                  rt(5) = akp%value(i,j,k)*f%value(i,j,k+1)
                  rt(6) = akm%value(i,j,k)*f%value(i,j,k-1)
                  rt(7) = -ap%value(i,j,k)*f%value(i,j,k)
                  rt(8) = con%value(i,j,k)
                  
                  ! Residual
                  res  = SUM(rt)
                  
                  DO n = 1,8
                     term = MAX(term,ABS(rt(n)))
                  END DO
                  
                  IF (ABS(res/term).GT.tmp) THEN
                     icon = 0
                  END IF
                  
                  blc  = blc + res
                  
               END IF
               
            END DO
         END DO
         
         denom = bl - ptz(k-1) * blm
         
         IF (ABS(denom/bl).LT.SMALL1) THEN
            denom = BIG
         END IF
         
         ptz(k) = blp / denom
         qtz(k) = (blc+blm*qtz(k-1))/denom
         
      END DO
      
      bl = 0.0
      
      DO k = objMalla%n2,2,-1
         
         bl = bl * ptz(k) + qtz(k)
         
         DO j = 2,objMalla%m2
            DO i = 2,objMalla%l2
               IF (ap%value(i,j,k).LT.BIG) THEN
                  f%value(i,j,k) = f%value(i,j,k) + bl
               ENDIF
            END DO
         END DO
         
      END DO
      
      DEALLOCATE(ptx,qtx)
      DEALLOCATE(pty,qty)
      DEALLOCATE(ptz,qtz)
      
      !WRITE(*,*) "Block Correction Scheme Performed in Z"
      
   END SUBROUTINE blockCorrectionScheme
   
   SUBROUTINE tdma3D(f)
      
      USE typeArreglo3D
      USE Coeficientes
      
      INTEGER          :: mm2,mm,ll2,ll,nn2,nn
      TYPE(arreglo3D)  :: f
      
      REAL(nP)        :: denom,temp
      
      REAL(nP), ALLOCATABLE :: ptx(:),qtx(:)
      REAL(nP), ALLOCATABLE :: pty(:),qty(:)
      REAL(nP), ALLOCATABLE :: ptz(:),qtz(:)
      
      ALLOCATE(ptx(ni),qtx(ni))
      ALLOCATE(pty(nj),qty(nj))
      ALLOCATE(ptz(nk),qtz(nk))
      
      mm2 = 2*objMalla%m2
      mm  = mm2 - 2
      
      ll2 = 2*objMalla%l2
      ll  = ll2 - 2
      
      nn2 = 2*objMalla%n2
      nn = nn2 - 2
      
      ! Aplicacion del TDMA en el plano XY con movimiento en Z
      DO kk = 2, nn
         
         k = MIN(kk,nn2-kk)
         
         ! Primera parte direccion X
         DO jj = 2, mm
            
            j = MIN(jj,mm2-jj)
            
            ptx(1) = 0.0
            qtx(1) = 0.0
            
            ! Forward Sustitution
            DO i = 2, objMalla%l2
               denom  = ap%value(i,j,k)-ptx(i-1)*aim%value(i,j,k)
               ptx(i) = aip%value(i,j,k) / denom
               
               temp   = con%value(i,j,k)
               temp   = temp + ajp%value(i,j,k) * f%value(i,j+1,k) &
               + ajm%value(i,j,k) * f%value(i,j-1,k)
               temp   = temp + akp%value(i,j,k) * f%value(i,j,k+1) &
               + akm%value(i,j,k) * f%value(i,j,k-1)
               
               qtx(i) = (temp + aim%value(i,j,k)*qtx(i-1)) / denom
            END DO
            
            ! Backward Sustitution
            DO i = objMalla%l2,2,-1
               f%value(i,j,k) = f%value(i+1,j,k) * ptx(i) + qtx(i)
            END DO
            
         END DO
         
         ! Segunda parte direccion Y
         DO ii = 2, ll
            
            i = MIN(ii,ll2-ii)
            
            pty(1) = 0.0
            qty(1) = 0.0
            
            ! Forward Sustitution
            DO j = 2, objMalla%m2
               denom  = ap%value(i,j,k)-pty(j-1)*ajm%value(i,j,k)
               pty(j) = ajp%value(i,j,k) / denom
               
               temp   = con%value(i,j,k)
               temp   = temp + aip%value(i,j,k) * f%value(i+1,j,k) &
               + aim%value(i,j,k) * f%value(i-1,j,k)
               temp   = temp + akp%value(i,j,k) * f%value(i,j,k+1) &
               + akm%value(i,j,k) * f%value(i,j,k-1)
               qty(j) = (temp + ajm%value(i,j,k)*qty(j-1)) / denom
            END DO
            
            ! Backward Sustitution
            DO j = objMalla%m2,2,-1
               f%value(i,j,k) = f%value(i,j+1,k) * pty(j) + qty(j)
            END DO
            
         END DO
         
      END DO
      
      ! Aplicacion del TDMA en el plano XZ con movimiento en Y
      DO jj = 2, mm
         
         j = MIN(jj,mm2-jj)
         
         DO ii = 2, ll
            
            i = MIN(ii,ll2-ii)
            
            ptz(1) = 0.0
            qtz(1) = 0.0
            
            ! Forward Sustitution
            DO k = 2, objMalla%n2
               denom  = ap%value(i,j,k)-ptz(k-1)*akm%value(i,j,k)
               ptz(k) = akp%value(i,j,k) / denom
               
               temp   = con%value(i,j,k)
               temp   = temp + aip%value(i,j,k) * f%value(i+1,j,k) &
               + aim%value(i,j,k) * f%value(i-1,j,k)
               temp   = temp + ajp%value(i,j,k) * f%value(i,j+1,k) &
               + ajm%value(i,j,k) * f%value(i,j-1,k)
               qtz(k) = (temp + akm%value(i,j,k)*qtz(k-1)) / denom
            END DO
            
            ! Backward Sustitution
            DO k = objMalla%n2,2,-1
               f%value(i,j,k) = f%value(i,j,k+1) * ptz(k) + qtz(k)
            END DO
            
         END DO
         
      END DO
      
      DEALLOCATE(ptx,qtx)
      DEALLOCATE(pty,qty)
      DEALLOCATE(ptz,qtz)
      
   END SUBROUTINE tdma3D
   
   SUBROUTINE useItpack2(f)
      
      USE typeArreglo3D
      USE Coeficientes
      
      REAL(nP)          :: vectorB((ni-2)*(nj-2)*(nk-2))
      
      TYPE (arreglo3D)  :: f
      
      ! Construye el arreglo AA & vectorB
      vectorB = 0.0
      n = 0
      DO nFil = 1, numNodes
         
         i = getCoord(nFil)%i
         j = getCoord(nFil)%j
         k = getCoord(nFil)%k
         
         n = n + 1
         
         AA(n) = ap%value(i,j,k)
         
         vectorB(nFil) = vectorB(nFil) + con%value(i,j,k)
         
         DO nVecino = 1, 6
            
            IF (getCoord(nFil)%nodeVec(nVecino).NE.0) THEN
               
               n = n + 1
               
               ! nVecino = 1 -- Cara E
               ! nVecino = 2 -- Cara W
               ! nVecino = 3 -- Cara N
               ! nVecino = 4 -- Cara S
               ! nVecino = 5 -- Cara T
               ! nVecino = 6 -- Cara B
               
               SELECT CASE (nVecino)
               CASE (1)
                  !matrixA(nFil,nPos) =
                  AA(n) = -aip%value(i,j,k)
               CASE (2)
                  !matrixA(nFil,nPos) = -aim%value(i,j,k)
                  AA(n) = -aim%value(i,j,k)
               CASE (3)
                  !matrixA(nFil,nPos) = -ajp%value(i,j,k)
                  AA(n) = -ajp%value(i,j,k)
               CASE (4)
                  !matrixA(nFil,nPos) = -ajm%value(i,j,k)
                  AA(n) = -ajm%value(i,j,k)
               CASE (5)
                  !matrixA(nFil,nPos) = -akp%value(i,j,k)
                  AA(n) = -akp%value(i,j,k)
               CASE (6)
                  !matrixA(nFil,nPos) = -akm%value(i,j,k)
                  AA(n) = -akm%value(i,j,k)
               END SELECT
               
            END IF
            
         END DO
         
      END DO
      
      itMax = 1000
      
      CALL setEnviroment(AA,JA,IA,vectorB,numNodes,nodesActTotal,itMax)
      
      ! Actualiza el valor de la variable f
      DO nFil = 1,numNodes
         
         i = getCoord(nFil)%i
         j = getCoord(nFil)%j
         k = getCoord(nFil)%k
         
         f%value(i,j,k) = vectorB(nFil)
         
      END DO
      
   END SUBROUTINE useItpack2
   
   SUBROUTINE setEnviroment(A,JA,IA,RHS,N,NZ,ITMAX)
      
      PARAMETER (ICTE = 12)
      
      INTEGER :: N,NZ,ITMAX
      INTEGER :: IA(N+1)
      INTEGER :: JA(NZ)
      INTEGER :: IPARM(ICTE)
      INTEGER :: IWKSP(3*N)
      DOUBLE PRECISION, ALLOCATABLE :: WKSP(:)
      DOUBLE PRECISION :: A(NZ), RHS(N), U(2*(3*N)),  RPARM(ICTE)
      
      INTENT (IN)      :: IA,JA,A
      INTENT (IN OUT)  :: RHS
      
      !Arreglos necesarios para el MULTIGRID
      
      INTEGER, allocatable 	:: JJA(:)
      INTEGER, allocatable 	:: IIA(:)
      REAL(nP), allocatable	:: UU(:)
      REAL(nP), allocatable	:: AAA(:)
      REAL(nP), allocatable	:: FF(:)
      REAL(nP), allocatable	:: IG(:)
      
      level = 0
      idgts = 0
      
      SELECT CASE (indAlgSol)
         
      CASE (2)
         
         NW = 5*N
         ALLOCATE(WKSP(NW))
         
         CALL DFAULT (IPARM, RPARM)
         IPARM(1) = ITMAX
         IPARM(2) = LEVEL
         IPARM(5) = 1
         IPARM(12) = IDGTS
         U = RHS
         CALL SSORSI (N,IA,JA,A,RHS,U,IWKSP,NW,WKSP,IPARM,RPARM,IERR)
         
      CASE (3)
         
         !---------------------------------------------------------------------
         
         allocate(JJA(2*(3*NZ+5*N)))
         allocate(IIA(2*(3*N)))
         allocate(AAA(2*(3*NZ+5*N)))
         allocate(UU(2*(3*N)))
         allocate(FF(2*(3*N)))
         
         IIA=0
         DO i=1,N+1
            IIA(i)=IA(i)
         END DO
         JJA=0
         AAA=0
         DO i=1,NZ
            JJA(i)=JA(i)
            AAA(i)=A(i)
         END DO
         FF=0
         DO i=1,N
            FF(i)=RHS(i)
         END DO
         !---------------------------------------------------------------------
         
         NDA  = SIZE(AAA)
         NDIA = SIZE(IIA)
         NDJA = SIZE(JJA)
         NDU  = SIZE(UU)
         NDF  = SIZE(FF)
         ALLOCATE(IG(INT(5.4*NDF)))
         NDIG = SIZE(IG)
         
         !------------ SWITCHES ---------------------------------------------
         
         MATRIX = 22
         ISWTCH = 4
         IOUT   = 10
         IPRINT = 10404
         IERR   = 0
         
         !------------------------------------------------------------------
         ECG1   = 0.0D0
         ECG2   = 0.25D0
         EWT2   = 0.35D0
         NWT    = 2
         NTR    = 0
         IERR   = 0
         IG     = 0
         U      =0.0D0
         
         LEVELX = 25
         IFIRST = 13
         NCYC   = 1036!8!10
         EPS    = 1.D-5
         MADAPT = 27
         NRD    = 1131
         NSOLCO = 110
         NRU    = 1131
         TIME   = 0.0D0
         
         
         CALL AMG1R5(AAA,IIA,JJA,UU,FF,IG,NDA,NDIA,NDJA,NDU,NDF,NDIG,N,MATRIX,&
         ISWTCH,IOUT,IPRINT,LEVELX,IFIRST,NCYC,EPS,MADAPT,NRD,NSOLCO,NRU,&
         ECG1,ECG2,EWT2,NWT,NTR,IERR)
         
      CASE (4)
         
         NW = 2*N!+ITMAX*4
         allocate(WKSP(NW))
         
         call DFAULT (IPARM, RPARM)
         IPARM(1) = ITMAX
         IPARM(2) = LEVEL
         IPARM(5) = 1
         IPARM(12) = IDGTS
         U = RHS
         call JSI (N,IA,JA,A,RHS,U,IWKSP,NW,WKSP,IPARM,RPARM,IERR)
         
      END SELECT
      
      ! Actualiza la solucion
      SELECT CASE (indAlgSol)
         
      CASE (2,4)
         RHS = U
         DEALLOCATE(WKSP)
      CASE (3)
         RHS = UU
         DEALLOCATE (IG)
         DEALLOCATE(JJA)
         DEALLOCATE(IIA)
         DEALLOCATE(AAA)
         DEALLOCATE(UU)
         DEALLOCATE(FF)
         
      END SELECT
      
   END SUBROUTINE setEnviroment
   
END MODULE solvers


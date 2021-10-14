!     Last change:  EF    3 Feb 2005   11:17 am

MODULE PRESION
   !Realiza la correcion de presion siguiendo los pasos
   !del algoritmo SIMPLE, y aplicando la interpolacion
   !de Rhie-Chow para calcular el valor de la velocidad
   !en las caras
   
   USE typePrecision
   USE typeArreglo3D
   USE typeMalla
   
   REAL(nP)        :: ssum,smax
   TYPE(arreglo3D) :: pp,p,dpx,dpy,dpz,dup,dvp,dwp
   TYPE(arreglo3D) :: F1,F2,F3
   
   contains
   
   SUBROUTINE setPresion()
      !Dimensiona y nombra todos los arreglos utilizados
      
      USE Propiedades
      USE LeerData
      
      ALLOCATE(pp%value(ni,nj,nk))
      ALLOCATE(p%value(ni,nj,nk))
      
      ALLOCATE(dpx%value(ni,nj,nk))
      ALLOCATE(dpy%value(ni,nj,nk))
      ALLOCATE(dpz%value(ni,nj,nk))
      ALLOCATE(dup%value(ni,nj,nk))
      ALLOCATE(dvp%value(ni,nj,nk))
      ALLOCATE(dwp%value(ni,nj,nk))
      
      ALLOCATE(F1%value(ni-1,nj-1,nk-1))
      ALLOCATE(F2%value(ni-1,nj-1,nk-1))
      ALLOCATE(F3%value(ni-1,nj-1,nk-1))
      
      pp%title  = "Correcci�n de Presi�n"
      p%title  = "Presion"
      dpx%title = "Ca�da de Presi�n en X"
      dpy%title = "Ca�da de Presi�n en Y"
      dpz%title = "Ca�da de Presi�n en Z"
      dup%title = "Coeficiente 1/apu"
      dvp%title = "Coeficiente 1/apv"
      dwp%title = "Coeficiente 1/apw"
      F1%title  = "Flujo de Masa en X"
      F2%title  = "Flujo de Masa en Y"
      F3%title  = "Flujo de Masa en Z"
      
      !Se inicializan todos los valores en cero
      dpx%value(:,:,:) = 0.0D0
      dpy%value(:,:,:) = 0.0D0
      dpz%value(:,:,:) = 0.0D0
      pp%value(:,:,:)  = 0.0D0
      p%value(:,:,:)   = 0.0D0
      dup%value(:,:,:) = 0.0D0
      dvp%value(:,:,:) = 0.0D0
      dwp%value(:,:,:) = 0.0D0
      
      F1%value(:,:,:) = 0.0D0
      F2%value(:,:,:) = 0.0D0
      F3%value(:,:,:) = 0.0D0
      
      if (DatPrev) then
         DO i = 1,objMalla%l1
            DO j = 1,objMalla%m1
               DO k = 1,objMalla%n1
                  
                  p%value(i,j,k) = temp1%value(i,j,k)
                  
               END DO
            END DO
         END DO
      end if
      
      DEALLOCATE(temp1%value)
      
   END SUBROUTINE setPresion
   
   SUBROUTINE calP(u,v,w)
      
      USE Coeficientes
      USE solvers
      USE Propiedades
      USE typeBorde
      USE LeerData
      
      integer  :: i,j,k,indiceDir
      real(nP) :: rhoFace, dpf, upro, dpxel, uface, dpxe
      real(nP) :: vol
      real(nP) :: PPE,PPW,PPN,PPS,PPT,PPB
      
      TYPE (arreglo3D) :: u,v,w
      
      INTENT (IN OUT) :: u,v,w
      
      ! C�lculo de P y P�
      ! Inicializa todos los coeficientes a cero
      CALL initCoef()
      DO indiceDir = 1, 3
         !Coeficientes para la ecuaci�n de correcci�n de presi�n P�
         DO i = 2,objMalla%l2
            DO j = 2,objMalla%m2
               DO k = 2,objMalla%n2
                  
                  SELECT CASE (indiceDir)
                     
                  CASE (1)
                     
                     IF (ni.gt.3.and.i.LE.objMalla%l2-1) THEN
                        
                        !Determinaci�n de la densidad en la cara e
                        rhoFace = InterLin(rho%value(i,j,k),rho%value(i+1,j,k),i,1)
                        
                        !Determinaci�n de los gradientes de presi�n interpolados en la cara e
                        dpxel = InterLin(dpx%value(i,j,k),dpx%value(i+1,j,k),i,1)
                        
                        !Determinaci�n del gradiente de presi�n en la cara de la celda
                        dpxe  = (p%value(i+1,j,k) - p%value(i,j,k))/(objMalla%x(i+1) - objMalla%x(i))
                        
                        !Interpolaci�n de 1/ap en la cara e de la celda
                        dpf   = InterLin(dup%value(i,j,k),dup%value(i+1,j,k),i,1)
                        !Evita la correcion de velocidades en las partes solidas de la estructura
                        IF (visc%value(i,j,k).EQ.visc%value(i+1,j,k).AND.visc%value(i,j,k).LT.1.0D+25) THEN
                           !Determinaci�n de la velocidad u* en la cara de la celda e
                           upro  = InterLin(u%value(i,j,k),u%value(i+1,j,k),i,1)
                           
                           !Volumen alrededor de la cara
                           vol   = objMalla%areaX(j,k) * (objMalla%x(i+1) - objMalla%x(i))
                           
                           !Correcci�n de la velocidad en la cara de la celda e
                           uface = upro - vol * dpf * (dpxe - dpxel)
                        else
                           uface = 0.0D0
                        END if
                        
                        !Correcci�n del flujo m�sico
                        F1%value(i,j,k)    = rhoface * objMalla%areaX(j,k) * uface
                        
                        !C�lculo de los coeficientes
                        aip%value(i,j,k)   = rhoface*dpf * objMalla%areaX(j,k) * objMalla%areaX(j,k)
                        aim%value(i+1,j,k) = aip%value(i,j,k)
                        
                     END IF
                     
                  CASE (2)
                     
                     IF (nj.gt.3.and.j.LE.objMalla%m2-1) THEN
                        
                        !Determinaci�n de la densidad en la cara n
                        rhoFace = InterLin(rho%value(i,j,k),rho%value(i,j+1,k),j,2)
                        
                        !Determinaci�n de los gradientes de presi�n interpolados en la cara n
                        dpxel = InterLin(dpy%value(i,j,k),dpy%value(i,j+1,k),j,2)
                        
                        !Determinaci�n del gradiente de presi�n en la cara de la celda
                        dpxe  = (p%value(i,j+1,k) - p%value(i,j,k))/(objMalla%y(j+1) - objMalla%y(j))
                        
                        !Interpolaci�n de An/ap en la cara de la celda
                        dpf   = InterLin(dvp%value(i,j,k),dvp%value(i,j+1,k),j,2)
                        !Evita la correcion de velocidades en las partes solidas de la estructura
                        IF (visc%value(i,j,k).EQ.visc%value(i,j+1,k).AND.visc%value(i,j,k).LT.1.0D+25) THEN
                           !Determinaci�n de la velocidad v* en la cara de la celda n
                           upro  = InterLin(v%value(i,j,k),v%value(i,j+1,k),j,2)
                           
                           !Volumen alrededor de la cara
                           vol   = objMalla%areaY(i,k) * (objMalla%y(j+1) - objMalla%y(j))
                           
                           !Correcci�n de la velocidad en la cara de la celda n
                           uface = upro - vol * dpf * (dpxe - dpxel)
                        else
                           uface = 0.0D0
                        END if
                        
                        !Correcci�n del flujo m�sico
                        F2%value(i,j,k) = rhoface * objMalla%areaY(i,k) * uface
                        
                        !C�lculo de los coeficientes
                        ajp%value(i,j,k)   = rhoface*dpf*objMalla%areaY(i,k) * objMalla%areaY(i,k)
                        ajm%value(i,j+1,k) = ajp%value(i,j,k)
                        
                     END IF
                     
                  CASE (3)
                     IF (nk.gt.3.and.k.LE.objMalla%n2-1) THEN
                        
                        !Determinaci�n de la densidad en la cara t
                        rhoFace = InterLin(rho%value(i,j,k),rho%value(i,j,k+1),k,3)
                        
                        !Determinaci�n de los gradientes de presi�n interpolados en la cara t
                        dpxel = InterLin(dpz%value(i,j,k),dpz%value(i,j,k+1),k,3)
                        
                        !Determinaci�n del gradiente de presi�n en la cara de la celda
                        dpxe  = (p%value(i,j,k+1) - p%value(i,j,k))/(objMalla%z(k+1) - objMalla%z(k))
                        !Interpolaci�n de At/ap en la cara de la celda
                        dpf   = InterLin(dwp%value(i,j,k),dwp%value(i,j,k+1),k,3)
                        !Evita la correcion de velocidades en las partes solidas de la estructura
                        IF (visc%value(i,j,k).EQ.visc%value(i,j,k+1).AND.visc%value(i,j,k).LT.1.0D+25) THEN
                           !Determinaci�n de la velocidad w* en la cara de la celda t
                           upro  = InterLin(w%value(i,j,k),w%value(i,j,k+1),k,3)
                           
                           !Volumen alrededor de la cara
                           vol   = objMalla%areaZ(i,j) * (objMalla%z(k+1) - objMalla%z(k))
                           
                           !Correcci�n de la velocidad en la cara de la celda t
                           uface = upro - vol * dpf * (dpxe - dpxel)
                        else
                           uface = 0.0D0
                        END IF
                        !Correcci�n del flujo m�sico
                        F3%value(i,j,k) = rhoface*objMalla%areaZ(i,j)*uface
                        
                        !C�lculo de los coeficientes de p'
                        akp%value(i,j,k)   = rhoface*dpf*objMalla%areaZ(i,j) * objMalla%areaZ(i,j)
                        akm%value(i,j,k+1) = akp%value(i,j,k)
                        
                     END IF
                     
                  END SELECT
                  
               END DO
            END DO
         END DO
         
      END DO
      
      smax = 0.0D0
      ssum = 0.0D0
      
      DO i=2,objMalla%l2
         DO j=2,objMalla%m2
            DO k=2,objMalla%n2
               
               !Determina el t�rmino ap de la ecuaci�n de p'
               ap%value(i,j,k)  = aip%value(i,j,k) + aim%value(i,j,k) + ajp%value(i,j,k) +  &
               ajm%value(i,j,k) + akp%value(i,j,k) + akm%value(i,j,k)
               
               !Determina el t�rmino fuente de la ecuaci�n de p'
               con%value(i,j,k) = F1%value(i-1,j,k) - F1%value(i,j,k) + F2%value(i,j-1,k) - &
               F2%value(i,j,k) + F3%value(i,j,k-1) - F3%value(i,j,k)
               
               !Determina el desbalance de masa global
               ssum  = ssum + con%value(i,j,k)
               
               !Determina el desbalance de masa m�ximo
               smax  = MAX(smax,ABS(con%value(i,j,k)))
               
            END DO
         END DO
      END DO
      
      !Calcula el nodo que se encuentra en el centro del dominio
      nii = ceiling(ni/2.0D0)
      njj = ceiling(nj/2.0D0)
      nkk = ceiling(nk/2.0D0)
      
      ap%value(nii,njj,nkk)  = 1.0D20
      con%value(nii,njj,nkk) = 0.0
      
      ! Resuelve el sistema de Ecuaciones
      CALL solveEq(pp)
      
      !Extrapola el valor de pp a los bordes
      CALL ppbound()
      
      ! Correcci�n de velocidades en los nodos y de los flujos
      DO i = 2,objMalla%l2
         DO j = 2,objMalla%m2
            DO k = 2,objMalla%n2
               
               ! Actualiza los flujo en las caras
               ! Flujo en la direcci�n X!
               if (i.lt.objMalla%l2) F1%value(i,j,k) = F1%value(i,j,k) - aip%value(i,j,k) &
               * (pp%value(i+1,j,k) - pp%value(i,j,k))
               ! Flujo en la direcci�n Y!
               if (j.lt.objMalla%m2) F2%value(i,j,k) = F2%value(i,j,k) - ajp%value(i,j,k) &
               * (pp%value(i,j+1,k) - pp%value(i,j,k))
               ! Flujo en la direcci�n Z!
               if (k.lt.objMalla%n2) F3%value(i,j,k) = F3%value(i,j,k) - akp%value(i,j,k) &
               * (pp%value(i,j,k+1) - pp%value(i,j,k))
               
               ! Correcion de velocidades
               ! Determina el valor de la presion corregida en las caras
               PPE = InterLin(pp%value(i,j,k),pp%value(i+1,j,k),i,1)
               PPW = InterLin(pp%value(i-1,j,k),pp%value(i,j,k),i-1,1)
               PPN = InterLin(pp%value(i,j,k),pp%value(i,j+1,k),j,2)
               PPS = InterLin(pp%value(i,j-1,k),pp%value(i,j,k),j-1,2)
               PPT = InterLin(pp%value(i,j,k),pp%value(i,j,k+1),k,3)
               PPB = InterLin(pp%value(i,j,k-1),pp%value(i,j,k),k-1,3)
               
               ! Actualiza el valor de las velocidades ( U = U* + du * (pe - pw))
               ! Velocidad U
               u%value(i,j,k) = u%value(i,j,k) - dup%value(i,j,k) * objMalla%areaX(j,k) * (PPE - PPW)
               
               ! Velocidad V
               v%value(i,j,k) = v%value(i,j,k) - dvp%value(i,j,k) * objMalla%areaY(i,k) * (PPN - PPS)
               
               ! Velocidad W
               w%value(i,j,k) = w%value(i,j,k) - dwp%value(i,j,k) * objMalla%areaZ(i,j) * (PPT - PPB)
               
               ! Actualiza el valor de presion (P = P + relax * P�). Los valores
               ! de P' son establecidos con respecto al VC de referencia nii,njj,nkk
               p%value(i,j,k) = p%value(i,j,k) + relax%value(4) * (pp%value(i,j,k) - pp%value(nii,njj,nkk))
               
            END DO
         END DO
      END DO
      
      !Extrapola el valor de presion hacia los bordes
      CALL pbound()
      
      !Calcula la caida de presion
      CALL PDrop()
      
      !Determina el valor de la presion en las esquinas
      CALL setEdgesValues(p)
      
   END SUBROUTINE calP
   
   SUBROUTINE pbound()
      
      USE typeBorde
      
      ! Extrapola la P hacia los bordes
      
      REAL (nP) :: wIJK, wLMN
      
      ! Direcci�n X:
      nDir  = 1
      nf    = 1
      do j = 2,objMalla%m2
         do k = 2,objMalla%n2
            
            if (get_BC_I1(j,k,nDir,nf).eq.WALL) then
               wijk = objMalla%fx(3)
            else
               wijk = 0.0D0
            end if
            
            i = 2
            p%value(1,j,k)    = p%value(i,j,k) + (p%value(i,j,k)-p%value(i+1,j,k)) * wijk
            
            if (get_BC_L1(j,k,nDir,nf).eq.WALL) then
               wlmn = objMalla%fxm(objMalla%l2)
            else
               wlmn = 0.0D0
            end if
            
            i = objMalla%l2
            p%value(objMalla%l1,j,k)   = p%value(i,j,k) + (p%value(i,j,k)-p%value(i-1,j,k)) * wlmn
            
         end do
      end do
      
      ! Direcci�n Y:
      nDir  = 2
      nf    = 2
      do i = 2,objMalla%l2
         do k = 2,objMalla%n2
            
            if (get_BC_J1(i,k,nDir,nf).eq.WALL) then
               wIJK = objMalla%fy(3)
            else
               wIJK = 0.0D0
            end if
            
            j = 2
            p%value(i,1,k)    = p%value(i,j,k) + (p%value(i,j,k)-p%value(i,j+1,k)) * wIJK
            
            if (get_BC_M1(i,k,nDir,nf).eq.WALL) then
               wLMN = objMalla%fym(objMalla%m2)
            else
               wLMN = 0.0D0
            end if
            
            j = objMalla%m2
            p%value(i,objMalla%m1,k)   = p%value(i,j,k) + (p%value(i,j,k)-p%value(i,j-1,k)) * wLMN
            
         end do
      end do
      
      ! direcci�n Z:
      nDir  = 3
      nf    = 3
      do i = 2,objMalla%l2
         do j = 2,objMalla%m2
            
            if (get_BC_K1(i,j,nDir,nf).eq.WALL) then
               wijk = objMalla%fz(3)
            else
               wijk = 0.0D0
            end if
            
            k = 2
            p%value(i,j,1)    = p%value(i,j,k) + (p%value(i,j,k)-p%value(i,j,k+1)) * wijk
            
            if (get_BC_N1(i,j,nDir,nf).eq.WALL) then
               wlmn = objMalla%fzm(objMalla%n2)
            else
               wlmn = 0.0D0
            end if
            
            k = objMalla%n2
            p%value(i,j,objMalla%n1)   = p%value(i,j,k) + (p%value(i,j,k)-p%value(i,j,k-1)) * wlmn
            
         end do
      end do
      
   END SUBROUTINE pbound
   
   SUBROUTINE ppbound()
      
      ! Extrapola la PP hacia los bordes
      
      ! Borde Este y Oeste
      
      DO j = 2,objMalla%m2
         DO k = 2,objMalla%n2
            
            i = 2
            pp%value(1,j,k) = pp%value(i,j,k)
            
            i = objMalla%l2
            pp%value(objMalla%l1,j,k) = pp%value(i,j,k)
            
         END DO
      END DO
      
      ! Borde Norte y Sur
      
      DO i = 2,objMalla%l2
         DO k = 2,objMalla%n2
            
            j = 2
            pp%value(i,1,k) = pp%value(i,j,k)
            
            j = objMalla%m2
            pp%value(i,objMalla%m1,k) = pp%value(i,j,k)
            
         END DO
      END DO
      
      ! Borde Top y Bottom
      
      DO i = 2,objMalla%l2
         DO j = 2,objMalla%m2
            
            k = 2
            pp%value(i,j,1) = pp%value(i,j,k)
            
            k = objMalla%n2
            pp%value(i,j,objMalla%n1) = pp%value(i,j,k)
            
         END DO
      END DO
      
   END SUBROUTINE ppbound
   
   SUBROUTINE PDrop()
      
      REAL(nP)          :: pe,pw
      
      ! C�lcula la ca�da de presi�n incluida en el t�rmino fuente
      
      DO i = 2,objMalla%l2
         DO j = 2,objMalla%m2
            DO k = 2,objMalla%n2
               
               ! Direcci�n X
               pe = InterLin(p%value(i,j,k),p%value(i+1,j,k),i,1)
               pw = InterLin(p%value(i-1,j,k),p%value(i,j,k),i-1,1)
               dpx%value(i,j,k) = (pe - pw)/objMalla%dx(i)
               
               ! Direcci�n Y
               pe = InterLin(p%value(i,j,k),p%value(i,j+1,k),j,2)
               pw = InterLin(p%value(i,j-1,k),p%value(i,j,k),j-1,2)
               dpy%value(i,j,k) = (pe - pw)/objMalla%dy(j)
               
               ! Direcci�n Z
               pe = InterLin(p%value(i,j,k),p%value(i,j,k+1),k,3)
               pw = InterLin(p%value(i,j,k-1),p%value(i,j,k),k-1,3)
               dpz%value(i,j,k) = (pe - pw)/objMalla%dz(k)
               
            END DO
         END DO
      END DO
      
   END SUBROUTINE Pdrop
   
   SUBROUTINE SetEdgesValues(f)
      ! Extrapola el valor de las variables en las esquinas del domino
      
      real(nP)                      	 :: tmp
      integer                       	 :: i,j,k
      type(arreglo3D), intent(in out) :: f
      
      tmp = 0.5D0
      ! Bordes Especiales    j <--|
      ! Direccion Z            ___| i
      i = 1
      j = 1
      f%value(i,j,:) = tmp * (f%value(i+1,j,:) + f%value(i,j+1,:))
      
      i = objMalla%l1
      j = 1
      f%value(i,j,:) = tmp * (f%value(i-1,j,:) + f%value(i,j+1,:))
      
      i = objMalla%l1
      j = objMalla%m1
      f%value(i,j,:) = tmp * (f%value(i-1,j,:) + f%value(i,j-1,:))
      
      i = 1
      j = objMalla%m1
      f%value(i,j,:) = tmp * (f%value(i+1,j,:) + f%value(i,j-1,:))
      
      ! Bordes Especiales    k <--|
      ! Direccion X            ___| j
      j = 1
      k = 1
      f%value(:,j,k) = tmp * (f%value(:,j+1,k) + f%value(:,j,k+1))
      
      j = objMalla%m1
      k = 1
      f%value(:,j,k) = tmp * (f%value(:,j-1,k) + f%value(:,j,k+1))
      
      j = 1
      k = objMalla%n1
      f%value(:,j,k) = tmp * (f%value(:,j+1,k) + f%value(:,j,k-1))
      
      j = objMalla%m1
      k = objMalla%n1
      f%value(:,j,k) = tmp * (f%value(:,j-1,k) + f%value(:,j,k-1))
      
      ! Bordes Especiales    k <--|
      ! Direccion Y            ___| i
      i = 1
      k = 1
      f%value(i,:,k) = tmp * (f%value(i+1,:,k) + f%value(i,:,k+1))
      
      i = objMalla%l1
      k = 1
      f%value(i,:,k) = tmp * (f%value(i-1,:,k) + f%value(i,:,k+1))
      
      i = 1
      k = objMalla%n1
      f%value(i,:,k) = tmp * (f%value(i+1,:,k) + f%value(i,:,k-1))
      
      i = objMalla%l1
      k = objMalla%n1
      f%value(i,:,k) = tmp * (f%value(i-1,:,k) + f%value(i,:,k-1))
      
   END SUBROUTINE SetEdgesValues
   
END MODULE PRESION

!     Last change:  EF    3 Feb 2005   11:17 am

MODULE Propiedades
   !Contiene las declaraciones de las propiedades
   !utilizadas en el programa (densidad, conductividad,
   !viscosidad), GamX toma el valor de Conductividad o
   !Viscosidad dependiendo de que variable se resuelva
   USE typeArreglo3D
   
   TYPE(arreglo3D) :: rho,gamX,gamT,visc,Cp,temp1
   
   CONTAINS
   
   SUBROUTINE setProp()
      !Dimensiona las propiedades de acuedo al tama�o de
      !la malla
      USE typeMalla
      
      ALLOCATE(rho%value(ni,nj,nk))
      ALLOCATE(gamT%value(ni,nj,nk))
      ALLOCATE(visc%value(ni,nj,nk))
      ALLOCATE(Cp%value(ni,nj,nk))
      ALLOCATE(temp1%value(ni,nj,nk))
      
      temp1%value(:,:,:) = 0.0D0
      visc%title = "Viscosidad"
      
   END SUBROUTINE setProp
   
END MODULE Propiedades

MODULE typeBorde
   !Este modulo contiene las condiciones de borde
   !para la velocidad y la temperatura
   USE typePrecision
   
   ! Posibilidades para las condiciones de borde de la velocidad
   INTEGER, parameter 	:: WALL    	= 1
   INTEGER, parameter 	:: INFLOW  	= 2
   INTEGER, parameter 	:: OUTFLOW 	= 3
   
   INTEGER, allocatable	:: bc_I1(:,:)
   INTEGER, allocatable	:: bc_J1(:,:)
   INTEGER, allocatable	:: bc_K1(:,:)
   INTEGER, allocatable	:: bc_L1(:,:)
   INTEGER, allocatable	:: bc_M1(:,:)
   INTEGER, allocatable	:: bc_N1(:,:)
   
   !Define los kbc de las condiciones de borde de velocidad
   character(len = 8)   :: kbc_I1, kbc_L1, kbc_J1, kbc_M1, kbc_K1, kbc_N1
   
   !Declara las variables correspondientes a los flujos de calor para
   !FLXC (Termino independiente) y FLXP (termino dependiente de la
   !variable temperatura, velocidad u otra), para cada pared
   TYPE BoundFlux2D
   REAL(nP), POINTER :: CI1(:,:),CJ1(:,:),CK1(:,:)
   REAL(nP), POINTER :: CL1(:,:),CM1(:,:),CN1(:,:)
   REAL(nP), POINTER :: PI1(:,:),PJ1(:,:),PK1(:,:)
   REAL(nP), POINTER :: PL1(:,:),PM1(:,:),PN1(:,:)
   END TYPE BoundFlux2D
   
   TYPE(BoundFlux2D) :: FLX
   
   !Condiciones de borde para la temperatura
   TYPE kbcIn
   INTEGER, POINTER :: I1(:,:),J1(:,:),K1(:,:)
   INTEGER, POINTER :: L1(:,:),M1(:,:),N1(:,:)
   END TYPE kbcIn
   
   TYPE(kbcIn)       :: kbc
   
   !FLujo de calor para cada pared volumen de control
   !No entra en los calculos, solo provee informacion
   !sobre la salida o entrada de flujo en el dominio
   !TYPE flujo
   !  REAL(nP), POINTER :: I1(:,:,:),J1(:,:,:),K1(:,:,:)
   !  REAL(nP), POINTER :: L1(:,:,:),M1(:,:,:),N1(:,:,:)
   !END TYPE flujo
   
   !TYPE(flujo) :: flux
   
   CONTAINS
   
   !SUBROUTINE setFlux()
   !Dimensiona el flujo para cada variable
   !Para temperatura: Flujo de calor
   !Para velocidad: Esfuerzos debido a fuerzas cortantes
   !   USE typeMalla
   !   USE readData
   
   !   ALLOCATE(flux%I1(nj,nk,nfmax))
   !   ALLOCATE(flux%J1(ni,nk,nfmax))
   !   ALLOCATE(flux%K1(ni,nj,nfmax))
   !   ALLOCATE(flux%L1(nj,nk,nfmax))
   !   ALLOCATE(flux%M1(ni,nk,nfmax))
   !   ALLOCATE(flux%N1(ni,nj,nfmax))
   
   !END SUBROUTINE setFlux
   
   SUBROUTINE setBC()
      !Asigna los valores a las condiciones de contorno de
      !velocidad
      USE typeMalla
      
      ! Nuevas condiciones de contorno
      allocate (bc_I1(nj,nk))
      allocate (bc_J1(ni,nk))
      allocate (bc_K1(ni,nj))
      allocate (bc_L1(nj,nk))
      allocate (bc_M1(ni,nk))
      allocate (bc_N1(ni,nj))
      
      ! Copia los valores de las condiciones de contorno
      
      bc_I1 = convertBC(kbc_I1)
      bc_J1 = convertBC(kbc_J1)
      bc_K1 = convertBC(kbc_K1)
      
      bc_L1 = convertBC(kbc_L1)
      bc_M1 = convertBC(kbc_M1)
      bc_N1 = convertBC(kbc_N1)
      
   END SUBROUTINE setBC
   
   INTEGER FUNCTION convertBC(bc_index)
   !Convierte el valor de texto en un valor numerico
   character(len = 8), intent(in)	:: bc_index
   IF (bc_index.eq."WALL") then
      convertBC = WALL
   else IF (bc_index.eq."INFLOW") then
      convertBC = INFLOW
   else IF (bc_index.eq."OUTFLOW") then
      convertBC = OUTFLOW
   END IF
END FUNCTION convertBC

! Utilidades para obtener las condiciones de borde en X
INTEGER FUNCTION get_BC_I1(j,k,nDir,nf)
INTEGER, intent (in)	:: j,k,nDir,nf
! Verifica que variable se est� resolviendo
SELECT CASE (nf)
   ! Velocidad U
CASE (1) 
   get_BC_I1 = 1
   IF (bc_I1(j,k).eq.WALL) get_BC_I1 = 2
   RETURN
   ! Velocidad V
CASE (2) 
   get_BC_I1 = 2
   IF (bc_I1(j,k).eq.WALL) get_BC_I1 = 1
   RETURN
   ! Velocidad W
CASE (3)
   get_BC_I1 = 2
   IF (bc_I1(j,k).eq.WALL) get_BC_I1 = 1
   RETURN
END SELECT
END FUNCTION get_BC_I1

INTEGER FUNCTION get_BC_L1(j,k,nDir,nf)
INTEGER, intent (in)	:: j,k,nDir,nf
! Verifica que variable se est� resolviendo
SELECT CASE (nf)
   ! Velocidad U
CASE (1) 
   get_BC_L1 = 1
   IF (bc_L1(j,k).eq.WALL) get_BC_L1 = 2
   RETURN
   ! Velocidad V
CASE (2) 
   get_BC_L1 = 2
   IF (bc_L1(j,k).eq.WALL) get_BC_L1 = 1
   RETURN
   ! Velocidad W
CASE (3)
   get_BC_L1 = 2
   IF (bc_L1(j,k).eq.WALL) get_BC_L1 = 1
   RETURN
END SELECT
END FUNCTION get_BC_L1


! Utilidades para obtener las condiciones de borde en Y
INTEGER FUNCTION get_BC_J1(i,k,nDir,nf)
INTEGER, intent (in)	:: i,k,nDir,nf
! Verifica que variable se est� resolviendo
SELECT CASE (nf)
   ! Velocidad U
CASE (1) 
   get_BC_J1 = 2
   IF (bc_J1(i,k).eq.WALL) get_BC_J1 = 1
   RETURN
   ! Velocidad V
CASE (2) 
   get_BC_J1 = 1
   IF (bc_J1(i,k).eq.WALL) get_BC_J1 = 2
   RETURN
   ! Velocidad W
CASE (3)
   get_BC_J1 = 2
   IF (bc_J1(i,k).eq.WALL) get_BC_J1 = 1
   RETURN
END SELECT
END FUNCTION get_BC_J1

INTEGER FUNCTION get_BC_M1(i,k,nDir,nf)
INTEGER, intent (in)	:: i,k,nDir,nf
! Verifica que variable se est� resolviendo
SELECT CASE (nf)
   ! Velocidad U
CASE (1) 
   get_BC_M1 = 2
   IF (bc_M1(i,k).eq.WALL) get_BC_M1 = 1
   RETURN
   ! Velocidad V
CASE (2) 
   get_BC_M1 = 1
   IF (bc_M1(i,k).eq.WALL) get_BC_M1 = 2
   RETURN
   ! Velocidad W
CASE (3)
   get_BC_M1 = 2
   IF (bc_M1(i,k).eq.WALL) get_BC_M1 = 1
   RETURN
END SELECT
END FUNCTION get_BC_M1

! Utilidades para obtener las condiciones de borde en Z
INTEGER FUNCTION get_BC_K1(i,j,nDir,nf)
INTEGER, intent (in)	:: i,j,nDir,nf
! Verifica que variable se est� resolviendo
SELECT CASE (nf)
   ! Velocidad U
CASE (1) 
   get_BC_K1 = 2
   IF (bc_K1(i,j).eq.WALL) get_BC_K1 = 1
   RETURN
   ! Velocidad V
CASE (2) 
   get_BC_K1 = 2
   IF (bc_K1(i,j).eq.WALL) get_BC_K1 = 1
   RETURN
   ! Velocidad W
CASE (3)
   get_BC_K1 = 1
   IF (bc_K1(i,j).eq.WALL) get_BC_K1 = 2
   RETURN
END SELECT
END FUNCTION get_BC_K1

INTEGER FUNCTION get_BC_N1(i,j,nDir,nf)
INTEGER, intent (in)	:: i,j,nDir,nf
! Verifica que variable se est� resolviendo
SELECT CASE (nf)
   ! Velocidad U
CASE (1) 
   get_BC_N1 = 2
   IF (bc_N1(i,j).eq.WALL) get_BC_N1 = 1
   RETURN
   ! Velocidad V
CASE (2) 
   get_BC_N1 = 2
   IF (bc_N1(i,j).eq.WALL) get_BC_N1 = 1
   RETURN
   ! Velocidad W
CASE (3)
   get_BC_N1 = 1
   IF (bc_N1(i,j).eq.WALL) get_BC_N1 = 2
   RETURN
END SELECT
END FUNCTION get_BC_N1

END MODULE typeBorde

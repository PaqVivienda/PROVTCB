!     Last change:  EF    3 Feb 2005   11:17 am

MODULE Propiedades
   !Contiene las declaraciones de las propiedades
   !utilizadas en el programa (densidad, conductividad,
   !viscosidad), GamX toma el valor de Conductividad o
   !Viscosidad dependiendo de que variable se resuelva
   USE typeArreglo3D
   
   TYPE(arreglo3D) :: rho,gamX,gamT,visc,Cp,temp1
   
   contains
   
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
   integer, parameter 	:: WALL    	= 1
   integer, parameter 	:: INFLOW  	= 2
   integer, parameter 	:: OUTFLOW 	= 3
   
   integer, allocatable	:: bc_I1(:,:)
   integer, allocatable	:: bc_J1(:,:)
   integer, allocatable	:: bc_K1(:,:)
   integer, allocatable	:: bc_L1(:,:)
   integer, allocatable	:: bc_M1(:,:)
   integer, allocatable	:: bc_N1(:,:)
   
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
   
   contains
   
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
   
   integer function convertBC(bc_index)
   !Convierte el valor de texto en un valor numerico
   character(len = 8), intent(in)	:: bc_index
   if (bc_index.eq."WALL") then
      convertBC = WALL
   else if (bc_index.eq."INFLOW") then
      convertBC = INFLOW
   else if (bc_index.eq."OUTFLOW") then
      convertBC = OUTFLOW
   end if
end function convertBC

! Utilidades para obtener las condiciones de borde en X
integer function get_BC_I1(j,k,nDir,nf)
integer, intent (in)	:: j,k,nDir,nf
! Verifica que variable se est� resolviendo
select case (nf)
   ! Velocidad U
case (1) 
   get_BC_I1 = 1
   if (bc_I1(j,k).eq.WALL) get_BC_I1 = 2
   return
   ! Velocidad V
case (2) 
   get_BC_I1 = 2
   if (bc_I1(j,k).eq.WALL) get_BC_I1 = 1
   return
   ! Velocidad W
case (3)
   get_BC_I1 = 2
   if (bc_I1(j,k).eq.WALL) get_BC_I1 = 1
   return
end select
end function get_BC_I1

integer function get_BC_L1(j,k,nDir,nf)
integer, intent (in)	:: j,k,nDir,nf
! Verifica que variable se est� resolviendo
select case (nf)
   ! Velocidad U
case (1) 
   get_BC_L1 = 1
   if (bc_L1(j,k).eq.WALL) get_BC_L1 = 2
   return
   ! Velocidad V
case (2) 
   get_BC_L1 = 2
   if (bc_L1(j,k).eq.WALL) get_BC_L1 = 1
   return
   ! Velocidad W
case (3)
   get_BC_L1 = 2
   if (bc_L1(j,k).eq.WALL) get_BC_L1 = 1
   return
end select
end function get_BC_L1


! Utilidades para obtener las condiciones de borde en Y
integer function get_BC_J1(i,k,nDir,nf)
integer, intent (in)	:: i,k,nDir,nf
! Verifica que variable se est� resolviendo
select case (nf)
   ! Velocidad U
case (1) 
   get_BC_J1 = 2
   if (bc_J1(i,k).eq.WALL) get_BC_J1 = 1
   return
   ! Velocidad V
case (2) 
   get_BC_J1 = 1
   if (bc_J1(i,k).eq.WALL) get_BC_J1 = 2
   return
   ! Velocidad W
case (3)
   get_BC_J1 = 2
   if (bc_J1(i,k).eq.WALL) get_BC_J1 = 1
   return
end select
end function get_BC_J1

integer function get_BC_M1(i,k,nDir,nf)
integer, intent (in)	:: i,k,nDir,nf
! Verifica que variable se est� resolviendo
select case (nf)
   ! Velocidad U
case (1) 
   get_BC_M1 = 2
   if (bc_M1(i,k).eq.WALL) get_BC_M1 = 1
   return
   ! Velocidad V
case (2) 
   get_BC_M1 = 1
   if (bc_M1(i,k).eq.WALL) get_BC_M1 = 2
   return
   ! Velocidad W
case (3)
   get_BC_M1 = 2
   if (bc_M1(i,k).eq.WALL) get_BC_M1 = 1
   return
end select
end function get_BC_M1

! Utilidades para obtener las condiciones de borde en Z
integer function get_BC_K1(i,j,nDir,nf)
integer, intent (in)	:: i,j,nDir,nf
! Verifica que variable se est� resolviendo
select case (nf)
   ! Velocidad U
case (1) 
   get_BC_K1 = 2
   if (bc_K1(i,j).eq.WALL) get_BC_K1 = 1
   return
   ! Velocidad V
case (2) 
   get_BC_K1 = 2
   if (bc_K1(i,j).eq.WALL) get_BC_K1 = 1
   return
   ! Velocidad W
case (3)
   get_BC_K1 = 1
   if (bc_K1(i,j).eq.WALL) get_BC_K1 = 2
   return
end select
end function get_BC_K1

integer function get_BC_N1(i,j,nDir,nf)
integer, intent (in)	:: i,j,nDir,nf
! Verifica que variable se est� resolviendo
select case (nf)
   ! Velocidad U
case (1) 
   get_BC_N1 = 2
   if (bc_N1(i,j).eq.WALL) get_BC_N1 = 1
   return
   ! Velocidad V
case (2) 
   get_BC_N1 = 2
   if (bc_N1(i,j).eq.WALL) get_BC_N1 = 1
   return
   ! Velocidad W
case (3)
   get_BC_N1 = 1
   if (bc_N1(i,j).eq.WALL) get_BC_N1 = 2
   return
END select
end function get_BC_N1

END MODULE typeBorde

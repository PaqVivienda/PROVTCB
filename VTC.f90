   !     Last change:  EF    2 Feb 2005   11:52 am
    PROGRAM VTC
      
      ! ----------------------------------------
      ! Programa VTC
      ! Programa para el estudio de viviendas
      ! t�rmicamente confortables
      ! Proyecto LUZ-FONACIT
      ! Facultad de Ingenieria
      ! Universidad del Zulia
      ! Maracaibo - Venezuela
      ! ----------------------------------------
      ! Realizado por: Ing. Javier E. Romero V.
      ! Modificado por: Ing. Javier V. Goicochea P.
      ! All Rigths Reserved
      ! Copyright 2002-2003
      !-----------------------------------------
      
      USE typePrecision
      USE UVWT
      USE presion
      USE output
      USE Propiedades
      USE LeerData
      
      !Declarar las variables que no estan definidas en ningun modulo
      INTEGER  :: iter,numIter,numTimeStep
      REAL(nP) :: timeStep,runTime,Tol,TIMER,div
      
      ! Llama a la subrutina DatosEntrada y lee los nodos y caras
      ! de la malla, para establecer el tama�o de los arreglos.
      ! La subroutina regresa los valores de ni,nj,nk.
      ! Tambien el archivo de salida y regresa los valores
      ! correspondientes a:
      ! 1) Dimensiones de la malla (x,y,z,xu,yv,zw)
      ! 2) Intervalos de Tiempo
      ! 3) Datos de la Simulacion
      CALL DatosEntrada(numIter,runTime,timeStep)
      
      numTimeStep = 0.0D0
      Tol = 1.0E-12         !Define el criterio de convergencia
      TIMER = 0.0D0         !Inicia el tiempo en cero
      
      !Dimensiona e inicializa todos los arreglos utilizados en el programa
      CALL InicioProb(TIMER)
      
      ! Inicializa el campo de velocidades en las caras (flujo en las caras)
      CALL initVelFace()
      div = 0.0D0
      DO WHILE (TIMER.LE.runTime)
         
         iter = 0
         smax = 1.0E20
         
         DO WHILE ((smax.GE.Tol).AND.(iter.LT.numIter))
            ! Imprime las variables en pantalla
            CALL OutScreen(iter,numTimeStep,timer)
            
            ! Construye el lazo para la soluci�n de u, v, w
            IF (lcaluvw) CALL calUVW(timer,TimeStep)
            
            ! Calculo de la temperatura
            IF (lcalT)  CALL calT(TimeStep,timer)
            
            ! Siguiente iteracion
            iter = iter + 1
            
         END DO
         
         CALL OutScreen(iter,numTimeStep,timer)
         
         IF (TIMER/60.ge.div) THEN
            CALL writeOutputFiles(u,v,w,T,visc,p,Timer) !Escribe los archivos de salida en texto y grafico
            div = div + 360                             !cada 360 minutos
         END IF
         
         !Siguiente intervalo de tiempo
         numTimeStep = numTimeStep + 1
         TIMER = TIMER + timestep
         
      END DO
      
      ! Escribe archivos de salida para el graficador
      CALL writeOutputFiles(u,v,w,T,visc,p,Timer)
      
      OPEN (14,FILE= 'fort.14') ;  close (14,STATUS='delete') !Elimina el archivo Temporal
      
      PRINT *,("")
      PRINT *, ""
      CALL system("PAUSE")
      
   END PROGRAM VTC
   
   SUBROUTINE InicioProb(TIMER)
      
      USE Coeficientes
      USE typeBorde
      USE Propiedades
      USE UVWT
      USE presion
      USE typeMalla
      USE solvers
      USE adapt
      
      REAL(nP) :: TIMER
      !Hace la llamada para las distintas subrutinas necesarias para establecer tanto
      !Condiciones Iniciales(KBC, etc),  como  variables  para
      !hacer calculos referentes a CM, Vel, etc.
      
      CALL setCoef()             !Inicializa coeficientes
      CALL setProp()             !Inicializa propiedades
      CALL setUVWT()             !Inicializa velocidades y temperatura
      CALL setPresion()          !Presiones y flujos masicos
      CALL setGridDimensions()   !Dimensiones de la malla
      !CALL setFlux()
      
      ! Crea los arreglos IA y JA necesarios para emplear los algoritmos de
      ! soluci�n de la libreria ITpack
      CALL setupSolver()
      
      ! Regresa todos los valores correspondientes de la Malla
      ! dx,dy,dz - areaX,areaY,areaZ - vol - Factores de Relajacion
      CALL makeGrid()
      
      ! Llama a la subrutina BEGIN del Adapt, Lee el archivo de
      ! las propiedades Densidad, Gamma, Gammat, Viscosidad
      CALL begin(u,v,w,t,TIMER)
      
      ! Inicializa las condiciones de borde de la velocidad
      CALL setBC()
      
   END SUBROUTINE InicioProb
   
   
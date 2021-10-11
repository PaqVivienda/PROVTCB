!     Last change:  EF   10 Feb 2005   11:03 am

MODULE LeerData
!Realiza la lectura de todos los datos iniciales necesarios
!para la corrida del programa

   USE typePrecision

   INTEGER              :: IndAlgSol
   CHARACTER(LEN = 120) :: title,nomRad,dir,plotf,nomPrev

   REAL(nP)     :: relaxUVW,relaxP,relaxT  ! Factores de relajacion

   INTEGER      :: Imon,Jmon,Kmon          ! Punto de Monitoreo

   LOGICAL      :: lcaluvw,lcalP,lcalT,DatPrev

   TYPE arreglo1D
     SEQUENCE
     REAL(nP), POINTER :: value(:)
   END TYPE

   TYPE(arreglo1D) :: relax

   CONTAINS

      SUBROUTINE DatosEntrada(numIter,runTime,timeStep)

         USE typeMalla

         CHARACTER(LEN = 120)  :: TEMP,nomDatos
         integer :: aa

  	 ! RUNOPTIONS
  	 INTEGER              :: numIter           !Numero de iteraciones por paso de
                                                   !tiempo
  	 REAL(nP)             :: timeStep,runTime  !Paso de tiempo y Tiempo total

  	 INTENT (IN OUT) :: numIter,runTime,timeStep

  	 ! GRID
  	 INTEGER               :: ncx,ncy,ncz      ! numero de caras
  	 INTEGER               :: ncvx,ncvy,ncvz   ! numero de volumenes de control

         DatPrev = .False.

         5 FORMAT (A100)

         !Lee los nombres de los archivos que contienen los datos
         !iniciales y los datos de radiación. Se tiene que incluir
         !el directorio en el que se encuentra cada archivo.
         OPEN (UNIT = 1, FILE = "nombre", STATUS = "OLD")
            READ (1,5) nomDatos              !Datos iniciales (Geometria, factores de relajacion)
            READ (1,5) nomRad                !Propiedades (Densidad, Cp, K) y datos de radiacion
            READ (1,*) temp
            !En caso de que se lea un archivo de una corrida previa
               if (temp.eq."s") then
                 DatPrev = .True.
               end if
            if (DatPrev) READ (1,5) nomPrev  !Lee el nombre del archivo de salida a ser leido
                                             !posteriormente
	 CLOSE(1)

         !Determina el nombre del directorio donde se encuentran
         !los archivos de entrada
         aa = 0
         do i = LEN_TRIM(NomDatos),0,-1
            if (aa.eq.0) then
               IF (nomDatos(i:i).eq."\") then
               DIR = nomDatos(1:i)
               aa = 1
               END IF
            end if
         end do

         !Lectura de los datos iniciales de la simulacion
	 OPEN (UNIT = 1, FILE = nomDatos)
            READ (1,*)
            READ (1,*) title
            READ (1,*)
            READ (1,*) numIter,IndAlgSol
            READ (1,*)
            READ (1,*) plotf
            READ (1,*)
            READ (1,*) runtime,timeStep
            READ (1,*)
            READ (1,*)
            READ (1,*) lcalUVW,lcalP,lcalT
            READ (1,*)
            READ (1,*)
            READ (1,*) relaxUVW,relaxP,relaxT
            READ (1,*)
            READ (1,*)
            READ (1,*) Imon,Jmon,Kmon

  	    !Transforma de minutos a segundos el tiempo
            !total de la corrida y el paso de tiempo
  	    runTime = runTime * 60
  	    timeStep = timeStep * 60

            !Lectura de la geometria de la malla
            !LECTURA DE LAS CARAS DE VOLUMENES DE CONTROL
    	    READ (1,*)
    	    READ (1,*) TEMP,ncx
    	    READ (1,*) TEMP,ncy
    	    READ (1,*) TEMP,ncz
    	    ALLOCATE(objMalla%xu(ncx),objMalla%yv(ncy),objMalla%zw(ncz))
    	    READ (1,*) (objmalla%xu(I),I=1,ncx)
    	    READ (1,*) (objmalla%yv(I),I=1,ncy)
    	    READ (1,*) (objmalla%zw(I),I=1,ncz)
  	    !LECTURA DE LOS NODOS
    	    READ (1,*)
    	    READ (1,*) TEMP,ncvx
    	    READ (1,*) TEMP,ncvy
    	    READ (1,*) TEMP,ncvz
    	    ni = ncvx + 2
    	    nj = ncvy + 2
    	    nk = ncvz + 2
    	    ALLOCATE(objMalla%x(ni),objMalla%y(nj),objMalla%z(nk))
    	    READ (1,*) (objmalla%x(I+1),I=1,ncvx)
    	    READ (1,*) (objmalla%y(I+1),I=1,ncvy)
    	    READ (1,*) (objmalla%z(I+1),I=1,ncvz)
	 CLOSE(1)
  !======================================

  	 ! nf = 1 -----> Velocidad U
  	 ! nf = 2 -----> Velocidad V
  	 ! nf = 3 -----> Velocidad W
  	 ! nf = 4 -----> Presión
  	 ! nf = 5 -----> Temperatura

  	 ! Dimensiona dinámicamente a relax
   	 ALLOCATE(relax%value(nfmax))

  	 ! Factores de relajación para cada variable
   	 relax%value(1:3)     = relaxUVW
  	 relax%value(4)       = relaxP
   	 relax%value(5)       = relaxT

         !Si el directorio tiene espacios en blanco, se cambian
         !temporalmente por *, para evitar errores en la corrida
         do aa=1,LEN_TRIM(DIR)
            if (DIR(aa:aa).eq."") then
              dir(aa:aa)="*"
            end if
         end do

     END SUBROUTINE DatosEntrada

END MODULE LeerData

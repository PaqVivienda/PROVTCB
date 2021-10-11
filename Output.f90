!     Last change:  DAN  19 Jul 2018    1:58 pm

MODULE OUTPUT
!Contiene todas las subrutinas correspondientes a la
!salida de datos

   USE LeerData
   USE typeMalla
   USE typeArreglo3D

   CHARACTER (LEN = 120) :: plotf1,nomArchivo1

   contains

      SUBROUTINE OutScreen(iter,numTimeStep,timer)
      !Subrutina para la salida de datos por pantalla

         USE typePrecision
         USE PRESION
         USE UVWT
         USE ADAPT

	 INTEGER         :: iter,numTimeStep
	 REAL(nP)        :: timer

	 INTENT(IN)  ::  iter,numTimeStep,timer

         IF (iter.EQ.0.AND.numTimeStep.EQ.0.0D0) THEN

         OPEN (unit= 4,FILE= 'Results.dat',STATUS='UNKNOWN')


      	    WRITE(*,*) "                              		      "
      	    WRITE(*,*) " -----------------------------------------        "
  	    WRITE(*,*) " Programa ProVTC                                  "
  	    WRITE(*,*) " Procesador para el estudio de Viviendas            "
  	    WRITE(*,*) " Termicamente Confortables                        "
  	    WRITE(*,*) " Proyecto LUZ-FONACIT                             "
  	    WRITE(*,*) " Facultad de Ingenieria                           "
  	    WRITE(*,*) " Universidad del Zulia                            "
  	    WRITE(*,*) " Maracaibo - Venezuela                            "
  	    WRITE(*,*) " -----------------------------------------        "
  	    WRITE(*,*) " Realizado  por: Ing. Javier E. Romero V.          "
            WRITE(*,*) " Modificado por: Ing. Javier V. Goicochea P.      "
            WRITE(*,*) " All Rigths Reserved                              "
  	    WRITE(*,*) " Copyright 2002-2004                              "
  	    WRITE(*,*) " -----------------------------------------        "
            WRITE(*,*) "                           			      "
      	    WRITE(*,*) " 	      RESULTADOS DEL ProVTC                "
      	    WRITE(*,*) "                                       	      "
            
            WRITE(4,*) "                              		      "
      	    WRITE(4,*) " -----------------------------------------        "
  	    WRITE(4,*) " Programa ProVTC                                  "
  	    WRITE(4,*) " Procesador para el estudio de Viviendas            "
  	    WRITE(4,*) " Termicamente Confortables                        "
  	    WRITE(4,*) " Proyecto LUZ-FONACIT                             "
  	    WRITE(4,*) " Facultad de Ingenieria                           "
  	    WRITE(4,*) " Universidad del Zulia                            "
  	    WRITE(4,*) " Maracaibo - Venezuela                            "
  	    WRITE(4,*) " -----------------------------------------        "
  	    WRITE(4,*) " Realizado  por: Ing. Javier E. Romero V.          "
            WRITE(4,*) " Modificado por: Ing. Javier V. Goicochea P.      "
            WRITE(4,*) " All Rigths Reserved                              "
  	    WRITE(4,*) " Copyright 2002-2004                              "
  	    WRITE(4,*) " -----------------------------------------        "
            WRITE(4,*) "                           			      "
      	    WRITE(4,*) " 	      RESULTADOS DEL ProVTC                "
      	    WRITE(4,*) "                                       	      "

         END IF

         IF (iter.EQ.0) THEN

            WRITE(*,*)  '					 '
            WRITE(*,*)  'Pasos en el Tiempo =',numTimeStep,'  Tiempo de simulacion =',CEILING(TIMER/60),'min.'
            WRITE(*,*)  '					 '
            WRITE(*,*) 'Temperatura Ambiente =', Tamb 
            WRITE(*,*)  '					 '
            WRITE(*,*) 'Imonitoreo=',Imon,',Jmonitoreo=',Jmon,',Kmonitoreo=',Kmon
            WRITE(*,*) '					                            '
            WRITE(*,*) 'Iter       velU     velV   	velW       Temp      smax       ssum'
            WRITE(*,*) '					                            '

            WRITE(4,*)  '					 '
            WRITE(4,*)  'Pasos en el Tiempo =',numTimeStep,'  Tiempo de simulacion =',CEILING(TIMER/60),'min.'
            WRITE(4,*)  '					 '
            WRITE(4,*) 'Temperatura Ambiente =', Tamb 
            WRITE(4,*)  '					 '
            WRITE(4,*) 'Imonitoreo=',Imon,',Jmonitoreo=',Jmon,',Kmonitoreo=',Kmon
            WRITE(4,*) '					                            '
            WRITE(4,*) 'Iter       velU     velV   	velW       Temp      smax       ssum'
            WRITE(4,*) '					                            '

         END IF

         IF (iter.EQ.1.OR.MOD(iter,1).EQ.0) THEN

  	     WRITE(*,10) iter,u%value(Imon,Jmon,kmon),v%value(Imon,Jmon,kmon),&
  	        	 w%value(Imon,Jmon,kmon),t%value(Imon,Jmon,kmon),smax,ssum

             WRITE(4,10) iter,u%value(Imon,Jmon,kmon),v%value(Imon,Jmon,kmon),&
  	        	 w%value(Imon,Jmon,kmon),t%value(Imon,Jmon,kmon),smax,ssum

         END IF

  10 FORMAT(1X,I4,2x,1P8E11.3,4X,1P8E11.3,4X,1P8E11.3,4X,1P8E11.3,4X,1P8E11.3,4X,1P8E11.3)

         RETURN

         CLOSE(4)

      END SUBROUTINE OutScreen

      subroutine writeOutputFiles(u,v,w,T,visc,p,Timer)

         type(arreglo3D), intent(in)	:: u,v,w,T,visc,p
         type(arreglo3D)    		:: f(nfmax)
         REAL(nP), INTENT (IN)          :: TIMER
         CHARACTER (LEN=6)              :: TEMP1
         INTEGER                        :: aa

         write (14,*) INT(TIMER/3600)


         !Escribe el nombre de los archivos de salida en un archivo temporal
         OPEN (14,FILE = "fort.14")
            REWIND (14)
            READ (14,*) temp1
            write (14,*) TRIM(DIR),TRIM(PLOTF),"_",TRIM(temp1),"_hrs.tec"
            BACKSPACE (14)
            READ (14,*) plotf1
            write (14,*) TRIM(DIR),TRIM(nomArchivo),"_",TRIM(temp1),"_hrs.txt"
            BACKSPACE (14)
            READ (14,*) nomArchivo1
            REWIND (14)
         CLOSE (14)

         !Vuelve a colocar los espacios en blanco en el nombre de los archivos
         DO aa=1,LEN_TRIM(nomArchivo1)
            if (nomArchivo1(aa:aa).eq."*") then
               nomArchivo1(aa:aa) = ""
            end if
         END DO

         DO aa=1,LEN_TRIM(plotf1)
            if (plotf1(aa:aa).eq."*") then
               plotf1(aa:aa) = ""
            end if
         END DO

         ! Copia los valores de las variables a imprimir
         f(1)  = u
         f(2)  = v
         f(3)  = w
         f(4)  = p
         f(5)  = T
         f(6)  = visc

         ! Llama al modulo Grid y escribe la malla en el archivo
         CALL writeMalla()

         ! Escribe el campo de solución para u,v,w,T
         call writeFileArreglos3D(f)

         ! Escribe el campo de solución para u,v,w,T en el archivo plotf
         call plotFileArreglos(f,timer)

      END SUBROUTINE writeOutputFiles

      subroutine writeMalla()
      !Escribe los valores de los nodos en cada direccion en
      !archivo de texto externo
         integer  :: iunit,iend,ibeg,i,jend,jbeg,j,kend,kbeg,k

         iunit = 10

         open (iunit,file = nomArchivo1)

         ! Crea la salida para la Malla
         1 FORMAT(/1x,6(1h*),3x,"Malla             ",3x,6(1h*)/9x,20(1h-))

         iend=0

         2 FORMAT(/,'I =',2x,8(i4,5x))
         4 FORMAT('X =',1p8e9.2)

         6 FORMAT(/,'J =',2x,8(i4,5x))
         8 FORMAT('Y =',1p8e9.2)

         10 FORMAT(/,'K =',2x,8(i4,5x))
         12 FORMAT('Z =',1p8e9.2)

         DO WHILE (iend.NE.objMalla%l1)
            ibeg=iend+1
            iend=iend+8
            iend=MIN(iend,objMalla%l1)
            WRITE(iunit,2) (i,i=ibeg,iend)
            WRITE(iunit,4) (objMalla%x(i),i=ibeg,iend)
         END DO

         DO WHILE (jend.NE.objMalla%m1)
            jbeg=jend+1
            jend=jend+8
            jend=MIN(jend,objMalla%m1)
            WRITE(iunit,6) (j,j=jbeg,jend)
            WRITE(iunit,8) (objMalla%y(j),j=jbeg,jend)
         END DO

         DO WHILE (kend.NE.objMalla%n1)
            kbeg=kend+1
            kend=kend+8
            kend=MIN(kend,objMalla%n1)
            WRITE(iunit,10) (k,k=kbeg,kend)
            WRITE(iunit,12) (objMalla%z(k),k=kbeg,kend)
         END DO
      
      end subroutine writeMalla

      subroutine writeFileArreglos3D(f)
      !Crea el archivo de datos de texto con los valores de todas las
      !variables estudiadas (f(nfmax)), para cada nodo

         type(arreglo3D), intent(in)   :: f(nfmax)
         integer                       :: i,j,k,nf,iunit,ibeg,iend

         iunit = 10

         open (iunit,file = nomArchivo1)

         ! crea la salida para las variables dependientes f
         14 format(//1x,6(1h*),3x,a24,3x,6(1h*)/9x,26(1h-))
         16 format(/'  i =',i6,6i9)
         18 format(1x,i2,3x,1p7e9.2)
         20 format(/'  k =',i6)
         22 format(1x,20a,1x,20(1h-))

         do nf=1,nfmax

            write (iunit,14) f(nf)%title

            do k = 1,nk

               ibeg = 1
               iend = ibeg + 6

               write(iunit,20) k

               if (ni.gt.iend) then

		            do while (iend.ne.ni)

                     iend = ibeg + 6
                     iend = min(iend,ni)

                     write(iunit,16) (i,i=ibeg,iend)
                     write(iunit,*) ' j ='

                     do j=nj,1,-1
                        write(iunit,18) j,(f(nf)%value(i,j,k),i=ibeg,iend)
                     end do

                     ibeg = iend + 1

                  end do

               else

                  iend = ibeg + 6
                  iend = min(iend,ni)

                  write(iunit,16) (i,i=ibeg,iend)
                  write(iunit,*) ' j ='

		            do j=nj,1,-1
                     write(iunit,18) j,(f(nf)%value(i,j,k),i=ibeg,iend)
                  end do

                  ibeg = iend + 1

               end if

            end do
         end do

         close (iunit)

      end subroutine writeFileArreglos3D

      SUBROUTINE plotFileArreglos(f,timer)
      !Crea el archivo de salida grafico

         type(arreglo3D), intent(in)	:: f(nfmax)
         REAL(nP) :: a,b,c,d

         nZone = 1

         iunit      = 10

         4 FORMAT('ZONE T = ZONA1, I=',I3,3X, 'J=',I3,3X,'K=',I3)
         5 FORMAT('  C=BLACK, F=POINT,')
         6 FORMAT(3(E16.6),E16.6,E16.6,E16.6,E16.6,2(E16.6))
         8 FORMAT(3(E16.6),E16.6,E16.6,E16.6,E16.5,2(E16.6))
        10 FORMAT(3(E16.6),E16.6,E16.6,E16.5,E16.6,2(E16.6))
        12 FORMAT(3(E16.6),E16.6,E16.5,E16.6,E16.6,2(E16.6))
        14 FORMAT(3(E16.6),E16.5,E16.6,E16.6,E16.6,2(E16.6))
        16 FORMAT(3(E16.6),E16.6,E16.6,E16.5,E16.5,2(E16.6))
        18 FORMAT(3(E16.6),E16.6,E16.5,E16.6,E16.5,2(E16.6))
        20 FORMAT(3(E16.6),E16.5,E16.6,E16.6,E16.5,2(E16.6))
        22 FORMAT(3(E16.6),E16.6,E16.5,E16.5,E16.6,2(E16.6))
        24 FORMAT(3(E16.6),E16.5,E16.6,E16.5,E16.6,2(E16.6))
        26 FORMAT(3(E16.6),E16.5,E16.5,E16.6,E16.6,2(E16.6))
        28 FORMAT(3(E16.6),E16.6,E16.5,E16.5,E16.5,2(E16.6))
        30 FORMAT(3(E16.6),E16.5,E16.6,E16.5,E16.5,2(E16.6))
        32 FORMAT(3(E16.6),E16.5,E16.5,E16.6,E16.5,2(E16.6))
        34 FORMAT(3(E16.6),E16.5,E16.5,E16.5,E16.6,2(E16.6))
        36 FORMAT(3(E16.6),E16.5,E16.5,E16.5,E16.5,2(E16.6))

         IF (nZone.EQ.1) THEN
            OPEN (iunit,FILE= plotf1, ACCESS='SEQUENTIAL',FORM='FORMATTED')
            WRITE (iunit,*) 'TITLE = "',title,'"'
	    WRITE (iunit,*) 'VARIABLES = "X"'
            WRITE (iunit,*) ' "Y           "'
	    WRITE (iunit,*) ' "Z           "'
	    WRITE (iunit,*) ' "Velocidad U "'
	    WRITE (iunit,*) ' "Velocidad V "'
	    WRITE (iunit,*) ' "Velocidad W "'
	    WRITE (iunit,*) ' "Presión     "'
	    WRITE (iunit,*) ' "Temperatura "'
	    WRITE (iunit,*) ' "Material    "'

          ELSE
      	    OPEN (iunit,FILE= plotf1, ACCESS='SEQUENTIAL',FORM='FORMATTED',STATUS='UNKNOWN' &
                                     ,POSITION='APPEND')
            REWIND(100)
          END IF

          WRITE (iunit,4) ni, nj, nk
          WRITE (iunit,5)

 	  DO k = 1,nk
   	  DO j = 1,nj
	  DO i = 1,ni

             a = f(1)%value(i,j,k)
             b = f(2)%value(i,j,k)
             c = f(3)%value(i,j,k)
             d = f(4)%value(i,j,k)

             IF (a.ge.0.0.and.b.ge.0.0.and.c.ge.0.0.and.d.ge.0.0) THEN
             WRITE (iunit,6) objMalla%x(i),objMalla%y(j),objMalla%z(k),(f(nf)%value(i,j,k),nf=1,5),tipo%value(i,j,k)
             ELSE IF (a.ge.0.0.and.b.ge.0.0.and.c.ge.0.0.and.d.lt.0.0) THEN
             WRITE (iunit,8) objMalla%x(i),objMalla%y(j),objMalla%z(k),(f(nf)%value(i,j,k),nf=1,5),tipo%value(i,j,k)
             ELSE IF (a.ge.0.0.and.b.ge.0.0.and.c.lt.0.0.and.d.ge.0.0) THEN
             WRITE (iunit,10) objMalla%x(i),objMalla%y(j),objMalla%z(k),(f(nf)%value(i,j,k),nf=1,5),tipo%value(i,j,k)
             ELSE IF (a.ge.0.0.and.b.lt.0.0.and.c.ge.0.0.and.d.ge.0.0) THEN
             WRITE (iunit,12) objMalla%x(i),objMalla%y(j),objMalla%z(k),(f(nf)%value(i,j,k),nf=1,5),tipo%value(i,j,k)
             ELSE IF (a.lt.0.0.and.b.ge.0.0.and.c.ge.0.0.and.d.ge.0.0) THEN
             WRITE (iunit,14) objMalla%x(i),objMalla%y(j),objMalla%z(k),(f(nf)%value(i,j,k),nf=1,5),tipo%value(i,j,k)
             ELSE IF (a.ge.0.0.and.b.ge.0.0.and.c.lt.0.0.and.d.lt.0.0) THEN
             WRITE (iunit,16) objMalla%x(i),objMalla%y(j),objMalla%z(k),(f(nf)%value(i,j,k),nf=1,5),tipo%value(i,j,k)
             ELSE IF (a.ge.0.0.and.b.lt.0.0.and.c.ge.0.0.and.d.lt.0.0) THEN
             WRITE (iunit,18) objMalla%x(i),objMalla%y(j),objMalla%z(k),(f(nf)%value(i,j,k),nf=1,5),tipo%value(i,j,k)
             ELSE IF (a.lt.0.0.and.b.ge.0.0.and.c.ge.0.0.and.d.lt.0.0) THEN
             WRITE (iunit,20) objMalla%x(i),objMalla%y(j),objMalla%z(k),(f(nf)%value(i,j,k),nf=1,5),tipo%value(i,j,k)
             ELSE IF (a.ge.0.0.and.b.lt.0.0.and.c.lt.0.0.and.d.ge.0.0) THEN
             WRITE (iunit,22) objMalla%x(i),objMalla%y(j),objMalla%z(k),(f(nf)%value(i,j,k),nf=1,5),tipo%value(i,j,k)
             ELSE IF (a.lt.0.0.and.b.ge.0.0.and.c.lt.0.0.and.d.ge.0.0) THEN
             WRITE (iunit,24) objMalla%x(i),objMalla%y(j),objMalla%z(k),(f(nf)%value(i,j,k),nf=1,5),tipo%value(i,j,k)
             ELSE IF (a.lt.0.0.and.b.lt.0.0.and.c.ge.0.0.and.d.ge.0.0) THEN
             WRITE (iunit,26) objMalla%x(i),objMalla%y(j),objMalla%z(k),(f(nf)%value(i,j,k),nf=1,5),tipo%value(i,j,k)
             ELSE IF (a.ge.0.0.and.b.lt.0.0.and.c.lt.0.0.and.d.lt.0.0) THEN
             WRITE (iunit,28) objMalla%x(i),objMalla%y(j),objMalla%z(k),(f(nf)%value(i,j,k),nf=1,5),tipo%value(i,j,k)
             ELSE IF (a.lt.0.0.and.b.ge.0.0.and.c.lt.0.0.and.d.lt.0.0) THEN
             WRITE (iunit,30) objMalla%x(i),objMalla%y(j),objMalla%z(k),(f(nf)%value(i,j,k),nf=1,5),tipo%value(i,j,k)
             ELSE IF (a.lt.0.0.and.b.lt.0.0.and.c.ge.0.0.and.d.lt.0.0) THEN
             WRITE (iunit,32) objMalla%x(i),objMalla%y(j),objMalla%z(k),(f(nf)%value(i,j,k),nf=1,5),tipo%value(i,j,k)
             ELSE IF (a.lt.0.0.and.b.lt.0.0.and.c.lt.0.0.and.d.ge.0.0) THEN
             WRITE (iunit,34) objMalla%x(i),objMalla%y(j),objMalla%z(k),(f(nf)%value(i,j,k),nf=1,5),tipo%value(i,j,k)
             ELSE IF (a.lt.0.0.and.b.lt.0.0.and.c.lt.0.0.and.d.lt.0.0) THEN
             WRITE (iunit,36) objMalla%x(i),objMalla%y(j),objMalla%z(k),(f(nf)%value(i,j,k),nf=1,5),tipo%value(i,j,k)
             END IF

   	  END DO
   	  END DO
   	  END DO

          write (iunit,*)
          write (iunit,*) TIMER

   	  CLOSE(iunit)

      END SUBROUTINE plotFileArreglos

END MODULE OUTPUT

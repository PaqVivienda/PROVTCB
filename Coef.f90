!     Last change:  LSC  26 Nov 2004   11:56 am

MODULE typeArreglo3D
	!Este modulo se usa cada vez que sea necesario declarar
	!variables que tengan una estructura tridimensional
	USE typePrecision
	
	TYPE arreglo3D
	SEQUENCE
	character(len = 20)  :: title
	REAL(nP), POINTER    :: value(:,:,:)
	END TYPE arreglo3D
	
	TYPE arreglo4D
	SEQUENCE
	REAL(nP), POINTER    :: value(:,:,:)
	END TYPE arreglo4D
	
	TYPE(arreglo4D) :: Tipo
	
END MODULE typeArreglo3D

MODULE Coeficientes
	!Contiene todos los coeficientes de la ecuacion
	!discretizada para la solucion de todas las variables
	!dependientes.
	
	USE typeArreglo3D
	
	TYPE(arreglo3D) :: ap,aip,aim,ajp,ajm,akp,akm,con
	
	contains
	
	SUBROUTINE setCoef()
		!Dimensiona todos los coeficientes
		
		USE typeMalla
		
		ALLOCATE(ap%value(ni,nj,nk))
		ALLOCATE(aip%value(ni,nj,nk))
		ALLOCATE(aim%value(ni,nj,nk))
		ALLOCATE(ajp%value(ni,nj,nk))
		ALLOCATE(ajm%value(ni,nj,nk))
		ALLOCATE(akp%value(ni,nj,nk))
		ALLOCATE(akm%value(ni,nj,nk))
		ALLOCATE(con%value(ni,nj,nk))
		
	END SUBROUTINE setCoef
	
	SUBROUTINE initCoef()
		!Inicializa todos los coeficientes a cero
		!para eliminar todos los valores que tenian
		!almacenados previamente
		
		ap%value  = 0.0D0
		con%value = 0.0D0
		aip%value = 0.0D0
		aim%value = 0.0D0
		ajp%value = 0.0D0
		ajm%value = 0.0D0
		akp%value = 0.0D0
		akm%value = 0.0D0
		
	END SUBROUTINE initCoef
	
	SUBROUTINE initCoefBordes()
		!Inicializa los coeficientes de los nodos adyacentes,
		!en todos los bordes del dominio, porque para las
		!tres direcciones de la velocidad, el valor de los
		!coeficientes en el interior del dominio va a ser igual
		!pero en los bordes se pueden presentar cambios.
		
		USE typeMalla
		
		INTEGER :: i,j,k
		
		ap%value  = 0.0D0
		con%value = 0.0D0
		
		do i = 1, ni, ni-1
			aip%value(i,:,:) = 0.0D0
			aim%value(i,:,:) = 0.0D0
			ajp%value(i,:,:) = 0.0D0
			ajm%value(i,:,:) = 0.0D0
			akp%value(i,:,:) = 0.0D0
			akm%value(i,:,:) = 0.0D0
		end do
		
		do j = 1, nj, nj-1
			aip%value(:,j,:) = 0.0D0
			aim%value(:,j,:) = 0.0D0
			ajp%value(:,j,:) = 0.0D0
			ajm%value(:,j,:) = 0.0D0
			akp%value(:,j,:) = 0.0D0
			akm%value(:,j,:) = 0.0D0
		end do
		
		do k = 1, nk, nk-1
			aip%value(:,:,k) = 0.0D0
			aim%value(:,:,k) = 0.0D0
			ajp%value(:,:,k) = 0.0D0
			ajm%value(:,:,k) = 0.0D0
			akp%value(:,:,k) = 0.0D0
			akm%value(:,:,k) = 0.0D0
		end do
		
	END SUBROUTINE initCoefBordes
	
END MODULE Coeficientes

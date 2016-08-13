!
!	调用GSS2.0求解稀疏矩阵
!	稀疏矩阵采用压缩列格式存储，下标从1开始,注意调用GSS时下标转化的过程
!
!	压缩列格式采用{ptr,ind,val}记录稀疏矩阵的非零元。
!	对于N阶有nnz个非零元的矩阵来说:
!	ind, val的长度是nnz。按列顺序记录每个非零元的行标和数值。
!	ptr的长度是N+1。ptr[i]记录第i列第一个非零元的位置,最后一个元素ptr[N+1]=nnz+1。
!		这样第i列的长度就是   ptr[i+1]-ptr[i]。
!
!	例如对于如下3阶的矩阵
!		1.0  0.0  5.0
!		0.0  3.0  6.0
!		2.0  4.0  7.0
!   则 	ptr[4]={1,3,5,8}
!		ind[7]={1,3,2,3,1,2,3}
!		val[7]={1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0}
!   
!   Copyright (c) 2006-present 陈英时(YingShiChen).    All Rights Reserved. 
!	Any problem,contact gsp@grusoft.com
!
	
	MODULE GSS_6_INTERFACE

      implicit none    
      INTERFACE 
		subroutine GRUS_init_heapinfo( );
 		!DEC$ ATTRIBUTES ALIAS:'_GRUS_init_heapinfo' :: grus_init_heapinfo
		!DEC$ ATTRIBUTES C :: GRUS_init_heapinfo
		end subroutine !GRUS_init_heapinfo 

		subroutine GRUS_exit_heapinfo( );
 		!DEC$ ATTRIBUTES ALIAS:'_GRUS_exit_heapinfo' :: grus_exit_heapinfo
		!DEC$ ATTRIBUTES C :: GRUS_exit_heapinfo
		end subroutine !GRUS_exit_heapinfo 
				
		function GSS_init_id(nRow,nCol,ptr,ind,val,m_type,setting)
		!DEC$ ATTRIBUTES ALIAS:'_GSS_init_id' :: GSS_init_id
		!DEC$ ATTRIBUTES C :: GSS_init_id
			integer GSS_init_id
			integer nRow,nCol,ptr(*),ind(*),m_type
			double precision val(*),setting(32)
		end function !GSS_init_id

		function GSS_symbol_id( nRow,nCol,ptr,ind,val )
		!DEC$ ATTRIBUTES ALIAS:'_GSS_symbol_id' :: GSS_symbol_id
		!DEC$ ATTRIBUTES C :: GSS_symbol_id
			integer GSS_symbol_id
			integer nRow,nCol,ptr(*),ind(*)
			double precision val(*)
		end function !GSS_symbol_udi

		function GSS_numeric_id( nRow,nCol,ptr,ind,val,hGSS )
		!DEC$ ATTRIBUTES ALIAS:'_GSS_numeric_id' :: GSS_numeric_id
		!DEC$ ATTRIBUTES C :: GSS_numeric_id
			integer hGSS,GSS_numeric_id
			integer nCol,ptr(*),ind(*)
			double precision val(*)
		end function !GSS_numeric_udi

		function GSS_solve_id( hGSS,nRow,nCol,ptr,ind,val,x )
		!DEC$ ATTRIBUTES ALIAS:'_GSS_solve_id' :: GSS_solve_id
		!DEC$ ATTRIBUTES C :: GSS_solve_id
			integer hGSS,GSS_solve_id,nRow,nCol,ptr(*),ind(*)
			double precision val(*),x(*)
		end function !GSS_solve_id

		function GSS_clear_id( hGSS )
		!DEC$ ATTRIBUTES ALIAS:'_GSS_clear_id' :: GSS_clear_id
		!DEC$ ATTRIBUTES C :: GSS_clear_id
			integer hGSS,GSS_clear_id
		end function !GSS_clear_id
	END    INTERFACE 
	end MODULE GSS_6_INTERFACE

	subroutine GSS_demo_6( dim,ptr,ind,val,rhs,m_type )	
! Copyright (c) 2006-present 陈英时(YingShiChen)
	use GSS_6_INTERFACE			!必须声明	use GSS_6_INTERFACE
	implicit none    
	integer dim,nnz,ptr(*),ind(*),loop,ret
	double precision val(*),rhs(*),start,setting(32)

	double precision  t_symbol,t_numeric,t_solve,x(:),err,rhs_norm
	allocatable:: x
	integer hGSS			!指向求解器
	integer i,j,r,strategy,nIterRefine,info,m_type
!	clock_t start;
	
	loop=3
	t_numeric=0.0;		t_solve=0.0
	nnz = ptr(dim+1)-1 
	write( *,* ), "dim=",dim,"nnz=",nnz,"m_type=",m_type
!/*************调用求解器******************/
	strategy=0;		info=0
	hGSS = 0
!下标转化
	do i = 1,dim+1
		ptr(i)=ptr(i)-1
	end do
	do i = 1,nnz
		ind(i)=ind(i)-1
	enddo	
	setting=0.0
	ret = GSS_init_id( dim,dim,ptr,ind,val,m_type,setting )  !初始化
	if( ret/=0 )	then
		write(*,*) "        ERROR at init GSS solver. ERROR CODE",ret
		goto 100;
	endif
	start=SECNDS(0.0) 
	hGSS = GSS_symbol_id( dim,dim,ptr,ind,val )	!符号分析
	t_symbol = SECNDS(start)
	write( *,* ),"symbol time=",t_symbol 
	if( hGSS==0 )	then
		write(*,*) "        ERROR at symbol."
		goto 100;
	endif
	start=SECNDS(0.0) 
	do i = 1,loop
		ret = GSS_numeric_id( dim,dim,ptr,ind,val,hGSS )  !数值分析
		if( ret /= 0 )	then
			write(*,*) "        ERROR at numeric. ERROR CODE",ret
			hGSS=0;			!必须归零,GSS已自动释放内存
			goto 100		
		endif
	end do 
	t_numeric = SECNDS(start)/loop
	write( *,* ),"numeric time=",t_numeric 
	
	allocate( x(dim) )
	start = SECNDS(0.0) 
	do i = 1,loop
		x(1:dim)=rhs(1:dim)
		call GSS_solve_id( hGSS,dim,dim,ptr,ind,val,x )  !回代求解
	end do
	t_solve = SECNDS(start)/loop
	write( *,* ), "solve time=",t_solve 
	call GSS_clear_id( hGSS )		!释放占用的内存

100	continue
	!下标复原
	do i = 1,dim+1
		ptr(i)=ptr(i)+1
	enddo
	do i = 1,nnz
		ind(i)=ind(i)+1
	enddo
	
	do i=1,dim
	rhs(i)=x(i)
	enddo
      deallocate( x )
!/*************||Ax-b||******************/
!	if( hGSS/=0 )	then
!		rhs_norm=0.0;			err=0.0
!		do i = 1,dim	
!			rhs_norm = rhs_norm+(rhs(i)*rhs(i))
!		enddo
!		rhs_norm = sqrt(rhs_norm);
!		do i = 1,dim	
!			do j = ptr(i), ptr(i+1)-1
!				r = ind(j)
!				rhs(r) = rhs(r)-val(j)*x(i)
!				if( (m_type==11 .or. m_type==12) .and. r/=i )then
!					rhs(i) = rhs(i)-val(j)*x(r)
!				endif
!			enddo		
!		end do
!		do i = 1,dim	
!			err = err+(rhs(i)*rhs(i))
!		end do
!		err = sqrt(err)	
!		write( *,* ), "Residual |Ax-b|=",err,"|b|=",rhs_norm 			
!		deallocate( x )
!	endif

	end subroutine	!GSS_demo_6

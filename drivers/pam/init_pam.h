
      vpd=0.4d0
      ed0=0.d0
      ep0=0.d0
      open(10,file="inputPAM.in")
      read(10,nml=pamvar)
      close(10)
!     Process command line variable change:
      if(n_command_arg/=0)then
         do i=1,n_command_arg
            call get_command_argument(i,arg_buffer)
            pos      = scan(arg_buffer,"=")
            nml_name = arg_buffer(1:pos-1);call s_cap(nml_name)
            nml_value= arg_buffer(pos+1:)
            select case(nml_name)
         case("VPD");read(nml_value,*)vpd
         case("ED0");read(nml_value,*)ed0
         case("EP0");read(nml_value,*)ep0
         case default
            print*,"No corresponging variable in NML"
         end select
      enddo
      endif
      gmu=xmu  
      if((ed0-ep0) > 0.d0)gzero=0.5*(ep0+ed0+sqrt((ep0-ed0)**2 + 4*Vpd**2))
      if((ed0-ep0) < 0.d0)gzero=0.5*(ep0+ed0-sqrt((ep0-ed0)**2 + 4*Vpd**2))
      if((ed0-ep0) /=0.d0)xmu=gmu+gzero !true ED chemical potential
      write(*,*)'shift mu to (from) = ',xmu,'(',gmu,')'
      write(*,*)'shift is           = ',gzero

!+-------------------------------------------------------------------+
!PROGRAM  : 
!TYPE     : Subroutine
!PURPOSE  : 
!+-------------------------------------------------------------------+
subroutine step_Gkret(ik,istep)
  integer :: ik,istep
  integer :: itp

  !Get the first collision integral I(0,0)
  if(istep==1)call Ik_ret(istep-1)
  do itp=0,istep-1
     Gkret(istep,itp)=Udelta(ik,istep-1)*Gkret(istep-1,itp)-&
          Vdelta(ik,istep-1)*Ikret(istep-1,itp)
     !print*,"Step",Gkret(istep,itp),istep,itp
  enddo
  Gkret(istep,istep)=-xi
  !print*,"Step",Gkret(istep,istep),istep,istep
  return
end subroutine step_Gkret

!+-------------------------------------------------------------------+
!PROGRAM  : 
!TYPE     : Subroutine
!PURPOSE  : 
!+-------------------------------------------------------------------+
subroutine restep_Gkret(ik,istep)
  integer :: ik,istep
  integer :: itp
  complex(8) :: avIkret
  call Ik_ret(istep) !build the new collision integral at istep
  do itp=0,istep-1
     avIkret=Ikret(istep,itp)+Ikret(istep-1,itp)
     avIkret=avIkret/2.d0
     Gkret(istep,itp)=Udelta(ik,istep-1)*Gkret(istep-1,itp)-&
          Vdelta(ik,istep-1)*avIkret
     !print*,"ReStep",Gkret(istep,itp),istep,itp
  enddo
  Gkret(istep,istep)=-xi!;print*,"ReStep",Gkret(istep,istep),istep,istep
  !print*,'==================================='
  return
end subroutine restep_Gkret

!+-------------------------------------------------------------------+
!PROGRAM  : 
!TYPE     : Function
!PURPOSE  : 
!+-------------------------------------------------------------------+
subroutine Ik_ret(istep)
  complex(8) :: I1
  integer    :: i,itp,istep
  !print*,"Ik_ret:"
  do itp=0,istep  ! itp = t_j (the index of this function)
     I1=zero
     do i=itp,istep !itp=t'_j istep=t_i
        !print*,Gkret(istep,i),Sret(i,itp)
        I1=I1+Gkret(istep,i)*Sret(i,itp)*dt
     enddo
     !print*,"Ir,t,t`",I1,istep,itp
     Ikret(istep,itp)=I1
  enddo
  if(istep==0)Ikret(0,0)=zero
  !print*,"= = = = = = = = = = ="
  return
end subroutine Ik_ret

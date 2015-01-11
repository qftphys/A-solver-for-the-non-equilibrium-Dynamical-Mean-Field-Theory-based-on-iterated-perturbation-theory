!+------------------------------------------------------------+
!PROGRAM  : 
!TYPE     : Subroutine
!PURPOSE  : 
!+------------------------------------------------------------+
subroutine step_Gkgtr(ik,istep)
  integer :: ik,istep
  integer :: it
  !if(istep==1)call Ik_gtr(istep-1)
  do it=0,istep-1
     Gkgtr(istep,it)=Udelta(ik,istep-1)*Gkgtr(istep-1,it)-&
          Vdelta(ik,istep-1)*Ikgtr(istep-1,it)
     !print*,"Step",Gkgtr(it,istep),it,istep
  enddo
  !Gkgtr(istep,istep-1)=(-xi+Gkless(istep-1,istep-1))*Udelta(ik,istep-1)+&
  !Ikgtr(istep-1,istep-1)*Vdelta(ik,istep-1)
  Gkgtr(istep,istep)=(Gkless(istep,istep)-xi)!*Udelta(ik,istep-1)+&
!       Ikgtr(istep-1,istep-1)*Vdelta(ik,istep-1)
  return
end subroutine step_Gkgtr

!+------------------------------------------------------------+
!PROGRAM  : 
!TYPE     : Subroutine
!PURPOSE  : 
!+------------------------------------------------------------+
subroutine restep_Gkgtr(ik,istep)
  integer :: ik,istep
  integer :: it
  complex(8) :: avIkgtr
  !call Ik_gtr(istep) !build the new "layer" at istep+1  
  do it=0,istep-1
     avIkgtr=Ikgtr(istep,it)+oldIkgtr(istep-1,it)
     avIkgtr=avIkgtr/2.d0
     Gkgtr(istep,it)=Udelta(ik,istep-1)*Gkgtr(istep-1,it)- &
          Vdelta(ik,istep-1)*avIkgtr
     !print*,"ReStep",Gkless(it,istep),it,istep
  enddo
  !Gkgtr(istep,istep)=Gkless(istep,istep)-xi
  !print*,"ReStep",Gkret(istep,istep),istep,istep
  !print*,'==================================='
  return
end subroutine restep_Gkgtr

!+------------------------------------------------------------+
!PROGRAM  : 
!TYPE     : Function
!PURPOSE  : 
!+------------------------------------------------------------+
subroutine Ik_gtr(istep)
  integer :: i,it,istep
  complex(8) :: I1,I2,Ib
  do it=0,istep
     I1=zero;I2=zero;Ib=zero
     do i=0,it-1
        I1=I1+Sret(istep,i)*Gkgtr(i,it)*dt
     enddo
     if(it==0)I1=zero
     do i=0,istep-1
        I2=I2+Sgtr(istep,i)*Gkadv(i,it)*dt
     enddo
     if(istep==0)I2=zero
     do i=0,Ltau-1
        Ib=Ib+Slceil(istep,i)*Gkrceil(i,it)*dtau
     enddo
     !print*,"Ir,t,t`",I1+I2-xi*Ib,it,istep
     Ikgtr(istep,it)=I1+I2-xi*Ib
  enddo
  !if(istep==0)Ikless(0,0)=-xi*Ib!zero
  return
end subroutine Ik_gtr



module buoyancy_mod
  implicit none

contains

  subroutine buoyancy(ncrms,tkelebuoy)
    use vars
    use params
    implicit none
    integer, intent(in) :: ncrms
    integer i,j,k,kb,icrm
    real(crm_rknd) betu, betd,coef
    real T_rho_0(ncrms,nzm), T_rho(ncrms,nzm)
    real(crm_rknd), dimension(ncrms,nzm), intent(inout) ::  tkelebuoy

    coef = 1./float(nx*ny)

    do k=1,nz
      do icrm=1,ncrms
       tkelebuoy(icrm,k)=0.
       T_rho_0(icrm,k)=0.
      end do
    end do

    !$acc parallel loop gang vector collapse(4) async(asyncid)
    do k=2,nzm
      do j=1,ny
        do i=1,nx
          do icrm = 1 , ncrms
            kb=k-1
            betu=adz(icrm,kb)/(adz(icrm,k)+adz(icrm,kb))
            betd=adz(icrm,k)/(adz(icrm,k)+adz(icrm,kb))

            dwdt(icrm,i,j,k,na)=dwdt(icrm,i,j,k,na) +  &
            bet(icrm,k)*betu* &
            ( tabs0(icrm,k)*(epsv*(qv(icrm,i,j,k)-qv0(icrm,k))-(qcl(icrm,i,j,k)+qci(icrm,i,j,k)-qn0(icrm,k)+qpl(icrm,i,j,k)+qpi(icrm,i,j,k)-qp0(icrm,k))) &
            +(tabs(icrm,i,j,k)-tabs0(icrm,k))*(1.+epsv*qv0(icrm,k)-qn0(icrm,k)-qp0(icrm,k)) ) &
            + bet(icrm,kb)*betd* &
            ( tabs0(icrm,kb)*(epsv*(qv(icrm,i,j,kb)-qv0(icrm,kb))-(qcl(icrm,i,j,kb)+qci(icrm,i,j,kb)-qn0(icrm,kb)+qpl(icrm,i,j,kb)+qpi(icrm,i,j,kb)-qp0(icrm,kb))) &
            +(tabs(icrm,i,j,kb)-tabs0(icrm,kb))*(1.+epsv*qv0(icrm,kb)-qn0(icrm,kb)-qp0(icrm,kb)) )
          end do ! i
        end do ! j
      end do ! k
    end do ! k

    do k=2,nzm
      do j=1,ny
        do i=1,nx
          do icrm=1,ncrms
            T_rho_0(icrm,k) =  T_rho_0(icrm,k) + tabs0(icrm,k)*(1.+epsv*qv0(icrm,k)-qn0(icrm,k)-qp0(icrm,k))
          end do ! i
        end do ! j
      end do ! k
    end do ! k


    do k=2,nzm
      do j=1,ny
        do i=1,nx
          do icrm=1,ncrms
            T_rho_0(icrm,k) =  T_rho_0(icrm,k) * coef
            T_rho(icrm,k) = tabs(icrm,i,j,k)*(1.+epsv*qv(icrm,i,j,k)-(qcl(icrm,i,j,k)+qci(icrm,i,j,k)+qpl(icrm,i,j,k)+qpi(icrm,i,j,k)))
            tkelebuoy(icrm,k) = tkelebuoy(icrm,k) + (T_rho(icrm,k)-T_rho_0(icrm,k))*dwdt(icrm,i,j,k,na)
          end do ! i
        end do ! j
      end do ! k
    end do ! k

    do k=1,nz
      do icrm=1,ncrms
       tkelebuoy(icrm,k)= 9.81*tkelebuoy(icrm,k)*coef/T_rho_0(icrm,k)
      end do
    end do


  end subroutine buoyancy

end module buoyancy_mod

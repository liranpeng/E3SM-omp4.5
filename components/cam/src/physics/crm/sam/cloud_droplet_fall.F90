
subroutine cloud_droplet_fall(ncrms)


! Sedimentation of cloud droplets

use vars
use microphysics, only: micro_field, index_water_vapor, &
     Nc0, sigmag,  rho_water, &
     qtot_sed, prec_accum,do_chunked_energy_budgets
!use micro_params
use params

implicit none
integer, intent(in) :: ncrms
integer :: icrm
integer i,j,k, kb, kc, iqcl
integer, dimension(ncrms) :: return_flag
real coef,dqcl,lat_heat,vt_liq, coef_cl
real omnu, omnc, omnd, qclu, qclc, qcld, tmp_theta, tmp_phi
real fz(ncrms,nx,ny,nz)
integer, allocatable :: kmax(:)
integer, allocatable :: kmin(:)

allocate( kmax(ncrms) )
allocate( kmin(ncrms) )


kmax(:)=0
kmin(:)=nzm+1

do k = 1,nzm
 do j = 1, ny
  do i = 1, nx
    do icrm = 1 , ncrms
      if(qcl(icrm,i,j,k).gt.0.) then
        ! find range of vertical levels with cloud liquid
        kmin(icrm) = min(kmin(icrm),k)
        kmax(icrm) = max(kmax(icrm),k)
      end if
    end do
  end do
 end do
end do

return_flag(:) = 0
! Do not compute sedimentation if no cloud liquid is present
do icrm = 1 , ncrms
  if(kmax(icrm).lt.kmin(icrm)) then
    return_flag(icrm) = 1
  end if
end do

fz(:,:,:,:) = 0.
!
! Take into account sedimentation of cloud water which may be important for stratocumulus case.
! Parameterization of sedimentation rate is taken from GCSS WG1 DYCOMS2_RF2 case, and base on
! Rogers and Yau, 1989
coef_cl = 1.19e8*(3./(4.*3.1415*rho_water*Nc0*1.e6))**(2./3.)

! Compute cloud ice flux (using flux limited advection scheme, as in
! chapter 6 of Finite Volume Methods for Hyperbolic Problems by R.J.
! LeVeque, Cambridge University Press, 2002). 
do k = max(1,min(kmin)-1),max(kmax)
  ! Set up indices for x-y planes above and below current plane.
   kc = min(nzm,k+1)
   kb = max(1,k-1)
   do j = 1,ny
    do i = 1,nx
      do icrm = 1 , ncrms  
         ! CFL number based on grid spacing interpolated to interface i,j,k-1/2
         coef = dtn/(0.5*(adz(icrm,kb)+adz(icrm,k))*dz(icrm))
         ! Compute cloud liquid density in this cell and the ones above/below.
         ! Since cloud liquid is falling, the above cell is u (upwind),
         ! this cell is c (center) and the one below is d (downwind). 

         qclu = rho(icrm,kc)*qcl(icrm,i,j,kc)
         qclc = rho(icrm,k) *qcl(icrm,i,j,k) 
         qcld = rho(icrm,kb)*qcl(icrm,i,j,kb) 

         ! From Rogers and Yau.  Leading coefficient (assumed to be uniform in space)
         !   is computed above.  Depends on (rho*qcl)^(2/3).  Small offset of 1.e-12
         !   prevents issues with raising zero to a fractional power.
         vt_liq = coef_cl*(qclc+1.e-12)**(2./3.) &
              *exp(5.*log( sigmag(qclc) )**2)

         ! Use MC flux limiter in computation of flux correction.
         ! (MC = monotonized centered difference).
         if (qclc.eq.qcld) then
            tmp_phi = 0.
         else
            tmp_theta = (qclu-qclc)/(qclc-qcld)
            tmp_phi = max(0.,min(0.5*(1.+tmp_theta),2.,2.*tmp_theta))
         end if

         ! Compute limited flux.
         ! Since falling cloud liquid is a 1D advection problem, this
         ! flux-limited advection scheme is monotonic.
         fz(icrm,i,j,k) = -vt_liq*(qclc - 0.5*(1.-coef*vt_liq)*tmp_phi*(qclc-qcld))
      end do
   end do
  end do
end do
fz(:,:,:,nz) = 0.

! This only works for schemes that advect water vapor and cloud liquid mass
!   together as a single species.
iqcl = index_water_vapor
do k=max(1,min(kmin)-2),max(kmax)
  do j=1,ny
    do i=1,nx
      do icrm = 1 , ncrms 
         coef=dtn/(dz(icrm)*adz(icrm,k)*rho(icrm,k)) 
         ! The cloud liquid increment is the difference of the fluxes.
         dqcl=coef*(fz(icrm,i,j,k)-fz(icrm,i,j,k+1))
         ! Add this increment to both non-precipitating and total water.
         micro_field(icrm,i,j,k,iqcl)  = micro_field(icrm,i,j,k,iqcl)  + dqcl
         ! Include this effect in the total moisture budget.
         qifall(icrm,k) = qifall(icrm,k) + dqcl ! combine cloud liquid and ice sedimentation in qifall
         precflux(icrm,k) = precflux(icrm,k) - fz(icrm,i,j,k)*dtn/dz(icrm)

         ! The latent heat flux induced by the falling cloud liquid enters
         ! the liquid-ice static energy budget in the same way as the
         ! precipitation.  Note: use latent heat of sublimation. 
         lat_heat  = fac_cond*dqcl
         ! Add divergence of latent heat flux to liquid-ice static energy.
         t(icrm,i,j,k)  = t(icrm,i,j,k)  - lat_heat
         ! Add divergence to liquid-ice static energy budget.
         tlatqi(icrm,k) = tlatqi(icrm,k) - lat_heat ! combined cloud liquid and ice sedimentation in tlatqi
      end do
   end do
  end do
end do

if(do_chunked_energy_budgets) then
  do k=max(1,min(kmin)-2),max(kmax)
    do j=1,ny
      do i=1,nx
        do icrm = 1 , ncrms
        !bloss: save sedimentation tendencies for energy budgets in mse.f90
         coef=dtn/(dz(icrm)*adz(icrm,k)*rho(icrm,k))
        ! The cloud liquid increment is the difference of the fluxes.
         dqcl=coef*(fz(icrm,i,j,k)-fz(icrm,i,j,k+1))
         qtot_sed(icrm,i,j,k) = qtot_sed(icrm,i,j,k) + dqcl
       end do
     end do
   end do
  end do
!bloss(TODO): Include cloud drop sedimentation contribution to surface precipitation in 
!    MSE Budget outputs.  This cloud be important for fog...
!!$   prec_accum(1:nx,1:ny) = prec_accum(1:nx,1:ny) &
!!$        - dtn*fz(1:nx,1:ny,1)
 end if

deallocate( kmax )
deallocate( kmin )

end subroutine cloud_droplet_fall


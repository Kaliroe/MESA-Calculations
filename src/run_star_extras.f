! ***********************************************************************
!
!   Copyright (C) 2012  Bill Paxton
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful, 
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************
 
      module run_star_extras 

      use star_lib
      use star_def
      use const_def
      use const_def
      use chem_def
      use num_lib
      use crlibm_lib
      use binary_def
      
      implicit none
      
      integer :: time0, time1, clock_rate
      real(dp), parameter :: expected_runtime = 120 ! minutes

      integer, parameter :: restart_info_alloc = 1
      integer, parameter :: restart_info_get = 2
      integer, parameter :: restart_info_put = 3
      
      
      contains
      
      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         ! this is the place to set any procedure pointers you want to change
         ! e.g., other_wind, other_mixing, other_energy  (see star_data.inc)
         s% other_wind => brott_wind
         
         s% extras_startup => extras_startup
         s% extras_check_model => extras_check_model
         s% extras_start_step => extras_start_step
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns  

         s% job% warn_run_star_extras=.false.

      end subroutine extras_controls
      
!FUNCTIONS ADDED FOR CALCULATING THE MASS TRANSFER AT THE OUTER LAGRANGIAN POINT

      real(dp) function Phi_func(b,x,y,z) result(Phi_eq)
         real(dp), intent(in) :: x, y, z
         real(dp) :: q, x_cm, r, r2
         type(binary_info), pointer :: b
         include 'formats.inc'

         q = b% m(b% d_i)/b% m(b% a_i)
         x_cm = 1/(1+q)
         r = sqrt(x*x+y*y+z*z)
         r2 = sqrt(((1-x)*(1-x))+y*y+z*z)
         Phi_eq = -q/(r*(1+q)) -1/((q+1)*r2)&
                  -(1/2)*((x-x_cm)*(x-x_cm)+y*y)
      end function Phi_func

      real(dp) function dPhidx(b,x,y,z) result(derivative)
         real(dp), intent(in) :: x, y, z
         real(dp) :: q, x_cm, r, r2
         type(binary_info), pointer :: b
         include 'formats.inc'

         q = b% m(b% d_i)/b% m(b% a_i)
         x_cm = 1/(1+q)
         r = sqrt(x*x+y*y+z*z)
         r2 = sqrt(((1-x)*(1-x))+y*y+z*z)
         derivative = q*x/((1+q)*r*r*r) -(1-x)/((1+q)*r2*r2) -(x-x_cm)

      end function dPhidx

      real(dp) function ddPhidy(b,x,y,z) result(ddpdy)
         real(dp), intent(in) :: x, y, z
         real(dp) :: x_cm, r, q, r2
         type(binary_info), pointer :: b
         include 'formats.inc'

         q = b% m(b% d_i)/b% m(b% a_i)
         r = sqrt(x*x+y*y+z*z)
         r2 = sqrt(((1-x)*(1-x))+y*y+z*z)
         ddpdy = q/((1+q)*r*r*r) - 3*q*y*y/((1+q)*r*r*r*r*r) + &
            1/((1+q)*r2*r2*r2)- 3*y*y/((1+q)*r2*r2*r2*r2*r2) - 1
      end function ddPhidy

      real(dp) function ddPhidz(b,x,y,z) result(ddpdz)
         real(dp), intent(in) :: x, y, z
         real(dp) :: q, x_cm, r, r2
         type(binary_info), pointer :: b
         include 'formats.inc'

         q = b% m(b% d_i)/b% m(b% a_i)
         r = sqrt(x*x+y*y+z*z)
         r2 = sqrt(((1-x)*(1-x))+y*y+z*z)
         ddpdz = q/((1+q)*r*r*r) - 3*q*z*z/((1+q)*r*r*r*r*r) + &
            1/((1+q)*r2*r2*r2) - 3*z*z/((1+q)*r2*r2*r2*r2*r2)
      end function ddPhidz

      real(dp) function ddddPhidy(b,x,y,z) result(ddddpdy)
         real(dp), intent(in) :: x, y, z
         real(dp) :: q, x_cm, r, r2
         type(binary_info), pointer :: b
         include 'formats.inc'

         q = b% m(b% d_i)/b% m(b% a_i)
         x_cm = 1/(1+q)
         r = sqrt(x*x+y*y+z*z)
         r2 = sqrt(((1-x)*(1-x))+y*y+z*z)
         ddddpdy = -(1/(1+q))*(q*(3/(r*r*r*r*r) - 15*y*y/(r*r*r*r*r*r*r) + 6/(r*r*r*r*r) &
            - 30*y*y/(r*r*r*r*r*r*r) - 45*y*y/(r*r*r*r*r*r*r) + 105*y*y*y*y/(r*r*r*r*r*r*r*r*r))&
            + 3/(r2*r2*r2*r2*r2) - 15*y*y/(r2*r2*r2*r2*r2*r2*r2) + 6/(r2*r2*r2*r2*r2*r2*r2)&
            - 30*y*y/(r2*r2*r2*r2*r2*r2*r2) - 45*y*y/(r2*r2*r2*r2*r2*r2*r2)&
            + 105*y*y*y*y/(r2*r2*r2*r2*r2*r2*r2*r2*r2))
      end function ddddPhidy

      real(dp) function ddddPhidz(b,x,y,z) result(ddddpdz)
         real(dp), intent(in) :: x, y, z
         real(dp) :: q, x_cm, r, r2
         type(binary_info), pointer :: b
         include 'formats.inc'

         q = b% m(b% d_i)/b% m(b% a_i)
         x_cm = 1/(1+q)
         r = sqrt(x*x+y*y+z*z)
         r2 = sqrt(((1-x)*(1-x))+y*y+z*z)
         ddddpdz = -(1/(1+q))*(q*(3/(r*r*r*r*r) - 15*z*z/(r*r*r*r*r*r*r) + 6/(r*r*r*r*r) &
            - 30*z*z/(r*r*r*r*r*r*r) - 45*z*z/(r*r*r*r*r*r*r) + 105*z*z*z*z/(r*r*r*r*r*r*r*r*r))&
            + 3/(r2*r2*r2*r2*r2) - 15*z*z/(r2*r2*r2*r2*r2*r2*r2) + 6/(r2*r2*r2*r2*r2*r2*r2)&
            - 30*z*z/(r2*r2*r2*r2*r2*r2*r2) - 45*z*z/(r2*r2*r2*r2*r2*r2*r2)&
            + 105*z*z*z*z/(r2*r2*r2*r2*r2*r2*r2*r2*r2))
      end function ddddPhidz

      real(dp) function ddddPhidyz(b,x,y,z) result(ddddpdyz)
         real(dp), intent(in) :: x, y, z
         real(dp) :: q, x_cm, r, r2
         type(binary_info), pointer :: b
         include 'formats.inc'

         q = b% m(b% d_i)/b% m(b% a_i)
         x_cm = 1/(1+q)
         r = sqrt(x*x+y*y+z*z)
         r2 = sqrt(((1-x)*(1-x))+y*y+z*z)
         ddddpdyz = (1/(1+q))*(q*(-3/(r*r*r*r*r) + 15*z*z/(r*r*r*r*r*r*r) + 15*y*y/&
            (r*r*r*r*r*r*r) - 75*z*z*y*y/(r*r*r*r*r*r*r)) - 3/(r2*r2*r2*r2*r2)&
            + 15*z*z/(r2*r2*r2*r2*r2*r2*r2) + 15*y*y/(r2*r2*r2*r2*r2*r2*r2)&
            - 105*z*z*y*y/(r2*r2*r2*r2*r2*r2*r2*r2*r2))
      end function ddddPhidyz

      real(dp) function find_L1(b) result(L1)
         real(dp) :: limit, tolerance, x, upper_bound, lower_bound, dPhi_new
         type(binary_info), pointer :: b
         include 'formats.inc'

         upper_bound = 1d0
         lower_bound = 0d0
         x = 0d0
         limit = abs(upper_bound-lower_bound)
         tolerance = 0.00001d0
  
         do while (limit > tolerance)
            x = (lower_bound+upper_bound)/2
            dPhi_new = dPhidx(b,x,0d0,0d0)
            if (dPhi_new > 0) then
               lower_bound = x
            else if (dPhi_new < 0) then
               upper_bound = x
            else
               exit
            end if
            limit = abs(upper_bound-lower_bound)
         end do
         L1 = (upper_bound + lower_bound)/2
      end function find_L1    

      real(dp) function find_L2(b) result(L2)
         real(dp) :: limit, tolerance, x, upper_bound, lower_bound, dPhi_new,q
         type(binary_info), pointer :: b
         include 'formats.inc'

         q = b% m(b% d_i)/b% m(b% a_i)

         if (q < 1) then
            upper_bound = 0d0
            lower_bound = -1d0
         end if

         if (q .GE. 1) then
            upper_bound = 2d0
            lower_bound = 1d0
         end if 

         x = 0d0

         limit = abs(upper_bound-lower_bound)
         tolerance = 0.00001d0
  
         do while (limit > tolerance)
            x = (lower_bound+upper_bound)/2
            dPhi_new = dPhidx(b,x,0d0,0d0)
            if (dPhi_new > 0) then
               lower_bound = x
            else if (dPhi_new < 0) then
               upper_bound = x
            else
               exit
            end if
            limit = abs(upper_bound-lower_bound)
         end do
         L2 = (upper_bound + lower_bound)/2
      end function find_L2

      real(dp) function find_L3(b) result(L3)
         real(dp) :: limit, tolerance, x, upper_bound, lower_bound, dPhi_new,q
         type(binary_info), pointer :: b
         include 'formats.inc'

         q = b% m(b% d_i)/b% m(b% a_i)

         if (q .GE. 1) then
            upper_bound = 0d0
            lower_bound = -1d0
         end if

         if (q < 1) then
            upper_bound = 2d0
            lower_bound = 1d0
         end if 

         x = 0d0

         limit = abs(upper_bound-lower_bound)
         tolerance = 0.00001d0
  
         do while (limit > tolerance)
            x = (lower_bound+upper_bound)/2
            dPhi_new = dPhidx(b,x,0d0,0d0)
            if (dPhi_new > 0) then
               lower_bound = x
            else if (dPhi_new < 0) then
               upper_bound = x
            else
               exit
            end if
            limit = abs(upper_bound-lower_bound)
         end do
         L3 = (upper_bound + lower_bound)/2
      end function find_L3

      subroutine brott_wind(id, Lsurf, Msurf, Rsurf, Tsurf, w, ierr)
         use star_def
         integer, intent(in) :: id
         real(dp), intent(in) :: Lsurf, Msurf, Rsurf, Tsurf ! surface values (cgs)
         ! NOTE: surface is outermost cell. not necessarily at photosphere.
         ! NOTE: don't assume that vars are set at this point.
         ! so if you want values other than those given as args,
         ! you should use values from s% xh(:,:) and s% xa(:,:) only.
         ! rather than things like s% Teff or s% lnT(:) which have not been set yet.
         real(dp), intent(out) :: w ! wind in units of Msun/year (value is >= 0)
         integer, intent(out) :: ierr

         integer :: h1, he4
         real(dp) :: Xs, Ys, Z_div_Z_solar, Teff_jump, alfa, L1, M1, R1, T1, &
            vink_wind, nieu_wind, hamann_wind, lowT_w, highT_w, Twindow
         type (star_info), pointer :: s
         type (binary_info), pointer :: b
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call binary_ptr(s% binary_id, b, ierr)
         if (ierr /= 0) return

         L1 = Lsurf
         M1 = Msurf
         R1 = Rsurf
         T1 = Tsurf

         h1 = s% net_iso(ih1)
         he4 = s% net_iso(ihe4)
         Xs = s% xa(h1,1)
         Ys = s% xa(he4,1)
         ! Z=0.017 is Z from Grevesse et al. 1996
         Z_div_Z_solar = s% Zbase/0.0142d0
         ! use Vink et al 2001, eqns 14 and 15 to set "jump" temperature
         Teff_jump = 1d3*(61.2d0 + 2.59d0*(-13.636d0 + 0.889d0*log10_cr(Z_div_Z_solar)))

         vink_wind = 0d0
         nieu_wind = 0d0
         hamann_wind = 0d0
         w = 0

         call eval_Vink_wind(vink_wind)
         call eval_Nieuwenhuijzen_wind(nieu_wind)
         call eval_Hamann_wind(hamann_wind)

         ! use 1/10 hamann
         hamann_wind = hamann_wind/10d0

         lowT_w = max(vink_wind, nieu_wind)

         alfa = 0d0
         if (Xs > 0.7d0) then
            alfa = 1d0
         else if (Xs > 0.4d0 .and. Xs < 0.7d0) then
            alfa = (Xs - 0.4d0)/0.3d0
         end if
         highT_w = alfa * vink_wind + (1d0-alfa) * hamann_wind

         ! have a 10% Teff_jump window to switch from the lowT to the highT wind
         Twindow = Teff_jump*0.10d0
         alfa = 0d0
         if (T1 < Teff_jump - Twindow/2d0) then
            alfa = 1d0
         else if (T1 > Teff_jump - Twindow/2d0 .and. T1 < Teff_jump + Twindow/2d0) then
            alfa = ((Teff_jump + Twindow/2d0)-T1)/Twindow
         end if
         w = alfa * lowT_w + (1d0-alfa) * highT_w

         ! further soften change in wind to avoid things going bad
         if (s% xtra28 /= 0) then
            if(abs(w) > abs(s% xtra28)*1.05) then
               w = s% xtra28*1.05
            else if(abs(w) < abs(s% xtra28)*0.95) then
               w = s% xtra28*0.95
            end if
         end if
         s% xtra28 = w

         ierr = 0

         contains

         subroutine eval_Vink_wind(w)
            real(dp), intent(inout) :: w
            real(dp) :: alfa, w1, w2, logMdot, dT, vinf_div_vesc

            ! alfa = 1 for hot side, = 0 for cool side
            if (T1 > 27500d0) then
               alfa = 1
            else if (T1 < 22500d0) then
               alfa = 0
            else
               dT = 100d0
               if (T1 > Teff_jump + dT) then
                  alfa = 1
               else if (T1 < Teff_jump - dT) then
                  alfa = 0
               else
                  alfa = (T1 - (Teff_jump - dT)) / (2*dT)
               end if
            end if
            
            if (alfa > 0) then ! eval hot side wind (eqn 24)
               vinf_div_vesc = 2.6d0 ! this is the hot side galactic value
               vinf_div_vesc = vinf_div_vesc*pow_cr(Z_div_Z_solar,0.13d0) ! corrected for Z
               logMdot = &
                  - 6.697d0 &
                  + 2.194d0*log10_cr(L1/Lsun/1d5) &
                  - 1.313d0*log10_cr(M1/Msun/30) &
                  - 1.226d0*log10_cr(vinf_div_vesc/2d0) &
                  + 0.933d0*log10_cr(T1/4d4) &
                  - 10.92d0*pow2(log10_cr(T1/4d4)) &
                  + 0.85d0*log10_cr(Z_div_Z_solar)
               w1 = exp10_cr(logMdot)
            else
               w1 = 0
            end if
            
            if (alfa < 1) then ! eval cool side wind (eqn 25)
               vinf_div_vesc = 1.3d0 ! this is the cool side galactic value
               vinf_div_vesc = vinf_div_vesc*pow_cr(Z_div_Z_solar,0.13d0) ! corrected for Z
               logMdot = &
                  - 6.688d0 &
                  + 2.210d0*log10_cr(L1/Lsun/1d5) &
                  - 1.339d0*log10_cr(M1/Msun/30) &
                  - 1.601d0*log10_cr(vinf_div_vesc/2d0) &
                  + 1.07d0*log10_cr(T1/2d4) &
                  + 0.85d0*log10_cr(Z_div_Z_solar)
               w2 = exp10_cr(logMdot)
            else
               w2 = 0
            end if
            
            w = alfa*w1 + (1 - alfa)*w2
            
         end subroutine eval_Vink_wind

         subroutine eval_Nieuwenhuijzen_wind(w)
            ! Nieuwenhuijzen, H.; de Jager, C. 1990, A&A, 231, 134 (eqn 2)
            real(dp), intent(out) :: w
            real(dp) :: log10w
            include 'formats'
            log10w = -14.02d0 &
                     +1.24d0*log10_cr(L1/Lsun) &
                     +0.16d0*log10_cr(M1/Msun) &
                     +0.81d0*log10_cr(R1/Rsun) &
                     +0.85d0*log10_cr(Z_div_Z_solar)
            w = exp10_cr(log10w)
         end subroutine eval_Nieuwenhuijzen_wind

         subroutine eval_Hamann_wind(w)
            ! Hamann, W.-R.; Koesterke, L.; Wessolowski, U. 1995, A&A, 299, 151
            real(dp), intent(out) :: w
            real(dp) :: log10w
            include 'formats'
            log10w = -11.95d0 &
                     +1.5d0*log10_cr(L1/Lsun) &
                     -2.85d0*Xs &
                     + 0.85d0*log10_cr(Z_div_Z_solar)
            w = exp10_cr(log10w)
         end subroutine eval_Hamann_wind

      end subroutine brott_wind
      
      integer function extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: restart_time, prev_time_used
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         if (.not. restart) then
            call system_clock(time0,clock_rate)
            call alloc_restart_info(s)
            s% xtra25 = 0d0
            s% xtra26 = 0d0
            s% xtra27 = 0d0
            s% scale_max_correction = 1d99
         else
            call unpack_restart_info(s)
            call system_clock(restart_time,clock_rate)
            prev_time_used = time1 - time0
            time1 = restart_time
            time0 = time1 - prev_time_used
            s% scale_max_correction = 0.1d0
         end if
         extras_startup = keep_going
      end function extras_startup
      

      integer function extras_check_model(id, id_extra)
         use chem_def
         integer, intent(in) :: id, id_extra
         type (star_info), pointer :: s
         type (binary_info), pointer :: b
         integer :: ierr,i
         logical :: is_ne_biggest
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call binary_ptr(s% binary_id, b, ierr)
         if (ierr /= 0) return
         extras_check_model = keep_going

      end function extras_check_model


      integer function how_many_extra_history_columns(id, id_extra)
         integer, intent(in) :: id, id_extra
         integer :: ierr
         type (star_info), pointer :: s
         type (binary_info), pointer :: b
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call binary_ptr(s% binary_id, b, ierr)
         if (ierr /= 0) return
         how_many_extra_history_columns = 2
      end function how_many_extra_history_columns
      
      
      subroutine data_for_extra_history_columns(id, id_extra, n, names, vals, ierr)
         integer, intent(in) :: id, id_extra, n
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         real(dp) :: dt, mu, I1, I2, a_min, a, ratio
         integer :: k, i , indexR
         type (star_info), pointer :: s
         type (binary_info), pointer :: b

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call binary_ptr(s% binary_id, b, ierr)
         if (ierr /= 0) return

         !if (n /= 1) then
         !   stop 'bad n for data_for_extra_history_columns'
         !end if
         dt = dble(time1 - time0) / clock_rate / 60
         names(1) = 'runtime_minutes'
         vals(1) = dt

         names(2) = 'Separation_Ratio'
         vals(2) = 0d0

         mu = b% m(1) *b% m(2) / (b% m(1) + b% m(2))

         I1 = 0
         do k=1, s% nz
            I1 = I1 + 0.6666666667 * (s% dm(k)) * (s% r(k))**2
         end do

         I2 = 0

         a_min = sqrt( 3*(I1+I2)/mu )

         a = b% separation

         ratio = a_min/a
         vals(2) = ratio !ratio
         write(*,*) "ratio is", ratio

      end subroutine data_for_extra_history_columns

      
      integer function how_many_extra_profile_columns(id, id_extra)
         integer, intent(in) :: id, id_extra
         how_many_extra_profile_columns = 0
      end function how_many_extra_profile_columns
      
      
      subroutine data_for_extra_profile_columns(id, id_extra, n, nz, names, vals, ierr)
         integer, intent(in) :: id, id_extra, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         integer :: k
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine data_for_extra_profile_columns
  

      integer function extras_start_step(id, id_extra)
         integer, intent(in) :: id, id_extra
         integer :: ierr
         type (star_info), pointer :: s
         type (binary_info), pointer :: b

         write(*,*) "in extras start step"

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call binary_ptr(s% binary_id, b, ierr)
         if (ierr /= 0) return

         extras_start_step = keep_going
      end function extras_start_step
  
      integer function extras_finish_step(id, id_extra)
         integer, intent(in) :: id, id_extra
         integer :: ierr, i, indexR
         real(dp) :: time, dt, stop_time, I1, I2, mu, a_min, a, ratio, radius, &
                     F3, G1, dP, q, rl3, q_temp, L3, A1, B1, start, q_old, F1_exact_1,&
                     Phi, dA, unit, Phi_last, F1_old, F3_old, G1_old, dP_old, mdot_thick,&
                     Phi_1, Phi_2, Phi_3, linear, C1, D1, E1, L1, L2, L
         type (star_info), pointer :: s
         type (binary_info), pointer :: b
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call binary_ptr(s% binary_id, b, ierr)
         if (ierr /= 0) return
         call system_clock(time1,clock_rate)
         call store_restart_info(s)

         if(.not. b% ignore_rlof_flag .and. abs(s% xtra25_old - s% center_h1) > 0.002) then
             s% dt_next = min(s% dt_next, s% dt * sqrt(s% min_timestep_factor))
             write(*,*) "reducing dt due to large change in central hydrogen"
         else if(.not. b% ignore_rlof_flag .and. abs(s% xtra26_old - s% center_he4) > 0.002) then
             s% dt_next = min(s% dt_next, s% dt * sqrt(s% min_timestep_factor))
             write(*,*) "reducing dt due to large change in central helium"
         else if(.not. b% ignore_rlof_flag .and. abs(s% xtra27_old - s% center_c12) > 0.002) then
             s% dt_next = min(s% dt_next, s% dt * sqrt(s% min_timestep_factor))
             write(*,*) "reducing dt due to large change in central carbon"
         end if

         s% xtra25 = s% center_h1
         s% xtra26 = s% center_he4
         s% xtra27 = s% center_c12

	 dt = (s% dt)/(365.25*24*60*60)
	 dt = log10(dt)
	 write(*,*) "dt =", dt

	 if (dt < -3.9428) then
	    write(*,*) "Reached dt minimum of an hour, terminating"
	    extras_finish_step = terminate
	    return
	 end if

	 time = dble(time1 - time0) / clock_rate / 60
	 write(*,*) "time passed is = ", time

	 stop_time = 1380
	 if (time > stop_time) then
	    write(*,*) "Reached runtime max of", stop_time ,"minutes, terminating"
            extras_finish_step = terminate
            return
         end if

         extras_finish_step = keep_going
      end function extras_finish_step
      
      
      subroutine extras_after_evolve(id, id_extra, ierr)
         integer, intent(in) :: id, id_extra
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         real(dp) :: dt
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         dt = dble(time1 - time0) / clock_rate / 60
         write(*,'(/,a50,2f18.6,99i10/)') 'runtime, retries, backups, steps', &
            dt, expected_runtime, s% num_retries, s% num_backups, s% model_number
      end subroutine extras_after_evolve
      
      
      ! routines for saving and restoring data so can do restarts

      
      subroutine alloc_restart_info(s)
         type (star_info), pointer :: s
         call move_restart_info(s,restart_info_alloc)
      end subroutine alloc_restart_info
      
      
      subroutine unpack_restart_info(s)
         type (star_info), pointer :: s
         call move_restart_info(s,restart_info_get)
      end subroutine unpack_restart_info
      
      
      subroutine store_restart_info(s)
         type (star_info), pointer :: s
         call move_restart_info(s,restart_info_put)
      end subroutine store_restart_info
      
      
      subroutine move_restart_info(s,op)
         type (star_info), pointer :: s
         integer, intent(in) :: op
         
         integer :: i, j, num_ints, num_dbls, ierr
         
         i = 0
         ! call move_int or move_flg 
         call move_int(time0)
         call move_int(time1)
         
         num_ints = i
         
         i = 0
         ! call move_dbl 
         
         num_dbls = i
         
         if (op /= restart_info_alloc) return
         if (num_ints == 0 .and. num_dbls == 0) return
         
         ierr = 0
         call star_alloc_extras(s% id, num_ints, num_dbls, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in star_alloc_extras'
            write(*,*) 'alloc_extras num_ints', num_ints
            write(*,*) 'alloc_extras num_dbls', num_dbls
            stop 1
         end if
         
         contains
         
         subroutine move_dbl(dbl)
            real(dp) :: dbl
            i = i+1
            select case (op)
            case (restart_info_get)
               dbl = s% extra_work(i)
            case (restart_info_put)
               s% extra_work(i) = dbl
            end select
         end subroutine move_dbl
         
         subroutine move_int(int)
            integer :: int
            include 'formats'
            i = i+1
            select case (op)
            case (restart_info_get)
               !write(*,3) 'restore int', i, s% extra_iwork(i)
               int = s% extra_iwork(i)
            case (restart_info_put)
               !write(*,3) 'save int', i, int
               s% extra_iwork(i) = int
            end select
         end subroutine move_int
         
         subroutine move_flg(flg)
            logical :: flg
            i = i+1
            select case (op)
            case (restart_info_get)
               flg = (s% extra_iwork(i) /= 0)
            case (restart_info_put)
               if (flg) then
                  s% extra_iwork(i) = 1
               else
                  s% extra_iwork(i) = 0
               end if
            end select
         end subroutine move_flg

      end subroutine move_restart_info
      
      


      end module run_star_extras
      

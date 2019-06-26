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
      module run_binary_extras 

      use star_lib
      use star_def
      use const_def
      use const_def
      use chem_def
      use num_lib
      use binary_def
      use binary_lib
      use crlibm_lib
      use utils_lib
      
      implicit none
      
      contains
      
      subroutine extras_binary_controls(binary_id, ierr)
         integer :: binary_id
         integer, intent(out) :: ierr
         type (binary_info), pointer :: b
         ierr = 0

         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if

         ! Set these function pinters to point to the functions you wish to use in
         ! your run_binary_extras. Any which are not set, default to a null_ version
         ! which does nothing.
         b% how_many_extra_binary_history_columns => how_many_extra_binary_history_columns
         b% data_for_extra_binary_history_columns => data_for_extra_binary_history_columns

         b% extras_binary_startup=> extras_binary_startup
         b% extras_binary_start_step=> extras_binary_start_step
         b% extras_binary_check_model=> extras_binary_check_model
         b% extras_binary_finish_step => extras_binary_finish_step
         b% extras_binary_after_evolve=> extras_binary_after_evolve
         
         b% other_rlo_mdot => my_rlo_mdot
         b% other_jdot_ml => my_jdot_ml
         b% use_other_rlo_mdot = .true.
         b% use_other_jdot_ml = .true.

         ! Once you have set the function pointers you want, then uncomment this (or set it in your star_job inlist)
         ! to disable the printed warning message,
          b% warn_binary_extra =.false.
         
      end subroutine extras_binary_controls

      subroutine my_jdot_ml(binary_id, ierr)
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr
         type (binary_info), pointer :: b
         real(dp) :: frac_L1, frac_outer, real_mdot, x_Louter, j1, j2
         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if

     !momentum gained from MT from winds
         !mass lost from vicinity of donor
         b% jdot_ml = (b% mdot_system_wind(b% d_i))*&
             (b% m(b% a_i)/(b% m(b% a_i)+b% m(b% d_i))*b% separation)**2*2*pi/b% period *&
             sqrt(1 - b% eccentricity**2)
         !mass lost from vicinity of accretor
         b% jdot_ml = b% jdot_ml + (b% mdot_system_wind(b% a_i))*&
             (b% m(b% d_i)/(b% m(b% a_i)+b% m(b% d_i))*b% separation)**2*2*pi/b% period *&
             sqrt(1 - b% eccentricity**2)

     !New part added to account for outer MT     
          call my_rlo_mdot(b% binary_id, real_mdot, ierr)

     !Separating MT from inner and outer
          frac_L1 = (b% mdot_thin+b% mdot_thick)/real_mdot
          frac_outer = 1-frac_L1

          write(*,*) "real vs guess MT", real_mdot/Msun*secyer, b% mtransfer_rate/Msun*secyer, &
             b% mdot_system_transfer(b% a_i)/Msun*secyer 
          write(*,*) "fractions of MT, inner vs outer", frac_L1, frac_outer

     !Picking L2 or L3 depending on q
          if(b% m(b% d_i)/b% m(b% a_i) > 1) then
             x_Louter = abs(find_L3(b))
          else
             x_Louter = abs(find_L2(b))
          end if

         if (abs(b% mtransfer_rate) > 0d0) then
            if (frac_outer > 0d0) then
              !mass lost from vicinity of donor, only from outer Lagrangian
               b% jdot_ml = b% jdot_ml + (b% mtransfer_rate*frac_outer)*&
                   ((b% m(b% a_i)/(b% m(b% a_i)+b% m(b% d_i))+x_Louter)*b% separation)**2*2*pi/b% period *&
                   sqrt(1 - b% eccentricity**2)
            end if
            if (frac_L1 > 0d0) then
              !mass lost from vicinity of accretor, only from L1 point
               b% jdot_ml = b% jdot_ml + (b% mtransfer_rate*frac_L1)*&
                   (b% m(b% d_i)/(b% m(b% a_i)+b% m(b% d_i))*b% separation)**2*2*pi/b% period *&
                   sqrt(1 - b% eccentricity**2)
            end if
          end if

      !Same as before but unused, only for error analysis
            !mass lost from vicinity of donor
            j1 = (b% mtransfer_rate*frac_outer)*&
                ((b% m(b% a_i)/(b% m(b% a_i)+b% m(b% d_i))+x_Louter)*b% separation)**2*2*pi/b% period *&
                sqrt(1 - b% eccentricity**2)
            !mass lost from vicinity of accretor
            j2 = (b% mtransfer_rate*frac_L1)*&
                (b% m(b% d_i)/(b% m(b% a_i)+b% m(b% d_i))*b% separation)**2*2*pi/b% period *&
                sqrt(1 - b% eccentricity**2)

         write(*,*) "Outer L position, and cm position", x_Louter, (b% m(b% a_i)/(b% m(b% a_i)+b% m(b% d_i)))
         write(*,*) "check jdot outer, jdot inner, and total jdot", j1, j2, b% jdot_ml
         write(*,*) "check sjdot outer, sjdot inner", j1/frac_outer, j2/frac_L1
      end subroutine my_jdot_ml

    subroutine my_rlo_mdot(binary_id, mdot, ierr)
         use const_def, only: dp
         integer, intent(in) :: binary_id
         real(dp), intent(out) :: mdot
         real(dp) :: mdot_outer
         integer, intent(out) :: ierr
         type (binary_info), pointer :: b
         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if

         mdot = 0d0

         if (b% mdot_scheme == "roche_lobe") then
            write(*,*) "mdot_scheme = roche_lobe not applicable for explicit scheme"
            write(*,*) "Not transfering mass"
            mdot = 0
            return
         else if (b% mdot_scheme /= "Ritter" .and. b% mdot_scheme /= "Kolb" .and. b% mdot_scheme /= "Arras") then
            write(*,*) "mdot_scheme = " , b% mdot_scheme , " not recognized"
            write(*,*) "Not transfering mass"
            mdot = 0
            return
         end if

         if (b% mdot_scheme == "Kolb" .and. b% eccentricity <= 0.0) then
            call get_info_for_ritter(b)
            mdot = b% mdot_thin
            call get_info_for_kolb(b)
            mdot = mdot + b% mdot_thick
            mdot_outer = get_info_for_kolb_outer(b)
            write(*,*) "check mdots", mdot/Msun*secyer, mdot_outer/Msun*secyer
            mdot = mdot + mdot_outer
            
         else if (b% mdot_scheme == "Kolb" .and. b% eccentricity > 0.0) then
            call get_info_for_ritter_eccentric(b)
            mdot = b% mdot_thin
            call get_info_for_kolb_eccentric(b)
            mdot = mdot + b% mdot_thick
            
         else if (b% mdot_scheme == "Ritter" .and. b% eccentricity <= 0.0) then
            call get_info_for_ritter(b)
            mdot = b% mdot_thin
            
         else if (b% mdot_scheme == "Ritter" .and. b% eccentricity > 0.0) then
            call get_info_for_ritter_eccentric(b)
            mdot = b% mdot_thin

         end if
            
         if (b% mdot_scheme == "Arras") then
            if (b% eccentricity > 0d0) &
               write(*,*) "mdot_scheme = Arras is not properly implemented for e>0"
            call get_info_for_arras(b)
            mdot = b% mdot_thin

         end if

      end subroutine my_rlo_mdot
         
      subroutine CE_rlo_mdot(binary_id, rlo_mdot, ierr)
         use const_def, only: dp
         integer, intent(in) :: binary_id
         real(dp), intent(out) :: rlo_mdot
         integer, intent(out) :: ierr
         type (binary_info), pointer :: b
         real(dp) :: exp_factor

         ierr = 0

         call binary_ptr(binary_id, b, ierr)

         write(*,*) "in CE_RLO"

         if (b% use_other_CE_rlo_mdot) then
            call b% other_CE_rlo_mdot(b% binary_id, rlo_mdot, ierr)
            return
         end if

         exp_factor = -log_cr(b% CE_mass_loss_rate_low/b% CE_mass_loss_rate_high)

         if (b% r(b% d_i)-b% rl(b% d_i) > 0d0) then
            rlo_mdot = -Msun/secyer*b% CE_mass_loss_rate_high
         else if (b% r(b% d_i)-b% rl(b% d_i) < -b% CE_rel_rlo_for_detachment*b% rl(b% d_i)) then
            rlo_mdot = -Msun/secyer*b% CE_mass_loss_rate_low
         else
            rlo_mdot = -Msun/secyer*b% CE_mass_loss_rate_high * &
               exp_cr(exp_factor*(b% r(b% d_i)-b% rl(b% d_i))/(b% rl(b% d_i)*b% CE_rel_rlo_for_detachment))
         end if
      end subroutine CE_rlo_mdot

      integer function how_many_extra_binary_history_columns(binary_id)
         use binary_def, only: binary_info
         integer, intent(in) :: binary_id
         how_many_extra_binary_history_columns = 6
      end function how_many_extra_binary_history_columns
      
      subroutine data_for_extra_binary_history_columns(binary_id, n, names, vals, ierr)
         use const_def, only: dp
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer, intent(in) :: n
         character (len=maxlen_binary_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         real(dp) :: beta, radius, mdot_outer, q, mdot_L1
         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if
                  
         names(1) = 'Outer_MT_flag'
         vals(1) = 0

         names(2) = 'Outer_MT_thick'
         vals(2) = 0d0

         names(3) = 'L1_MT_flag'
         vals(3) = 0

         names(4) = 'L1_MT_thick'
         vals(4) = 0d0

         names(5) = 'L1_MT_thin_flag'
         vals(5) = 0

         names(6) = 'L1_MT_thin'
         vals(6) = 0d0

         q = b% m(b% d_i)/b% m(b% a_i)

         !radius = 0.422d0 + 0.298d0*atan(1.06d0*log10_cr(q)+0.329d0)
         !radius = radius*b% separation
         radius = outer_radius(b)
         write(*,*) "radius" , radius
         if (.not. b% CE_flag) then

            call get_info_for_ritter(b)
            vals(6) = b% mdot_thin
     
            if (vals(6) .NE. 0d0) then
               vals(5) = 1
            else
               vals(5) = 0
            end if
     
           write(*,*) "thin flag" , vals(5)
           write(*,*) "thin value", vals(6)
 
           if(b% r(b% d_i)- b% rl(b% d_i) > 0.0d0) then
           vals(3) = 1
           else
              vals(3) = 0
           end if

           call get_info_for_kolb(b)
           mdot_L1 = b% mdot_thick
           vals(4) = mdot_L1/Msun*secyer

            if(b% r(b% d_i)- radius > 0.0d0) then
               vals(1) = 1
            else
               vals(1) = 0
            end if

            mdot_outer = get_info_for_kolb_outer(b)
            vals(2) = mdot_outer/Msun*secyer

            write(*,*) "Outer MT Flag is", vals(1)
            write(*,*) "Outer MT thick is", vals(2)

         end if

      end subroutine data_for_extra_binary_history_columns
      
      
      integer function extras_binary_startup(binary_id,restart,ierr)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr
         logical, intent(in) :: restart
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if
         
         !if (.not. restart) then
         !end if

         !if (restart) then
         !end if
         
         extras_binary_startup = keep_going
      end function  extras_binary_startup
      
      integer function extras_binary_start_step(binary_id,ierr)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr
         real(dp) :: center_c12, center_he4
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if
         
         extras_binary_start_step = keep_going

         if (.not. b% CE_flag) then
            b% other_rlo_mdot => my_rlo_mdot
            b% s_donor% mix_factor = 1d0
            b% s_donor% dxdt_nuc_factor = 1d0
            if (b% point_mass_i == 0) then
               b% s_accretor% mix_factor = 1d0
               b% s_accretor% dxdt_nuc_factor = 1d0
            end if
            b% keep_donor_fixed = .false.
         else
            b% other_rlo_mdot => CE_rlo_mdot
            b% s_donor% mix_factor = 1d0
            b% s_donor% dxdt_nuc_factor = 1d0
            if (b% point_mass_i == 0) then
               b% s_accretor% mix_factor = 1d0
               b% s_accretor% dxdt_nuc_factor = 1d0
            end if
            b% keep_donor_fixed = .true.
         end if


         if (b% point_mass_i /= 2) then
            if (b% CE_flag) then
               b% s2% delta_HR_limit = -1d0
               b% s2% delta_HR_hard_limit = -1d0
            else
               b% s2% delta_HR_limit = 0.005d0
               b% s2% delta_HR_hard_limit = 0.01d0
            end if
            if (b% s2% X(b% s2% nz) > 1d-4) then
               b% s2% burn_he_mix_region_logT = 1d10
            else
               b% s2% burn_he_mix_region_logT = 7.7d0
            end if
            b% s2% scale_max_correction = 0.1d0
            b% s2% relax_hard_limits_after_retry = .false.
            b% s2% relax_hard_limits_after_backup = .false.
            !if (abs(b% s2% mstar_dot) > 0d0) then
            !   b% s2% max_timestep = 1d-2*Msun/abs(b% s2% mstar_dot)
            !   !write(*,*) "check max timestep", b% s2% max_timestep/secyer
            !end if
            center_c12 = b% s2% xa(b% s2% net_iso(ic12),b% s2% nz)
            center_he4 = b% s2% xa(b% s2% net_iso(ihe4),b% s2% nz)
            !if (center_he4 < 1d-2) then
            if (center_he4 < 1d-2 .and. center_c12 < 1d-2) then
               if (b% point_mass_i == 0) then
                  write(*,*) "Secondary has depleted central carbon"
                  call binary_set_point_mass_i(b% binary_id, 2, ierr)
                  write(*,*) "check masses", b% m(1)/Msun, b% m(2)/Msun
                  b% limit_retention_by_mdot_edd = .true.
                  b% eq_initial_bh_mass = b% s2% m(1)
                  if (ierr /= 0) then
                     extras_binary_start_step = terminate
                     return
                  end if
               else
                  extras_binary_start_step = terminate
                  write(*,*) "Terminate due to secondary depleting carbon"
               end if
            end if
         end if

         if (b% point_mass_i /= 1) then
            if (b% CE_flag) then
               b% s1% delta_HR_limit = -1d0
               b% s1% delta_HR_hard_limit = -1d0
            else
               b% s1% delta_HR_limit = 0.005d0
               b% s1% delta_HR_hard_limit = 0.01d0
            end if
            if (b% s1% X(b% s1% nz) > 1d-4) then
               b% s1% burn_he_mix_region_logT = 1d10
            else
               b% s1% burn_he_mix_region_logT = 7.7d0
            end if
            b% s1% scale_max_correction = 0.1d0
            b% s1% relax_hard_limits_after_retry = .false.
            b% s1% relax_hard_limits_after_backup = .false.
            !if (abs(b% s1% mstar_dot) > 0d0) then
            !   b% s1% max_timestep = 1d-2*Msun/abs(b% s1% mstar_dot)
            !   !write(*,*) "check max timestep", b% s1% max_timestep/secyer
            !end if
            center_c12 = b% s1% xa(b% s1% net_iso(ic12),b% s1% nz)
            center_he4 = b% s1% xa(b% s1% net_iso(ihe4),b% s1% nz)
            !if (center_he4 < 1d-2) then
            if (center_he4 < 1d-2 .and. center_c12 < 1d-2) then
               if (b% point_mass_i == 0) then
                  write(*,*) "Primary has depleted central carbon"
                  call binary_set_point_mass_i(b% binary_id, 1, ierr)
                  b% limit_retention_by_mdot_edd = .true.
                  b% eq_initial_bh_mass = b% s1% m(1)
                  write(*,*) "check masses", b% m(1)/Msun, b% m(2)/Msun
                  if (ierr /= 0) return
               else
                  extras_binary_start_step = terminate
                  write(*,*) "Terminate due to primary depleting carbon"
               end if
            end if
         end if

         if (b% point_mass_i == 0) then
            b% limit_retention_by_mdot_edd = .false.
         else
            b% limit_retention_by_mdot_edd = .true.
         end if

      end function  extras_binary_start_step
      
      !Return either rety,backup,keep_going or terminate
      integer function extras_binary_check_model(binary_id)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer :: ierr
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if  
         extras_binary_check_model = keep_going
        
      end function extras_binary_check_model
      
      
      ! returns either keep_going or terminate.
      ! note: cannot request retry or backup; extras_check_model can do that.
      integer function extras_binary_finish_step(binary_id)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer :: ierr, star_id, i
         real(dp) :: q, mdot_limit_low, mdot_limit_high, &
            center_h1, center_h1_old, center_he4, center_he4_old
         logical :: is_ne_biggest

         extras_binary_finish_step = keep_going

         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if

         if (extras_binary_finish_step == terminate) then
            if (b% have_star_1) then
               call star_write_profile_info(b% s1% id, "LOGS1/prof_9FINAL.data", b% s1% id, ierr)
               if (ierr /= 0) return ! failure in profile
            end if
            if (b% have_star_2) then
               call star_write_profile_info(b% s2% id, "LOGS2/prof_9FINAL.data", b% s2% id, ierr)
               if (ierr /= 0) return ! failure in profile
            end if
         else
            if (b% have_star_1) then
               !additional profiles to be saved
               center_h1 = b% s1% xa(b% s1% net_iso(ih1),b% s1% nz)
               center_h1_old = b% s1% xa_old(b% s1% net_iso(ih1),b% s1% nz_old)
               center_he4 = b% s1% xa(b% s1% net_iso(ihe4),b% s1% nz)
               center_he4_old = b% s1% xa_old(b% s1% net_iso(ihe4), b% s1% nz_old)
               if (center_h1 < 0.5 .and. center_h1_old > 0.5) then
                  call star_write_profile_info(b% s1% id, "LOGS1/prof_2H50.data", b% s1% id, ierr)
                  if (ierr /= 0) return ! failure in profile
               else if (center_h1 < 0.25 .and. center_h1_old > 0.25) then
                  call star_write_profile_info(b% s1% id, "LOGS1/prof_3H25.data", b% s1% id, ierr)
                  if (ierr /= 0) return ! failure in profile
               else if (center_h1 < 1d-6 .and. center_h1_old > 1d-6) then
                  call star_write_profile_info(b% s1% id, "LOGS1/prof_4H00.data", b% s1% id, ierr)
                  if (ierr /= 0) return ! failure in profile
               else if (center_h1 < 1d-6) then
                  if (center_he4 < 0.75 .and. center_he4_old > 0.75) then
                     call star_write_profile_info(b% s1% id, "LOGS1/prof_5He75.data", b% s1% id, ierr)
                     if (ierr /= 0) return ! failure in profile
                  else if (center_he4 < 0.5 .and. center_he4_old > 0.5) then
                     call star_write_profile_info(b% s1% id, "LOGS1/prof_6He50.data", b% s1% id, ierr)
                     if (ierr /= 0) return ! failure in profile
                  else if (center_he4 < 0.25 .and. center_he4_old > 0.25) then
                     call star_write_profile_info(b% s1% id, "LOGS1/prof_7He25.data", b% s1% id, ierr)
                     if (ierr /= 0) return ! failure in profile
                  else if (center_he4 < 1d-5 .and. center_he4_old > 1d-5) then
                     !create the profile at Y<1d-5, to avoid profiles to be created every single
                     !timestep if we only evolve until helium depletion, defined at Y=1d-6
                     call star_write_profile_info(b% s1% id, "LOGS1/prof_8He00.data", b% s1% id, ierr)
                     if (ierr /= 0) return ! failure in profile
                  end if
               end if
            end if

            if (b% have_star_2) then
               center_h1 = b% s2% xa(b% s2% net_iso(ih1),b% s2% nz)
               center_h1_old = b% s2% xa_old(b% s2% net_iso(ih1),b% s2% nz_old)
               center_he4 = b% s2% xa(b% s2% net_iso(ihe4), b% s2% nz)
               center_he4_old = b% s2% xa_old(b% s2% net_iso(ihe4),b% s2% nz_old)
               if (center_h1 < 0.5 .and. center_h1_old > 0.5) then
                  call star_write_profile_info(b% s2% id, "LOGS2/prof_2H50.data", b% s2% id, ierr)
                  if (ierr /= 0) return ! failure in profile
               else if (center_h1 < 0.25 .and. center_h1_old > 0.25) then
                  call star_write_profile_info(b% s2% id, "LOGS2/prof_3H25.data", b% s2% id, ierr)
                  if (ierr /= 0) return ! failure in profile
               else if (center_h1 < 1d-6 .and. center_h1_old > 1d-6) then
                  call star_write_profile_info(b% s2% id, "LOGS2/prof_4H00.data", b% s2% id, ierr)
                  if (ierr /= 0) return ! failure in profile
               else if (center_h1 < 1d-6) then
                  if (center_he4 < 0.75 .and. center_he4_old > 0.75) then
                     call star_write_profile_info(b% s2% id, "LOGS2/prof_5He75.data", b% s2% id, ierr)
                     if (ierr /= 0) return ! failure in profile
                  else if (center_he4 < 0.5 .and. center_he4_old > 0.5) then
                     call star_write_profile_info(b% s2% id, "LOGS2/prof_6He50.data", b% s2% id, ierr)
                     if (ierr /= 0) return ! failure in profile
                  else if (center_he4 < 0.25 .and. center_he4_old > 0.25) then
                     call star_write_profile_info(b% s2% id, "LOGS2/prof_7He25.data", b% s2% id, ierr)
                     if (ierr /= 0) return ! failure in profile
                  else if (center_he4 < 1d-6 .and. center_he4_old > 1d-6) then
                     call star_write_profile_info(b% s2% id, "LOGS2/prof_8He00.data", b% s2% id, ierr)
                     if (ierr /= 0) return ! failure in profile
                  end if
               end if
            end if
         end if
         
      end function extras_binary_finish_step
      
      real(dp) function eval_rlobe(m1, m2, a) result(rlobe)
         real(dp), intent(in) :: m1, m2, a
         real(dp) :: q
         q = pow_cr(m1/m2,one_third)
      ! Roche lobe size for star of mass m1 with a
      ! companion of mass m2 at separation a, according to
      ! the approximation of Eggleton 1983, apj 268:368-369
         rlobe = a*0.49d0*q*q/(0.6d0*q*q + log1p_cr(q))
      end function eval_rlobe
      
      subroutine extras_binary_after_evolve(binary_id, ierr)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr
         integer :: iounit
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if      
 
      end subroutine extras_binary_after_evolve     

      !This is just copied from binary_mdot

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

      real(dp) function dPhidy(b,x,y,z) result(dpdy)
         real(dp), intent(in) :: x, y, z
         real(dp) :: q, x_cm, r, r2
         type(binary_info), pointer :: b
         include 'formats.inc'
         
         q = b% m(b% d_i)/b% m(b% a_i)
         r = sqrt(x*x+y*y+z*z)
         r2 = sqrt(((1-x)*(1-x))+y*y+z*z)
         dpdy = q*y/((1+q)*r*r*r) + y/((1+q)*r2*r2*r2) - y
      end function dPhidy

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
         limit = abs(upper_bound-lower_bound)/upper_bound
         tolerance = 0.000001d0
  
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
            limit = abs(upper_bound-lower_bound)/upper_bound
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
            limit = abs(upper_bound-lower_bound)/abs(lower_bound)
         end if

         if (q .GE. 1) then
            upper_bound = 2d0
            lower_bound = 1d0
            limit = abs(upper_bound-lower_bound)/abs(upper_bound)
         end if

         x = 0d0

         tolerance = 0.000001d0

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
 
            if (q < 1) then
               limit = abs(upper_bound-lower_bound)/abs(lower_bound)
            end if
 
            if (q .GE. 1) then
               limit = abs(upper_bound-lower_bound)/abs(upper_bound)
            end if

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
            limit = abs(upper_bound-lower_bound)/abs(lower_bound)
         end if

         if (q < 1) then
            upper_bound = 2d0
            lower_bound = 1d0
            limit = abs(upper_bound-lower_bound)/upper_bound
         end if

         x = 0d0

         tolerance = 0.000001d0

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

         if (q .GE. 1) then
            limit = abs(upper_bound-lower_bound)/abs(lower_bound)
         end if

         if (q < 1) then
            limit = abs(upper_bound-lower_bound)/upper_bound
         end if

         end do
         L3 = (upper_bound + lower_bound)/2
      end function find_L3

      real(dp) function outer_radius(b) result(radius)
         real(dp) :: q, a, b1, c, d, e, f, g, sigma, arg
         type(binary_info), pointer :: b
         include 'formats.inc'

         q = b% m(b% d_i)/b% m(b% a_i)
         q = log10_cr(q)
         a = 2.7412d0
         b1 = 0.4426d0
         c = 21.4669d0
         d = -0.478d0
         e = 0.8894d0
         f = 7.1311d0
         g = 12.1788d0
         sigma = c/(g+exp_cr(d*q))
         arg = ((q-b1)/sigma)*((q-b1)/sigma)
         radius = a/(1+arg)/(f+exp_cr(e*q))+1
         radius = b% rl(b% d_i) * radius

      end function outer_radius

      subroutine get_info_for_ritter(b)
         type(binary_info), pointer :: b
         real(dp) :: rho_exponent, F1, q, rho, p, grav, hp, v_th, rl3, q_temp
         include 'formats.inc'

         !--------------------- Optically thin MT rate -----------------------------------------------
         ! As described in H. Ritter 1988, A&A 202,93-100 and U. Kolb and H. Ritter 1990, A&A 236,385-392

         rho = b% s_donor% rho(1) ! density at surface in g/cm^3
         p = b% s_donor% p(1) ! pressure at surface in dynes/cm^2
         grav = b% s_donor% cgrav(1)*b% m(b% d_i)/(b% r(b% d_i))**2 ! local gravitational acceleration
         hp = p/(grav*rho) ! pressure scale height
         v_th = sqrt(kerg * b% s_donor% T(1) / (mp * b% s_donor% mu(1)))

         q = b% m(b% a_i)/b% m(b% d_i) ! Mass ratio, as defined in Ritter 1988
                                       ! (Kolb & Ritter 1990 use the opposite!)
         ! consider range of validity for F1, do not extrapolate! Eq. A9 of Ritter 1988
         q_temp = min(max(q,0.5d0),10d0)
         F1 = (1.23d0  + 0.5D0* log10_cr(q_temp))
         rl3 = (b% rl(b% d_i))*(b% rl(b% d_i))*(b% rl(b% d_i))
         b% mdot_thin0 = (2.0D0*pi/exp_cr(0.5d0)) * v_th*v_th*v_th * &
             rl3/(b% s_donor% cgrav(1)*b% m(b% d_i)) * rho * F1
         !Once again, do not extrapolate! Eq. (7) of Ritter 1988
         q_temp = min(max(q,0.04d0),20d0)
         if (q_temp < 1.0d0) then
            b% ritter_h = hp/( 0.954D0 + 0.025D0*log10_cr(q_temp) - 0.038D0*(log10_cr(q_temp))**2 )
         else
            b% ritter_h = hp/( 0.954D0 + 0.039D0*log10_cr(q_temp) + 0.114D0*(log10_cr(q_temp))**2 )
         end if

         b% ritter_exponent = (b% r(b% d_i)-b% rl(b% d_i))/b% ritter_h

         if (b% mdot_scheme == "Kolb") then
            if (b% ritter_exponent > 0) then
               b% mdot_thin = -b% mdot_thin0
            else
               b% mdot_thin = -b% mdot_thin0 * exp_cr(b% ritter_exponent)
            end if
         else
            b% mdot_thin = -b% mdot_thin0 * exp_cr(b% ritter_exponent)
         end if

      end subroutine get_info_for_ritter

      real(dp) function calculate_kolb_mdot_thick(b, indexR, rl_d) result(mdot_thick)
         real(dp), intent(in) :: rl_d
         integer, intent(in) :: indexR
         real(dp) :: F1, F3, G1, dP, q, rl3, q_temp, L1, A1, B1, start, q_old, F1_exact_1,&
                     Phi, dA, unit, Phi_last, F1_old, F3_old, G1_old, dP_old, mdot_thick_old,&
                     mdot_thick_alt, C1, D1, E1, linear, F1_exact_2, g_q1, g_q2, quad
         integer :: i
         type(binary_info), pointer :: b
         include 'formats.inc'

         !--------------------- Optically thick MT rate -----------------------------------------------
         
        ! write (*,*) "Starting optically thick mass transfer"

         L1 = find_L1(b)
         A1 = ddPhidy(b,L1, 0d0, 0d0)/2
         B1 = ddPhidz(b,L1, 0d0, 0d0)/2

         start = pi/(sqrt(A1*B1))

         C1 = ddddPhidy(b,L1, 0d0, 0d0)/24
         D1 = ddddPhidz(b,L1, 0d0, 0d0)/24
         E1 = ddddPhidyz(b,L1, 0d0, 0d0)/4

         linear = -2*pi*(3*D1*A1*A1 + A1*B1*E1 + 3*C1*B1*B1)/&
                    (8*sqrt(A1*A1*A1*A1*A1)*sqrt(B1*B1*B1*B1*B1))

	 mdot_thick = 0d0
         Phi = 0d0
         do i=indexR-1, 1, -1
            Phi = Phi + (b% s_donor% P(i+1)-b% s_donor% P(i))/(b% s_donor% rho(i))
            G1 = b% s_donor% gamma1(i)
            F3 = sqrt(G1) * pow_cr(2d0/(G1+1d0), (G1+1d0)/(2d0*G1-2d0))
            unit = b% separation/(b% s_donor% cgrav(i)*(b% m(b% a_i) + b% m(b% d_i)))
            dA = (b% separation * b% separation * b% separation)/&
               (b% s_donor% cgrav(i)*(b% m(b% a_i) + b% m(b% d_i)))*(start + linear*unit*Phi) 
            mdot_thick = mdot_thick + dA*F3*sqrt(kerg * b% s_donor% T(i) / &
               (mp * b% s_donor% mu(i)))*(b% s_donor% P(i+1)-b% s_donor% P(i))
         end do

         i = indexR
         ! only take a fraction of dP for last cell
         G1 = b% s_donor% gamma1(i)
         F3 = sqrt(G1) * pow_cr(2d0/(G1+1d0), (G1+1d0)/(2d0*G1-2d0))
         dP = (b% s_donor% r(indexR) - rl_d) / &
            (b% s_donor% r(indexR) - b% s_donor% r(indexR+1)) * (b% s_donor% P(i+1)-b% s_donor% P(i))
         Phi_last = dP/(b% s_donor% rho(i))
         unit = b% separation/(b% s_donor% cgrav(i)*(b% m(b% a_i) + b% m(b% d_i)))
         dA = (b% separation * b% separation * b% separation)/&
               (b% s_donor% cgrav(i)*(b% m(b% a_i) + b% m(b% d_i)))*(start + linear*unit*Phi_last)
         mdot_thick = mdot_thick + dA*F3*sqrt(kerg * b% s_donor% T(i) / (mp*b% s_donor% mu(i)))*dP

         mdot_thick = - mdot_thick
         write(*,*) "New calc mdot thick is", mdot_thick/Msun*secyer

     !Previous method of calculation using approximates.
         ! As described in Kolb and H. Ritter 1990, A&A 236,385-392

         ! compute integral in Eq. (A17 of Kolb & Ritter 1990)
         mdot_thick_old = 0d0
         do i=1,indexR-1
            G1_old = b% s_donor% gamma1(i)
            F3_old = sqrt(G1_old) * pow_cr(2d0/(G1_old+1d0), (G1_old+1d0)/(2d0*G1_old-2d0))
            mdot_thick_old = mdot_thick_old + F3_old*sqrt(kerg * b% s_donor% T(i) / &
               (mp * b% s_donor% mu(i)))*(b% s_donor% P(i+1)-b% s_donor% P(i))
         end do
         ! only take a fraction of dP for last cell 
         G1_old = b% s_donor% gamma1(i)
         F3_old = sqrt(G1_old) * pow_cr(2d0/(G1_old+1d0), (G1_old+1d0)/(2d0*G1_old-2d0))
         dP_old = (b% s_donor% r(indexR) - rl_d) / &
            (b% s_donor% r(indexR) - b% s_donor% r(indexR+1)) * (b% s_donor% P(i+1)-b% s_donor% P(i))
         mdot_thick_old = mdot_thick_old + F3_old*sqrt(kerg * b% s_donor% T(i) / (mp*b% s_donor% mu(i)))*dP_old

         q_old = b% m(b% a_i)/b% m(b% d_i) ! Mass ratio, as defined in Ritter 1988
                                       ! (Kolb & Ritter 1990 use the opposite!)
         ! consider range of validity for F1, do not extrapolate! Eq. A9 of Ritter 1988
         q_temp = min(max(q_old,0.5d0),10d0)
         F1_old = (1.23d0  + 0.5D0* log10_cr(q_temp))
        ! write(*,*) "approx F1 is" , F1_old
         g_q1 = 1/(L1*L1*L1) + q_old/((1-L1)*(1-L1)*(1-L1))
         F1_exact_1 = (b% separation/rl_d)*(b% separation/rl_d)*(b% separation/rl_d)&
                     /(sqrt(g_q1*(g_q1 - 1 - q_old)))
        ! write(*,*) "Kolb F1 is " , F1_exact_1
         F1_exact_1 = F1_exact_1/b% m(b% d_i)
         mdot_thick_alt = -2.0D0*pi*F1_exact_1*rl_d*rl_d*rl_d/(b% s_donor% cgrav(1))*mdot_thick_old
        ! write(*,*) "exact mdot thick is", mdot_thick_alt/Msun*secyer
         mdot_thick_old = -2.0D0*pi*F1_old*rl_d*rl_d*rl_d/(b% s_donor% cgrav(1)*b% m(b% d_i))*mdot_thick_old
         write(*,*) "original approx mdot thick is", mdot_thick_old/Msun*secyer
      
      end function calculate_kolb_mdot_thick
      
      subroutine get_info_for_kolb(b)
         type(binary_info), pointer :: b
         real(dp) :: F3, FF, G1, x_L1, q, g
         real(dp) :: mdot_thick0,  R_gas, dP, rl, s_div_rl
         integer :: i, indexR
         include 'formats.inc'

         !--------------------- Optically thick MT rate -----------------------------------------------
         ! As described in H. Ritter 1988, A&A 202,93-100 and U. Kolb and H. Ritter 1990, A&A 236,385-39

         ! First we need to find how deep inside the star the Roche lobe reaches. In other words the mesh point of the star at which R=R_RL
         b% mdot_thick = 0d0
         indexR=-1
         if(b% r(b% d_i)-b% rl(b% d_i) > 0.0d0) then
            i=1
            do while (b% s_donor% r(i) > b% rl(b% d_i))
               i=i+1
            end do
            
            if (i .eq. 1) then
               b% mdot_thick = 0d0
            else
               b% mdot_thick = calculate_kolb_mdot_thick(b, i-1, b% rl(b% d_i))
            end if
         end if

      end subroutine get_info_for_kolb

      subroutine get_info_for_ritter_outer(b)
         type(binary_info), pointer :: b
         real(dp) :: rho_exponent, F1, q, rho, p, grav, hp, v_th, rl3, q_temp, rl_d, g_q1, L, L2, L3
         include 'formats.inc'

         !--------------------- Optically thin MT rate -----------------------------------------------
         ! As described in H. Ritter 1988, A&A 202,93-100 and U. Kolb and H. Ritter 1990, A&A 236,385-392

         rho = b% s_donor% rho(1) ! density at surface in g/cm^3
         p = b% s_donor% p(1) ! pressure at surface in dynes/cm^2
         grav = b% s_donor% cgrav(1)*b% m(b% d_i)/(b% r(b% d_i))**2 ! local gravitational acceleration
         hp = p/(grav*rho) ! pressure scale height
         v_th = sqrt(kerg * b% s_donor% T(1) / (mp * b% s_donor% mu(1)))

       !This is the radius for the outer Lagrangian
         rl_d = outer_radius(b) !need radius fit

       !Added part
         q = b% m(b% a_i)/b% m(b% d_i)
         if (q < 1) then
            L3 = find_L3(b)
            L = L3
         else
            L2 = find_L2(b)
            L = L2
         end if

         g_q1 = 1/(L*L*L) + q/((1-L)*(1-L)*(1-L))
         F1 = (b% separation/rl_d)*(b% separation/rl_d)*(b% separation/rl_d)&
                     /(sqrt(g_q1*(g_q1 - 1 - q)))
         !F1 = F1/b% m(b% d_i)

       !change the radius --> done
         rl3 = (rl_d)*(rl_d)*(rl_d)
         b% mdot_thin0 = (2.0D0*pi/exp_cr(0.5d0)) * v_th*v_th*v_th * &
             rl3/(b% s_donor% cgrav(1)*b% m(b% d_i)) * rho * F1
         !Once again, do not extrapolate! Eq. (7) of Ritter 1988

       !needs to be changed
         q_temp = min(max(q,0.04d0),20d0)
         if (q_temp < 1.0d0) then
            b% ritter_h = hp/( 0.954D0 + 0.025D0*log10_cr(q_temp) - 0.038D0*(log10_cr(q_temp))**2 )
         else
            b% ritter_h = hp/( 0.954D0 + 0.039D0*log10_cr(q_temp) + 0.114D0*(log10_cr(q_temp))**2 )
         end if

         b% ritter_exponent = (b% r(b% d_i)-b% rl(b% d_i))/b% ritter_h

         if (b% mdot_scheme == "Kolb") then
            if (b% ritter_exponent > 0) then
               b% mdot_thin = -b% mdot_thin0
            else
               b% mdot_thin = -b% mdot_thin0 * exp_cr(b% ritter_exponent)
            end if
         else
            b% mdot_thin = -b% mdot_thin0 * exp_cr(b% ritter_exponent)
         end if

      end subroutine get_info_for_ritter_outer

      real(dp) function calculate_kolb_mdot_thick_outer(b, indexR, radius) result(mdot_thick)
         real(dp), intent(in) :: radius
         integer, intent(in) :: indexR
         real(dp) :: F1, F3, G1, dP, q, rl3, q_temp, L1, A1, B1, start, F1_exact_1,&
                     Phi, dA, unit, Phi_last, C1, D1, E1, linear, L, L3, L2 
         integer :: i
         type(binary_info), pointer :: b
         include 'formats.inc'

         !--------------------- Optically thick MT rate -----------------------------------------------
        
         q = b% m(b% a_i)/b% m(b% d_i)
         if (q < 1) then
            L3 = find_L3(b)
            L = L3
         else
            L2 = find_L2(b)
            L = L2
         end if

         A1 = ddPhidy(b,L, 0d0, 0d0)/2
         B1 = ddPhidz(b,L, 0d0, 0d0)/2
         C1 = ddddPhidy(b,L, 0d0, 0d0)/24
         D1 = ddddPhidz(b,L, 0d0, 0d0)/24
         E1 = ddddPhidyz(b,L, 0d0, 0d0)/4
  
         start = pi/(sqrt(A1*B1))
         linear = -2*pi*(3*D1*A1*A1 + A1*B1*E1 + 3*C1*B1*B1)/&
                       (8*sqrt(A1*A1*A1*A1*A1)*sqrt(B1*B1*B1*B1*B1))
  
         mdot_thick = 0d0
         Phi = 0d0
         do i=indexR-1, 1, -1
            Phi = Phi + (b% s_donor% P(i+1)-b% s_donor% P(i))/(b% s_donor% rho(i))
            G1 = b% s_donor% gamma1(i)
            F3 = sqrt(G1) * pow_cr(2d0/(G1+1d0), (G1+1d0)/(2d0*G1-2d0))
            unit = b% separation/(b% s_donor% cgrav(i)*(b% m(b% a_i) + b% m(b% d_i)))
            dA = (b% separation * b% separation * b% separation)/&
                 (b% s_donor% cgrav(i)*(b% m(b% a_i) + b% m(b% d_i)))*(start + linear*unit*Phi)
            mdot_thick = mdot_thick + dA*F3*sqrt(kerg * b% s_donor% T(i) / &
               (mp * b% s_donor% mu(i)))*(b% s_donor% P(i+1)-b% s_donor% P(i))
         end do

         i = indexR
         ! only take a fraction of dP for last cell
         G1 = b% s_donor% gamma1(i)
         F3 = sqrt(G1) * pow_cr(2d0/(G1+1d0), (G1+1d0)/(2d0*G1-2d0))
         dP = (b% s_donor% r(indexR) - radius) / &
            (b% s_donor% r(indexR) - b% s_donor% r(indexR+1)) * (b% s_donor% P(i+1)-b% s_donor% P(i))
         Phi_last = dP/(b% s_donor% rho(i))
         unit = b% separation/(b% s_donor% cgrav(i)*(b% m(b% a_i) + b% m(b% d_i)))
         dA = (b% separation * b% separation * b% separation)/&
               (b% s_donor% cgrav(i)*(b% m(b% a_i) + b% m(b% d_i)))*(start + linear*unit*Phi_last)
         mdot_thick = mdot_thick + dA*F3*sqrt(kerg * b% s_donor% T(i) / (mp*b% s_donor% mu(i)))*dP

         mdot_thick = - mdot_thick
         write(*,*) "outer mdot thick is", mdot_thick/Msun*secyer
 
      end function calculate_kolb_mdot_thick_outer
      
      real(dp) function get_info_for_kolb_outer(b) result(mdot_thick)
         type(binary_info), pointer :: b
         real(dp) :: F3, FF, G1, x_L1, q, g
         real(dp) :: R_gas, dP, rl, s_div_rl, radius_outer
         integer :: i, indexR
         include 'formats.inc'

         !--------------------- Optically thick MT rate --------------------------------------------
         ! As described in H. Ritter 1988, A&A 202,93-100 and U. Kolb and H. Ritter 1990, A&A 236,385-39

         ! First we need to find how deep inside the star the Roche lobe reaches. In other words the mesh point of the star at which R=R_RL

         q = b% m(b% d_i)/b% m(b% a_i)
         !radius_outer = 0.422d0 + 0.298d0*atan(1.06d0*log10_cr(q)+0.329d0)
         !radius_outer = radius_outer*b% separation

         radius_outer = outer_radius(b)

         mdot_thick = 0d0
         indexR=-1
         if(b% r(b% d_i)- radius_outer > 0.0d0) then
            i=1
            do while (b% s_donor% r(i) > radius_outer)
               i=i+1
            end do
            
            if (i .eq. 1) then
               mdot_thick = 0d0
            else
               mdot_thick = calculate_kolb_mdot_thick_outer(b, i-1, radius_outer)
            end if
         end if

      end function get_info_for_kolb_outer

      subroutine get_info_for_arras(b)
         type(binary_info), pointer :: b
         real(dp) :: q, rho, p, grav, hp, v_th
         real(dp) :: area,Asl,G,ma,md,mfac1,mfac2,my_mdot_thin,my_ritter_exponent,Omega,&
            phi,phiL1,q13,rfac,rhoL1,rv,rvL1,sep
         include 'formats.inc'

         !--------------------- Optically thin MT rate -----------------------------------------------
         ! Ritter 1988 but with better fits for the various formulas that work at extreme q

         rho = b% s_donor% rho(1) ! density at surface in g/cm^3
         p = b% s_donor% p(1) ! pressure at surface in dynes/cm^2
         grav = b% s_donor% cgrav(1)*b% m(b% d_i)/(b% r(b% d_i))**2 ! local gravitational acceleration
         hp = p/(grav*rho) ! pressure scale height
         v_th = sqrt(kerg * b% s_donor% T(1) / (mp * b% s_donor% mu(1)))

         q = b% m(b% a_i) / b% m(b% d_i)
         G = b% s_donor% cgrav(1)
         md = b% m(b% d_i)
         ma = b% m(b% a_i)
         sep = b% separation
         Omega = 2.d0*pi / b% period
         rvL1 = b% rl(b% d_i)
         rv = b% r(b% d_i)
         mfac1 = (md+ma)/md
         mfac2 = ( (md+ma)**2 + 3d0*ma*(md+ma) + 9d0*ma**2 ) / md**2
         rfac=rvL1/sep
         phiL1 = - G*md/rvL1 &
              * ( 1.d0 +  mfac1*pow3(rfac)/3.d0 + 4.d0*mfac2*pow6(rfac)/45.d0  )
              rfac=rv/sep
         phi = - G*md/rv &
              * ( 1.d0 +  mfac1*pow3(rfac)/3.d0 + 4.d0*mfac2*pow6(rfac)/45.d0  )
         my_ritter_exponent = - (phiL1-phi)/v_th**2
         rhoL1 = rho/sqrt(exp_cr(1.d0)) * exp_cr( my_ritter_exponent )
              q13=pow_cr(q,one_third)
              Asl = 4.d0 + 4.16d0/(-0.96d0 + q13 + 1.d0/q13)
         area = 2.d0 * pi * (v_th/Omega)**2 / sqrt( Asl*(Asl-1d0) )
         my_mdot_thin = - rhoL1 * v_th * area
         b% mdot_thin = my_mdot_thin

      end subroutine get_info_for_arras

      subroutine get_info_for_ritter_eccentric(b)
         type(binary_info), pointer :: b
         integer :: i
         real(dp) :: rho_exponent, F1, q, q_temp, rho, p, grav, hp, v_th, dm
         real(dp), DIMENSION(b% anomaly_steps):: mdot0, mdot, Erit, rl_d
         include 'formats.inc'
         
         ! Optically thin MT rate adapted for eccentric orbits 
         ! As described in H. Ritter 1988, A&A 202,93-100 and U. Kolb and H. Ritter 1990, A&A 236,385-392

         rho = b% s_donor% rho(1) ! density at surface in g/cm^3
         p = b% s_donor% p(1) ! pressure at surface in dynes/cm^2
         grav = b% s_donor% cgrav(1)*b% m(b% d_i)/(b% r(b% d_i))**2 ! local gravitational acceleration
         hp = p/(grav*rho) ! pressure scale height
         v_th = sqrt(kerg * b% s_donor% T(1) / (mp * b% s_donor% mu(1))) ! kerg = Boltzmann's constant
         
         ! phase dependant RL radius
         rl_d = b% rl(b% d_i) * (1d0 - b% eccentricity**2) / &
                (1 + b% eccentricity * cos(b% theta_co) )

         q = b% m(b% a_i)/b% m(b% d_i) ! Mass ratio, as defined in Ritter 1988
                                       ! (Kolb & Ritter 1990 use the opposite!)
         q_temp = min(max(q,0.5d0),10d0)
         F1 = (1.23d0  + 0.5D0* log10_cr(q_temp))

         mdot0 = (2.0D0*pi/exp_cr(0.5d0)) * pow3(v_th) * rl_d*rl_d*rl_d / &
             (b% s_donor% cgrav(1)*b% m(b% d_i)) * rho * F1   
             
         q_temp = min(max(q,0.04d0),20d0)
         if (q_temp < 1.0d0) then
            b% ritter_h = hp/( 0.954D0 + 0.025D0*log10_cr(q_temp) - 0.038D0*(log10_cr(q_temp))**2 )
         else
            b% ritter_h = hp/( 0.954D0 + 0.039D0*log10_cr(q_temp) + 0.114D0*(log10_cr(q_temp))**2 )
         end if

         Erit = (b% r(b% d_i)- rl_d) / b% ritter_h

         if (b% mdot_scheme == "Kolb") then
            do i = 1,b% anomaly_steps
               if (Erit(i) > 0) then
                  mdot(i) = -1 * mdot0(i)
               else
                  mdot(i) = -1 * mdot0(i) * exp(Erit(i))
               end if
            end do
         else
            mdot = -1 * mdot0 * exp(Erit)
         end if
         
         b% mdot_donor_theta = mdot
         
         !integrate to get total massloss
         dm = 0d0
         do i = 2,b% anomaly_steps ! trapezoidal integration
            dm = dm + 0.5d0 * (mdot(i-1) + mdot(i)) * (b% time_co(i) - b% time_co(i-1)) 
         end do
         
         b% mdot_thin = dm

      end subroutine get_info_for_ritter_eccentric
      
      subroutine get_info_for_kolb_eccentric(b)
         type(binary_info), pointer :: b
         real(dp) :: e, dm
         integer :: i, j
         real(dp), DIMENSION(b% anomaly_steps):: rl_d_i, mdot_thick_i
         include 'formats.inc'
         
         ! Optically thick MT rate adapted for eccentric orbits
         ! As described in H. Ritter 1988, A&A 202,93-100 and U. Kolb and H. Ritter 1990, A&A 236,385-392

         b% mdot_thick = 0d0
         e = b% eccentricity
         
         ! If the radius of the donor is smaller as the smallest RL radius,
         ! there is only atmospheric RLOF, thus return.
         if ( b% r(b% d_i) < b% rl(b% d_i) * (1-e**2)/(1+e) ) then
            return
         end if
         
         ! phase dependant RL radius
         rl_d_i = b% rl(b% d_i) * (1d0 - b% eccentricity**2) / &
                  (1 + b% eccentricity * cos(b% theta_co) )
         
         ! For each point in the orbit calculate mdot_thick 
         do i = 1,b% anomaly_steps
            ! find how deep in the star we are
            j=1
            do while (b% s_donor% r(j) > rl_d_i(i))
               j=j+1
            end do
            
            ! calculate mdot_thick
            if (j .eq. 1) then
               mdot_thick_i(i) = 0d0
            else
               mdot_thick_i(i) = calculate_kolb_mdot_thick(b, j-1, rl_d_i(i))
            end if
         end do
         
         b% mdot_donor_theta = b% mdot_donor_theta + mdot_thick_i
         
         ! Integrate mdot_thick over the orbit
         dm = 0d0
         do i = 2,b% anomaly_steps ! trapezoidal integration
            dm = dm + 0.5d0 * (mdot_thick_i(i-1) + mdot_thick_i(i)) * &
                              (b% time_co(i) - b% time_co(i-1)) 
         end do
         
         b% mdot_thick = dm
         
      end subroutine get_info_for_kolb_eccentric
      
      end module run_binary_extras

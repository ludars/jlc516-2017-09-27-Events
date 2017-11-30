WITH
    timeframe AS (
      SELECT
          MIN(start_date) AS start_date
        , MAX(end_date)   AS end_date
      FROM (
        SELECT
            atvfisc_code                   AS fy
          , TRUNC(atvfisc_start_date)      AS start_date
          , TRUNC(atvfisc_end_date)        AS end_date
          , ROW_NUMBER()
            OVER (
              PARTITION BY NULL
              ORDER BY atvfisc_code DESC ) AS r
        FROM
          atvfisc
        WHERE
          TRUNC(SYSDATE) > TRUNC(atvfisc_start_date))
      WHERE r <= 4)
  , event_comments AS (
    SELECT
      *
    FROM (
      SELECT
          gersubj_evnt_crn         AS evnt_crn
        , gersubj_func_code        AS func_code
        , gersubj_subj_code        AS subj_code
        , TO_CHAR(gerfcom_comment) AS text
      FROM
        gerfcom
        , gersubj
        , gtvsubj
      WHERE
        gersubj_subj_code IN (
          'EVTITL'
          , 'EDESC'
          , 'BUDGET'
          , 'PRIM'
          , 'ECOMM'
          , 'LEAD'
          , 'LEADPR')
        AND gerfcom_evnt_crn = gersubj_evnt_crn
        AND gerfcom_func_code = gersubj_func_code
        AND gerfcom_grp_seq_no = gersubj_grp_seq_no
        AND gersubj_subj_code = gtvsubj_code)
      PIVOT (LISTAGG(text, ';')
      WITHIN GROUP (
        ORDER BY NULL)
        FOR subj_code
        IN (
          'EVTITL' AS event_title
          , 'EDESC' AS event_desc
          , 'BUDGET' AS event_budget
          , 'PRIM' AS primary_region
          , 'ECOMM' AS event_comments
          , 'LEAD' AS univ_leadership_request
          , 'LEADPR' AS presidental_request)))
  , target_audience AS (
    SELECT
      *
    FROM (
      SELECT
          gerfcon_evnt_crn  AS evnt_crn
        , gerfcon_func_code AS func_code
        , gerfcon_targ_code AS targ_code
      FROM
        gerfcon
        , gtvsdax
      WHERE gtvsdax_internal_code = 'TARG_AUD'
            AND gerfcon_targ_code = gtvsdax_external_code)
      PIVOT (MAX('Y')
        FOR targ_code
        IN (
          'AFACC' AS affinity_accounting,
          'AFCHA' AS affinity_choral_arts,
          'AFFIN' AS affinity_financial_services,
          'AFHEA' AS affinity_healthcare,
          'AFLUG' AS affinity_lugala,
          'AFLUI' AS affinity_lu_innovators,
          'AFM97' AS affinity_marching_97,
          'AFRE' AS affinity_real_estate,
          'WRES' AS affinity_wrestling_club,
          'WSC' AS affinity_wall_street_council,
          'BDCAS' AS board_cas_advisory_council,
          'BDCBE' AS board_cbe_advisory_council,
          'BDCOE' AS board_coe_advisory_council,
          'BDGAC' AS board_greek_advisory_council,
          'BDPAR' AS board_parents_council_members,
          'BDRCE' AS board_rceas_advisory_council,
          'LBRD' AS board_luaa_board_of_directors,
          'LLC' AS board_leadership_councilllc,
          'TRUST' AS board_board_of_trustees,
          'YAC' AS board_young_alumni_council,
          'GRACO' AS greek_alpha_chi_omega,
          'GRACR' AS greek_alpha_chi_rho,
          'GRAEP' AS greek_alpha_epsilon_pi_fmr,
          'GRAGD' AS greek_alpha_gamma_delta,
          'GRALL' AS greek_all_members,
          'GRAOP' AS greek_alpha_omicron_pi,
          'GRAP' AS greek_alpha_phi,
          'GRAPA' AS greek_alpha_phi_alpha,
          'GRAPO' AS greek_alpha_phi_omega,
          'GRASP' AS greek_alpha_sigma_phi,
          'GRATO' AS greek_alpha_tau_omega,
          'GRBTP' AS greek_beta_theta_pi,
          'GRCO' AS greek_chi_omega,
          'GRCPH' AS greek_chi_phi,
          'GRCPS' AS greek_chi_psi,
          'GRDCH' AS greek_delta_chi,
          'GRDG' AS greek_delta_gamma,
          'GRDPH' AS greek_delta_phi,
          'GRDSP' AS greek_delta_sigma_phi,
          'GRDTD' AS greek_delta_tau_delta,
          'GRDUP' AS greek_delta_upsilon,
          'GRDZ' AS greek_delta_zeta,
          'GRGPB' AS greek_gamma_phi_beta,
          'GRKA' AS greek_kappa_alpha,
          'GRKAP' AS greek_kappa_alpha_psi,
          'GRKAT' AS greek_kappa_alpha_theta,
          'GRKD' AS greek_kappa_delta,
          'GRKS' AS greek_kappa_sigma,
          'GRLCA' AS greek_lambda_chi_alpha,
          'GRLSU' AS greek_lambda_sigma_upsilon,
          'GRLTA' AS greek_lambda_theta_alpha,
          'GRMSU' AS greek_mu_sigma_upsilon,
          'GROPP' AS greek_omega_psi_phi,
          'GRPBP' AS greek_pi_beta_phi,
          'GRPDT' AS greek_phi_delta_theta,
          'GRPGD' AS greek_phi_gamma_delta,
          'GRPHT' AS greek_phi_kappa_theta,
          'GRPKA' AS greek_pi_kappa_alpha,
          'GRPKP' AS greek_phi_kappa_psi_former,
          'GRPLP' AS greek_pi_lambda_phi,
          'GRPSK' AS greek_phi_sigma_kappa,
          'GRPSU' AS greek_psi_upsilon,
          'GRSAM' AS greek_sigma_alpha_mu,
          'GRSC' AS greek_sigma_chi,
          'GRSN' AS greek_sigma_nu,
          'GRSPE' AS greek_sigma_phi_epsilon,
          'GRSPH' AS greek_sigma_phi,
          'GRTCH' AS greek_theta_chi,
          'GRTDC' AS greek_theta_delta_chi,
          'GRTDP' AS greek_tau_delta_phi,
          'GRTEP' AS greek_tau_epsilon_phi,
          'GRTKP' AS greek_theta_kappa_phi,
          'GRTXI' AS greek_theta_xi,
          'GRZP' AS greek_zeta_psi,
          'GRZTA' AS greek_zeta_tau_alpha,
          '$1M +' AS general_one_million_donors,
          'ALMB' AS general_undergraduate_alumni,
          'APS' AS general_asa_packer_society,
          'GALMG' AS general_graduate_alumni,
          'GALUM' AS general_all_alumni,
          'GFAST' AS general_lu_faculty_and_staff,
          'GGLOB' AS general_global_village,
          'GLEHC' AS general_greater_lehigh_valley,
          'GSENR' AS general_current_seniors,
          'GTOWR' AS general_tower_society,
          'GVOLT' AS general_volunteers,
          'GYALM' AS general_young_alumni,
          'LGP' AS general_lgo_prospects,
          'PARN' AS general_parents,
          'PR06' AS general_prospects_rated_6_,
          'PR12' AS general_prospects_rated_12,
          'STDT' AS general_students,
          'ALSK' AS ak_alaska_club,
          'ARIZ' AS az_arizona_club,
          'ARKS' AS ar_arkansas_club,
          'ATLA' AS ga_atlanta_club,
          'AUST' AS tx_austinsan_antonio_tx_club,
          'BAL' AS affinity_balance,
          'BALT' AS md_maryland_club,
          'BAMA' AS al_alabama_club,
          'BOST' AS ma_boston_club,
          'CCMA' AS ma_cape_cod_club,
          'CEFL' AS fl_central_florida_club,
          'CENJ' AS nj_central_new_jersey,
          'CENY' AS ny_central_new_york_club,
          'CHAR' AS nc_charlotte_nc_club,
          'CHIC' AS il_chicago_club,
          'COLA' AS sc_greater_columbia_club,
          'COSC' AS sc_coastal_south_carolina_club,
          'CTVA' AS ct_connecticut_valley_club,
          'DAFW' AS tx_dallasft_worth_tx_club,
          'DELA' AS de_delaware_club,
          'FFCT' AS ct_fairfield_county_club,
          'GOCO' AS njny_gold_coast_club,
          'GUAM' AS gu_guam_club,
          'HAWI' AS hi_hawaii_club,
          'HOME' AS pa_home_club,
          'IDAH' AS id_idaho_club,
          'INDY' AS in_indiana_club,
          'INTL' AS international,
          'IOWA' AS ia_iowa_club,
          'KASA' AS ks_kansas_club,
          'KENT' AS ky_kentucky_club,
          'LAWY' AS affinity_lehigh_lawyers,
          'LGIS' AS ny_long_island_club,
          'LUMEC' AS affinity_lumeca,
          'MAIN' AS me_maine_club,
          'MICH' AS mi_michigan_club,
          'MIFL' AS fl_midflorida_club,
          'MILT' AS aaaeap_military_apo_fpo,
          'MISS' AS ms_mississippi_club,
          'MONT' AS mt_montana_club,
          'NDAK' AS nd_north_dakota_club,
          'NEBR' AS ne_nebraska_club,
          'NEPA' AS pa_northeast_pennsylvania_club,
          'NEVA' AS nv_nevada_club,
          'NEWH' AS nh_new_hampshire_club,
          'NJSH' AS nj_monmouth__ocean_co_club,
          'NMEX' AS nm_new_mexico_club,
          'NOFL' AS fl_northern_florida_club,
          'NOLA' AS la_new_orleans_club,
          'NONJ' AS nj_northern_new_jersey_club,
          'NONY' AS ny_northern_new_york_club,
          'NOOH' AS oh_northern_ohio_club,
          'NWIN' AS in_northwest_indiana_club,
          'NWPA' AS pa_northwest_pennsylvania_club,
          'NYNY' AS ny_new_york_club,
          'OHVA' AS oh_ohio_valley_club,
          'OKLA' AS ok_oklahoma_club,
          'ORCY' AS ca_orange_county_club,
          'PACN' AS orwa_pacific_northwest_club,
          'PAFL' AS fl_florida_panhandle,
          'PBGH' AS pa_pittsburgh_club,
          'PHIL' AS pa_philadelphia_club,
          'PRCO' AS pr_puerto_rico_club,
          'RDNC' AS nc_raleighdurham_club,
          'RHIS' AS ri_rhode_island_club,
          'ROCH' AS ny_rochester_club,
          'ROCK' AS co_rocky_mountain_club,
          'SAND' AS ca_san_diego_club,
          'SARA' AS fl_sarasota_club,
          'SDAK' AS sd_south_dakota_club,
          'SEPA' AS pa_southestern_pa_club,
          'SFBY' AS ca_san_francisco_bay_area_club,
          'SFLA' AS fl_southern_florida_club,
          'SOAZ' AS az_southern_arizona_club,
          'SOCA' AS ca_southern_california_club,
          'SONJ' AS nj_southern_new_jersey_club,
          'SONY' AS ny_southern_new_york_club,
          'STLU' AS mo_st_louis_club,
          'SUSQ' AS pa_lancastercentralyork_club,
          'SVHH' AS scga_savannahhilton_head_cb,
          'SWFL' AS fl_southwest_florida,
          'TBFL' AS fl_tampa_bay_club,
          'TCFL' AS fl_treasure_coast,
          'TENN' AS tn_tennessee_club,
          'TEXS' AS tx_houston_tx_club,
          'TWIN' AS mn_twin_cities_club,
          'UTAH' AS ut_utah_club,
          'VERT' AS vt_vermont_club,
          'VIRG' AS va_virginia_club,
          'VRIS' AS vi_virgin_islands_club,
          'WADC' AS dc_washington_d_c_club,
          'WCAR' AS scnc_western_carolinas_club,
          'WECH' AS ny_westchesterrockland_club,
          'WENY' AS ny_west_new_york_club,
          'WISC' AS wi_wisconsin_club,
          'WVIR' AS wv_west_virginia_club,
          'WYOM' AS wy_wyoming_club,
          'SPALL' AS sports_all_members,
          'SPBAB' AS sports_baseball,
          'SPBAK' AS sports_basketball,
          'SPCRC' AS sports_cross_country,
          'SPFLD' AS sports_field_hockey,
          'SPFTB' AS sports_football,
          'SPGLF' AS sports_golf,
          'SPLAC' AS sports_lacrosse,
          'SPROW' AS sports_rowing,
          'SPSOC' AS sports_soccer,
          'SPSOF' AS sports_softball,
          'SPSWD' AS sports_swimming__diving,
          'SPTEN' AS sports_tennis,
          'SPTRF' AS sports_track__field,
          'SPVOL' AS sports_volleyball,
          'SPWRS' AS sports_wrestling,
          'DONR' AS donor,
          'GSEDR' AS general_senior_alumni)))
  , virtual AS (
    SELECT
        gerfcon_evnt_crn  AS evnt_crn
      , gerfcon_func_code AS func_code
    FROM gerfcon
    WHERE gerfcon_targ_code = 'VIR')
SELECT
    gerattd_pidm                                                  AS attendee_pidm
  , slbevnt_crn                                                   AS evnt_crn
  , gebfunc_func_code                                             AS func_code
  , atvfisc_code                                                  AS fy
  , TO_CHAR(ssrmeet_start_date, 'MM/DD/YYYY')                     AS start_date
  , TO_CHAR(TO_DATE(ssrmeet_begin_time, 'HH24:MI'), 'HH12:MI AM') AS start_time
  , TO_CHAR(ssrmeet_end_date, 'MM/DD/YYYY')                       AS end_date
  , TO_CHAR(TO_DATE(ssrmeet_end_time, 'HH24:MI'), 'HH12:MI AM')   AS end_time
  , slbevnt_desc                                                  AS event_short_title
  , event_comments.event_title                                    AS event_long_title
  , event_comments.event_desc
  , gtvfunc_desc                                                  AS func_short_desc
  , gebfunc_long_description                                      AS func_long_desc
  , gebfunc_location                                              AS func_location
  , gtvrsvp_desc                                                  AS attendee_rsvp
  , gerattd_rsvp_date                                             AS attendee_rsvp_date
  , NVL(gerattd_attendance_ind, 'N')                              AS attended_ind
  , gtvfees_desc                                                  AS attendee_fees
  , gerattd_fees_desc                                             AS attendee_fees_desc
  , gerattd_fee_date                                              AS attendee_fees_date
  , gerattd_name_tag_name                                         AS attendee_name_tag_name
  , gerattd_place_card_name                                       AS attendee_place_card_name
  , gerattd_ticket_cnt                                            AS attendee_ticket_count
  , gtvmenu_desc                                                  AS attendee_menu
  , CASE
    WHEN gerattd_involve_ind = 'I' THEN 'Invited'
    WHEN gerattd_involve_ind = 'G' THEN 'Guest'
    ELSE 'Unknown'
    END                                                           AS attendee_involvement_type
  , gerattd_comment                                               AS attendee_comment
  , gerattd_number_of_guests                                      AS attendee_number_of_guests
  , gtvfsta_desc                                                  AS func_status
  , stvcoll_desc                                                  AS event_college
  , NVL(stvdept_desc_long, stvdept_desc)                          AS event_department
  , gtvsysi_desc                                                  AS event_owner
  , event_comments.event_budget
  , gebfunc_plan_attend                                           AS func_num_attend_max
  , gebfunc_brk_even_attend                                       AS func_num_attend_min
  , gebfunc_committee_name                                        AS func_num_dar_staff_requested
  , event_comments.univ_leadership_request
  , event_comments.presidental_request
  , evnt_etyp.stvetyp_desc                                        AS event_type
  , func_etyp.stvetyp_desc                                        AS func_type
  , gtvpurp_desc                                                  AS func_primary_purpose
  , gtvemph_desc                                                  AS func_secondary_purpose
  , event_comments.primary_region
  , NVL(target_audience.affinity_accounting, 'N')                 AS affinity_accounting
  , NVL(target_audience.affinity_balance, 'N')                    AS affinity_balance
  , NVL(target_audience.affinity_choral_arts, 'N')                AS affinity_choral_arts
  , NVL(target_audience.affinity_financial_services, 'N')         AS affinity_financial_services
  , NVL(target_audience.affinity_healthcare, 'N')                 AS affinity_healthcare
  , NVL(target_audience.affinity_lehigh_lawyers, 'N')             AS affinity_lehigh_lawyers
  , NVL(target_audience.affinity_lu_innovators, 'N')              AS affinity_lu_innovators
  , NVL(target_audience.affinity_lugala, 'N')                     AS affinity_lugala
  , NVL(target_audience.affinity_lumeca, 'N')                     AS affinity_lumeca
  , NVL(target_audience.affinity_marching_97, 'N')                AS affinity_marching_97
  , NVL(target_audience.affinity_real_estate, 'N')                AS affinity_real_estate
  , NVL(target_audience.affinity_wall_street_council, 'N')        AS affinity_wall_street_council
  , NVL(target_audience.affinity_wrestling_club, 'N')             AS affinity_wrestling_club
  , NVL(target_audience.board_board_of_trustees, 'N')             AS board_board_of_trustees
  , NVL(target_audience.board_cas_advisory_council, 'N')          AS board_cas_advisory_council
  , NVL(target_audience.board_cbe_advisory_council, 'N')          AS board_cbe_advisory_council
  , NVL(target_audience.board_coe_advisory_council, 'N')          AS board_coe_advisory_council
  , NVL(target_audience.board_greek_advisory_council, 'N')        AS board_greek_advisory_council
  , NVL(target_audience.board_leadership_councilllc, 'N')         AS board_leadership_councilllc
  , NVL(target_audience.board_luaa_board_of_directors, 'N')       AS board_luaa_board_of_directors
  , NVL(target_audience.board_parents_council_members, 'N')       AS board_parents_council_members
  , NVL(target_audience.board_rceas_advisory_council, 'N')        AS board_rceas_advisory_council
  , NVL(target_audience.board_young_alumni_council, 'N')          AS board_young_alumni_council
  , NVL(target_audience.donor, 'N')                               AS donor
  , NVL(target_audience.greek_all_members, 'N')                   AS greek_all_members
  , NVL(target_audience.greek_alpha_chi_omega, 'N')               AS greek_alpha_chi_omega
  , NVL(target_audience.greek_alpha_chi_rho, 'N')                 AS greek_alpha_chi_rho
  , NVL(target_audience.greek_alpha_epsilon_pi_fmr, 'N')          AS greek_alpha_epsilon_pi_fmr
  , NVL(target_audience.greek_alpha_gamma_delta, 'N')             AS greek_alpha_gamma_delta
  , NVL(target_audience.greek_alpha_omicron_pi, 'N')              AS greek_alpha_omicron_pi
  , NVL(target_audience.greek_alpha_phi, 'N')                     AS greek_alpha_phi
  , NVL(target_audience.greek_alpha_phi_alpha, 'N')               AS greek_alpha_phi_alpha
  , NVL(target_audience.greek_alpha_phi_omega, 'N')               AS greek_alpha_phi_omega
  , NVL(target_audience.greek_alpha_sigma_phi, 'N')               AS greek_alpha_sigma_phi
  , NVL(target_audience.greek_alpha_tau_omega, 'N')               AS greek_alpha_tau_omega
  , NVL(target_audience.greek_beta_theta_pi, 'N')                 AS greek_beta_theta_pi
  , NVL(target_audience.greek_chi_omega, 'N')                     AS greek_chi_omega
  , NVL(target_audience.greek_chi_phi, 'N')                       AS greek_chi_phi
  , NVL(target_audience.greek_chi_psi, 'N')                       AS greek_chi_psi
  , NVL(target_audience.greek_delta_chi, 'N')                     AS greek_delta_chi
  , NVL(target_audience.greek_delta_gamma, 'N')                   AS greek_delta_gamma
  , NVL(target_audience.greek_delta_phi, 'N')                     AS greek_delta_phi
  , NVL(target_audience.greek_delta_sigma_phi, 'N')               AS greek_delta_sigma_phi
  , NVL(target_audience.greek_delta_tau_delta, 'N')               AS greek_delta_tau_delta
  , NVL(target_audience.greek_delta_upsilon, 'N')                 AS greek_delta_upsilon
  , NVL(target_audience.greek_delta_zeta, 'N')                    AS greek_delta_zeta
  , NVL(target_audience.greek_gamma_phi_beta, 'N')                AS greek_gamma_phi_beta
  , NVL(target_audience.greek_kappa_alpha, 'N')                   AS greek_kappa_alpha
  , NVL(target_audience.greek_kappa_alpha_psi, 'N')               AS greek_kappa_alpha_psi
  , NVL(target_audience.greek_kappa_alpha_theta, 'N')             AS greek_kappa_alpha_theta
  , NVL(target_audience.greek_kappa_delta, 'N')                   AS greek_kappa_delta
  , NVL(target_audience.greek_kappa_sigma, 'N')                   AS greek_kappa_sigma
  , NVL(target_audience.greek_lambda_chi_alpha, 'N')              AS greek_lambda_chi_alpha
  , NVL(target_audience.greek_lambda_sigma_upsilon, 'N')          AS greek_lambda_sigma_upsilon
  , NVL(target_audience.greek_lambda_theta_alpha, 'N')            AS greek_lambda_theta_alpha
  , NVL(target_audience.greek_mu_sigma_upsilon, 'N')              AS greek_mu_sigma_upsilon
  , NVL(target_audience.greek_omega_psi_phi, 'N')                 AS greek_omega_psi_phi
  , NVL(target_audience.greek_phi_delta_theta, 'N')               AS greek_phi_delta_theta
  , NVL(target_audience.greek_phi_gamma_delta, 'N')               AS greek_phi_gamma_delta
  , NVL(target_audience.greek_phi_kappa_psi_former, 'N')          AS greek_phi_kappa_psi_former
  , NVL(target_audience.greek_phi_kappa_theta, 'N')               AS greek_phi_kappa_theta
  , NVL(target_audience.greek_phi_sigma_kappa, 'N')               AS greek_phi_sigma_kappa
  , NVL(target_audience.greek_pi_beta_phi, 'N')                   AS greek_pi_beta_phi
  , NVL(target_audience.greek_pi_kappa_alpha, 'N')                AS greek_pi_kappa_alpha
  , NVL(target_audience.greek_pi_lambda_phi, 'N')                 AS greek_pi_lambda_phi
  , NVL(target_audience.greek_psi_upsilon, 'N')                   AS greek_psi_upsilon
  , NVL(target_audience.greek_sigma_alpha_mu, 'N')                AS greek_sigma_alpha_mu
  , NVL(target_audience.greek_sigma_chi, 'N')                     AS greek_sigma_chi
  , NVL(target_audience.greek_sigma_nu, 'N')                      AS greek_sigma_nu
  , NVL(target_audience.greek_sigma_phi, 'N')                     AS greek_sigma_phi
  , NVL(target_audience.greek_sigma_phi_epsilon, 'N')             AS greek_sigma_phi_epsilon
  , NVL(target_audience.greek_tau_delta_phi, 'N')                 AS greek_tau_delta_phi
  , NVL(target_audience.greek_tau_epsilon_phi, 'N')               AS greek_tau_epsilon_phi
  , NVL(target_audience.greek_theta_chi, 'N')                     AS greek_theta_chi
  , NVL(target_audience.greek_theta_delta_chi, 'N')               AS greek_theta_delta_chi
  , NVL(target_audience.greek_theta_kappa_phi, 'N')               AS greek_theta_kappa_phi
  , NVL(target_audience.greek_theta_xi, 'N')                      AS greek_theta_xi
  , NVL(target_audience.greek_zeta_psi, 'N')                      AS greek_zeta_psi
  , NVL(target_audience.greek_zeta_tau_alpha, 'N')                AS greek_zeta_tau_alpha
  , NVL(target_audience.general_all_alumni, 'N')                  AS general_all_alumni
  , NVL(target_audience.general_asa_packer_society, 'N')          AS general_asa_packer_society
  , NVL(target_audience.general_current_seniors, 'N')             AS general_current_seniors
  , NVL(target_audience.general_global_village, 'N')              AS general_global_village
  , NVL(target_audience.general_graduate_alumni, 'N')             AS general_graduate_alumni
  , NVL(target_audience.general_greater_lehigh_valley, 'N')       AS general_greater_lehigh_valley
  , NVL(target_audience.general_lgo_prospects, 'N')               AS general_lgo_prospects
  , NVL(target_audience.general_lu_faculty_and_staff, 'N')        AS general_lu_faculty_and_staff
  , NVL(target_audience.general_one_million_donors, 'N')          AS general_one_million_donors
  , NVL(target_audience.general_parents, 'N')                     AS general_parents
  , NVL(target_audience.general_prospects_rated_12, 'N')          AS general_prospects_rated_12
  , NVL(target_audience.general_prospects_rated_6_, 'N')          AS general_prospects_rated_6_
  , NVL(target_audience.general_senior_alumni, 'N')               AS general_senior_alumni
  , NVL(target_audience.general_students, 'N')                    AS general_students
  , NVL(target_audience.general_tower_society, 'N')               AS general_tower_society
  , NVL(target_audience.general_undergraduate_alumni, 'N')        AS general_undergraduate_alumni
  , NVL(target_audience.general_volunteers, 'N')                  AS general_volunteers
  , NVL(target_audience.general_young_alumni, 'N')                AS general_young_alumni
  , NVL(target_audience.aaaeap_military_apo_fpo, 'N')             AS aaaeap_military_apo_fpo
  , NVL(target_audience.ak_alaska_club, 'N')                      AS ak_alaska_club
  , NVL(target_audience.al_alabama_club, 'N')                     AS al_alabama_club
  , NVL(target_audience.ar_arkansas_club, 'N')                    AS ar_arkansas_club
  , NVL(target_audience.az_arizona_club, 'N')                     AS az_arizona_club
  , NVL(target_audience.az_southern_arizona_club, 'N')            AS az_southern_arizona_club
  , NVL(target_audience.ca_orange_county_club, 'N')               AS ca_orange_county_club
  , NVL(target_audience.ca_san_diego_club, 'N')                   AS ca_san_diego_club
  , NVL(target_audience.ca_san_francisco_bay_area_club, 'N')      AS ca_san_francisco_bay_area_club
  , NVL(target_audience.ca_southern_california_club, 'N')         AS ca_southern_california_club
  , NVL(target_audience.co_rocky_mountain_club, 'N')              AS co_rocky_mountain_club
  , NVL(target_audience.ct_connecticut_valley_club, 'N')          AS ct_connecticut_valley_club
  , NVL(target_audience.ct_fairfield_county_club, 'N')            AS ct_fairfield_county_club
  , NVL(target_audience.dc_washington_d_c_club, 'N')              AS dc_washington_d_c_club
  , NVL(target_audience.de_delaware_club, 'N')                    AS de_delaware_club
  , NVL(target_audience.fl_central_florida_club, 'N')             AS fl_central_florida_club
  , NVL(target_audience.fl_florida_panhandle, 'N')                AS fl_florida_panhandle
  , NVL(target_audience.fl_midflorida_club, 'N')                  AS fl_midflorida_club
  , NVL(target_audience.fl_northern_florida_club, 'N')            AS fl_northern_florida_club
  , NVL(target_audience.fl_sarasota_club, 'N')                    AS fl_sarasota_club
  , NVL(target_audience.fl_southern_florida_club, 'N')            AS fl_southern_florida_club
  , NVL(target_audience.fl_southwest_florida, 'N')                AS fl_southwest_florida
  , NVL(target_audience.fl_tampa_bay_club, 'N')                   AS fl_tampa_bay_club
  , NVL(target_audience.fl_treasure_coast, 'N')                   AS fl_treasure_coast
  , NVL(target_audience.ga_atlanta_club, 'N')                     AS ga_atlanta_club
  , NVL(target_audience.gu_guam_club, 'N')                        AS gu_guam_club
  , NVL(target_audience.hi_hawaii_club, 'N')                      AS hi_hawaii_club
  , NVL(target_audience.ia_iowa_club, 'N')                        AS ia_iowa_club
  , NVL(target_audience.id_idaho_club, 'N')                       AS id_idaho_club
  , NVL(target_audience.il_chicago_club, 'N')                     AS il_chicago_club
  , NVL(target_audience.in_indiana_club, 'N')                     AS in_indiana_club
  , NVL(target_audience.in_northwest_indiana_club, 'N')           AS in_northwest_indiana_club
  , NVL(target_audience.international, 'N')                       AS international
  , NVL(target_audience.ks_kansas_club, 'N')                      AS ks_kansas_club
  , NVL(target_audience.ky_kentucky_club, 'N')                    AS ky_kentucky_club
  , NVL(target_audience.la_new_orleans_club, 'N')                 AS la_new_orleans_club
  , NVL(target_audience.ma_boston_club, 'N')                      AS ma_boston_club
  , NVL(target_audience.ma_cape_cod_club, 'N')                    AS ma_cape_cod_club
  , NVL(target_audience.md_maryland_club, 'N')                    AS md_maryland_club
  , NVL(target_audience.me_maine_club, 'N')                       AS me_maine_club
  , NVL(target_audience.mi_michigan_club, 'N')                    AS mi_michigan_club
  , NVL(target_audience.mn_twin_cities_club, 'N')                 AS mn_twin_cities_club
  , NVL(target_audience.mo_st_louis_club, 'N')                    AS mo_st_louis_club
  , NVL(target_audience.ms_mississippi_club, 'N')                 AS ms_mississippi_club
  , NVL(target_audience.mt_montana_club, 'N')                     AS mt_montana_club
  , NVL(target_audience.nc_charlotte_nc_club, 'N')                AS nc_charlotte_nc_club
  , NVL(target_audience.nc_raleighdurham_club, 'N')               AS nc_raleighdurham_club
  , NVL(target_audience.nd_north_dakota_club, 'N')                AS nd_north_dakota_club
  , NVL(target_audience.ne_nebraska_club, 'N')                    AS ne_nebraska_club
  , NVL(target_audience.nh_new_hampshire_club, 'N')               AS nh_new_hampshire_club
  , NVL(target_audience.nj_central_new_jersey, 'N')               AS nj_central_new_jersey
  , NVL(target_audience.nj_monmouth__ocean_co_club, 'N')          AS nj_monmouth__ocean_co_club
  , NVL(target_audience.nj_northern_new_jersey_club, 'N')         AS nj_northern_new_jersey_club
  , NVL(target_audience.nj_southern_new_jersey_club, 'N')         AS nj_southern_new_jersey_club
  , NVL(target_audience.njny_gold_coast_club, 'N')                AS njny_gold_coast_club
  , NVL(target_audience.nm_new_mexico_club, 'N')                  AS nm_new_mexico_club
  , NVL(target_audience.nv_nevada_club, 'N')                      AS nv_nevada_club
  , NVL(target_audience.ny_central_new_york_club, 'N')            AS ny_central_new_york_club
  , NVL(target_audience.ny_long_island_club, 'N')                 AS ny_long_island_club
  , NVL(target_audience.ny_new_york_club, 'N')                    AS ny_new_york_club
  , NVL(target_audience.ny_northern_new_york_club, 'N')           AS ny_northern_new_york_club
  , NVL(target_audience.ny_rochester_club, 'N')                   AS ny_rochester_club
  , NVL(target_audience.ny_southern_new_york_club, 'N')           AS ny_southern_new_york_club
  , NVL(target_audience.ny_west_new_york_club, 'N')               AS ny_west_new_york_club
  , NVL(target_audience.ny_westchesterrockland_club, 'N')         AS ny_westchesterrockland_club
  , NVL(target_audience.oh_northern_ohio_club, 'N')               AS oh_northern_ohio_club
  , NVL(target_audience.oh_ohio_valley_club, 'N')                 AS oh_ohio_valley_club
  , NVL(target_audience.ok_oklahoma_club, 'N')                    AS ok_oklahoma_club
  , NVL(target_audience.orwa_pacific_northwest_club, 'N')         AS orwa_pacific_northwest_club
  , NVL(target_audience.pa_home_club, 'N')                        AS pa_home_club
  , NVL(target_audience.pa_lancastercentralyork_club, 'N')        AS pa_lancastercentralyork_club
  , NVL(target_audience.pa_northeast_pennsylvania_club, 'N')      AS pa_northeast_pennsylvania_club
  , NVL(target_audience.pa_northwest_pennsylvania_club, 'N')      AS pa_northwest_pennsylvania_club
  , NVL(target_audience.pa_philadelphia_club, 'N')                AS pa_philadelphia_club
  , NVL(target_audience.pa_pittsburgh_club, 'N')                  AS pa_pittsburgh_club
  , NVL(target_audience.pa_southestern_pa_club, 'N')              AS pa_southestern_pa_club
  , NVL(target_audience.pr_puerto_rico_club, 'N')                 AS pr_puerto_rico_club
  , NVL(target_audience.ri_rhode_island_club, 'N')                AS ri_rhode_island_club
  , NVL(target_audience.sc_coastal_south_carolina_club, 'N')      AS sc_coastal_south_carolina_club
  , NVL(target_audience.sc_greater_columbia_club, 'N')            AS sc_greater_columbia_club
  , NVL(target_audience.scga_savannahhilton_head_cb, 'N')         AS scga_savannahhilton_head_cb
  , NVL(target_audience.scnc_western_carolinas_club, 'N')         AS scnc_western_carolinas_club
  , NVL(target_audience.sd_south_dakota_club, 'N')                AS sd_south_dakota_club
  , NVL(target_audience.tn_tennessee_club, 'N')                   AS tn_tennessee_club
  , NVL(target_audience.tx_austinsan_antonio_tx_club, 'N')        AS tx_austinsan_antonio_tx_club
  , NVL(target_audience.tx_dallasft_worth_tx_club, 'N')           AS tx_dallasft_worth_tx_club
  , NVL(target_audience.tx_houston_tx_club, 'N')                  AS tx_houston_tx_club
  , NVL(target_audience.ut_utah_club, 'N')                        AS ut_utah_club
  , NVL(target_audience.va_virginia_club, 'N')                    AS va_virginia_club
  , NVL(target_audience.vi_virgin_islands_club, 'N')              AS vi_virgin_islands_club
  , NVL(target_audience.vt_vermont_club, 'N')                     AS vt_vermont_club
  , NVL(target_audience.wi_wisconsin_club, 'N')                   AS wi_wisconsin_club
  , NVL(target_audience.wv_west_virginia_club, 'N')               AS wv_west_virginia_club
  , NVL(target_audience.wy_wyoming_club, 'N')                     AS wy_wyoming_club
  , NVL(target_audience.sports_all_members, 'N')                  AS sports_all_members
  , NVL(target_audience.sports_baseball, 'N')                     AS sports_baseball
  , NVL(target_audience.sports_basketball, 'N')                   AS sports_basketball
  , NVL(target_audience.sports_cross_country, 'N')                AS sports_cross_country
  , NVL(target_audience.sports_field_hockey, 'N')                 AS sports_field_hockey
  , NVL(target_audience.sports_football, 'N')                     AS sports_football
  , NVL(target_audience.sports_golf, 'N')                         AS sports_golf
  , NVL(target_audience.sports_lacrosse, 'N')                     AS sports_lacrosse
  , NVL(target_audience.sports_rowing, 'N')                       AS sports_rowing
  , NVL(target_audience.sports_soccer, 'N')                       AS sports_soccer
  , NVL(target_audience.sports_softball, 'N')                     AS sports_softball
  , NVL(target_audience.sports_swimming__diving, 'N')             AS sports_swimming__diving
  , NVL(target_audience.sports_tennis, 'N')                       AS sports_tennis
  , NVL(target_audience.sports_track__field, 'N')                 AS sports_track__field
  , NVL(target_audience.sports_volleyball, 'N')                   AS sports_volleyball
  , NVL(target_audience.sports_wrestling, 'N')                    AS sports_wrestling
  , NVL2(virtual.evnt_crn, 'Y', 'N')                              AS virtual
  , event_comments.event_comments
FROM
  slbevnt
  , gebfunc
  , ssrmeet
  , event_comments
  , target_audience
  , virtual
  , gerattd
  , timeframe
  , atvfisc
  , gtvfunc
  , stvetyp evnt_etyp
  , stvetyp func_etyp
  , stvcoll
  , stvdept
  , gtvsysi
  , gtvfsta
  , gtvpurp
  , gtvemph
  , gtvrsvp
  , gtvfees
  , gtvmenu
WHERE NVL(gerattd_attendance_ind, 'N') = 'Y'
      AND slbevnt_crn = gebfunc_evnt_crn
      AND gebfunc_evnt_crn = ssrmeet_crn
      AND gebfunc_func_code = ssrmeet_func_code
      AND gebfunc_func_code = gtvfunc_code
      AND gebfunc_evnt_crn = event_comments.evnt_crn (+)
      AND gebfunc_func_code = event_comments.func_code (+)
      AND gebfunc_evnt_crn = target_audience.evnt_crn (+)
      AND gebfunc_func_code = target_audience.func_code (+)
      AND gebfunc_evnt_crn = virtual.evnt_crn (+)
      AND gebfunc_func_code = virtual.func_code (+)
      AND gebfunc_evnt_crn = gerattd_evnt_crn
      AND gebfunc_func_code = gerattd_func_code
      AND TRUNC(ssrmeet_start_date) BETWEEN timeframe.start_date AND timeframe.end_date
      AND TRUNC(ssrmeet_start_date) BETWEEN TRUNC(atvfisc_start_date) AND TRUNC(atvfisc_end_date)
      AND slbevnt_etyp_code = evnt_etyp.stvetyp_code (+)
      AND gebfunc_etyp_code = func_etyp.stvetyp_code (+)
      AND slbevnt_coll_code = stvcoll_code (+)
      AND slbevnt_dept_code = stvdept_code (+)
      AND slbevnt_sysi_code = gtvsysi_code (+)
      AND gebfunc_fsta_code = gtvfsta_code (+)
      AND gebfunc_purp_code = gtvpurp_code (+)
      AND gebfunc_emph_code = gtvemph_code (+)
      AND gerattd_rsvp_code = gtvrsvp_code (+)
      AND gerattd_rsvp_code = gtvfees_code (+)
      AND gerattd_menu_code = gtvmenu_code (+)
ORDER BY
  fy
  , start_date
  , end_date
  , event_short_title
  , func_short_desc
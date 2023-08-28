; Author:      Magda Chudzinska, Sophie Smout, Bernie McConnell, SMRU and Jacob Nabe-Nielsen Aarhus University
; Date:        Created (started) 2018-02-02 in NetLogo 6.0.2
; Description:
; Patch size: 1000 x 1000 m
; time step: 900 s (15 min) but can be changed via slider
; The real world is 185 by 174 km 185 by 174 patches (St Andrews area)
; Model starts xxxxxx and finishes no later than xxxxx

;ask hauls [show (word xcor " " ycor)] ; to get coordinates
extensions
[
  gis                         ; for importing gis layers as background
  profiler                    ; for testing CPU use of various submodels
  time                        ; to have ticks saved as real time
  csv
  nw                          ; netwerk extension to set shorstest path
  matrix
  array
  table
]

globals
[
  ;;; LANDSCAPE RELATED GLOBALS
  dist                        ; file showing distance from the shore
  dist_HA31                    ; distance form each patch haul-out with id 1
  dist_HA32
  dist_HA42
  dist_HA52
  dist_HA62
  dist_HA72
  dist_HA82
  dist_HA92
  dist_HA30
;  dist_HA10
;  dist_HA11
;  dist_HA12
;  dist_HA13
;  dist_HA14
;  dist_HA15
;  dist_HA16
;  dist_HA17
;  dist_HA18
;  dist_HA19
;  dist_HA20
;  dist_HA21
;  dist_HA22
;  dist_HA23
;  dist_HA24
;  dist_HA25
;  dist_HA26
;  dist_HA27
;  dist_HA28
;  dist_HA29
  land                        ; raster (ascii) file with background land map (2 = land 0-1 = water (hsi))
  areas5km                    ; model domain with 5x5 km resolution, used for memory procedure
  areas25km
  hsi25km
  ho-list                     ; list of haul-out positions imported from text as agents, not background raster
  xllcorner                   ; x of lower left corner of the environment, read from ascii file
  yllcorner                   ; y of lower left corner of the environment, read from ascii file
  patch_size                  ; reasambling 1 x 1 km
  dayNumber                   ; consequite days (sunrise to sunrise) of the model
  start-time                  ; the time at which the model starts
  finish-time                 ; the time at the end of model
  current-time                ; the current simulation time
  com_count                   ; for debugging, checking if each movement possibility has a procedure assigned
  water-patches
  land-patches
  water-shore-patches
  shore-patches


  ;;; GENERAL MOVEMENT GLOBALS
  step_length                 ; this is 'normal' step length, defined by average speed, in absence of crw, some sort of 'normal, regular' seal move;  used for seals to avoid land and during some travelling modes
  ;time_step                  ; duration of each time step [seconds]
  ;ref-mem-decay-rate          ; memory decay rate, same rate applies to memorising haul-outs and good foraging spots, the unit is [decay/hour] so has to be adjusted to step length
  expindex                    ; used in combination with BehaviorSpace to define paramet combination

  ;; LAND AVOIDANCE PROCEDURE GLOBALS
  land-distance               ; distance used to calculate how much land is ahead of seals in land avoidance procedure [number of patches]

  ;;; CRW WALK basic movement defined by correlated random walk (CRW) and hsi, based partly on kinesis
  ; movement as used in mackerel model and on my equations

;  avg_Vsus                    ; average sustained swimming velocity, based on the observed value for speed for seals from East coast (data filtered from McClintock paper, see McClintock.R code) [m/s]
;  std_Vsus                    ; standard deviation of the sustained swimming velocity calculated as above [m/s]
  Vmin                        ; minimum swimming speed (bodylenghts/time step) calculated as above [m/s]
  Vmax                        ; maximum swimming speed - set to 2m/s (Bernie, pers conv) [m/s]
  avg_Vsus_corr               ; corrected for time step and patch size [patch/time step]
  std_Vsus_corr               ; corrected for time step and patch size [patch/time step]
  Vmin_corr                   ; corrected for time step and patch size [patch/time step]
  Vmax_corr                   ; corrected for time step and patch size [patch/time step]
;  sigK                        ; shape parameter for habitat suitability dependency ("functional responce") [unitless]
  hsi_opt                     ; optimal habitat, most preferable by seals, if set to 1 (so max hsi, sigK does not matter and seals will always behave differently in habitat close to hab_opt [unitless]
;  b                           ; parameter defining wiggliness of the movement, ranging from -1 to 1 with -1 resulting in zigzagging and 1 going in circles
;  imp_dist                    ; importance of distance to target
;  scale_alpha_k, lambda_rate  ; two parameters defining gamma distribution of step length

  ;; DIVING, FORAGING, HOUL-OUT

   ;DiveDuration                ; average dive duartion in seconds, set as slider but literature says it is about 3 minutes (Lesage 1999, Chudzinska master, Suryan 1998, Bjorge 1995) [s]
   ;search_rate                  ; search efficiency in m2/s so number of m2 a seal can 'scan' to search for fish during one second
;  a_time                       ; a and b are coeeficients defining relartion between something (time since last haul-out, blubber or distance) and probabily of hauling-out (1 / (1 + a*exp(b*something))
;  b_time
;  a_blub
;  b_blub
;  a_dist
;  b_dist
  a_prob
;  b_prob                        ; part of equation used for modelling probability of haulout with various stuff
;  b_prob2                        ; part of equation used for modelling probability of haulout with distance
  ;mean_durHO                    ; mean duration [h] of time hauling out, used in right skewed distribution
  ;multProbHaul                  ; multiplier used in the equation defining probability ofhaulout with distance to coast
; mean_g_per_fish               ; weighted mean based on obsrved diet
; mean_kJ_per_fish              ; weighted mean based on fish diet

  ;; MEMORY
  ; mem_level_passedBy_ho      ; memory level of haul-out sites which were remembered by being passed by, parametereised value
  ; ref-mem-decay-rate         ; decay rate of memorised patches, per time step, parameterised value
  ; haulOut_detection_distance ; distance at which seals are able to detect and remember passed by haul-out sites [km], parameterised value
  ; max_memory_day             ; dusration [days] when memory drops close to 0 (see Figure 5 in Trade)

  patches_km25_maxhsi          ; list of max hsi in each 25km patch and their correpsonding ids
  ;; PRODUCING OUTPUT FILES

  outputFileName_movement      ; NetLogo overwrites files everytime a new model is run, to avoid it I want a new file to be produced with date and time in the file name
  outputFileName_bodyCond      ; NetLogo overwrites files everytime a new model is run, to avoid it I want a new file to be produced with date and time in the file name
  outputFileName_crw           ; NetLogo overwrites files everytime a new model is run, to avoid it I want a new file to be produced with date and time in the file name
  outputFileName_spatial       ; NetLogo overwrites files everytime a new model is run, to avoid it I want a new file to be produced with date and time in the file name
  outputFileName_tripDurationExtend
  outputFileName_ho_list
  outputFileName_depletionStart
  outputFileName_depletionEnd
  outputFileName_hauloutCount
  ;; PARAMETERISATION
  param-data                  ; matrix with all the parameters' combination used in parameterisation

]

patches-own
[
  categ;                      ; land, water, shore (within 1 patch from land), haulout (patches directly underneet haul-outs)
  hsi                         ; habitat suitability index, the higher the better the habitat [unitless, has value between 0 and 1]
  hsi25km_max                 ; max value of hsi, calculated at the begining of the simulations, of patches in 25x25km2
  distance2shore              ; as the name says, distance of each patch to shore (land) [unit xy patches] - this is distance pre-calculated in R
  distance2shoreNL            ; this is the same distance as above but calculated within setup procedure, this slows down the setup A LOT but i keep it in the code to be used in case there is no distance file as input
  Nseals                      ; cumulative number of seals which visited this patch
  #Fish_m2                    ; number of fish per m2 (!), used in energy intake procedure and in the future in habitat depletion
  #FishTotal                  ; total number of fish per patch
  ;distance2HO_NL              ; distance from a given patch to all haulout sites
  distance_HO31                 ; distance from a given patch to haulout sites with id 1, calculated in R, as an input file
  distance_HO32
  distance_HO42
  distance_HO52
  distance_HO62
  distance_HO72
  distance_HO82
  distance_HO92
  distance_HO30
;  distance_HO10
;  distance_HO11
;  distance_HO12
;  distance_HO13
;  distance_HO14
;  distance_HO15
;  distance_HO16
;  distance_HO17
;  distance_HO18
;  distance_HO19
;  distance_HO20
;  distance_HO21
;  distance_HO22
;  distance_HO23
;  distance_HO24
;  distance_HO25
;  distance_HO26
;  distance_HO27
;  distance_HO28
;  distance_HO29
  km_5_id                   ; id of 5km resolution grid
  km_25_id
]

breed [hauls haul]            ; haul-out sites imported as text file, MUST be declared before seals as seals will be on top of haul-out sites
breed [seals seal]
breed [sps sp]                ; agents used to create a linked shore line used in calculating shortest path from seal position to a haul-out site where seal is currently moving


sps-own
[
  path                        ; list of sps on seal's way to the next haul-out
]


hauls-own
[
  ho_id
  ;; based on annual surveys during moutling: Chris and Callan
  ;; for east coast data are from 1996-2015 and for Wash 1996-2016
  mean96_16                   ; mean number of counted seal 1996-2016
  max96_16                    ; max number of counted seal 1996-2016
  min96_16                    ; min number of counted seal 1996-2016
  mean10_16                   ; mean number of counted seal 2010-2016
  max10_16                    ; max number of counted seal 2010-2016
  min10_16                    ; min number of counted seal 2010-2016
  PercAll                     ; percentage of all observed seals for each haul-out for entire survey time
  Perc10_16                   ; percentage of all observed seals for each haul-out for entire survey time for 2010-2015(6)
  numberOfModelledSeals       ; starting number of modelled seals hauling our on each haul-out based on PercAll or Perc10_16
  closest_sp                  ; the shortest path point closest to haul-out
  HONseals                    ; cumulative number of seals hauling out at each haul-out at a given tick
]

seals-own
[
  ;;;; BASIC BODY PARAMETERS

  Blength                     ; body length [cm]
  Tmass                       ; total body mass (not lean body mass)[g]
  LBM                         ; lean body mass [g]
  blubber%                    ; percent of total body mass consisting of blubber - KW
  ResMass                     ; mass of reserves (mainly blubber) Tmass-LBM [g]
  sex                         ; male or female
  stomachCap                  ; stomach capacity calculated from relationship established for harp seals (Christiansen 2004); stomach weight (g) = 3.155*body weight (kg) + 132.930; 100g of stomach = 0.8l = 0.8 kg [g]

  ;;;; LAND AVOIDANCE procedure
  is-land?                    ; used in land avoidance procedure
  head-current                ; finding suitable turning angle
  count-left                  ; deciding the land avoidance direction
  count-right                 ; deciding the land avoidance direction
  #turning-angle
  turning-counter             ; used in avoid land procedure
  avoidance-mode-right        ; used when seals avoid land to the right
  avoidance-mode-left         ; used when seals avoid land to the left
  avoid?                      ; true if seals are trying to avoid land, false if far from land, used for debugging
  my_path                     ; list of sps along which seal goes to the next haul-out site of this seal ancounter land on its way to haul-out
  my_path_outOfBay            ; list of sps along which seal goes out of Firth of Forth if it has not been eating for sometime

  ;;; CRW AND BIASED-CRW MOVEMENT

  prev_angle                  ; turning angle (not heading!) from the last step [degress from the last angle]
  speed                       ; seals' speed (step length during crw) [patches/time step]
  curr_turn_output            ; turning angle in a given time step - for debugging and TRACE only
  hsi_Im_on                   ; value of hsi of a patch seal is on - for debugging and TRACE only
  target                      ; the target point (either haul-out or good patch) towards seal head with biased-crw
  diff_heading                ; local variable saved as seal-own for debugging only, see description of the variable in crw procedure
  turn_hsi                    ; local variable saved as seal-own for debugging only, see description of the variable in crw procedure
  dist2target                 ; local variable saved as seal-own for debugging only, see description of the variable in crw procedure
  turn_bias                   ; local variable saved as seal-own for debugging only, see description of the variable in crw procedure
  turn_crw

  ;;; ENERGETICS
  activity                   ; foraging (F), short resting at sea (SRS) - to empty stomach, long resting at sea (LRS) to digest, haul-out (HO), land avoidance (LA), TR-HO - travelling from the moment seal decides to haul-out to actual houl-out site
  BMR                        ; based on equation by Kleiber et al (1975): BMR = ((70 * Tmass^0.75)*0.004184*time_step/(24*60) ; some papers said it should be lean body mass but in Kleibers it is total body mass which make sense as we would expect super fat seals to have different bmr than skinny seals [MJ/time_step]
  ee                         ; energy expenditure at a given time step [MJ/time_step]
  cum_ee                     ; cumulative energy expenditure, zerod if sth happens XXXXXXXX [MJ]
  daily_ee                   ; daily energy expenditure, zerod at the begininng of each day [MJ/day]
  ei                         ; energy intake at a given time step (calculated from g fish) [MJ/time_step]
  cum_ei                     ; cumulative energy intake, not sure if I will use for anything [MJ]
  daily_ei                   ; daily energy intake, zerod at the begininng of each day [MJ/day]
  DailyNetEnergy             ; net energy intake CALCULATED AT THE BEGGINING OF EACH DAY !!!!!! it would be more logic to do it at the end of each day but coding so the model knows which is the last tick of the day is CPU consuming so lets still to first step of each day [MJ/day]
  DaysWithNegativeDNE        ; a counter showing number of consequitive days with daily net energy <= 0. It is set to zero everytime seals have DNE >=0
  ; FORAGING

  #FishCaught                ; number of fish caught during a dive bout (=one time step)
  fishConsumed_g             ; gram of consumed fish, restarted after each resting event [g]
  TotalfishConsumed_g        ; gram of consumed fish over the model duration [g]
  fishConsumed_g_longResting ; gram of consumed fish, restarted after long-term resting , defines what happens after after seal consume > set  body weight [g]
  my_next_patch              ; patch with larger attraction index towards which seal head after leaving a haul-out site
  my_move_away_patch         ; in 'move away' foraging scenarios this is the patch seals move away from
  foraging_trip#             ; foraging trip is defined as movement between two consequitive haul-out events (including haul-out which proceed the trip). This parameter is a counter and increases by 1 everytime a seal starts a new trip. Used to output trip duration and extend

  ; RESTING MODEL
  ; this is first, simple model trying to establish whyand when seals rest (incl resting at sea)
  durationOfResting             ; once seals are in a resting activity they spend x time resting and not moving [number of time step to rest]

  ; HAUL-OUT BEHAVIOUR
  need2ho_TimeSinceLastHo          ; set to true if seals have to haul-out for skin maintanance [true or false]
  need2ho_BLUBBER_ALLOWS     ; set to true if seals are fat enough to haul-out (false if they should not haul-out because they are too skinny) [true or false]
  need2ho_DIGESTION_NEED     ; set to true if seals have to haul-out to digest, if not they rest at sea [true or false]
  my_next_ha                 ; seal's next haul-out site
  durationSinceLastHa        ; starts when seals leave a haul-out site and is reset when next haul-out starts. Needed in procedure mimicing a need a haul-out unrelated to digestion (skin?) [number of time steps]
  typeOfHa                   ; for debugging, to establish frequency of various haul-out triggers (skin, CONDITION ot digest)
  DurationOfDigestion        ; there is a marked difference between energy expenditure during haul-out and during digestion. If seals decide to digest while hauling out, their energy expendture must equal to ee of digestion not haul-out so the part of haulout needed for digestion hass ee = digestion, rest of the same haul-out ee=ee of haul out (seals which digest while hauling out usually spend longer time than it is neccessary just for digestion [number of time steps]
  distHaulOutDigestion       ; distance to haul-out at the moment seals decide to haul-out for digestion - for debugging only
  visited_ho_list            ; list of haul-out sites (their whos) in the order of visit, used to calculate haulout matrices
  previous_ho                ; my last visited haulout

  ; MEMORY RELATED PARAMETERS
  ; haul-out memory

  mem-haul-ids               ; stores memorised haul-out who
  memory-hauls-list          ; list of memory level of visited/memorised haulouts

  ; patch related memory
  patch-ids-5k               ; stores locations of visited 5X5KM AREAS to whcih a visited patch belongs
  patch-ei-5k                ; list of lists showing all ei obtained on a given square, items from this list decay if stored for longer than what it takes memory to decay to almost zero
  patch-area5id              ; ids of 5x5km squares
  ;patch-ei-5k-days           ; days when a given ei was memorised
  patch-ei-5k-memory         ; memory level of each entry of ei per square
  patch-memory               ; list of memory level of visited/memorised patches
  ;patch-ei-list              ; list of visited patches and their corresponding energy intake (nested list of lists)
  ;patch-hsi-list             ; list of visited patches and their corresponding hsi (nested list of lists)
  patch-ei-table-array       ; a table of visted square ids and corresponding ei
  patch-mem-table-array      ; a table of visted square ids and corresponding memory level
  patch-pxcor-table-array    ; a table of visted square ids and corresponding ei
  patch-pycor-table-array    ; a table of visted square ids and corresponding memory level
  patch-ei-burnin-table-array       ; a table of visted square ids and corresponding ei during burnin period
  ;patch-mem-burnin-table-array      ; a table of visted square ids and corresponding memory level during burnin period - it is not needed because memory does to decay in burn in period
  patch-pxcor-burnin-table-array    ; a table of visted square ids and corresponding ei during burnin period
  patch-pycor-burnin-table-array    ; a table of visted square ids and corresponding memory level during burnin period
]

;;;;;;;;;;;;;;;;;;;;; SETTING ENVIRONMENT, BACKGROUND MAPS AND LANDSCAPE ;;;;;;;;;;;;;;;;;;

to set_landscape
  ; all gis data should be loaded in one procedure, otherwise the projection is different
  ; I have to be sure that all my ascci file start from the same xll and yllcorner
  ;(in gis while processing from feature to raster and than from ratser to ascci I have to specify
  ;in 'environment' 'processing extent' the same extent as study area.
  ;All my ascii data must have the same number of columns and rows (for St Andrews 174 (174km cause one cell = 1000m) columns and 185 rows)
  ;and cell size must be specified to 1000. however in model settings it must be max-pxcor
  ;173 and pycor 184 cause first row/column is 0, not 1 in case of St Andrews
  clear-all
  set patch_size 1000
   set dist gis:load-dataset "Input/land_distance_EastCoast_noBay_large.asc" ; Increase map size to allow the seals to travel further from shore - map extension done by MC - KW
  set land gis:load-dataset "Input/hsi_greySeal_Matt_large.asc" ; hsi with values extrapolated to the little bays
  ;set land gis:load-dataset "Input/GrecianEtAl2018.asc" ; hsi with values extrapolated to the little bays
  set areas5km gis:load-dataset "Input/land_EastCoast_5km_NoBay_large.asc"
  set areas25km gis:load-dataset "Input/land_EastCoast_25km_NoBay_large.asc"
  set hsi25km gis:load-dataset "Input/Maxhsi_EastCoast_25km_NoBay_GreySeal_large.asc"
  set dist_HA31 gis:load-dataset "Input/ho_distance_noBay_id_GS_large31.asc"
  set dist_HA32 gis:load-dataset "Input/ho_distance_noBay_id_GS_large32.asc"
  set dist_HA42 gis:load-dataset "Input/ho_distance_noBay_id_GS_large42.asc"
  set dist_HA52 gis:load-dataset "Input/ho_distance_noBay_id_GS_large52.asc"
  set dist_HA62 gis:load-dataset "Input/ho_distance_noBay_id_GS_large62.asc"
  set dist_HA72 gis:load-dataset "Input/ho_distance_noBay_id_GS_large72.asc"
  set dist_HA82 gis:load-dataset "Input/ho_distance_noBay_id_GS_large82.asc"
  set dist_HA92 gis:load-dataset "Input/ho_distance_noBay_id_GS_large92.asc"
  set dist_HA30 gis:load-dataset "Input/ho_distance_noBay_id_GS_large30.asc"
;  set dist gis:load-dataset "Input/land_distance_EastCoast_noBay.asc"
;  set land gis:load-dataset "Input/hsi_greySeal_Matt.asc" ; hsi with values extrapolated to the little bays
  ;set land gis:load-dataset "Input/GrecianEtAl2018.asc" ; hsi with values extrapolated to the little bays
;  set areas5km gis:load-dataset "Input/land_EastCoast_5km_NoBay.asc"
;  set areas25km gis:load-dataset "Input/land_EastCoast_25km_NoBay.asc"
;  set hsi25km gis:load-dataset "Input/Maxhsi_EastCoast_25km_NoBay_GreySeal.asc"
;  set dist_HA31 gis:load-dataset "Input/ho_distance_noBay_id_GS_31.asc"
;  set dist_HA32 gis:load-dataset "Input/ho_distance_noBay_id_GS_32.asc"
;  set dist_HA42 gis:load-dataset "Input/ho_distance_noBay_id_GS_42.asc"
;  set dist_HA52 gis:load-dataset "Input/ho_distance_noBay_id_GS_52.asc"
;  set dist_HA62 gis:load-dataset "Input/ho_distance_noBay_id_GS_62.asc"
;  set dist_HA72 gis:load-dataset "Input/ho_distance_noBay_id_GS_72.asc"
;  set dist_HA82 gis:load-dataset "Input/ho_distance_noBay_id_GS_82.asc"
;  set dist_HA92 gis:load-dataset "Input/ho_distance_noBay_id_GS_92.asc"
;  set dist_HA30 gis:load-dataset "Input/ho_distance_noBay_id_GS_30.asc"
;  set dist_HA20 gis:load-dataset "Input/ho_distance_noBay_id_20.asc"
;  set dist_HA21 gis:load-dataset "Input/ho_distance_noBay_id_21.asc"
;  set dist_HA22 gis:load-dataset "Input/ho_distance_noBay_id_22.asc"
;  set dist_HA25 gis:load-dataset "Input/ho_distance_noBay_id_25.asc"
;  set dist_HA27 gis:load-dataset "Input/ho_distance_noBay_id_27.asc"
;  set dist_HA28 gis:load-dataset "Input/ho_distance_noBay_id_28.asc"
;  set dist_HA29 gis:load-dataset "Input/ho_distance_noBay_id_29.asc"

  ;gis:set-world-envelope gis:envelope-of land ;(gis:envelope-union-of (gis:envelope-of land) (gis:envelope-of haul-out)) ;I leave haulout code for now in case I decide to keep haulouts as raster
  ;resize-world 0 gis:width-of land 0 gis:height-of land
  gis:apply-raster land hsi
  gis:apply-raster areas5km km_5_id
  gis:apply-raster areas25km km_25_id
  gis:apply-raster hsi25km hsi25km_max
  gis:apply-raster dist distance2shore
  gis:apply-raster dist_HA31 distance_HO31
  gis:apply-raster dist_HA32 distance_HO32
  gis:apply-raster dist_HA42 distance_HO42
  gis:apply-raster dist_HA52 distance_HO52
  gis:apply-raster dist_HA62 distance_HO62
  gis:apply-raster dist_HA72 distance_HO72
  gis:apply-raster dist_HA82 distance_HO82
  gis:apply-raster dist_HA92 distance_HO92
  gis:apply-raster dist_HA30 distance_HO30
;  gis:apply-raster dist_HA20 distance_HO20
;  gis:apply-raster dist_HA21 distance_HO21
;  gis:apply-raster dist_HA22 distance_HO22
;  gis:apply-raster dist_HA25 distance_HO25
;  gis:apply-raster dist_HA27 distance_HO27
;  gis:apply-raster dist_HA28 distance_HO28
;  gis:apply-raster dist_HA29 distance_HO29

  ask patch 98 39 [set hsi 0.01] ;this is Isle of May, I have to remove it for now

    ask patches
  [
    ;let value_l gis:raster-sample land patch pxcor pycor ; I should not use this function as it does something weird to the cells along the coastline
    set distance2shore distance2shore / patch_size
    set distance_HO31 distance_HO31 / patch_size
    set distance_HO32 distance_HO32 / patch_size
    set distance_HO42 distance_HO42 / patch_size
    set distance_HO52 distance_HO52 / patch_size
    set distance_HO62 distance_HO62 / patch_size
    set distance_HO72 distance_HO72 / patch_size
    set distance_HO82 distance_HO82 / patch_size
    set distance_HO92 distance_HO92 / patch_size
    set distance_HO30 distance_HO30 / patch_size
;    set distance_HO20 distance_HO20 / patch_size
;    set distance_HO21 distance_HO21 / patch_size
;    set distance_HO22 distance_HO22 / patch_size
;    set distance_HO25 distance_HO24 / patch_size
;    set distance_HO27 distance_HO27 / patch_size
;    set distance_HO28 distance_HO28 / patch_size
;    set distance_HO29 distance_HO29 / patch_size


    if (Habitat = "GrecianEtAl2018")
    [
    if (hsi = 2)
    [
    set pcolor grey
    set categ "land"
    set hsi25km_max 0
    ]

    if (hsi = 0 or hsi < 0.01 or hsi > 2)
    [
    set pcolor blue
    set categ "water"
    set hsi 0.01
    ]

    if (hsi > 0 and hsi != 2)
    [
      set pcolor scale-color red hsi 1 0
      set categ "water"
    ]

    if (distance2shore <= 1 and categ = "water")
    [
    set categ "shore"
    ;set pcolor blue
    ]
    ]

    if (distance2shore <= 1 and categ = "water")
    [
    set categ "shore"
    ;set pcolor blue
    ]

    set #Fish_m2 hsi * #Fish_hsi_multiplier
    set #FishTotal round (#Fish_m2 * patch_size ^ 2)


  ]

    ;calculate_distance2shore_NetLogo ; use it only if there is no input file with distance to shore
  set water-patches patches with [categ = "water"]
  set land-patches patches with [categ = "land"]
  set water-shore-patches patches with [categ != "land"]
  set shore-patches patches with [categ = "shore"]
  ;calculate_meanHSI_of25km_grids
  set patches_km25_maxhsi table:make
  if (Large_scale_foraging = "Omniscience" or Large_scale_foraging = "OmniSwitch")
  [
    if (Habitat = "GrecianEtAl2018") [makingGlobalList_of_best_patches]
  ]

  reset-ticks

end


;to calculate_meanHSI_of25km_grids
;
;  ; this takes long time and there is something wrong so I will calculate in r for now
;  ; but in the future it must be calculated in here in cases we want to recalculate hsi more than once
;    ask water-patches
;  [
;    let patchHere min-one-of water-patches [distance myself]
;    let id_temp [km_25_id] of patchHere
;    set hsi25 mean [hsi] of patches with [km_25_id = id_temp]
;  ]
;
;end

to makingGlobalList_of_best_patches

  let ids ([km_25_id] of water-shore-patches)
  ;show ids
  foreach ids
  [
    x ->
    let maxHSI [hsi25km_max] of one-of water-shore-patches with [km_25_id = x]
    if maxHSI >= 2 [set maxHSI 0.05] ; due to my poor raster operation in r, some water patches have max hsi set to 2
    table:put patches_km25_maxhsi x maxHSI
  ]

  ; making subset of this table with let's say best 15%
  ; I use the same code as for calculating 15% patches in Burn in scenario

  let squareIDss table:keys patches_km25_maxhsi
  let maxhsi table:values patches_km25_maxhsi
  let temp_hsi_id []

  (foreach maxhsi squareIDss
    [
      [x y] ->
      set  temp_hsi_id lput (list x y) temp_hsi_id
    ]
  )
  ; then we sort these hsi from max to min

  let temp_hsi_id-sorted sort-nested-lists-from-highest temp_hsi_id
  ;show temp_hsi_id-sorted
  ;show length temp_hsi_id-sorted
  ; then we take x% of the patches with highest hsi
  ; because East coast is small and I only have 27 25x25 km patches, I take 90% of patches

  let list-length length temp_hsi_id-sorted
  let x_perc round ((percOfsavedPatchesOmniscience * list-length) / 100)
  if x_perc = 0 [set x_perc 1]

  let x%_temp_hsi_id-sorted sublist temp_hsi_id-sorted 0 x_perc

  ;show (word "15% " fifteen%_temp_ei_id-sorted)

  ; now I have to swap ids and ei so id are first
  let id_hsi_x% []

  (foreach x%_temp_hsi_id-sorted
    [
      [x] ->
      ;show x
      set  id_hsi_x% lput (list (item 1 x) (item 0 x)) id_hsi_x%
    ]
  )

  ;show (word "ids first "id_ei_15%)

  set patches_km25_maxhsi table:from-list id_hsi_x%
  show patches_km25_maxhsi
  show length id_hsi_x%
end


; below os the option to calculate distance to shore if we dont have an input file
to calculate_distance2shore_NetLogo
  ;ask patches with [categ = "water"]
  ask water-patches
  [
  ;let closest_land_patch min-one-of patches with [categ = "land"] [distance myself]
  let closest_land_patch min-one-of land-patches [distance myself]
  set distance2shoreNL distance closest_land_patch
 ]
end

; below is the option to calculate distance to all haulouts if we dont have an input file
;to calculate_distance2haulouts_NetLogo
;  ask patches with [categ = "water"]
;
; [
;  set distance2HO_NL []
;  set distance2HO_NL distance hauls in-radius 200
; ]
;end

to create-shortest-path-link

  ask patches with [categ = "shore"] [ sprout-sps 1 ]
  ask patches with [categ = "water" and distance2shore <= 10] [ sprout-sps 1 ]
  ask sps
  [
    set shape "circle"
    set color white
    set size 0.05
    let close_sp other sps in-radius 2
    create-links-with close_sp
]

  ask links [set color [255 255 255 25]] ; white color in rgb,last number defines transparency
  ;ask links [set color "grey"]
end


to load_haulouts

  ; values for East Coast
  set xllcorner 428292.82952926 ; define the lower left corner (from ascii files), only used if we import agent positions from text
  set yllcorner 6186935.0366855

 ; values for Wash
;  set xllcorner 681440 ; define the lower left corner (from ascii files), only used if we import agent positions from text
;  set yllcorner 5796687

  file-open "Input/HaulOutPoints_EastCoast_GreySeal.txt"
  while [ not file-at-end? ]
  [
    create-hauls 1
    [
      set ho_id file-read
      let x-coord_temp (file-read - xllcorner) / patch_size
      let y-coord_temp (file-read - yllcorner) / patch_size
      setxy x-coord_temp y-coord_temp
      set mean96_16 file-read
      set min96_16 file-read
      set max96_16 file-read
      set mean10_16 file-read
      set min10_16 file-read
      set max10_16 file-read
      set PercALl file-read
      set Perc10_16 file-read
      set size 3
      set color green
      set shape "square"
      move-to min-one-of sps [distance myself]
      set closest_sp min-one-of sps [distance myself]
      ;move-to min-one-of patches with [categ != "land"] [distance myself] ; haulout coordinates are based on 1*1 km grids and the coordinates are center of the grid. So to make sure that I dont have haul-out on land, I move them to water
      let patch_under_haul_out min-one-of patches [distance myself]
      ask patch_under_haul_out [set hsi 0 set categ "haulout" set #Fish_m2 0 set #FishTotal 0] ; patches on which haul-out sites are cannot have high hsi. Seals cannot forage on haul-out sites. It is maybe not super biologically correct but otherwise there are errors
    ]
  ]
  file-close

end


;;;;;;;;;;;;;;;;;;;;;  CREATING SEALS ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to create_seals
  ;reset-timer
    create-seals n_of_seals
  [
    ifelse Habitat = "Uniform" [set color black] [set color yellow]
    setxy  0 0 ; starting point does not matter cause I tell seals to move to a another haulout immidietaly
    ;move-to one-of hauls

    ; the who of various agents depands on the order of agents being created.
    ; I want my seals to move to haul-out with min who
    let HA_with_lowestWHO item 0 ([who] of hauls with [who = min [who] of hauls])
    move-to haul HA_with_lowestWHO

    set size 3
    set is-land? []
    set avoidance-mode-right false
    set avoidance-mode-left false
    set heading towards one-of other patches with [categ != "land"] ; because haul-out are in the middle of 5x5 km grid cell, they sometime end up on land and we want all houl-out site to be close to the shore, or at least not inland
    set speed random-normal avg_Vsus_corr std_Vsus_corr
    set prev_angle random-normal 0  std_dir
    ;if show_pd? [pd]
    let ran random 2
    ifelse ran = 0 [set sex "M"] [set sex "F"]
    ifelse sex = "M" [set Blength random-normal 189.1 7.85  set Tmass (1.3 * Blength) * 1000] ;  length distribution of adult (>5 years) males seals and length to body mass ratio - Twiss 2000  - KW
                     [set Blength random-normal 146.5 12.79 set Tmass (1.03 * Blength) * 1000] ;  length distribution of adult (>5 years) females seals and length to body mass ratio - Twiss 2000  - KW
    ;ifelse sex = "M" [set Blength random-normal 146 6.86   set Tmass (1.27 * Blength - 102.67) * 1000] ; length distribution of adult (>5 years) males seals - data from Aisla
     ;                [set Blength random-normal 138.6 5.7 set Tmass (0.86 * Blength - 49.67) * 1000] ; length distribution of adult (>5 years) females seals - data from Aisla
    ;set Tmass (1.13 * Blength - 84.5) * 1000      ; relationship between body length and total mass for adults (>5) in August - data from Aisla
    ;set Tmass (1.15 * Blength - 91.23) * 1000      ; relationship between body length and total mass for adults (>5) in sept october - data from Aisla + NL

    ;let blubber% random (32 - 22) + 23 ; starting reserve mass is between 23-32% (Markussen 1992, Sparling 2006, Beltran 2017, Moray Firth data with blubber% calculated based on Ryg's equation (1990) for P. hispida)
    ifelse sex = "M" [set blubber% random-normal 22.5 5.85] ; starting male reserve mass is between 22.5% +/- 5.85% (Beck 2003) - KW
                     [set blubber% random-normal 21.4 4.92] ; starting female reserve mass is between 21.4% +/- 4.92% (Beck 2003) - KW
    set ResMass (Tmass * blubber%) / 100
    set LBM (Tmass - ResMass)
    set BMR ((70 * (Tmass / 1000) ^ 0.75) * 0.004184 * (time_step / 60)) / (24 * 60) ; Tmass must be in kg and time in mins
    set stomachCap (3.155 * (Tmass / 1000) + 132.930) * 8 ; 0.8/100*1000; Christansen 2004, see description in seal-own attributes
    set fishConsumed_g 0                 ; seals start with empty stomach on a houl-out site at the beginning of day
    set TotalfishConsumed_g 0
    set fishConsumed_g_longResting 0
    set activity "HO"
    set ee 0
    set daily_ee 0
    set daily_ei 0
    set cum_ee 0
    set need2ho_DIGESTION_NEED false
    set need2ho_TimeSinceLastHo false
    set need2ho_BLUBBER_ALLOWS false
    set my_path []
    set my_path_outOfBay []
    set foraging_trip# 1
    set shape "seal" ; unneccessary fun with the icon shape - KW
    set size 10 ; increased icon size for visibility - KW
  ]

  distribute_seals_on_haulouts
  ;move-to haul 3124
  ; this has to be after seals are distributed over haul-out
  ask seals
  [

  ;let ha_im_on min-one-of hauls [distance myself]
  let all_hauls [who] of hauls

  set mem-haul-ids all_hauls;nearby_ha
  ; storing as list made from array to speed up, all haulout sites have low memory value, this value will be updated once a site is used or passed by
  let lengthMem length mem-haul-ids
  let temp_array array:from-list n-values lengthMem [0.20]
  set memory-hauls-list array:to-list temp_array

 ; the site I am on should have the highest memory value (0.99)
  let ha_im_on min-one-of hauls [distance myself]
  let who_ha_im_on [who] of ha_im_on

    if member? who_ha_im_on mem-haul-ids
    [
      let ha_im_on_pos position who_ha_im_on mem-haul-ids
      ;let temp_array2 array:from-list memory-hauls-list
      ;array:set temp_array2 ha_im_on_pos 0.99
      set memory-hauls-list replace-item ha_im_on_pos memory-hauls-list 0.99
      ;set memory-hauls-list array:to-list temp_array2
    ]


  ;show memory-hauls-list
  ;show memory-hauls-array
  ;set memory-hauls-list fput 0.99 memory-hauls-list
  ;set patch-ids []
  set patch-memory []
  ;set patch-hsi-list []
  ;set patch-ei-list []
  set my_next_patch  "none"
  set my_move_away_patch "none"
  set my_next_ha ha_im_on ; this has no influence on seal movement, it is used to output coordinates of the starting haul-out for calculating trip duration and extend
  set visited_ho_list []
  set visited_ho_list lput [who] of ha_im_on visited_ho_list
  set previous_ho ha_im_on
  set patch-ei-table-array table:make
  set patch-mem-table-array table:make
  set patch-pxcor-table-array table:make
  set patch-pycor-table-array table:make

  set patch-ei-burnin-table-array table:make
  set patch-pxcor-burnin-table-array table:make
  set patch-pycor-burnin-table-array table:make

  ;print table:length patch-mem-table-array
  ;set patch-ei-5k-days []
  set patch-ids-5k    []
  set patch-ei-5k []
  set patch-area5id   []
  set patch-ei-5k-memory []

  ]
  ;show timer
end

  ; seals should not be randomly distributed over the haul-out but according to proportion observed in reality. So most popular haulout sites should
  ; have most seals at the start
  ; procedure distributing seals at the start of the model over haul-outs in accordance to number of observed seals on various haul-outs


to distribute_seals_on_haulouts

  ;let list_ha n-values count hauls [i -> i]

  ask seals
  [
    carefully
    [
      loop
      [
        let countSeals count seals-here
        let whoOf_ho_im_on [who] of hauls-here ;ho-Im-on
        ifelse countSeals > [numberOfModelledSeals] of hauls-here
        [move-to one-of hauls with [who = (whoOf_ho_im_on + 1)]]
        [stop]
      ]
    ]
    [
      move-to one-of hauls
    ]
  ]

    ask hauls
  [
    set HONseals count seals-here
  ]


end

;;;;;;;;;;;;;;;;;;;;; LAND AVOIDANCE, COPIED FROM SAIMAA SEAL MODEL ;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;   Liukonnen et al, 2018, Ecological Modelling 368 (2018) 321–335


to  see-if-land-ahead

  ; the list can be extended if measuring land at step length is not enough but shorter/longer distancees are needed
  ; but the more to calculate the slower the model so I need to find a compromise
  ; in Saimaa model short distances where added, in our case these have to be dostances > step_length/patch_size
  let check-land-distances
  (list
  ;(step_length / patch_size)
  (step_length  * 2) ;/ patch_size
  (step_length  * 3) ;/ patch_size
  (step_length  * 4) ;/ patch_size
  )
  set is-land? []

  foreach check-land-distances
  [

x ->
    carefully
     [
      if ([categ = "land"] of patch-ahead x)
      [
       set is-land? fput true is-land?
       ;set pcolor green
       ;print word [who] of seals "land ahead"
      ]
     ]


    [
      set is-land? fput true is-land? ; originally should be true, if I change it to false nothing happens anyway
      ;print word [who] of seals "debug land nobody"
    ]
  ]


end

;----------------------------------------------------
to avoid-land-decision

 ;; AVOIDANCE MODE DECISION
 ;; Count how much is land in each side of the individual and avoid land where less

  set head-current heading

  ;; right side
  set count-right 0
  set heading (head-current + 45)
  set count-right count patches in-cone land-distance 90 with [categ = "land"]

  ;; left side
  set count-left 0
  set heading (head-current - 45)
  set count-left count patches in-cone land-distance 90 with [categ = "land"]

  ;; if equal amount of land patches, random selection
  if count-left = count-right
  [
    ifelse random-float 1 > 0.5
    [
      set count-right count-right + 1
      stop
    ]
    [
      set count-left count-left + 1
    ]
  ]

  ifelse count-left < count-right
  [
    ; choose left direction
    set heading head-current
      set avoidance-mode-right false
      set avoidance-mode-left true
    avoid-land-left
    stop
    ;print(word [who] of seals with [avoidance-mode-left = true] "avoiding left")
  ]
  [
    ; choose right direction
    set heading head-current
      set avoidance-mode-right true
      set avoidance-mode-left false
    avoid-land-right
    ;print(word [who] of seals with [avoidance-mode-right = true] "avoiding left")
  ]


end
;
;;;------------------------------
;
to avoid-land-left
;; After making a decision to left OR in left land avoidance mode
;; Adjusts the turning angle for avoiding the land from left side

  set heading head-current
  see-if-land-ahead
  set turning-counter 0
  ifelse empty? is-land?
  ;; If no land ahead and the seal wants to avoid the land from the left side,
  ;; it needs to turn "right" from current heading (clockwise). (makes the seal to follow the shoreline of an island)
  [
    let head-prev head-current
    loop
    [
      turn-clockwise
      ;print (word [who] of seals "turning clockwise")
      see-if-land-ahead
      if not empty? is-land?
      [
        ;; When land encountered while turning, the seal knows that no more turning is needed
        ;; and it can choose the turning angle as the previous one
        set heading head-prev
        fd step_length; / patch_size
        stop
      ]
      if turning-counter = 36
      [
        ;; The step lenght is sometimes too small to encounter land.
        ;; If this happens, the seal can simply take a step towards target
        ;; as it went around 360 degrees already
        ;;let patch_Notland min-one-of patches with [categ != "land"] [distance myself]
        let patch_Notland min-one-of water-shore-patches [distance myself]
        set heading towards patch_Notland
        ;print (word [who] of seals with [turning-counter = 36] "1turning counter = 36")
        move-to patch_Notland ; sometimes seals move on land during crw procedure and get stuck there, this is to 'force' them to move back to water
        ;fd step_length; / patch_size
        stop
      ]

      set head-prev head-current
    ]
  ]
  ;; If there is land ahead and the seal wants to avoid the land from the left side,
  ;; it needs to turn left from current heading.
  [
    loop
    [
      turn-anticlockwise
      ;print (word [who] of seals "turning anticlockwise")
      see-if-land-ahead
      if empty? is-land?
      [
        ;; When land is not encountered while turning, the seal knows that no more turning is needed
        fd step_length; / patch_size
        stop
      ]

      if turning-counter = 36
      [
        ;; The step lenght is sometimes too large and encounters always land.
        ;; If this happens, the seal can decrease the step lenght and start trying again
        ;print (word [who] of seals with [turning-counter = 36] "2turning counter = 36")
        set turning-counter 0
        ;fd step_length * 0.5
        ;set step_length step_length * 0.5
        ;move-to min-one-of patches with [categ != "land"] [distance myself] ; sometimes seals move on land during crw procedure and get stuck there, this is to 'force' them to move back to water
        move-to min-one-of water-shore-patches [distance myself] ; sometimes seals move on land during crw procedure and get stuck there, this is to 'force' them to move back to water
      ]
    ]
  ]

end
;
;;------------------------------
;
to avoid-land-right
;; After making a decision to right OR in right land avoidance mode (possible in movement modes 1 and 3, targeted moving)
;; Adjusts the turning angle for avoiding the land from the right side

  set heading head-current
  see-if-land-ahead
  set turning-counter 0
  ifelse empty? is-land?
  ;; If no land ahead and the seal wants to avoid the land from the right side,
  ;; it needs to turn "left" from current heading. (makes the seal to follow the shoreline of an island)
  [
    let head-prev head-current
    loop
    [
      turn-anticlockwise
      ;print (word [who] of seals "turning anticlockwise")
      see-if-land-ahead
      if not empty? is-land?
      [
        ;; When land encountered while turning, the seal knows that no more turning is needed
        ;; and it can choose the turning angle as the previous one
        set heading head-prev
        fd step_length; / patch_size
        stop
      ]

      if turning-counter = 36 ; 360 degrees
      [
        ;; The step lenght is sometimes too small to encounter land.
        ;; If this happens, the seal can simply take a step towards target
        ;; as it went around 360 degrees already
        ;;let patch_land min-one-of patches with [categ != "land"] [distance myself]
        let patch_land min-one-of water-shore-patches [distance myself]
        set heading towards patch_land
        ;set heading towards target
        ;print (word [who] of seals with [turning-counter = 36] "3turning counter = 36")
        move-to patch_land ; sometimes seals move on land during crw procedure and get stuck there, this is to 'force' them to move back to water
        ;fd step_length; / patch_size

        stop
      ]
      set head-prev head-current
    ]
  ]
  ;; If there is land ahead and the seal wants to avoid the land from the right side,
  ;; it needs to turn right from current heading.
  [
    loop
    [
      turn-clockwise
      ;print (word [who] of seals "turning clockwise")
      see-if-land-ahead

      if empty? is-land?
      [
        ;; When land is not encountered while turning, the seal knows that no more turning is needed
        fd step_length; / patch_size
        stop
      ]

      if turning-counter = 36
      [
        ;; The step lenght is sometimes too large and encounters always land.
        ;; If this happens, the seal can decrease the step lenght and start trying again
        ;print (word [who] of seals with [turning-counter = 36] "4turning counter = 36")
        set turning-counter 0
        ;;move-to min-one-of patches with [categ != "land"] [distance myself] ; sometimes seals move on land during crw procedure and get stuck there, this is to 'force' them to move back to water
        move-to min-one-of water-shore-patches [distance myself] ; sometimes seals move on land during crw procedure and get stuck there, this is to 'force' them to move back to water
        ;fd step_length * 0.5
        ;set step_length step_length * 0.5
      ]
    ]
  ]

end

;-----------------------------------------------------------------------------

to turn-clockwise
  ifelse head-current < 350
      [
        set head-current head-current + 10
        ;; increase of 10 degrees to the right
      ]
      [
        set head-current (head-current - 350)
        ;; if the degrees are already more than 350, need to start the circle from 0
      ]
      set turning-counter turning-counter + 1
      set heading head-current

end
;-----------------------------------------------------------------------------

to turn-anticlockwise

  ifelse head-current >= 10
      [
        set head-current head-current - 10
        ;; increase of 10 degrees to the left
      ]
      [
        set head-current 360 - (10 - head-current)
        ;; if the degrees are already less than 10, need to start the circle from 360
      ]
      set turning-counter turning-counter + 1
      set heading head-current

end
;;----------------------------------------------------------------------------

; below procedure is used when seals are on their way to haul-out but have to avoid land in between
to go2_houl_out_along_shortest_path

    let sp-Im_on min-one-of sps [distance myself]
    ;let my_next_ha_temp haul my_next_ha
    let target_sp [closest_sp] of my_next_ha
    ;show target_sp
    ;let target_sp min-one-of sps [distance myself]
    ;let target_sp sps with [who = target_sp] ;one-of sps in-radius 1 my_next_ha
    ask sp-Im_on [set path nw:turtles-on-path-to target_sp]
    set my_path [path] of sp-Im_on
    ;show my_path
end

to go_outOfBay_along_shortest_path

    let sp-Im_on min-one-of sps [distance myself]
    ;let my_next_ha_temp haul my_next_ha

    let my_outOfBay_patch patch 105 34 ;  this is a patch just outside Firth of Forth
    if ycor > 52 [set my_outOfBay_patch patch 82 69];  this is a patch just outside Eden Bay
    let target_sp 1 ; just to create target sp
    ask my_outOfBay_patch [set target_sp min-one-of sps [distance myself]]
    ask sp-Im_on [set path nw:turtles-on-path-to target_sp]
    set my_path_outOfBay [path] of sp-Im_on
    ;show my_path
end

;;;;;;;;;;;;;;;;;;;;;; SETUP PROCEDURES ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


to setups
  ;ca

  ;; GENERAL SETUP
  ;testing_b_prob
  set time_step time_step ; time step in the model can vary, we agreed on 15-20 min
  set dayNumber 1
  set monthsOfSimulation monthsOfSimulation
  if (Large_scale_foraging = "BurnInTime" ) [set monthsOfSimulation monthsOfSimulation + monthsOfBurnIn] ; in burn-in scenarios, first two months of simulations are times for seal to learn the area
  ;parameterisation_energyIntake_geert
  ;parameterisation_memories
  ;parameterisation_global_6param


  ;; LANDSCAPE RELATED SETUPS
  set_landscape
  create-shortest-path-link
  load_haulouts
  ;ask patches with [categ = "water"] [set hsi random-float (0.8 - 0.6) + 0.6 set pcolor white] ; for testing bias only

  ;set expindex "none"
  ;; MOVEMENT/CRW RELATED SETUP
  ;set b -0.20
  set Vmin 0.0002
  set Vmax 2
  set step_length (avg_Vsus * time_step / patch_size) ; this is 'normal' step length, defined by average speed, in absence of crw, some sort of 'normal, regular' seal move
  set land-distance (step_length / patch_size) ; for now for land avoidance procedure seals will count how much land is ahead with one time_step
  set avg_Vsus_corr (avg_Vsus * time_step / patch_size)   ; same as step_length. Unneccessary duplicated but I leave it for now
  set std_Vsus_corr (std_Vsus * time_step / patch_size)
  set Vmin_corr (Vmin * time_step / patch_size)
  set Vmax_corr (Vmax * time_step / patch_size)
  set hsi_opt 1
  set ref-mem-decay-rate ref-mem-decay-rate ; with this value = 0.21 as in Saimaa model, seals need approx 4.6 days to forget a haul-out site
  set haulOut_detection_distance haulOut_detection_distance / patch_size * 1000
;  set a_dist 1000
;  set b_dist -0.15
;  set a_time 1000
;  set b_time -0.15
;  set a_blub 1000
;  set b_blub -0.15
  set a_prob 1000

  ;; HAUL-OUT RELATED SETUP

  ask hauls
  [
    set numberOfModelledSeals round (perc10_16 * n_of_seals)
  ]

  create_seals
  let date (remove "-" (substring date-and-time 16 23))
  let time (remove "." remove ":" remove " " (substring date-and-time 0 8))


  ; final-many replicates
   let hb "?"
   ;let dodatek "uniform004"
   ifelse hab_depletion [set hb "Depl"] [set hb "NoDepl"]
   let BSN behaviorspace-run-number ; file names have to be unique. I thought date and time would so the job but sometimes simulations might have been started at the same time and the simulation never ended (at least I think this was the reason)
   set outputFileName_movement (word "Output/Movement.csv")
   set outputFileName_spatial (word "Output/SpatialP2.csv")
   set outputFileName_ho_list (word "Output/HaulOut.csv")
   set outputFileName_depletionEnd (word "Output//DepletionEnd.csv")
   set outputFileName_depletionStart (word "Output/DepletionStart.csv")
  ;; Initialize the simulation time
  set start-time time:create "2018/04/30 00:00"
  ;set finish-time time:create "2018/11/01 00:00"
  set finish-time time:create end-time
  ;set start-time time:show time "yyyy/MM/dd HH:mm"
  set current-time time:anchor-to-ticks start-time time_step "second"

 end


;;;;;;;;;;;;;;;;;; MOVING PROCEDURE ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


to vectorised_crw

  ;;;;; turning angle defined by the relationship with HSI and CRW
  let closest_patch min-one-of patches [distance myself]
  let hsi_on [hsi] of closest_patch
  let expHSI exp((-0.5)*(((hsi_on - hsi_opt) / sigK) ^ 2))
  set turn_hsi random-normal 0  (std_dir * expHSI)
  ifelse turn_hsi > 180
  [set turn_hsi turn_hsi - 360] [set turn_hsi turn_hsi]
  ifelse turn_hsi < -180
  [set turn_hsi turn_hsi + 360] [set turn_hsi turn_hsi]

  set turn_crw (b * prev_angle + turn_hsi)

  ;;;;; turning angle defined by the bias towards a target

  ifelse activity = "TR-HO"
  [set target my_next_ha ]; if target != "none" [ask my_next_ha [set color black]]]
  [set target my_next_patch ]; if target != "none" [ask my_next_patch [ask neighbors [set pcolor blue]]]]

  if target = closest_patch [set target one-of neighbors]; show "te dziwy crw"] ; sometimes when seals are on the target the diff_heading cannot be calculated and the model crashes
  if target != "none" [set diff_heading subtract-headings towards target heading] ; smallest difference between seal current heading and heading towards the target

  ;;;;; step_length = speed - not related to previous step length and not a resultand of vector sum
  ;;;;; related to HSI. Observed data show that speed has ln distribution
  ; NetLogo does not have build-in function for log normal distribution so I follow description by Steve Railbeck

  ; old log normal distribution
;  let beta ln(1 + ((std_Vsus_corr ^ 2) / (avg_Vsus_corr ^ 2)))
;  let em ln (avg_Vsus_corr) - (beta / 2)
;  let es sqrt(beta)
;  let rSpeed exp (random-normal (em * expHSI) es)
  ;show rSpeed

;  while [ rSpeed > Vmax_corr] [set rSpeed exp (random-normal (em * expHSI) es)]
;  while [rSpeed < Vmin_corr] [set rSpeed Vmin_corr]

  ;;;; gamma distribution
  let rSpeed random-gamma (scale_alpha_k / expHSI) lambda_rate
  while [ rSpeed > Vmax_corr] [set rSpeed random-gamma (scale_alpha_k) lambda_rate ];show (word "speed=" rSpeed)]
  while [rSpeed < Vmin_corr] [set rSpeed Vmin_corr]

  set speed rSpeed
  ;show (word rSpeed "after loop")

  ;;;;; resultant turning angle - it is not a sum of three turning angles but trigonometric sum,
  ;;;;; I assume that vectors turn_crw and hsi_on have the same length = 1 but length of vector turn_bias changes with hsi
  ;;;;; so in good habitats (closer to 1) this vector will be shorter and therefore 'pull' less towards the target
  ;;;;; in bed habitat (close to zero), the length will be close to 1 and therefore have the same 'pulling' force as the other two vectors

  ;;;;; solution from jacob

  let old_heading heading
  rt turn_crw
  let crw_dx dx * speed
  let crw_dy dy * speed
  set heading old_heading

  let bias_dx 0
  let bias_dy 0

  ifelse target != "none"
  [
  set diff_heading subtract-headings towards target heading ; smallest difference between seal current heading and heading towards the target
  face target
  set dist2target (distance target)
  set bias_dx (dx * speed) / (dist2target * expHSI * imp_dist)
  set bias_dy (dy * speed) / (dist2target * expHSI * imp_dist)
  set heading old_heading
    ;show (word "bx "bias_dx " by " bias_dy)
  ]
  [
    if (Large_scale_foraging = "MoveAway")
      [
    if my_move_away_patch != "none"
    [
      set diff_heading subtract-headings towards my_move_away_patch heading ; smallest difference between seal current heading and heading towards the target
      face my_move_away_patch
      rt (random (270 - 100) + 100) ; this is the away part, it cannot be strictly 180 because they get stuck in the bays
      set dist2target (distance my_move_away_patch)
      set bias_dx (dx * speed) / (expHSI)
      set bias_dy (dy * speed) / (expHSI)
      set heading old_heading
      ask my_move_away_patch [set pcolor red]
      ;show "I am moving away"
    ]
    ]

  ]


  let new_dx crw_dx + bias_dx
  let new_dy crw_dy + bias_dy

  facexy ( new_dx + xcor ) ( new_dy + ycor )
  let res_turn subtract-headings old_heading heading
  ifelse res_turn > 180
  [set res_turn res_turn - 360] [set res_turn res_turn]
  ifelse res_turn < -180
  [set res_turn res_turn + 360] [set res_turn res_turn]

  set prev_angle res_turn
  fd speed

end


to go

   ;if time:is-after current-time finish-time or time:is-equal current-time finish-time


    if ticks >= round(2881 * monthsOfSimulation)
  [
    if Large_scale_foraging = "BurnInTime" [set monthsOfSimulation monthsOfSimulation - monthsOfBurnIn]
    let date (remove "-" (substring date-and-time 16 23))
    let time (remove "." remove ":" remove " " (substring date-and-time 0 8))
    ;let Export_world_FileName (word "Output/Global/Interface/" expindex "-" date time "interface")

    ;print general profiler
;    profiler:stop
;    print profiler:report
;    profiler:reset
;    show ( word " Finished in " timer " sec ")
    if output? [output_spatial stop]
    ;export-world (word Export_world_FileName ".csv") ; this only shows the 'state' of the model at the last tick but also saves all the model settings in one file
    ;export-interface (word Export_world_FileName ".png")



    stop
  ]

  if (ticks = round(2881 * monthsOfBurnIn + monthsOfSimulation) and (Large_scale_foraging = "BurnInTime"))
  [
    if output? [output_spatial]
  ]


  ; date should be reseted when burn in is over - i cannot make it work, fucjk it

;  if (Large_scale_foraging = "burn-in time")
;  [if ticks = round(2881 * monthsOfBurnIn)
;  [ show "wtf" reset-timer set current-time time:anchor-to-ticks start-time time_step "second"]]


 ; profiler:start
  if ticks < 1 [reset-timer]
  ;show timer
  if time:get "hour" current-time = 8 and time:get "minute" current-time = 0
  [
   set dayNumber dayNumber + 1
  ]

  ask seals
  [
    movement
    ;show mem-haul-ids
    ;show memory-hauls-list

;    ifelse ((length patch-memory) = (length patch-ei-list)); = (length patch-hsi-list))
;    [show "ok"]
;    [show "dupa"]
    ;show mem-haul-ids
    ;show visited_ho_list

  ]; show patch-ids] ;fishConsumed_g show TotalfishConsumed_g show fishConsumed_g_longResting]

  ;ask patches with [categ = "water"]
  ask water-patches
  [
    if any? seals-here [set Nseals Nseals + count seals-here]
  ]

  if ticks > 0
  [
    ask hauls
  [
    set HONseals count seals-here
  ]
  ]

;  if ticks >= 0
;  [
   if output? and Large_scale_foraging != "BurnInTime" [output_file_csv_extension] ;output-tripExtentDuration
    ;if output? and Large_scale_foraging ! = "burn-in time"[output_crw]
   if output? and Large_scale_foraging = "BurnInTime" and ticks >= round (2881 * monthsOfBurnIn) [output_file_csv_extension] ;output-tripExtentDuration
   plotting
;  ]
  tick
end



to movement

  if show_pd? [pd]

  ;;; these are procedures which happen every time step, regardless seal bahaviour, position, state etc
  energy_expenditure
  ;net_energy ; tu moze polukac czy warto zmienic
  ; seals don't calculate net energy and therefore don't gain or loose energy during burn-in period
  ; memory of patches also does not decyas during burn-in time
      if (ticks = round(2881 * monthsOfBurnIn) and (Large_scale_foraging = "BurnInTime"))
    [
      calc-attraction-of-patchesAfterBurnIn-5km-table-array
       ; to get seal distribution at the end of burn in
;      show "I am recalculating my next patch"
;      show (word "15 % " patch-ei-burnin-table-array)
;      show (word "My next patch is " my_next_patch)
    ]

  if time:get "hour" current-time = 8 and time:get "minute" current-time = 0
  [
    ifelse (Large_scale_foraging = "BurnInTime")
    [
      ifelse (ticks >= round(2881 * monthsOfBurnIn))
      [
        decay-patch-memory-table-array
        net_energy
        ;show "should be this"
      ]
      [
        net_energy
        ;show "This cannot be2"
      ]

    ]
    [
      ;if time:get "hour" current-time = 8 and time:get "minute" current-time = 0
      decay-patch-memory-table-array  ;decay-patch-memory
      net_energy
      ;show "This cannot be"
    ]
  ]

;  if time:get "hour" current-time = 8 and time:get "minute" current-time = 0
;  [decay-patch-memory-table-array] ;decay-patch-memory
  ;;;

  ;; this little ifelse loop can be written everytime seals got to a haul-out
  ;; but it would mean writing the same code over and over again so it is here even it it looks out of place
  ifelse activity = "HO"
  [
  set durationSinceLastHa 0
  ;remember-haul-outs-Im-currently-hauling-out
  ]
  [
  set durationSinceLastHa durationSinceLastHa + 1
  see-if-land-ahead ; seals don't have to check if there is land ahead when they haul-out
  ;remember-patches ; seals don't have to remember patches when they haul-out
  if durationOfResting = 0 [remember-haul-outs-Im-passing-by remember-patches-5km-table-array] ; remember-patches-5km seals dont have remember patches or haul-out when they rest at sea, this is to safe CPU
  ]

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; HERE STARTS 'MAIN' MOVEMENT PROCEDURE. THE MAIN DIFFERENCES IN MOVEMENT OF SEALS ARE RELATED TO WHETHER SEALS
  ;; ARE ABOUT TO REST/HAUL-OUT AND WHETHER SEALS AVOID LAND OR NOT
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ifelse durationOfResting = 0 ; (if duration of resting is > 0 it means seals are about to rest/haul-out or are resting/hauling out)
  [
    ; this is what happens when seals do not rest (they forage and/or avoid land)

    ifelse empty? is-land?

  ; this is what happens when seals do not avoid land during non-resting time (they forage:)
    [
      set avoidance-mode-left false
      set avoidance-mode-right false
      set avoid? false

      set activity "F"
      set DurationOfDigestion 0 ; just to reset
      ;energy_intake
      energy_intake_poisson
      ;crw
      ;vectorised_crw

      ; this procedure ensures that seals do not get stuck in Firth of Forth
      ifelse empty? my_path_outOfBay
      [
        vectorised_crw
        ;set color yellow
      ]
      [
        let my_next_sp first my_path_outOfBay
        move-to my_next_sp
        set my_path_outOfBay remove-item 0 my_path_outOfBay
        set color red
        face my_next_sp
        ;show "Im moving outside bay 1"

      ]
      set typeOfHa "No" ; just to reset
      should-I-rest
      ;set com_count com_count + 1
    ]

  ; and when they avoid land in normal land procedure, outside resting time

    [
      ; this is normal land avoidance procedure outisde resting time
      fd step_length ; step_length not related to crw, just to avoid land
      avoid-land-decision
      ;energy_intake
      energy_intake_poisson
      set avoid? true
      set activity "LA"
      ;show "Im avoiding land - normal la"
      ;set com_count com_count + 1
    ]
  ]

  [
    ; this is what happens if seals are on their way to haul-out site or rest at sea or haul-out (are about to rest in any form)

    ; this is what happens if they have to avoid land on their way to haul-out site
    if activity = "LA"
    [
     ;energy_intake
      energy_intake_poisson
      ifelse empty? my_path
      [
        let distance2haul distance my_next_ha
        ifelse distance2haul <= (2 * avg_Vsus_corr) ; if seals are super close to haul-out they should just go there
        [
          move-to my_next_ha
          set previous_ho my_next_ha
          set activity "HO"
          set visited_ho_list lput [who] of  my_next_ha visited_ho_list
          set foraging_trip# foraging_trip# + 1
          set DurationOfDigestion round(((47.3 + 12.72 * (100 * fishConsumed_g_longResting / Tmass)) * 60) / time_step) ; has to be in seconds
          set fishConsumed_g_longResting 0
          set fishConsumed_g 0

          ;ifelse not empty? patch-memory [calc-attraction-of-patches-5km calc-attraction-of-patches-5km-table-array] [set my_next_patch "none"] ; to musze zmienic po array
          ifelse (Patch_memory = "on") [calc-attraction-of-patches-5km-table-array] [set my_next_patch "none"] ; just to be able to switch the memory off

          ; in burn-in scenarios, seals move with crw for the duration of burn-in. They don't calculate attraction of patches but they still memorise them
          ifelse (Large_scale_foraging = "BurnInTime")
          [
            ifelse (ticks < 2881 * monthsOfBurnIn)
            [set my_next_patch "none"]; show "I move as crw"]
            [calc-attraction-of-patches-5km-table-array]
          ]
          [
            calc-attraction-of-patches-5km-table-array
          ]

          remember-haul-outs-Im-currently-hauling-out
          ;set patch-hsi-list []
          ;show "I started hauling out - from LA"
          ;set com_count com_count + 1
        ]
        [
          go2_houl_out_along_shortest_path
          ;show "setting shortest path"
          ;set com_count com_count + 1
        ]
      ]
      [
        let my_next_sp first my_path
        move-to my_next_sp
        set my_path remove-item 0 my_path
        ;show my_path
        ;show "shortest path to ho"
        ;set com_count com_count + 1
      ]

    ]

    ; and this is how they get to the haul-out site but do not avoid land
    ; they still evaluate if land is ahead and if there is they avoid land (go one level up in the loop)

    if activity = "TR-HO"

    [
      let patch_im_on patch-here
      ;if [categ] of patch_im_on = "land" [set patch_im_on min-one-of patches with [categ = "shore"] [distance myself]] ; if seals are close to shore, somwtimes patch-here is set to land

      let dist2land [distance2shore] of patch_im_on
      let distance2haul distance my_next_ha

      ifelse dist2land < 5 and distance2haul <= avg_Vsus_corr
      [
        move-to my_next_ha
        set activity "HO"
        set previous_ho my_next_ha
        set foraging_trip# foraging_trip# + 1
        set visited_ho_list lput [who] of my_next_ha visited_ho_list
        set DurationOfDigestion round(((47.3 + 12.72 * (100 * fishConsumed_g_longResting / Tmass)) * 60) / time_step) ; has to be in seconds
        set fishConsumed_g_longResting 0
        set fishConsumed_g 0
        ;ifelse not empty? patch-memory [calc-attraction-of-patches-5km calc-attraction-of-patches-5km-table-array] [set my_next_patch "none"]
        ifelse (Patch_memory = "on") [calc-attraction-of-patches-5km-table-array] [set my_next_patch "none"] ; just to be able to switch the memory off

        ; in burn-in scenarios, seals move with crw for the duration of burn-in. They don't calculate attraction of patches but they still memorise them
        ifelse (Large_scale_foraging = "BurnInTime")
        [
          ifelse (ticks < 2881 * monthsOfBurnIn)
          [set my_next_patch "none"]; show "I move as crw"]
          [calc-attraction-of-patches-5km-table-array]
        ]
        [
          calc-attraction-of-patches-5km-table-array
        ]

        remember-haul-outs-Im-currently-hauling-out
        ;set patch-hsi-list []
        ;show "I started hauling out - no LA"
        ;set com_count com_count + 1

      ]
      [
        ;energy_intake
        energy_intake_poisson
        ifelse empty? is-land?

        [
          set target my_next_ha
          ;crw
          ;vectorised_crw

          ; this procedure ensures that seals do not get stuck in Firth of Forth
          ifelse empty? my_path_outOfBay
          [
            vectorised_crw
            ifelse Habitat = "Uniform" [set color black] [set color yellow]
          ]
          [
            let my_next_sp first my_path_outOfBay
            move-to my_next_sp
            set my_path_outOfBay remove-item 0 my_path_outOfBay
            set color red
            ;show "Im moving outside bay 1"

          ]
          ;show "biased-crw on my way to ha"
          ;show towards my_next_ha
          ;ask haul my_next_ha [set color black]
          ;set com_count com_count + 1
        ]

        [
          ;fd step_length ; step_length not related to crw, just to avoid land
          avoid-land-decision
          set avoid? true
          set activity "LA"
          ;show "Im avoiding land on my way to ha - the shortest path procedure should be next"
        ]

      ]

    ]

    ; and this if they have to rest at sea

    if activity = "SRS" or activity = "LRS"
    [

      set durationOfResting durationOfResting - 1
      set DurationOfDigestion DurationOfDigestion - 1
      if DurationOfDigestion < 0 or DurationOfDigestion > durationOfResting [set DurationOfDigestion 0]
      fd 0 ; seals do not move when resting
      ;if activity = "SRS" [show "Im resting at sea - short"]
      ;if activity = "LRS" [show "Im resting at sea - long"]
      ;set com_count com_count + 1

    ]

    if activity = "HO"
    [
      set durationOfResting durationOfResting - 1
      set DurationOfDigestion DurationOfDigestion - 1
      if DurationOfDigestion < 0 or DurationOfDigestion > durationOfResting [set DurationOfDigestion 0]
      fd 0 ; seals do not move when resting
      ;show "Im hauling out"
      ;set com_count com_count + 1
    ]

  ]

  if ResMass / Tmass <= 0.05 [show "I died as I am too skinny" die] ; 5% taken from Beltran, search for proper reference there

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; HERE ENDS 'MAIN' MOVEMENT PROCEDURE

end

to should-I-rest

;; resting model 1 and 2 only differ in terms of how much fish a seal gets at each dive

  ; ###################################

  ; there are three mechanisms triggering resting, but at the moment we use first two:
  ;1: eating more than 80-100% of stomach volume (we assume 1l of stomach = 1 kg of fish) and then resting 45-60 min, mimicing large-scale resting from Ramasco 2014

  if fishConsumed_g >= (random (stomachCap - (0.8 * stomachCap)) + (0.8 * stomachCap))
  [
    set activity "SRS"
    set fishConsumed_g 0
    set durationOfResting  random (round((60 * 60) / time_step) - round((45 * 60) / time_step) + 1) + round((45 * 60) / time_step) ; seals should rest around 45min - 1h mimicking long- term resting as in Ramasco 2015
    ;show (word "resting > stomach " durationOfResting)
    ;show "only fishConsumed reset"
 ]

  ;2: food consumption since last large resting is more than 7-14% of seal's body weight (large-scale rest) - (Renouf 1991)
  ; the length of resting is proportional to how much they ate (sparling 2007, manual quick and dirty regression from Fig 2b)
  ; if they have to rest more than minimum time observed haul-out (which is 10 min based on analysis of
  ; East coast seals so less than one time step but leave it if it changes in the future),, they go to houl out if within certain distance to a haul-out, otheriwse they stay and rest at sea
  ; if the haul-out, they may have the time spend at haul-out increased in compariosn to time they would otherwise spent resting at sea
  ; this is following the idea that seals would haul-out for additional reason than just digestion
  ; but the time to haul-out cannot be shorter than time otherwise spent resting at sea

  ;let MaxPercBodyWeight (random-float (14 - 7) + 7) ; seals can consume max 7-14% of their body weight daily (Renouf 1991, Rosen 2004, Ryan pers comm) after this they have to rest (long rest)
  let MaxPercBodyWeight (random-float (6 - 4) + 4) ; seals can consume max 4-6% of their body weight daily () after this they have to rest (long rest) - KW
  if fishConsumed_g_longResting >= ((MaxPercBodyWeight * Tmass) / 100)
  [
    set activity "LRS"
    set fishConsumed_g 0
    set DurationOfDigestion round(((47.3 + 12.72 * (100 * fishConsumed_g_longResting / Tmass)) * 60) / time_step) ; has to be in seconds
    set durationOfResting  DurationOfDigestion  ; seals should rest proportionally to what they ate
    ;show (word "resting > % body mass " durationOfResting " " fishConsumed_g_longResting)
    set fishConsumed_g_longResting 0 ;
    ;show "Long resting and fish consumed reset"

  ]
  should-I-haul-out

  ; 3: it is after dark regarding how much seals ate during the day (ok,it has to be more than 0, otherwie they will rest the entire night) (sleeping or longer digestion breaks)
  ; the length of resting is proportional to how much they ate (sparling 2007, manual quick anddirty regression from Fig 2b)
  ;if they have to rest more than minum time observed haul-out (which is 10 min based on analysis of
  ; East coast seals so less than one time step but leave it if it changes in the future), they go to houl out if within 1km to a haul-out, otheriwse they stay and rest at sea
  ; if the y haul-out, they may have the time spend at haul-out increased in compariosn to time they would otherwise spent resting at sea
  ; this is following the idea that seals would haul-out for additional reason than just digestion
  ; but the time to haul-out cannot be shorter than time otherwise spent resting at sea

;    if ticks_since_sunset = 1 and fishConsumed_g_longResting > 0
;  [
;    set activity "LRS"
;    set fishConsumed_g 0
;    let TimeNeeded2Digest (47.3 + 12.72 * (100 * fishConsumed_g_longResting / Tmass)) * 60 ; has to be in seconds
;    set durationOfResting  round(TimeNeeded2Digest / time_step)  ; seals should rest proportionally to what they ate
;    show (word "resting night " durationOfResting " " fishConsumed_g_longResting)
;    set fishConsumed_g_longResting 0
;
;  ]


end


to should-I-haul-out ; in previous model versions there was an older procedure called should-I-haul-out, this is a new version but even if the old one is deleted from this code I keep the name

;Seals haul-out either to digest or to maintain their skin. Whether seals digest at the
;haul-out site (prob_HO_DIGESTION_NEED, need2ho_DIGESTION_NEED) or at sea depends on the distance
;between place where the seal is at the moment it has to digest
;and its next houl-out (calc-attraction-of-hauls defines what is the next haul-out).
;The distance / probability of haul-out relationship is defined by exponential function
; digestion related haul-out always take place regarding skin and blubber condition of seals

; the probability of hauling out increase with time since last haul-out meaning that the longer the seals are at sea, the worse the condition of their skin is
; however, if seals are in bad condition (low blubber %), they are more likely to remain foraging regardless their skin condition (idea taken from fur seal IBM)



  if activity = "F" or activity = "LA" ; when seals already rest (HO, LRS, SRS) they should not reevaluate need of re-resting
  [
  ; #1 calculating probability of haul-out in relation to time since last haul-out
  let prob_HO_TimeSinceLastHo (1 / (0.8 + (5 * exp(-0.7 * (durationSinceLastHa / 328)) * ((durationSinceLastHa / 2) + 1)))) ; increase probability of haul out given time to allow seals to travel for durations up to 855 hrs (per Janeke's Data) - KW
  let bin1 binomial 1 prob_HO_TimeSinceLastHo
  ;let rand1 random-float 1.0001
  ;let rand1 0 ; for debugging
  ;ifelse prob_HO_TimeSinceLastHo >= rand1 [set need2ho_TimeSinceLastHo FALSE] [set need2ho_TimeSinceLastHo TRUE ]  ;print (word precision prob_HO_TimeSinceLastHo 2 " " precision rand1 2 " I should ha skin")
    ifelse bin1 = 1 [set need2ho_TimeSinceLastHo TRUE ] [set need2ho_TimeSinceLastHo FALSE ]  ;print (word precision prob_HO_TimeSinceLastHo 2 " " precision rand1 2 " I should ha skin")
  ;show (word "HO_Time " precision prob_HO_TimeSinceLastHo 3 " " durationSinceLastHa " " bin1)

  ; #2 calculating probability of haul-out in relation to blubber %
  let prob_HO__Blubber (1 / (1 + a_prob * exp(b_prob * ((ResMass * 100) / Tmass))))
  let bin2 binomial 1 prob_HO__Blubber
  ;show (word "HO_Blubber " precision prob_HO__Blubber 3 " " round(((ResMass * 100) / Tmass)) " " bin2)
;  let rand2 random-float 1.0001
;  ifelse prob_HO__Blubber >= rand2 [set need2ho_BLUBBER_ALLOWS FALSE ] [set need2ho_BLUBBER_ALLOWS TRUE] ;print (word precision prob_HO__Blubber 2 " " precision rand2 2 " I should ha I am too fat")
  ifelse bin2 = 1 [set need2ho_BLUBBER_ALLOWS TRUE ] [set need2ho_BLUBBER_ALLOWS FALSE] ;print (word precision prob_HO__Blubber 2 " " precision rand2 2 " I should ha I am too fat")
  ]

  ; #3 calculating probability of haul-out in relation to distance to coast if seal has to digest
  if activity = "LRS"
  [
    calc-attraction-of-hauls
    ;show "attraction because of LRS"
    ;set my_next_ha min-one-of hauls [distance myself]
    let dist_closest_haul-out distance my_next_ha
    set distHaulOutDigestion dist_closest_haul-out ; for debugging only
    ;let prob_HO_DIGESTION_NEED (0.81 ^ dist_closest_haul-out) ; exponential decrease
    ; wersja ze wszystkie trzy grafy maja ten sama formule
    ;let prob_HO_DIGESTION_NEED (1 / (1 + a_prob * 0.0001 * exp(-1 * b_prob * dist_closest_haul-out * 2.4))) ; s-shaped decrease, na oko, jak powyzej
     ; wersja ze ten graf ma inna formule niz time i blubber
    let prob_HO_DIGESTION_NEED 0
    ifelse b_prob >= -0.15
    [set prob_HO_DIGESTION_NEED (1 / (1 + a_prob * 0.0001 * 10 * exp(-1 * b_prob * dist_closest_haul-out * 2.4)))]
    [set prob_HO_DIGESTION_NEED (1 / (1 + a_prob * 0.0001 * exp(-1 * b_prob * dist_closest_haul-out * 2.4)))]
    ; s-shaped decrease, na oko, jak powyzej
    ;let prob_HO_DIGESTION_NEED (1 / (1 + a_prob * 0.0001 * multProbHaul * exp(-1 * b_prob2 * dist_closest_haul-out * 2.4))) ; s-shaped decrease, na oko, jak powyzej
    let bin3 binomial 1 prob_HO_DIGESTION_NEED
    ;show (word "HO_DIGESTION_NEED " precision prob_HO_DIGESTION_NEED 3 " " round(dist_closest_haul-out) " " bin3)
    ;let rand3 random-float 1.00001
    ;let rand3 0 ; this is for debugging, just to see what happens if they haul-out
    ;ifelse prob_HO_DIGESTION_NEED >= rand3 [set need2ho_DIGESTION_NEED FALSE ] [set need2ho_DIGESTION_NEED TRUE ]  ;print (word precision prob_HO_DIGESTION_NEED 2" " precision rand3 2 " I should ha digestion-too close to shore"); print (word precision prob_HO_DIGESTION_NEED 2" " precision rand3 2 " I rest at sea, too far from ha")
    ifelse bin3 = 1 [set need2ho_DIGESTION_NEED TRUE ] [set need2ho_DIGESTION_NEED FALSE]; show "Chyba powinnam rest at sea-long"]  ;print (word precision prob_HO_DIGESTION_NEED 2" " precision rand3 2 " I should ha digestion-too close to shore"); print (word precision prob_HO_DIGESTION_NEED 2" " precision rand3 2 " I rest at sea, too far from ha")
  ]


  ; ################# START DIGESTION + SKIN + BLUBBER #######################

  ifelse need2ho_DIGESTION_NEED = TRUE
  [
    set durationOfResting DurationOfDigestion  ; here durationOfResting = time needed for digestion
    should-I-haul-out-subprocedure
    ;show (word "I am going to " my_next_ha " for digestion, distance allowed" )
    set typeOfHa "DIG"

    stop
    ;print (word need2ho_DIGESTION_NEED " stop in first if")
  ]

  [

    if need2ho_TimeSinceLastHo = TRUE
    [

      ifelse need2ho_BLUBBER_ALLOWS = TRUE
      [
       should-I-haul-out-subprocedure
       ;show (word "Long time since last houlout, I am fat, I am going to " my_next_ha)
       set typeOfHa "SKIN_FAT"
      ]

      [
      ;show "My skin is bad but I am in poor condition and have to forage"
      ;show (word "HO bo czas, not digestion")
      ;show "My skin is bad but I am in poor condition - this should rarely happen now"
      ]

    ]
  ;show (word need2ho_DIGESTION_NEED " nie chce na HO")
  ]

  ; ################# END DIGESTION + SKIN + BLUBBER #######################
end

to should-I-haul-out-subprocedure

; this part of the code has to be repeated several times in should-I-houl-out procedure. In order to reduce the lines, I have it as seperate procedure

  set activity "TR-HO"

  ;for parameterisation
  let haul-out-duration random-normal ln(mean_durHO) ln(3.6)
  set haul-out-duration round(((exp haul-out-duration) * 60 * 60) / time_step)
  ; sometimes the duration gets ridicously long so I cut it at 30h (see fig. 18 in TRACE)
  if haul-out-duration >= (30 * 60 * 60 / time_step) [set haul-out-duration (30 * 60 * 60 / time_step)]; show "My haul-out was super long"]

  ;;;(cunningham 2009, Thompson 1990, Ramasco 2014)
  ;let haul-out-duration round(random-normal round((4.77 * 60 * 60) / time_step) round((3.6 * 60 * 60) / time_step))
  ;let haul-out-duration random-normal ln(4.77) ln(3.6)
  ;set haul-out-duration round(((exp haul-out-duration) * 60 * 60) / time_step)

  ;;; based on data by Isla and Paul incorporating close to shore resting
  ;;; the data are right skewed so I draw the number from normal distribution around log(mean) and log(sd) and then backtransfrom it
  ;let haul-out-duration random-normal 2.21 1.06
  ;set haul-out-duration round(((exp haul-out-duration) * 60 * 60) / time_step)

  calc-attraction-of-hauls
  ;show "attract becuase of TR-HO"
  ;set my_next_ha min-one-of hauls [distance myself]


  set need2ho_TimeSinceLastHo FALSE
  set need2ho_BLUBBER_ALLOWS FALSE
  set need2ho_DIGESTION_NEED FALSE
  ; the time to haul-out cannot be shorter than time otherwise spent resting at sea
  ifelse haul-out-duration > durationOfResting
  [
    set durationOfResting haul-out-duration
  ]
  [
    set durationOfResting durationOfResting
  ]
end

; ################################ MEMORY PROCEDURES ##########################
; #############################################################################

;; ............................................................................
;; HAUL-OUT MEMORY PROCEDURES

; seals remember each haul-out side on which they hauled out (remember-haul-outs-Im-currently-hauling-out)
; they also remember each haul-out passed by distance xxxx, regardless if they hauled out on it or not (remember-haul-outs-Im-passing-by)
; in order for seals to not to calculate distance to haul-out each time step,they only calculate this distance if are within xxxx distance from shore
; but before they remember this passed one, they have to 'pass' 0,1 probability test. The new site is given memory value a bit smaller then 0.99

; IN THE FUTURE WE MAY ADD EXTRA HAUL-OUT QUALITY ASPECTS INFLUENCING PERCEPTION/MEMORY LIKE NUMBER OF OTHER SEALS HAULINGOUT OR AVAILABILITY RELATED TO TIDES

to remember-haul-outs-Im-currently-hauling-out

      let ha_im_on min-one-of hauls [distance myself]
      let who_ha_im_on [who] of ha_im_on
      ;ask ha_im_on [set color black]

  if member? who_ha_im_on mem-haul-ids
      [
        let replace-pos position who_ha_im_on mem-haul-ids
        set memory-hauls-list replace-item replace-pos memory-hauls-list 0.99
        ;set mem-haul-ids remove-item remove-pos mem-haul-ids
      ]

      ;let max-memory-hauls 0.99
      ;set memory-hauls-list fput max-memory-hauls memory-hauls-list
      ;set mem-haul-ids fput [who] of ha_im_on mem-haul-ids
      ;show (word "remIon " memory-hauls-list)
end

to remember-haul-outs-Im-passing-by


  let patch_im_on patch-here
  ;if [categ] of patch_im_on = "land" [set patch_im_on min-one-of patches with [categ = "shore"] [distance myself]] ; if seals are close to shore, somwtimes patch-here is set to land
  let dist2coast [distance2shore] of patch_im_on
  ifelse dist2coast <= haulOut_detection_distance
  [
    ;show "remember-haul-outs-Im-passing-by"
    ;; before any further operations/calculations are done, we first have 0,1 test draw from uniform distribution
    ;; so if the test fails the further code is not implemented saving time

    let rand random 2
    ;show rand
    if rand = 1
    [
      ;let paased_ha hauls in-radius haulOut_detection_distance
      let paased_ha hauls with [distance myself <= haulOut_detection_distance]
      ;show [who] of paased_ha
      if any? paased_ha;hauls in-radius haulOut_detection_distance
      [
      let close_ha one-of paased_ha
      let who_close_ha [who] of close_ha
      ;show who_close_ha
      ;ask close_ha [set color black]
      ;show "I remember this ha"

        if member? who_close_ha mem-haul-ids
        [
          let replace-pos position who_close_ha mem-haul-ids

          ifelse item replace-pos memory-hauls-list > mem_level_passedBy_ho
          [
            set memory-hauls-list replace-item replace-pos memory-hauls-list (item replace-pos memory-hauls-list)
          ]
          [
            set memory-hauls-list replace-item replace-pos memory-hauls-list mem_level_passedBy_ho
          ]

        ]

      ]
    ]
  ]
  [
    stop
  ]
;show mem-haul-ids
;show memory-hauls-list

end


to calc-attraction-of-hauls

  ; this procedure is copied from Saimaa model (Liukkonen et al 2018) and Nabe-Nielsen et al 2013
  ;; calculate distances to each haul out sites from the memory list and calculates attraction based on this distance and memory value of each haul-out
  ;; NOW SEALS' HAUL-OUT MEMORY DOES NOT DECAY SO THE ATRRACTIVENESS IS ONLY BASED ON DISTANCE AND THE PRESET MEMORY VALUES (0.99 FOR ALREADY VISITED HA AND 0.70 FOR PASSED BY)

  ;; for seals from St Andrews this should not be a big deal becuase there is not too many but in areas where mamory list can get long and we won't model many seals, this procedure may, however, really slow down the model
  let dist-to-hauls []
  (foreach mem-haul-ids
    [
      x ->
      set dist-to-hauls lput (precision (distance one-of hauls with [who = x]) 1) dist-to-hauls
      ;;resolution scale transforms the distance to meters, although it does not matter, it can also be in number of patches
    ]
  )

  ;; calculate attraction values
  let profits-list []
  (foreach memory-hauls-list dist-to-hauls
    [
      [x y] ->
      if y = 0 [set y 0.0001]
      set  profits-list lput (precision (x / y) 3)  profits-list
    ]
  )
 ;; atrractiveness is first expressed in %, than multplied by a random number between 0 and 1 and next haul-out site is then the one with lowest attractivenss (same idea as for choosing next target patch)

   let profits-list-% [] ; attractivensss list expressed in %
   let sum_attr sum profits-list

    (foreach profits-list
    [
      x ->
      set profits-list-% lput (precision (100 * x / sum_attr) 3) profits-list-%
    ]
    )

  let profits-list-%_rand [] ; multplying the atrraction list expressed in % by a random number
  let c 0
  (foreach profits-list-%
    [
      x ->
      let randn random-float 1
      ;show randn
      set profits-list-%_rand lput (precision (x * randn) 3) profits-list-%_rand
      ;set profits-list-%_rand lput (precision ((item c profits-list-%) * randn) 2) profits-list-%_rand
      ;set c c + 1
    ]
  )

  let max_att max profits-list-%_rand
  ;show sum attract_list_%_rand
  let max-pos (position max_att profits-list-%_rand)


  ;show (word "ids " mem-haul-ids)
  ;show (word "mem " memory-hauls-list)
  ;show (word "dist " dist-to-hauls)
  ;show (word "attract% " profits-list-%)
  ;show (word "attract %*rand " profits-list-%_rand)
  ;show (word "min is " max_att " at " max-pos)
  ;show visited_ho_list

  let next_who item max-pos mem-haul-ids

  if not empty? mem-haul-ids ; WAS IFELSE BEFORE
  [
  set my_next_ha one-of hauls with [who = next_who]
  ;show (word "my next ha " my_next_ha)
  ]



end

;; ......................................................................................................
;; PATCH - MEMORY PROCEDURES


to decay-patch-memory-table-array

;  show "decaying"
; this procedure is copied from Saimaa model (Liukkonen et al 2018) and Nabe-Nielsen et al 2013
; ref-mem-decay-rate has to be paramaterised and at the moment set to 0.21 (decay rate per hour) as in Saimaa paper

; decay of patch (square) memory
  let seeIfEmpty table:length patch-mem-table-array
  ifelse seeIfEmpty > 0
  [
    ; decaying memory of patches
    let mems table:values patch-mem-table-array
    let squareIDs table:keys patch-mem-table-array
;    show mems
;    show squareIDs

    (foreach mems squareIDs
      [
        [i y]->
        ;increment-table-value patch-mem-table-array x ref-mem-decay-rate
        let newMem (precision (i - ((ref-mem-decay-rate * 96) * (1 - i) * i)) 3)
        table:put patch-mem-table-array y newMem
      ]
    )
  remove_forgotten_patches-5km-table-array
  ]
  [
    stop
  ]




end


to remove_forgotten_patches-5km-table-array

  ; this is only to save CPU as patches with low memory value are not likely to be chosen as next patch to go
  ; remove patches with low memory level

    let mems table:values patch-mem-table-array
    let squareIDs table:keys patch-mem-table-array
    (foreach mems squareIDs
      [
        [i y]->
        ifelse i < 0.01
        [
          table:remove patch-mem-table-array y
          table:remove patch-ei-table-array y
        table:remove patch-pxcor-table-array y
        table:remove patch-pycor-table-array y
        ]
        [
          stop
        ]

      ]
    )


end




to remember-patches-5km-table-array


 ; SEALS DO NOT REMEMBER EACH PATCH BUT A 5X5 KM AREA TO WHCIH THIS PATCH BELONGS
 ; SEALS REMEMBER THE ENERGY INTAKE OF VISITED AREA (which are not land)
 ; FOR REDUCING CPU, ONLY PATCHES WITH HSI > 0.01 ARE STORED, LATER ON IF LIST OF VISITED PATCHES IS EMPTY (SEALS PREVIOUSLY ONLY VISITED CRAP PLACES)
 ; SEALS LEAVE HAUL-OUT WITH CRW AND NOT BIASED CRW TOWARDS A PATCH
 ; IF SEALS ALREADY VISITED THIS AREA BEFORE, ITS MEMORY IS UPDATED AND THE MEMORISED EI IS THE MEAN VALUE OF THE
 ; CURRENT EI AND EI FROM THE LAST VISIT

;patch-all-table-array

      let patch_im_on patch-here
      ;if [categ] of patch_im_on = "land" [set patch_im_on min-one-of patches with [categ = "shore"] [distance myself]] ; if seals are close to shore, somwtimes patch-here is set to land
      if [categ] of patch_im_on = "land" [set patch_im_on min-one-of shore-patches [distance myself]] ; if seals are close to shore, somwtimes patch-here is set to land
      let area5_im_on [km_5_id] of patch_im_on
      let hsi_patch_im_on [hsi] of patch_im_on
      let x [pxcor] of patch_im_on
      let y [pycor] of patch_im_on
      ;show hsi_patch_im_on

  ifelse hsi_patch_im_on > 0.01
  [
    let prevEi table:has-key? patch-ei-table-array area5_im_on
    ifelse prevEi = true
    [
      let prevEi1 table:get patch-ei-table-array area5_im_on
      let prevEi2 fput ei prevEi1
      set prevEi prevEi2
    ]
    [
      set prevEi (list ei)
    ]
    ;print prevEi
    ;if any? prevEi [show prevEi]
    table:put patch-ei-table-array area5_im_on prevEi
    table:put patch-mem-table-array area5_im_on 0.99
    table:put patch-pxcor-table-array area5_im_on x
    table:put patch-pycor-table-array area5_im_on y
  ]
  [
    stop
  ]

;    show patch-ei-table-array
;    show patch-mem-table-array
;    show patch-pxcor-table-array
;    show patch-pycor-table-array
    ;show patch_im_on
    ;show area5_im_on
end


to calc-attraction-of-patches-5km-table-array

   ; attractiveness of patches is calculated based on equation by Jacob (Nabe-Nielsen et al 2013) which takes distance, memory and energy intake
   ; distance is an input. Energy intake is the mean intake from a given square from the last x days. This x = length of memory decay when seals almost forget stuff (with r=0.005 it is a about 16 days) n(max-memory-day)
  ;show "zonk?"
  let seeIfEmpty table:length patch-mem-table-array
  ifelse seeIfEmpty > 0
  [
  let startOne timer
  let ha_id item 0 ([ho_id] of hauls-here)
  let mems table:values patch-mem-table-array
  let squareIDs table:keys patch-mem-table-array
  let eis table:values patch-ei-table-array
  let xpacza table:values patch-pxcor-table-array
  let ypacza table:values patch-pycor-table-array
  let attract_tbl table:make
;  show mems
;  show squareIDs
;  show eis
  let i 0
  ;(foreach mems squareIDs eis
  (foreach mems xpacza ypacza eis squareIDs
    [

      [x x1 y2 z id] ->

      let dist_temp (word "[distance_HO" ha_id "] of" patch x1 y2)
      ;let pacz one-of water-shore-patches with [km_5_id = y]
      ;let dist_temp (word "[distance_HO" ha_id "] of" pacz)
      let dist_temp2 runresult dist_temp
      if dist_temp2 < 1 [set dist_temp2 1]
      ; based on remembered intake rate
      let attract ((x * (mean z)) / dist_temp2)
      table:put attract_tbl id attract
    ]
  )

    ;let endOne timer
    ;let FirstPart endOne - startOne
    ;show (word "First part " FirstPart)
  ;show (word "attr" attract_tbl)
  ;show patch-ei-list


  ;; seal choose their next patch to go based on Jacob's approach: first a random number is generated (between 0 and 1) for each of the patch, then
  ;; it is multiplied by attractiveness expressed as %. Seal takes the patch with highest number after multiplication
  ;let startTwo timer
  let attract table:values attract_tbl
  let sum_attr sum attract
  if sum_attr = 0 [set sum_attr 0.000001] ; seals which spend time in araeas where ei is almost always zero, may end up with atrractiveness = 0 and we cannot divide by 0

  let attract_tbl_% table:make
    (foreach squareIDs attract
    [
      [x y] ->
      let attract% (100 * y / sum_attr) * random-float 1
      table:put attract_tbl_% x attract%
    ]
    )
  ;show attract_list_%
  ;show attract_tbl_%
  let attract_rand table:values attract_tbl_%
  let max_att max attract_rand
  ;show max_att
  let max-pos (position max_att attract_rand)
  let SqID item max-pos squareIDs
  let patch_im_on patch-here

  ;set my_next_patch one-of water-shore-patches with [km_5_id = SqID] ; and [self != patch_im_on]);other patch_im_on with [km_5_id = SqID];and one-of water-shore-patches with [pxcor != [pxcor] of patch-here] one-of water-shore-patches with [km_5_id = SqID]
  set my_next_patch one-of water-shore-patches with [km_5_id = SqID and self != patch_im_on]
    if [pxcor] of my_next_patch > 105 and xcor < 80 [go_outOfBay_along_shortest_path]; show "out patch"]; show "Zaraz wybywam z zatoki"]
    if DaysWithNegativeDNE >= NdaysWithNegativeNEB and xcor < 80 [go_outOfBay_along_shortest_path]; show "out hungry"]; show "Zaraz wybywam z zatoki"]
  ;show [km_5_id] of my_next_patch_temp
;   show (word "attract% " attract_list_%_rand)
;   show (word "max is " max_att " at " max-pos)
;   show (word "patch-ids " patch-ids-5k)
    ;let endTwo timer
    ;let SecPart endTwo - startTwo
    ;show (word "Sec part " SecPart)
  ]
  [
    set my_next_patch "none"
  ]





end

to calc-attraction-of-patchesAfterBurnIn-5km-table-array

  ; after months of burn in, seals store location and id of patches with 15% best energy intake 'learned' during burn in time
  ; they are stored in new tables which after the burn in time stay fixed (unchanged)
  ; after, seals choose one of these 15% at random and put it to their normal memory tables and wipe out the rest,
  ; this one is also set as their target patch - this all happens once at the end of burn in

  ; I first calculates mean ei for each visited square (bo tam jest kilka wartosci)
  ; to do it I make a nested lists with ei and ids of squares as I have a code written for nested lists
  ;show (word "stara przed" patch-ei-table-array)
  let squareIDss table:keys patch-ei-table-array
  ;show squareIDs
  let eiss table:values patch-ei-table-array
  ;show eis
  let temp_ei_id []

  (foreach eiss squareIDss
    [
      [x y] ->
      set  temp_ei_id lput (list (mean x) y) temp_ei_id
    ]
  )
  ;show (word "nowa " temp_ei_id)
  ;show (word "stara po " patch-ei-table-array)

  ; then we sort these ei from max to min

  let temp_ei_id-sorted sort-nested-lists-from-highest temp_ei_id
  ;show (word "sorted " temp_ei_id-sorted)

  ; then we take 15% of the patches with highest ei

  let list-length length temp_ei_id-sorted
  let x_perc round ((percOfsavedPatchesFromBurnIn * list-length) / 100)
  if x_perc = 0 [set x_perc 1]

  let fifteen%_temp_ei_id-sorted sublist temp_ei_id-sorted 0 x_perc

  ;show (word "15% " fifteen%_temp_ei_id-sorted)

  ; now I have to swap ids and ei so id are first
  let id_ei_15% []

  (foreach fifteen%_temp_ei_id-sorted
    [
      [x] ->
      ;show x
      set  id_ei_15% lput (list (item 1 x) (item 0 x)) id_ei_15%
    ]
  )

  ;show (word "ids first "id_ei_15%)

   set patch-ei-burnin-table-array table:from-list id_ei_15%

  ; now I have to match x and y of each id from 15% table

  let ids (table:keys patch-ei-burnin-table-array)

  foreach ids
  [
    [x]->
    let xs table:get patch-pxcor-table-array x
    let ys table:get patch-pycor-table-array x
    table:put patch-pxcor-burnin-table-array x xs
    table:put patch-pycor-burnin-table-array x ys
  ]

;  show (word "ids first "id_ei_15%)
;  show (word "ei " patch-ei-table-array)
;  show (word "x " patch-pxcor-table-array)
;  show (word "y " patch-pycor-table-array)
;  show (word "ei15% " patch-ei-burnin-table-array)
;  show (word "x15% " patch-pxcor-burnin-table-array)
;  show (word "y15% " patch-pycor-burnin-table-array)

  ; seals should now choose a random patch from these 15% and set it as their next patch
  ; than clear all the memories patches from the burn in period and replace it with the newly chosen patch

  let patch_rand15 table:keys patch-ei-burnin-table-array
  let randSqID one-of patch_rand15
  let patch_im_on patch-here
  set my_next_patch one-of water-shore-patches with [km_5_id = randSqID and self != patch_im_on]

  let eis table:get patch-ei-burnin-table-array randSqID
  let xs table:get patch-pxcor-table-array randSqID
  let ys table:get patch-pycor-table-array randSqID

  table:clear patch-ei-table-array
  table:clear patch-mem-table-array
  table:clear patch-pxcor-table-array
  table:clear patch-pycor-table-array

  table:put patch-ei-table-array randSqID (list eis)
  table:put patch-mem-table-array randSqID 0.99
  table:put patch-pxcor-table-array randSqID xs
  table:put patch-pycor-table-array randSqID ys

end
;;;;;;;;;;;;;;;;;;;;;;;;;;;; ENERGETICS ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to energy_expenditure ; all in MJ per time step
  let seasonal_rmr_multiplier 2.08 ; Sparling 2006 and Chapter 4,8 (table 4.1) 1.95 0.26 -
  ;this is mean value over the entire year, here we introduce seasonal changes in RMR

  if activity = "HO" or activity = "SRS"
  [
     set ee (random-normal seasonal_rmr_multiplier 0.13) * BMR ; Sparling 2006 and Chapter 4,8 (table 4.1) 1.95 0.26
    ;- this is mean value over the entore year, here we introduce seasonal changes in RMR

    ; post moulting/breeding 2.08 0.13
    ; moulting 2.28 0.24
    ; pre-breeding 1.76 0.23
    ; breeding 1.86 0.22
  ]

  ; there is a marked difference between energy expenditure during haul-out and during digestion.
  ; If seals decide to digest while hauling out, their energy expendture must equal to ee of digestion
  ; not haul-out so the part of haulout needed for digestion hass ee = digestion,
  ; rest of the same haul-out ee=ee of haul out (seals which digest while hauling out usually
  ; spend longer time than it is neccessary just for digestion

  if activity = "LRS" or DurationOfDigestion > 0 ;
  [
    set ee (random-normal 3.5 1.2) * seasonal_rmr_multiplier * BMR; Sparling 2007 - digestion, has to be as multiplier of RMR
  ]


    if activity = "F" or activity = "LA" or activity = "TR-HO"
  [
    set ee (random-normal 1.8 0.45) * BMR ; Sparling 2004 says 1.7, but it was set in artificial consitions so it is likley to be higher in the field, I set it to 1.8 (Dave T pers. comm)
  ]

  set daily_ee (daily_ee + ee)
  set cum_ee (cum_ee + ee)
end

to energy_intake_poisson

  ; this is energy intake procedure based on Geert's example. The main reason for trying this approach is to
  ; reduce stochasticy and make the procedure more straightforward and widely applicable

  ; first a seal would 'calculate' number of fish encounter per bout (= one time step) and this is related to search rate and number of fish per patch, which is
  ; related to hsi of this patch (see landcape load)

  let patch_im_on patch-here
  ;if [categ] of patch_im_on = "land" [set patch_im_on min-one-of patches with [categ = "shore"] [distance myself]] ; if seals are close to shore, somwtimes patch-here is set to land
  if [categ] of patch_im_on = "land" [set patch_im_on min-one-of shore-patches [distance myself]] ; if seals are close to shore, somwtimes patch-here is set to land
  let #FishPatch [#Fish_m2] of patch_im_on
  ;show (word "per patch " #FishPatch " on p " patch_im_on)
  let #FishEncounteredPerBout (#FishPatch * search_rate );* (time_step / 60)) ; [#fish/m2] * [m2/min] * [s]
  ;show (word "Av " #FishEncounteredPerBout)
  ; the actual realized number of caught fish is subject of stochasticity and assumed to follow a poisson distribution

  ;let #FishCaught_temp random-poisson #FishEncounteredPerBout
  ifelse ((ResMass / Tmass) * 100) > MeanResMass%
  [
    ; this procedures is based on DEPONS with Cara's modifications and defines the reduction of food by fat animals
    let fatOver ((((ResMass / Tmass) * 100) - MeanResMass% )/(MaxResMass% - MeanResMass%))
    ;show fatOver
    let eiModifier exp (- 5 * fatOver)
    ;show eiModifier
    ;set #FishCaught round(#FishCaught_temp * eiModifier)
    set #FishEncounteredPerBout round(#FishEncounteredPerBout * eiModifier)
    ;show (word "Red from " #FishCaught_temp " to "#FishCaught)
  ]
  [
    ;set #FishCaught #FishCaught_temp
    set #FishEncounteredPerBout #FishEncounteredPerBout
  ]

  set #FishCaught random-poisson #FishEncounteredPerBout

  ; if burn-in scenario, there should not be habitat depletion for the duration of burn in
  ifelse (Large_scale_foraging = "BurnInTime" )
  [
    if ticks >= round(2881 * monthsOfBurnIn)
    [
      if #FishCaught > 0 and hab_depletion [habitat_depletion_and_HSIchange]
    ]

  ]
  [if #FishCaught > 0 and hab_depletion [habitat_depletion_and_HSIchange]]


  ; the translation from number of fish into gram and kJ of fish is based on two approaches
  ; once the number of fish is set, the gram and kJ of each = weighted mean of observed
  let eaten_grams (#FishCaught * mean_g_per_fish)
  set fishConsumed_g fishConsumed_g + (eaten_grams)
  ;show (word "Caught " #FishCaught " of " eaten_grams)
  set TotalfishConsumed_g TotalfishConsumed_g + (eaten_grams)
  set fishConsumed_g_longResting fishConsumed_g_longResting + (eaten_grams)
  set ei (eaten_grams * mean_kJ_per_gOffish * 0.0010) ; 0.0010 is conversion to MJ); * 0.80) ; 80% digestive efficiency, at the moment inlcuded in energy expenditure related to digestion
  set cum_ei cum_ei + ei
  set daily_ei daily_ei + ei

end


to net_energy

  ; calculates net energy at the beginning of each day (set as 8:00 in the morning) and turn into grams of blubber if net energy > 0
  ; for adult seals we assume that 100% OF ENERGY EXCESS IS TURNED INTO BLUBBER
  ; and currently for simplicity we assume that if net energy < 0, all energy comes from fat
  ; (CHANGE LATER IS WE FOLLOW BELTRAN ASSUMPTION THAT 53-60% COMES FROM FAT AND 40-47% FROM PROTEIN)
  ; efficiency of fat synthesis (deposition) is 0.74 - 0.90 (Malavear 2002 in Table 2 in Beltran) and energy content of 1g of fat is 0.0394 MJ (Blaxter 1989 in Table 2 in Beltran)
  ; efficiency of fat catabolism is 80% (Barboza 2008 in Beltran)

    set DailyNetEnergy (daily_ei - daily_ee) ; 0.80 digestion efficiency
    set daily_ee 0
    set daily_ei 0
    ifelse (Large_scale_foraging = "BurnInTime")
    [
      if (ticks >= round(2881 * monthsOfBurnIn))
      [
      ;show "I loose/gain weight"
        ifelse DailyNetEnergy > 0
        [
          let fatGain (DailyNetEnergy * (random-float (0.90 - 0.74) + 0.74)) / 0.0394
          set Tmass Tmass + fatGain
          set ResMass ResMass + fatGain
          set DaysWithNegativeDNE 0
          ifelse ((ResMass / Tmass) * 100) > MaxResMass% [set color orange] [ifelse Habitat = "Uniform" [set color black] [set color yellow]]
          ;stop
        ]

        ; FOR NOW WE ASSUME THAT 100% OF ENERGY IS TAKEN FROM FAT, MAY HAVE TO BE CHANGED, SEE WHAT I WROTE FEW LINES ABOVE
        [
          let fatLoss (DailyNetEnergy / (0.0394 * 0.80)) ; this will have 'minus' sign
          set Tmass Tmass + fatLoss
          set ResMass ResMass + fatLoss
          ifelse ((ResMass / Tmass) * 100) > MaxResMass% [set color orange] [set color blue]
          set DaysWithNegativeDNE DaysWithNegativeDNE  + 1
          ;set color blue
        ]

        ; as BMR is a function of Tmass is has to be recalculated
        ;set BMR ((70 * (Tmass / 1000)^ 0.75) * 0.004184 * (time_step / 60)) / (24 * 60) ; Tmass must be in kg and time in mins
        ;set BMR ((70 * (LBM / 1000)^ 0.75) * 0.004184 * (time_step / 60)) / (24 * 60)
        switch_to_large_scale_movement
      ]
    ]
    [
    ifelse DailyNetEnergy > 0
    [
      let fatGain (DailyNetEnergy * (random-float (0.90 - 0.74) + 0.74)) / 0.0394
      set Tmass Tmass + fatGain
      set ResMass ResMass + fatGain
      ;ifelse ((ResMass / Tmass) * 100) > MaxResMass% [set color orange] [ifelse Habitat = "Uniform" [set color black] [set color yellow]]
      set DaysWithNegativeDNE 0
      ;set color yellow
      ;show (word "Im eating yellow ok, my target patch is " my_next_patch " my bad is " my_move_away_patch)
      ;stop
    ]

    ; FOR NOW WE ASSUME THAT 100% OF ENERGY IS TAKEN FROM FAT, MAY HAVE TO BE CHANGED, SEE WHAT I WROTE FEW LINES ABOVE
    [
      let fatLoss (DailyNetEnergy / (0.0394 * 0.80)) ; this will have 'minus' sign
      set Tmass Tmass + fatLoss
      set ResMass ResMass + fatLoss
      set DaysWithNegativeDNE DaysWithNegativeDNE  + 1
      ifelse ((ResMass / Tmass) * 100) > MaxResMass% [set color orange] [set color blue]
      ;show (word "Im eating blue ok, my target patch is " my_next_patch " my bad is " my_move_away_patch)

    ]
      switch_to_large_scale_movement
      ; as BMR is a function of Tmass is has to be recalculated
      ;set BMR ((70 * (Tmass / 1000)^ 0.75) * 0.004184 * (time_step / 60)) / (24 * 60) ; Tmass must be in kg and time in mins
      ;set BMR ((70 * (LBM / 1000)^ 0.75) * 0.004184 * (time_step / 60)) / (24 * 60)

      ;show "Co cholera"
    ;]
    ]


end

; here is procedure what happens of seals are loosing energy for X number of consequtive days
to switch_to_large_scale_movement

  ; the first option results in seals 'clearing out' their patch memory if they keep loosing energy for NdaysWithNegativeNEB, so they basically go back to CRW
  if (Large_scale_foraging = "Switch2CRW")
  [
  ifelse DaysWithNegativeDNE >= NdaysWithNegativeNEB
    [
      set my_next_patch  "none"
      set patch-ei-table-array table:make
      set patch-mem-table-array table:make
      set patch-pxcor-table-array table:make
      set patch-pycor-table-array table:make
      set color black
          ]
    [
      stop
    ]
  ]

  ; if seals don't get enough fish for 7 consequtive days, they choose a random patch from their 15% best 'collected' during burn in, set ot as their target patch and than clear the
  ; memory of the remaining patches. 15% list stays intacked
  if (Large_scale_foraging = "BurnInTime" )
  [
    ifelse DaysWithNegativeDNE >= NdaysWithNegativeNEB
    [
      let patch_rand15 table:keys patch-ei-burnin-table-array
      let randSqID one-of patch_rand15
      let patch_im_on patch-here
      set my_next_patch one-of water-shore-patches with [km_5_id = randSqID and self != patch_im_on]

      let eis table:get patch-ei-burnin-table-array randSqID
      let xs [pxcor] of my_next_patch
      let ys [pycor] of my_next_patch

      table:clear patch-ei-table-array
      table:clear patch-mem-table-array
      table:clear patch-pxcor-table-array
      table:clear patch-pycor-table-array

      table:put patch-ei-table-array randSqID (list eis)
      table:put patch-mem-table-array randSqID 0.99
      table:put patch-pxcor-table-array randSqID xs
      table:put patch-pycor-table-array randSqID ys
      ;show "I switched to latrge scale movement"
    ]
    [
      stop
    ]
  ]

 if (Large_scale_foraging = "MoveAway")
  [
  ifelse DaysWithNegativeDNE >= NdaysWithNegativeNEB
    [
      if DaysWithNegativeDNE = NdaysWithNegativeNEB [set my_move_away_patch my_next_patch]; show (word "I just started LSM, my away patch is " my_move_away_patch " my old target patch is "my_next_patch)]
      set my_next_patch  "none"
      set patch-ei-table-array table:make
      set patch-mem-table-array table:make
      set patch-pxcor-table-array table:make
      set patch-pycor-table-array table:make
      set color black
      ;show "Im moving away"
    ]
    [
      set my_move_away_patch "none"
      ;show (word "Im eating blue ok, my target patch is " my_next_patch " my bad is " my_move_away_patch)
      ;stop
    ]
  ]


  ; if seals don't get enough fish for 7 consequtive days, they choose a random patch from their best patches based on 90% of best patches in the area, set ot as their target patch and than clear the
  ; memory of the remaining patches. 90% list stays intacked
  if (Large_scale_foraging = "Omniscience" )
  [
    ifelse DaysWithNegativeDNE >= NdaysWithNegativeNEB
    [
      set color black
      let patch_rand table:keys patches_km25_maxhsi
      let randSqID25 one-of patch_rand
      let patch_im_on patch-here
      set my_next_patch one-of water-shore-patches with [km_25_id = randSqID25 and self != patch_im_on]
      ;show (word "ID is " randSqID25 " Patch test is " my_next_patch " target is " target)
      ; seals now calculate potential ei on this chosen patch, as they do for energy intake
      let this_maxhsi  table:get patches_km25_maxhsi randSqID25
      let #Fish_m2_maxHSI this_maxhsi * #Fish_hsi_multiplier
      let #FishEncounteredPerBout (#Fish_m2_maxHSI * search_rate )
      let #FishCaughtMaxHSI random-poisson #FishEncounteredPerBout
      ; the translation from number of fish into gram and kJ of fish is based on two approaches
      ; once the number of fish is set, the gram and kJ of each = weighted mean of observed
      let eaten_grams (#FishCaughtMaxHSI * mean_g_per_fish)
      let eiMaxHSI (eaten_grams * mean_kJ_per_gOffish * 0.0010) ; 0.0010 is conversion to MJ); * 0.80) ; 80% digestive efficiency, at the moment inlcuded in energy expenditure related to digestion
      ;show (word "jem " eiMaxHSI)
      let xs [pxcor] of my_next_patch
      let ys [pycor] of my_next_patch
      let id5 [km_5_id] of my_next_patch

      table:clear patch-ei-table-array
      table:clear patch-mem-table-array
      table:clear patch-pxcor-table-array
      table:clear patch-pycor-table-array

      table:put patch-ei-table-array id5 (list eiMaxHSI)
      table:put patch-mem-table-array id5 0.99
      table:put patch-pxcor-table-array id5 xs
      table:put patch-pycor-table-array id5 ys

;      show patch-ei-table-array
;      show patch-mem-table-array
;      show patch-pxcor-table-array
;      show patch-pycor-table-array
      ;show "I switched to large scale movement"
    ]
    [
      stop
    ]
  ]

  if (Large_scale_foraging = "OmniSwitch" )
  [
    let whichLSF random 1.1
    ifelse whichLSF = 0 ; if zero than as Omniscience, if 1 as Switch2crw
    [
      ifelse DaysWithNegativeDNE >= NdaysWithNegativeNEB
      [
        set color black
        let patch_rand table:keys patches_km25_maxhsi
        let randSqID25 one-of patch_rand
        let patch_im_on patch-here
        set my_next_patch one-of water-shore-patches with [km_25_id = randSqID25 and self != patch_im_on]
        ;show (word "ID is " randSqID25 " Patch test is " my_next_patch " target is " target)
        ; seals now calculate potential ei on this chosen patch, as they do for energy intake
        let this_maxhsi  table:get patches_km25_maxhsi randSqID25
        let #Fish_m2_maxHSI this_maxhsi * #Fish_hsi_multiplier
        let #FishEncounteredPerBout (#Fish_m2_maxHSI * search_rate )
        let #FishCaughtMaxHSI random-poisson #FishEncounteredPerBout
        ; the translation from number of fish into gram and kJ of fish is based on two approaches
        ; once the number of fish is set, the gram and kJ of each = weighted mean of observed
        let eaten_grams (#FishCaughtMaxHSI * mean_g_per_fish)
        let eiMaxHSI (eaten_grams * mean_kJ_per_gOffish * 0.0010) ; 0.0010 is conversion to MJ); * 0.80) ; 80% digestive efficiency, at the moment inlcuded in energy expenditure related to digestion
                                                                  ;show (word "jem " eiMaxHSI)
        let xs [pxcor] of my_next_patch
        let ys [pycor] of my_next_patch
        let id5 [km_5_id] of my_next_patch

        table:clear patch-ei-table-array
        table:clear patch-mem-table-array
        table:clear patch-pxcor-table-array
        table:clear patch-pycor-table-array

        table:put patch-ei-table-array id5 (list eiMaxHSI)
        table:put patch-mem-table-array id5 0.99
        table:put patch-pxcor-table-array id5 xs
        table:put patch-pycor-table-array id5 ys

        ;      show patch-ei-table-array
        ;      show patch-mem-table-array
        ;      show patch-pxcor-table-array
        ;      show patch-pycor-table-array
        ;show "I switched to large scale movement"
      ]
      [
        stop
      ]
    ]
    [
      ifelse DaysWithNegativeDNE >= NdaysWithNegativeNEB
      [
        set my_next_patch  "none"
        set patch-ei-table-array table:make
        set patch-mem-table-array table:make
        set patch-pxcor-table-array table:make
        set patch-pycor-table-array table:make
        set color black
      ]
      [
        stop
      ]
    ]
  ]

end

to habitat_depletion_and_HSIchange

  ;show "habitat depletion"
  let patch_im_on patch-here
  ;if [categ] of patch_im_on = "land" [set patch_im_on min-one-of patches with [categ = "shore"] [distance myself]] ; if seals are close to shore, somwtimes patch-here is set to land
  if [categ] of patch_im_on = "land" [set patch_im_on min-one-of shore-patches [distance myself]] ; if seals are close to shore, somwtimes patch-here is set to land

  ;show [#FishTotal] of patch-here
  let #FishTotal_temp ([#FishTotal] of patch_im_on - #FishCaught)
  ask patch_im_on
  [
  set #FishTotal #FishTotal_temp
  if #FishTotal < 0 [set #FishTotal 0]
  set #Fish_m2 (#FishTotal / patch_size ^ 2)
  set hsi (#Fish_m2 / #Fish_hsi_multiplier)
  ;show #FishTotal
  ]

end



;;;;;;;;;;;;;;;;;;; DEBUGGING ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; 1. Debugging land avoidance procedure - If there is land ahead in step length or smaller, seal calculates the amount of land
; on both sides and chooses the direction where there is less land. Count-left and count-right is outputted abd checked whether seals choose direction towards less land

to debugging_land_avoidance
  file-open "Debugging/land_avoidance_RandL_count.csv"
;  ask seals
;  [
;  file-type ticks                  file-type ","
;  file-type is-land?               file-type ","
;  file-type count-right            file-type ","
;  file-type count-left             file-type ","
;  file-type avoidance-mode-right   file-type ","
;  file-print avoidance-mode-left
;  ]

  ; same code as above just less to write (I think so)
    ask seals
  [
  file-type (word ticks ",")
  file-type (word avoid? ",")
  file-type (word count-right ",")
  file-type (word count-left ",")
  file-type (word avoidance-mode-right ",")
  file-print (word avoidance-mode-left ",")
  ]
  file-close
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; OUTPUT FILES ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to output_file_csv_extension


  file-open outputFileName_movement
  ask seals
    [

      file-print csv:to-row (list
        expindex
        ticks
        who
        dayNumber
        current-time
        activity
        fishConsumed_g
        TotalfishConsumed_g
        daily_ee
        daily_ei
        dailyNetEnergy
        DaysWithNegativeDNE
        typeOfHa
        fishConsumed_g_longResting
        distHaulOutDigestion
        Blength
        Tmass
        LBM
        ResMass
        sex
        xcor
        ycor
        [km_25_id] of patch-here
        [km_5_id] of patch-here
        foraging_trip#
        durationOfresting
        previous_ho
      )

    ]
  file-close




end

to output_spatial
  file-open outputFileName_spatial
  ask patches
  [
    file-print csv:to-row (list
      expindex
      ticks
      pxcor
      pycor
      Nseals
    )
  ]
  file-close


  file-open outputFileName_ho_list
  ask seals
  [
    file-print csv:to-row (list
      expindex
      ticks
      who
      visited_ho_list
      )
  ]
  file-close

    file-open outputFileName_depletionEnd
  ask patches
  [
    file-print csv:to-row (list
      expindex
      ticks
      pxcor
      pycor
      hsi
      #fish_m2
    )
  ]
  file-close

end

 to output-tripExtentDuration
    file-open outputFileName_tripDurationExtend
  ask seals
  [
    file-print csv:to-row (list
      expindex
     ticks
      who
      dayNumber
      current-time
      xcor
      ycor
      foraging_trip#
      durationOfresting
      activity
      previous_ho
    )
  ]
  file-close
end

to output_crw

  file-open outputFileName_crw
  ask seals
    [

      file-print csv:to-row (list
        ticks
        who
        xcor
        ycor
        b
        sigK
        imp_dist
        turn_crw
        turn_hsi
        turn_bias
        dist2target
        prev_angle
        hsi_Im_on
        heading
        diff_heading
        speed
      )
    ]

  file-close
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; PARAMETERISATION ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;to parameterisation_energyIntake_geert
;  file-open "Input/Param-energIntakeGeert.txt"      ;; Import matrix for parameterisation
;  let param-matrix file-read
;  set param-data matrix:from-row-list param-matrix
;  file-close
;
;  set #Fish_hsi_multiplier (matrix:get param-data expindex 0)
;  set search_rate          (matrix:get param-data expindex 1)
;
;end

;to testing_b_prob
;  file-open "Input/Param_newb_prob2.txt"      ;; Import matrix for parameterisation
;  let param-matrix file-read
;  set param-data matrix:from-row-list param-matrix
;  file-close
;
;  ;set Large_scale_foraging (matrix:get param-data expindex 0)
;  set hab_depletion          (matrix:get param-data expindex 0)
;  set b_prob (matrix:get param-data expindex 1)
;
;end
;
;to parameterisation_memories
;  file-open "Input/Param-memories.txt"      ;; Import matrix for parameterisation
;  let param-matrix file-read
;  set param-data matrix:from-row-list param-matrix
;  file-close
;
;  set ref-mem-decay-rate            (matrix:get param-data expindex 0)
;  set haulOut_detection_distance    (matrix:get param-data expindex 1)
;  set mem_level_passedBy_ho         (matrix:get param-data expindex 2)
;
;end
;
;to parameterisation_global
;  ;file-open "Input/Param-global.txt"      ;; Import matrix for parameterisation
;  ;file-open "Input/Param-global-goodDur.txt"      ;; Import matrix for parameterisation
;  file-open "Input/Parameteryzcja#2.txt"      ;; Import matrix for parameterisation
;  let param-matrix file-read
;  set param-data matrix:from-row-list param-matrix
;  file-close
;
;  set ref-mem-decay-rate            (matrix:get param-data expindex 0)
;  set haulOut_detection_distance    (matrix:get param-data expindex 1)
;  set mem_level_passedBy_ho         (matrix:get param-data expindex 2)
;  set mean_durHO                    (matrix:get param-data expindex 3)
;  set b_prob                        (matrix:get param-data expindex 4)
;  ; tylko w Parameteryzcja#2.txt
;  set b_prob2                       (matrix:get param-data expindex 5)
;  set multProbHaul                  (matrix:get param-data expindex 6)
;
;end

;to parameterisation_global_6param
;  file-open "Input/Param-global-6param_150comb.txt"      ;; Import matrix for parameterisation
;  let param-matrix file-read
;  set param-data matrix:from-row-list param-matrix
;  file-close
;
;  set #Fish_hsi_multiplier          (matrix:get param-data expindex 0)
;  set search_rate                   (matrix:get param-data expindex 1)
;  set ref-mem-decay-rate            (matrix:get param-data expindex 2)
;  set haulOut_detection_distance    (matrix:get param-data expindex 3)
;  set mean_durHO                    (matrix:get param-data expindex 4)
;  set b_prob                        (matrix:get param-data expindex 5)
;  ; tylko w Parameteryzcja#2.txt
;;  set b_prob2                       (matrix:get param-data expindex 5)
;;  set multProbHaul                  (matrix:get param-data expindex 6)
;
;end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;; PLOTTING ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to plotting

  ; First update the line graph of mean similarity
  set-current-plot "Mean blubber %"
  plot (mean [ResMass] of seals / mean [Tmass] of seals) * 100

;  set-current-plot "Energy expenditure [MJ]"
;  if time:get "hour" current-time = 7 and time:get "minute" current-time = 45
;  [plot mean [daily_ee] of seals]

  set-current-plot "Days with negative energy"
  histogram [DaysWithNegativeDNE] of seals
;  ; Histogram of the similarity of each patch with each of its neighbors
;  ; This requires a list of all the similarity values
;  ; ("histogram [([similarity-with myself] of neighbors4] of patches" does not
;  ; work because it creates a list of lists instead of a big list of similarity values.)
;  let a-similarity-list (list)
;  ask patches
;  [
;    set a-similarity-list sentence a-similarity-list ([similarity-with myself] of neighbors4)
;  ]
;  set-current-plot "Similarity histogram"
;  histogram a-similarity-list

end


;;;;;;;;;;;;;;;;;;; REPORTERS (equivalents of functions in R) ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


to-report heading-to-angle [ h ]
  report (90 - h) mod 360
end

to-report arctan [x]
  report asin (x / sqrt(1 + x * x))
end

to-report  log-normal [mu sigma]

; to use it I just have to write log-normal 1 2 and it will show log-normal value with mu = 1 and sigma = 2
  let beta ln (1 + ((sigma ^ 2) / (mu ^ 2)))


  let x exp (random-normal (ln (mu) - (beta / 2)) sqrt beta)


  report x


end

to-report cumsum [ x ]
  let ll x
  let tempsum 0
  let tempsum_list []
  foreach ll
  [
    a ->
    set tempsum (tempsum + a )
    set tempsum_list lput tempsum tempsum_list
  ]

  report tempsum_list
end

to-report FindFish  [ r ] ; r -random number between 0 and 1

  ;foreach [CumSum100] [a -> show (r > a)]
  ;report Fish
  let perc (list 3.8 8.3 18.4 44.6 3.6 8.4 4.4)
  let csum cumsum (perc)
  let csum100 map [ a -> a / 100] csum
  ;let r random-float 1
  let tr []
  ;let tr map
  foreach csum100
  [c -> ;c > r; show (word d " -> " round d)
    ifelse c < r
    [set tr lput 1 tr]
    [set tr lput 0 tr]
  ]
  ;csum100
;  ;let tr map [z -> r > z] csum100
  let sumtr sum tr
  report sumtr

end


;The binomial reporter essentially runs a bunch of Bernoulli trials (is a U(0,1) number less than
; specified probability of success 'p') and then sums up the cases where this
;is true, which is  Binomial distribution (defined by  choice of 'n'
;trials each with probability of success 'p').

to-report binomial [n p]
    let bin_ct 0
    repeat n [
        if random-float 1 < p [set bin_ct bin_ct + 1]
    ]
    report bin_ct
end


to-report sort-nested-lists-from-highest [l]
  report sort-by [[list1 list2] -> first list1 > first list2] l
end

to-report sort-nested-lists-from-highest-idFirst [l]
  report sort-by [[list2 list1] -> first list1 > first list2] l
end

to-report sort-nested-lists-from-lowest [l]
  report sort-by [[list1 list2] -> first list1 < first list2] l
end

to-report numbers_in_list_below_sth [threshold lista]

  let indices n-values (length lista) [i -> i]
  ;show indices
  report filter [x -> item x lista <= threshold] indices
end
@#$#@#$#@
GRAPHICS-WINDOW
224
10
998
518
-1
-1
2.7
1
10
1
1
1
0
0
0
1
0
283
0
184
1
1
1
ticks
20.0

BUTTON
8
10
63
43
NIL
setups
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SLIDER
8
219
180
252
n_of_seals
n_of_seals
0
3683
22.0
1
1
NIL
HORIZONTAL

SLIDER
823
572
995
605
avg_Vsus
avg_Vsus
0
2
0.29
0.01
1
m/s
HORIZONTAL

SLIDER
823
608
995
641
std_Vsus
std_Vsus
0
2
0.33
0.01
1
m/s
HORIZONTAL

SLIDER
823
535
995
568
std_dir
std_dir
0
180
58.0
1
1
NIL
HORIZONTAL

SWITCH
104
50
199
83
show_pd?
show_pd?
0
1
-1000

SLIDER
824
719
996
752
time_step
time_step
60
1200
900.0
60
1
s
HORIZONTAL

MONITOR
1013
489
1118
534
Current-time
time:show current-time \"yyyy/MM/dd HH:mm\"
17
1
11

MONITOR
1135
489
1212
534
NIL
dayNumber
17
1
11

PLOT
1010
10
1210
160
Mean blubber %
Ticks
%
0.0
4000.0
0.0
50.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" ""

PLOT
1011
171
1211
321
Energy expenditure [MJ]
ticks
Daily EE
0.0
40.0
0.0
20.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" ""

BUTTON
67
10
130
43
NIL
go
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
134
10
197
43
NIL
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SWITCH
7
50
97
83
output?
output?
0
1
-1000

INPUTBOX
520
655
584
715
sigK
0.52
1
0
Number

INPUTBOX
592
655
655
715
b
-0.7
1
0
Number

INPUTBOX
345
523
456
583
#Fish_hsi_multiplier
3.0
1
0
Number

INPUTBOX
460
524
533
584
search_rate
13.176
1
0
Number

SWITCH
8
259
160
292
hab_depletion
hab_depletion
0
1
-1000

INPUTBOX
665
590
818
650
end-time
2018/07/31 00:00
1
0
String

INPUTBOX
446
656
511
716
imp_dist
0.015
1
0
Number

INPUTBOX
539
524
637
584
mean_g_per_fish
23.1
1
0
Number

INPUTBOX
226
525
333
585
mean_kJ_per_gOffish
5.4
1
0
Number

SLIDER
823
643
995
676
scale_alpha_k
scale_alpha_k
0
1
0.6
0.01
1
NIL
HORIZONTAL

SLIDER
824
683
996
716
lambda_rate
lambda_rate
0
3
1.88
0.01
1
NIL
HORIZONTAL

TEXTBOX
830
522
980
540
CRW parameters
11
0.0
1

INPUTBOX
225
591
380
651
mem_level_passedBy_ho
0.5
1
0
Number

INPUTBOX
388
590
498
650
ref-mem-decay-rate
0.009
1
0
Number

INPUTBOX
506
590
661
650
haulOut_detection_distance
2.07
1
0
Number

INPUTBOX
740
526
816
586
b_prob
-0.343
1
0
Number

INPUTBOX
645
526
733
586
mean_durHO
7.372
1
0
Number

CHOOSER
9
296
152
341
Large_scale_foraging
Large_scale_foraging
"Switch2CRW" "BurnInTime" "MoveAway" "Omniscience" "Basic" "OmniSwitch"
5

PLOT
1012
331
1212
481
Days with negative energy
NIL
NIL
0.0
20.0
0.0
30.0
true
false
"" ""
PENS
"default" 1.0 1 -16777216 true "" ""

CHOOSER
10
351
153
396
Patch_memory
Patch_memory
"on" "off"
0

INPUTBOX
9
653
135
713
NdaysWithNegativeNEB
7.0
1
0
Number

INPUTBOX
9
400
124
460
monthsOfSimulation
3.0
1
0
Number

MONITOR
1015
541
1072
586
#seals
count seals
17
1
11

INPUTBOX
9
464
105
524
monthsOfBurnIn
0.0
1
0
Number

INPUTBOX
9
527
164
587
percOfsavedPatchesFromBurnIn
15.0
1
0
Number

INPUTBOX
9
590
164
650
percOfsavedPatchesOmniscience
90.0
1
0
Number

INPUTBOX
225
657
327
717
MeanResMass%
22.0
1
0
Number

INPUTBOX
335
657
437
717
MaxResMass%
45.0
1
0
Number

CHOOSER
5
90
145
135
Habitat
Habitat
"GrecianEtAl2018"
0

@#$#@#$#@
## WHAT IS IT?

(a general understanding of what the model is trying to show or explain)

## HOW IT WORKS

(what rules the agents use to create the overall behavior of the model)

## HOW TO USE IT

(how to use the model, including a description of each of the items in the Interface tab)

## THINGS TO NOTICE

(suggested things for the user to notice while running the model)

## THINGS TO TRY

(suggested things for the user to try to do (move sliders, switches, etc.) with the model)

## EXTENDING THE MODEL

(suggested things to add or change in the Code tab to make the model more complicated, detailed, accurate, etc.)

## NETLOGO FEATURES

(interesting or unusual features of NetLogo that the model uses, particularly in the Code tab; or where workarounds were needed for missing features)

## RELATED MODELS

(models in the NetLogo Models Library and elsewhere which are of related interest)

## CREDITS AND REFERENCES

(a reference to the model's URL on the web if it has one, as well as any other necessary credits, citations, and links)
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

seal
true
0
Polygon -7500403 true true 30 120 30 120 150 120 195 90 225 75 255 75 270 90 285 90 285 105 240 120 210 165 195 180 180 195 135 195 90 180 30 165 30 120 30 120
Polygon -7500403 true true 30 120 21 87 15 86 0 120 15 150 0 180 13 214 20 212 30 165
Polygon -1 true false 150 150 195 165 150 180 121 180 91 174 105 135
Circle -16777216 true false 240 86 10

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.2.0
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="Param-EnergyIntake" repetitions="1" runMetricsEveryStep="true">
    <setup>setups</setup>
    <go>go</go>
    <timeLimit steps="2978"/>
    <metric>ticks</metric>
    <steppedValueSet variable="expindex" first="11" step="1" last="11"/>
  </experiment>
  <experiment name="Param-mamories" repetitions="1" runMetricsEveryStep="true">
    <setup>setups</setup>
    <go>go</go>
    <timeLimit steps="2978"/>
    <metric>ticks</metric>
    <steppedValueSet variable="expindex" first="0" step="1" last="64"/>
  </experiment>
  <experiment name="Param-global" repetitions="1" runMetricsEveryStep="true">
    <setup>setups</setup>
    <go>go</go>
    <metric>ticks</metric>
    <steppedValueSet variable="expindex" first="95" step="1" last="99"/>
  </experiment>
  <experiment name="experiment" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="haulOut_detection_distance">
      <value value="3.9827731230000003"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="search_rate">
      <value value="0.017"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="show_pd?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean_durHO">
      <value value="2.090821431"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="avg_Vsus">
      <value value="0.29"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="#Fish_hsi_multiplier">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="time_step">
      <value value="900"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="hab_depletion">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ref-mem-decay-rate">
      <value value="0.004351249"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="b">
      <value value="-0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="end-time">
      <value value="&quot;2018/10/31 00:00&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mem_level_passedBy_ho">
      <value value="0.822854609"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="b_prob">
      <value value="-0.243764576"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std_dir">
      <value value="60"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="expindex">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n_of_seals">
      <value value="350"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std_Vsus">
      <value value="0.33"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sigK">
      <value value="0.54"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="Param-global-test" repetitions="1" runMetricsEveryStep="true">
    <setup>setups</setup>
    <go>go</go>
    <metric>ticks</metric>
    <steppedValueSet variable="expindex" first="47" step="1" last="55"/>
  </experiment>
  <experiment name="experiment" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="haulOut_detection_distance">
      <value value="2.235341894"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="search_rate">
      <value value="0.017"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="show_pd?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean_durHO">
      <value value="9.994918348"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="avg_Vsus">
      <value value="0.29"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="#Fish_hsi_multiplier">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="time_step">
      <value value="900"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="b">
      <value value="-0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ref-mem-decay-rate">
      <value value="0.007362162"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="hab_depletion">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="end-time">
      <value value="&quot;2018/10/31 00:00&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mem_level_passedBy_ho">
      <value value="0.502516316"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="b_prob">
      <value value="-0.234681495"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std_dir">
      <value value="60"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="expindex">
      <value value="74"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n_of_seals">
      <value value="350"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std_Vsus">
      <value value="0.33"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sigK">
      <value value="0.54"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="Param2-distVspropbOfHaulOut" repetitions="1" runMetricsEveryStep="true">
    <setup>setups</setup>
    <go>go</go>
    <metric>ticks</metric>
    <steppedValueSet variable="expindex" first="0" step="1" last="3"/>
  </experiment>
  <experiment name="Final-index43_many replicates" repetitions="5" runMetricsEveryStep="true">
    <setup>setups</setup>
    <go>go</go>
    <timeLimit steps="100"/>
    <metric>ticks</metric>
    <steppedValueSet variable="expindex" first="0" step="1" last="0"/>
  </experiment>
  <experiment name="Param-global-6param" repetitions="1" runMetricsEveryStep="true">
    <setup>setups</setup>
    <go>go</go>
    <metric>ticks</metric>
    <steppedValueSet variable="expindex" first="62" step="1" last="62"/>
  </experiment>
  <experiment name="Final_index127" repetitions="24" runMetricsEveryStep="true">
    <setup>setups</setup>
    <go>go</go>
    <metric>ticks</metric>
  </experiment>
  <experiment name="b_prob" repetitions="2" runMetricsEveryStep="true">
    <setup>setups</setup>
    <go>go</go>
    <metric>ticks</metric>
    <enumeratedValueSet variable="hab_depletion">
      <value value="false"/>
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="b_prob">
      <value value="-0.16"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Large_scale_foraging">
      <value value="&quot;Basic&quot;"/>
      <value value="&quot;Omniscience&quot;"/>
      <value value="&quot;BurnInTime&quot;"/>
      <value value="&quot;Switch2crw&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="Repeated runs" repetitions="4" runMetricsEveryStep="true">
    <setup>setups</setup>
    <go>go</go>
    <timeLimit steps="8644"/>
    <metric>ticks</metric>
    <enumeratedValueSet variable="lambda_rate">
      <value value="1.88"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="haulOut_detection_distance">
      <value value="2.07"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="percOfsavedPatchesFromBurnIn">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="search_rate">
      <value value="13.176"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="monthsOfSimulation">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="#Fish_hsi_multiplier">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="MeanResMass%">
      <value value="28"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dodatek">
      <value value="&quot;corrNeed2HO&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="time_step">
      <value value="900"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Habitat">
      <value value="&quot;James&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ref-mem-decay-rate">
      <value value="0.009"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="b">
      <value value="-0.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="hab_depletion">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean_g_per_fish">
      <value value="48"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale_alpha_k">
      <value value="0.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="NdaysWithNegativeNEB">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="b_prob">
      <value value="-0.343"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n_of_seals">
      <value value="350"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Patch_memory">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Large_scale_foraging">
      <value value="&quot;OmniSwitch&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="show_pd?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean_kJ_per_gOffish">
      <value value="4.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="avg_Vsus">
      <value value="0.29"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean_durHO">
      <value value="7.372"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="MaxResMass%">
      <value value="45"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="end-time">
      <value value="&quot;2018/10/31 00:00&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mem_level_passedBy_ho">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="imp_dist">
      <value value="0.015"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std_dir">
      <value value="58"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="percOfsavedPatchesOmniscience">
      <value value="90"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std_Vsus">
      <value value="0.33"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="monthsOfBurnIn">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sigK">
      <value value="0.52"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="4" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="lambda_rate">
      <value value="1.88"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="haulOut_detection_distance">
      <value value="2.07"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="percOfsavedPatchesFromBurnIn">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="search_rate">
      <value value="13.176"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="monthsOfSimulation">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="#Fish_hsi_multiplier">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dodatek">
      <value value="&quot;corrNeed2HO&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="MeanResMass%">
      <value value="28"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="time_step">
      <value value="900"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="hab_depletion">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ref-mem-decay-rate">
      <value value="0.009"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Habitat">
      <value value="&quot;James&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="b">
      <value value="-0.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean_g_per_fish">
      <value value="48"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale_alpha_k">
      <value value="0.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="NdaysWithNegativeNEB">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="b_prob">
      <value value="-0.343"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n_of_seals">
      <value value="350"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Patch_memory">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Large_scale_foraging">
      <value value="&quot;OmniSwitch&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="show_pd?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean_kJ_per_gOffish">
      <value value="4.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean_durHO">
      <value value="7.372"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="MaxResMass%">
      <value value="45"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="avg_Vsus">
      <value value="0.29"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="end-time">
      <value value="&quot;2018/10/31 00:00&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mem_level_passedBy_ho">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="imp_dist">
      <value value="0.015"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="percOfsavedPatchesOmniscience">
      <value value="90"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std_dir">
      <value value="58"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="std_Vsus">
      <value value="0.33"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="monthsOfBurnIn">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sigK">
      <value value="0.52"/>
    </enumeratedValueSet>
  </experiment>
</experiments>
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@

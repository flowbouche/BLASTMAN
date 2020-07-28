; Developed by: Adam J Bouché for a Master of Science Thesis in Forest Ecosystems and Society at Oregon State University
; 2017 - 2020

;; _________________________TO DO ASAP_______________________________
;  - Implement counts of trees infected by root, sc, pf, hn (current and cumulative)
;  - Add patch properties to "patches-own"?
;  - make sure all "set-" and "reset-" procedures for patches and cells include all variables
;  - add to data: export stand data - worth it? Implemented already?
;  - check attraction thresholds for infected trees
;  - why are there errors in n and % trees killed by different factors (e.g., % > 100)?
;  - add attraction for recently thinned and harvested stands at initialization
;  - fix root radius use
;  - Can the model work with relative paths?

;; _________________________TO ADD _______________________________
;  - Greater flexibility in silvicultural treatments (set densities and auto-calculate % removed for thinning) (maybe a third management approach "custom")

;__________________________________________________________________________________________________________________________________________________
;__________________________________________________________________________________________________________________________________________________
; MODEL INITIALIZATION
;__________________________________________________________________________________________________________________________________________________
; Included files: ring globals defined in "ring_globals_330.nls" (to keep code clean)
; file to set ring globals with coordenates is too large to open the model
__includes [
  "supporting_files/basic_functions_1.10.0.0.nls"
  "supporting_files/insect-attraction_1.10.0.0.nls"
  "supporting_files/probability_parameters_1.10.0.0.nls"
  "supporting_files/ring_lists_1.10.0.0.nls"
  "supporting_files/spread-infection_Go_1.10.0.0.nls"
  ;"supporting_files/under_development_1.9.2.0.nls"
  ;"supporting_files/WIP_attraction_function_dvlpment_1.9.2.0.nls"  ; TEMP: Development of attraction function
]

; Import extensions
extensions [
  gis
  r
  ;profiler
  ;palette
]

; Define global variables
globals [
  ; GENERAL, GIS & files
  base-setup?              ; Whether base has been set up
  run-id-no                ; unique run id number
  run-id                   ; unique run id
  model-dir                ; model directory
  ;run-dir                 ; run directory ; check delete?
  run-file                 ; path to run file
  progress-file            ; file to keep track of model progress
  stand-data-file          ; file for stand-level data output
  tree-data-file           ; file for data output
  cell-coords-file         ; file will cell distances
  test-timer               ; Useful for a timer that needs to be used in a global context
  current-date-and-time    ; Used to store an abbreviated timestamp of current date and time
  intercell-distance       ; Varb to set distance between cells/patch dimensions (m)
  cell-area                ; The area represented by each patch/cell (m^2)

  ; GIS, OUTPUTS
  projection               ; projection used for output files

  ; LANDSCAPE GENERATOR / MANAGEMENT
  seed-radius               ; used to make the stands more evenly sized
  ; Stand lists: (lists of stand-id numbers, NOT agentsets of cells with those mgmt types)
  stand-list                ; a list of all stand-id numbers (1 to total number of stands)
  short-rotation-stand-list      ; list of stand-id #s of 'short-rotation' stands
  long-rotation-stand-list      ; list of stand-id #s of 'long-rotation stands'
  set-aside-stand-list      ; list of stand-id #s of 'set-aside stands'
  regen-list-short-rotation      ; list of short-rotation stands to be regenerated
  regen-list-long-rotation      ; list of long-rotation stands to be regenerated
  thin-list-PCT             ; list of thinned stands (for attraction)
  thin-list-CT1             ; list of thinned stands (for attraction)
  thin-list-CT2             ; list of thinned stands (for attraction)
  plantation-road-density   ; density of roads in plantations in meters per hectare (34.2 m/ha, from John Sessions, personal communication)
  ; set-aside-road-density    ; ADD
  thinning-occurred?        ; bool to track whether a thinning occurred
  harvest-occurred?         ; bool to track whether a harvest occurred

  ; GLOBAL AGENTSETS
  trees_current              ; an agentset of all current trees
  trees_potential            ; an agentset of all cells than will be trees during a model run
  edges                      ; an agentset of all cells/patches at a stand edge
  trees-to-thin              ; agentset: all trees to be thinned
  trees-to-harvest           ; agentset: all trees to be harvested

  ; INFECTION MODEL          ; CHECK & delete those not used (if any)
  trees_newly-infected       ; an agentset of all trees that are newly infected ; check - in reset-globals
  trees_losing-infection     ; an agentset of all trees that have lost their infection (due to mortality and inoculum viability decay ; add to globals to reset
  trees_infected             ; an agentset of all infected trees
  cells-projecting           ; an agentset of all cells projecting infection or attraction
  empty-agentset             ; an empty agentset for clearing

  ; Maximum dispersal distances (measured in number of cells, not meters), used for generating cell-distance-file
  ; and calculated based on the distances set by max-cell-distance-___-m (which is in meters)
  ; FIX: fix these values later
  ; max-root-distance_cell         ; Max cell distance transmission for roots DELETE
  max-root-distance_long-rotation_cell       ; Max cell distance transmission for roots (long-rotation stands == 3 cell dists)
  max-root-distance_short-rotation_cell       ; Max cell distance transmission for roots (short-rotation stands == 4 cell dists)
  max-insect-distance_SC_cell      ; Max cell dispersal for S. carinatus
  max-insect-distance_PF_cell      ; Max cell dispersal for P. fasciatus
  max-insect-distance_HN_cell      ; Max cell dispersal for H. nigrinus
  max-inf-spread-distance_cell     ; Max cell distance used for spread (max of dispersal distances in cells
  max-attraction-distance_cell     ; max distance for attraction of insects (in cell units): 10 cells = 15 m / 49.2 ft

  ; BACKGROUND PROBABILITY OF insect DISPERSAL
  n-stands-projecting-insect-dispersal_background   ; the number of stands projecting background dispersal
                                                    ; occurs when the % infection is very high (for BSRD)
  prob-insect-dispersal_background_sum            ; the probability of insect dispersal to project based on the number of stands with dispersal spillover

  ;max-distance-used_cell       ; Max distance used for distance rings
  bsrd-inf_cumul               ; Count of all bsrd infections during model run (cumulative)
  bsrd-mort_cumul              ; Count of all bsrd mortalities during model run (cumulative)
  trees_cumul                  ; count of all trees ever existing in the model
  bsrd-inf_short-rotation_cumul           ; Count of all bsrd infections during model run in short-rotation stands (cumulative)
  bsrd-mort_short-rotation_cumul          ; Count of all bsrd mortalities during model run in short-rotation stands (cumulative)
  bsrd-inf_long-rotation_cumul           ; Count of all bsrd infections during model run in long-rotation stands (cumulative)
  bsrd-mort_long-rotation_cumul          ; Count of all bsrd mortalities during model run in long-rotation stands (cumulative)
  trees_short-rotation_cumul                  ; count of all trees ever existing in the model
  trees_long-rotation_cumul                  ; count of all trees ever existing in the model

  ; for checking parameters
  parameter-values

  ; Define colors for easy visual understanding of the model
  color_empty              ; Non-tree cells
  color_tree-live-notinf   ; Live trees: not infected
  color_tree-live-inf      ; Live trees: infected
  color_tree-dead-notinf   ; Dead trees: not infected
  color_tree-dead-inf      ; Dead trees: infected
  color_road               ; Color for road cells
  color_patches            ; Background color

;  errorFile ; delete
  kill-run? ; kills the run if a major error occurs
]

; Patch properties (used to quicky produce raster files of the landscape)
patches-own [
  p_cell             ; bool:     Whether or not a patch has a cell; used to avoid having >1 cell per patch
  p_seed             ; bin:      Whether or not a patch serves as a "seed" for "growing" a stand (0 = NO, 1 = YES)
  p_near-seed        ; bin:      Whether or not near a seed (to make stand sizes more equal)
  p_edge             ; bin:      Whether or not an edge (0: no, 1: yes)                                                            ; USED? CHECK
  p_road             ; bin:      Whether or not a road
  p_stand-id         ; int:      unique ID # for each stand
  p_mgmt             ; int:      management regime (1: short-rotation, 2: long-rotation, 3: old growth)
  p_tree             ; bin:      whether or not patch contains tree/stump (alive or dead)
  p_age              ; int:      tree age
  p_alive            ; bin:      If patch contains tree, whether patch has a live tree (1) or dead tree/stump (0). Default: 0
  p_inf              ; bin:      whether or not a tree is infected (0: no, 1: yes)
  p_inf-root         ; bin:      whether infected via roots
  p_inf-SC           ; bin:      whether infected via insect SC
  p_inf-PF           ; bin:      whether infected via insect PF
  p_inf-HN           ; bin:      whether infected via insect HN
  p_mort-cause       ; int:      cause of mortality: 1 == thin, 2 == harvest, 3 == BSRD
  p_p-inf_root       ; float:    probability of infection from root transmission
  p_p-inf_SC         ; float:    probability of infection from insect S. carinatus
  p_p-inf_PF         ; float:    probability of infection from insect P. fasciatus
  p_p-inf_HN         ; float:    probability of infection from insect H. nigrinus

  ;patch attraction properties
  p_attr-road        ; real:     Attraction factor when near a road
  p_attr-mgmt
  p_attr-mgmt_source
]

; Create a class ("breed") of turtles to serve as the hexagon cells of the model
breed [ cells cell ]

; Hexagon cell properties - see patch properties for info
cells-own [
  ;cell-nbors            ; agentset: The six neighboring hexagons ; DEACTIVATED
  ; STAND SETUP
  seed?                 ; bool:     Whether or not cell serves as a "seed" for "growing" a stand
  stand-id              ; int:      Unique ID # for each stand
  near-seed?            ; bool:     Whether or not cells are near a "seed", used to generate similarly sized stands
  edge?                 ; bool:     Whether or not an edge ; CHECK: Necessary?
  road?                 ; bool:     Whether or not a road
  mgmt                  ; str:      Management regime (i.e., short-rotation, long-rotation, set-aside)

  ; TREE STATUS
  potential-tree?       ; bOOL: used so that project-inf and attr affect those cells that are trees or will become trees post-regen
  tree?                 ; bool:     whether or not cell contains tree/stump (alive or dead)
  alive?                ; bool:     whether or not the tree is alive
  age                   ; int:      cell/stand age
  ; Infection
  infected?             ; bool:     whether or not a tree is infected
  newly-infected?       ; bool:     if a tree is infected during this time step (clear before next time step) ; check
  losing-infection?     ; bool:     if a tree is losing its infection during this time step ; check
  inf-initial?          ; bool:     whether tree became infection in initial setup
  inf-root?             ; bool:     whether tree became infected by roots
  inf-SC?               ; bool:     whether tree became infected by SC
  inf-PF?               ; bool:     whether tree became infected by PF
  inf-HN?               ; bool:     whether tree became infected by HN
  time-since-inf        ; int:      Time since infection, used to create a 1-2 year lag between infection and mortality
  mort-cause            ; str:      cause of mortality for a tree (used for harvest/regeneration and bsrd mort)
  time-since-mort       ; int:      Time since mortality, used to eliminate the viability of dead trees as a source of inoculum after 1-2 years

  age-bsrd-inf-suscept  ; real:    age-based infection susceptibility
  age-bsrd-mort-suscept ; real:    age-based mortality susceptibility

  ; PROBABILITY OF INFECTION
  at-risk?                     ; bool:     whether the cell is at risk to become infected: is (a) alive or recently dead AND (b) has yet to be infected by all mechs

  ; CHECK necessary?
  ; set-default-cell-properties: YES
  ; reset-cell-properties:       NO

  ; ROOTS
  ; Probability of root infection and projecting probability of root contact
  prob-inf-root                ; real:     probability of infection by root transmission
  prob-root-transm          ; one-of prob-root-transmission_list
  prob-root-contact_union   ; Calculated as the SUM of all prob-root-contact for individual trees
  any-inf-root-nbors?       ; used to track whether probability of root contact has changed and needs to be recalculated
  inf-root-nbors_ages_r1    ; used to store a list of ages of infected neighbors within the distance necessary for root contact and transmission
  inf-root-nbors_ages_r2    ; used to store a list of ages of infected neighbors within the distance necessary for root contact and transmission
  inf-root-nbors_ages_r3    ; used to store a list of ages of infected neighbors within the distance necessary for root contact and transmission
  inf-root-nbors_ages_r4    ; used to store a list of ages of infected neighbors within the distance necessary for root contact and transmission
  inf-root-nbors_ids_r1     ; used to store a list of ages of infected neighbors within the distance necessary for root contact and transmission
  inf-root-nbors_ids_r2     ; used to store a list of ages of infected neighbors within the distance necessary for root contact and transmission
  inf-root-nbors_ids_r3     ; used to store a list of ages of infected neighbors within the distance necessary for root contact and transmission
  inf-root-nbors_ids_r4     ; used to store a list of ages of infected neighbors within the distance necessary for root contact and transmission

  ; Probability of insect infection and projecting probability of insect dispersal
  prob-inf-SC                  ; real:     probability of infection by S. carinatus
  prob-inf-PF                  ; real:     probability of infection by P. fasciatus
  prob-inf-HN                  ; real:     probability of infection by H. nigrinus
  ;  recalc-insect?                      ; bool:     whether or not to recalculate the cell's infection probability (BETA) ; check delete
  ; base prob for insect infection = ( prob-insect-transmission * prob-insect-infested ). Set from the time the tree is initialized
  prob-inf_SC_base          ; (( one-of prob-insect-infested_SC_list ) * ( one-of prob-insect-transmission_SC_list ))
  prob-inf_PF_base          ; (( one-of prob-insect-infested_PF_list ) * ( one-of prob-insect-transmission_PF_list ))
  prob-inf_HN_base          ; (( one-of prob-insect-infested_HN_list ) * ( one-of prob-insect-transmission_HN_list ))

  prob-insect-disp_sum_SC             ; Summation of prob insect dispersal for S. carinatus (sum of prob dispersal from all infected trees)
  prob-insect-disp_sum_PF             ; Summation of prob insect dispersal for P. fasciatus '                                             '
  prob-insect-disp_sum_HN             ; Summation of prob insect dispersal for H. nigrinus  '                                             '
  n-infected-trees-disp-radius_SC     ; Count of the number of infected trees in the dispersal radius for SC
  n-infected-trees-disp-radius_PF     ; Count of the number of infected trees in the dispersal radius for PF - CHECK needed for all 3 or just SC and PF_HN
  n-infected-trees-disp-radius_HN     ; Count of the number of infected trees in the dispersal radius for PF - CHECK needed for all 3 or just SC and PF_HN
  ; CHECK new probabilities for union?

  ; ATTRACTION
  prob-insect-wound            ; list of real: probability of wounding based on the state of the tree (live/dead, inf/non-inf)
  attr-road                    ; list of real: factor by which wounding is increased when near a road
  attr-mgmt                    ; list of real: factor by which wounding is increased when in/near a thinned/harvested stand
  attr-mgmt_source             ; int:  the stand id of the stand that was thinned or managed to create the attraction
  stand-w-thin-disturbance?    ; bool: whether the stand has been recently thinned
  stand-w-harv-disturbance?    ; bool: whether the stand has been recently harvested
  losing-attr-dead?            ; bool: if a dead tree is losing its attractiveness to insects (after death and decay, no longer suitable host for pathogen)
  new-bsrd-mortality?          ; bool: if a tree has just died
  ; Yet to implement - CHECK ????
  attr-inf                     ; list of real: factor by which wounding is increased when near an infected tree
  attr-inf_source              ; list: cell IDs (who) of trees in attraction radius providing infection attraction
  attr-dead                    ; list of real: factor by which wounding is increased when near a dead tree
  attr-dead_source             ; list: cell IDs (who) of trees in attraction radius providing dead attraction

  ; DELETE
  ; test
  ; tree-count
  ; tree-count_decay
  ; CHECK TO SEE IF USED/NECESSARY
  ;stand-nbor-id    ; int:      Neighboring stand id ; ADD: Use later to get stand neighbors   ; check useD?
]

;__________________________________________________________________________________________________________________________________________________________
; USEFUL REPORTERS - USED IN PROGRESS FILE

; REPORTING WHOLE STANDS
to-report stand [ stand-id-no ]                                ; Reports an agentset of all cells with a given stand ID #
  report ( cells with [ stand-id = stand-id-no ] )
end

to-report stand-patches [ stand-id-no ]
  report ( patches with [ p_stand-id = stand-id-no ] )         ; Reports an agentset of all patches with a given stand ID #
end

to-report trees                                                ; Reports an agentset of all trees (cells)
  ; !!! Keep this one reporter as-is because it's used to set the trees_current agentset
  report ( cells with [ tree? ] )
end

to-report live-trees                                           ; Reports an agentset of live trees (infected an uninfected)
  report ( trees_current with [ alive? ] )
end

to-report live-infected-trees                                  ; reports an agentset of live, infected trees
  report ( trees_current with [ infected? and alive? ] )
end

to-report dead-trees                                           ; reports an agentset of dead trees (infected and uninfected)
  report ( trees_current with [ not alive? ] )
end

to-report dead-infected-trees                                  ; reports an agentset of dead, infected trees
  report ( trees_current with [ infected? and not alive? ] )
end

to-report infected-trees                                       ; reports an agentset of infected trees (alive and dead)
  report ( trees_current with [ infected? ] )
end

to-report newly-infected-trees                                 ; agentset of trees losing infection (some may already be non-trees)
  report ( trees_current with [ newly-infected? ] )
end

to-report losing-infection-trees                               ; agentset of trees losing infection (some may already be non-trees)
  report ( trees_potential with [ losing-infection? ] )
end


to-report root-infected-trees                                  ; reports an agentset of trees infected by roots (alive and dead)
  report ( trees_current with [ inf-root? ] )
end

to-report SC-infected-trees                                    ; reports an agentset of trees infected by SC (alive and dead)
  report ( trees_current with [ inf-SC? ] )
end

to-report PF-infected-trees                                    ; reports an agentset of trees infected by PF (alive and dead)
  report ( trees_current with [ inf-PF? ] )
end

to-report HN-infected-trees                                    ; reports an agentset of trees infected by HN (alive and dead)
  report ( trees_current with [ inf-HN? ] )
end

to-report susceptible-trees                                    ; reports an agentset of trees susceptible to infection (all live trees and dead trees that died less than 1 year ago)
  report ( trees_current with [ alive? or ( not alive? and time-since-mort <= dead-host-viability-duration_min ) ] )
end

to-report bsrd-killed-trees                                    ; reports an agentset of trees killed by BSRD
  report ( trees_current with [ mort-cause = "bsrd" ] )
end

to-report thinned-trees                                        ; reports an agentset of trees killed by thinning
  report ( trees_current with [ mort-cause = "thin" ] )
end

to-report harvested-trees                                      ; reports an agentset of trees killed by harvest
  report ( trees_current with [ mort-cause = "harvest" ] )
end

to reset-globals
  reset-ticks
  ; CHECK - Are all globals reset? are all globals necessary?
  ;detect-computer               ; detect which computer is being used, then set the model directory accordingly; FIX removed because it was taking minutes to do this...
  if comp = "asus-big"           [ set model-dir "D:/AdamB/Documents/Masters/Thesis/Modeling/Models/BSRD_NetLogo/BSRD-full-model/" ]         ; Adam's big ASUS
  if comp = "asus-small"         [ set model-dir "C:/Users/ASUS/Documents/Masters/Thesis/Modeling/Models/BSRD_NetLogo/BSRD-full-model/" ]    ; Adam's ASUS Zenbook
  if comp = "law"                [ set model-dir "/home/bouchea/" ]                                                                          ; LAW cluster

  set run-id "r0"            ; unique run id number
  set run-file ""            ; path to run file
  set progress-file   ""     ; file to keep track of model progress
  set stand-data-file ""     ; file for stand-level data output
  set tree-data-file  ""     ; file to store prob-inf data at the cell scale
  set projection ( word "supporting_files/projections/wgs_84.prj" )     ; GIS projection file
  set test-timer 0
  set current-date-and-time ""               ; Abbreviated timestamp of current date and time
  set intercell-distance 1.524               ; Global to set distance between cells (1.524 meter resolution) = patch dimensions
  set cell-area ( intercell-distance ^ 2 )   ; Area of each patch (m^2)

  ; LANDSCAPE GENERATOR / MANAGEMENT
  set seed-radius 1
  ; Stand lists
  set stand-list                      []
  set short-rotation-stand-list       []
  set long-rotation-stand-list        []
  set set-aside-stand-list            []
  set regen-list-short-rotation       []
  set regen-list-long-rotation        []
  set thinning-occurred?              false       ; bool to track whether a thinning occurred
  set harvest-occurred?               false       ; bool to track whether a harvest occurred

  ; BACKGROUND PROBABILITY OF DISPERSAL FOR INSECTS (PF & HN)
  set n-stands-projecting-insect-dispersal_background 0
  set prob-insect-dispersal_background_sum          0

  ; Lists to track management attraction decay
  set thin-0yr-ago    []
  set thin-1yr-ago    []
  set thin-2yr-ago    []
  set thin-3yr-ago    []
  ;set thin-4yr-ago    []
  set harv-0yr-ago    []
  set harv-1yr-ago    []
  set harv-2yr-ago    []
  set harv-3yr-ago    []
  ;set harv-4yr-ago    []

  set plantation-road-density         34.93                               ; density of roads in meters per hectare (34.2 m/ha, from John Sessions, personal communication)
  ;set set-aside-road-density          14.92                               ; density of roads in meters per hectare (14.92 m/ha)
  set empty-agentset                  cells with [ tree? and not tree? ]  ; An impossible combination
  set trees-to-thin                   empty-agentset                      ; agentset: all trees to be thinned
  set trees-to-harvest                empty-agentset                      ; agentset: all trees to be thinned
  set edges                           empty-agentset                      ; an agentset of all hex cells considered to be at a stand edge ; REMOVE?
  set trees_current                   empty-agentset                      ; an agentset of all trees
  set trees_potential                 empty-agentset                      ; an agentset of all cells than will be trees during a model run
  set trees_newly-infected            empty-agentset                      ; an agentset of all trees that are newly infected
  set trees_losing-infection          empty-agentset                      ; an agentset of all trees that lost infection (dead trees)
  ;set trees_losing-dead-attraction    empty-agentset                      ; an agentset of all trees that lost infection (dead trees)
  set trees_infected                  empty-agentset                      ; an agentset of all trees that are infected ; check necessary
  set cells-projecting                empty-agentset

  ; INFECTION MODEL
  ; PROB PARAM LISTS
  reset-infection-parameters
  ; Maximum dispersal distances (measured in cells, not meters), used for generating cell-distance-file
  set-max-cell-distances

  ; FIX - why do these differ from live inf counts?
  set bsrd-inf_cumul                     0
  set bsrd-mort_cumul                    0
  set bsrd-inf_short-rotation_cumul      0
  set bsrd-mort_short-rotation_cumul     0
  set bsrd-inf_long-rotation_cumul       0
  set bsrd-mort_long-rotation_cumul      0
  set trees_cumul                        0
  set trees_short-rotation_cumul         0         ; count of all trees ever existing in the model
  set trees_long-rotation_cumul          0         ; count of all trees ever existing in the model

  ; KILL RUN - kills the run if an error occurs that invalidate the results but the model doesn't crash
  set kill-run? false
end


; Set the default/initial value of cell variables
to set-default-cell-properties ; with the exception of some properties that need to be set separately (i.e., cell-nbors)
  if visualize? [ set color color_empty ]       ; set cells to default (empty) cell color
  ;set cell-nbors          [] ; deactivated
  ; STAND SETUP
  set seed?                  false
  set near-seed?             false
  set edge?                  false
  set stand-id                0
  set mgmt                   ""
  set road?                  false

  ; STATUS
  set tree?                  false
  set alive?                 false
  set age                     0
  ; Infection
  set infected?              false
  set newly-infected?        false
  set losing-infection?      false
  set inf-initial?           false
  set inf-root?              false
  set inf-SC?                false
  set inf-PF?                false
  set inf-HN?                false
  set time-since-inf          0
  set age-bsrd-inf-suscept    0.0
  set age-bsrd-mort-suscept   0.0
  set prob-inf-root           0.0
  set prob-inf-SC             0.0
  set prob-inf-PF             0.0
  set prob-inf-HN             0.0
  ; Attraction
  set prob-insect-wound       0.0
  set attr-road               1
  set attr-mgmt            [  1 ]
  set attr-mgmt_source     [ -1 ] ; impossible cell ID value
  ifelse attr-inf-dead_version = "proportional" [
    set attr-inf             1
    set attr-inf_source      0
  ]
  [
    set attr-inf             [  1 ]
    set attr-inf_source      [ -1 ] ; impossible cell ID value
  ]
  ; check move into previous
  set attr-dead                        [  1 ]
  set attr-dead_source                 [ -1 ] ; impossible cell ID value
  set stand-w-thin-disturbance?        false ; bool: whether the stand has been recently thinned
  set stand-w-harv-disturbance?        false ; bool: whether the stand has been recently harvested
  set losing-attr-dead?                false
  set new-bsrd-mortality?              false
  set mort-cause                       ""
  set time-since-mort                  0

  set prob-root-transm                 0
  set prob-root-contact_union          0
  set prob-inf_SC_base                 0
  set prob-inf_PF_base                 0
  set prob-inf_HN_base                 0
  set prob-insect-disp_sum_SC          0
  set prob-insect-disp_sum_PF          0
  set prob-insect-disp_sum_HN          0
  set n-infected-trees-disp-radius_SC  0
  set n-infected-trees-disp-radius_PF  0
  set n-infected-trees-disp-radius_HN  0
  set at-risk?                         false
  set potential-tree?                  false

  ; For calculating probability of root contact
  set inf-root-nbors_ages_r1 []   ; used to store a list of ages of infected neighbors within the distance necessary for root contact and transmission
  set inf-root-nbors_ages_r2 []   ; used to store a list of ages of infected neighbors within the distance necessary for root contact and transmission
  set inf-root-nbors_ages_r3 []   ; used to store a list of ages of infected neighbors within the distance necessary for root contact and transmission
  set inf-root-nbors_ages_r4 []   ; used to store a list of ages of infected neighbors within the distance necessary for root contact and transmission
  set inf-root-nbors_ids_r1  []   ; used to store a list of ages of infected neighbors within the distance necessary for root contact and transmission
  set inf-root-nbors_ids_r2  []   ; used to store a list of ages of infected neighbors within the distance necessary for root contact and transmission
  set inf-root-nbors_ids_r3  []   ; used to store a list of ages of infected neighbors within the distance necessary for root contact and transmission
  set inf-root-nbors_ids_r4  []   ; used to store a list of ages of infected neighbors within the distance necessary for root contact and transmission
  set any-inf-root-nbors? false

end

; RESET CELL PROPERTIES: Turns trees into empty cells (e.g., after tree has died and disappeared)
to reset-cell-properties ; clears cell properties except for stand id, mgmt, road, and attraction ; CHECK - does this include everything?
  if visualize? [ set color color_empty ]
  set tree?                         false               ; true/false:  whether or not the cell has a tree
  set alive?                        false               ; true/false:  whether or not the cell is alive
  set age                           0                   ; integer:     the age of the cell
  set infected?                     false               ; true/false:  whether or not the cell is infected
  set newly-infected?               false
  set inf-initial?                  false
  set inf-root?                     false
  set inf-SC?                       false
  set inf-PF?                       false
  set inf-HN?                       false
  set mort-cause                    ""                  ; list:        cause of mortality (e.g., thin, harvest, bsrd)
  set prob-inf-SC                   0.0                 ; real:        probability of infection by S. carinatus
  set prob-inf-PF                   0.0                 ; real:        probability of infection by P. fasciatus
  set prob-inf-HN                   0.0                 ; real:        probability of infection by H. nigrinus
  set time-since-inf                0                   ; int:         used to create a 1-2 year lag between infection and mortality
  if mort-cause != "harvest" [ set time-since-mort 0 ]  ; int:         used to create a lag between mortality and infection loss, attraction, etc. Not used for harvest to allow for regen
  set at-risk?                     false                ;
;  set recalc-insect?              false               ; bool:     whether of not there's any need to recalculate prob-inf because prob-exposure has changed

  ; check - should attr be reset here in any circumstance? Or not nec?

  ; DON'T RESET
  ; attr-mgmt, road, dead, or inf
  ; inf-root-nbor age or id lists
  ; losing-infection? losing-attr-dead? or new-bsrd-mortality?
  ; prob-insect-disp_sum for any insect, nor n-infected-trees-disp-radius for any insect
end

; Set the default/initial value of patch variables
to set-default-patch-properties
  if visualize? [ set pcolor color_patches ]
  set p_cell     false
  set p_seed       0
  set p_near-seed  0
  set p_edge       0
  set p_road       0
  set p_stand-id   0
  set p_mgmt       0
  set p_tree       0
  set p_age        0
  set p_alive      0
  set p_inf        0
  set p_inf-root   0
  set p_inf-SC     0
  set p_inf-PF     0
  set p_inf-HN     0
  set p_mort-cause 0

  set p_p-inf_root 0.0
  set p_p-inf_SC   0.0
  set p_p-inf_PF   0.0
  set p_p-inf_HN   0.0
  ;patch attraction properties
  set p_attr-road  0.0
  set p_attr-mgmt  0.0
  set p_attr-mgmt_source 0.0
end

; SET DEFAULT VALUES FOR USER-DEFINED VARIABLES
; CHECK: Does this affect headless runs?
to startup
  no-display
  ;detect-computer
  set-patch-size                           1E-10
  set base                                "patch"
  set visualize?                           false
  set limit-ticks?                         true
  set max-ticks                            160
  set config                              "random"
  set n-stands                               1
  set pct-short-rotation                    50
  set pct-long-rotation                     50
  set pct-set-aside                          0
  set rotation:short-rotation               37
  set rotation:long-rotation                80
  set replant-delay                          1
  set spread-setup?                        true
  set age-inf-suscept                      "age-indep"
  set dead-attr-duration_min                 0 ; default 2 yr
  set dead-attr-duration_max                 0 ; default 4 yr
  set dead-host-viability-duration_min       0 ; default 1 yr
  set dead-host-viability-duration_max       0 ; default 2 yr
  set pct-initial-infection                  0.5
  set reinitiate-infection?                false
  set initial-infection                    "generate-initial-infection-pct"
  set max-root-distance_m                    8
  set max-insect-distance_SC_m               8
  set max-insect-distance_PF_m               8
  set max-insect-distance_HN_m               8
  set max-attraction-distance_m             15
  set sensitivity-analysis-param           "none"
  set sensitivity-analysis-multiplier        1.0
  set track-run?                           false
  set export-tree-data?                    false
  set tree-data-export-freq                300
  set export-rasters?                      false
  set raster-export-frequency              300
end

to detect-computer
  ; Detect computer based on folders present
  ifelse ( file-exists?  "D:/AdamB"  )                  [ set comp "asus-big"    ]
  [
    if ( file-exists?  "C:/Users/ASUS"  )               [ set comp "asus-small"  ]
    if ( file-exists?  "/home/bouchea/law.txt"  )       [ set comp "law" ]
  ]
end

to create-run-id ; used in generate-landscape to create a unique run identifier based on the date/time and configuration
  create-datetimestamp
  set run-id-no (word "r" ( random 100000 ) )
  set run-id (word run-id-no "_" current-date-and-time "_dm-" max-pxcor "-" max-pycor "_" config "_" n-stands "_L-I-" pct-short-rotation "-E-" pct-long-rotation "_" "SNS-" sensitivity-analysis-param "-mult-" sensitivity-analysis-multiplier  )
end

to create-run-files
  ; Detect computer
  ;detect-computer
  ; Change model directory to reflect comp being used
  if comp = "asus-big"       [ set model-dir "D:/AdamB/Documents/Masters/Thesis/Modeling/Models/BSRD_NetLogo/BSRD-full-model/" ]         ; Adam's big ASUS
  if comp = "asus-small"     [ set model-dir "C:/Users/ASUS/Documents/Masters/Thesis/Modeling/Models/BSRD_NetLogo/BSRD-full-model/" ]    ; Adam's ASUS Zenbook
  if comp = "law"            [ set model-dir "/home/bouchea/models/" ]                                                                   ; LAW cluster

  set run-file ( word "output/" run-id "_run.txt" ) ; Set the path for the run info
  file-open  run-file
  write-run-details
  file-close
  ; Print out run details in log file for reference
  file-open run-file
  while [ not file-at-end? ] [print file-read-line ]
  file-close

  ; initialize progress file with column names
  initialize-progress-file

  ; initialize stand data file with column names
  if export-stand-data? [ initialize-stand-data-file ]
end

to write-run-details
  file-print ( word "Run ID: " run-id )
  file-print "----BASICS----\n"
  file-print ( word "Start: " date-and-time )
  file-print ( word "Time stamp: " current-date-and-time )
  file-print ( word "Projection: " projection )

  file-print "\n----RUN INFO----\n"
  file-print ( word "comp: " comp )
  file-print ( word "set-world-size? " set-world-size? )
  if set-world-size? = "preset" [
    file-print ( word "world-size: " world-size )
  ]
  if set-world-size? = "automatic" [
    file-print ( word "stands:                  "           n-stands                           " stands"   )
    file-print ( word "mean stand size (acres): "           stand-size_mean_acres              " acres"   )
    file-print ( word "mean stand size (ha):    "         ( stand-size_mean_acres * 0.404686 ) " ha"      )
    file-print ( word "World size:      "                 ( stand-size_mean_acres * n-stands ) " acre(s)" )
  ]
  file-print ( word "n-stands: "                          n-stands )
  file-print ( word "World dim - x_min: "                 min-pxcor )
  file-print ( word "World dim - x_max: "                 max-pxcor )
  file-print ( word "World dim - y_min: "                 min-pycor )
  file-print ( word "World dim - y_max: "                 max-pycor )
  file-print ( word "limit-ticks? "                       limit-ticks? )
  file-print ( word "max-ticks: "                         max-ticks )
  file-print ( word "visualize? "                         visualize? )
  file-print ( word "track-run? "                         track-run? )
  file-print ( word "export-tree-data? "                  export-tree-data? )
  file-print ( word "tree-data-export-freq: "             tree-data-export-freq )
  file-print ( word "export-rasters? "                    export-rasters? )
  file-print ( word "raster-export-frequency: "           raster-export-frequency )
  file-print ( word "base: "                              base )
  file-print ( word "config: "                            config )
  file-print ( word "n-blocks: "                          n-blocks )
  file-print ( word "pct-short-rotation: "                     pct-short-rotation )
  file-print ( word "pct-long-rotation: "                     pct-long-rotation )
  file-print ( word "pct-set-aside: "                     pct-set-aside )
  file-print ( word "rotation:short-rotation: "                rotation:short-rotation )
  file-print ( word "rotation:long-rotation: "                rotation:long-rotation )
  file-print ( word "thinning? "                          thinning? )
  file-print ( word "harvest? "                           harvest? )
  file-print ( word "replant-delay: "                     replant-delay )
  file-print ( word "spread-setup? "                      spread-setup? )
  file-print ( word "spread-infection-version? "          spread-infection-version? )
  file-print ( word "initial-infection: "                 initial-infection )
  file-print ( word "pct-initial-infection: "             pct-initial-infection )
  file-print ( word "reinitiate-infection? "              reinitiate-infection? )
  file-print ( word "sensitivity-analysis-param: "        sensitivity-analysis-param )
  file-print ( word "sensitivity-analysis-multiplier: "   sensitivity-analysis-multiplier )
  file-print ( word "max-attraction-distance_m: "         max-attraction-distance_m )
  file-print ( word "max-root-distance_m: "               max-root-distance_m )
  file-print ( word "max-insect-distance_SC_m: "          max-insect-distance_SC_m )
  file-print ( word "max-insect-distance_PF_m: "          max-insect-distance_PF_m )
  file-print ( word "max-insect-distance_HN_m: "          max-insect-distance_HN_m )
  file-print ( word "background-prob-insect-dispersal? " background-prob-insect-dispersal? )
  file-print ( word "insect-disp_background_stand-pct-inf-spillover-threshold: " insect-disp_background_stand-pct-inf-spillover-threshold )
  file-print ( word "prob-insect-disp_background: "      prob-insect-disp_background )
  file-print ( word "age-inf-suscept: "                  age-inf-suscept )
  file-print ( word "export-jsons?   "                   export-jsons?  )
  file-print ( word "export-potential-trees? : "         export-potential-trees?  )
  file-print ( word "debug-mode?   "                     debug-mode?  )
  file-print ( "----- Adjustable parameter settings -----" )
  file-print ( word "prc_conservative?"                  prc_conservative? )
  file-print ( word "prt-exclude-data? "                 prt-exclude-data? )
  file-print ( word "prt-discount-by-contact-type? "     prt-discount-by-contact-type? )
  file-print ( word "new-attr? "                         new-attr? )
  file-print ( word "attr-inf? "                         attr-inf? )
  file-print ( word "attr-dead? "                        attr-dead? )
  file-print ( word "attr-mgmt_conservative?"            attr-mgmt_conservative? )
  file-print ( word "attr-inf-dead_version: "            attr-inf-dead_version )
  file-print ( word "inf-center-attr-threshold: "        inf-center-attr-threshold )
  file-print ( word "inf-center-attr-threshold_pct: "    inf-center-attr-threshold_pct )
  file-print ( word "dead-attr-duration (min): "         dead-attr-duration_min )                        ; data indicates: 2   yr
  file-print ( word "dead-attr-duration (max): "         dead-attr-duration_max )                        ; data indicates: 4   yr
  file-print ( word "dead-host-viability-duration_min: " dead-host-viability-duration_min )              ; data indicates: 0.5 yr
  file-print ( word "dead-host-viability-duration_max: " dead-host-viability-duration_max )              ; data indicates: 2   yr
end

to complete-run-file
  if run-file != 0 and run-file != "" [
    file-open run-file
    file-print ( word "End: " date-and-time )
    file-close
  ]
  print ( word "End: " date-and-time )
end
;__________________________________________________________________________________________________________________________________________________
;__________________________________________________________________________________________________________________________________________________
;__________________________________________________________________________________________________________________________________________________
;__________________________________________________________________________________________________________________________________________________
; MODEL SETUP
;__________________________________________________________________________________________________________________________________________________
;__________________________________________________________________________________________________________________________________________________
;; LANDSCAPE GENERATION PROCEDURES:
;__________________________________________________________________________________________________________________________________________________
; AGE DISTRIBUTIONS

; AGE DISTRIBUTION (TO DRAW STAND AGES FROM)
to-report age-dist:short-rotation
  report ( random rotation:short-rotation + 1 )
end

to-report age-dist:long-rotation
  report ( random rotation:long-rotation + 1 )
end

;__________________________________________________________________________________________________________________________________________________
; LANDSCAPE SETUP

; SETUP THE BASE FOR THE WORLD
to setup-base
  ; Clear the model
  no-display                                            ; Turn off display: keeps the setup running faster
  clear-all                                             ; Clears everything
  r:clear                                               ; clear the R environent
  reset-timer
  reset-ticks                                           ; Reset the model ticks (time steps) to 0
  reset-globals                                         ; Set all globals to default values, including globals that will not change
  if ( patch-size != 1E-10 ) [ set-patch-size 1E-10 ]   ; Set the patch size very small to avoid running out of memory w/big worlds
  if set-world-size? != "manual" [ set-world-size ]     ; Set the dimentions of the world
  if visualize?      [ setup-colors set-default-shape cells "hex" ]                   ; Set default colors and those that represent cell states
  ; set default shape for cells. The hexagon shape must be imported if this is a new
  ; installation of NetLogo. This can be done from the menu bar >> Tools >> Turtle Shapes Editor
  set-default-shape cells "hex"

  ; Set GIS settings
  gis:load-coordinate-system ( projection )             ; Set the coordinate system (default WGS 1984)
  gis:set-world-envelope [ 0 4503.42 0 4503.42 ]        ; World dimensions in meters ( 2955 cells * 1.524 m/cell = 4503.42 m per side)

  ask patches [ set-default-patch-properties ] ; Set default patch properties

  if visualize? [
    color-world
    display ; Turn the display back on when the landscape is finished.
  ]
  print ( word "Base generated successfully @ " date-and-time "; Duration: " timer )
  set base-setup? true
end

; BASE SETUP SUB-PROCEDURES - these are model procedures used within both landscape generation procedures

; SET WORLD SIZE
to set-world-size ; DELETE: relax the world size constraint
  if set-world-size? = "pre-set" [
    if world-size =     "10 x 10" [ resize-world 0    9 0    9   set-patch-size 1E-10 ]
    if world-size =     "26 x 26" [ resize-world 0   25 0   25   set-patch-size 1E-10 ]
    if world-size =     "50 x 50" [ resize-world 0   49 0   49   set-patch-size 1E-10 ]
    if world-size =   "100 x 100" [ resize-world 0   99 0   99   set-patch-size 1E-10 ]
    if world-size =   "300 x 300" [ resize-world 0  299 0  299   set-patch-size 1E-10 ]
    if world-size =   "500 x 500" [ resize-world 0  499 0  499   set-patch-size 1E-10 ]
    if world-size = "1000 x 1000" [ resize-world 0  999 0  999   set-patch-size 1E-10 ]
    if world-size = "2952 x 2952" [ resize-world 0 2951 0 2951   set-patch-size 1E-10 ] ; 100, 50-acre stands
    if world-size = "3616 x 3616" [ resize-world 0 3615 0 3615   set-patch-size 1E-10 ] ; 150, 50-acre stands
    if world-size = "4176 x 4176" [ resize-world 0 4175 0 4175   set-patch-size 1E-10 ] ; 200, 50-acre stands
  ]
  if set-world-size? = "automatic" [
    if config = "blocks" [ user-message "Blocks requires manual or pre-set setup, sorry!" set kill-run? true stop ]
    let stand-size_mean_ha      ( stand-size_mean_acres * 0.404686 )                      ; convert acres to ha
    let world-size_ha           ( stand-size_mean_ha * n-stands )                         ; multiply by n-stands to get world size
    let world-size_m2           ( world-size_ha * 10000 )                                 ; convert ha to m^2
    let world-dim_m             ( sqrt world-size_m2 )                                    ; calc dimensions of a square with this area
    let world-dim_patches       ( ceiling ( world-dim_m / 1.524 ) )                       ; calc n patches necessary for this world (rounding up)
    if world-dim_patches mod 2 = 1 [ set world-dim_patches ( world-dim_patches + 1 ) ]    ; avoid having a misalignment at world edges
    resize-world 0 ( world-dim_patches - 1 ) 0 ( world-dim_patches - 1 )                  ; Resize the world

    ; Output values for documentation
    print ( word "Mean stand size: " stand-size_mean_acres                     " acre(s)"   )
    print ( word "Mean stand size: " stand-size_mean_ha                        " ha"        )
    print ( word "World size:      " (stand-size_mean_acres * n-stands )       " acre(s)"   )
    print ( word "World size:      " stand-size_mean_ha                        " ha"        )
    print ( word "World size:      " world-size_m2                             " m^2"       )
    print ( word "world dim:       " world-dim_m " x " world-dim_m             " m"         )
    print ( word "world dim:       " world-dim_patches " x " world-dim_patches " patch(es)" )
  ]
end

to check-world-size
  let stand-size_mean_ha      ( stand-size_mean_acres * 0.404686 )                      ; convert acres to ha
  let world-size_ha           ( stand-size_mean_ha * n-stands )                         ; multiply by n-stands to get world size
  let world-size_m2           ( world-size_ha * 10000 )                                 ; convert ha to m^2
  let world-dim_m             ( sqrt world-size_m2 )                                    ; calc dimensions of a square with this area
  let world-dim_patches       ( ceiling ( world-dim_m / 1.524 ) )                       ; calc n patches necessary for this world (rounding up)
  if world-dim_patches mod 2 = 1 [ set world-dim_patches ( world-dim_patches + 1 ) ]    ; avoid having a misalignment at world edges

  print ( word "Mean stand size: " stand-size_mean_acres                     " acre(s)"   )
  print ( word "Mean stand size: " stand-size_mean_ha                        " ha"        )
  print ( word "World size:      " (stand-size_mean_acres * n-stands )       " acre(s)"   )
  print ( word "World size:      " stand-size_mean_ha                        " ha"        )
  print ( word "World size:      " world-size_m2                             " m^2"       )
  print ( word "world dim:       " world-dim_m " x " world-dim_m             " m"         )
  print ( word "world dim:       " world-dim_patches " x " world-dim_patches " patch(es)" )
end

;_____________________________________________________________________________________________________________
; Clear the landscape with the appropriate process without clearing ring globals (if already read)
to reset
  no-display
  ifelse ( base-setup? = true ) [ reset-globals ] [ setup-base print "base being set up from clear landscape" ] ; If base set up: reset globals, If not: setup-base
  reset-ticks      ; Reset tick count
  ask patches [ set-default-patch-properties ]   ; Reset all patch properties
  if any? cells [ ask cells [ die ] ]                          ; Deletes all cells
  if visualize? [ setup-plots update-plots display ]
end

; "Setup" function (once the cell grid has been created). Generates the landscape to user specifications.
to setup
  no-display
  ;____________________________________________________
  ; BEGIN SETUP PROCEDURE
  let approx-world-size_acres ( ( ( world-width * 1.524 ) ^ 2 ) / 10000 / 0.404686 )
  let desired-world-size_acres ( n-stands * stand-size_mean_acres )

  carefully [
    set approx-world-size_acres ( ( ( world-width * 1.524 ) ^ 2 ) / 10000 / 0.404686 )
    set desired-world-size_acres ( n-stands * stand-size_mean_acres )

    if ( base-setup? != true or
      ( set-world-size? = "pre-set" and not member? ( word world-width) world-size ) OR
      ( set-world-size? = "automatic" AND ( ( approx-world-size_acres / desired-world-size_acres ) < 0.97 OR ( approx-world-size_acres / desired-world-size_acres ) > 1.03 ) ) ; make sure world size is within 3% of what it should be
    ) [ setup-base print "base being set up from generate landscape" ]                             ; Setup base if it has not been set up
    if ( any? patches with [ p_stand-id != 0 ] or any? cells ) [ reset ]                           ; Resets globals and removes stand from previous landscape

    ; CHECK TO MAKE SURE THE SETTINGS WILL ALLOW A LANDSCAPE TO BE SET UP
    if config = "random" OR config = "clustered" AND n-stands <= 0  [ user-message "Please set n-stands to a number > 0." set kill-run? true stop ]
    if config = "random" OR config = "clustered" AND n-stands >= ( 10 * count patches ) [ user-message "Too many stands for a landscape of this size." set kill-run? true stop ]
    let total-of-mgmt-proportions ( pct-short-rotation + pct-long-rotation + pct-set-aside )
    if total-of-mgmt-proportions != 100 [ user-message "Percentages of forest management types do not add to 100%. Try again!" print "Percentages of forest management types do not add to 100%. Try again!" set kill-run? true stop ]
  ] [ print "ERROR: setup -> initial" print error-message ]
  reset-timer reset-ticks

  ;____________________________________________________
  ; CREATE RUN ID
  create-run-id
  ;____________________________________________________

  ; SET UP FOR SPREADING INFECTION
  set-max-cell-distances                      ; Find the max rings necessary based on spread distance settings
  print "set-max-cell-distances: successful"
  setup-infection-parameters                  ; setup lists that store param values; includes attraction setting
  print "setup-infection-parameters: successful"
  if ( spread-setup? = true ) [ set-ring-lists_R ]       ; NEW! Now more efficient
  print "set-ring-lists_R: successful"
  ;____________________________________________________

  ; SET UP THE LANDSCAPE
  ; Create the stands in the landscape
  generate-stands
  print "generate-stands: successful"
  ; CRATE AN INITIAL INFECTION
  run initial-infection                             ; set the agentset of infected trees
  print "generate initial-infection: successful"
  set trees_infected ( infected-trees )
  set trees_newly-infected ( newly-infected-trees )
  ; Landscape setup is complete...

  ;____________________________________________________

  ; Set up for Go spread infection
  create-spread-infection-directory
  print "Directories set up for Go"
  if ( spread-setup? and spread-infection-version? = "Go" and
    ( ( max-insect-distance_SC_m + max-insect-distance_SC_m + max-insect-distance_SC_m ) > 0 )
    ) [
    write-spread-infection-ouputs_initial
    print "Initial go outputs written"
  ]

  ;____________________________________________________

  ; RUN TRACKING AND DATA EXPORT: Keep track of run progress and export data
  if track-run? [ create-run-files ]
  let duration timer
  if run-file != 0 and run-file != "" [
    file-open  run-file
    file-print ( word "Run ID: " run-id )
    file-print ( word "The landscape has been generated successfully @" date-and-time "; Duration: " duration )
    file-print (word "% short-rotation: " pct-short-rotation "; % Exensive:  " pct-long-rotation )
    file-close-all
  ]
  print (word "Run ID: " run-id)
  print ( word "The landscape has been generated successfully @" date-and-time "; Duration: " duration )
  print (word "% short-rotation: " pct-short-rotation "; % Exensive:  " pct-long-rotation )

  if export-rasters? [ export-stand-map ] ; Export a raster of management, age, etc

  if visualize? [
    color-world
    setup-plots update-plots display
  ]
  if track-run? [  write-settings-file ]

  ; DELETE
;  set errorFile ( word "fixThis_" random 100 ".csv" )
;  if file-exists? errorFile [ file-delete errorFile ]
;  file-open errorFile
;  file-print "tick,who,tree,alive,inf,mortCause,newlyInfSet,losingInfSet,procedure"
;  print (word "errorFile: " errorFile )
;  file-close
end

to generate-stands
  if config = "random" [
    setup-random-stands    ; Creates blocks in the landscape that become stands
    assign-mgmt            ; Stands are assigned a management class
    create-initial-trees   ; Trees are created based on the management and age of the stand
  ]
  ;if config = "blocks"    [ ; complete later
  ;if config = "clustered" [ ; complete later
  setup-roads

  ; Realign cells in the world
  ;ask cells [ ifelse (pxcor mod 2 = 0) [ set ycor ( ycor - 0.25 ) ] [ set ycor ( ycor + 0.25 ) ] ]    ; if even column (based on pxcor), shift down one half coordinate, else realign cells in world (gets cells off the boundary between patches)
  ask cells [ if (pxcor mod 2 = 0) [ set ycor ( ycor - 0.5 ) ] ]    ; if even column (based on pxcor), shift down one half coordinate, else realign cells in world (gets cells off the boundary between patches)
  set trees_potential ( cells with [ potential-tree? ] )
  set trees_current ( trees )
  print ( word "Count trees: " count trees_current )

  print "generate-stands complete"
end

; RANDOM STAND LANDSCAPE
; Stand origins are randomly seeded with some constraints to make them evenly-sized, and the spatial arrangement of
; management classes is random
to setup-random-stands
  ifelse n-stands > 1 [
    setup-seeds        ;; create the seed cells from which the stands will "grow"
    grow-seeds         ;; "grow" the stands from the "seeds"
  ]
  [
    set stand-list [ 1 ]
    ask patches [ set p_stand-id 1 ]
  ]
end

to setup-seeds
  ; Create a list of numbers from 1 to the total number of stands as set by
  ; "n-stands" (increment of 1)
  set stand-list ( range 1 ( n-stands + 1) )

  set seed-radius  round ( 1.2 * sqrt ( ( ( max-pxcor * max-pycor ) / n-stands ) / pi ) ) ; Create a radius in which no other patches will serve as stand seeds
  ; Calculate total number of patches per stand, calculate an appropriate radius

  ; Select a set of random cells of size "n-stands", assign one to seed each stand
  foreach stand-list [
    this-stand ->
    ask one-of patches with [ p_near-seed = 0 ] [
      set p_stand-id this-stand  ; set the cell's stand id
      set p_seed 1               ; change property to indicate that this cell is a seed
      ask patches in-radius seed-radius [ set p_near-seed 1 ] ; makes a minimum radius for stand-building, making stands more even. If this number is too high relative to
      ; world size and stand number, space for new stand seeds will run out! FIX
    ]
  ]
end

to grow-seeds
  ; To fill the landscape: while there are any patches with no stand id (stand-id = 0), have each seed (and later, the stand patches)
  ; set all unassigned neighbors to their stands
  while [ any? patches with [ p_stand-id = 0 ] ] [
    ;; check if setup looks the same when ask concurrent is used
    ask patches with [ p_stand-id != 0 ] [
      ask neighbors with [ p_stand-id = 0 ] [
        set p_stand-id [ p_stand-id ] of myself
      ]
    ]
  ]
end

to assign-mgmt
  let remaining-stands stand-list ; Create a temporary variable 'remaining stands' that is a copy of the stand list
  set short-rotation-stand-list ( n-of ( pct-short-rotation / 100 * length ( stand-list ) ) remaining-stands ) ; subset of remaining stands to be short-rotation stands
  foreach short-rotation-stand-list [ this-stand ->
    ask stand-patches this-stand [ set p_mgmt 1 ]             ; set mgmt
    set remaining-stands remove this-stand remaining-stands   ; remove from remaining stand list
  ]  ; Repeat for long-rotation stands and old growth
  set long-rotation-stand-list ( n-of ( pct-long-rotation / 100 * length ( stand-list ) ) remaining-stands )
  foreach long-rotation-stand-list [ this-stand ->
    ask stand-patches this-stand [ set p_mgmt 2 ]    ; set mgmt
    set remaining-stands remove this-stand remaining-stands   ; remove from remaining stand list
  ]
  set set-aside-stand-list ( n-of ( pct-set-aside / 100 * length ( stand-list ) ) remaining-stands )
  foreach set-aside-stand-list [ this-stand ->
    ask stand-patches this-stand [ set p_mgmt 3 ]
    set remaining-stands remove this-stand remaining-stands
  ]
  foreach remaining-stands [ this-stand -> ; Randomly assign all unassigned stands one of the appropriate mgmt options based on the settings
    ; This occurs when the proportions of each forest type do not work 'perfectly' with the number of stands (i.e., when the proportion * n-stands is not an whole number)
    let mgmt-options [] ; empty list of possible management classes (based on the user-defined percentages)
    if pct-short-rotation > 0 [ set mgmt-options fput 1 mgmt-options ] ; add the appropriate mgmt classes to the
    if pct-long-rotation > 0 [ set mgmt-options fput 2 mgmt-options ]
    if pct-set-aside > 0 [ set mgmt-options fput 3 mgmt-options ]
    let mgmt-to-set ( one-of mgmt-options )                              ; choose 1 mgmt type
    ask stand-patches this-stand [ set p_mgmt mgmt-to-set ]                        ; set it
    ; assign to a list based on what mgmt was randomly chosen
    if mgmt-to-set = 1 [ set short-rotation-stand-list lput this-stand short-rotation-stand-list ]
    if mgmt-to-set = 2 [ set long-rotation-stand-list lput this-stand long-rotation-stand-list ]
    if mgmt-to-set = 3 [ set set-aside-stand-list lput this-stand set-aside-stand-list ]
    set remaining-stands remove this-stand remaining-stands ; update remaining stand list
  ]
end

to create-initial-trees
  ; Setup short-rotation stands
  foreach short-rotation-stand-list [ this-stand ->
    let age-to-set age-dist:short-rotation                                    ; chooses an age from short-rotation age dist (so that all trees in the stand have the same age)
    ask ( stand-patches this-stand ) [
      ; set every 4.572 (15 ft) as a tree (every 3 cells)
      if pxcor mod 3 = 0 AND pycor mod 3 = 0                 [ create-tree age-to-set ]                                         ; create a tree
      ; create cells where trees will be (offset placement) after harvest, allowing stumps to persist
      if ( pxcor + 2 ) mod 3 = 0 AND ( pycor + 2 ) mod 3 = 0 [
        if not p_cell [ set p_cell true sprout-cells 1 [ set-default-cell-properties ] ; Create cells at all places where they will be needed
          ask cells-here [
            set potential-tree? true
            if p_edge = 1 [ set edge? true ]
            set stand-id p_stand-id
            if p_mgmt = 1 [ set mgmt "short-rotation" ]
            if p_mgmt = 2 [ set mgmt "long-rotation" ]
            if p_mgmt = 3 [ set mgmt "set-aside" ]
            ; set attr-road p_attr-road ; set road attractiveness based on patch ; check - not anymore?
          ]
        ]
      ]
    ]
  ]
  ; Setup long-rotation stands
  ; ADD CHECK FIX: Attraction if recent thinning and dead trees?
  foreach long-rotation-stand-list [ this-stand ->
    let age-to-set age-dist:long-rotation                       ; chooses age from long-rotation age dist (so all trees in a stand have same age)
    ask ( stand-patches this-stand ) with [ ( pxcor mod 2 = 0 AND pycor mod 2 = 0 ) OR  ( ( pxcor + 1 ) mod 2 = 0 AND ( pycor + 1 ) mod 2 = 0 ) ] [
      if not p_cell [ set p_cell true sprout-cells 1 [ set-default-cell-properties ] ; Create cells at all places where they will be needed
        ask cells-here [
          set potential-tree? true
          if p_edge = 1 [ set edge? true ]
          set stand-id p_stand-id
          if p_mgmt = 1 [ set mgmt "short-rotation" ]
          if p_mgmt = 2 [ set mgmt "long-rotation" ]
          if p_mgmt = 3 [ set mgmt "set-aside" ]
        ]
      ]
    ]
    ; Initialize the stand only where the coordinates meet the following criteria
    ask ( stand-patches this-stand ) with [ pxcor mod 2 = 0 AND pycor mod 2 = 0 ] [             ; 3.048 m / 10 ft / 2 cell spacing
      ; The number of trees to generate depends on the age, which determines how many trees have been thinned and therefore the density of trees
      if ( age-to-set <= 15 ) [ create-tree age-to-set ]                                                          ; if before the first thinning, all patches at the spacing sprout trees
      if ( age-to-set >  15 and age-to-set <= 35 and random-float 1 < 0.688705234 ) [ create-tree age-to-set ]    ; if after the first but before the second thinning, a smaller proportion of patches make trees
      if ( age-to-set >  35 and age-to-set <= 55 and random-float 1 < 0.367309458 ) [ create-tree age-to-set ]    ; if after the second but before the third thinning, a smaller proportion of patches make trees
      if ( age-to-set >  55 and random-float 1 < 0.229568411 ) [ create-tree age-to-set ]                         ; if after the third thinning, a smaller proportion of patches make trees
    ]
  ]
end

; CREATE NEW TREES - Used in "create-initial-trees" and "regenerate-stands" to create new trees - They start at a defined age
to create-tree [ this-tree-age ]
  ; ADD PROB-insect-WOUND
  if ( not p_cell ) [ set p_cell true set p_tree 1 sprout-cells 1 [ set-default-cell-properties ] ] ; Shows this patch has a cell and thus cannot make more
  ask cells-here [
    set tree? true
    set alive? true
    set age this-tree-age
    ifelse age-inf-suscept != "age-indep" [ set age-bsrd-inf-suscept  ( one-of ( item ( this-tree-age - 1 ) age-bsrd-inf-suscept_list ) ) ]   ; assign a susceptibility from the distribution based on age
      [ set age-bsrd-inf-suscept  ( one-of age-bsrd-inf-suscept_list ) ]
    set prob-insect-wound ( one-of prob-insect-wound_live-nonInf_list )          ; Set prob-wound

    set  prob-root-transm ( one-of prob-root-transmission_list )
    set  prob-inf_SC_base (( one-of prob-insect-infested_SC_list ) * ( one-of prob-insect-transmission_SC_list ))
    set  prob-inf_PF_base (( one-of prob-insect-infested_PF_list ) * ( one-of prob-insect-transmission_PF_list ))
    set  prob-inf_HN_base (( one-of prob-insect-infested_HN_list ) * ( one-of prob-insect-transmission_HN_list ))

    ; CHECK: necessary?
    set mort-cause              ""
    set time-since-mort         0


    ; Set stand id and mgmt based on patch
    if stand-id = 0 [ set stand-id p_stand-id ]
    if mgmt = "" [
      if p_mgmt = 1 [ set mgmt "short-rotation" ]
      if p_mgmt = 2 [ set mgmt "long-rotation" ]
      if p_mgmt = 3 [ set mgmt "set-aside" ]
    ]
    set potential-tree? true
    set trees_cumul ( trees_cumul + 1) ; add to cumulative tree count
    if mgmt = "short-rotation" [ set trees_short-rotation_cumul ( trees_short-rotation_cumul + 1) ] ; add to cumulative tree count
    if mgmt = "long-rotation" [ set trees_long-rotation_cumul ( trees_long-rotation_cumul + 1) ] ; add to cumulative tree count
  ]
end

;_____________________________________________________________________________________________________________
; ROAD PROCEDURES

; CREATE ROADS IN THE LANDSCAPE PROPORTIONAL TO THE AREA
to setup-roads
  ; IF THERE IS ONLY ONE STAND, THERE ARE NO EDGES. THEREFORE, CREATE THE APPROPRIATE NUMBER OF ROAD CELLS BASED ON AREA AND DISTRIBUTE RANDOMLY
  if ( config = "random" or config = "clustered" ) and n-stands = 1 [
    if any? patches with [ p_mgmt = 1 or p_mgmt = 2 ] [
      let area-plantations ( count patches * cell-area / 10000 )                                       ; Calculate the area of plantations in ha
      let n-road-cells ( area-plantations * plantation-road-density / intercell-distance )             ; Calculate number of road cells needed given the area of plantations
      ask n-of n-road-cells patches with [ p_tree = 0 ] [ create-road_patch ]                                ; Have that number of patches make roads
    ]
  ]

  ; IF THERE ARE MULTIPLE STANDS, TURN SOME EDGES INTO ROADS AT A PROPORTION BASED ON THE AREA OF PLANTATIONS IN THE LANDSCAPE
  if config = "random" or config = "clustered" and n-stands > 1 [
    ; Assign edge patches & cells to have p_edge = 1 & edge? = true
    ; CHECK: Assign all edge patches to a global variable (to make certain calculations easier)
    set edges ( patches with [ any? neighbors with [ p_stand-id != [ p_stand-id ] of myself ] ] )                ; Edges are patches at the edge of a stand (neighbors have diff stand id)
    ask edges [                                                                                                                 ; Ask edge patches
      set p_edge 1 if visualize? [ set pcolor [ pcolor ] of self + 1.25 ]                                                       ; designate as edges w/ p_edge property
      if not p_cell [ set p_cell true sprout-cells 1 [ set-default-cell-properties ] ]                                          ; create an edge cell to spread attr and make roads (if not already there)
      ask cells-here [ set edge? true set stand-id p_stand-id                                                                   ; assign that cell to have propeties of the patch
        if p_mgmt = 1 [ set mgmt "short-rotation" ] if p_mgmt = 2 [ set mgmt "long-rotation" ] if p_mgmt = 3 [ set mgmt "set-aside" ]
      ]
    ]
    set edges ( edges with [ p_tree = 0 ] ) ; reduce edges to patches that could have roads (i.e. those that don't have trees)
    ; Set roads for plantations; total plantation area = number of plantation patches * patch area
    if ( length short-rotation-stand-list > 0 ) OR ( length long-rotation-stand-list > 0 ) [                    ; if there are any short-rotation or long-rotation stands
      let area-plantations ( count patches * cell-area / 10000 )                                         ; Calculate the area of plantations in ha
      let n-necessary-roads ( area-plantations * plantation-road-density / intercell-distance )               ; Calculate number of road cells needed given the area of plantations
      let n-potential-roads ( count edges with [ p_tree = 0 AND ( p_mgmt = 1 or p_mgmt = 2 ) ] )         ; Calculate the number of cells available to become roads (those without trees)
      let potential-road-deficit_initial ( n-potential-roads - n-necessary-roads )                       ; Calculate the number needed in case there aren't sufficient patches available to make roads
      let potential-road-deficit potential-road-deficit_initial                                          ; Calculate the number needed in case there aren't sufficient patches available to make roads
      if  potential-road-deficit <= 0 [                                                                  ; If there's a road deficit (i.e. there are less edge cells than necessary to be roads):
        print (word "Potential road deficit: " potential-road-deficit ) ; delete
        while [ potential-road-deficit <= 0 ] [                                                            ; While the deficit remains
          ask edges with [ p_mgmt = 1 or p_mgmt = 2 ] [                                                      ; Ask edge patches with plantation management and no trees
            ask neighbors with [ p_tree = 0 ] [ if ( p_mgmt = 1 or p_mgmt = 2 ) [ set p_edge 1 ] ]           ; to make their nbors edges too
          ]
          set edges ( patches with [ p_edge = 1 and p_tree = 0 ] )                                           ; Set edge agentset again
          set n-potential-roads ( count edges with [ not p_cell AND ( p_mgmt = 1 or p_mgmt = 2 ) ] )         ; Recalculate the number of patches available to make roads
          set potential-road-deficit ( n-potential-roads - n-necessary-roads )                               ; Realculate the number needed in case there aren't sufficient patches available to make roads
        ]
        ask n-of ( ( abs potential-road-deficit_initial ) + 1 ) edges with [ not p_cell ] [                    ; Add as many cells as you need to fill the deficit of potential road cells
          set p_cell true sprout-cells 1 [ set-default-cell-properties ]                                           ; make them edges
          ask cells-here [ set edge? true set stand-id p_stand-id
            if p_mgmt = 1 [ set mgmt "short-rotation" ] if p_mgmt = 2 [ set mgmt "long-rotation" ] if p_mgmt = 3 [ set mgmt "set-aside" ]
          ]
        ]
      ]
      ask ( n-of n-necessary-roads cells with [ edge? AND ( mgmt = "short-rotation" or mgmt = "long-rotation" ) and not tree? ] ) [ set road? true ]        ; Have that number of cells become roads
    ]
  ]
  if spread-setup? [
    if any? cells with [ road? ] [ add-attr-road ]
  ]
end

; CREATE ROADS - Used in "setup-roads"
to create-road_patch ; For use by patches only
  if not p_cell [ set p_cell true sprout-cells 1 [ set-default-cell-properties ] ]   ; Create a cell if there's not one
  set p_road 1  ; Designate self as a road & having a cell
  ask cells-here [ set road? true set edge? true ]
end


to remove-edges ; check necessary? should I remove roads too?
  if ( any? cells   with [ edge?  ] ) [ ask cells   [ set edge?  false set road?  false ] ]
  if ( any? patches with [ p_edge ] ) [ ask patches [ set p_edge false set p_road false ] ]
  set edges empty-agentset
end


;_____________________________________________________________________________________________________________
; GENERATE INITIAL INFECTION

to generate-initial-infection-pct
  ; Create an initial infection (based on initial-infection-proportion slider), with trees set to be infected
  ; based on the proportion of trees they account for to avoid always having higher density of initial infections
  ; in the more densely planted stands. With this approach, % initial infected in short-rotation = % initial infected in long-rotation
  carefully [
    ifelse any? trees_current [
      let prop-trees_long-rotation (( count trees_current with [ mgmt = "long-rotation" ] ) / ( count trees_current ))
      let prop-trees_short-rotation (( count trees_current with [ mgmt = "short-rotation" ] ) / ( count trees_current ))
      let n-initial-infections ( count trees_current ) * ( pct-initial-infection / 100 )

      ; Set initial infections in long-rotation stands
      ask n-of ( ceiling ( n-initial-infections * prop-trees_long-rotation ) ) trees_current with [ mgmt = "long-rotation" ] [
        set infected? true set newly-infected? true set inf-initial? true
        set prob-insect-wound ( one-of prob-insect-wound_live-inf_list )
        set age-bsrd-mort-suscept ( one-of ( item ( age - 1 ) age-bsrd-mort-suscept_year-of-infection_list ) )
      ]
      ; Set initial infections in short-rotation stands
      ask n-of ( ceiling ( n-initial-infections * prop-trees_short-rotation ) ) trees_current with [ mgmt = "short-rotation" ] [
        set infected? true set newly-infected? true set inf-initial? true
        set prob-insect-wound ( one-of prob-insect-wound_live-inf_list )
        set age-bsrd-mort-suscept ( one-of ( item ( age - 1 ) age-bsrd-mort-suscept_year-of-infection_list ) )
      ]
    ] [ print "PROBLEM: setup -> generate-initial-infection-pct -- No \"trees_current\"! Trees were either not generated or the agentset was never set." ]
  ] [ print "ERROR: generate-initial-infection-pct" print error-message set kill-run? true set kill-run? true stop ]
end

to generate-initial-infection-alt
  ask patch ( ceiling max-pxcor * ( 1 / 4 ) ) ( ceiling max-pycor * ( 1 / 4 ) ) [
    ask one-of cells in-radius 20 with [ tree? and alive? ] [ set infected? true set newly-infected? true set inf-initial? true set prob-insect-wound ( one-of prob-insect-wound_live-inf_list )
      set age-bsrd-mort-suscept ( one-of ( item ( age - 1 ) age-bsrd-mort-suscept_year-of-infection_list ) )
    ]
  ]
  ask patch ( ceiling max-pxcor * ( 1 / 4 ) ) ( ceiling max-pycor * ( 3 / 4 ) ) [
    ask one-of cells in-radius 20 with [ tree? and alive? ] [ set infected? true set newly-infected? true set inf-initial? true set prob-insect-wound ( one-of prob-insect-wound_live-inf_list )
      set age-bsrd-mort-suscept ( one-of ( item ( age - 1 ) age-bsrd-mort-suscept_year-of-infection_list ) )
    ]
  ]
  ask patch ( ceiling max-pxcor * ( 3 / 4 ) ) ( ceiling max-pycor * ( 1 / 4 ) ) [
    ask one-of cells in-radius 20 with [ tree? and alive? ] [ set infected? true set newly-infected? true set inf-initial? true set prob-insect-wound ( one-of prob-insect-wound_live-inf_list )
      set age-bsrd-mort-suscept ( one-of ( item ( age - 1 ) age-bsrd-mort-suscept_year-of-infection_list ) )
    ]
  ]
  ask patch ( ceiling max-pxcor * ( 3 / 4 ) ) ( ceiling max-pycor * ( 3 / 4 ) ) [
    ask one-of cells in-radius 20 with [ tree? and alive? ] [ set infected? true set newly-infected? true set inf-initial? true set prob-insect-wound ( one-of prob-insect-wound_live-inf_list )
      set age-bsrd-mort-suscept ( one-of ( item ( age - 1 ) age-bsrd-mort-suscept_year-of-infection_list ) )
    ]
  ]
  ;ask trees_infected [ apply attraction ]        ; ADD CHECK
end

to generate-initial-infection-center
  ask patch ( ceiling max-pxcor * ( 1 / 2 ) ) ( ceiling max-pycor * ( 1 / 2 ) ) [
    ask one-of cells in-radius 10 with [ tree? and alive? ] [ set infected? true set newly-infected? true set inf-initial? true set prob-insect-wound ( one-of prob-insect-wound_live-inf_list )
      set age-bsrd-mort-suscept ( one-of ( item ( age - 1 ) age-bsrd-mort-suscept_year-of-infection_list ) )
    ]
  ]
  ;ask trees_infected [ apply attraction ]        ; ADD CHECK
end


to reinitiate-infection ; If infection dies out, start a new infection? CHECK FIX does this make sense
  carefully [
    let drew-age-0? false
    ifelse any? trees_current [
      if not any? trees_current with [ infected? ] [
        ; Create an initial infection (based on initial-infection-proportion slider), with trees set to be infected
        ; based on the proportion of trees they account for to avoid always having higher density of initial infections
        ; in the more densely planted stands. With this approach, % initial infected in short-rotation = % initial infected in long-rotation
        print "infection lost... reinitiating based on \"initial-infection\""

        let prop-trees_long-rotation (( count trees_current with [ mgmt = "long-rotation" ] ) / ( count trees_current ))
        let prop-trees_short-rotation (( count trees_current with [ mgmt = "short-rotation" ] ) / ( count trees_current ))
        let n-initial-infections ( count trees_current ) * ( pct-initial-infection / 100  * 0.1 )

        ; Set initial infections in long-rotation stands
        ask n-of ( ceiling ( n-initial-infections * prop-trees_long-rotation ) ) trees_current with [ mgmt = "long-rotation" ] [
          set infected? true set newly-infected? true set inf-initial? true
          set prob-insect-wound ( one-of prob-insect-wound_live-inf_list )
          ; In case you draw age 0 (shouldn't happen)
          ifelse age > 0 [
            set age-bsrd-mort-suscept ( one-of ( item ( age - 1 ) age-bsrd-mort-suscept_year-of-infection_list ) )
          ]
          [ set age-bsrd-mort-suscept ( one-of ( item ( age ) age-bsrd-mort-suscept_year-of-infection_list ) ) set drew-age-0? true ]
        ]
        ; Set initial infections in short-rotation stands
        ask n-of ( ceiling ( n-initial-infections * prop-trees_short-rotation ) ) trees_current with [ mgmt = "short-rotation" ] [
          set infected? true set newly-infected? true set inf-initial? true
          set prob-insect-wound ( one-of prob-insect-wound_live-inf_list )
          ifelse age > 0 [
            set age-bsrd-mort-suscept ( one-of ( item ( age - 1 ) age-bsrd-mort-suscept_year-of-infection_list ) )
          ]
          [ set age-bsrd-mort-suscept ( one-of ( item ( age ) age-bsrd-mort-suscept_year-of-infection_list ) ) set drew-age-0? true ]
        ]
      ]
    ][ print "No trees_current. Could not reinitiate infection" ]
    if drew-age-0? [ print "PROBLEM: tried to assign age-0 agent in trees_current to be an infection. Likely not a tree." ]
  ][ print "ERROR go -> reinitiate infection" print error-message set kill-run? true set kill-run? true stop ]
end

;_____________________________________________________________________________________________________________
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
;_____________________________________________________________________________________________________________
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
;_____________________________________________________________________________________________________________
; GO PROCEDURE
;_____________________________________________________________________________________________________________
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
;_____________________________________________________________________________________________________________
to go ; PROCEDURE THAT RUNS EVERY TICK
  ;_____________________________________________________________________________________
  ; GO SETUP
  ;_____________________________________________________________________________________
  ; Stop the model once you reach the maximum number of ticks (failsafe)
  no-display
  if ( limit-ticks? AND ( ticks > max-ticks ) ) [
    if run-file != 0 and run-file != "" [
      file-open  run-file
      file-print (word max-ticks " time steps reached!")
      file-print ( word "End: " date-and-time )
      file-close
    ]
    print (word max-ticks " time steps reached!")
    print ( word "End: " date-and-time )
    set kill-run? true stop
  ]

  if do-kill-runs? [
    if kill-run? [
      print "///////////////////////////////////////////////////////////////////////////////////////////////"
      print "------------------------ FATAL ERROR: kill-run triggered. Run aborted. ------------------------"
      print "///////////////////////////////////////////////////////////////////////////////////////////////"
      stop
    ]
  ]

  reset-timer ; for keeping track of time ; delete

  print (word "\n---------------- Tick " ticks " ----------------\n" ) ; for easier visual parsing of output

  ;_____________________________________________________________________________________

  ; SPREAD INFECTION - Calculate infection probabilities and have cells become infected
  carefully [
    if ( spread-infection-version? != "none" ) [
      project-infection
    ]
  ] [ print "ERROR go -> project-infection" print error-message set export-rasters? false set export-tree-data? false set kill-run? true stop ]

  ;_____________________________________________________________________________________

  ; STAND MANAGEMENT PROCEDURES
  carefully [
    manage-stands
  ] [ print "ERROR go -> manage-stands" print error-message set export-rasters? false set export-tree-data? false set kill-run? true stop  ]


  ;_____________________________________________________________________________________

  ; TIME-PASSES: Various time-related model processes - trees age, gain or lose infection susceptibility, die, old dead trees and stumps disappear
  carefully [
    time-passes
  ] [ print "ERROR go -> time-passes" print error-message set export-rasters? false set export-tree-data? false set kill-run? true stop ]
  ;_____________________________________________________________________________________

  ; EXPORT DATA AND PRINT TIMER
  ; Export data: This happens after time passes to account for new deaths
  ; Progress file tracks the number of current and cumulative infections and mortalities at each tick
  carefully [
    if track-run? [
      ; ensure the progress file exists and is written to
      if progress-file != 0 and progress-file != "" [ write-to-progress-file ]
    ]
  ] [ print "ERROR go -> writing progress file" print error-message set kill-run? true stop ]

  ; Data file tracks the states of each tree at the frequency determined by the user
  carefully [
    if ( export-tree-data? and ticks mod tree-data-export-freq = 0 ) [ export-tree-data ]
  ] [ print "ERROR go -> setting export-tree-data" print error-message set kill-run? true stop ]
  carefully [
    if ( export-stand-data? and ticks mod stand-data-export-freq = 0 ) [ export-stand-data ]
  ] [ print "ERROR go -> setting export-stand-data" print error-message set kill-run? true stop ]

  ; Raster files track the states of each cell tree at the frequency determined by the user
  carefully [
    if ( export-rasters? AND ( ticks mod raster-export-frequency = 0 ) ) [ export-inf-map ]
  ] [ print "ERROR go -> setting trees_current" print error-message set kill-run? true stop ]

  if visualize? [
    color-world display
    ask newly-infected-trees [ set color magenta - 1 ]
    ask losing-infection-trees [ set color violet + 1 ]
    ask trees_potential with [ new-bsrd-mortality? ] [ set color orange ]
    ask trees_potential with [ losing-attr-dead? ] [ set color yellow ]
  ]

  let go-duration timer
  print (word "Tick " ticks " took: " go-duration " seconds." )
  if track-run? [
    file-open run-file
    file-print (word "Tick " ticks " took: " go-duration " seconds." )
    file-close-all
  ]





;  ; test whether project attraction worked
;  print "_________________________________________"
;  print "did project attraction work?"
;  print "---------------"

;  print (word "are there any trees with attr-inf? " ( any? trees with [ length attr-inf > 1 ] ) )
;  print (word "are there any infected trees? " any? infected-trees )
;  print (word "are there any trees losing infection? " any? losing-infection-trees )
;  print "---------------"
;  print (word "are there any trees with attr-dead? " ( any? trees with [ length attr-dead > 1 ] ) )
;  print (word "are there any bsrd-killed trees? " ( any? trees with [ not alive? and mort-cause = "bsrd" ] ) )
;  print (word "are there any trees losing attr-dead? " ( any? trees with [ losing-attr-dead? ] ) )
;  print "_________________________________________"
;  print "_________________________________________"

  ;____________________________________________________________________________
  ; FINAL ACTIONS:
  ; Set trees_current
  carefully [
    set trees_current ( trees )   ; set an agentset for all current trees
  ] [ print "ERROR go -> setting trees_current" print error-message set kill-run? true stop ]

  ; If infection lost, reinitiate
  if reinitiate-infection? [ reinitiate-infection ] ; ADD IN RANDOM INFECTION IF THERE ARE NONE - CHECK FIX does this make sense?

  ; Increase ticks (to next time step)
  tick

  if limit-ticks? and ticks = max-ticks [ delete-program-directory ]
end



;////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
; PROJECT INFECTION - Core of the infection spread for the model
;////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

to project-infection
  ; SET AGENTSETS FOR TREES NEWLY-INFECTED and LOSING-INFECTION (the cells projecting infection)
  carefully [
    ;print ( word "Finding trees spreading infection..." )
    set trees_newly-infected ( newly-infected-trees )
    print ( word "Number of trees_newly-infected: " ( count trees_newly-infected ) )
    set trees_losing-infection ( losing-infection-trees )
    print ( word "Number of trees_losing-infection: " ( count trees_losing-infection ) )
  ] [ print "Error project-infection -> set tree agentsets" print error-message set kill-run? true stop ]

  ;_______________________________________________________________
  ; PROJECT PROBABILITY OF INSECT DISPERSAL PROCEDURE
  ; if insect dispersal
  if ( ( max-insect-distance_SC_m + max-insect-distance_PF_m + max-insect-distance_HN_m ) > 0 ) [ ; Check that insect dispersal is active
    ; output the trees_newly-infected and trees_losing-infection to JSON
    export-newly-infected-trees-json
    export-losing-infection-trees-json
    ;print "Exporting newly-infected and losing-infection files: complete"

    ; run the program and load the data back into NetLogo
    run-spread-infection_Go
    ;print "Running spread-infection program: complete"

    ; use the data to assign probabilities
    ifelse not probability-correction? [ apply-solution_Go ] [ apply-solution_Go_correction ] ; DELETE once tested
    ;print "Applying solution from spread-infection: complete"

    ; Background probability of dispersal
    ; ADD: (4) Background projection
    ; FIX: CHECK: This way of doing things results in trees within NO DISPERSAL RADIUS having prob-insect-disp > than some WITHIN a dispersal radius
    if background-prob-insect-dispersal? [
      calculate-insect-dispersal_background
    ]
    ; FIX HERE
  ]

  ;_______________________________________________________________
  ; PROJECT ROOT CONTACT: PROJECT TREE AGES AND IDS

  ; Difference between "proportional" attraction vs. "median" vs. "max":
  ;   - "proportional" approach (**used for thesis runs**) uses
  ;   - "median" and "max" appraches use the respective calue

  ifelse new-attr? [
    ifelse attr-inf-dead_version = "proportional" [
      project_root-infection_attraction-proportional
    ]
    [ project_root-infection_attraction ]
  ]
  [ project-root-infection ]

  ;_______________________________________________________________

  ; Clear agentsets
  set trees_newly-infected     empty-agentset
  set trees_losing-infection   empty-agentset
end

;_______________________________________________________________________________________________
; PROJECT ROOT INFECTION PROBABILITIES
to project_root-infection_attraction
  ; trees spreading root infection (newly-infected and losing-infection) defined previously in main project-infection procedure

  ; Newly-infected add their who and age to their neighbors' properties
  let attr-distance-range ( range 1 ( max-attraction-distance_cell + 1 ) )               ; create a list for the range of distances
;  print (word "attr-distance-range " attr-distance-range )

  ; ID all of the trees that are projecting infection and/or attraction
  ifelse attr-dead? [
    set cells-projecting ( turtle-set
      trees_current with   [ newly-infected? ]               ; a) newly infected:     add to prob root contact (+) and infection attraction (+)
      trees_potential with [ losing-infection? ]             ; b) losing infection:   remove prob root contact (-) AND infection attraction (-)
      trees_potential with [ new-bsrd-mortality? ]           ; c) newly dead (from BSRD):                 add to dead attraction (attr-dead) (+)
      trees_potential with [ losing-attr-dead? ]             ; d) dead losing attraction and being reset: remove dead attraction (attr-dead) (-)
    )
  ]
  [   ; NO DEAD ATTRACTION
    set cells-projecting ( turtle-set
      trees_current with   [ newly-infected? ]               ; a) newly infected:     add to prob root contact (+) and infection attraction (+)
      trees_potential with [ losing-infection? ]             ; b) losing infection:   remove prob root contact (-) AND infection attraction (-)
    )
  ]

  ; If there are any trees spreading root infection
  carefully [
    ifelse is-agentset? cells-projecting and any? cells-projecting [
      let source-is-newly-infected?      false
      let source-is-losing-infection?    false
      let source-is-new-bsrd-mortality?  false
      let source-is-losing-attr-dead?    false
      ask cells-projecting [
        ; Set temp variables to track source tree properties
        let age_source-tree age
        let id_source-tree who
        ifelse newly-infected?     [ set source-is-newly-infected?        true ] [ set source-is-newly-infected?        false ]
        ifelse losing-infection?   [ set source-is-losing-infection?      true ] [ set source-is-losing-infection?      false ]
        ifelse new-bsrd-mortality? [ set source-is-new-bsrd-mortality?    true ] [ set source-is-new-bsrd-mortality?    false ]
        ifelse losing-attr-dead?   [ set source-is-losing-attr-dead?      true ] [ set source-is-losing-attr-dead?      false ]

        ; Root spread with distance based on management and stand position
        let max-root-contact-distance_source 0 ; the number of cells the tree can have root contact based on management spacing and whether an edge
        ifelse pct-short-rotation > 0 [
          ; if there are short-rotation stands, allow long-rotation on their borders to have root contact
          ifelse ( ( edge? and mgmt = "long-rotation" ) or mgmt = "short-rotation" ) [
            set max-root-contact-distance_source max-root-distance_short-rotation_cell   ; short-rotation and long-rotation edges
          ]
          [ set max-root-contact-distance_source max-root-distance_long-rotation_cell ] ; long-rotation interior
        ]
        [ set max-root-contact-distance_source max-root-distance_long-rotation_cell ]   ; long-rotation

        ; Loop through root contact distances
        foreach attr-distance-range [ current-ring-distance ->
          ; Calculate the current distance (in meters)
          let current-distance_m ( current-ring-distance * intercell-distance )

          ; for each potential tree in each ring
          ask cells at-points ( item ( current-ring-distance - 1 ) ring-list_attraction ) [
            if potential-tree? [
              ;------------------------------------------------------------------------------------------
              ; PROJECT ROOT CONTACT INFECTIONS AND INF-ATTR - if the tree is newly infected
              if source-is-newly-infected? [
                ; PROJECT ROOT CONTACT (+) (add) - if within the distance for root contact
                if max-root-distance_m > 0 AND current-ring-distance <= max-root-contact-distance_source [
                  set any-inf-root-nbors? true    ; track that there is some root infection probability that needs to be accounted for
                  carefully [
                    ; add the age and who of each source tree to the relevant distance ring of all neighbors
                    if current-ring-distance = 1 [
                      set inf-root-nbors_ages_r1 ( lput age_source-tree inf-root-nbors_ages_r1 )
                      set inf-root-nbors_ids_r1  ( lput id_source-tree inf-root-nbors_ids_r1 )
                    ]
                    if current-ring-distance = 2  [
                      set inf-root-nbors_ages_r2 ( lput age_source-tree inf-root-nbors_ages_r2 )
                      set inf-root-nbors_ids_r2  ( lput id_source-tree inf-root-nbors_ids_r2 )
                    ]
                    if current-ring-distance = 3 [
                      set inf-root-nbors_ages_r3 ( lput age_source-tree inf-root-nbors_ages_r3 )
                      set inf-root-nbors_ids_r3  ( lput id_source-tree inf-root-nbors_ids_r3 )
                    ]
                    if current-ring-distance = 4 [
                      set inf-root-nbors_ages_r4 ( lput age_source-tree inf-root-nbors_ages_r4 )
                      set inf-root-nbors_ids_r4  ( lput id_source-tree inf-root-nbors_ids_r4 )
                    ]
                  ] [ print "ERROR project-infection -> project_root-infection_attraction -> trees_newly-infected -> project root infection" print error-message set kill-run? true stop  ]
                ]
                ;----------------------------------------------------------------------------------------
                ; PROJECT ATTR-INF (+) (add) - attraction due to infection
                if attr-inf? [
                  ; (i) Draw an attraction value from the list and run through distance decay equation,
                  carefully [
                    let attr-inf_to-set 0   ; to store value
                    ; Draw infection attraction based on whether the tree is alive or dead (potential trees that are not trees treated as alive because they will be regenerated)
                    ; account for distance decay
                    ifelse (not tree? or alive? ) [ set attr-inf_to-set ( attraction-distance-decay_FUNCTION ( one-of attr-effect_inf-vs-live_list ) ( current-distance_m ) ) ]
                    [ set attr-inf_to-set ( attraction-distance-decay_FUNCTION ( one-of attr-effect_inf-vs-dead_list ) ( current-distance_m ) ) ]
                    ; (ii) Add to the attr-inf list for the target tree/cell
                    set attr-inf        ( lput attr-inf_to-set attr-inf )          ; set target cell inf attraction by adding to end of list
                    set attr-inf_source ( lput id_source-tree attr-inf_source )    ; set the mgmt attraction source stand to their own stand id
                  ] [ print "ERROR project-infection -> project_root-infection_attraction -> trees_newly-infected -> project attraction (infection)" print error-message set kill-run? true stop ]
                ]
              ]
              ;------------------------------------------------------------------------------------------
              ; PROJECT INFECTION LOSS AND REMOVE INFECTION ATTRACTION - if source is losing-infection
              if source-is-losing-infection? [
                carefully [
                  ; remove the age and who of each source tree from the relevant distance ring of all potential tree neighbors
                  ; PROJECT ROOT INFECTION (-) (remove) - if within the distance for root contact
                  if max-root-distance_m > 0 AND current-ring-distance <= max-root-contact-distance_source [
                    if current-ring-distance = 1 [
                      ; print ( word "entered ring " current-ring-distance )
                      let value-position ( position id_source-tree inf-root-nbors_ids_r1 )                  ; Find the position of the tree losing infection in the root nbor lists
                      ; print ( word "value position? " value-position )
                      if ( value-position != false ) [
                        ; print "entered remove loop"
                        set inf-root-nbors_ages_r1 ( remove-item value-position inf-root-nbors_ages_r1 ) ; remove the age at that position
                        set inf-root-nbors_ids_r1  ( remove id_source-tree inf-root-nbors_ids_r1 )       ; remove the who of that tree losing infection
                      ]
                    ]
                    if current-ring-distance = 2 [
                      ; print ( word "entered ring " current-ring-distance )
                      let value-position ( position id_source-tree inf-root-nbors_ids_r2 )                  ; Find the position of the tree losing infection in the root nbor lists
                      ; print ( word "value posit? " value-position )
                      if ( value-position != false ) [
                        ; print "entered remove loop"
                        set inf-root-nbors_ages_r2 ( remove-item value-position inf-root-nbors_ages_r2 ) ; remove the age at that position
                        set inf-root-nbors_ids_r2  ( remove id_source-tree inf-root-nbors_ids_r2 )       ; remove the who of that tree losing infection
                      ]
                    ]
                    if current-ring-distance = 3 [
                      ; print ( word "entered ring " current-ring-distance )
                      let value-position ( position id_source-tree inf-root-nbors_ids_r3 )                  ; Find the position of the tree losing infection in the root nbor lists
                      ; print ( word "value posit? " value-position )
                      if ( value-position != false ) [
                        ; print "entered remove loop"
                        set inf-root-nbors_ages_r3 ( remove-item value-position inf-root-nbors_ages_r3 ) ; remove the age at that position
                        set inf-root-nbors_ids_r3  ( remove id_source-tree inf-root-nbors_ids_r3 )       ; remove the who of that tree losing infection
                      ]
                    ]
                    if current-ring-distance = 4 [
                      ; print ( word "entered ring " current-ring-distance )
                      let value-position ( position id_source-tree inf-root-nbors_ids_r4 )                  ; Find the position of the tree losing infection in the root nbor lists
                      ; print ( word "value posit? " value-position )
                      if ( value-position != false ) [
                        ; print "entered remove loop"
                        set inf-root-nbors_ages_r4 ( remove-item value-position inf-root-nbors_ages_r4 ) ; remove the age at that position (can't use "remove" function because it would remove all)
                        set inf-root-nbors_ids_r4  ( remove id_source-tree inf-root-nbors_ids_r4 )       ; remove the who of that tree losing infection
                      ]
                    ]
                  ]
                ] [ print "ERROR project-infection -> project_root-infection_attraction -> trees_losing-infection -> project root infection" print error-message set kill-run? true stop ]
                ;----------------------------------------------------------------------------------------
                ; PROJECT INFECTION ATTRACTION (-) (remove) - infection due to attraction
                if attr-inf? [
                  carefully [
                    let value-position_attr-inf ( position id_source-tree attr-inf_source )              ; Find the position of the tree losing infection in the attr-inf list
                    ifelse ( value-position_attr-inf != false ) [                                        ; If the source id is found in the list
                      set attr-inf        ( remove-item value-position_attr-inf attr-inf        )        ; remove the age at that position
                      set attr-inf_source ( remove-item value-position_attr-inf attr-inf_source )        ; remove the who of that tree losing infection
                    ] [ print (word "PROBLEM in project_root-infection_attraction -> trees_losing-infection -> project attraction (infection): for target tree **" who "**, value of source **" id_source-tree "** not found" ) ]
                  ] [ print "ERROR project-infection -> project_root-infection_attraction -> trees_losing-infection -> project attraction (infection)" print error-message set kill-run? true stop ]
                ]
              ]
              ;------------------------------------------------------------------------------------------
              ; PROJECT DEAD ATTRACTION (+) (add)    - if tree is new bsrd-caused mortality (thin- and harvest-killed trees dealt with separately, at the stand-scale)
              if attr-dead? [
                if source-is-new-bsrd-mortality? [
                  carefully [
                    ; draw an attraction value from the list and run through distance decay equation,
                    let attr-dead_to-set ( attraction-distance-decay_FUNCTION ( one-of attr-effect_dead_list ) ( current-distance_m ) )
                    set attr-dead        ( lput attr-dead_to-set attr-dead )        ; set target cell dead attraction by adding to the list
                    set attr-dead_source ( lput id_source-tree attr-dead_source )   ; set the mgmt attraction source stand to their own stand id
                  ] [ print "ERROR project-infection -> project_root-infection_attraction -> new-bsrd-mortality -> project attraction (dead)" print error-message set kill-run? true stop ]
                ]
                ;------------------------------------------------------------------------------------------
                ; PROJECT DEAD ATTRACTION (-) (remove) - if dead tree is losing attraction and being reset
                if source-is-losing-attr-dead? [
                  carefully [
                    let value-position_attr-dead ( position id_source-tree attr-dead_source )               ; Find the position of the tree losing infection in the attr-dead list
                    ; print ( word "source: " id_source-tree ) print (word "target: " print who)
                    if ( value-position_attr-dead != false ) [                                              ; If the source id is found in the list
                      set attr-dead        ( remove-item value-position_attr-dead attr-dead )               ; remove the age at that position
                      set attr-dead_source ( remove-item value-position_attr-dead attr-dead_source )        ; remove the who of that tree losing infection
                    ]
                  ] [ print "ERROR project-infection -> project_root-infection_attraction -> losing-attr-dead -> project attraction (infection)" print error-message ]
                ]
              ]
            ]
          ]
        ]
      ]
    ] [ print "------------------>" print "NO CELLS PROJECTING" print "------------------>" ]
  ] [ print "Error in project-infection -> project-root-infection -> ask cells-projecting" print error-message set kill-run? true stop ]

  ask cells-projecting [
    ; RESET PROPERTIES OF PROJECTORS
    set newly-infected?      false              ; remove newly-infected status
    set losing-infection?    false              ; remove losing-infection status
    set new-bsrd-mortality?  false              ; remove newly-dead status
    set losing-attr-dead?    false              ; remove losing-dead-attr status
  ]

  ifelse any? ( cells-projecting with [ newly-infected? OR losing-infection? OR new-bsrd-mortality? OR losing-attr-dead? ] ) [
    print "PROBLEM !!! cell status used for projection WAS NOT CLEARED!!!!"
    user-message "PROBLEM !!! cell status used for projection WAS NOT CLEARED!!!!"
  ]
  [ print "All projecting properties correctly cleared" ]
  ; Clear the cells-projecting agentset
  set cells-projecting empty-agentset
end

to project_root-infection_attraction-proportional ; NEWEST. USED FOR THESIS RUNS.
  ; trees spreading root infection (newly-infected and losing-infection) defined previously in main project-infection procedure

  ; Newly-infected add their who and age to their neighbors' properties
  let attr-distance-range ( range 1 ( max-attraction-distance_cell + 1 ) )               ; create a list for the range of distances
;  print (word "attr-distance-range " attr-distance-range )

  ; ID all of the trees that are projecting infection and/or attraction
  ifelse attr-dead? [
    set cells-projecting ( turtle-set
      trees_current   with [ newly-infected? ]               ; a) newly infected:     add to prob root contact (+) and infection attraction (+)
      trees_potential with [ losing-infection? ]             ; b) losing infection:   remove prob root contact (-) AND infection attraction (-)
      trees_potential with [ new-bsrd-mortality? ]           ; c) newly dead (from BSRD):                 add to dead attraction (attr-dead) (+)
      trees_potential with [ losing-attr-dead? ]             ; d) dead losing attraction and being reset: remove dead attraction (attr-dead) (-)
    )
  ]
  [   ; NO DEAD ATTRACTION
    set cells-projecting ( turtle-set
      trees_current with   [ newly-infected? ]               ; a) newly infected:     add to prob root contact (+) and infection attraction (+)
      trees_potential with [ losing-infection? ]             ; b) losing infection:   remove prob root contact (-) AND infection attraction (-)
    )
  ]

  ; If there are any trees spreading root infection
  carefully [
    ifelse is-agentset? cells-projecting and any? cells-projecting [
      let source-is-newly-infected?      false
      let source-is-losing-infection?    false
      let source-is-new-bsrd-mortality?  false
      let source-is-losing-attr-dead?    false
      ask cells-projecting [
        ; Set temp variables to track source tree properties
        let age_source-tree age
        let id_source-tree who
        ifelse newly-infected?     [ set source-is-newly-infected?        true ] [ set source-is-newly-infected?        false ]
        ifelse losing-infection?   [ set source-is-losing-infection?      true ] [ set source-is-losing-infection?      false ]
        ifelse new-bsrd-mortality? [ set source-is-new-bsrd-mortality?    true ] [ set source-is-new-bsrd-mortality?    false ]
        ifelse losing-attr-dead?   [ set source-is-losing-attr-dead?      true ] [ set source-is-losing-attr-dead?      false ]

        ; Root spread with distance based on management and stand position
        let max-root-contact-distance_source 0 ; the number of cells the tree can have root contact based on management spacing and whether an edge
        ifelse pct-short-rotation > 0 [
          ; if there are short-rotation stands, allow long-rotation on their borders to have root contact
          ifelse ( ( edge? and mgmt = "long-rotation" ) or mgmt = "short-rotation" ) [
            set max-root-contact-distance_source max-root-distance_short-rotation_cell   ; short-rotation and long-rotation edges
          ]
          [ set max-root-contact-distance_source max-root-distance_long-rotation_cell ] ; long-rotation interior
        ]
        [ set max-root-contact-distance_source max-root-distance_long-rotation_cell ]   ; long-rotation

        ; Loop through root contact distances
        foreach attr-distance-range [ current-ring-distance ->
          ; Calculate the current distance (in meters)
          let current-distance_m ( current-ring-distance * intercell-distance )

          ; for each potential tree in each ring
          ask cells at-points ( item ( current-ring-distance - 1 ) ring-list_attraction ) [
            if potential-tree? [
              ;------------------------------------------------------------------------------------------
              ; PROJECT ROOT CONTACT INFECTIONS AND INF-ATTR - if the tree is newly infected
              if source-is-newly-infected? [
                ; PROJECT ROOT CONTACT (+) (add) - if within the distance for root contact
                if max-root-distance_m > 0 AND current-ring-distance <= max-root-contact-distance_source [
                  set any-inf-root-nbors? true    ; track that there is some root infection probability that needs to be accounted for
                  carefully [
                    ; add the age and who of each source tree to the relevant distance ring of all neighbors
                    if current-ring-distance = 1 [
                      set inf-root-nbors_ages_r1 ( lput age_source-tree inf-root-nbors_ages_r1 )
                      set inf-root-nbors_ids_r1  ( lput id_source-tree inf-root-nbors_ids_r1 )
                    ]
                    if current-ring-distance = 2  [
                      set inf-root-nbors_ages_r2 ( lput age_source-tree inf-root-nbors_ages_r2 )
                      set inf-root-nbors_ids_r2  ( lput id_source-tree inf-root-nbors_ids_r2 )
                    ]
                    if current-ring-distance = 3 [
                      set inf-root-nbors_ages_r3 ( lput age_source-tree inf-root-nbors_ages_r3 )
                      set inf-root-nbors_ids_r3  ( lput id_source-tree inf-root-nbors_ids_r3 )
                    ]
                    if current-ring-distance = 4 [
                      set inf-root-nbors_ages_r4 ( lput age_source-tree inf-root-nbors_ages_r4 )
                      set inf-root-nbors_ids_r4  ( lput id_source-tree inf-root-nbors_ids_r4 )
                    ]
                  ] [ print "ERROR project-infection -> project_root-infection_attraction-proportional -> trees_newly-infected -> project root infection" print error-message set kill-run? true stop ]
                ]
                ;----------------------------------------------------------------------------------------
                ; PROJECT ATTR-INF (+) (add) - attraction due to infection
                if attr-inf? [   ; If attraction to infected trees is activated...
                  carefully [
                    ; Draw an weighted tree value from the list based on distance and add to the attr-inf_source property
                    set attr-inf_source ( attr-inf_source + ( item ( current-ring-distance - 1 ) attr-effect_tree-weight-over-distance_list ) )
                  ] [ print "ERROR project-infection -> project_root-infection_attraction-proportional -> trees_newly-infected -> project attraction (infection)" print error-message set kill-run? true stop ]
                ]
              ]
              ;------------------------------------------------------------------------------------------
              ; PROJECT INFECTION LOSS AND REMOVE INFECTION ATTRACTION - if source is losing-infection
              if source-is-losing-infection? [
                carefully [
                  ; remove the age and who of each source tree from the relevant distance ring of all potential tree neighbors
                  ; PROJECT ROOT INFECTION (-) (remove) - if within the distance for root contact
                  if max-root-distance_m > 0 AND current-ring-distance <= max-root-contact-distance_source [
                    if current-ring-distance = 1 [
                      ; print ( word "entered ring " current-ring-distance )
                      let value-position ( position id_source-tree inf-root-nbors_ids_r1 )                  ; Find the position of the tree losing infection in the root nbor lists
                      ; print ( word "value position? " value-position )
                      if ( value-position != false ) [
                        ; print "entered remove loop"
                        set inf-root-nbors_ages_r1 ( remove-item value-position inf-root-nbors_ages_r1 ) ; remove the age at that position
                        set inf-root-nbors_ids_r1  ( remove id_source-tree inf-root-nbors_ids_r1 )       ; remove the who of that tree losing infection
                      ]
                    ]
                    if current-ring-distance = 2 [
                      ; print ( word "entered ring " current-ring-distance )
                      let value-position ( position id_source-tree inf-root-nbors_ids_r2 )                  ; Find the position of the tree losing infection in the root nbor lists
                      ; print ( word "value posit? " value-position )
                      if ( value-position != false ) [
                        ; print "entered remove loop"
                        set inf-root-nbors_ages_r2 ( remove-item value-position inf-root-nbors_ages_r2 ) ; remove the age at that position
                        set inf-root-nbors_ids_r2  ( remove id_source-tree inf-root-nbors_ids_r2 )       ; remove the who of that tree losing infection
                      ]
                    ]
                    if current-ring-distance = 3 [
                      ; print ( word "entered ring " current-ring-distance )
                      let value-position ( position id_source-tree inf-root-nbors_ids_r3 )                  ; Find the position of the tree losing infection in the root nbor lists
                      ; print ( word "value posit? " value-position )
                      if ( value-position != false ) [
                        ; print "entered remove loop"
                        set inf-root-nbors_ages_r3 ( remove-item value-position inf-root-nbors_ages_r3 ) ; remove the age at that position
                        set inf-root-nbors_ids_r3  ( remove id_source-tree inf-root-nbors_ids_r3 )       ; remove the who of that tree losing infection
                      ]
                    ]
                    if current-ring-distance = 4 [
                      ; print ( word "entered ring " current-ring-distance )
                      let value-position ( position id_source-tree inf-root-nbors_ids_r4 )                  ; Find the position of the tree losing infection in the root nbor lists
                      ; print ( word "value posit? " value-position )
                      if ( value-position != false ) [
                        ; print "entered remove loop"
                        set inf-root-nbors_ages_r4 ( remove-item value-position inf-root-nbors_ages_r4 ) ; remove the age at that position (can't use "remove" function because it would remove all)
                        set inf-root-nbors_ids_r4  ( remove id_source-tree inf-root-nbors_ids_r4 )       ; remove the who of that tree losing infection
                      ]
                    ]
                  ]
                ] [ print "ERROR project-infection -> project_root-infection_attraction-proportional -> trees_losing-infection -> project root infection" print error-message set kill-run? true stop ]
                ;----------------------------------------------------------------------------------------
                ; PROJECT INFECTION ATTRACTION (-) (remove) - infection due to attraction
                if attr-inf? [
                  carefully [
                    ; Draw an weighted tree value from the list based on distance and subtract from attr-inf_source
                    set attr-inf_source ( attr-inf_source - ( item ( current-ring-distance - 1 ) attr-effect_tree-weight-over-distance_list ) )
                  ] [ print "ERROR project-infection -> project_root-infection_attraction-proportional -> trees_losing-infection -> project attraction (infection)" print error-message set kill-run? true stop ]
                ]
              ]
              ;------------------------------------------------------------------------------------------
              ; PROJECT DEAD ATTRACTION (+) (add)    - if tree is new bsrd-caused mortality (thin- and harvest-killed trees dealt with separately, at the stand-scale)
              if attr-dead? [
                if source-is-new-bsrd-mortality? [
                  carefully [
                    ; draw an attraction value from the list and run through distance decay equation,
                    let attr-dead_to-set ( attraction-distance-decay_FUNCTION ( one-of attr-effect_dead_list ) ( current-distance_m ) )
                    set attr-dead        ( lput attr-dead_to-set attr-dead )        ; set target cell dead attraction by adding to the list
                    set attr-dead_source ( lput id_source-tree attr-dead_source )   ; set the mgmt attraction source stand to their own stand id
                  ] [ print "ERROR project-infection -> project_root-infection_attraction-proportional -> new-bsrd-mortality -> project attraction (dead)" print error-message set kill-run? true stop ]
                ]
                ;------------------------------------------------------------------------------------------
                ; PROJECT DEAD ATTRACTION (-) (remove) - if dead tree is losing attraction and being reset
                if source-is-losing-attr-dead? [
                  carefully [
                    let value-position_attr-dead ( position id_source-tree attr-dead_source )               ; Find the position of the tree losing infection in the attr-dead list
                    ; print ( word "source: " id_source-tree ) print (word "target: " print who)
                    if ( value-position_attr-dead != false ) [                                              ; If the source id is found in the list
                      set attr-dead        ( remove-item value-position_attr-dead attr-dead )               ; remove the age at that position
                      set attr-dead_source ( remove-item value-position_attr-dead attr-dead_source )        ; remove the who of that tree losing infection
                    ]
                  ] [ print "ERROR project-infection -> project_root-infection_attraction-proportional -> losing-attr-dead -> project attraction (infection)" print error-message ]
                ]
              ]
            ]
          ]
        ]
      ]
    ] [ print "------------------>" print "NO CELLS PROJECTING" print "------------------>" ]
  ] [ print "Error in project-infection -> project-root-infection -> ask cells-projecting" print error-message set kill-run? true stop ]

  ask cells-projecting [
    ; RESET PROPERTIES OF PROJECTORS
    set newly-infected?      false              ; remove newly-infected status
    set losing-infection?    false              ; remove losing-infection status
    set new-bsrd-mortality?  false              ; remove newly-dead status
    set losing-attr-dead?    false              ; remove losing-dead-attr status
  ]

  ifelse any? ( cells-projecting with [ newly-infected? OR losing-infection? OR new-bsrd-mortality? OR losing-attr-dead? ] ) [
    print "PROBLEM !!! cell status used for projection WAS NOT CLEARED!!!!"
  ]
  [ print "All projecting properties correctly cleared" ]
  ; Clear the cells-projecting agentset
  set cells-projecting empty-agentset
end

to project-root-infection
  ; trees spreading root infection (newly-infected and losing-infection) defined previously in main project-infection procedure

  ; Newly-infected add their who and age to their neighbors' properties
  let root-contact-distance-range 0
  let distance-range_int ( range 1 ( 4 + 1 ) )               ; create a list for the range of distances
  let distance-range_ext ( range 1 ( 3 + 1 ) )               ; create a list for the range of distances (tight spacing prevents root contact beyond 3 cells)

  ; If there are any trees spreading root infection
  carefully [
    if is-agentset? trees_newly-infected and any? trees_newly-infected [
      ask trees_newly-infected [
        let age_source-tree age
        let id_source-tree who
        ; Root spread distance based on management
        if edge? or mgmt = "short-rotation"      [ set root-contact-distance-range distance-range_int ] ; short-rotation and edges
        if mgmt = "long-rotation" and not edge? [ set root-contact-distance-range distance-range_ext ]
        foreach root-contact-distance-range [ current-ring-distance ->
          ask cells at-points ( item ( current-ring-distance - 1 ) ring-list_root ) [
            if potential-tree? [
              set any-inf-root-nbors? true    ; just acknowledging that there is some root infection probability that needs to be accounted for
              ; add the age and who of each source tree to the relevant distance ring of all neighbors
              if current-ring-distance = 1 [
                set inf-root-nbors_ages_r1 ( lput age_source-tree inf-root-nbors_ages_r1 )
                set inf-root-nbors_ids_r1  ( lput id_source-tree inf-root-nbors_ids_r1 )
              ]
              if current-ring-distance = 2  [
                set inf-root-nbors_ages_r2 ( lput age_source-tree inf-root-nbors_ages_r2 )
                set inf-root-nbors_ids_r2  ( lput id_source-tree inf-root-nbors_ids_r2 )
              ]
              if current-ring-distance = 3 [
                set inf-root-nbors_ages_r3 ( lput age_source-tree inf-root-nbors_ages_r3 )
                set inf-root-nbors_ids_r3  ( lput id_source-tree inf-root-nbors_ids_r3 )
              ]
              if current-ring-distance = 4 [
                set inf-root-nbors_ages_r4 ( lput age_source-tree inf-root-nbors_ages_r4 )
                set inf-root-nbors_ids_r4  ( lput id_source-tree inf-root-nbors_ids_r4 )
              ]
            ]
          ]
        ]
        set newly-infected? false                ; remove newly-infected status
      ]
    ]
  ] [ print "ERROR project-infection -> project-root-infection -> ask trees_newly-infected" print error-message set kill-run? true stop ]

  ; Losing-infection remove their who and age from their neighbors' properties
  carefully [
    if is-agentset? trees_losing-infection and any? trees_losing-infection [
      print "entered project-root -> losing infection"
      ask trees_losing-infection [
        let id_source-tree who
        ;      print ( word "cell losing: " id_source-tree )
        ; Root spread distance based on management
        if edge? or mgmt = "short-rotation"      [ set root-contact-distance-range distance-range_int ] ; short-rotation and edges
        if mgmt = "long-rotation" and not edge? [ set root-contact-distance-range distance-range_ext ]

        ; loop through rings
        foreach root-contact-distance-range [ current-ring-distance ->
          ; ask cells in each ring
          ask cells at-points ( item ( current-ring-distance - 1 ) ring-list_root ) [
  ;          print ( word "ring " current-ring-distance )
            if potential-tree? [
              ; remove the age and who of each source tree from the relevant distance ring of all potential tree neighbors
              if current-ring-distance = 1 [
  ;              print ( word "entered ring " current-ring-distance )
                let value-position ( position id_source-tree inf-root-nbors_ids_r1 )                  ; Find the position of the tree losing infection in the root nbor lists
  ;              print ( word "value posit? " value-position )
                if ( value-position != false ) [
  ;                print "entered remove loop"
                  set inf-root-nbors_ages_r1 ( remove-item value-position inf-root-nbors_ages_r1 ) ; remove the age at that position
                  set inf-root-nbors_ids_r1  ( remove id_source-tree inf-root-nbors_ids_r1 )       ; remove the who of that tree losing infection
                ]
              ]
              if current-ring-distance = 2 [
  ;              print ( word "entered ring " current-ring-distance )
                let value-position ( position id_source-tree inf-root-nbors_ids_r2 )                  ; Find the position of the tree losing infection in the root nbor lists
  ;              print ( word "value posit? " value-position )
                if ( value-position != false ) [
  ;                print "entered remove loop"
                  set inf-root-nbors_ages_r2 ( remove-item value-position inf-root-nbors_ages_r2 ) ; remove the age at that position
                  set inf-root-nbors_ids_r2  ( remove id_source-tree inf-root-nbors_ids_r2 )       ; remove the who of that tree losing infection
                ]
              ]
              if current-ring-distance = 3 [
  ;              print ( word "entered ring " current-ring-distance )
                let value-position ( position id_source-tree inf-root-nbors_ids_r3 )                  ; Find the position of the tree losing infection in the root nbor lists
  ;              print ( word "value posit? " value-position )
                if ( value-position != false ) [
  ;                print "entered remove loop"
                  set inf-root-nbors_ages_r3 ( remove-item value-position inf-root-nbors_ages_r3 ) ; remove the age at that position
                  set inf-root-nbors_ids_r3  ( remove id_source-tree inf-root-nbors_ids_r3 )       ; remove the who of that tree losing infection
                ]
              ]
              if current-ring-distance = 4 [
  ;              print ( word "entered ring " current-ring-distance )
                let value-position ( position id_source-tree inf-root-nbors_ids_r4 )                  ; Find the position of the tree losing infection in the root nbor lists
  ;              print ( word "value posit? " value-position )
                if ( value-position != false ) [
  ;                print "entered remove loop"
                  set inf-root-nbors_ages_r4 ( remove-item value-position inf-root-nbors_ages_r4 ) ; remove the age at that position (can't use "remove" function because it would remove all)
                  set inf-root-nbors_ids_r4  ( remove id_source-tree inf-root-nbors_ids_r4 )       ; remove the who of that tree losing infection
                ]
              ]
              if length ( sentence inf-root-nbors_ids_r1 inf-root-nbors_ids_r2 inf-root-nbors_ids_r3 inf-root-nbors_ids_r4 ) = 0 [
                set any-inf-root-nbors? false    ; just acknowledging that there is some root infection probability that needs to be accounted for
              ]
            ]
          ]
        ]
        set losing-infection? false              ; remove losing-infection status
      ]
    ]
  ] [ print "Error in project-infection -> project-root-infection -> ask trees_losing-infection" print error-message set kill-run? true stop ]
end

;_______________________________________________________________________________________________
; BACKGROUND PROBABILITY OF insect DISPERSAL ; Check fix test

to calculate-insect-dispersal_background ; EDIT to remove minimum values
  carefully [
    ; Reset the number of stands projecting background dispersal
    set n-stands-projecting-insect-dispersal_background 0
    set prob-insect-dispersal_background_sum 0
    ; Go through each stand
    foreach stand-list [ this-stand-id ->
      ; count the number of trees in the stand
      let stand-trees ( trees_current with [ stand-id = this-stand-id ] )
      ; calculate the proportion of trees infected in the stand
      let pct-trees-infected ( ( count stand-trees with [ infected? ] ) / ( count stand-trees ) * 100 )
      ; if the % of trees infected meets or exceeds the limit for dispersal spillover
      if pct-trees-infected >= insect-disp_background_stand-pct-inf-spillover-threshold [
        set n-stands-projecting-insect-dispersal_background ( n-stands-projecting-insect-dispersal_background + 1 )
      ]
    ]
    ; Calculate the cumulative prob-insect-dispersal
    ifelse n-stands-projecting-insect-dispersal_background > 0 [
      set prob-insect-dispersal_background_sum ( n-stands-projecting-insect-dispersal_background * prob-insect-disp_background )
    ] [ set prob-insect-dispersal_background_sum prob-insect-disp_background ] ; always a minimum
  ] [ print "ERROR in project-infection -> calculate-insect-dispersal_background" print error-message ]

end

;////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

; TIME-PASSES: Various time-related model processes Trees age, gain or lose infection susceptibility, die, old dead trees and stumps disappear
to time-passes
  ; TO COMPLETE
  ;  - Combine where appropriate, minimizing computing to find the relevant cells to apply the processes to
  ;  - Check the order makes sense
  ;  - Check that inf loss and attr loss are applied for the appropriate amount of steps
  ;      - if infection can be lost immediately, that's not realistic b/c there should be at least 7 months of infectiousness
  ;      - Since the temp. res. of the model is 1 yr, this can't be adequately captured.

  ; Begin procedure to update the status of all trees

  ; To calculate cumulative prob root contact using the lists - use the max distance to account for trees near other management
  let root-contact-distance-range ( range 1 ( 4 + 1 ) )

  ; REMOVE PROB INFECTION CAUSED BY ERRORS WITH FLOATING POINT MATH
  ask trees_potential [
    if prob-insect-disp_sum_SC < 0.000001 [ set prob-insect-disp_sum_SC 0 ]
    if prob-insect-disp_sum_PF < 0.000001 [ set prob-insect-disp_sum_PF 0 ]
    if prob-insect-disp_sum_HN < 0.000001 [ set prob-insect-disp_sum_HN 0 ]
  ]


;  file-open errorFile ;delete
  ; MAIN LOOP --> Update status for all current trees
  ask trees_current [
    ; Setup variables to track whether the cell WAS or BECAME infected (used throughout)
    let previously-infected?    infected?     ; reset
    let became-infected?        false         ; reset
    carefully [
      ; (1) INCREASE COUNTER: MORTALITY
      ;     (this goes before infected tree death so that trees dying at this step will not have their time-since-mort counter increased)
      ifelse not alive? [                                    ; Ask trees that are dead
        set time-since-mort ( time-since-mort + 1 )          ; add 1 to their mortality counter - trees killed by thinning or harvest will have time-since-mort = 1 for the same tick in which
                                                             ; they died because the attraction due to the mgmt action that killed them has already been accounted for in this tick

        ; (2) DEAD CELLS ARE RESET AND DESIGNATED TO LOSE ATTRACTION, REDUCE NBOR ATTR, AND DECREASE NBOR PROB-INF (IF INFECTED)
        ; IF dead at least 2 years, tree can lose attraction and be reset; after 4 years (dead-attr-duration), all are reset
        ifelse ( ( ( time-since-mort >= dead-attr-duration_max ) OR
          ( time-since-mort >= dead-attr-duration_min and random-float 1 < prob-attr-dead-loss ) ) AND
          time-since-mort > dead-host-viability-duration_max ; check NEW
        )
        [
          if infected? or newly-infected? [
            set losing-infection? true
;            file-print ( word ticks "," who "," tree? "," alive? "," infected? "," mort-cause ",false,true,time-passes -> step 2" ) ; delete
          ]                      ; note whether it was infected
          if ( mort-cause = "bsrd" ) [
            set losing-attr-dead? true  ; designate as losing insect attraction but only when not in a managed setting (in thinned or harvested stands, attraction loss and reset cell properties managed by remove-attr-mgmt_by-stand)
          ]
          reset-cell-properties       ; reset cell properties only for dead by bsrd, NOT management (separate proceduree)
        ]
        ; ELSE, if not being reset...
        [
          carefully [
            ; (3) PREVIOUSLY DEAD & INF TREES LOSE INFECTION
            ; Check how long ago the tree died relative to inoculum-viability-loss (the amount of time after tree/host death before (a) the pathogen can't survive
            ; or (b) the pathogen can't colonize the host [7-24 months / 0.5-2 years])
            ; IF time-since-mort is >= inoc-viab-loss, THEN lose infection; otherwise, THERE'S A probabilistic loss of infection
            if infected? [
              ; IF the dead tree cannot maintain viable pathogen inoculum (based on time-since-mortality), it loses infection
              ifelse ( ( time-since-mort >= dead-host-viability-duration_min AND ( random 2 = 0 ) ) OR ( time-since-mort >= dead-host-viability-duration_max ) ) [
                ;OR (  AND ( random 2 = 0 ) ) ; fix later random ( max - min + 1 ), more flex ; what?
                ;( time-since-mort > dead-host-viability-duration_min ) AND ( random-float 2 = 0 )
                ;                file-print ( word ticks "," who "," tree? "," alive? "," infected? "," mort-cause ",false,true,time-passes -> step 3" ) ; delete
                set losing-infection? true
                set infected? false
                set inf-root? false
                set inf-SC? false
                set inf-PF? false
                set inf-HN? false
                set inf-initial? false
              ]
              ; ELSE as long as not losing infection...
              ; (4) INCREASE COUNTER: INFECTION
              ;     (this goes BEFORE infection so that trees infected during this step will not have their time-since-inf counter increased)
              ;     (this goes AFTER INCREASE COUNTER: MORTALITY, DEAD CELL RESET, and INFECTION LOSS to identify trees eligible to continue being infected)
              ; (4a) for dead trees
              [ set time-since-inf  ( time-since-inf + 1 ) ]      ;   INCREASE INFECTION COUNTER
            ]
          ] [  print "ERROR: time-passes -> Steps 3-4" print error-message set kill-run? true stop ]
        ]
      ]
      [ ; (4b) for live trees
        set time-since-inf  ( time-since-inf  + 1 )                                                                ; INCREASE INFECTION COUNTER
        set age-bsrd-mort-suscept ( one-of ( item ( age - 1 ) age-bsrd-mort-suscept_years-after-infection_list ) ) ; UPDATE PROBABILITY OF MORTALITY (IF time-since-inf = 1, otherwise, already updated)
      ]
    ] [ print "ERROR: time-passes -> Steps 1-4" print error-message set kill-run? true stop ]


    if tree? [ ; IF STILL A TREE (NOT RESET EARLIER IN PROCEDURE) - Check if tree? before continuing, since some trees have now been reset to non-trees
      ; (5) PROBABILITY OF INFECTION CALCULATION - Clear prob dispersal that are too low
      set prob-root-contact_union 0          ; clear prob root contact from last tick (always changes because of aging and tree growth)
      set prob-inf-root 0


      ; (6) RECALCULATE PROBABILITIES OF INFECTION, INFECT TREES
      ; Clear prob-root-contact and prob-inf-root for all trees (needs to be recalculated to account for changing age) ; check
      ; (A) RECALCULATE PROBABILITY OF INFECTION
      ; (i) PROB ROOT INFECTION
      carefully [
        if max-root-distance_m > 0 [                                                                                ; if root transmission is on (non-zero distance for transmission)
          if any-inf-root-nbors? [                                                                                  ; if the tree has infected root neighbors. CHECK Redundant?
            let inf-root-nbor-age-lists ( list                                                                      ; create a nested list with all root-neighbor age lists (to iterate through)
              inf-root-nbors_ages_r1 inf-root-nbors_ages_r2 inf-root-nbors_ages_r3 inf-root-nbors_ages_r4
            )
            ;print ( word "age: " age )
            ;type "nested age list " print inf-root-nbor-age-lists
            ifelse length ( sentence inf-root-nbor-age-lists ) != 0 [                                               ; if there are any values in the root neighbors age lists, proceed with calculation, else
                                                                                                                    ; no values in the root neighbors age lists, set prc and prob root infection to 0
              foreach root-contact-distance-range [ this-distance ->                                                ; for each distance in the root radius
                let current-age-list ( item ( this-distance - 1 ) inf-root-nbor-age-lists )                             ; get the infected nbor age list for that distance

                ; Calcuate the prob root contact with the inf nbor trees in that age list by mapping the draw prc function to the neighbor age list and summing,
                let prc_current-list []            ; store in temp variable prc_current-list
                carefully [  set prc_current-list ( calculate-prob-root-contact age current-age-list this-distance ) ]
                [ print "Error during calc-prob-root-contact/draw-prob-root-contact" print error-message set kill-run? true stop ]
                set prob-root-contact_union ( prob-root-contact_union + prc_current-list ) ; add to the prc for this tree
              ]
              ; Calculate prob infection via roots
              set prob-inf-root ( age-bsrd-inf-suscept * prob-root-contact_union * prob-root-transm )
            ]
            [ ; ELSE there are no infected root nbors - set any infected root neighbors to false to indicate lack of inf root nbors and prevent calculation
              set any-inf-root-nbors? false
            ]
          ]
        ]
      ] [ print "ERROR: time-passes -> Step 7.A.i" print error-message set kill-run? true stop ]

      carefully [
        if ( max-insect-distance_SC_m + max-insect-distance_PF_m + max-insect-distance_HN_m ) > 0 OR background-prob-insect-dispersal? [  ; if insect dispersal is active, recalculate prob insect infection
          ; PROB. insect INFECTION

          ; Attraction is applied to prob insect wounding
          ; (i) initialize temp variables
          let prob-insect-wound_calc 0 ; Initialize a temporary variable to modify for probability of insect wounding
          let attr-mgmt_to-use 0
          let attr-dead_to-use 0
          let attr-inf_to-use  0

          ; (ii) Set management, inf, and dead attraction based on version used
          ; For the "median" version, the median is used for all three
          if attr-inf-dead_version = "median" OR attr-inf-dead_version = "proportional" [
            set attr-mgmt_to-use ( median attr-mgmt )
            set attr-dead_to-use ( median attr-dead )
            ; For the "proportional" version, the median is used for attr-mgmt and attr-dead
            if attr-inf-dead_version = "median" [
              set attr-inf_to-use  ( median attr-inf  )
            ]
          ]
          if attr-inf-dead_version = "max" [
            set attr-mgmt_to-use ( max  attr-mgmt )
            set attr-dead_to-use ( max  attr-dead )
            set attr-inf_to-use  ( max  attr-inf  )
          ]

          ; CONSERVATIVE/PROPORTIONAL: weight the infected neighbor trees based on distance and set an attr-inf proportional to the amount of inf trees
          ; in the attraction radius versus the max possible, with a threshold for what is considered an infection center (with the full range
          ; of attr-inf values)
          if attr-inf-dead_version = "proportional" [
            ; draw an attraction factor from the appropriate list
            ifelse alive? [ set attr-inf_to-use ( one-of attr-effect_inf-vs-live_list ) ]
            [ set attr-inf_to-use ( one-of attr-effect_inf-vs-dead_list ) ]
            ; scale based on number and distance of inf trees in surroundings
            let proportion-of-inf-nbors-relative-to-threshold ( min ( list ( attr-inf_source / inf-center-attr-threshold_n-source-trees ) 1 ) ) ; prop relative to threshold - max = 1
            set attr-inf_to-use ( max ( list ( attr-inf_to-use * proportion-of-inf-nbors-relative-to-threshold ) 1 ) )                          ; attr-inf effect min = 1
            set attr-inf attr-inf_to-use
          ]

          ; Combine attraction effects and apply to prob wound
          ifelse ( stand-w-thin-disturbance? or stand-w-harv-disturbance? ) [ ; If the cell is in a thinned or harvested stand, don't apply dead attraction, which would be redundant
            set prob-insect-wound_calc ( prob-insect-wound * attr-road * attr-mgmt_to-use * attr-inf_to-use )
          ]
          [ set prob-insect-wound_calc ( prob-insect-wound * attr-road * attr-mgmt_to-use * attr-dead_to-use * attr-inf_to-use ) ]
          ; Prob insect wounding should not exceed 1
          set prob-insect-wound_calc ( min ( list prob-insect-wound_calc 1 ) )

          ; Set probability of insect infection by each insect, should not exceed 1
          set prob-inf-SC ( min ( list ( age-bsrd-inf-suscept * prob-inf_SC_base * prob-insect-disp_sum_SC * prob-insect-wound_calc ) 1 ) )                                              ; limit to 1
          set prob-inf-PF ( min ( list ( age-bsrd-inf-suscept * prob-inf_PF_base * ( prob-insect-disp_sum_PF + prob-insect-dispersal_background_sum ) * prob-insect-wound_calc ) 1 ) ) ; limit to 1
          set prob-inf-HN ( min ( list ( age-bsrd-inf-suscept * prob-inf_HN_base * ( prob-insect-disp_sum_HN + prob-insect-dispersal_background_sum ) * prob-insect-wound_calc ) 1 ) ) ; limit to 1
        ]
      ][ print "Error in time-passes -> Step 7.A.ii: Prob insect inf calc" print error-message set kill-run? true stop ]

      ; (B) DETERMINE WHETHER INFECTION OCCURS & INFECT TREES
      ;    (i) IDENTIFY TREES THAT CAN BECOME INFECTED (by setting "at-risk?" to TRUE)
      ; If the tree has been infected by all four mechanisms OR has not been dead for too long to be infected
      ifelse ( ( infected? AND inf-root? AND inf-SC? AND inf-PF? AND inf-HN? ) OR losing-infection? OR losing-attr-dead? ) [
        set at-risk? false
      ] ; ELSE at risk to become infected
      [ set at-risk? true ]
      ;    (ii) Infect trees
      carefully [
        if at-risk? [          ; Only check whether infection occurs if at-risk
          if infected? [ set previously-infected? true set newly-infected? false ]   ; If the cell was infected before, make sure that is tracked and that the tree is not labeled as newly-inf

          ; (1) ROOT INFECTION
          if not inf-root? [                                                         ; IF not already infected by roots
            if ( random-float 1 < prob-inf-root )  [                                 ; infection occurs probabilistically based on prob-inf-root
              set inf-root? true set became-infected? true                           ; set infected true for that mechanism
            ]
          ]
          ; (2) insect INFECTION
          ;  (a) Stremnius carinatus (SC)
          if not inf-SC? [                                                           ; IF not already infected by SC, check whether infection happened
            if ( random-float 1 < prob-inf-SC )  [                                   ; infection occurs probabilistically based on prob-inf-SC
              set inf-SC? true set became-infected? true                             ; set infected true for that mechanism
            ]
          ]
          ;   (b) Pissodes fasciatus (PF)
          if not inf-PF? [                                                           ; IF not already infected by PF, check whether infection happened
            if ( random-float 1 < prob-inf-PF )  [                                   ; infection occurs probabilistically based on prob-inf-PF
              set inf-PF? true set became-infected? true                             ; set infected true for that mechanism
            ]
          ]
          ;   (c) Hylastes nigrinus (HN)
          if not inf-HN? [                                                           ; IF not already infected by HN, check whether infection happened
            if ( random-float 1 < prob-inf-HN   )  [                                 ; infection occurs probabilistically based on prob-inf-HN
              set inf-HN? true set became-infected? true                             ; set infected true for that mechanism
            ]
          ]
          ; IF the tree was NOT INFECTED and BECAME INFECTED,
          if ( not previously-infected? ) and became-infected? [
            set infected? true                                                                   ; set as infected
            set newly-infected? true                                                             ; designate as newly infected                                                                    ;
            ifelse alive? [ set prob-insect-wound ( one-of prob-insect-wound_live-inf_list ) ]   ; update prob-insect-wound to refect current cell state if alive...
            [ set prob-insect-wound ( one-of prob-insect-wound_dead-inf_list ) ]                 ; or dead
            set bsrd-inf_cumul ( bsrd-inf_cumul + 1 )                                            ; add to total infection count
            if mgmt = "short-rotation" [ set bsrd-inf_short-rotation_cumul ( bsrd-inf_short-rotation_cumul + 1 ) ]          ; add to total infection count for short-rotation
            if mgmt = "long-rotation" [ set bsrd-inf_long-rotation_cumul ( bsrd-inf_long-rotation_cumul + 1 ) ]          ; and long-rotation stands
            set age-bsrd-mort-suscept  ( one-of ( item ( age - 1 ) age-bsrd-mort-suscept_year-of-infection_list ) ) ; change bsrd mortality susceptibility to reflect newly infected status
          ]
        ]
      ] [ print "Error in time-passes -> Step 7.B: Infect trees" print error-message set kill-run? true stop ]


      if alive? [
        ; (7) LIVE, INFECTED TREES DIE (including NEWLY-INFECTED, which can occasionally die in the year of infection, especially when seedlings)
        if infected? [
          carefully [
            if ( ( random-float 1 ) < age-bsrd-mort-suscept ) [
              set alive? false set mort-cause "bsrd"                                         ; Kill the tree if YES and set BSRD as cause of death
              set new-bsrd-mortality? true
              set prob-insect-wound ( one-of prob-insect-wound_dead-inf_list )               ; Change prob wound to that of a dead tree
              set bsrd-mort_cumul ( bsrd-mort_cumul + 1 )                                    ; Increase the cumulative count of mortalities caused by bsrd overall...
              if mgmt = "short-rotation" [ set bsrd-mort_short-rotation_cumul ( bsrd-mort_short-rotation_cumul + 1 ) ]  ; and for each management practice
              if mgmt = "long-rotation" [ set bsrd-mort_long-rotation_cumul ( bsrd-mort_long-rotation_cumul + 1 ) ]
            ]
          ][ print "ERROR time-passes -> Step 8: Live, infected trees die" print error-message set kill-run? true stop ]

          ; (8) NEWLY-INFECTED DEAD TREES LOSE INFECTION (PREVIOUSLY INFECTED TREES ALREADY ACCOUNTED FOR IN STEPS 3 & 4)
          carefully [
            if newly-infected? [
              ; IF the dead tree cannot maintain viable pathogen inoculum (based on time-since-mortality), it loses infection
              if ( ( time-since-mort >= dead-host-viability-duration_min AND ( random 2 = 0 ) ) OR ( time-since-mort >= dead-host-viability-duration_max ) ) [
                set newly-infected? false
                set infected?       false
                set inf-initial?    false
                set inf-root?       false
                set inf-SC?         false
                set inf-PF?         false
                set inf-HN?         false
                ;set losing-infection? false ; no longer newly infected to avoid spreading then removing prob inf
              ]
            ]
          ] [ print "ERROR time-passes -> Step 9: Newly-infected dead trees lose infection" print error-message ]
        ]
        ; (9) AGING: LIVE TREES get older and have their bsrd susceptibilities change based on their age
        carefully [
          set age ( age + 1 )                                               ; Live trees set age to age + 1
          if age-inf-suscept != "age-indep" [
            set age-bsrd-inf-suscept ( one-of ( item ( age - 1 ) age-bsrd-inf-suscept_list ) )     ; Draw a new age/inf-suscept from the list
          ]
          if infected? [                               ; If the tree is infected, set a susceptibility to mortality value based on age and time-since-infection
            ifelse time-since-inf = 0 [ set age-bsrd-mort-suscept ( one-of ( item ( age - 1 ) age-bsrd-mort-suscept_year-of-infection_list ) ) ]
            [ set age-bsrd-mort-suscept ( one-of ( item ( age - 1 ) age-bsrd-mort-suscept_years-after-infection_list ) ) ]
          ]
        ] [ print "ERROR time-passes -> Step 10: aging" print error-message ]
      ]
      ; (10) Root neighbor ages increase to account for ages of infected neighbors
      carefully [
        if any-inf-root-nbors? [    ; If the tree has infected root neighbors, update ages (if any) in each list for the next tick
          if ( length inf-root-nbors_ages_r1 ) > 0 [ set inf-root-nbors_ages_r1 ( map [ i -> i + 1 ] inf-root-nbors_ages_r1 ) ] ; increment all ages by 1
          if ( length inf-root-nbors_ages_r2 ) > 0 [ set inf-root-nbors_ages_r2 ( map [ i -> i + 1 ] inf-root-nbors_ages_r2 ) ]
          if ( length inf-root-nbors_ages_r3 ) > 0 [ set inf-root-nbors_ages_r3 ( map [ i -> i + 1 ] inf-root-nbors_ages_r3 ) ]
          if ( length inf-root-nbors_ages_r4 ) > 0 [ set inf-root-nbors_ages_r4 ( map [ i -> i + 1 ] inf-root-nbors_ages_r4 ) ]
        ]
      ] [ print "ERROR time-passes -> Step 11: Root neighbor ages increase" print error-message ]
    ]
  ]


  ; (11) MANAGED LOSE ATTRACTION
  ; REMOVE ATTR-MGMT for cells after appropriate amount of time (2-4 years influencing infection)
  carefully [
    remove-attr-mgmt_by-stand
  ] [ print "ERROR: time-passes -> Step 1: remove-attr-mgmt_by-stand" print error-message set kill-run? true stop ]
;  file-close ;delete
end

;---------------------------------------------------------------- NOT USED BUT RELEVANT FOR DESIGN ----------------------------------------------------------------

;to recalculate-prob-inf
;  ; CHECK: multiplying CUMULATIVE dispersal probability by ATTR
;  ; CHECK: add together attr before multiplying?)
;  set prob-inf-root ( age-bsrd-inf-suscept * prob-root-contact_union )
;  set prob-inf-SC   ( age-bsrd-inf-suscept * prob-inf_SC_base * prob-insect-wound * prob-insect-disp_sum_SC * ( attr-road ) * ( max attr-mgmt ) * attr-inf * attr-dead )
;  set prob-inf-PF   ( age-bsrd-inf-suscept * prob-inf_PF_base * prob-insect-wound * prob-insect-disp_sum_PF * ( attr-road ) * ( max attr-mgmt ) * attr-inf * attr-dead )
;  set prob-inf-HN   ( age-bsrd-inf-suscept * prob-inf_HN_base * prob-insect-wound * prob-insect-disp_sum_HN * ( attr-road ) * ( max attr-mgmt ) * attr-inf * attr-dead )
;end


;_____________________________________________________________________________________________________________
;_____________________________________________________________________________________________________________
;_____________________________________________________________________________________________________________
; MANAGEMENT
;_____________________________________________________________________________________________________________

to manage-stands  ; Perform the three management activities (thinning, harvest, and regeneration) on appropriate stands (based on their age and management history)

  ; Reset trackers for harvest and thinning
  set thinning-occurred?        false       ; bool to track whether a thinning occurred
  set harvest-occurred?         false       ; bool to track whether a harvest occurred

  ; THIN STANDS
  if ( length long-rotation-stand-list > 0 ) and thinning? [                                                                               ; If any stands with long-rotation management
    set trees-to-thin ( trees_current with [ mgmt = "long-rotation" AND ( age = 15 or age = 35 or age = 55 ) AND alive? ] )  ; Find the trees to thin
    if any? trees-to-thin [ thin-stands ]                                                                                ; Thin them
  ]
  ; HARVEST STANDS
  if harvest? [
    set trees-to-harvest ( trees_current with [ ( ( mgmt = "short-rotation" AND age = rotation:short-rotation ) OR ( mgmt = "long-rotation" AND age = rotation:long-rotation ) ) AND alive? ] ) ; If any trees to harvest
    if any? trees-to-harvest [ harvest-stands ] ; Harvest them
  ]

  ; REGENERATE STANDS
  if ( length regen-list-short-rotation > 0 ) or ( length regen-list-long-rotation > 0 ) [  ; IF any stands that have been harvested
    regenerate-stands                                                              ; replant them
  ]

end


; THIN long-rotationLY MANAGED STANDS (with % removed based on age)
to thin-stands  ; TEST - CHECK - IN PROGRESS
  ; For increased flexibility in thinning:
  ; IMPLEMENT: globals for age-pct = 15
  ; IMPLEMENT: globals for age-ct1 = 35
  ; IMPLEMENT: globals for age-ct2 = 55
  ; IMPLEMENT: globals for prop-rm-pct = 0.311294766
  ; IMPLEMENT: globals for prop-rm-ct1 = 0.466666667
  ; IMPLEMENT: globals for prop-rm-ct2 = 0.375000000

  ; Create temp lists for stands
  let thinned-stands_pct_list []
  let thinned-stands_ct1_list []
  let thinned-stands_ct2_list []

  ; If there are long-rotationly managed stands
  ; Perform PCT thinning on 15-year-old trees in long-rotation mgmt stands
  let trees-to-thin_pct ( trees-to-thin with [ age = 15 ] )                                    ; subset trees to thin
  if any? trees-to-thin_pct [                                                                  ; if there are any of pct age
    ask trees-to-thin_pct [
      if ( random-float 1 < 0.311294766 ) [
        set alive? false set mort-cause "thin"                                                 ; have a proportion of them cut
        ifelse infected? [ set prob-insect-wound ( one-of prob-insect-wound_dead-inf_list ) ]       ; Set prob-wound
        [ set prob-insect-wound ( one-of prob-insect-wound_dead-nonInf_list ) ]
      ]
    ]
    set thinned-stands_pct_list ( remove-duplicates ( [ stand-id ] of trees-to-thin_pct ) )    ; make a list of stands thinned w/pct
    set thinning-occurred? true  ; mark that a thinning occurred this year
  ]

  ; Perform CT1 thinning on 35-year-old trees in long-rotation mgmt stands
  let trees-to-thin_ct1 ( trees-to-thin with [ age = 35 ] )                                    ; subset trees to thin
  if any? trees-to-thin_ct1 [                                                                  ; if there are any of ct1 age
    ask trees-to-thin_ct1 [
      if ( random-float 1 < 0.466666667 ) [
        set alive? false set mort-cause "thin"                                                 ; have a proportion of them cut
        ifelse infected? [ set prob-insect-wound ( one-of prob-insect-wound_dead-inf_list ) ]       ; Set prob-wound
        [ set prob-insect-wound ( one-of prob-insect-wound_dead-nonInf_list ) ]
      ]
    ]
    set thinned-stands_ct1_list ( remove-duplicates ( [ stand-id ] of trees-to-thin_ct1 ) )    ; make a list of stands thinned w/ct1
    set thinning-occurred? true  ; mark that a thinning occurred this year
  ]

  ; Perform CT2 thinning on 55-year-old trees in long-rotation mgmt stands
  let trees-to-thin_ct2 ( trees-to-thin with [ age = 55 ] )                                    ; subset trees to thin
  if any? trees-to-thin_ct2 [                                                                  ; if there are any of ct2 age
    ask trees-to-thin_ct2 [
      if ( random-float 1 < 0.375 ) [
        set alive? false set mort-cause "thin"                                                 ; have a proportion of them cut
        ifelse infected? [ set prob-insect-wound ( one-of prob-insect-wound_dead-inf_list ) ]       ; Set prob-wound
        [ set prob-insect-wound ( one-of prob-insect-wound_dead-nonInf_list ) ]
      ]
    ]
    set thinned-stands_ct2_list ( remove-duplicates ( [ stand-id ] of trees-to-thin_ct2 ) )    ; make a list of stands thinned w/ct2
    set thinning-occurred? true  ; mark that a thinning occurred this year
  ]

  ; Set list of stands thinned 1 yr ago to the stands thinned in the current year
  set thin-0yr-ago remove-duplicates ( sentence thinned-stands_pct_list thinned-stands_ct1_list thinned-stands_ct2_list )

  ; Create attraction if spread-infection is active
  if spread-setup? and max-attraction-distance_m > 0 [ add-attr-mgmt_thin thinned-stands_pct_list thinned-stands_ct1_list thinned-stands_ct2_list ]

  ; Clear agentset of trees to be thinned
  set trees-to-thin empty-agentset
  print ( word "stands thinned: " thin-0yr-ago )
end

to harvest-stands ; TEST - CHECK - IN PROGRESS - NOTHING CHANGED SO FAR
  ; Harvest the trees at rotation age
  ask trees-to-harvest [
    set alive? false set mort-cause "harvest"                                                                  ; Kill trees, set mort cause
    set time-since-mort 0                                                                                      ; Reset time-since-mort (just in case)
    ifelse infected? [ set prob-insect-wound ( one-of prob-insect-wound_dead-inf_list ) ]                      ; Set prob-wound
    [ set prob-insect-wound ( one-of prob-insect-wound_dead-nonInf_list ) ]
  ]
  ; create a list from that agentset of stands to regenerate for each management class
  let regen-list-short-rotation_to-add ( remove-duplicates ( [ stand-id ] of trees-to-harvest with [ mgmt = "short-rotation" ] ) ) ; Get a list of stand-ids of trees to harvest
  set regen-list-short-rotation ( sentence regen-list-short-rotation regen-list-short-rotation_to-add )                                 ; Add to regen list so these stands can be regenerated
  let regen-list-long-rotation_to-add ( remove-duplicates ( [ stand-id ] of trees-to-harvest with [ mgmt = "long-rotation" ] ) ) ; Get a list of stand-ids of trees to harvest
  set regen-list-long-rotation ( sentence regen-list-long-rotation regen-list-long-rotation_to-add )                                 ; Add to regen list so these stands can be regenerated
  set harv-0yr-ago ( sentence regen-list-short-rotation_to-add regen-list-long-rotation_to-add ) ; Recently harvested added to the list
  print ( word "stands harvested: " harv-0yr-ago )
  ; ADD ATTRACTION TO HARVESTED STANDS
  if ( spread-setup? and max-attraction-distance_m > 0 ) [ add-attr-mgmt_harvest regen-list-short-rotation_to-add regen-list-long-rotation_to-add ]

  set harvest-occurred? true  ; mark that a harvest occurred this year
  ; clear agentset of trees to harvest
  set trees-to-harvest empty-agentset
end


to regenerate-stands ; REGENERATE STANDS AFTER HARVEST
  ; FIX: For replant-delays >= 4 years or even for very short dead-attr-duration, this won't work because trees will be set to non-trees and time-since-mort will stop increasing. Not a problem now but maybe later.
  ; Check the conditions used to regen stands and whether or not
  foreach regen-list-short-rotation [ this-stand ->                                                                              ; FOR EACH stand on the regen list
    if ( any? ( ( stand this-stand ) with [ tree? and ( mort-cause = "harvest" AND time-since-mort >= replant-delay ) ] ) ) [   ; IF there are any ready to replant (based on cause of morality and replant delay)

      ifelse not any? ( ( stand this-stand ) with [ tree? AND pxcor mod 3 = 0 AND pycor mod 3 = 0 ] ) [                     ; THEN, IF there are no trees on patches with coord divisible by 3
        ask ( stand-patches this-stand ) with [ ( pxcor ) mod 3 = 0 AND ( pycor ) mod 3 = 0 ]         [ create-tree 1 ]     ; regenerate trees on those patches
      ]                                                                                                                     ; ELSE IF those patches have (dead) trees
      [ ask ( stand-patches this-stand ) with [ ( pxcor + 2 ) mod 3 = 0 AND ( pycor + 2 ) mod 3 = 0 ] [ create-tree 1 ] ]   ; offset the regeneration
      set regen-list-short-rotation ( remove this-stand regen-list-short-rotation )                                                   ; Remove the stand from the regen list
      ask ( stand this-stand ) [ set time-since-mort 0 ]                                                                    ; Reset time-since-mort
      print ( word "stand " this-stand " regenerated." )
    ]
  ]
  foreach regen-list-long-rotation [ this-stand ->                                                                              ; FOR EACH stand on the regen list
    if ( any? ( ( stand this-stand) with [ tree? and mort-cause = "harvest" AND time-since-mort >= replant-delay ] ) ) [    ; IF there are any ready to replant (based on cause of morality and replant delay)
      ifelse not any? ( ( stand this-stand ) with [ tree? AND pxcor mod 2 = 0 AND pycor mod 2 = 0 ] ) [                     ; THEN, IF there are no trees on patches with coord divisible by 2
        ask ( stand-patches this-stand ) with [ ( pxcor ) mod 2 = 0 AND ( pycor ) mod 2 = 0 ] [ create-tree 1 ]             ; regenerate trees on those patches
      ]                                                                                                                     ; ELSE IF those patches have (dead) trees
      [ ask ( stand-patches this-stand ) with [ ( pxcor + 1 ) mod 2 = 0 AND ( pycor + 1 ) mod 2 = 0 ] [ create-tree 1 ] ]           ; offset the regeneration
      set regen-list-long-rotation ( remove this-stand regen-list-long-rotation )                                                   ; Remove the stand from the regen list
      ask ( stand this-stand ) [ set time-since-mort 0 ]                                                                    ; Reset time-since-mort
      print ( word "stand " this-stand " regenerated." )
    ]
  ]
  ; GLOBAL AGENTSETS
  carefully [
    set trees_current ( trees )   ; set an agentset for all current trees
  ] [ print "ERROR regenerate-stands -> setting trees_current" print error-message ]
end

;_____________________________________________________________________________________________________________
; COLOR PROCEDURES - used to color landscape in a way that visually represents stand and cell characteristics
;_____________________________________________________________________________________________________________

; Define colors for easy visual understanding of the model
to setup-colors
  set  color_patches                9.2     ; background patch color
  set  color_empty                  9.7     ; empty (non-tree) cells
  set  color_tree-live-notinf      66.0     ; non-infected trees
  set  color_tree-live-inf         15.5     ; infected trees
  set  color_tree-dead-notinf      36.0     ; dead, non-infected trees
  set  color_tree-dead-inf         18.5     ; dead, infected trees
  set  color_road                   6.0     ; cells that are roads
end

; Colors all cells based on their states (e.g., tree vs. empty, alive vs. dead,
; infected vs. not-infected)
to color-world
  ask patches [ set pcolor color_patches ]
  ask cells [
    if ( not tree? )                               [ set color color_empty ]
    if ( alive? AND tree? AND not infected? )      [ set color color_tree-live-notinf ]
    if ( alive? AND infected? AND tree? )          [ set color color_tree-live-inf ]
    if ( tree? AND not alive? AND not infected? )  [ set color color_tree-dead-notinf ]
    if ( infected? AND tree? AND not alive? )      [ set color color_tree-dead-inf ]
    if ( road? )                                   [ set color color_road ]
  ]
end

to color-cells ; Colors cells to "color-world" colors without changing patch colors
  ask cells [
    if ( not tree? )                               [ set color color_empty ]
    if ( alive? AND tree? AND not infected? )      [ set color color_tree-live-notinf ]
    if ( alive? AND infected? AND tree? )          [ set color color_tree-live-inf ]
    if ( tree? AND not alive? AND not infected? )  [ set color color_tree-dead-notinf ]
    if ( infected? AND tree? AND not alive? )      [ set color color_tree-dead-inf ]
    if ( road? )                                   [ set color color_road ]
  ]
end

; COLOR STANDS BY AGE
to color-by-age ; The darker the shade of red, the older the stand
  color-world
  let min-age min [ age ] of cells
  let max-age max [ age ] of cells
  ifelse max-age != min-age [ ask cells with [ tree? and alive? ] [ set color scale-color red age ( max-age * 1.1 ) ( min-age * 0.9 ) ] ] [
    ask cells with [ tree? and alive? ] [ set color red + 2 ] ]
  ; if any? cells with [ edge? ] [ ask edges [ set color [ color ] of self + 1.25 ] ] ; REMOVE?
end

; COLOR STANDS BY STAND ID #
; Stand colors are set sequentially based on stand ID #, from darkest to lightest shade of blue
to color-by-stand-id ; UPDATED - incorporate
  ask patches [
    if config = "random" OR config = "clustered" [ set pcolor scale-color blue p_stand-id ( ( n-stands ) + 1 ) 1 ]
    if config = "blocks"                                [ set pcolor scale-color blue p_stand-id ( ( n-blocks ) + 1 ) 1 ]
  ]
  ask cells [ set color ( ( [ pcolor ] of patch-here ) + 1 ) ]
  ;if any? cells with [ edge? ] [ ask edge-cells [ set color [ color ] of self + 1.25 ] ] ; remove?
end

; Set colors to unique management classes
to color-by-mgmt
  ask patches [
    if p_mgmt = 1  [ set pcolor green ]
    if p_mgmt = 2  [ set pcolor blue ]
    if p_mgmt = 3  [ set pcolor brown ]
  ]
  ask cells [ set color ( ( [ pcolor ] of patch-here ) + 1 ) ]
end

; Set colors to unique combinations of management and stand-id
to color-by-mgmt-&-stand-id
  ; short-rotation
  ask patches with [ p_mgmt = 1 ] [
    if config = "random" OR config = "clustered" [
      set pcolor scale-color green p_stand-id ( max short-rotation-stand-list + 50 ) ( min short-rotation-stand-list - 50)
    ]
    if config = "blocks" [ set pcolor scale-color green p_stand-id ( max short-rotation-stand-list + 50 ) ( min short-rotation-stand-list - 50 )
    ]
  ]
  ; long-rotation
  ask patches with [ p_mgmt = 2 ] [
    if config = "random" OR config = "clustered" [
      set pcolor scale-color blue p_stand-id ( max long-rotation-stand-list + 50 ) ( min long-rotation-stand-list - 50)
    ]
    if config = "blocks" [ set pcolor scale-color blue p_stand-id ( max long-rotation-stand-list + 50 ) ( min long-rotation-stand-list - 50 )
    ]
  ]
  ; SET-ASIDE
  ask cells with [ p_mgmt = 3 ] [
    if config = "random" OR config = "clustered" [
      set pcolor scale-color brown p_stand-id ( max set-aside-stand-list + 50 ) ( min set-aside-stand-list - 50 )
    ]
    if config = "blocks" [
      set pcolor scale-color brown p_stand-id ( max set-aside-stand-list + 50 ) ( min set-aside-stand-list - 50 )
    ]
  ]
  ask cells [ set color ( ( [ pcolor ] of patch-here ) + 1 ) ]
end



; Color based on inf suscept
to color-by-inf-susc
  ask cells with [ not tree? ] [ set color color_empty ]
  if any? trees_current [
    ask trees_current [
      set color scale-color red age-bsrd-inf-suscept 1 0
    ]
  ]
end

to color-by-root
  if any? trees_current [
    let max-value ( max [ prob-inf-root ] of trees_current )
    ask trees_current [
      ;set color scale-color red prob-inf-root 1 0
      ifelse prob-inf-root = 0 [ set color color_patches ] [
        set color scale-color red prob-inf-root max-value 0.0000000000001
      ]
    ]
  ]
end

to color-by-SC
  if any? trees_current [
    let max-value ( max [ prob-inf-SC ] of trees_current )
    ask trees_current [
      ;set color scale-color red prob-inf-SC 1 0
      set color scale-color red prob-inf-SC max-value -0.001
    ]
  ]
end

to color-by-PF
  if any? trees_current [
    let max-value ( max [ prob-inf-PF ] of trees_current )
    ask trees_current [
      ;set color scale-color red prob-inf-PF 1 0
      set color scale-color red prob-inf-PF max-value -0.001
    ]
  ]
end

to color-by-HN
  if any? trees_current [
    let max-value ( max [ prob-inf-HN ] of trees_current )
    ask trees_current [
      ;set color scale-color red prob-inf-HN 1 0
      set color scale-color red prob-inf-HN max-value -0.001
    ]
  ]
end
@#$#@#$#@
GRAPHICS-WINDOW
346
48
354
57
-1
-1
1.0E-10
1
10
1
1
1
0
1
1
1
0
25
0
25
0
0
1
ticks
30.0

CHOOSER
6
12
144
57
comp
comp
"asus-big" "asus-small" "law" "remote" "cgrb"
0

CHOOSER
7
65
145
110
set-world-size?
set-world-size?
"manual" "pre-set" "automatic"
1

CHOOSER
155
62
253
107
world-size
world-size
"10 x 10" "26 x 26" "50 x 50" "100 x 100" "300 x 300" "500 x 500" "1000 x 1000" "2400 x 2400" "2952 x 2952"
1

SLIDER
4
116
158
149
stand-size_mean_acres
stand-size_mean_acres
0
100
8.0
1
1
NIL
HORIZONTAL

BUTTON
163
117
258
150
check-world-size
check-world-size
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

CHOOSER
6
160
100
205
base
base
"patch"
0

CHOOSER
112
165
250
210
config
config
"random" "blocks" "clustered"
0

CHOOSER
8
210
106
255
n-blocks
n-blocks
1 2 4 9 16 25 36 49 64 81 100 121 144 169 196 225 256 289 324 361 400
0

SLIDER
113
221
246
254
n-stands
n-stands
0
300
1.0
1
1
NIL
HORIZONTAL

SWITCH
9
264
120
297
limit-ticks?
limit-ticks?
0
1
-1000

SLIDER
126
263
251
296
max-ticks
max-ticks
0
1000
160.0
1
1
yr
HORIZONTAL

SWITCH
10
302
118
335
visualize?
visualize?
1
1
-1000

SWITCH
124
303
236
336
track-run?
track-run?
1
1
-1000

BUTTON
281
10
381
43
NIL
setup-colors
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
388
10
458
43
NIL
display
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
463
10
543
43
color-world
no-display color-world display
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
547
10
619
43
NIL
color-cells
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
632
9
743
42
NIL
color-by-mgmt
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
748
10
874
43
NIL
color-by-stand-id
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
880
10
1021
43
NIL
color-by-mgmt-&-stand-id
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
1025
10
1132
43
clearCellColor
ask cells [ set color color_empty ]
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
1166
76
1230
109
resize
ask cells [ set size 1 ]
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
1233
76
1297
109
HL-inf
ask infected-trees [ set color red set size 2 ]\n
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
1302
76
1394
109
HL-new-inf
ask cells with [ newly-infected? ] [\nset color magenta set size 3]
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
1396
76
1496
109
HL-losing-inf
ask cells with [ losing-infection? ] [\nset color violet set size 3 ]
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

CHOOSER
1266
115
1458
160
sensitivity-analysis-param
sensitivity-analysis-param
"none" "age-infection-susceptibility" "age-bsrd-mortality-susceptibility" "prob-root-contact" "prob-root-transmission" "prob-vector-dispersal" "prob-vector-dispersal_SC" "prob-vector-dispersal_PF" "prob-vector-dispersal_HN" "prob-vector-transmission" "prob-vector-transmission_SC" "prob-vector-transmission_PF" "prob-vector-transmission_HN" "prob-vector-infested" "prob-vector-infested_SC" "prob-vector-infested_PF" "prob-vector-infested_HN" "prob-vector-wound" "attr-effect_mgmt" "attr-effect_road" "attr-effect_inf" "attr-effect_dead"
0

SLIDER
1266
163
1457
196
sensitivity-analysis-multiplier
sensitivity-analysis-multiplier
0
2
1.0
0.025
1
NIL
HORIZONTAL

SLIDER
1262
200
1466
233
max-attraction-distance_m
max-attraction-distance_m
0
30
15.0
1
1
m
HORIZONTAL

SLIDER
1263
239
1441
272
max-root-distance_m
max-root-distance_m
0
8
8.0
0.5
1
m
HORIZONTAL

SLIDER
1263
278
1478
311
max-insect-distance_SC_m
max-insect-distance_SC_m
0
125
8.0
1
1
m
HORIZONTAL

SLIDER
1266
318
1479
351
max-insect-distance_PF_m
max-insect-distance_PF_m
0
500
8.0
1
1
m
HORIZONTAL

SLIDER
1266
359
1480
392
max-insect-distance_HN_m
max-insect-distance_HN_m
0
500
8.0
1
1
m
HORIZONTAL

SWITCH
1272
402
1521
435
background-prob-insect-dispersal?
background-prob-insect-dispersal?
0
1
-1000

SLIDER
1272
443
1524
476
prob-insect-disp_background
prob-insect-disp_background
0
0.05
0.015
0.0005
1
NIL
HORIZONTAL

CHOOSER
1273
609
1411
654
age-inf-suscept
age-inf-suscept
"age-indep" "age-dep-1" "age-dep-2"
0

SLIDER
1415
659
1573
692
dead-attr-duration_max
dead-attr-duration_max
0
4
0.0
1
1
yr
HORIZONTAL

SLIDER
1273
526
1482
559
dead-host-viability-duration_min
dead-host-viability-duration_min
0
2
0.0
1
1
yr
HORIZONTAL

SLIDER
1272
568
1485
601
dead-host-viability-duration_max
dead-host-viability-duration_max
0
2
0.0
1
1
yr
HORIZONTAL

SWITCH
1223
740
1333
773
new-attr?
new-attr?
0
1
-1000

SWITCH
8
343
148
376
export-tree-data?
export-tree-data?
1
1
-1000

SWITCH
10
426
131
459
export-rasters?
export-rasters?
1
1
-1000

SLIDER
7
382
155
415
tree-data-export-freq
tree-data-export-freq
0
500
300.0
1
1
NIL
HORIZONTAL

SLIDER
139
425
300
458
raster-export-frequency
raster-export-frequency
0
100
300.0
1
1
NIL
HORIZONTAL

BUTTON
277
60
340
93
NIL
setup
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
1229
15
1316
49
setup-base
no-display setup-base display
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
1589
25
1649
58
reset
no-display reset display
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
1409
20
1472
53
go
no-display go display
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
274
103
337
136
go
no-display go display
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
1619
868
1791
901
pct-set-aside
pct-set-aside
0
100
0.0
1
1
%
HORIZONTAL

SLIDER
9
462
148
495
pct-short-rotation
pct-short-rotation
0
100
50.0
1
1
%
HORIZONTAL

SLIDER
152
462
285
495
pct-long-rotation
pct-long-rotation
0
100
50.0
1
1
%
HORIZONTAL

SLIDER
16
585
120
618
replant-delay
replant-delay
0
5
1.0
1
1
yr
HORIZONTAL

CHOOSER
26
767
191
812
spread-infection-version?
spread-infection-version?
"none" "Go"
1

CHOOSER
14
676
191
721
initial-infection
initial-infection
"generate-initial-infection-pct" "generate-initial-infection-alt" "generate-initial-infection-center"
0

SLIDER
211
777
351
810
pct-initial-infection
pct-initial-infection
0
100
0.5
0.5
1
%
HORIZONTAL

SWITCH
146
629
288
662
reinitiate-infection?
reinitiate-infection?
1
1
-1000

SLIDER
17
508
201
541
rotation:short-rotation
rotation:short-rotation
0
37
37.0
1
1
yr
HORIZONTAL

SLIDER
19
545
197
578
rotation:long-rotation
rotation:long-rotation
0
80
80.0
1
1
yr
HORIZONTAL

BUTTON
1169
162
1250
195
p-disp-pf/hn
no-display\nlet the-max ( max [ prob-insect-disp_sum_PF ] of trees ) \nask trees [ ifelse prob-insect-disp_sum_PF = 0 [ set color white ]\n  [ set color scale-color red prob-insect-disp_sum_PF the-max 0 ]\n  if infected? [ set color blue ]\n]\ndisplay
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
1169
118
1250
151
p-disp-sc
no-display\nlet the-max ( max [ prob-insect-disp_sum_SC ] of trees ) \nask trees [ ifelse prob-insect-disp_sum_SC = 0 [ set color white ]\n  [ set color scale-color red prob-insect-disp_sum_SC the-max 0 ]\n  if infected? [ set color blue ]\n]\ndisplay
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
1169
203
1253
236
pi SC base
ask patches [ set pcolor color_patches ]\nlet minVal ( min [ prob-inf_sc_base ] of trees_current )\nlet maxVal ( max [ prob-inf_sc_base ] of trees_current )\nask cells   [ set  color scale-color red prob-inf_SC_base maxVal minVal ]\nprint (word \"min: \" minVal)\nprint (word \"max: \" maxVal)
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
1172
242
1254
275
col p-v-base
no-display ask patches [ set pcolor white ]\nask cells [ set color white\nif tree? [ set color grey\n  if prob-inf_SC_base  = 0  [ set color red ]\n  if prob-inf_PF_base  = 0  [ set color yellow ]\n  if prob-inf_HN_base  = 0  [ set color blue + 2 ]\n  if prob-inf_SC_base  = 0 AND prob-inf_PF_base = 0 [ set color orange + 1 ]\n  if prob-inf_SC_base  = 0 AND prob-inf_HN_base = 0 [ set color violet - 2]\n  if prob-inf_PF_base  = 0 AND prob-inf_HN_base = 0 [ set color green ]\n  if prob-inf_SC_base  = 0 AND prob-inf_PF_base = 0 AND prob-inf_HN_base = 0 [ set color black ]\n]\n]\ndisplay
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
1052
248
1137
281
p-rt-cntct
no-display\nask cells [ set color color_patches ]\nlet the-max ( max [ prob-root-contact_union ] of trees ) \nprint the-max\nask trees [ set color scale-color red prob-root-contact_union the-max 0\n]\ndisplay
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
1016
782
1106
815
fill-cell-grid
no-display\nask patches [\n  if not p_cell [\n    sprout-cells 1 [\n      set-default-cell-properties\n      if pxcor mod 2 = 0 [\n        set ycor ycor - 0.5\n      ]\n    ]\n  ]\n]\ndisplay
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
1476
20
1539
53
NIL
stop
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SWITCH
13
728
145
761
export-jsons?
export-jsons?
0
1
-1000

BUTTON
1084
615
1262
648
all-possible-potential-trees
ask cells [\n  if ( ( pxcor mod 2 = 0 AND pycor mod 2 = 0 ) OR \n       ( ( pxcor + 1 ) mod 2 = 0 AND ( pycor + 1 ) mod 2 = 0 ) OR \n       ( pxcor mod 3 = 0 AND pycor mod 3 = 0 ) OR \n       ( ( pxcor + 2 ) mod 3 = 0 AND ( pycor + 2 ) mod 3 = 0 )\n     ) [ set potential-tree? true set color 38 ]\n]\nset trees_potential ( cells with [ potential-tree? ] )
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SWITCH
21
629
137
662
spread-setup?
spread-setup?
0
1
-1000

SWITCH
150
728
314
761
export-potential-trees?
export-potential-trees?
0
1
-1000

SWITCH
1495
200
1626
233
debug-mode?
debug-mode?
1
1
-1000

BUTTON
1056
290
1131
323
in root rad
ask trees_potential [ \n  if length ( sentence inf-root-nbors_ids_r1 inf-root-nbors_ids_r2 inf-root-nbors_ids_r3 inf-root-nbors_ids_r4 ) > 0 [\n    set color sky\n  ]\n  if tree? and infected? [ set color red ]\n  if losing-infection? [ set color violet ]\n  if newly-infected? [ set color magenta ]\n]\n
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
1266
483
1646
516
insect-disp_background_stand-pct-inf-spillover-threshold
insect-disp_background_stand-pct-inf-spillover-threshold
0
100
6.0
1
1
%
HORIZONTAL

SWITCH
150
342
291
375
export-stand-data?
export-stand-data?
1
1
-1000

BUTTON
1172
363
1252
396
p-inf-root
no-display\nask cells [ set color color_patches ]\nlet the-max ( max [ prob-inf-root ] of trees ) \nprint the-max\nask trees [ set color scale-color red prob-root-contact_union the-max 0\n]\ndisplay
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
159
382
303
415
stand-data-export-freq
stand-data-export-freq
0
100
0.0
1
1
NIL
HORIZONTAL

SWITCH
1529
248
1634
281
thinning?
thinning?
0
1
-1000

BUTTON
1176
452
1256
485
w/attr-inf
no-display\nask trees_potential [\n  ifelse is-list? attr-inf [\n    if length attr-inf > 1 [ set color pink ]\n  ] [\n    if attr-inf > 1 [ set color pink ]\n  ]\n  if infected? [ set color red ]\n  if losing-infection? [ set color violet ]\n  if infected? and losing-infection? [ set color black ]\n]\ndisplay
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
1176
489
1256
522
w/attr-dead
ask trees_potential [\n if length attr-dead > 1 [ set color blue ]\n ;set color scale-color red ( max attr-dead ) 5 1\n if tree? and not alive? and mort-cause = \"bsrd\" and not new-bsrd-mortality? [\n   set color black\n ]\n]
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
1149
652
1218
685
new-inf
ask cells with [ newly-infected? ] [ set color magenta]\n
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
1053
659
1145
692
losing-inf
ask cells with [ losing-infection? ] [ set color violet]\n
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
1693
233
1788
266
cells-proj
ask cells-projecting [ set color blue ]
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
1053
697
1146
730
losing-dead
ask cells with [ losing-attr-dead? ] [ set color yellow - 1]\n
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
1149
689
1219
722
newly-dead
ask cells with [ new-bsrd-mortality? ] [ set color orange]\n
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
1686
182
1821
216
wInfRootNbors
ask trees_potential [\n  let inf-nbor-lists ( sentence inf-root-nbors_ages_r1 inf-root-nbors_ages_r2 inf-root-nbors_ages_r3 inf-root-nbors_ages_r4 )\n  if length inf-nbor-lists > 0 [\n    set color blue\n  ]\n  if infected? [ set color red ]\n]
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
1236
659
1408
692
dead-attr-duration_min
dead-attr-duration_min
0
4
0.0
1
1
yr
HORIZONTAL

BUTTON
1166
412
1256
445
col_attr-mgmt
ask trees_potential [\n  set color scale-color red ( max attr-mgmt ) 3.5 1\n]
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
1090
532
1266
565
NIL
project_root-infection_attraction
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
1089
452
1170
485
max attr-inf
let maxValue 0\nifelse attr-inf-dead_version = \"proportional\" [\n  set maxValue ( max [ attr-inf ] of trees_current )\n] [ set maxValue ( max attr-effect_inf-vs-live_list ) ]\nask trees_potential [\n  if attr-inf-dead_version = \"proportional\" [\n    set color scale-color red attr-inf maxValue 1\n  ]\n  if attr-inf-dead_version = \"mean\" [\n    set color scale-color red ( mean attr-inf ) maxValue 1\n  ]\n  if attr-inf-dead_version = \"max\" [\n    set color scale-color red ( mean attr-inf ) maxValue 1\n  ]\n  if infected? [ set color blue ]   \n]
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
1089
491
1171
524
max attr-dead
let maxValue ( max attr-effect_dead_list )\nask trees_potential [\n  set color scale-color red ( max attr-dead ) maxValue 1\n  ;set color scale-color red ( one-of attr-dead ) maxValue 1\n if tree? and not alive? and mort-cause = \"bsrd\" and not new-bsrd-mortality? [\n   set color violet\n ]\n]
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
1069
143
1142
176
  col_tsi
ask trees_current [\n  set color scale-color red time-since-inf 5 0\n  if not infected? [ set color sky ]\n]
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
1069
193
1142
226
col_tsm
ask trees_current [\n  set color scale-color red time-since-mort 5 0\n  if alive? [ set color sky ]\n  if alive? and time-since-mort > 0 [ set color violet ]\n]
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SWITCH
1599
80
1736
113
prt-exclude-data?
prt-exclude-data?
0
1
-1000

BUTTON
1505
602
1594
635
NIL
clear-plot
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SWITCH
1532
296
1636
329
harvest?
harvest?
0
1
-1000

SWITCH
1602
121
1825
154
prt-discount-by-contact-type?
prt-discount-by-contact-type?
0
1
-1000

SWITCH
1352
743
1455
776
attr-inf?
attr-inf?
0
1
-1000

SWITCH
1466
743
1580
776
attr-dead?
attr-dead?
0
1
-1000

SWITCH
1395
703
1590
736
attr-mgmt_conservative?
attr-mgmt_conservative?
0
1
-1000

SWITCH
1229
703
1386
736
prc_conservative?
prc_conservative?
0
1
-1000

SWITCH
400
823
503
856
R-gc?
R-gc?
1
1
-1000

SWITCH
24
822
198
855
R-clear-environment?
R-clear-environment?
0
1
-1000

SWITCH
209
823
393
856
set-Go-solution_direct?
set-Go-solution_direct?
1
1
-1000

SLIDER
1382
783
1576
816
inf-center-attr-threshold_pct
inf-center-attr-threshold_pct
0
100
30.0
1
1
%
HORIZONTAL

CHOOSER
1223
838
1384
883
inf-center-attr-threshold
inf-center-attr-threshold
"high" "medium" "low"
1

CHOOSER
1219
782
1362
827
attr-inf-dead_version
attr-inf-dead_version
"proportional" "median" "max"
0

BUTTON
1063
346
1138
379
color-by-HN
no-display\ncolor-by-HN\ndisplay
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
1043
580
1143
613
col_attr-source
ifelse attr-inf-dead_version = \"proportional\" [\n  let maxValue 20\n  ask trees_potential [\n    set color scale-color red attr-inf_source maxValue 0\n    if infected? [ set color blue ]   \n  ]\n  print \"---- attr-inf_source stats for current trees ----\"\n  print ( word \"Min: \"    ( min    [ attr-inf_source ] of trees_current ) )\n  print ( word \"Median: \" ( median [ attr-inf_source ] of trees_current ) )\n  print ( word \"Mean: \"   ( mean   [ attr-inf_source ] of trees_current ) )\n  print ( word \"Max: \"    ( max    [ attr-inf_source ] of trees_current ) )\n]\n[ user-message \"only works with attr-inf-dead_version = proportional\" ]\n
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
1023
839
1154
872
print attr-inf stats
print \"---- attr-inf stats for current trees ----\"\n  print ( word \"Min: \"    ( min    [ attr-inf ] of trees_current ) )\n  print ( word \"Median: \" ( median [ attr-inf ] of trees_current ) )\n  print ( word \"Mean: \"   ( mean   [ attr-inf ] of trees_current ) )\n  print ( word \"Max: \"    ( max    [ attr-inf ] of trees_current ) )
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

MONITOR
1066
58
1144
103
Age (mean)
mean [ age ] of trees_current
17
1
11

SWITCH
212
683
333
716
do-kill-runs?
do-kill-runs?
1
1
-1000

SWITCH
1663
400
1842
433
probability-correction?
probability-correction?
1
1
-1000

@#$#@#$#@
# Black-stain root disease landscape simulator model
#### Created by Adam J. Bouché under the guidance of Major Professors Dr. Klaus Puettmann and Dr. David Shaw

**Thesis in partial fulfillment of the requirements for a Master of Science in Forest Ecosystems and Society**
*Department of Forest Ecosystems and Society*
*College of Forestry*
*Oregon State University*
#####
# DELETE:::::
#### ASAP
  - RESOLVE:
    - prob-mort VS time-since-mort
    - prob-inf/age-inf-susc VS time-since-mort
    - infectiousness VS time-since-mort
    - attractiveness VS time-since-mort


#### KNOWN ISSUES TO IMPROVE:
  - For small stand sizes (1 stand), roads randomly distributed throughout stands making all areas attractive b/c of roads. "make-one-big-road" designed to fix this but isn't great
  - Some world sizes will be too small for a given number of stands and an error will be thrown. However, this will not happen at any size remotely resembling the real world.

#### TO CLEAN UP:
  - All other WIP code
  - All "check" "fix" "deactivate" and "delete" in code (walkthrough)
  - "make big road" buttons/sliders/code


### VERSION GUIDE:
#### v.1.10.0.0
Language changes in code:
  - "vector" => "insect"
  - "cumul" => "union" (only for probability of root contact)
  - "cumul" => "sum" (only for probability of insect dispersal)
  

#### v.1.9.3.0
  - Modify attr-inf, normalizing drawn value to be more conservative: data measured attractiveness of infected trees in clustered infection centers. Now, infection attraction can be scaled to infection centers.

#### v.1.9.2.0
  - Fixed probability of mortality list errors
  - Fixed time-since-infection counter functioning
  - Added dead-attr-duration min and max to add flexibility
    - Adjusted remove attr-mgmt accordingly
  - Combined project-root-infection loops (from 2 -> 1 )
  - Fixed the major error: previously I was setting the values for prob-vector-dispersal_cumul TO BE the outputs of the Go program. What I needed to do was modify the present prob-vector-dispersal_cumul values of each tree by adding the Go output.
  - Finally added attr-inf and attr-dead coming from new infections and new bsrd mortalities
  - Added flexibility to parameter values to allow for liberal/conservative estimates.

#### v.1.9.0.1
  - projecting root contact implemented using infected root neighbor age lists to avoid repeating projection from trees with sustained infection
  - projecting vector dispersal using the Go program fully implemented
  - ongoing cleanup work
  - changed "inf-cause" from a list to five booleans for "initial", "root", "sc", "pf", and "hn"
  - added "empty-agentset" to replace "nobody" for clearing agentsets to fix non-agentset error

#### v.1.8.3.3
  - clean model base (in case corrupted)
  - spread with a single loop as an alterative to see if it prevents the core dumps

#### v.1.8.3.0
  - project-infection_simpleFix to account for the fact that prob-root-contact_cumulative was not changing as trees got older
  - generate-initial-infection_pct was changed so that each management class had the same % of trees initially infected (because before, there were more infected trees in extensive because of the fact that extensive stands have higher planting density and therefore more trees)
  - project-infection_simpleFix_PROFILING - to chop it up and see what's happening in the profiler

#### v.1.8.2.0
Many changes since the last version.
  - project-infection w/ fix for changing prob root contact with age (was constant before)
  - resize world in different ways

#### v.1.8.1.1
Many changes since the last version.
  - project-infection
  - time-passes_NEW

##### NEW MODEL APPROACH - REDUCE MEMORY
  - implement newly-infected? using was* and became*
  - REDUCE overall number of CALCULATIONS (especially in R). This means not repeating so many times:
    - limit-prob-values (how to get around this?)
    - all other R-based
  - have cells that lose infection reduce infection of neighbors (project negative infection)
  - with age-bsrd-inf-susc, you can strip that completely out of the spread process and just use it at the end when calculating prob-inf
  - near the end, calculate SPECIFIC prob-infs first then add to get the total prob inf
  - make prob-exposure functions work
  - add a multiplier so that if dead the probability of infection decreases... add a multiplier to the age-inf-susc (although if this is the case, you would have to recalculate prob-infs as well by first changing susceptibility and then multiplying by exposure

#### v.1.8.0.1
  - Incorporating attraction
  - Removing old useless procedures (0.1)

#### v.1.7.1.1
  - Making the switch for new spread-infection
  - Reordering events
  - "spread-infection" changed from concatenating ring with ring # to calling coordinates from a global "ring-list"

#### v.1.7.1.0
  - this is where I begin renovating spread-infection for the last time, making it faster
  - repeate-use r variables declared in the global environment (so that r:clearLocal can be used)
  - Added prob-exposure properties, added susceptibility property
  - more reporters for specific groups of trees based on status (e.g., to-report root-infected-trees)

##### v.1.7.0.2
  - patch-based landscape brought into primary model code, old cell-based code moved to new supporting file for cell-based landscape
  - tiny tweaks for cluster runs
  - profile-go procedure to do exactly as it sounds

##### v.1.7.0.1
  - tiny tweaks for cluster runs

##### v.1.7.0.0
  - Fully implemented patch-based world with cells only used for trees (to reduce memory use and speed up spread)
  - General cleaning and organizing, removing redundant and outdated code, variables
  - Finally made "config" variables lower case
  - create-datetimestamp revamped to "YYYYMMDD_HHMMXM" format

##### v.1.6.2.3
Implementing patch-based world with cells only used for trees (to reduce memory use and speed up spread)
  - General cleaning and organizing, removing redundant and outdated code
  - Moved infection probability code to the file "probability_parameters.nls"
  - Removed "grow-seeds-alt"
  - 
  - 

##### v.1.6.2.2
Probably more changes than this, but:
  - Added "gc?" to turn of R/Java garbage collector that might be causing a GC error when running multiple runs in parallel in an experiment
  - Added "raster-export-frequency" to change OTF how often rasters can be exported (rather than have it hard coded)
  - Updated export-data-csv (data-file output) and progress-file output to make them more useful when combining data from different runs (including world size, etc)
  - Added more details to run file as well to help figure out what I'm looking at while testing the model

##### v.1.6.1_agentsets
  - Test of using agentsets to improve efficiency
  - Fixed the issue of management assignment to inappropriate management classes in "assign-mgmt" procedure
  - "max-ticks" changes to user-defined varb on a slider
  - removed "lowest-value" and "highest-value" redundant functions (min and max already exist...)
  - Move some functions to "basic_functions.nls"
  - Moved WIP to separate, included .nls files ("WIP_oldgrowth_development.nls", "under_development.nls")
  - Revamped "time-passes" procedure order and functions
  - Added "visualize?" switch to deactivate display and coloring when not necessary
  - Added "limit-ticks?" switch and "max-ticks" slider to make adjustment and use of tick limits easier for model runs

##### v.1.6.0
  - __includes added to keep ring globals in a separate .nls file
  - Begin implementing agentsets
  - "make-one-big-road"

##### v.1.5.6
  - spread-infection - Previously, calculations had been made with all surrounding cells rather than just with trees because of an error in parentheses'
  - changed "-alt" functions for root transmission to be primary

##### v.1.5.X
  - Lots of changes and fixes to parameters, and efforts to avoid issues with R by minimizing use of that extension
  - Always enough roads

##### v.1.5.0
  - Parameter reboot: roots
  - working on old-growth development (OG-development-3, 3-alt)

##### v.1.4.1
  - working on old-growth development (OG-development-3, 3-alt)

##### v.1.4.0
  - incorporation of ring globals to speed up spread processes (grow-seeds)
  - fixed road cell availability for landscapes with few edges

##### v.1.3.9
  - updates to profiler procedure
  - developing alternatives to ring file for spread processes

##### v.1.3.8
  - "old-growth" changed to "set-aside"
  - only live trees were able to become infected before
  - "limit-values" function to make sure probabilities are never >1
  - file outputs reorganized by landscape configuration
  - spread-infection improved to be faster (changed from "ask trees at-points with [...]" to "ask cells at-points [ if (...) ]")
  - create-datetimestamp updated to have full month code (to avoid confusing e.g., January, June, and July)

##### v.1.3.7
  - Added "attraction-test" as prototype for attraction
  - Added Behaviorspace tests to run on the cluster

##### v.1.3.5
  - Seed radius changed from user-defined to hard-coded
  - set-edges procedure removed, edge-width global removed (no longer necessary)

##### v.1.3.4
  - Infection model implemented (with old versions of parameters, need to be updated)
  - Stumps can only be infected up to 1 year from mortality

##### v.1.2.2
  -  Made efficient thinning, harvest, regeneration

##### v.1.2.1
  - Added thinnings
  - Fixed extensive stands so that stands with age > thinning age have trees removed
  - Fixed roads (for single-stand landscapes with no edges)

##### v.1.2.0
  - Added roads in the landscape at realistic proportions for each stand type
  - Changed management proportions from percents to whole numbers because of errors in summing to 1 with floating point math
  - Old-growth age fixed

### PIPE-DREAM EXTENSIONS

##### MANAGEMENT
  - Old-growth/set-aside age distribution (finding the right balance of inputs/outputs to maintain the age distribution)
  - More flexible thinning


## WHAT IS IT?

Developing a method for setting up landscapes composed of forest stands for a landscape-scale model.

## QUICK USE GUIDE
  1. Choose desired _**Configuration**_ using the drop-down menu:
  _**Random**_ – Randomly places stand “seeds” in the landscape and “grows” them outwards to create randomly distributed stands.

  2. Choose the dimensions of the world using the _**world-size**_ drop-down menu.
  _Note: The largest landscapes take some time to set up (potentially up to 5 minutes, depending on your computer), so it’s recommended to use the smaller world sizes for demonstration purposes._

  3. If generating a _**random**_ landscape: choose the number of stands with the _**n-stands**_ slider.
If generating a _**block**_ landscape: choose the number of stands with the _**n-blocks**_ drop down list.

  4. Choose the mean and standard deviation for generating a stand age distribution and tree spacing distribution (_**stand-age-mean**_, _**stand-age-sd**_, _**spacing-mean**_, _**spacing-sd**_).

  5. Click the _**Generate landscape**_ button.

  6. To generate stand edges, choose the cell width for edges by setting the _**edge-width**_ slider, then click the _**Set edges**_ button. The color of the edge cells will change to highlight the cells considered to be edges. To remove stand edges, use the _**Remove edges**_ button.

_Note: Coloration_

  - _**Set colors (ID)**_ - Once the landscape has been generated, stands are colored sequentially from the lightest to darkest shade of blue based on their stand ID number. If you change the color scheme (e.g., with the _**Set colors (age)**_ button), you can return to the stand-ID color scheme with the _**Set colors (ID) button**_.

  - _**Set colors (age)**_ - Stands are colored sequentially from the lightest to darkest shade of red based on their relative ages.


## HOW IT WORKS
Each landscape is composed of a grid of hexagonal “cells”, the smallest units in the model.

Both processes begin by clearing the landscape and generating a list of stand ID numbers based on the total number of stands to be generated in the landscape.

#### Stand property distributions

The model currently works as described below. It could easily be adjusted to set the ages and rotation lengths differently (e.g., at a user-defined proportion of the landscape).

A stand age distribution is generated based on an approximated random normal distribution with user-defined mean (Slider: _**stand-age-mean**_) and standard deviation (Slider: _**stand-age-sd**_). Stand ages are drawn from this distribution and assigned to each of the stands. 

A stand density distribution is generated based on an approximated random normal distribution is also generated for stand density in order to account for the difference in root transmission probability with changing distance between trees. The mean and standard deviation is set using the _**spacing-mean**_ and _**spacing-sd**_ sliders.

#### Random stands
  - _**Configuration**_ must be set to _Random_, and the other settings are chosen by the user.
  - When _**Generate landscape**_ is selected, the world is cleared, the world size is set, and the age and spacing distributions are generated. Then, an empty hexagon cell grid is produced.
  - The _**setup-seeds**_ procedure randomly selects _**n-stands**_ cells (as determined by the _**n-stands**_ slider) to serve as “seeds” for each stand. Each of these cells are assigned each a unique stand ID number and an age and spacing from the respective distributions.
  - The _**grow-seeds**_ procedure causes each cell with a stand ID number (starting with the “seeds”) to "grow" outward by assigning its properties to neighbors. Each stand’s expansion stops where another stand is reached. Stands can continue to grow until all cells have become part of a stand.
  - Set edge changes the "edge?" property to "true" for all patches within a certain, user-defined edge distance (edge-width slider). The color of these patches is changed to highlight their status as edges.

#### Block stands
  - _**Configuration**_ must be set to _Blocks_, and the other settings are chosen by the user.
  - When _**Generate landscape**_ is selected, the world is cleared, the world size is set, and the age and spacing distributions are generated. Then, an empty hexagon cell grid is produced.
  - This procedure then subdivides the landscape into evenly sized, square stands in a grid based on the number of blocks desired (_**n-blocks**_). This grid has sqrt(n-blocks) rows and sqrt(n-blocks) columns. Patches are then assigned

## Forest cover classes

#### A) SILVICULTURAL SYSTEM: Intensive
  - Manager class:	Industrial/investment
  - Land use:	Production
  - Forest type:	Plantation
  - Composition:	Douglas-fir monoculture
  - Age structure:	Even-aged
  - Rotation length:	30
  - Thinning?	No
  - Thinning age:	NA
  - Spacing (initial):	13
  - Density (initial):	257.75 tpa
  - Density (residual):	NA

#### B) SILVICULTURAL SYSTEM: Extensive
  - Manager class:	Small woodland owner
  - Land use:	Production
  - Forest type:	Plantation
  - Composition:	Douglas-fir monoculture
  - Age structure:	Even-aged
  - Rotation length:	80
  - Thinning?	Yes (3)
  - Thinning age:	15 (PCT), 35, 55
  - Spacing (initial):	10
  - Density (initial):	43560 tpa
  - Density (residual):	300, 160, 100 tpa


## FUTURE CHANGES & ADDITIONS
#### Planned additions
  - Create two age distributions (younger and older stands)
  - Assigning the rotation length
  - Create roads
  - Create "barriers" of old forest in the landscape
  - create vector attraction (another neighborset?)

## HOW TO USE IT

(how to use the model, including a description of each of the items in the Interface tab)

## THINGS TO NOTICE

(suggested things for the user to notice while running the model)

## THINGS TO TRY

(suggested things for the user to try to do (move sliders, switches, etc.) with the model)

## EXTENDING THE MODEL

  - Reimplment block stands
  - Isolate thinning vs. rotation age to see what the main driver of infection is

#### C) SILVICULTURAL SYSTEM: Passive
  - Management class: Set-aside
  - Land use: Conservation
  - Forest type: Old growth
  - Composition: Douglas-fir mixture
  - Age structure: Old growth
  - Thinning?	No
  - Thinning age:	NA
  - Spacing (initial):	____
  - Density (initial):	____
  - Density (residual):	NA


## NETLOGO FEATURES

(interesting or unusual features of NetLogo that the model uses, particularly in the Code tab; or where workarounds were needed for missing features)

## RELATED MODELS

(models in the NetLogo Models Library and elsewhere which are of related interest)


## CREDITS AND REFERENCES

Adam Bouché developed this model for his Master of Science thesis in the Department of Forest Ecosystems and Society in the College of Forestry at Oregon State University with major professors Dr. Klaus Puettmann and Dr. David Shaw.

This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License.

The correct citation for this model is:
Bouché, A. J., Klaus J. Puettmann, and David C. Shaw. Black stain root disease landscape model
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

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

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

hex
false
0
Polygon -7500403 true true 0 150 75 30 225 30 300 150 225 270 75 270

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
NetLogo 6.0.4
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="test_noThin" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>delete-program-directory
r:stop</final>
    <timeLimit steps="5"/>
    <enumeratedValueSet variable="comp">
      <value value="&quot;law&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-world-size?">
      <value value="&quot;automatic&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-stands">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stand-size_mean_acres">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pct-intensive">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pct-extensive">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="world-size">
      <value value="&quot;50 x 50&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="base">
      <value value="&quot;patch&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="config">
      <value value="&quot;random&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-blocks">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="limit-ticks?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-ticks">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="visualize?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="track-run?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="export-stand-data?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stand-data-export-freq">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="export-tree-data?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="tree-data-export-freq">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="export-rasters?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="raster-export-frequency">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="export-jsons?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="export-potential-trees?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-Go-solution_direct?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="R-clear-environment?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="R-gc?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug-mode?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pct-set-aside">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rotation:intensive">
      <value value="37"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rotation:extensive">
      <value value="80"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="thinning?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="harvest?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="replant-delay">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="spread-setup?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="spread-infection-version?">
      <value value="&quot;Go&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-infection">
      <value value="&quot;generate-initial-infection-pct&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pct-initial-infection">
      <value value="0.6125"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="reinitiate-infection?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-attraction-distance_m">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-root-distance_m">
      <value value="6.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-vector-distance_SC_m">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-vector-distance_PF_m">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-vector-distance_HN_m">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="background-prob-vector-dispersal?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prob-vector-disp_background">
      <value value="0.015"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="age-inf-suscept">
      <value value="&quot;age-indep&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sensitivity-analysis-param">
      <value value="&quot;none&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sensitivity-analysis-multiplier">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dead-host-viability-duration_min">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dead-host-viability-duration_max">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="new-attr?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="attr-inf?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="attr-dead?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prt-exclude-data?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="attr-inf-dead_version">
      <value value="&quot;proportional&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="attr-mgmt_conservative?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inf-center-attr-threshold">
      <value value="&quot;high&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prc_conservative?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prt-discount-by-contact-type?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vector-disp_background_stand-pct-inf-spillover-threshold">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dead-attr-duration_min">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dead-attr-duration_max">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inf-center-attr-threshold_pct">
      <value value="10"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="full_st-60_stSz-30_I-25_E-75_initInf-0.6125_intermediate_law" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>delete-program-directory
r:stop</final>
    <timeLimit steps="301"/>
    <enumeratedValueSet variable="comp">
      <value value="&quot;law&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-world-size?">
      <value value="&quot;automatic&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-stands">
      <value value="60"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stand-size_mean_acres">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pct-intensive">
      <value value="25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pct-extensive">
      <value value="75"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="world-size">
      <value value="&quot;50 x 50&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="base">
      <value value="&quot;patch&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="config">
      <value value="&quot;random&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-blocks">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="limit-ticks?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-ticks">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="visualize?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="track-run?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="export-stand-data?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stand-data-export-freq">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="export-tree-data?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="tree-data-export-freq">
      <value value="25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="export-rasters?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="raster-export-frequency">
      <value value="25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="export-jsons?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="export-potential-trees?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-Go-solution_direct?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="R-clear-environment?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="R-gc?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug-mode?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pct-set-aside">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rotation:intensive">
      <value value="37"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rotation:extensive">
      <value value="80"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="thinning?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="harvest?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="replant-delay">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="spread-setup?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="spread-infection-version?">
      <value value="&quot;Go&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-infection">
      <value value="&quot;generate-initial-infection-pct&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pct-initial-infection">
      <value value="0.6125"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="reinitiate-infection?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-attraction-distance_m">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-root-distance_m">
      <value value="6.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-vector-distance_SC_m">
      <value value="125"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-vector-distance_PF_m">
      <value value="350"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-vector-distance_HN_m">
      <value value="350"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="background-prob-vector-dispersal?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prob-vector-disp_background">
      <value value="0.015"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="age-inf-suscept">
      <value value="&quot;age-indep&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sensitivity-analysis-param">
      <value value="&quot;none&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sensitivity-analysis-multiplier">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dead-host-viability-duration_min">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dead-host-viability-duration_max">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="new-attr?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="attr-inf?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="attr-dead?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prt-exclude-data?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="attr-inf-dead_version">
      <value value="&quot;proportional&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="attr-mgmt_conservative?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inf-center-attr-threshold">
      <value value="&quot;high&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prc_conservative?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prt-discount-by-contact-type?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vector-disp_background_stand-pct-inf-spillover-threshold">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dead-attr-duration_min">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dead-attr-duration_max">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inf-center-attr-threshold_pct">
      <value value="10"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="full_st-60_stSz-30_I-50_E-50_initInf-0.6125_intermediate_noThin_law" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>delete-program-directory
r:stop</final>
    <timeLimit steps="301"/>
    <enumeratedValueSet variable="comp">
      <value value="&quot;law&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-world-size?">
      <value value="&quot;automatic&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-stands">
      <value value="60"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stand-size_mean_acres">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pct-intensive">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pct-extensive">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="world-size">
      <value value="&quot;50 x 50&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="base">
      <value value="&quot;patch&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="config">
      <value value="&quot;random&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-blocks">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="limit-ticks?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-ticks">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="visualize?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="track-run?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="export-stand-data?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stand-data-export-freq">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="export-tree-data?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="tree-data-export-freq">
      <value value="25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="export-rasters?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="raster-export-frequency">
      <value value="25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="export-jsons?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="export-potential-trees?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-Go-solution_direct?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="R-clear-environment?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="R-gc?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug-mode?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pct-set-aside">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rotation:intensive">
      <value value="37"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rotation:extensive">
      <value value="80"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="thinning?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="harvest?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="replant-delay">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="spread-setup?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="spread-infection-version?">
      <value value="&quot;Go&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-infection">
      <value value="&quot;generate-initial-infection-pct&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pct-initial-infection">
      <value value="0.6125"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="reinitiate-infection?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-attraction-distance_m">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-root-distance_m">
      <value value="6.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-vector-distance_SC_m">
      <value value="125"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-vector-distance_PF_m">
      <value value="350"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-vector-distance_HN_m">
      <value value="350"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="background-prob-vector-dispersal?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prob-vector-disp_background">
      <value value="0.015"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="age-inf-suscept">
      <value value="&quot;age-indep&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sensitivity-analysis-param">
      <value value="&quot;none&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sensitivity-analysis-multiplier">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dead-host-viability-duration_min">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dead-host-viability-duration_max">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="new-attr?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="attr-inf?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="attr-dead?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prt-exclude-data?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="attr-inf-dead_version">
      <value value="&quot;proportional&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="attr-mgmt_conservative?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inf-center-attr-threshold">
      <value value="&quot;high&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prc_conservative?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prt-discount-by-contact-type?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vector-disp_background_stand-pct-inf-spillover-threshold">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dead-attr-duration_min">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dead-attr-duration_max">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inf-center-attr-threshold_pct">
      <value value="10"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="full_st-60_stSz-30_I-0_E-100_initInf-0.6125_root_noThin_law" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>delete-program-directory
r:stop</final>
    <timeLimit steps="301"/>
    <enumeratedValueSet variable="comp">
      <value value="&quot;law&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-world-size?">
      <value value="&quot;automatic&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-stands">
      <value value="60"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stand-size_mean_acres">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pct-intensive">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pct-extensive">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="world-size">
      <value value="&quot;50 x 50&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="base">
      <value value="&quot;patch&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="config">
      <value value="&quot;random&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-blocks">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="limit-ticks?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-ticks">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="visualize?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="track-run?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="export-stand-data?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stand-data-export-freq">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="export-tree-data?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="tree-data-export-freq">
      <value value="25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="export-rasters?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="raster-export-frequency">
      <value value="25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="export-jsons?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="export-potential-trees?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-Go-solution_direct?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="R-clear-environment?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="R-gc?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug-mode?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pct-set-aside">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rotation:intensive">
      <value value="37"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rotation:extensive">
      <value value="80"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="thinning?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="harvest?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="replant-delay">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="spread-setup?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="spread-infection-version?">
      <value value="&quot;Go&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-infection">
      <value value="&quot;generate-initial-infection-pct&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pct-initial-infection">
      <value value="0.6125"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="reinitiate-infection?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-attraction-distance_m">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-root-distance_m">
      <value value="6.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-vector-distance_SC_m">
      <value value="125"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-vector-distance_PF_m">
      <value value="350"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-vector-distance_HN_m">
      <value value="350"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="background-prob-vector-dispersal?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prob-vector-disp_background">
      <value value="0.015"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="age-inf-suscept">
      <value value="&quot;age-indep&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sensitivity-analysis-param">
      <value value="&quot;none&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sensitivity-analysis-multiplier">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dead-host-viability-duration_min">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dead-host-viability-duration_max">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="new-attr?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="attr-inf?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="attr-dead?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prt-exclude-data?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="attr-inf-dead_version">
      <value value="&quot;proportional&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="attr-mgmt_conservative?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inf-center-attr-threshold">
      <value value="&quot;high&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prc_conservative?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prt-discount-by-contact-type?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vector-disp_background_stand-pct-inf-spillover-threshold">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dead-attr-duration_min">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dead-attr-duration_max">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inf-center-attr-threshold_pct">
      <value value="15"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="full_st-60_stSz-30_I-50_E-50_initInf-0.6125_root_noThin_law" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>delete-program-directory
r:stop</final>
    <timeLimit steps="301"/>
    <enumeratedValueSet variable="comp">
      <value value="&quot;law&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-world-size?">
      <value value="&quot;automatic&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-stands">
      <value value="60"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stand-size_mean_acres">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pct-intensive">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pct-extensive">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="world-size">
      <value value="&quot;50 x 50&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="base">
      <value value="&quot;patch&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="config">
      <value value="&quot;random&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-blocks">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="limit-ticks?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-ticks">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="visualize?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="track-run?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="export-stand-data?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stand-data-export-freq">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="export-tree-data?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="tree-data-export-freq">
      <value value="25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="export-rasters?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="raster-export-frequency">
      <value value="25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="export-jsons?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="export-potential-trees?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-Go-solution_direct?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="R-clear-environment?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="R-gc?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug-mode?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pct-set-aside">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rotation:intensive">
      <value value="37"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rotation:extensive">
      <value value="80"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="thinning?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="harvest?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="replant-delay">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="spread-setup?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="spread-infection-version?">
      <value value="&quot;Go&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-infection">
      <value value="&quot;generate-initial-infection-pct&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pct-initial-infection">
      <value value="0.6125"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="reinitiate-infection?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-attraction-distance_m">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-root-distance_m">
      <value value="6.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-vector-distance_SC_m">
      <value value="125"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-vector-distance_PF_m">
      <value value="350"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-vector-distance_HN_m">
      <value value="350"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="background-prob-vector-dispersal?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prob-vector-disp_background">
      <value value="0.015"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="age-inf-suscept">
      <value value="&quot;age-indep&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sensitivity-analysis-param">
      <value value="&quot;none&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sensitivity-analysis-multiplier">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dead-host-viability-duration_min">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dead-host-viability-duration_max">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="new-attr?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="attr-inf?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="attr-dead?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prt-exclude-data?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="attr-inf-dead_version">
      <value value="&quot;proportional&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="attr-mgmt_conservative?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inf-center-attr-threshold">
      <value value="&quot;high&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prc_conservative?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prt-discount-by-contact-type?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vector-disp_background_stand-pct-inf-spillover-threshold">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dead-attr-duration_min">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dead-attr-duration_max">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inf-center-attr-threshold_pct">
      <value value="15"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="full_st-60_stSz-30_I-0_E-100_initInf-0.6125_vector_noThin_law" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>delete-program-directory
r:stop</final>
    <timeLimit steps="301"/>
    <enumeratedValueSet variable="comp">
      <value value="&quot;law&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-world-size?">
      <value value="&quot;automatic&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-stands">
      <value value="60"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stand-size_mean_acres">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pct-intensive">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pct-extensive">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="world-size">
      <value value="&quot;50 x 50&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="base">
      <value value="&quot;patch&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="config">
      <value value="&quot;random&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-blocks">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="limit-ticks?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-ticks">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="visualize?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="track-run?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="export-stand-data?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stand-data-export-freq">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="export-tree-data?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="tree-data-export-freq">
      <value value="25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="export-rasters?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="raster-export-frequency">
      <value value="25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="export-jsons?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="export-potential-trees?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-Go-solution_direct?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="R-clear-environment?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="R-gc?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug-mode?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pct-set-aside">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rotation:intensive">
      <value value="37"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rotation:extensive">
      <value value="80"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="thinning?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="harvest?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="replant-delay">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="spread-setup?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="spread-infection-version?">
      <value value="&quot;Go&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-infection">
      <value value="&quot;generate-initial-infection-pct&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pct-initial-infection">
      <value value="0.6125"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="reinitiate-infection?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-attraction-distance_m">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-root-distance_m">
      <value value="6.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-vector-distance_SC_m">
      <value value="125"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-vector-distance_PF_m">
      <value value="350"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-vector-distance_HN_m">
      <value value="350"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="background-prob-vector-dispersal?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prob-vector-disp_background">
      <value value="0.015"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="age-inf-suscept">
      <value value="&quot;age-indep&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sensitivity-analysis-param">
      <value value="&quot;none&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sensitivity-analysis-multiplier">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dead-host-viability-duration_min">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dead-host-viability-duration_max">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="new-attr?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="attr-inf?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="attr-dead?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prt-exclude-data?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="attr-inf-dead_version">
      <value value="&quot;proportional&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="attr-mgmt_conservative?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inf-center-attr-threshold">
      <value value="&quot;high&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prc_conservative?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prt-discount-by-contact-type?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vector-disp_background_stand-pct-inf-spillover-threshold">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dead-attr-duration_min">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dead-attr-duration_max">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inf-center-attr-threshold_pct">
      <value value="3"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="full_st-60_stSz-30_I-50_E-50_initInf-0.6125_vector_noThin_law" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>delete-program-directory
r:stop</final>
    <timeLimit steps="301"/>
    <enumeratedValueSet variable="comp">
      <value value="&quot;law&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-world-size?">
      <value value="&quot;automatic&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-stands">
      <value value="60"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stand-size_mean_acres">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pct-intensive">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pct-extensive">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="world-size">
      <value value="&quot;50 x 50&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="base">
      <value value="&quot;patch&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="config">
      <value value="&quot;random&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-blocks">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="limit-ticks?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-ticks">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="visualize?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="track-run?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="export-stand-data?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stand-data-export-freq">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="export-tree-data?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="tree-data-export-freq">
      <value value="25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="export-rasters?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="raster-export-frequency">
      <value value="25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="export-jsons?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="export-potential-trees?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-Go-solution_direct?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="R-clear-environment?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="R-gc?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug-mode?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pct-set-aside">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rotation:intensive">
      <value value="37"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rotation:extensive">
      <value value="80"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="thinning?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="harvest?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="replant-delay">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="spread-setup?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="spread-infection-version?">
      <value value="&quot;Go&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-infection">
      <value value="&quot;generate-initial-infection-pct&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pct-initial-infection">
      <value value="0.6125"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="reinitiate-infection?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-attraction-distance_m">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-root-distance_m">
      <value value="6.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-vector-distance_SC_m">
      <value value="125"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-vector-distance_PF_m">
      <value value="350"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-vector-distance_HN_m">
      <value value="350"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="background-prob-vector-dispersal?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prob-vector-disp_background">
      <value value="0.015"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="age-inf-suscept">
      <value value="&quot;age-indep&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sensitivity-analysis-param">
      <value value="&quot;none&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sensitivity-analysis-multiplier">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dead-host-viability-duration_min">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dead-host-viability-duration_max">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="new-attr?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="attr-inf?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="attr-dead?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prt-exclude-data?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="attr-inf-dead_version">
      <value value="&quot;proportional&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="attr-mgmt_conservative?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inf-center-attr-threshold">
      <value value="&quot;high&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prc_conservative?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prt-discount-by-contact-type?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vector-disp_background_stand-pct-inf-spillover-threshold">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dead-attr-duration_min">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dead-attr-duration_max">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inf-center-attr-threshold_pct">
      <value value="3"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="full_st-60_stSz-30_I-25_E-75_initInf-0.6125_vector_law" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>delete-program-directory
r:stop</final>
    <timeLimit steps="301"/>
    <enumeratedValueSet variable="comp">
      <value value="&quot;law&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-world-size?">
      <value value="&quot;automatic&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-stands">
      <value value="60"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stand-size_mean_acres">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pct-intensive">
      <value value="25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pct-extensive">
      <value value="75"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="world-size">
      <value value="&quot;50 x 50&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="base">
      <value value="&quot;patch&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="config">
      <value value="&quot;random&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-blocks">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="limit-ticks?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-ticks">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="visualize?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="track-run?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="export-stand-data?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stand-data-export-freq">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="export-tree-data?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="tree-data-export-freq">
      <value value="25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="export-rasters?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="raster-export-frequency">
      <value value="25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="export-jsons?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="export-potential-trees?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-Go-solution_direct?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="R-clear-environment?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="R-gc?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug-mode?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pct-set-aside">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rotation:intensive">
      <value value="37"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rotation:extensive">
      <value value="80"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="thinning?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="harvest?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="replant-delay">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="spread-setup?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="spread-infection-version?">
      <value value="&quot;Go&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-infection">
      <value value="&quot;generate-initial-infection-pct&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pct-initial-infection">
      <value value="0.6125"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="reinitiate-infection?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-attraction-distance_m">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-root-distance_m">
      <value value="6.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-vector-distance_SC_m">
      <value value="125"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-vector-distance_PF_m">
      <value value="350"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-vector-distance_HN_m">
      <value value="350"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="background-prob-vector-dispersal?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prob-vector-disp_background">
      <value value="0.015"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="age-inf-suscept">
      <value value="&quot;age-indep&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sensitivity-analysis-param">
      <value value="&quot;none&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sensitivity-analysis-multiplier">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dead-host-viability-duration_min">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dead-host-viability-duration_max">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="new-attr?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="attr-inf?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="attr-dead?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prt-exclude-data?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="attr-inf-dead_version">
      <value value="&quot;proportional&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="attr-mgmt_conservative?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inf-center-attr-threshold">
      <value value="&quot;high&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prc_conservative?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prt-discount-by-contact-type?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vector-disp_background_stand-pct-inf-spillover-threshold">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dead-attr-duration_min">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dead-attr-duration_max">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inf-center-attr-threshold_pct">
      <value value="3"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="full_st-60_stSz-30_I-75_E-25_initInf-0.6125_vector_law" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>delete-program-directory
r:stop</final>
    <timeLimit steps="301"/>
    <enumeratedValueSet variable="comp">
      <value value="&quot;law&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-world-size?">
      <value value="&quot;automatic&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-stands">
      <value value="60"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stand-size_mean_acres">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pct-intensive">
      <value value="75"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pct-extensive">
      <value value="25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="world-size">
      <value value="&quot;50 x 50&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="base">
      <value value="&quot;patch&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="config">
      <value value="&quot;random&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-blocks">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="limit-ticks?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-ticks">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="visualize?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="track-run?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="export-stand-data?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stand-data-export-freq">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="export-tree-data?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="tree-data-export-freq">
      <value value="25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="export-rasters?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="raster-export-frequency">
      <value value="25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="export-jsons?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="export-potential-trees?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-Go-solution_direct?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="R-clear-environment?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="R-gc?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug-mode?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pct-set-aside">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rotation:intensive">
      <value value="37"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rotation:extensive">
      <value value="80"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="thinning?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="harvest?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="replant-delay">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="spread-setup?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="spread-infection-version?">
      <value value="&quot;Go&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-infection">
      <value value="&quot;generate-initial-infection-pct&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pct-initial-infection">
      <value value="0.6125"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="reinitiate-infection?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-attraction-distance_m">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-root-distance_m">
      <value value="6.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-vector-distance_SC_m">
      <value value="125"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-vector-distance_PF_m">
      <value value="350"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-vector-distance_HN_m">
      <value value="350"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="background-prob-vector-dispersal?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prob-vector-disp_background">
      <value value="0.015"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="age-inf-suscept">
      <value value="&quot;age-indep&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sensitivity-analysis-param">
      <value value="&quot;none&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sensitivity-analysis-multiplier">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dead-host-viability-duration_min">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dead-host-viability-duration_max">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="new-attr?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="attr-inf?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="attr-dead?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prt-exclude-data?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="attr-inf-dead_version">
      <value value="&quot;proportional&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="attr-mgmt_conservative?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inf-center-attr-threshold">
      <value value="&quot;high&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prc_conservative?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prt-discount-by-contact-type?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vector-disp_background_stand-pct-inf-spillover-threshold">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dead-attr-duration_min">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dead-attr-duration_max">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inf-center-attr-threshold_pct">
      <value value="3"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="full_st-60_stSz-30_I-25_E-75_initInf-0.6125_root_law" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>delete-program-directory
r:stop</final>
    <timeLimit steps="301"/>
    <enumeratedValueSet variable="comp">
      <value value="&quot;law&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-world-size?">
      <value value="&quot;automatic&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-stands">
      <value value="60"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stand-size_mean_acres">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pct-intensive">
      <value value="25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pct-extensive">
      <value value="75"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="world-size">
      <value value="&quot;50 x 50&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="base">
      <value value="&quot;patch&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="config">
      <value value="&quot;random&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-blocks">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="limit-ticks?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-ticks">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="visualize?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="track-run?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="export-stand-data?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stand-data-export-freq">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="export-tree-data?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="tree-data-export-freq">
      <value value="25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="export-rasters?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="raster-export-frequency">
      <value value="25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="export-jsons?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="export-potential-trees?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-Go-solution_direct?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="R-clear-environment?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="R-gc?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug-mode?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pct-set-aside">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rotation:intensive">
      <value value="37"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rotation:extensive">
      <value value="80"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="thinning?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="harvest?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="replant-delay">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="spread-setup?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="spread-infection-version?">
      <value value="&quot;Go&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-infection">
      <value value="&quot;generate-initial-infection-pct&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pct-initial-infection">
      <value value="0.6125"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="reinitiate-infection?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-attraction-distance_m">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-root-distance_m">
      <value value="6.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-vector-distance_SC_m">
      <value value="125"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-vector-distance_PF_m">
      <value value="350"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-vector-distance_HN_m">
      <value value="350"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="background-prob-vector-dispersal?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prob-vector-disp_background">
      <value value="0.015"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="age-inf-suscept">
      <value value="&quot;age-indep&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sensitivity-analysis-param">
      <value value="&quot;none&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sensitivity-analysis-multiplier">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dead-host-viability-duration_min">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dead-host-viability-duration_max">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="new-attr?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="attr-inf?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="attr-dead?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prt-exclude-data?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="attr-inf-dead_version">
      <value value="&quot;proportional&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="attr-mgmt_conservative?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inf-center-attr-threshold">
      <value value="&quot;high&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prc_conservative?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prt-discount-by-contact-type?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vector-disp_background_stand-pct-inf-spillover-threshold">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dead-attr-duration_min">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dead-attr-duration_max">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inf-center-attr-threshold_pct">
      <value value="15"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="full_st-60_stSz-30_I-75_E-25_initInf-0.6125_root_law" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>delete-program-directory
r:stop</final>
    <timeLimit steps="301"/>
    <enumeratedValueSet variable="comp">
      <value value="&quot;law&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-world-size?">
      <value value="&quot;automatic&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-stands">
      <value value="60"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stand-size_mean_acres">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pct-intensive">
      <value value="75"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pct-extensive">
      <value value="25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="world-size">
      <value value="&quot;50 x 50&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="base">
      <value value="&quot;patch&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="config">
      <value value="&quot;random&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-blocks">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="limit-ticks?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-ticks">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="visualize?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="track-run?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="export-stand-data?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stand-data-export-freq">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="export-tree-data?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="tree-data-export-freq">
      <value value="25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="export-rasters?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="raster-export-frequency">
      <value value="25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="export-jsons?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="export-potential-trees?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-Go-solution_direct?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="R-clear-environment?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="R-gc?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug-mode?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pct-set-aside">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rotation:intensive">
      <value value="37"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rotation:extensive">
      <value value="80"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="thinning?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="harvest?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="replant-delay">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="spread-setup?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="spread-infection-version?">
      <value value="&quot;Go&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-infection">
      <value value="&quot;generate-initial-infection-pct&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pct-initial-infection">
      <value value="0.6125"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="reinitiate-infection?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-attraction-distance_m">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-root-distance_m">
      <value value="6.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-vector-distance_SC_m">
      <value value="125"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-vector-distance_PF_m">
      <value value="350"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-vector-distance_HN_m">
      <value value="350"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="background-prob-vector-dispersal?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prob-vector-disp_background">
      <value value="0.015"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="age-inf-suscept">
      <value value="&quot;age-indep&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sensitivity-analysis-param">
      <value value="&quot;none&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sensitivity-analysis-multiplier">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dead-host-viability-duration_min">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dead-host-viability-duration_max">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="new-attr?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="attr-inf?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="attr-dead?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prt-exclude-data?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="attr-inf-dead_version">
      <value value="&quot;proportional&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="attr-mgmt_conservative?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inf-center-attr-threshold">
      <value value="&quot;high&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prc_conservative?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prt-discount-by-contact-type?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vector-disp_background_stand-pct-inf-spillover-threshold">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dead-attr-duration_min">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dead-attr-duration_max">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inf-center-attr-threshold_pct">
      <value value="15"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="full_st-60_stSz-30_I-0_E-100_initInf-0.6125_intermediate_noThin_law" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>delete-program-directory
r:stop</final>
    <timeLimit steps="301"/>
    <enumeratedValueSet variable="comp">
      <value value="&quot;law&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-world-size?">
      <value value="&quot;automatic&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-stands">
      <value value="60"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stand-size_mean_acres">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pct-intensive">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pct-extensive">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="world-size">
      <value value="&quot;50 x 50&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="base">
      <value value="&quot;patch&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="config">
      <value value="&quot;random&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-blocks">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="limit-ticks?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-ticks">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="visualize?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="track-run?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="export-stand-data?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stand-data-export-freq">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="export-tree-data?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="tree-data-export-freq">
      <value value="25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="export-rasters?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="raster-export-frequency">
      <value value="25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="export-jsons?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="export-potential-trees?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-Go-solution_direct?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="R-clear-environment?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="R-gc?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug-mode?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pct-set-aside">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rotation:intensive">
      <value value="37"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rotation:extensive">
      <value value="80"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="thinning?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="harvest?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="replant-delay">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="spread-setup?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="spread-infection-version?">
      <value value="&quot;Go&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-infection">
      <value value="&quot;generate-initial-infection-pct&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pct-initial-infection">
      <value value="0.6125"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="reinitiate-infection?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-attraction-distance_m">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-root-distance_m">
      <value value="6.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-vector-distance_SC_m">
      <value value="125"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-vector-distance_PF_m">
      <value value="350"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-vector-distance_HN_m">
      <value value="350"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="background-prob-vector-dispersal?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prob-vector-disp_background">
      <value value="0.015"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="age-inf-suscept">
      <value value="&quot;age-indep&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sensitivity-analysis-param">
      <value value="&quot;none&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sensitivity-analysis-multiplier">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dead-host-viability-duration_min">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dead-host-viability-duration_max">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="new-attr?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="attr-inf?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="attr-dead?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prt-exclude-data?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="attr-inf-dead_version">
      <value value="&quot;proportional&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="attr-mgmt_conservative?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inf-center-attr-threshold">
      <value value="&quot;high&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prc_conservative?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prt-discount-by-contact-type?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vector-disp_background_stand-pct-inf-spillover-threshold">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dead-attr-duration_min">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dead-attr-duration_max">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inf-center-attr-threshold_pct">
      <value value="10"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="full_st-60_stSz-30_I-75_E-25_initInf-0.6125_intermediate_law" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>delete-program-directory
r:stop</final>
    <timeLimit steps="301"/>
    <enumeratedValueSet variable="comp">
      <value value="&quot;law&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-world-size?">
      <value value="&quot;automatic&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-stands">
      <value value="60"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stand-size_mean_acres">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pct-intensive">
      <value value="75"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pct-extensive">
      <value value="25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="world-size">
      <value value="&quot;50 x 50&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="base">
      <value value="&quot;patch&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="config">
      <value value="&quot;random&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-blocks">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="limit-ticks?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-ticks">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="visualize?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="track-run?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="export-stand-data?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stand-data-export-freq">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="export-tree-data?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="tree-data-export-freq">
      <value value="25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="export-rasters?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="raster-export-frequency">
      <value value="25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="export-jsons?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="export-potential-trees?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-Go-solution_direct?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="R-clear-environment?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="R-gc?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug-mode?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pct-set-aside">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rotation:intensive">
      <value value="37"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rotation:extensive">
      <value value="80"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="thinning?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="harvest?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="replant-delay">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="spread-setup?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="spread-infection-version?">
      <value value="&quot;Go&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-infection">
      <value value="&quot;generate-initial-infection-pct&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pct-initial-infection">
      <value value="0.6125"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="reinitiate-infection?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-attraction-distance_m">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-root-distance_m">
      <value value="6.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-vector-distance_SC_m">
      <value value="125"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-vector-distance_PF_m">
      <value value="350"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-vector-distance_HN_m">
      <value value="350"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="background-prob-vector-dispersal?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prob-vector-disp_background">
      <value value="0.015"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="age-inf-suscept">
      <value value="&quot;age-indep&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sensitivity-analysis-param">
      <value value="&quot;none&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sensitivity-analysis-multiplier">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dead-host-viability-duration_min">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dead-host-viability-duration_max">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="new-attr?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="attr-inf?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="attr-dead?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prt-exclude-data?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="attr-inf-dead_version">
      <value value="&quot;proportional&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="attr-mgmt_conservative?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inf-center-attr-threshold">
      <value value="&quot;high&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prc_conservative?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prt-discount-by-contact-type?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vector-disp_background_stand-pct-inf-spillover-threshold">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dead-attr-duration_min">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dead-attr-duration_max">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inf-center-attr-threshold_pct">
      <value value="10"/>
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

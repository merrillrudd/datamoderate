#V3.30.13.00-trans;_2019_03_09;_Stock_Synthesis_by_Richard_Methot_(NOAA)_using_ADMB_12.0
#Stock Synthesis (SS) is a work of the U.S. Government and is not subject to copyright protection in the United States.
#Foreign copyrights may apply. See copyright.txt for more information.
#_user_support_available_at:NMFS.Stock.Synthesis@noaa.gov
#_user_info_available_at:https://vlab.ncep.noaa.gov/group/stock-synthesis
#_data_and_control_files: codOM.dat // codOM.ctl
0  # 0 means do not read wtatage.ss; 1 means read and use wtatage.ss and also read and use growth parameters
1  #_N_Growth_Patterns
1 #_N_platoons_Within_GrowthPattern
#_Cond 1 #_Morph_between/within_stdev_ratio (no read if N_morphs=1)
#_Cond  1 #vector_Morphdist_(-1_in_first_val_gives_normal_approx)
#
2 # recr_dist_method for parameters:  2=main effects for GP, Settle timing, Area; 3=each Settle entity; 4=none, only when N_GP*Nsettle*pop==1
1 # not yet implemented; Future usage: Spawner-Recruitment: 1=global; 2=by area
1 #  number of recruitment settlement assignments
0 # unused option
#GPattern month  area  age (for each settlement assignment)
 1 1 1 0
#
#_Cond 0 # N_movement_definitions goes here if Nareas > 1
#_Cond 1.0 # first age that moves (real age at begin of season, not integer) also cond on do_migration>0
#_Cond 1 1 1 2 4 10 # example move definition for seas=1, morph=1, source=1 dest=2, age1=4, age2=10
#
0 #_Nblock_Patterns
#_blocks_per_pattern
# begin and end years of blocks
#
# controls for all timevary parameters
1 #_env/block/dev_adjust_method for all time-vary parms (1=warn relative to base parm bounds; 3=no bound check)
#
# AUTOGEN
0 0 0 0 0 # autogen: 1st element for biology, 2nd for SR, 3rd for Q, 4th reserved, 5th for selex
# where: 0 = autogen all time-varying parms; 1 = read each time-varying parm line; 2 = read then autogen if parm min==-12345
#
#_Available timevary codes
#_Block types: 0: P_block=P_base*exp(TVP); 1: P_block=P_base+TVP; 2: P_block=TVP; 3: P_block=P_block(-1) + TVP
#_Block_trends: -1: trend bounded by base parm min-max and parms in transformed units (beware); -2: endtrend and infl_year direct values; -3: end and infl as fraction of base range
#_EnvLinks:  1: P(y)=P_base*exp(TVP*env(y));  2: P(y)=P_base+TVP*env(y);  3: null;  4: P(y)=2.0/(1.0+exp(-TVP1*env(y) - TVP2))
#_DevLinks:  1: P(y)*=exp(dev(y)*dev_se;  2: P(y)+=dev(y)*dev_se;  3: random walk;  4: zero-reverting random walk with rho;  21-24 keep last dev for rest of years
#
#
#
# setup for M, growth, maturity, fecundity, recruitment distibution, movement
#
0 #_natM_type:_0=1Parm; 1=N_breakpoints;_2=Lorenzen;_3=agespecific;_4=agespec_withseasinterpolate
  #_no additional input for selected M option; read 1P per morph
#
1 # GrowthModel: 1=vonBert with L1&L2; 2=Richards with L1&L2; 3=age_specific_K_incr; 4=age_specific_K_decr; 5=age_specific_K_each; 6=NA; 7=NA; 8=growth cessation
0 #_Age(post-settlement)_for_L1;linear growth below this
999 #_Growth_Age_for_L2 (999 to use as Linf)
-999 #_exponential decay for growth above maxage (value should approx initial Z; -999 replicates 3.24; -998 to not allow growth above maxage)
0  #_placeholder for future growth feature
#
0 #_SD_add_to_LAA (set to 0.1 for SS2 V1.x compatibility)
0 #_CV_Growth_Pattern:  0 CV=f(LAA); 1 CV=F(A); 2 SD=F(LAA); 3 SD=F(A); 4 logSD=F(A)
#
1 #_maturity_option:  1=length logistic; 2=age logistic; 3=read age-maturity matrix by growth_pattern; 4=read age-fecundity; 5=disabled; 6=read length-maturity
0 #_First_Mature_Age
1 #_fecundity option:(1)eggs=Wt*(a+b*Wt);(2)eggs=a*L^b;(3)eggs=a*Wt^b; (4)eggs=a+b*L; (5)eggs=a+b*W
0 #_hermaphroditism option:  0=none; 1=female-to-male age-specific fxn; -1=male-to-female age-specific fxn
1 #_parameter_offset_approach (1=none, 2= M, G, CV_G as offset from female-GP1, 3=like SS2 V1.x)
#
#_growth_parms
#_ LO HI INIT PRIOR PR_SD PR_type PHASE env_var&link dev_link dev_minyr dev_maxyr dev_PH Block Block_Fxn
# Sex: 1  BioPattern: 1  NatMort
 0.01 1.8 0.18 -1.715 0.5 3 -3 0 0 0 0 0 0 0 # NatM_p_1_Fem_GP_1
# Sex: 1  BioPattern: 1  Growth
 1 80 14.25 14.25 0.2 0 -2 0 0 0 0 0 0 0 # L_at_Amin_Fem_GP_1
 25 250 55 55 0.2 0 -5 0 0 0 0 0 0 0 # L_at_Amax_Fem_GP_1
 0.01 2 0.30 0.30 0.8 0 -2 0 0 0 0 0 0 0 # VonBert_K_Fem_GP_1
 -0.01 0.5 0.1 0.1 0.8 0 -3 0 0 0 0 0 0 0 # CV_young_Fem_GP_1
 0.01 0.5 0.1 0.1 0.8 0 -5 0 0 0 0 0 0 0 # CV_old_Fem_GP_1
# Sex: 1  BioPattern: 1  WtLen
 0 3 6.8e-06 6.8e-06 0 0 -1 0 0 0 0 0 0 0 # Wtlen_1_Fem
 2.5 3.5 3.101 3.101 0.2 0 -3 0 0 0 0 0 0 0 # Wtlen_2_Fem
# Sex: 1  BioPattern: 1  Maturity&Fecundity
 10 50 36.3 36.3 0 0 -3 0 0 0 0 0 0 0 # Mat50%_Fem
 -2 2 -0.276 0 0 0 -3 0 0 0 0 0 0 0 # Mat_slope_Fem
 -3 3 1 0 0 0 -3 0 0 0 0 0 0 0 # Eggs/kg_inter_Fem
 -3 4 0 0 0 0 -3 0 0 0 0 0 0 0 # Eggs/kg_slope_wt_Fem
# Hermaphroditism
#  Recruitment Distribution
 -4 4 0 0 0 0 -4 0 0 0 0 0 0 0 # RecrDist_GP_1
 -4 4 0 0 0 0 -4 0 0 0 0 0 0 0 # RecrDist_Area_1
 -4 4 0 0 0 0 -4 0 0 0 0 0 0 0 # RecrDist_timing_1
#  Cohort growth dev base
 -4 4 1 0 0 0 -4 0 0 0 0 0 0 0 # CohortGrowDev
#  Movement
#  Age Error from parameters
#  catch multiplier
#  fraction female, by GP
 0.000001 0.999999 0.5 0.5  0.5 0 -99 0 0 0 0 0 0 0 # FracFemale_GP_1
#
#_no timevary MG parameters
#
#_seasonal_effects_on_biology_parms
 0 0 0 0 0 0 0 0 0 0 #_femwtlen1,femwtlen2,mat1,mat2,fec1,fec2,Malewtlen1,malewtlen2,L1,K
#_ LO HI INIT PRIOR PR_SD PR_type PHASE
#_Cond -2 2 0 0 -1 99 -2 #_placeholder when no seasonal MG parameters
#
3 #_Spawner-Recruitment; Options: 2=Ricker; 3=std_B-H; 4=SCAA; 5=Hockey; 6=B-H_flattop; 7=survival_3Parm; 8=Shepherd_3Parm; 9=RickerPower_3parm
0  # 0/1 to use steepness in initial equ recruitment calculation
0  #  future feature:  0/1 to make realized sigmaR a function of SR curvature
#_          LO            HI          INIT         PRIOR         PR_SD       PR_type      PHASE    env-var    use_dev   dev_mnyr   dev_mxyr     dev_PH      Block    Blk_Fxn #  parm_name
            5            20          10          10            10             0          1          0          0          0          0          0          0          0 # SR_LN(R0)
           0.2             1          0.7           0.7          0.05             0         -4          0          0          0          0          0          0          0 # SR_BH_steep
             0             2           0.8           0.8           0.8             0         -5          0          0          0          0          0          0          0 # SR_sigmaR
            -5             5             0             0             1             0         -4          0          0          0          0          0          0          0 # SR_regime
             0             0             0             0             0             0        -99          0          0          0          0          0          0          0 # SR_autocorr
#_no timevary SR parameters
1 #do_recdev:  0=none; 1=devvector (R=F(SSB)+dev); 2=deviations (R=F(SSB)+dev); 3=deviations (R=R0*dev; dev2=R-f(SSB)); 4=like 3 with sum(dev2) adding penalty
1 # first year of main recr_devs; early devs can preceed this era
100 # last year of main recr_devs; forecast devs start in following year
-2 #_recdev phase
1 # (0/1) to read 13 advanced options
 -29 #_recdev_early_start (0=none; neg value makes relative to recdev_start)
 -4 #_recdev_early_phase
 0 #_forecast_recruitment phase (incl. late recr) (0 value resets to maxphase+1)
 1 #_lambda for Fcast_recr_like occurring before endyr+1
 1 #_last_yr_nobias_adj_in_MPD; begin of ramp
 1 #_first_yr_fullbias_adj_in_MPD; begin of plateau
 100 #_last_yr_fullbias_adj_in_MPD
 100 #_end_yr_for_ramp_in_MPD (can be in forecast to shape ramp, but SS sets bias_adj to 0.0 for fcast yrs)
 1 #_max_bias_adj_in_MPD (-1 to override ramp and set biasadj=1.0 for all estimated recdevs)
 0 #_period of cycles in recruitment (N parms read below)
 -5 #min rec_dev
 5 #max rec_dev
140 #_read_recdevs
-29 -0.44838051724177
-28 -0.184141991586624
-27 1.2469666513193
-26 0.0564067131396608
-25 0.103430188128757
-24 1.37205198950662
-23 0.368732964791362
-22 -1.01204898768523
-21 -0.549482281514821
-20 -0.356529576079966
-19 0.979265437951569
-18 0.287851061645891
-17 0.320617160475242
-16 0.0885461727560958
-15 -0.44467290780326
-14 1.42953050944246
-13 0.398280382583392
-12 -1.57329372530371
-11 0.561084721250948
-10 -0.378233126182347
-9 -0.854258964789476
-8 -0.174379931726636
-7 -0.820803558645792
-6 -0.583112983432912
-5 -0.500031414279405
-4 -1.34935464859393
-3 0.67022963559562
-2 0.122698494269212
-1 -0.910509549609558
0 1.00305193685594
1 0.341171377181451
2 -0.236057186393817
3 0.716100528836018
4 0.702506790026434
5 0.65726486530999
6 0.550912203280073
7 0.443134122830071
8 -0.0495293684613773
9 -0.244770130991933
10 -0.304376800809906
11 -0.55576558313641
12 -0.166333822415679
13 -1.01231708125461
14 1.73516477227081
15 0.966369598643992
16 -0.898486866562679
17 -0.322307868239261
18 -0.373324282898575
19 0.623972094669054
20 -0.0666952531774634
21 0.202654811195804
22 -0.0228374042789624
23 -0.0342963658330529
24 1.09488182721157
25 -0.180616788527414
26 1.21317648354363
27 -1.23900224338418
28 0.467690999708855
29 0.099083395075691
30 0.172753254995178
31 0.303711586207906
32 -0.401858762487442
33 -0.266565906935536
34 -0.814860306485671
35 -0.857432981180462
36 0.242822913123406
37 0.358567822903541
38 0.0424033813844033
39 0.73781397430379
40 1.64006774850172
41 -0.392824932845228
42 -1.84733510051265
43 0.804590819569805
44 -0.567360610065914
45 -0.550406893173886
46 0.820457095757359
47 -0.227818405640807
48 -0.976574169803629
49 0.14504278379932
50 -0.111113089951236
51 0.00461134871990955
52 0.308224320901064
53 -0.296528025433927
54 0.515501238815066
55 -0.176389249455001
56 0.265425571132558
57 0.877471210519478
58 0.348145192667042
59 -0.260745268424981
60 0.919046094760875
61 0.794803084769696
62 0.438717567606456
63 0.190985388089153
64 -0.502324860831497
65 1.08852195882401
66 -0.480207669717702
67 1.74986639441326
68 1.22608850094815
69 -0.188560287280382
70 -0.821136720245425
71 -0.568325250959441
72 0.205506967325224
73 -0.197353502769899
74 -0.278034079518187
75 -0.761294853812013
76 -0.0360221798471362
77 -0.627923575565661
78 -1.33435354927051
79 -0.30418121623021
80 0.735197287248613
81 -0.460277570086714
82 0.486371457780027
83 -1.29430616663133
84 -0.0444495724196316
85 0.41552576315477
86 0.240922689733372
87 0.0845409553191547
88 -0.512564806644301
89 -0.679763476826866
90 -0.819303032483931
91 0.0941172776801007
92 -0.757979691347842
93 -0.392445954960535
94 -0.204873753758598
95 1.47508960418577
96 -0.521559921356367
97 0.188309257827886
98 0.0623686796509686
99 -0.769485307304103
100 -0.0570464688988789
101 1.15564068673868
102 0.361203242463372
103 0.0329863375943519
104 -0.3379974658717
105 -1.64259777723241
106 0.90506977073134
107 -1.16851205673986
108 0.591958008701868
109 1.52728285537399
110 -1.15511452877744
#Fishing Mortality info
0.3 # F ballpark
-2001 # F ballpark year (neg value to disable)
2 # F_Method:  1=Pope; 2=instan. F; 3=hybrid (hybrid is recommended)
4  # max F or harvest rate, depends on F_Method
#
#
#
0 1 100 # overall start F value; overall phase; N detailed inputs to read
#Fleet Yr Seas F_value se phase (for detailed setup of F_Method=2; -Yr to fill remaining years)
1 1 1 0.0017 0.1 1
1 2 1 0.00198333333333333 0.1 1
1 3 1 0.00226666666666667 0.1 1
1 4 1 0.00255 0.1 1
1 5 1 0.00283333333333333 0.1 1
1 6 1 0.00311666666666667 0.1 1
1 7 1 0.0034 0.1 1
1 8 1 0.00368333333333333 0.1 1
1 9 1 0.00396666666666667 0.1 1
1 10 1 0.00425 0.1 1
1 11 1 0.00453333333333333 0.1 1
1 12 1 0.00481666666666667 0.1 1
1 13 1 0.0051 0.1 1
1 14 1 0.00538333333333333 0.1 1
1 15 1 0.00566666666666667 0.1 1
1 16 1 0.00595 0.1 1
1 17 1 0.00623333333333333 0.1 1
1 18 1 0.00651666666666667 0.1 1
1 19 1 0.0068 0.1 1
1 20 1 0.00708333333333333 0.1 1
1 21 1 0.00736666666666667 0.1 1
1 22 1 0.00765 0.1 1
1 23 1 0.00793333333333333 0.1 1
1 24 1 0.00821666666666667 0.1 1
1 25 1 0.0085 0.1 1
1 26 1 0.0085 0.1 1
1 27 1 0.0199310344827586 0.1 1
1 28 1 0.0313620689655172 0.1 1
1 29 1 0.0427931034482759 0.1 1
1 30 1 0.0542241379310345 0.1 1
1 31 1 0.0656551724137931 0.1 1
1 32 1 0.0770862068965517 0.1 1
1 33 1 0.0885172413793104 0.1 1
1 34 1 0.099948275862069 0.1 1
1 35 1 0.111379310344828 0.1 1
1 36 1 0.122810344827586 0.1 1
1 37 1 0.134241379310345 0.1 1
1 38 1 0.145672413793103 0.1 1
1 39 1 0.157103448275862 0.1 1
1 40 1 0.168534482758621 0.1 1
1 41 1 0.179965517241379 0.1 1
1 42 1 0.191396551724138 0.1 1
1 43 1 0.202827586206897 0.1 1
1 44 1 0.214258620689655 0.1 1
1 45 1 0.225689655172414 0.1 1
1 46 1 0.237120689655172 0.1 1
1 47 1 0.248551724137931 0.1 1
1 48 1 0.25998275862069 0.1 1
1 49 1 0.271413793103448 0.1 1
1 50 1 0.282844827586207 0.1 1
1 51 1 0.294275862068965 0.1 1
1 52 1 0.305706896551724 0.1 1
1 53 1 0.317137931034483 0.1 1
1 54 1 0.328568965517241 0.1 1
1 55 1 0.34 0.1 1
1 56 1 0.34 0.1 1
1 57 1 0.34 0.1 1
1 58 1 0.34 0.1 1
1 59 1 0.34 0.1 1
1 60 1 0.34 0.1 1
1 61 1 0.34 0.1 1
1 62 1 0.327473684210526 0.1 1
1 63 1 0.314947368421053 0.1 1
1 64 1 0.302421052631579 0.1 1
1 65 1 0.289894736842105 0.1 1
1 66 1 0.277368421052632 0.1 1
1 67 1 0.264842105263158 0.1 1
1 68 1 0.252315789473684 0.1 1
1 69 1 0.239789473684211 0.1 1
1 70 1 0.227263157894737 0.1 1
1 71 1 0.214736842105263 0.1 1
1 72 1 0.20221052631579 0.1 1
1 73 1 0.189684210526316 0.1 1
1 74 1 0.177157894736842 0.1 1
1 75 1 0.164631578947368 0.1 1
1 76 1 0.152105263157895 0.1 1
1 77 1 0.139578947368421 0.1 1
1 78 1 0.127052631578947 0.1 1
1 79 1 0.114526315789474 0.1 1
1 80 1 0.102 0.1 1
1 81 1 0.102 0.1 1
1 82 1 0.102 0.1 1
1 83 1 0.102 0.1 1
1 84 1 0.102 0.1 1
1 85 1 0.102 0.1 1
1 86 1 0.102 0.1 1
1 87 1 0.102 0.1 1
1 88 1 0.102 0.1 1
1 89 1 0.102 0.1 1
1 90 1 0.102 0.1 1
1 91 1 0.102 0.1 1
1 92 1 0.102 0.1 1
1 93 1 0.102 0.1 1
1 94 1 0.102 0.1 1
1 95 1 0.102 0.1 1
1 96 1 0.102 0.1 1
1 97 1 0.102 0.1 1
1 98 1 0.102 0.1 1
1 99 1 0.102 0.1 1
1 100 1 0.102 0.1 1
#_Q_setup for fleets with cpue or survey data
#_1:  fleet number
#_2:  link type: (1=simple q, 1 parm; 2=mirror simple q, 1 mirrored parm; 3=q and power, 2 parm; 4=mirror with offset, 2 parm)
#_3:  extra input for link, i.e. mirror fleet# or dev index number
#_4:  0/1 to select extra sd parameter
#_5:  0/1 for biasadj or not
#_6:  0/1 to float
#_   fleet      link link_info  extra_se   biasadj     float  #  fleetname
         1         1         0         0         0         0  #  Fishery
         2         1         0         0         0         0  #  Survey
-9999 0 0 0 0 0
#
#_Q_parms(if_any);Qunits_are_ln(q)
#_          LO            HI          INIT         PRIOR         PR_SD       PR_type      PHASE    env-var    use_dev   dev_mnyr   dev_mxyr     dev_PH      Block    Blk_Fxn  #  parm_name
            -3             3             0             0            99             0         -5          0          0          0          0          0          0          0  #  LnQ_base_Fishery(1)
            -3             3             0             0            99             0         -5          0          0          0          0          0          0          0  #  LnQ_base_Survey(2)
#_no timevary Q parameters
#
#_size_selex_patterns
#Pattern:_0; parm=0; selex=1.0 for all sizes
#Pattern:_1; parm=2; logistic; with 95% width specification
#Pattern:_5; parm=2; mirror another size selex; PARMS pick the min-max bin to mirror
#Pattern:_15; parm=0; mirror another age or length selex
#Pattern:_6; parm=2+special; non-parm len selex
#Pattern:_43; parm=2+special+2;  like 6, with 2 additional param for scaling (average over bin range)
#Pattern:_8; parm=8; New doublelogistic with smooth transitions and constant above Linf option
#Pattern:_9; parm=6; simple 4-parm double logistic with starting length; parm 5 is first length; parm 6=1 does desc as offset
#Pattern:_21; parm=2+special; non-parm len selex, read as pairs of size, then selex
#Pattern:_22; parm=4; double_normal as in CASAL
#Pattern:_23; parm=6; double_normal where final value is directly equal to sp(6) so can be >1.0
#Pattern:_24; parm=6; double_normal with sel(minL) and sel(maxL), using joiners
#Pattern:_25; parm=3; exponential-logistic in size
#Pattern:_27; parm=3+special; cubic spline
#Pattern:_42; parm=2+special+3; // like 27, with 2 additional param for scaling (average over bin range)
#_discard_options:_0=none;_1=define_retention;_2=retention&mortality;_3=all_discarded_dead;_4=define_dome-shaped_retention
#_Pattern Discard Male Special
 24 0 0 0 # 1 Fishery
 24 0 0 0 # 2 Survey
#
#_age_selex_patterns
#Pattern:_0; parm=0; selex=1.0 for ages 0 to maxage
#Pattern:_10; parm=0; selex=1.0 for ages 1 to maxage
#Pattern:_11; parm=2; selex=1.0  for specified min-max age
#Pattern:_12; parm=2; age logistic
#Pattern:_13; parm=8; age double logistic
#Pattern:_14; parm=nages+1; age empirical
#Pattern:_15; parm=0; mirror another age or length selex
#Pattern:_16; parm=2; Coleraine - Gaussian
#Pattern:_17; parm=nages+1; empirical as random walk  N parameters to read can be overridden by setting special to non-zero
#Pattern:_41; parm=2+nages+1; // like 17, with 2 additional param for scaling (average over bin range)
#Pattern:_18; parm=8; double logistic - smooth transition
#Pattern:_19; parm=6; simple 4-parm double logistic with starting age
#Pattern:_20; parm=6; double_normal,using joiners
#Pattern:_26; parm=3; exponential-logistic in age
#Pattern:_27; parm=3+special; cubic spline in age
#Pattern:_42; parm=2+special+3; // cubic spline; with 2 additional param for scaling (average over bin range)
#_Pattern Discard Male Special
 11 0 0 0 # 1 Fishery
 11 0 0 0 # 2 Survey
#
#_          LO            HI          INIT         PRIOR         PR_SD       PR_type      PHASE    env-var    use_dev   dev_mnyr   dev_mxyr     dev_PH      Block    Blk_Fxn  #  parm_name
# 1   Fishery LenSelex
            20           199          42          42          0.05             0          2          0          0          0          0        0.5          0          0  #  SizeSel_P1_Fishery(1)
            -5             3            -3            -3          0.05             0        -99          0          0          0          0        0.5          0          0  #  SizeSel_P2_Fishery(1)
            -4            12           5.1           5.1          0.05             0          3          0          0          0          0        0.5          0          0  #  SizeSel_P3_Fishery(1)
            -2            16            15            15          0.05             0        -99          0          0          0          0        0.5          0          0  #  SizeSel_P4_Fishery(1)
           -15             5          -999          -999          0.05             0        -99          0          0          0          0        0.5          0          0  #  SizeSel_P5_Fishery(1)
            -5          1000           999           999          0.05             0        -99          0          0          0          0        0.5          0          0  #  SizeSel_P6_Fishery(1)
# 2   Survey LenSelex
            20           199          42          42          0.05             0          2          0          0          0          0        0.5          0          0  #  SizeSel_P1_Fishery(1)
            -5             3            -3            -3          0.05             0        -99          0          0          0          0        0.5          0          0  #  SizeSel_P2_Fishery(1)
            -4            12           5.1           5.1          0.05             0          3          0          0          0          0        0.5          0          0  #  SizeSel_P3_Fishery(1)
            -2            16            15            15          0.05             0        -99          0          0          0          0        0.5          0          0  #  SizeSel_P4_Fishery(1)
           -15             5          -999          -999          0.05             0        -99          0          0          0          0        0.5          0          0  #  SizeSel_P5_Fishery(1)
            -5          1000           999           999          0.05             0        -99          0          0          0          0        0.5          0          0  #  SizeSel_P6_Fishery(1)
# 1   Fishery AgeSelex
             0             1           0.1           0.1            99             0         -3          0          0          0          0        0.5          0          0  #  AgeSel_P1_Fishery(1)
             0           101           100           100            99             0         -3          0          0          0          0        0.5          0          0  #  AgeSel_P2_Fishery(1)
# 2   Survey AgeSelex
             0             1           0.1           0.1            99             0         -3          0          0          0          0        0.5          0          0  #  AgeSel_P1_Survey(2)
             0           101           100           100            99             0         -3          0          0          0          0        0.5          0          0  #  AgeSel_P2_Survey(2)
#_no timevary selex parameters
#
0   #  use 2D_AR1 selectivity(0/1):  experimental feature
#_no 2D_AR1 selex offset used
#
# Tag loss and Tag reporting parameters go next
0  # TG_custom:  0=no read; 1=read if tags exist
#_Cond -6 6 1 1 2 0.01 -4 0 0 0 0 0 0 0  #_placeholder if no parameters
#
# no timevary parameters
#
#
# Input variance adjustments factors:
 #_1=add_to_survey_CV
 #_2=add_to_discard_stddev
 #_3=add_to_bodywt_CV
 #_4=mult_by_lencomp_N
 #_5=mult_by_agecomp_N
 #_6=mult_by_size-at-age_N
 #_7=mult_by_generalized_sizecomp
#_Factor  Fleet  Value
 -9999   1    0  # terminator
#
4 #_maxlambdaphase
1 #_sd_offset; must be 1 if any growthCV, sigmaR, or survey extraSD is an estimated parameter
# read 0 changes to default Lambdas (default value is 1.0)
# Like_comp codes:  1=surv; 2=disc; 3=mnwt; 4=length; 5=age; 6=SizeFreq; 7=sizeage; 8=catch; 9=init_equ_catch;
# 10=recrdev; 11=parm_prior; 12=parm_dev; 13=CrashPen; 14=Morphcomp; 15=Tag-comp; 16=Tag-negbin; 17=F_ballpark; 18=initEQregime
#like_comp fleet  phase  value  sizefreq_method
-9999  1  1  1  1  #  terminator
#
# lambdas (for info only; columns are phases)
#  0 0 0 0 #_CPUE/survey:_1
#  1 1 1 1 #_CPUE/survey:_2
#  1 1 1 1 #_CPUE/survey:_3
#  1 1 1 1 #_lencomp:_1
#  1 1 1 1 #_lencomp:_2
#  0 0 0 0 #_lencomp:_3
#  1 1 1 1 #_agecomp:_1
#  1 1 1 1 #_agecomp:_2
#  0 0 0 0 #_agecomp:_3
#  1 1 1 1 #_init_equ_catch
#  1 1 1 1 #_recruitments
#  1 1 1 1 #_parameter-priors
#  1 1 1 1 #_parameter-dev-vectors
#  1 1 1 1 #_crashPenLambda
#  0 0 0 0 # F_ballpark_lambda
0 # (0/1) read specs for more stddev reporting
 # 0 0 0 0 0 0 0 0 0 # placeholder for # selex_fleet, 1=len/2=age/3=both, year, N selex bins, 0 or Growth pattern, N growth ages, 0 or NatAge_area(-1 for all), NatAge_yr, N Natages
 # placeholder for vector of selex bins to be reported
 # placeholder for vector of growth ages to be reported
 # placeholder for vector of NatAges ages to be reported
999


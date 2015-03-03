;;RENAME TO "DEVICE"_ROUTINES I.E. D3D_ROUTINES AND RENAME FILE ACCORDINGLY
PRO templete_routines,inputs,grid,$     ;;INPUT: INPUTS AND GRID POINTS DO NOT CHANGE
					   nbi,$ 			;;OUTPUT: NEUTRAL BEAM INJECTION INFO STRUCTURE
					   chords,$ 		;;OUTPUT: CHORDS INFO STRUCTURE
					   profiles,$		;;OUTPUT: PROFILES STRUCTURE
					   equil,$			;;OUTPUT: MAGNETIC GRID STRUCTURE
					   err				;;OUTPUT: ERROR STATUS ERR=1 == SOMETHING WENT WRONG


	;;IN THIS SECTION YOU CAN USE WHATEVER ROUTINES
	;;YOU WANT SO LONG AS YOU DEFINE THE OUTPUT STRUCTURES
	;;CONTAIN AT LEAST THE FOLLOWING TAGS

	;;	IDL> help,chords
	;;	** Structure <1d447c48>, 11 tags, length=728, data length=724, refs=1:
	;;	   NCHAN           LONG                11
	;;	   DIAG            STRING    'OBLIQUE'
    ;;     CHAN_ID         LONG      Array[11]
	;;	   ULOS            DOUBLE    Array[11]
	;;	   VLOS            DOUBLE    Array[11]
	;;	   WLOS            DOUBLE    Array[11]
	;;	   ULENS           DOUBLE    Array[11]
	;;	   VLENS           DOUBLE    Array[11]
	;;	   WLENS           DOUBLE    Array[11]
	;;	   SIGMA_PI_RATIO  DOUBLE    Array[11]
	;;	   RA              FLOAT     Array[11]
	;;	   RD              FLOAT     Array[11]
	;;	   H               FLOAT     Array[11]

	;;	IDL> help,equil
	;;	** Structure <1d474638>, 10 tags, length=6636160, data length=6636138, refs=1:
	;;	   RHO2D           FLOAT     Array[40, 60]
	;;	   BR              DOUBLE    Array[40, 60]
	;;	   BPHI            DOUBLE    Array[40, 60]
	;;	   BW              DOUBLE    Array[40, 60]
	;;	   ER              DOUBLE    Array[40, 60]
	;;	   EW              DOUBLE    Array[40, 60]
	;;	   ERR             INT              0

	;;	IDL> help,profiles
	;;	** Structure <1d475698>, 7 tags, length=5816, data length=5810, refs=1:
	;;	   RHO             DOUBLE    Array[121]
	;;	   TE              DOUBLE    Array[121]
	;;	   TI              DOUBLE    Array[121]
	;;	   VTOR            DOUBLE    Array[121]
	;;	   DENE            DOUBLE    Array[121]
	;;	   ZEFF            DOUBLE    Array[121]
	;;	   ERR             INT              0

	;;	IDL> help,nbi
	;;	** Structure <1d475af8>, 13 tags, length=168, data length=168, refs=1:
	;;	   EINJ            DOUBLE           80.775734
	;;	   PINJ            DOUBLE           2.4117758
	;;	   FULL            DOUBLE          0.54850105
	;;	   HALF            DOUBLE          0.28972649
	;;	   THIRD           DOUBLE          0.16177245
	;;	   UVW_SRC         DOUBLE    Array[3]
	;;	   UVW_POS         DOUBLE    Array[3]
	;;	   BMWIDRA         DOUBLE           6.0000000
	;;	   BMWIDZA         DOUBLE           24.000000
	;;	   DIVY            DOUBLE    Array[3]
	;;	   DIVZ            DOUBLE    Array[3]
	;;	   FOCY            DOUBLE           999999.90
	;;	   FOCZ            DOUBLE           1000.0000

	;;FOR CONVINIENCE HERE ARE THE MINIMUM STRUCTURE DEFINITIONS
	equil={rho2d:rho2d,$    	   			;;INTERP. GRID IN MAGNETIC FLUX COORDINATES (RHO)
		   br:br,$					   		;;R MAGNETIC FIELD COMPONENT AT GRID POINTS
		   bphi:bphi,$					   	;;PHI MAGNETIC FIELD COMPONENT AT GRID POINTS
		   bw:bw,$					   		;;W MAGNETIC FIELD COMPONENT AT GRID POINTS
		   er:er,$							;;R ELECTRIC FIELD COMPONENT AT GRID POINTS
		   ew:ew }							;;W ELECTRIC FIELD COMPONENT AT GRID POINTS

	nbi={einj:einj,$				   		;;BEAM INJECTION ENERGY [keV]
		 pinj:pinj,$				   		;;BEAM INJECTION POWER  [MW]
		 full:full,$				   		;;FULL BEAM FRACTION
		 half:half,$				   		;;HALF BEAM FRACTION
		 third:third,$				   		;;THIRD BEAM FRACTION
		 uvw_src:uvw_src,$			   		;;POSITION OF BEAM SOURCE IN MACHINE COORDINATES [cm]
		 uvw_pos:uvw_pos,$			   		;;BEAM CROSSOVER POINT IN MACHINE COORDINATES [cm]
		 bmwidra:bmwidra,$			   		;;HORIZONTAL BEAM WIDTH [cm]
		 bmwidza:mbwidza,$			   		;;VERTICAL BEAM WIDTH   [cm]
		 focy:focy,$				   		;;HORIZONTAL FOCAL LENGTH [cm]
	     focz:focz,$						;;VERTICAL FOCAL LENGTH [cm]
		 divy:divy,$				   		;;HORIZONTAL BEAM DIVERGENCE [rad]
		 divz:divz }				   		;;VERTICAL BEAM DIVERGENCE [rad]

  	chords={sigma_pi_ratio:sigma_pi_ratio,$	;;RATIO OF SIGMA LINES TO PI LINES  (0 IF NPA)
		 nchan:nchan,$				  		;;NUMBER OF CHANNELS
         chan_id:chan_id,$                  ;;CHANNEL ID (0 FOR FIDA,1 FOR NPA)
		 umid:umid,$						;;U POS. OF WHERE CHORD CROSSES MIDPLANE [cm]
		 vmid:vmid,$						;;V POS. OF WHERE CHORD CROSSES MIDPLANE [cm]
         wmid:wmid,$						;;W POS. OF WHERE CHORD CROSSES MIDPLANE [cm]
		 ulens:ulens,$						;;U POS. OF LENS/APERTURE [cm]
		 vlens:vlens,$						;;V POS. OF LENS/APERTURE [cm]
		 wlens:wlens,$						;;W POS. OF LENS/APERTURE [cm]
		 ra:ra,$				            ;;RADIUS OF APERTURE [cm] (0 IF FIDA)
		 rd:rd,$				            ;;RADIUS OF DETECTOR [cm] (0 IF FIDA)
		 h:h}		                        ;;SEPERATION BETWEEN DETECTOR AND APERTURE [cm] (0 IF FIDA)

	profiles={time:time,$					;;SHOT TIME
			  rho:rho,$						;;RHO VALUES
			  ti:ti,$						;;ION TEMPERATURE [eV]
			  omega:omega,$					;;TORODIAL ANGULAR VELOCITY [rad/s]
			  te:te,$						;;ELECTRON TEMPERATURE [eV]
			  dene:dene,$					;;ELECTRON DENSITY [m^-3]
			  zeff:zeff}					;;ZEFF
END

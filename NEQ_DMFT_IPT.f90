MODULE NEQ_DMFT_IPT
  USE NEQ_CONTOUR
  USE NEQ_CONTOUR_GF
  USE NEQ_INPUT_VARS
  USE ELECTRIC_FIELD
  USE THERMOSTAT
  USE NEQ_IPT
  USE NEQ_AUX_FUNX
  USE NEQ_MEASURE

  !KELDYSH CONTOUR PARAMETERS 
  !=========================================================  
  type(kb_contour_params),public        :: cc_params 

END MODULE NEQ_DMFT_IPT

MODULE NEQ_VARS_GLOBAL
  USE NEQ_CONTOUR
  implicit none

  !OTHER VARIABLES THAT MAY EVENTUALLY BE NEEDED IN THE FUTURE
  !SHOULD GO HERE. CFR. WITH THE ED CODE CASE 

  !KELDYSH CONTOUR PARAMETERS 
  !=========================================================  
  type(kb_contour_params)                :: cc_params 


  !ELECTRIC FIELD
  !=========================================================  
  real(8),dimension(2)                   :: Ak,Ek         !Electric field vector potential and vector

end module NEQ_VARS_GLOBAL


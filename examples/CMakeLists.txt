#  
# Include and configure the channelflow examples
#
# This file is a part of channelflow version 2.0, https://channelflow.ch .
# License is GNU GPL version 2 or later: ./LICENSE
#
set(example_APPS
  couette
  writestuff
  channel
   evaluateSecondDeriv
  evaluateVectorfield
evaluate_second_order_term	
test_linear_operator
SolveEq
  )

foreach (program ${example_APPS})
    install_channelflow_application(${program} examples)
endforeach ()

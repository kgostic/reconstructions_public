for(yy in 2000:2017){
  for(bb in 1998:2017){
    if(bb>yy){
      break
    }
    print(c(yy,bb))
  }
}

b_yr = 2000
i_yr = 2005
nn = length(b_yr:min(b_yr+17, i_yr))
for(exp1 in 1:nn){
  mm = nn - exp1
  if(mm>0){
    for(exp2 in (1+exp1):(mm+exp1)){
      print(c(exp1,exp2))
    }
  }
}

master.g1p_00 == master.g2p_11
round(master.g1p_00,10) == round(master.g2p_11,10)

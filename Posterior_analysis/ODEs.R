Cisplatinvivo<- function(timepoint, state, parameters) {
  with(as.list(c(state, parameters)), {
    dPlasma = KidneyPt * k_KidPlas - Plasma * (k_PlasKid + k_ePlas);
    dKidneyPt = AccuPt * k_AccuKid + Plasma * k_PlasKid - KidneyPt * (k_KidPlas + k_KidAccu + k_eKid);
    dAccuPt = KidneyPt * k_KidAccu - AccuPt * k_AccuKid;
    dDD = k_kidDD*scale*(KidneyPt + AccuPt)/(hillDD + scale*(KidneyPt + AccuPt)) - d_DD*DD;
    dCD = (maxdeath*DD^h1/(k_hillnec^h1 + DD^h1))*(1-CD) - d_CD*CD*INF;
    dINF =  k_CDINF*CD - d_INF*INF;   
    dKF =  (k_INFKF*INF^(p+1)/(hillKF^(p+1) + INF^(p+1))) * (1-KF);
    dp = -d_p * p;
    list(c(dPlasma, dKidneyPt, dAccuPt, dDD, dCD, dINF, dKF,dp))
  })
}

VivoPK2 = function(timepoint,state,parameters){
  with(as.list(c(state, parameters)),{ 
    
    dPlasma= KidneyPt * k_KidPlas - Plasma * (k_PlasKid + k_ePlas)
    dKidneyPt = Plasma * k_PlasKid - KidneyPt * (k_KidPlas + k_eKid)
    list(c(dPlasma,dKidneyPt))})}

VivoPK3 <- function(timepoint, state, parameters) {
  with(as.list(c(state, parameters)), {
    dPlasma = KidneyPt * k_KidPlas - Plasma * (k_PlasKid + k_ePlas);
    dKidneyPt = AccuPt * k_AccuKid + Plasma * k_PlasKid - KidneyPt * (k_KidPlas + k_KidAccu + k_eKid);
    dAccuPt = KidneyPt * k_KidAccu - AccuPt * k_AccuKid;
    list(c(dPlasma, dKidneyPt, dAccuPt))
  })
}

VitroPK <- function(timepoint, state, parameters) {
  with(as.list(c(state, parameters)), {
    dQapi = Ncell * (Fouta * Qcell/(Ncell*Vcell) - Fina * Qapi/Vapi - Kmet * Qcell/(Ncell*Vcell)) # change in time in amount cisplatin in apical medium in ug/s
    dQbas = Ncell * (Foutb * Qcell/(Ncell*Vcell) - Finb * Qbas/Vbas); # change in time in amount cisplatin in basolateral medium in ug/s
    dQcell = Ncell * (Fina * Qapi/Vapi + Finb * Qbas/Vbas - (Fouta + Foutb) * Qcell/(Ncell*Vcell)) # change in amount of cisplatin in cells in ug/s
    list(c(dQapi, dQbas, dQcell))#
  })
}

Model7 <- function(timepoint, state, parameters) {
  with(as.list(c(state, parameters)), {
    dQapi = Ncell* (Fouta * Qcell/(Ncell*Vcell) - Fina * Qapi/Vapi- Kmet * Qcell/(Ncell*Vcell)) # change in time in amount cisplatin in apical medium in ug/s
    dQbas = Ncell * (Foutb * Qcell/(Ncell*Vcell) - Finb * Qbas/Vbas); # change in time in amount cisplatin in basolateral medium in ug/s
    dQcell = Ncell * (Fina * Qapi/Vapi + Finb * Qbas/Vbas - (Fouta + Foutb) * Qcell/(Ncell*Vcell)) # change in amount of cisplatin in cells in ug/s
    dQinter = k_inter * (Qcell - Qinter)
    dDD = - degDD * DD + (k_cDD * Qinter)/(hillDD + Qinter)
    dNEC = (maxdeath * DD^h/(k_hillnec^h + DD ^h))*(p - NEC) #
    list(c(dQapi, dQbas, dQcell, dQinter, dDD, dNEC))#
  })
}

setModelParameters_vivo <- function(model, draws_cis, k) {
  # Initialize empty lists for parameters and initial state
  pars <- list()
  inistate <- list()
  
  pars <- c(k_kidDD = draws_cis[[k,"k_kidDD"]],
          d_DD = draws_cis[[k,"d_DD"]],
            maxdeath= draws_cis[[k,"maxdeath"]],
            hillDD = draws_cis[[k,"hillDD"]],
            k_hillnec = draws_cis[[k,"k_hillnec"]],
            d_CD = draws_cis[[k,"d_CD"]],
            h1 = draws_cis[[k,"h1"]],
            k_CDINF = draws_cis[[k,"k_CDINF"]],
            d_INF = draws_cis[[k,"d_INF"]],
            k_INFKF = draws_cis[[k,"k_INFKF"]],
            hillKF = draws_cis[[k,"hillKF"]],
            d_p = 0.005,
            scale = PK_pars[["scale", "mean"]],
            k_PlasKid = PK_pars[["k_PlasKid", "mean"]],
            k_ePlas = PK_pars[["k_ePlas", "mean"]],
            k_KidPlas = PK_pars[["k_KidPlas", "mean"]],
            k_eKid = PK_pars[["k_eKid", "mean"]],
            k_KidAccu = PK_pars[["k_KidAccu", "mean"]],
            k_AccuKid = PK_pars[["k_AccuKid", "mean"]])
   inistate <- c(Plasma = PK_pars[["Plasma0", "mean"]], KidneyPt = 0, AccuPt = 0, 
              DD = 0, CD = 0, INF = draws_cis[[k,"INF0"]], KF = 0, p = 10)
   return(list(pars = pars, inistate = inistate))}

setModelParameters_PK3 = function(model, draws, k){    
  pars <- c(scale = draws[[k,"scale"]],
          k_PlasKid = draws[[k,"k_PlasKid"]],
          k_ePlas = draws[[k,"k_ePlas"]],
          k_KidPlas = draws[[k,"k_KidPlas"]],
          k_eKid = draws[[k,"k_eKid"]],
          k_KidAccu = draws[[k,"k_KidAccu"]],
          k_AccuKid = draws[[k,"k_AccuKid"]])

 
# Return a list containing both the parameters and initial state
return(list(pars = pars))
}
    
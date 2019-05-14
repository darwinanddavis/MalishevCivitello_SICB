## DEB IBM 

### Files required  
# "DEB_IBM.R"
# "DEB_INF_GUTS_IBM.nlogo"
# "FullStarve_shrink_dilute_damage3.Rda"
# "IndividualModel_IBM3.c"
# "IndividualModel_IBM3.so" # Mac OSX. generated from C
# "IndividualModel_IBM3.o" # Mac OSX. generated from C  
# "IndividualModel_IBM.dll" # Windows. generated from C  

# version 1.3 (post sicb ver with molluscicide)
# DEB_INF_GUTS_IBM_1.1.nlogo
# IndividualModel_IBM3.c (generates .o and .so files)  
# FullStarve_shrink_dilute_damage3.Rda

# 18-4-19 (v. 1.2)
# fixed snail control 

#################################  Running NetLogo in Mac ##################################
# if using Mac OSX El Capitan + 
# if using mac, set gui = F to run model in RStudio. 
# if you want to view the netlogo model in mac (with gui), set mac_gui = 1 

mac_gui <- 0

if(mac_gui==1){
  install.packages('JGR',,'http://www.rforge.net/')
  library(JGR)
  JGR::JGR()
}

# run java test  
install.packages("RCurl"); require(RCurl)
script <- getURL("https://raw.githubusercontent.com/darwinanddavis/SchistoIBM/master/mac/java_test.R", ssl.verifypeer = FALSE)
eval(parse(text = script))
# check rJava version  
.jinit()
.jcall("java/lang/System", "S", "getProperty", "java.runtime.version")
# get latest Java/Oracle version: https://www.oracle.com/technetwork/java/javase/downloads/index-jsp-138363.html  

# :three: [GCC compiler in R (unconfirmed)](https://stackoverflow.com/questions/1616983/building-r-packages-using-alternate-gcc)
# [Running Netlogo 6.0.+](https://github.com/NetLogo/NetLogo/issues/1282)

####################################  set user inputs ####################################### 
# Search "@netlogo" for netlogo code in file  

# set user outputs
mac <- 0 # mac or windows system? 1 = mac, 0 = windows 
gui <- 0 # display the gui? 1 = yes, 0 = no
pck <- 1 # if not already, install rnetlogo and rjava from source? 1 = yes, 0 = already installed 
save_to_file <- 0 # 1 = save simulation outputs to local dir, 0 = plot in current R session
mcmcplot <- 0 # 1 = save mcmc plots to dir
traceplot <- 0 # 1 = include traceplot in mcmc plots? intensive!!!

# set dir paths (for "/" for both Windows and Mac)
wd <- "<path to working dir>" # set working directory  
ver_nl <-"<netlogo version, eg. 6.0.4>" # type in Netlogo version. found in local dir. 
ver_gcc <-"<gcc version, e.g. 4.6.3>" # NULL # type in gcc version (if known). leave as "NULL" if unknown   
nl.path <- "<path to nelogo program>" # set path to Netlogo program location

# define starting conditions for simulation model @netlogo
n.ticks <- 120 # set number of days to simulate
day <- 1 # number of days to run simulation  
resource_type <- "detritus" # set resource type as "algae" or "detritus"
resources <- "event" # set resources: "cyclical" or "event"

####################################  set model paths #######################################
# load files
deb_samps <- "FullStarve_shrink_dilute_damage3.Rda"
deb_compile <- "IndividualModel_IBM3"

setwd(wd)
nl.model <- list.files(pattern="*.nlogo") ;nl.model # Netlogo model
if(mac==1){
  nl.path <- paste0(nl.path,"/NetLogo ",ver_nl,"/Java/"); cat("Mac path:",nl.path)
}else{
  nl.path <- paste0(nl.path,"/NetLogo ",ver_nl,"/app/"); cat("Windows path:",nl.path)
}
model.path <- paste0(wd,"/"); model.path # set path to Netlogo model   

####################################  load packages #######################################
# if already loaded, uninstall RNetlogo and rJava
if(pck==1){
  p<-c("rJava", "RNetLogo"); remove.packages(p)
  # then install rJava and RNetLogo from source
  if(mac==1){
    install.packages("rJava", repos = "https://cran.r-project.org/", type="source"); library(rJava)
    install.packages("RNetLogo", repos = "https://cran.r-project.org/", type="source"); library(RNetLogo)
  }
}
require(RNetLogo)
require(rJava)

# check pck versions
installed.packages()["RNetLogo","Version"] 
installed.packages()["rJava","Version"]

# install relevant packages   
packages <- c("Matrix","deSolve","mvtnorm","LaplacesDemon","coda","adaptMCMC","sp","RNetLogo","ggplot2","RCurl","RColorBrewer","Interpol.T","lubridate","ggExtra","tidyr","ggthemes","reshape2","pse","sensitivity","beepr")  
if(require(packages)){
  install.packages(packages,dependencies = T)
}
# load annoying packages manually because they're stubborn 
if(mac==0){
  install.packages("RNetLogo")
  install.packages("RCurl")
  install.packages("Interpol.T")
  install.packages("lubridate")
  install.packages("tidyr")
  install.packages("ggthemes")
  install.packages("ggExtra")
  install.packages("beepr")
}
ppp <- lapply(packages,require,character.only=T)
if(any(ppp==F)){cbind(packages,ppp);cat("\n\n\n ---> Check packages are loaded properly <--- \n\n\n")}

cs <- list() # diagnostics list for checking NAs in create snails command  

# load plot function 
script <- getURL("https://raw.githubusercontent.com/darwinanddavis/plot_it/master/plot_it.R", ssl.verifypeer = FALSE)
eval(parse(text = script))
display.brewer.all()
# Set global plotting parameters
cat("plot_it( \n0 for presentation, 1 for manuscript, \nset colour for background, \nset colour palette 1. use 'display.brewer.all()', \nset colour palette 2. use 'display.brewer.all()', \nset alpha for colour transperancy, \nset font style \n)")
plot_it(1,"blue","YlOrRd","Greens",1,"mono") # set plot function params       
plot_it_gg("white") # same as above for ggplot     

################################  compile packages and load files ###################################

### Install rtools and gcc for using C code and coda package 
#### https://cran.r-project.org/bin/macosx/tools/

# define paths for gcc compiler 
if(mac==1){ #### Mac OSX
  rtools <- "/usr/local/clang6/bin"
  gcc <- paste0("usr/local/clang6/gcc-",ver_gcc,"/bin")
  # Mac OSX
}else{ #### Windows 
  rtools <- "C:\\Rtools\\bin"
  gcc <- paste0("C:\\Rtools\\gcc-",ver_gcc,"\\bin")
}

#### point to path on comp to access rtools and gcc for C compiler  
path <- strsplit(Sys.getenv("PATH"), ";")[[1]]
new_path <- c(rtools, gcc, path)
new_path <- new_path[!duplicated(tolower(new_path))]
Sys.setenv(PATH = paste(new_path, collapse = ";"))

if(mac==1){
  # dyn.unload("IndividualModel_IBM3.so") # unLoad .so (Mac OSX
  system(paste0("R CMD SHLIB ",deb_compile,".c")) # generates .o and .so files 
  dyn.load(paste0(deb_compile,".so")) # Load .so (Mac OSX)
}else{
  # compile model from C definition
  # dyn.unload(paste0(deb_compile,".dll")) # unload dll (Windows only)
  system(paste0("R CMD SHLIB ",deb_compile,".c"))
  dyn.load(paste0(deb_compile,".dll"))# Load dll (Windows only)
}

####################################  load deb params #######################################

# load DEB starvation model parameters and create mcmc (and convert mcmc chain to coda format)
samps = readRDS(deb_samps)
samps <- as.mcmc(samps[, c("iM", "k", "M", "EM", "Fh", "muD", "DR", "fe", "yRP",
                           "ph", "yPE", "iPM", "eh", "mP", "alpha", "yEF", "LM",
                           "kd", "z", "kk", "hb", "theta", "mR", "yVE", "yEF2",
                           "sd.LI1", "sd.LU1", "sd.EI1", "sd.EU1", "sd.W1",  "sd.LI2",
                           "sd.LU2", "sd.EI2", "sd.EU2", "sd.W2", "gammaH", "gammaP", "lpost")])


# ---------- summarise and plot estimated params
svar <- "M" # select variable 
sampsvar <- samps[,svar] # pull from mcmc
summary(sampsvar) # get mean, sd, se, and quantiles for each input variable  

den <- density(sampsvar) # get AUC
densplot(sampsvar, show.obs = F,type="n") # density estimate of each variable
polygon(den, col=adjustcolor(colv,alpha=0.5),border=colv) # fill AUC 

plot(sampsvar,trace=T,density=T,col=colv) # traceplot (below) and density plot (above)
# intensive  
traceplot(sampsvar,smooth=T,type="l",lwd=0.3,xlim=c(0,length(sampsvar)),col=colv[2],xlab=paste0("Iterations"),ylab=paste0("Sampled values"),main=paste0("Sampled values over iterations for ",svar)) # iterations vs sampled valued per variable 

if(mcmcplot==1){
  par(mfrow=c(1,1))
  plotlist <- list()
  pdf("mcmc_vars.pdf",onefile = T,paper="a4")
  for(i in colnames(samps)){
    par(bty="n", las = 1)
    if(traceplot==1){
      traceplot(sampsvar,smooth=T,type="l",xlim=c(0,length(sampsvar)),col=colv[2],xlab=paste0("Iterations"),ylab=paste0("Sampled values"),main=paste0("Sampled values over iterations for ",svar)) # iterations vs sampled valued per variable
    }
    svar <- i # select variable 
    sampsvar <- samps[,svar] # pull from mcmc
    den <- density(sampsvar) # get AUC
    densplot(sampsvar, show.obs = F,type="n",main=paste0("Density estimate of ",i)) # density estimate of each variable
    polygon(den, col=adjustcolor(colv,alpha=0.5),border=colv) # fill AUC 
  }
  dev.off()
} # end mcmcplot
# ----------

# get the best fit DEB parameters to match the data (using mcmc)
pars = as.vector(data.frame(samps)[max(which(data.frame(samps)$lpost >= max(data.frame(samps)$lpost) -0.001)),])
pars["Fh"] = 2 # f_scaled (for v.1.1)
pars["ENV"] = 500 # Units: L
pars["r"] = 1   # Units: day-1
pars["step"] = 1  # Units: day
pars["epsilon"] = 20 # Units: L host-1, day-1 (Rounded estimate from Civitello and Rohr)
pars["sigma"] = 0.5 
pars["m_M"] = 1   # Units: day-1
pars["m_Z"] = 1   # Units: day-1
pars["M_in"] = 10
pars["K"] = 5
pars["Det"] = 0.1 # Units mg C/L-1 d-1 (detritus) 
pars["yRP"] = 0.824 # new  0.0471 * (1/0.4) * 7 [expected shedding output per week * (snails shed 40% of their total cercariae during 9-11 AM.) * 7 days] 

####################################  solve deb eqs #######################################
# display list of param definitions
DEB = function(step, Food, L, e, D, RH, P, RP, DAM, HAZ, iM, k, M, EM, 
               Fh, muD, DR, yRP, ph, yPE, iPM, eh, mP, alpha, yEF,
               LM, kd, z, kk, hb, theta, mR, yVE, ENV, Lp, SAtotal, r, K, Det){
  # starting conditions 
  initials = c(Food=Food, L=L, e=e, D=D, RH=RH, P=P, RP=RP, DAM=DAM, HAZ=HAZ)
  # deb parameters
  parameters = c(iM, k, M, EM, Fh, muD, DR, yRP, ph, yPE, iPM,
                 eh, mP, alpha, yEF, LM, kd, z, kk, hb, theta, mR, yVE, ENV, Lp, SAtotal, r, K, Det)
  # estimate starting deb conditions using fitted params by solving ode's
  ## return survival and host shell length  
  DEBstep <- lsoda(initials, c(0,step), func = "derivs", dllname = deb_compile, 
                   initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"), maxsteps=500000,
                   as.numeric(parameters),  rtol=1e-6, atol=1e-6, hmax=1)
  DEBstep[2, 2:12] # 12 = survival
} # end deb model

### deb output for each timestep 
result = DEB(step=1, Food=5, L=10, e=0.9, D=as.numeric(pars["DR"]), RH=0, P=0, RP=0, DAM=0, HAZ=0, iM=pars["iM"], k=pars["k"], M=pars["M"], EM=pars["EM"], 
             Fh=pars["Fh"], muD=pars["muD"], DR=pars["DR"], yRP=pars["yRP"], ph=pars["ph"], yPE=pars["yPE"], iPM=pars["iPM"], eh=pars["eh"],
             mP=pars["mP"], alpha=pars["alpha"], yEF=pars["yEF"], LM=pars["LM"], kd=pars["kd"], z=pars["z"], kk=pars["kk"], hb=pars["hb"],
             theta=pars["theta"], mR=pars["mR"], yVE=pars["yVE"], ENV=pars["ENV"], Lp=10,SAtotal=7007.822, r=pars["r"], K=pars["K"], Det=pars["Det"]) 

### Exposure submodel
# pass the deb state vars into infection model 
Infection = function(snail.stats, miracidia, parameters){
  # Parameters
  epsilon = as.numeric(parameters["epsilon"])
  sigma = as.numeric(parameters["sigma"])
  ENV = as.numeric(parameters["ENV"])
  m_M = as.numeric(parameters["m_M"])
  step = as.numeric(parameters["step"])
  
  # Later calculations depend on exposure probabilities
  exp.rates = epsilon/ENV*(snail.stats[,"L"]>0) # This is just to get uniform exposure rates
  sum.exp.rates = sum(exp.rates)
  
  # Probabilities for fate of miracidia
  ## Still in water
  P.left.in.water = exp(-(m_M+sum(exp.rates))*step)
  ## Infect a snail
  P.infects.this.snail = (1 - P.left.in.water)*(sigma*exp.rates/(m_M+sum.exp.rates)) 
  ## Die in water or fail to infect
  P.dead = (1 - P.left.in.water)*(m_M/(m_M+sum.exp.rates)) + sum((1 - P.left.in.water)*((1-sigma)*exp.rates/(m_M+sum.exp.rates)))
  
  prob.vector = c(P.infects.this.snail, P.left.in.water, P.dead)
  
  # Multinomial outcome from number of miracidia in env based on their survival probability
  rmultinom(n=1, size=miracidia, prob=prob.vector)
  #sum(P.left.in.water, P.invades.this.snail, P.dead)
} # end infection model 

### update all the snails @netlogo
update.snails = function(who, new.L, new.e, new.D, new.RH, new.P, new.RP, new.DAM, new.HAZ, new.LG){
  paste("ask snail", who, 
        "[set L", new.L,
        "set ee", new.e,
        "set D", new.D,
        "set RH", new.RH,
        "set P", new.P,        
        "set RPP", new.RP,       
        "set DAM", new.DAM,
        "set HAZ", new.HAZ,
        "set LG", new.LG, # new max length
        "]")
} # end host update

#Example update
#paste(mapply(update.snails, who=snail.stats[,"who"], new.L=L, new.e=e, new.D=D, new.RH=RH, new.P=P, new.RP=RP, new.DAM=DAM, new.HAZ=HAZ), collapse=" ")

geterrmessage() # check if there were any error messages  

###########################################################################################
####################################  load netlogo ######################################## 
###########################################################################################
# @netlogo

# working NLStart in RStudio. works with gui=F (2018/09/24)
if(gui==0){
  NLStart(nl.path,gui=F,nl.jarname = paste0("netlogo-",ver_nl,".jar")) # open netlogo without a gui  
}else{
  NLStart(nl.path,nl.jarname = paste0("netlogo-",ver_nl,".jar")) # open netlogo
}

NLLoadModel(paste0(model.path,nl.model),nl.obj=NULL) # load model  
# if java.lang error persists on Mac, try copying all .jar files from the 'Java' folder where Netlogo is installed into the main Netlogo folder   	

################################################################################################
####################################  start netlogo sim ######################################## 
################################################################################################
resource_type="detritus" # set resource type
resources="event" # set resource cycles

# algae params
rg_pars <- c(0.1,0.25,0.5) # resource growth rates (rs)
alpha_pars <- c(0,0.25,0.5,0.75,1) # amplitude of resources (alphas)
rho_pars <- c(1,seq(10,120,10)) # periodicity of resources (rhos)

# detritus params
detr_pars <- c(0.1,0.25,0.5) # detritus input (mg L^-1 day^-1)

# mortality params 
hb_pars <- c(0.001, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1); hb_pars

# set type of resource input @netlogo
set_resource_type<-function(resource_type){ # set resource input in env  
  if(resource_type == "detritus"){NLCommand("set resource_type \"detritus\" ")}else{NLCommand("set resource_type \"algae\" ")}}
set_resource_type(resource_type) # set resource type: "detritus" or "algae"  @netlogo

# set type of resource dynamics @netlogo
set_resources<-function(resources){ # set resource input in env  
  if (resources == "cyclical"){NLCommand("set resources \"cyclical\" ")}else{NLCommand("set resources \"event\" ")}}
set_resources(resources) # set resources: "cyclical" or "event"  @netlogo
cat("\nResource type = ",resource_type,"\nResources = ",resources)

if(save_to_file==1){pdf(paste0(wd,"/master_sim.pdf"),onefile=T,paper="a4")}

me_pars <- 0 # set molluscicide events to 0 
if(resource_type=="detritus"){
  detr_pars <- detr_pars; alpha_pars <- 0; rho_pars <- 1; rg_pars <- 0;cat("detritus input = ",detr_pars,"\nrg = ",rg_pars,"\nhb = ", hb_pars)
}else{detr_pars <- 0;cat("detritus input = ", detr_pars,"\nrg = ",rg_pars)}
# set resource to cycle or be constant
if(resource_type=="algae"){
  if(resources=="cyclical"){alpha_pars <- alpha_pars; rho_pars <- rho_pars ; rg_pars <- rg_pars;cat("alphas = ",alpha_pars,"\nrhos = ",rho_pars,"\nrgs = ",rg_pars)
  }else{alpha_pars <- 0; rho_pars <- 1; rg_pars <- rg_pars;cat("\nalphas = ",alpha_pars,"\nrhos = ",rho_pars,"\nrg = ",rg_pars,"\nhb = ", hb_pars)}
}

snail_control = 0 # run molluscicide sims? 
Env_G = numeric() # create empty environment vector 
day <- 1 # reset sim days 

# individual outputs
cerc_list <- list() # cercariae   
food_list <- list() # food in env 
juv_list <- list() # juvenile hosts
adult_list <- list() # adult hosts 
infec_list <- list() # infected hosts
infec_shed_list <- list() # infected shedding hosts
hl_list <- list() # host length
pmass_list <- list() # parasite biomass 
host_biomass_list <- list() # host biomass list 

# master outputs
cerc_master <- list() # master list for cerc density (Env_Z) 
food_master <- list() # master list for food dynamics (Env_F) 
juv_master <- list() # master list for total host pop () 
adult_master <- list() # master list for total host pop () 
infec_master <- list() # master list for infected host pop () 
infec_shed_master <- list() # master list for infected shedding host pop
hl_master <- list() # master list for host length
pmass_master <- list() # master list for parasite biomass 
host_biomass_master <- list() # master list for host biomass

# define plot window
plot.matrix <- matrix(c(length(alpha_pars),length(rho_pars)))
par(mfrow=plot.matrix)

if(snail_control==1){me_pars <- me_pars; hb_pars <- 0.001}else{me_pars <- 1000000; hb_pars <- hb_pars}; cat("Snail control will occur every ",max(me_pars)/length(me_pars)," days \n Mortality is ",me_im) 
if(resource_type=="algae"){detr_pars <- 0; algae <- rg_pars}else{detr_pars <- detr_pars; rg_pars <- 0}
cat("algae:",rg_pars,"\ndetritus:",detr_pars,"\nrho:",0,"\nalpha:",alpha_pars,"\nmortality (if not mollusciciding):",hb_pars,"\nmolluscicide days:",me_pars, "\nmolluscicide impact: ",me_im)

# file name to save results
global_output_fh = paste0(wd,"/global_output_",resource_type,"_",resources,"_me_test.R"); global_output_fh

####################################  start netlogo sim ######################################## 
for(hb in hb_pars){
  for(detr in detr_pars){ # loop through detritus inputs 
    for(alpha in alpha_pars){ # loop through alphas (amplitude in food cycle)
      for(rho in rho_pars){ # loop through rhos (periodicity of food cycle)
        for(rg in rg_pars){ # loop through rgs (food growth rates)
          for(me in me_pars){ # loop through mes (molluscicide events)
            for(me_im in me_im_pars){ # loop through mes (molluscicide events)
              NLCommand("setup")
              day <- 1 # reset days 
              for(t in 1:n.ticks){ # start nl sim  @netlogo
                
                # define food dynamics for cyclical algal (logistic food growth equation) or detritus food sources
                alpha <- alpha # amplitude of resources
                rho <- rho  # periodicity (time range of resource cycles)  
                rg <- rg # resource growth rate 
                rg_t <- rg + alpha * rg * sin(2 * pi * t/rho) # equilibrium cyclical resource dynamics (19-12-18)
                pars["r"] <- rg_t # set resource growth rate 
                pars["Det"] <- detr # Units mg C/L-1 d-1 (detritus)
                # snail control (me events)
                if(snail_control==1){
                  if(day==me){
                    hb <- me_im
                    cat("---------- DEAD SNAILS \ndays = ", day,"\nnumber of snails = ", length(snail.stats$L))
                  }else{
                    hb <- hb}
                }else{hb <- hb}
                
                # set environment variables 
                environment = as.numeric(NLGetAgentSet(c("F", "M", "Z", "G"), "patches")) # calc food, free miracidia, cercariae released, and eggs, per patch
                
                # ----------------------------------- if there are hosts 
                if(NLReport("count snails") > 0){
                  # set host variables
                  snail.stats = NLGetAgentSet(c("who", "L", "ee", "D", "RH", "P", "RPP", "DAM", "HAZ", "LG"), "snails")
                  N.snails = length(snail.stats[,"L"])
                  
                  # Infect snails
                  Infection.step = as.vector(Infection(snail.stats, environment[2], pars)) # Who gets infected
                  snail.stats[which(Infection.step[1:N.snails] > 0),"P"] = snail.stats[which(Infection.step[1:N.snails] > 0),"P"] + 2.85e-5 # add biomass of one miracidia
                  
                  # Update DEBS, HAZ=0 so survival probs are calculated for the current day
                  snail.update = t(mapply(DEB, L=snail.stats[,2], e=snail.stats[,3], D=snail.stats[,4], RH=snail.stats[,5],
                                          P=snail.stats[,6], RP=snail.stats[,7], DAM=snail.stats[,8], Lp=snail.stats[,10],# Food=environment[1]*(snail.stats[,2]^2)/sum(snail.stats[,2]^2), # update food availability per snail 
                                          MoreArgs = list(step=1, HAZ=0, Food=environment[1],# constant food available (23-1-19)
                                                          iM=pars["iM"], k=pars["k"], M=pars["M"], EM=pars["EM"], Fh=pars["Fh"], 
                                                          muD=pars["muD"],
                                                          DR=pars["DR"], yRP=pars["yRP"], ph=pars["ph"], yPE=pars["yPE"], iPM=pars["iPM"], eh=pars["eh"],
                                                          mP=pars["mP"], alpha=pars["alpha"], yEF=pars["yEF"], LM=pars["LM"], kd=pars["kd"], z=pars["z"], 
                                                          kk=pars["kk"], 
                                                          hb=hb,
                                                          theta=pars["theta"], mR=pars["mR"], yVE=pars["yVE"], SAtotal= sum(snail.stats[,2]^2), 
                                                          ENV=pars["ENV"], r=pars["r"], K=pars["K"], 
                                                          Det=pars["Det"]))) # detritus (Det) defined in C file
                  
                  L = snail.update[,"L"] # host structural length
                  e = snail.update[,"e"] # host scaled reserve density    
                  D = snail.update[,"D"] # host development 
                  RH = snail.update[,"RH"] # host energy to reproduction buffer  
                  DAM = snail.update[,"DAM"] # host damage from starvation  
                  HAZ = snail.update[,"HAZ"] # host hazard rate from starvation   
                  LG = snail.update[,"LG"] # host shell length  
                  P = snail.update[,"P"] # parasite mass (sum within host)
                  RP = snail.update[,"RP"] # parasite reproductive buffer  
                  
                  # ingestion = environment[1] - sum(snail.update[,"Food"]) # food intake by host from environment (for v.1.1)
                  chi <- pars["M"]/(1 + pars["EM"]) # length to volume conversion factor for getting biomass
                  host_biomass <- sum(chi * L^3) # get total host biomass
                  
                  hl_list[t] <- mean(L) # get average host lengths per model step 
                  pmass_list[t] <- sum(P) # get total parasite mass per model step
                  host_biomass_list <- host_biomass # get host mass per model step
                  
                  Eggs = floor(RH/0.015)  # Figure out how many (whole) eggs are released  
                  # if(day==me){Eggs <- Eggs[1:round(0.1*length(Eggs))]} # kill off 90% of snail eggs in water with molluscicide event  
                  RH = RH %% 0.015        # Remove released cercariae from the buffer
                  Cercs = floor(RP/4e-5)  # Figure out how many (whole) cercs are released
                  RP = RP %% 4e-5         # Remove released cercariae from buffer
                  Eggs = as.integer(Eggs); Cercs = as.integer(Cercs) 
                  
                  # Update environment 
                  Env_M = as.numeric(Infection.step[N.snails + 1] + pars["M_in"]) # total miracidia density 
                  Env_Z = as.numeric(environment[3]*exp(-pars["m_Z"]*pars["step"]) + sum(Cercs)/pars["ENV"]) # total cerc density
                  Env_G = as.integer(Env_G) # set pop density outputs to integer to pass into Env_G and rbinom func
                  ifelse(day==me,Env_G[day] <- max(0, 0.1*sum(Eggs),na.rm=T),Env_G[day] <- max(0, sum(Eggs),na.rm=T)) # kill off 90% of snail eggs in water with molluscicide event 
                  Env_G[day] <- max(0, sum(Eggs),na.rm=T)
                  
                  Env_G[is.na(Env_G)] <- 0 # turn NAs to 0 to feed into rbinom function  
                  # Env_F = max(0.001, as.numeric(pars["K"]*environment[1]/(environment[1] + (pars["K"] - environment[1])*exp(-pars["r"]*pars["step"])) - ingestion)) # Analytical soln to logistic - ingestion (alphas [1,100]) (original r growth equation)
                  # Env_F = max(0.001, as.numeric(pars["K"]*environment[1]/(environment[1] + (pars["K"] - environment[1])*exp(-rg_t*pars["step"])) - ingestion)) # Analytical soln to logistic - ingestion with equilibrium resource growth wave (rg_t) (alphas [0,1]) (for v.1.1)     
                  # F = F * exp(- r + alpha * r * sin(2 * pi * t/rho) * s) * (1 - F/K) - f(i_{M} * sum(L^2) # v. 1.2 algae and detritus with cyclical algal growth
                  # r_t <- pars["r"] + alpha * pars["r"] * sin(2 * pi * t/rho) # equilibrium resource dynamics (static)
                  Env_F = max(0.001, snail.update[1,"Food"]) # algal or detritus food sources (for v.1.2)
                  # Command back to NL @netlogo
                  snail.commands = paste(mapply(update.snails, who=snail.stats[,"who"], new.L=L, new.e=e, new.D=D, new.RH=RH, new.P=P, new.RP=RP, new.DAM=DAM, new.HAZ=HAZ, new.LG=LG), collapse=" ")
                  NLCommand(snail.commands) 
                  
                }else{ # ----------------------------------- if there are no hosts 
                  
                  Env_M = as.numeric(environment[2]*exp(-pars["m_M"]*pars["step"]) + pars["M_in"]) # total miracidia density 
                  Env_Z = as.numeric(environment[3]*exp(-pars["m_Z"]*pars["step"])) # total cerc density
                  Env_G = as.integer(Env_G) # set pop density outputs to integer to pass into Env_G and rbinom func
                  Env_G[day] <- 0
                  Env_F = ifelse(pars["Det"] == 0, as.numeric(pars["K"]*environment[1]/(environment[1] + (pars["K"] - environment[1])*exp(-pars["r"]*pars["step"]))), as.numeric(environment[1] + pars["Det"]))
                } # end no hosts
                
                # update patch variables 
                NLCommand("ask patch 0 0 [set F", Env_F, "set M", Env_M, "set Z", Env_Z, "set G", Env_G[day], "]")
                
                # kill snail eggs with molluscicide event  
                if(day > 10){
                  if(snail_control==1){ 
                    if(day==me){
                      # Env_G[day - 10, day-1] <- round(0.1*Env_G[day - 10, day-1]) # mollusciciding kills eggs
                      create_snails <- rbinom(n=1, size=Env_G[day - 10], prob=0.1)
                    }else{create_snails <- rbinom(n=1, size=Env_G[day - 10], prob=0.5)}
                  }else{create_snails <- rbinom(n=1, size=Env_G[day - 10], prob=0.5)}
                  NLCommand("create-snails ", create_snails, "[set L 0.75 set ee 0.9 set D 0 set RH 0 set P 0 set RPP 0 set DAM 0 set HAZ 0 set LG 0.75]")
                } # end create snails
                
                NLCommand("go") # run @netlogo sim steps
                #cs[t] <- rbinom(n=1, size=Env_G[day - 10], prob=0.5) # list to check 'create snails' output doesn't produce NAs
                day = day + 1 
                # results outputs
                cerc_list[t] <- Env_Z # get cercariae density 
                food_list[t] <- Env_F # get food growth
                juv_list[t] <- length(which(snail.stats$RH==0)) # get juvenile hosts
                adult_list[t] <- length(which(snail.stats$RH>0)) # get adult hosts
                infec_list[t] <- length(which(snail.stats$P>0)) # get just infected hosts
                infec_shed_list[t] <- length(which(snail.stats$RP>0)) # get infected hosts that are shedding
              } # --------------------------------------- end nl sim
              # save individual outputs 
              cerc_list <- as.numeric(cerc_list) 
              food_list <- as.numeric(food_list)
              juv_list <- as.numeric(juv_list)
              adult_list <- as.numeric(adult_list)
              infec_list <- as.numeric(infec_list)
              infec_shed_list <- as.numeric(infec_shed_list)
              hl_list <- as.numeric(hl_list)
              pmass_list <- as.numeric(pmass_list)
              host_biomass_list <- as.numeric(host_biomass_list)
              
              # save master outputs 
              cerc_master[[length(cerc_master)+1]] <- cerc_list # cerc master list
              food_master[[length(food_master)+1]] <- food_list # food master list
              juv_master[[length(juv_master)+1]] <- juv_list # juv pop master list
              adult_master[[length(adult_master)+1]] <- adult_list # adult pop master list
              infec_master[[length(infec_master)+1]] <- infec_list # infected host pop master list
              infec_shed_master[[length(infec_shed_master)+1]] <- infec_shed_list # infected shedding host pop master list
              hl_master[[length(hl_master)+1]] <- hl_list # host length master
              pmass_master[[length(pmass_master)+1]] <- pmass_list # parasite mass master
              host_biomass_master[[length(host_biomass_master)+1]] <- host_biomass_list # host biomass master
              # day_master[[length(day_master)+1]] <- day_list
              ### plot outputs 
              #   plot(cerc_list,type="l",las=1,bty="n",ylim=c(0,do.call(max,cerc_master)),col=round(do.call(max,cerc_master)),
              # 	main=paste0("alpha = ",alpha, "; rho = ", rho, "; r = ", rg),ylab="Cercariae density",xlab="Days") 
              # paste(expression("alpha = ",alpha, "; rho = ", rho, "; r = ", rg))
              #   text(which(cerc_list==max(cerc_list)),max(cerc_list),paste0("a= ",alpha," \n p= ",rho)#,col=max(cerc_list),
              #        )
              #abline(h=which(cerc_list==max(cerc_list)),type=3,col=round(do.call(max,cerc_master))) # draw line at max peak
              if(save_to_file==1){dev.off()}
            } # ----end me_im
          } # --------------- end mes
        } # ------------------------------ end rgs
      } # --------------------------------------------- end rhos
    } # ----------------------------------------------------------- end alphas
  } # ------------------------------------------------------------------------- end detritus
} # end hb pars 
####################################  end netlogo sim ######################################## 

# results output 
# save sim results to dir 
str(list(cerc_master,food_master,juv_master, adult_master,infec_master,infec_shed_master,hl_master,pmass_master,host_biomass_master))  
global_output <- list(cerc_master,food_master,juv_master, adult_master,infec_master,infec_shed_master,hl_master,pmass_master,host_biomass_master) 
# global_output <- cerc_master # save just cercs
saveRDS(global_output,global_output_fh) # save to dir 
cat("Output saved in ", global_output_fh)
# read in saved sim results
cat("order = cerc, food, juv, adult, infected, infected shedding, mean host length, mean parasite mass, summed host biomass")

# ------------------------- plot individual outputs -------------------------
mm_ = readRDS(global_output_fh)

#------------------------- plot all sim results in one window -------------------------
sim_type <- resource_type
require(gridExtra)
gspl <- list()
ttl_list <- c("cerc","food", "juv", "adult", "infec", "infec (shed)", "host L", "parasite mass","host biomass")
if(sim_type=="algae"){legend_pars <- rg_pars}; if(sim_type=="detritus"){legend_pars <- detr_pars};  

# choose sim to plot
global_sim_plot <- mm_
for(g in 1:length(global_sim_plot)){
  par(bty="n", las = 1)
  mm <- global_sim_plot[[g]]
  y_m <- melt(mm);y_m
  gspl[[g]] <- ggplot() +
    geom_line(data = y_m, aes(x = rep.int(1:length(mm_[[1]][[1]]),max(L1)) , y = value, group = L1, colour=factor(L1)), ) +
    theme_tufte() +
    theme(legend.position = "right") +
    scale_color_manual(name=paste0(sim_type," days"), # stupid ggplot legend @ggplotlegend  # http://r-statistics.co/Complete-Ggplot2-Tutorial-Part2-Customizing-Theme-With-R-Code.html
                       labels = legend_pars,
                       values = unique(y_m$L1)) +
    geom_vline(xintercept = me_pars, linetype="dotted", # add me event lines  
               color = "gray", size=0.5) +
    labs(title=ttl_list[g],x="",y="") +
    if(g==length(global_sim_plot)-1){
      labs(x="Time") 
    }else{
      theme(legend.position="none")
    } 
}
do.call(grid.arrange,gspl) # plot in one window    
gspl[1] # plot individual outputs
global_output_fh

#################################################################################################
##########################################  end body ############################################ 
#################################################################################################

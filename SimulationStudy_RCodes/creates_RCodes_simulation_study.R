# setwd('')
library(tidyverse)
n <- 1400 ## number of repeats of simulation for all file
raw.dir = "/hpc/home/ma521/SimulationStudy_RCodes/"
code.dir = "/cwork/ma521/HPImage/"

for (i in 1:n){
  fnm<-ifelse(i%/%100==0,"sim_himage13.R",
              ifelse(i%/%100==1,"sim_himage14.R",
                     ifelse(i%/%100==2,"sim_himage15.R",
                            ifelse(i%/%100==3,"sim_himage20.R",
                                   ifelse(i%/%100==4,"sim_himage25.R",
                                          ifelse(i%/%100==5,"sim_himage35.R",
                                                 ifelse(i%/%100==6,"sim_himage60.R",
                                                        ifelse(i%/%100==7,"sim_hpath13.R",
                                                               ifelse(i%/%100==8,"sim_hpath14.R",
                                                                      ifelse(i%/%100==9,"sim_hpath15.R",
                                                                             ifelse(i%/%100==10,"sim_hpath20.R",
                                                                                    ifelse(i%/%100==11,"sim_hpath25.R",
                                                                                           ifelse(i%/%100==12,"sim_hpath35.R",
                                                                                                  "sim_hpath60.R")))))))))))))
  if(i==100){
    fnm<-"sim_himage13.R"
  } 
  
  if(i==200){
    fnm<-"sim_himage14.R"
  } 
  
  if(i==300){
    fnm<-"sim_himage15.R"
  } 

  if(i==400){
    fnm<-"sim_himage20.R"
  } 
  
  if(i==500){
    fnm<-"sim_himage25.R"
  } 
  
  if(i==600){
    fnm<-"sim_himage35.R"
  } 
  
  if(i==700){
    fnm<-"sim_himage60.R"
  } 
  
  if(i==800){
    fnm<-"sim_hpath13.R"
  } 
  
  if(i==900){
    fnm<-"sim_hpath14.R"
  } 
  
  if(i==1000){
    fnm<-"sim_hpath15.R"
  } 
  
  if(i==1100){
    fnm<-"sim_hpath20.R"
  } 
  
  if(i==1200){
    fnm<-"sim_hpath25.R"
  } 
  
  if(i==1300){
    fnm<-"sim_hpath35.R"
  } 
  
  if(i==1400){
    fnm<-"sim_hpath60.R"
  } 
  
  source_code <- paste0("source('", raw.dir, fnm,"')")
  source_code <- as.data.frame(source_code) %>% mutate(V1=source_code) %>% dplyr::select(V1)
  seed_code <- paste0("iter <- ", ifelse(i%%100==0,100,i%%100))
  seed_code <- as.data.frame(seed_code) %>% mutate(V1=seed_code) %>% dplyr::select(V1)
  code_file <- paste0(code.dir, i, '.R')
  R_code <- rbind(seed_code, source_code)
  write_delim(R_code, code_file, col_names = F, delim = '')
}

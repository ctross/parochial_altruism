############### Plots of ind-level stuff

N = ncol(model_dat_Coast$Taking)
give = model_dat_Coast$Taking 
eth = model_dat_Coast$Indigenous
friends = model_dat_Coast$Friends

give[74,74]=21 #fix data error

simp = array(NA,c(N,5))

for(i in 1:N){
 scrap = rep(0,5)
    for(j in 1:N){
 scrap[1] = scrap[1] + ifelse(j==i, give[i,j], 0 ) 

 scrap[2] = scrap[2] + ifelse((eth[j]==eth[i]) & (i != j), give[i,j], 0 ) 
 scrap[3] = scrap[3] + ifelse((eth[j]!=eth[i]) & (i != j), give[i,j], 0 ) 

 scrap[4] = scrap[4] + ifelse((eth[j]==eth[i]) & (i != j), friends[i,j] + friends[j,i], 0 ) 
 scrap[5] = scrap[5] + ifelse((eth[j]!=eth[i]) & (i != j), friends[i,j] + friends[j,i], 0 ) 
               }
 scrap = scrap/sum(scrap)
if(eth[i]==0){
simp[i,] = scrap / c(1, (1-0.2413793), 0.2413793 ,1,1)
} else{
simp[i,] = scrap / c(1, 0.2413793,(1-0.2413793) ,1,1)
}
}


df2c = data.frame(Ethnicity=c("Afrocolombian", "Embera")[eth+1], Self=(simp[,1]), InGroup=(simp[,2]), OutGroup=(simp[,3]),
                  InGroupFriends=(simp[,4]), OutGroupFriends=(simp[,5]),NoFood=model_dat_Coast$NoFood, Site="Coast")

############### Plots of ind-level stuff

N = ncol(model_dat_Inland$Taking)
give = model_dat_Inland$Taking 
eth = model_dat_Inland$Indigenous
friends = model_dat_Inland$Friends


simp = array(NA,c(N,5))

for(i in 1:N){
 scrap = rep(0,5)
    for(j in 1:N){
 scrap[1] = scrap[1] + ifelse(j==i, give[i,j], 0 ) 

 scrap[2] = scrap[2] + ifelse((eth[j]==eth[i]) & (i != j), give[i,j], 0 ) 
 scrap[3] = scrap[3] + ifelse((eth[j]!=eth[i]) & (i != j), give[i,j], 0 ) 

 scrap[4] = scrap[4] + ifelse((eth[j]!=eth[i]) & (i != j), friends[i,j] + friends[j,i], 0 ) 
 scrap[5] = scrap[5] + ifelse((eth[j]!=eth[i]) & (i != j), friends[i,j] + friends[j,i], 0 ) 
               }
 scrap = scrap/sum(scrap)
if(eth[i]==0){
simp[i,] = scrap / c(1, (1-0.1390728), 0.1390728 ,1,1)
} else{
simp[i,] = scrap / c(1, 0.1390728,(1-0.1390728) ,1,1)
}
}



df2i = data.frame(Ethnicity=c("Afrocolombian", "Embera")[eth+1], Self=(simp[,1]), InGroup=(simp[,2]), OutGroup=(simp[,3]),
                  InGroupFriends=(simp[,4]), OutGroupFriends=(simp[,5]),NoFood=model_dat_Inland$NoFood, Site="Inland")

df2 = rbind(df2c,df2i)

df2$InterTie=factor(ifelse(df2$OutGroupFriends>0,1,0))


p <- ggplot(df2, aes(InGroup, OutGroup, color = Ethnicity))+
  geom_jitter(width = 0.025, height = 0.025,size=1.25) + scale_colour_manual("Ethnicity: ", values = c("royalblue4", "orange3"))+
   geom_abline(intercept = 0, slope = 1, linetype=3) + facet_grid(.~Site) +theme(legend.position="bottom", legend.box = "horizontal")+
   labs(y="Leaving rate for out-group", x="Leaving rate for in-group") + theme(strip.text.x = element_text(size=14,face="bold"), 
     strip.text.y = element_text(size=14,face="bold"),axis.text=element_text(size=12),axis.title.y=element_text(size=14,
     face="bold"),axis.title.x=element_text(size=14,
     face="bold"))+theme(strip.text.y = element_text(angle = 360)) + theme(panel.spacing = unit(1, "lines")) +
     theme(legend.title=element_text(size=14), legend.text=element_text(size=14))
p

   
ggsave("Exploitation_Individual_Level.pdf",p,height=4, width=7)

                  

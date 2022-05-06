color <- ifelse(model_dat_Inland$Indigenous == 0, "royalblue4", "orange3")
                                                   
# Plots
gg_SC_Friends <- ggnet2(network(model_dat_Inland$Friends, directed = TRUE),  mode = "fruchtermanreingold", layout.par = list(cell.jitter = 1.175), label = FALSE, size=2, color=color)
gg_SC_Exchange <- ggnet2(network(model_dat_Inland$Transfer, directed = TRUE),  mode = "fruchtermanreingold", layout.par = list(cell.jitter = 1.175), label = FALSE, size=2, color=color)

gg_SC_Give <- ggnet2(network(model_dat_Inland$Giving, directed = TRUE),  mode = "fruchtermanreingold", layout.par = list(cell.jitter = 1.175), label = FALSE, size=2, color=color)
gg_SC_Take <- ggnet2(network(model_dat_Inland$Taking, directed = TRUE),  mode = "fruchtermanreingold", layout.par = list(cell.jitter = 1.175), label = FALSE, size=2, color=color)
gg_SC_Puni <- ggnet2(network(model_dat_Inland$Reducing, directed = TRUE),  mode = "fruchtermanreingold", layout.par = list(cell.jitter = 1.175), label = FALSE, size=2, color=color)



color <- ifelse(model_dat_Coast$Indigenous == 0, "royalblue4", "orange3")
                                                   
# Plots
gg_BS_Friends <- ggnet2(network(model_dat_Coast$Friends, directed = TRUE),  mode = "fruchtermanreingold", layout.par = list(cell.jitter = 1.175), label = FALSE, size=2, color=color)
gg_BS_Exchange <- ggnet2(network(model_dat_Coast$Transfer, directed = TRUE),  mode = "fruchtermanreingold", layout.par = list(cell.jitter = 1.175), label = FALSE, size=2, color=color)

gg_BS_Give <- ggnet2(network(model_dat_Coast$Giving, directed = TRUE),  mode = "fruchtermanreingold", layout.par = list(cell.jitter = 1.175), label = FALSE, size=2, color=color)
gg_BS_Take <- ggnet2(network(model_dat_Coast$Taking, directed = TRUE),  mode = "fruchtermanreingold", layout.par = list(cell.jitter = 1.175), label = FALSE, size=2, color=color)
gg_BS_Puni <- ggnet2(network(model_dat_Coast$Reducing, directed = TRUE),  mode = "fruchtermanreingold", layout.par = list(cell.jitter = 1.175), label = FALSE, size=2, color=color)


ggsave("Friends_coast.pdf",gg_BS_Friends,width=2.5,height=2.5)
ggsave("Transfers_coast.pdf",gg_BS_Exchange,width=2.5,height=2.5)
ggsave("Give_coast.pdf",gg_BS_Give,width=2.5,height=2.5)
ggsave("Take_coast.pdf",gg_BS_Take,width=2.5,height=2.5)
ggsave("Punish_coast.pdf",gg_BS_Puni,width=2.5,height=2.5)

ggsave("Friends_inland.pdf",gg_SC_Friends,width=2.5,height=2.5)
ggsave("Transfers_inland.pdf",gg_SC_Exchange,width=2.5,height=2.5)
ggsave("Give_inland.pdf",gg_SC_Give,width=2.5,height=2.5)
ggsave("Take_inland.pdf",gg_SC_Take,width=2.5,height=2.5)
ggsave("Punish_inland.pdf",gg_SC_Puni,width=2.5,height=2.5)




nullstuff=data.frame(label="Friendships",x=0,y=0)
ftext=ggplot() + geom_text(data = nullstuff, aes(x = x, y = y,label=label),angle=90)+
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())
  ggsave("Friends.pdf",ftext,width=0.5,height=2.5)

  nullstuff=data.frame(label="Sharing",x=0,y=0)
ftext=ggplot() + geom_text(data = nullstuff, aes(x = x, y = y,label=label),angle=90)+
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())
  ggsave("Sharing.pdf",ftext,width=0.5,height=2.5)

  nullstuff=data.frame(label="Allocation",x=0,y=0)
ftext=ggplot() + geom_text(data = nullstuff, aes(x = x, y = y,label=label),angle=90)+
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())
  ggsave("Allocation.pdf",ftext,width=0.5,height=2.5)

  nullstuff=data.frame(label="Exploitation",x=0,y=0)
ftext=ggplot() + geom_text(data = nullstuff, aes(x = x, y = y,label=label),angle=90)+
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())
  ggsave("Exploitation.pdf",ftext,width=0.5,height=2.5)

  nullstuff=data.frame(label="Reduction",x=0,y=0)
ftext=ggplot() + geom_text(data = nullstuff, aes(x = x, y = y,label=label),angle=90)+
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())
  ggsave("Reduction.pdf",ftext,width=0.5,height=2.5)

  nullstuff=data.frame(label="Coastal",x=0,y=0)
ftext=ggplot() + geom_text(data = nullstuff, aes(x = x, y = y,label=label),angle=0)+
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())
  ggsave("Coastal.pdf",ftext,width=2.5,height=0.5)

    nullstuff=data.frame(label="Inland",x=0,y=0)
ftext=ggplot() + geom_text(data = nullstuff, aes(x = x, y = y,label=label),angle=0)+
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())
  ggsave("Inland.pdf",ftext,width=2.5,height=0.5)
  
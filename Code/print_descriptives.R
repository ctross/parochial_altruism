############# Descriptives
# Giving
dgive <- rowSums(model_dat_Inland$Giving)-diag(model_dat_Inland$Giving)
dgive <- dgive[which(rowSums(model_dat_Inland$Giving) != 0)]

print("Give inland")
print(paste("mean=",mean(dgive)))
print(paste("median=",median(dgive)))
print(paste("sd=",sd(dgive)))
print(paste("max=",max(dgive)))
print(paste("min=",min(dgive)))


dgive <- rowSums(model_dat_Coast$Giving)-diag(model_dat_Coast$Giving)
dgive <- dgive[which(rowSums(model_dat_Coast$Giving) != 0)]

print("Give coast")
print(paste("mean=",mean(dgive)))
print(paste("median=",median(dgive)))
print(paste("sd=",sd(dgive)))
print(paste("max=",max(dgive)))
print(paste("min=",min(dgive)))

# Leaving
dgive <- rowSums(model_dat_Inland$Taking)-diag(model_dat_Inland$Taking)
dgive <- dgive[which(rowSums(model_dat_Inland$Taking) != 0)]

print("Leave inland")
print(paste("mean=",mean(dgive)/2))
print(paste("median=",median(dgive)/2))
print(paste("sd=",sd(dgive)/2))
print(paste("max=",max(dgive)/2))
print(paste("min=",min(dgive)/2))


dgive <- rowSums(model_dat_Coast$Taking)-diag(model_dat_Coast$Taking)
dgive <- dgive[which(rowSums(model_dat_Coast$Taking) != 0)]

print("Leave coast")
print(paste("mean=",mean(dgive)/2))
print(paste("median=",median(dgive)/2))
print(paste("sd=",sd(dgive)/2))
print(paste("max=",max(dgive)/2))
print(paste("min=",min(dgive)/2))

# Punishing
dgive <- rowSums(model_dat_Inland$Reducing)-diag(model_dat_Inland$Reducing)
dgive <- dgive[which(rowSums(model_dat_Inland$Reducing) != 0)]

print("Reduce inland")
print(paste("mean=",mean(dgive)))
print(paste("median=",median(dgive)))
print(paste("sd=",sd(dgive)))
print(paste("max=",max(dgive)))
print(paste("min=",min(dgive)))


dgive <- rowSums(model_dat_Coast$Reducing)-diag(model_dat_Coast$Reducing)
dgive <- dgive[which(rowSums(model_dat_Coast$Reducing) != 0)]

print("Reduce coast")
print(paste("mean=",mean(dgive)))
print(paste("median=",median(dgive)))
print(paste("sd=",sd(dgive)))
print(paste("max=",max(dgive)))
print(paste("min=",min(dgive)))

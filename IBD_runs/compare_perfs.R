perfs0 = read.table("perf0.csv", h = T, sep = ",")
perfs0$X = NULL
perfs0
#colMeans(perfs0)

perfs1 = read.table("perf1.csv", h = T, sep = ",")
perfs1$X = NULL
perfs1

print(rbind(apply(perfs0,2,median),apply(perfs1,2,median)))

df_diff = as.data.frame(t(rbind(apply(perfs0,2,median),apply(perfs1,2,median))))
it1 = sum(df_diff$V2-df_diff$V1 >= 0)/nrow(df_diff)*100
print(it1)

print(apply(df_diff, 2, median))

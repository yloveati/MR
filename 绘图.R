res<- read.csv(file = "res.csv", header = T)

mrResult=mr(res)

mrTab=generate_odds_ratios(mrResult)
write.csv(mrTab, file="table.MRresult.csv", row.names=F)
heterTab=mr_heterogeneity(res)
write.csv(heterTab, file="table.heterogeneity.csv", row.names=F)
pleioTab=mr_pleiotropy_test(res)
write.csv(pleioTab, file="table.pleiotropy.csv", row.names=F)


#绘制散点图
mr_scatter_plot(mr_results=mrResult, dat=res)
#森林图
res_single=mr_singlesnp(res)
mr_forest_plot(res_single)
#漏斗图
mr_funnel_plot(singlesnp_results = res_single)
#留一法敏感性分析
mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(res))

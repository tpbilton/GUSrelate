

library(GUSrelate)

ra = readRA("~/Academia/Otago_University/GUSrelate/test/simDS.vcf.ra.tab")
grm = makeGRM(ra, samfile = "~/Academia/Otago_University/GUSrelate/test/test_sim_ploid.csv")

grm$addSampleInfo("~/Academia/Otago_University/GUSrelate/test/test_sim_addSamInfo.csv")

grm$HWEtest(nThreads = 3)
grm$computeGRM(name = "test", ep=0)

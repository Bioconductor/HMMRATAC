# .jcall("HMMR_ATAC.Main_HMMR_Driver", "V", "main", .jarray(list(), "java/lang/String"))
hmmr <- rJava::new(J("HMMR_ATAC.Main_HMMR_Driver"))
hmmr$main(.jarray(list(), "java/lang/String"))

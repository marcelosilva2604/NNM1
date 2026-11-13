args <- commandArgs(trailingOnly=TRUE)
src <- args[1]; csv <- args[2]; out <- args[3]
wrapmessage <- function(...) invisible(NULL)              # só mensagens do pacote
ex <- parse(file.path(src))
for (e in ex) {                                            # avalia só threshold e evalues.RR do pacote
  if (is.call(e) && length(e) >= 3 && as.character(e[[1]]) %in% c("=","<-") &&
      as.character(e[[2]]) %in% c("threshold","evalues.RR")) eval(e, globalenv())
}
d <- read.csv(csv, stringsAsFactors=FALSE)
res <- do.call(rbind, lapply(seq_len(nrow(d)), function(i) {
  r <- suppressMessages(evalues.RR(est=d$irr[i], lo=d$lo95[i], hi=d$hi95[i]))
  data.frame(term=d$term[i], irr=d$irr[i], E_point_R=r["E-values","point"],
             E_ci_R=ifelse(d$irr[i] > 1, r["E-values","lower"], r["E-values","upper"]))
}))
thr <- do.call(rbind, lapply(c(1.10,1.25), function(x){ r <- evalues.RR(est=x); data.frame(term=paste0("residual RR ",x), irr=x, E_point_R=r["E-values","point"], E_ci_R=NA)}))
res <- rbind(res, thr)
write.csv(res, out, row.names=FALSE)
print(res, digits=6)

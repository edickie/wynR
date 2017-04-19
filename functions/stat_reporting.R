report_pvalue <- function(pvalue) {
  if (pvalue > 0.01) {
    pstr <- sprintf('%3.2f',pvalue)
  } else if (pvalue < 0.01 & pvalue > 0.001) {
    pstr <- sprintf('%4.3f', pvalue) 
  } else if (pvalue < 0.001 & pvalue > 0.0001) {
    pstr <- sprintf('%5.4f', pvalue)
  } else if (pvalue < 0.0001) {
    pstr <- sprintf("%2.1e", pvalue)
  }
  return(pstr)
}

## a little function to format t statistics for reporting in the text
report_tstat <- function(tstatistic, df, pvalue) {
  pstr <- report_pvalue(pvalue)
  tstring <- sprintf("t(%i)=%3.2f, p=%s",df,tstatistic, pstr)
  return(tstring)
}
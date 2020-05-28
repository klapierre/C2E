### Subsetting only particular CDR and other studies
###
### Author: Kim Komatsu
###         Kevin Wilcox
###
### Last updated: April 23, 2020


divTrt1 <- div%>%
  
  filter(pulse==0, plot_mani>0, project_name!="e001"&project_name!="e002")


#removing a subset of CDR e001 and e002 treatments to prevent the majority of data being from CDR; keeping lowest, highest, and 10 gm-2 N additions (levels most comparable to other studies)

divCDRe001<-div%>%
  
  filter(site_code=="CDR"&treatment==1|treatment==6|treatment==8|treatment==9,plot_mani>0)

divCDRe002<-div%>%
  
  filter(site_code=="CDR"&treatment=='1_f_u_n'|treatment=='6_f_u_n'|treatment=='8_f_u_n'|treatment=='9_f_u_n',plot_mani>0)
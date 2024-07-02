### P-val Matrix Adjustment Function
##### Used in the Absorption_MS_Analysis Drafts. Included in draft 1, removed in draft 2.
rm_adjusted.corr = function(rmc_mat,df, p.adjust = FALSE, p.adjust.method = "none", threshold = 0.05, cor.method = "repeated"){
  # rmc_mat is the result from a call to function rmcorr_mat from rmcorr library
  threshold = threshold
  
  rmc_mat = rmc_mat
  r.mat = rmc_mat$matrix
  
  if (p.adjust==TRUE){
    p.vals <- p.adjust(rmc_mat$summary$p.vals, 
                       method = p.adjust.method,
                       n = length(rmc_mat$summary$p.vals))
  } else {
    p.vals = rmc_mat$summary$p.vals
  }
  
  cor_df<-df
  
  # initiate empty matrix for the p-values only
  MAT = matrix(NA, nrow = ncol(cor_df), ncol = ncol(cor_df))
  diag(MAT) = 0
  
  rownames(MAT) = names(cor_df)
  colnames(MAT) = names(cor_df)
  
  # starting to convert p.list to p.values matrix
  for(ind in 1:(length(p.vals))){
    var1 = rmc_mat$summary$measure1[ind] # var1
    var2 = rmc_mat$summary$measure2[ind] # var1
    p.value = p.vals[ind] # p value
    MAT[var1, var2] = p.value
    MAT[var2, var1] = p.value
  }
  
  # At this point, MAT has the p.values in a matrix
  
  # subset only the coefficients with p values < 0.05 (or than the specified threshold)
  
  subset = ifelse(MAT < threshold, r.mat, NA)
  rownames(subset) = names(cor_df)
  colnames(subset) = names(cor_df)
  
  output = list(adj.p.values = MAT, threshold.r = subset)
  return(output)
}

### Correlation Chart Function
##### Used in the Absorption_MS_Analysis Drafts. Included in draft 1, removed in draft 2.
rm_correlation.matrix = function(rm_adj_corr, colours_list, Title, Subtitle){
  # rm_adj_corr is the result of a call to rm_adjusted.corr
  pval = rm_adj_corr$adj.p.values  # shows the p.values in a matrix format
  thresR = rm_adj_corr$threshold.r  # shows only corr (r) values where the p value <0.05
  
  # convert to triangle
  thresR[lower.tri(thresR)]<- NA
  pval[lower.tri(pval)]<- NA
  #create melted values
  melted_cormat<-melt(thresR,value.name = "r")
  melted_pmat<-melt(pval, value.name = "p")
  corPlotData<-full_join(melted_cormat, melted_pmat, by = c("Var1", "Var2"))
  corPlotData$r[corPlotData$r==1]<-NA
  corPlotData.c<-na.omit(corPlotData)
  corPlotData.c$stars<-cut(corPlotData.c$p, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))  # Create column of significance
  corPlotData.c$r<-round(corPlotData.c$r,2)
  
  # Create a ggheatmap
  ggheatmap <- ggplot(corPlotData.c, aes(Var2, Var1, fill = r))+
    geom_tile(color = "white")+
    scale_fill_gradientn(
      colours = colours_list,
      limits = c(-1,1)
    )+
    theme_minimal()+ 
    theme(axis.text.x = element_text(vjust = 1, hjust = 1, angle = 45, size = 6))+
    theme(axis.text.y = element_text(size = 6))+
    labs(title = Title,
         subtitle = Subtitle)+
    coord_fixed()
  
  print(ggheatmap)
  
  ggheatmap2<-ggheatmap+
    geom_text(aes(Var2, Var1, label = paste(stars,r, sep ="\n")), color = "black", size = 2) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 10),
      plot.subtitle = element_text(hjust = 0.5, size = 8),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank(),
      legend.justification = c(1, 0),
      legend.position = c(0.4, 0.8),
      legend.direction = "horizontal")+
    guides(fill = guide_colorbar(barwidth = 5, barheight = 0.4,
                                 title.position = "top", title.hjust = 0.5, 
                                 title.theme = element_text(size = 5),
                                 label.theme = element_text(size = 5)))
  
  print(ggheatmap2)
}

neighbours_LR_moved<-function(pt_id){
  #data = spatial_dat
  #seat = data$Y[data$Pt_ID == pt_id]
  #row = data$X[data$Pt_ID == pt_id]
  seat = Y
  row = X
  if (seat>=1 & seat<=12){
    if (-6.5<=row& row<=6.5){ 
      neighbour_L = seat-1
      neighbour_R = seat+1
      L_moved = as.numeric(data%>%filter(X ==row & Y ==neighbour_L)%>%select(moved_Beethoven))
      R_moved = as.numeric(data%>%filter(X ==row & Y ==neighbour_R)%>%select(moved_Beethoven))
      average = mean(c(L_moved,R_moved))
      print(average)
    }
  }
}


adj.cor = function(df, p.adjust = FALSE, p.adjust.method = "none", threshold = 1, cor.method = "kendall"){
  cor_test = rcor.test(df, p.adjust = p.adjust, p.adjust.method = p.adjust.method, method = cor.method)
  r.mat = cor_test$cor.mat     # matrix of coefficient values
  p.list = cor_test$p.values   # p.list will be 3 columns
  
  # initiate empty matrix for the p values only
  MAT = matrix(, nrow = ncol(df), ncol = ncol(df))
  diag(MAT) = 0
  
  # starting to convert p.list to p.values matrix
  for(ind in 1:(length(p.list[,1]))){
    var1 = p.list[ind, 1] # var1
    var2 = p.list[ind, 2] # var2
    p.value = p.list[ind, 3] # p value
    MAT[var1, var2] = p.value
    MAT[var2, var1] = p.value
  }
  rownames(MAT) = names(df)
  colnames(MAT) = names(df)
  # At this point, MAT has the p.values in a matrix
  
  # subset only the coefficients with p values < 0.05 (or threshold)
  subset = ifelse(MAT < threshold, r.mat, NA)
  rownames(subset) = names(df)
  colnames(subset) = names(df)
  
  output = list(adj.p.values = MAT, threshold.r = subset)
  return(output)
}


lm_eqn <- function(df){
  # from https://stackoverflow.com/questions/7549694/add-regression-line-equation-and-r2-on-graph
  # and originally from: https://groups.google.com/g/ggplot2/c/1TgH-kG5XMA
  m <- lm(y ~ x, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 2)))
  as.character(as.expression(eq));
}

#Define gppr_theme() function

#library(extrafont)
#font_import()
theme_DSQ <- function(){ 
  font <-"arial"# "Nunito"   #assign font family up front
  
  theme_bw() %+replace%    #replace elements we want to change
    
    theme(
      
      #grid elements
      #panel.grid.major = element_blank(),    #strip major gridlines
      #panel.grid.minor = element_blank(),    #strip minor gridlines
      #axis.ticks = element_blank(),          #strip axis ticks
      
      #since theme_minimal() already strips axis lines, 
      #we don't need to do that again
      
      #text elements
      plot.title = element_text(             #title
        family = font,            #set font family
        size = 20,                #set font size
        face = 'bold',            #bold typeface
        hjust = 0,                #left align
        vjust = 2),               #raise slightly
      
      plot.subtitle = element_text(          #subtitle
        family = font,            #font family
        size = 14),               #font size
      
      plot.caption = element_text(           #caption
        family = font,            #font family
        size = 9,                 #font size
        hjust = 1),               #right align
      
      axis.title = element_text(             #axis titles
        family = font,            #font family
        size = 10),               #font size
      
      axis.text = element_text(              #axis text
        family = font,            #axis famuly
        size = 11),                #font size
      
      axis.text.x = element_text(            #margin for axis text
        margin=margin(5, b = 10)),
      
      #strip.text.x = element_text(size = 14)
      
      #since the legend often requires manual tweaking 
      #based on plot content, don't define it here
    )
}

theme_DSQ_SMPC <- function(){ 
  font <-"TT Arial"# "Nunito"   #assign font family up front
  
  theme_bw() %+replace%    #replace elements we want to change
    
    theme(
      
      #grid elements
      #panel.grid.major = element_blank(),    #strip major gridlines
      #panel.grid.minor = element_blank(),    #strip minor gridlines
      #axis.ticks = element_blank(),          #strip axis ticks
      
      #since theme_minimal() already strips axis lines, 
      #we don't need to do that again
      
      #text elements
      plot.title = element_text(             #title
        family = font,            #set font family
        size = 20,                #set font size
        face = 'bold',            #bold typeface
        hjust = 0,                #left align
        vjust = 2),               #raise slightly
      
      plot.subtitle = element_text(          #subtitle
        family = font,            #font family
        size = 14),               #font size
      
      plot.caption = element_text(           #caption
        family = font,            #font family
        size = 9,                 #font size
        hjust = 1),               #right align
      
      axis.title = element_text(             #axis titles
        family = font,            #font family
        size = 13),               #font size
      
      axis.text = element_text(              #axis text
        family = font,            #axis famuly
        size = 13),                #font size
      
      axis.text.x = element_text(            #margin for axis text
        margin=margin(5, b = 10)),
      
      strip.text.x = element_text(size = 14),
      
      legend.text=element_text(size=13),
      legend.title=element_text(size=13)
      
      #since the legend often requires manual tweaking 
      #based on plot content, don't define it here
    )
}

chart.correlation = function(corData, colours_list, Title, Subtitle){
  # create corrected p-values
  correlation = adj.cor(corData, p.adjust = TRUE, p.adjust.method = "BH", threshold = 0.05, cor.method = "kendall")
  pval = correlation$adj.p.values  # shows the p.values in a matrix format
  thresR = correlation$threshold.r  # shows only the p value <0.05
  
  # convert to triangle
  thresR[lower.tri(thresR)]<- NA
  pval[lower.tri(pval)]<- NA
  #create melted values
  melted_cormat<-melt(thresR,value.name = "r")
  melted_pmat<-melt(pval, value.name = "p")
  corPlotData<-full_join(melted_cormat, melted_pmat, by = c("Var1", "Var2"))
  corPlotData$r[corPlotData$r==1]<-NA
  corPlotData.c<-na.omit(corPlotData)
  corPlotData.c$stars<-cut(corPlotData.c$p, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))  # Create column of significance
  corPlotData.c$r<-round(corPlotData.c$r,2)
  
  # Create a ggheatmap
  ggheatmap <- ggplot(corPlotData.c, aes(Var2, Var1, fill = r))+
    geom_tile(color = "white")+
    scale_fill_gradientn(
      colours = colours_list,
      limits = c(-1,1)
    )+
    theme_minimal()+ 
    theme(axis.text.x = element_text(vjust = 1, hjust = 1, angle = 45, size = 6))+
    theme(axis.text.y = element_text(size = 6))+
    labs(title = Title,
         subtitle = Subtitle)+
    coord_fixed()
  
  print(ggheatmap)
  
  ggheatmap2<-ggheatmap+
    geom_text(aes(Var2, Var1, label = paste(stars,r, sep ="\n")), color = "black", size = 2) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 10),
      plot.subtitle = element_text(hjust = 0.5, size = 8),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank(),
      legend.justification = c(1, 0),
      legend.position = c(0.6, 0.8),
      legend.direction = "horizontal")+
    guides(fill = guide_colorbar(barwidth = 5, barheight = 0.4,
                                 title.position = "top", title.hjust = 0.5, 
                                 title.theme = element_text(size = 5),
                                 label.theme = element_text(size = 5)))
  
  print(ggheatmap2)
}

chart.correlation.spearman = function(corData, colours_list, Title, Subtitle){
  # create corrected p-values
  correlation = adj.cor(corData, p.adjust = TRUE, p.adjust.method = "BH", threshold = 0.05, cor.method = "spearman")
  pval = correlation$adj.p.values  # shows the p.values in a matrix format
  thresR = correlation$threshold.r  # shows only the p value <0.05
  
  # convert to triangle
  thresR[lower.tri(thresR)]<- NA
  pval[lower.tri(pval)]<- NA
  #create melted values
  melted_cormat<-melt(thresR,value.name = "r")
  melted_pmat<-melt(pval, value.name = "p")
  corPlotData<-full_join(melted_cormat, melted_pmat, by = c("Var1", "Var2"))
  corPlotData$r[corPlotData$r==1]<-NA
  corPlotData.c<-na.omit(corPlotData)
  corPlotData.c$stars<-cut(corPlotData.c$p, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))  # Create column of significance
  corPlotData.c$r<-round(corPlotData.c$r,2)
  
  # Create a ggheatmap
  ggheatmap <- ggplot(corPlotData.c, aes(Var2, Var1, fill = r))+
    geom_tile(color = "white")+
    scale_fill_gradientn(
      colours = colours_list,
      limits = c(-1,1)
    )+
    theme_minimal()+ 
    theme(axis.text.x = element_text(vjust = 1, hjust = 1, angle = 45, size = 6))+
    theme(axis.text.y = element_text(size = 6))+
    labs(title = Title,
         subtitle = Subtitle)+
    coord_fixed()
  
  print(ggheatmap)
  
  ggheatmap2<-ggheatmap+
    geom_text(aes(Var2, Var1, label = paste(stars,r, sep ="\n")), color = "black", size = 2) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 10),
      plot.subtitle = element_text(hjust = 0.5, size = 8),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank(),
      legend.justification = c(1, 0),
      legend.position = c(0.6, 0.8),
      legend.direction = "horizontal")+
    guides(fill = guide_colorbar(barwidth = 5, barheight = 0.4,
                                 title.position = "top", title.hjust = 0.5, 
                                 title.theme = element_text(size = 5),
                                 label.theme = element_text(size = 5)))
  
  print(ggheatmap2)
}

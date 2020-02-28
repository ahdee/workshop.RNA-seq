library ( heatmap3)
library(d3heatmap)
library (GO.db)
do.histo <- function (result.df){
  
  q.fc <- cut( result.df$logFC, breaks=c(-10,-1.5, -1, 0, 1, 1.5, 10))
  q.fc  <-  data.frame ( table ( q.fc ) ) 
  
  
  q.fdr <-  cut( result.df$adj.P.Val , breaks=c(0.05, .1, .2, .3, .5, .6, 1))
  q.fdr <-  data.frame ( table ( q.fdr) ) 
  q.p <-  cut( result.df$P.Value, breaks=c(0,0.01, .05, 1))
  q.p <-  data.frame ( table ( q.p ) ) 
  
  h1 <- ggplot(data=result.df, aes(result.df$logFC)) +
    geom_histogram(
      binwidth = .2,
      col="grey", 
      fill="#6e9be5", 
      alpha = 1
    ) + 
    labs(title="logFC distrubtion") +
    labs(x="logFC", y="Count")
  
  
  h2 <- ggplot(data=result.df, aes(result.df$P.Value)) +
    geom_histogram(
      binwidth = .01,
      col="grey", 
      fill="#efd37f", 
      alpha = 1
    ) + 
    labs(title="p.value distrubtion") +
    labs(x="p.value", y="Count") +  
    geom_vline(xintercept=.05, color = "red") +
    geom_text(aes(x=.01, label=" p < .05", y=200), colour="red", angle=0 )
  
  h3 <- ggplot(data=result.df, aes(result.df$adj.P.Val)) +
    geom_histogram(
      binwidth = .01,
      col="grey", 
      fill="#bae8a2", 
      alpha = 1
    ) + 
    labs(title="FDR distrubtion") +
    labs(x="corrected p.value", y="Count") +  
    geom_vline(xintercept=.05, color = "red") +
    geom_text(aes(x=.01, label=" p < .05", y=200), colour="red", angle=0 )
  
  return (list(logfc=h1, p.value = h2, fdr=h3, q.fc=q.fc, q.p=q.p, q.fdr=q.fdr ))
}

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


plot.post <- function (result.gene, g1,g2, new.key, r.sub,exp.group="group",exp.this="exp", normal.this="control", sample.id="sample.id",p.val=.05, fdr = .05, fold_thres = 1.5, top10 = NA, GENE_SYMBOL= "GENE_SYMBOL", title1="") {  
  
  #### 
  normal.color <- "#9BB3DF"
  exp.color <- "#DCB335"
  new.key$colorCodes <- normal.color
  new.key[new.key[[exp.group]] == exp.this, ]$colorCodes <- exp.color
  ###
  
  
  d.filter <- result.gene[result.gene$P.Value < p.val & result.gene$adj.P.Val < fdr & (result.gene$logFC > fold_thres | result.gene$logFC < -fold_thres )  ,]
  
  mds.semi <-plotMDS( d.filter[,!colnames(d.filter) %in% r.sub  ]  )
  
  mds.semi <- data.frame(x=mds.semi$x, y=mds.semi$y)
  
  mds.gene <- ggplot( mds.semi , aes(x=x, y=y, shape= new.key[,g1], col=new.key[,g2] ))  +
    geom_point(size=8) +
    
    theme_bw() +
    theme(legend.position="bottom", legend.title=element_blank(), legend.key = element_blank()) +
    labs(x="Dimension 1", y="Dimension 2") + ggtitle(paste0(title1, " MDS Supervised") )   
  
  mds.gene
  
  ### heatmap 
  
  my_palette <- colorRampPalette(c("blue", "#f2f2f2", "#fc8600"))(n = 75)
  par(mar=c(0,0,0,0))
  # bottom, left, top and right 
  
  hm= d3heatmap(as.matrix( d.filter[,!colnames(d.filter) %in% rsub  ]  ), labRow=d.filter$gene, 
                ColSideColors = this.color, distfun = function(x) dist(x,method = 'euclidean'), 
                hclustfun= function(x) hclust(x,method = 'ward.D2')  ) 
  
  
  my_palette <- colorRampPalette(c("blue", "#f2f2f2", "#fc8600"))(n = 75)
  par(mar=c(10,0,0,0))
  # bottom, left, top and right 
  heatmap3(as.matrix( d.filter[,!colnames(d.filter) %in% rsub  ]  ), 
           Rowv=F, dendrogram="column", density.info="none", 
           trace="none", key=TRUE, keysize=1.5, labRow=NA, cexCol=1.0, scale="row", 
           col = my_palette, showRowDendro=TRUE,   margins=c(8,2),
           main= paste0(title1, " Heatmap" ), ColSideColor = new.key$colorCodes, ColSideLabs = "")
  legend("right",legend=c(normal.this,exp.this),fill=c(normal.color,exp.color) )
  
  hm2 <- recordPlot()
  plot.new() ## clean up device
  
  
  
  
  
  #### Scatterplot here 
  
  # scatter plot now
  ## prepping scatter plot   
  
  result.gene$class <- 'no-change'
  result.gene[ result.gene$P.Value < p.val & (result.gene$logFC > fold_thres) & result.gene$adj.P.Val < fdr, ]$class <- 'up'
  result.gene[ result.gene$P.Value < p.val & (result.gene$logFC < -fold_thres) & result.gene$adj.P.Val < fdr, ]$class  <- 'down'
  
  
  ##############
  
  
  overexp_col <- '#fc8600'
  underexp_col <- 'blue'
  nochange_col <- '#F7F6F5'
  
  
  c <- normal.this
  e <- exp.this
  
  # name of average 
  ave.control <- paste0(c,".ave")
  ave.exp <- paste0(e,".ave")
  ave.exp = gsub("-",".", ave.exp)
  
  # get all the colname from key sample.id needs to match that of data
  id.control <- as.character ( new.key[new.key[[exp.group]] %in% normal.this, sample.id] )
  id.exp <- as.character ( new.key[new.key[[exp.group]] %in% exp.this, sample.id] )
  
  result.gene[[ave.exp]] <- rowMeans(result.gene[, id.exp  ] )
  result.gene[[ave.control]] <- rowMeans(result.gene[, id.control ] )
  
  result.gene$AVE <- rowMeans( result.gene[, c(  id.control , id.exp )] )
  
  
  
  class <- "class"
  result.gene <- result.gene[order(result.gene$adj.P.Val, result.gene$P.Value),]
  
  if ( is.na (top10) ){
    top10 <- head (result.gene[ !result.gene$class == "no-change" & !result.gene[[GENE_SYMBOL]] == "",  GENE_SYMBOL ],10 )
    top10 <- unique (top10)
  }
  
  label.this <- result.gene[ result.gene[[GENE_SYMBOL]] %in% top10, ]
  
  scatterplot1 <- ggplot(result.gene, aes_string(ave.control ,ave.exp, col="class"))  + 
    scale_color_manual(values=c("up"=overexp_col, "no-change"=nochange_col, "down"=underexp_col), breaks=c("up", "down", "no-change"), labels=c(paste0("up in ",e,"(",nrow(result.gene [result.gene $class=='up',]),")"), paste0("down in ", e ,"(",nrow(result.gene[result.gene $class=='down',]),")"), "no-change")) +
    
    
    geom_point() +
    theme_bw() +
    theme(legend.position="bottom", legend.title=element_blank(), legend.key = element_blank()) +
    labs(x=paste0(c ," average (log2 fpkm)" ), y= paste0(e ," average (log2 fpkm)" ) ) + 
    ggtitle(paste0(title1, "Scatter.plot ")) +
    geom_text_repel(
      data = label.this, 
      aes(
        label = factor ( label.this[[GENE_SYMBOL]])
      ), 
      
      point.padding = unit(.55, "lines"),
      box.padding = unit(2.25, "lines"),
      nudge_y = 0.5,
      fontface=2, 
      color="black"
    )
  
  ### 
  
  
  
  
  ma_p_h <- ggplot(result.gene , aes(x=AVE,y=logFC, col=factor(class) ) )+ 
    scale_color_manual(values=c("up"=overexp_col, "no-change"=nochange_col, "down"=underexp_col), breaks=c("up", "down", "no-change"), labels=c(paste0("up in ",e,"(",nrow(result.gene [result.gene $class=='up',]),")"), paste0("down in ", e ,"(",nrow(result.gene[result.gene $class=='down',]),")"), "no-change")) +
    scale_size_manual(guide="none", values=c("N"=1, "Y"=5)) +
    scale_alpha_manual(guide="none", values=c("N"=0.3, "Y"=1)) +
    #scale_shape_manual(breaks='Y',label='KRAS',values=c("Y"=17,"N"=16))+
    geom_point() +
    theme_bw() +
    theme(legend.position="bottom", legend.title=element_blank(), legend.key = element_blank()) +
    labs(x="Ave Expression", y="FC-log2") +  
    ggtitle(paste0(title1, " MA ")) +
    geom_text_repel(
      data = label.this, 
      aes(
        label = factor ( factor ( label.this[[GENE_SYMBOL]])  )
      ), 
      
      point.padding = unit(2.55, "lines"),
      box.padding = unit(2.25, "lines"),
      nudge_y = 0.5,
      fontface=1, 
      color="black"
    ) 
  
  ### volcano 
  
  
  volcano.g <- ggplot(result.gene , aes(logFC, -log10(P.Value))) +
    geom_point(aes(col=class)) +
    scale_color_manual(values=c(underexp_col, 'black',overexp_col) , name="sig", breaks = c('up','down','no-change'), labels = c(paste0("up in ",e,"(",nrow(result.gene [result.gene $class=='up',]),")"), paste0("down in ", e ,"(",nrow(result.gene[result.gene $class=='down',]),")"), "no-change") ) + 
    ggtitle(paste0(title1, " Volcano.plot ")) +
    geom_text_repel(
      data = label.this, 
      aes(
        label = factor ( label.this[[GENE_SYMBOL]])
      ), 
      fontface=2, 
      color="red",
      
      point.padding = unit(.25, "lines"),
      box.padding = unit(.25, "lines"),
      nudge_y = 0.1
    ) +theme_bw(base_size = 20)  +  theme(legend.text=element_text(size=15)) + guides(color = guide_legend(override.aes = list(size=10)))
  
  
  pm <- list ( 
    mds = mds.gene, 
    hm = hm, hmstat = hm2, 
    scatter = scatterplot1, 
    ma = ma_p_h, 
    volcano = volcano.g, 
    data = result.gene
  )
  
  return (pm)
  
}


get.go <- function (g,mmid,g.name="GENE_SYMBOL", m.name="external_gene_name", e.name="entrezgene", species="Hs",bg, fdr=.05){
  # hg38 must be defined 
  # data frame must contain class to define up and down
  
  g <- merge(g,mmid, by.x=g.name, by.y=m.name )
  
  tm.this <- function (x){
    x$Term <- strtrim(x$Term,60)
    return (x)
  }
  
  # new version will let user decide what is up or down instead of relying on "class"
  #up <- unique ( g[g$class == "up",e.name] )
  #down <-unique (  g[g$class == "down",e.name] )
  up <- unique ( g[g$logFC > 0 ,e.name] )
  down <-unique (  g[g$logFC < 0 ,e.name] )
  
  # bg is now force to have input so you need to provide a background! 
  bg = merge(data.frame(bg=bg),mmid, by.x="bg", by.y=m.name )
  bg <- unique ( mmid[ , e.name] )
  
  
  
  go <- goana(list(Up=up, Down=down),
              species=species,
              universe=bg, FDR=fdr)
  
  
  
  
  # note that I commented out fdr, reason is because fdr already screened above
  
  
  bp <- tm.this ( topGO(go, 'BP', number=Inf) )
  #bp$P.Upfdr = p.adjust(bp$P.Up , n = nrow( bp) ) 
  #bp$P.Downfdr = p.adjust(bp$P.Down , n = nrow( bp) )
  bp = bp[bp$P.Down < .05 | bp$P.Up < .05, ]
  #bp = bp[bp$P.Downfdr < fdr | bp$P.Upfdr < fdr, ]
  bp$P.Up = ifelse ( bp$P.Up < .05, "*", "ns")
  bp$P.Down = ifelse ( bp$P.Down < .05, "*", "ns")
  
  
  
  cc <- tm.this ( topGO(go, 'CC', number=Inf) )
  #cc$P.Upfdr = p.adjust(cc$P.Up , n = nrow( cc) ) 
  #cc$P.Downfdr = p.adjust(cc$P.Down , n = nrow( cc) )
  cc = cc[cc$P.Down < .05 | cc$P.Up < .05, ]
  #cc = cc[cc$P.Downfdr < fdr | cc$P.Upfdr < fdr, ]
  cc$P.Up = ifelse ( cc$P.Up < .05, "*", "ns")
  cc$P.Down = ifelse ( cc$P.Down < .05, "*", "ns")
  
  
  mf <- tm.this ( topGO(go, 'MF', number=Inf) )
  #mf$P.Upfdr = p.adjust(mf$P.Up , n = nrow( mf) ) 
  #mf$P.Downfdr = p.adjust(mf$P.Down , n = nrow( mf) )
  mf = mf[mf$P.Down < .05 | mf$P.Up < .05, ]
  #mf = mf[mf$P.Downfdr < fdr | mf$P.Upfdr < fdr, ]
  mf$P.Up = ifelse ( mf$P.Up < .05, "*", "ns")
  mf$P.Down = ifelse ( mf$P.Down < .05, "*", "ns")
  
  kegg <- kegga(list(Up=up, Down=down),
                species=species,
                universe=bg, FDR=fdr)
  
  kegg = kegg[(kegg$P.Up < .05 | kegg$P.Down < .05 ), ]
  keggUP = kegg[order(kegg$P.Up), ]
  keggDown = kegg[order(kegg$P.Down), ]
  kegg$p = ifelse ( kegg$P.Up < kegg$P.Down, kegg$P.Up, kegg$P.Down)
  # the last one makes sure that p value is not specific to up or down rather it capture both and you can subset later
  keggn = kegg[ order(kegg$p), ]
  keggn$P.Up = ifelse ( keggn$P.Up < .05, "*", "ns")
  keggn$P.Down = ifelse ( keggn$P.Down < .05, "*", "ns")
  
  list ( bp=bp, cc=cc, mf=mf, kegg=keggn #, 
         #k.up = k.up, k.down = k.down
  )
  
}   
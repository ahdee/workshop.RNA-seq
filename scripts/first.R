# R data structure
## think of it a series of "objects"
# Most basic is the atomic vector 


## there are three types here 

1 # numeric 
"cat" # character 
TRUE # logical 

# you can store a single object or multiple  
a=10 
b=c(2,3,3, 10, 10000)

# you can print and do stuff to them 
b[1:3]
b + 12 


# you can name your vector 
names ( b ) = c("a","b","c","d","e")
b["e"]
b["e"] + 12



# logical vectors can be build from comparison 
a > b 
a >= b 
a == b 
a != b
a %in% b 

# you can use the results of your logical vector to subset things 

b[ a > b  ]
# you can apply logical operations 

b[ a>b & b<3] # AND
b[ b > 100 | b<3] # OR 
b[ ! ( b > 100 | b<3 )] # NEGATION !

# you cannot coerce a numeric type into a character and vice versa 
c = c ( 'cat', 'dog')
c + a
as.numeric ( c )

# sometimes things can be deceptive 
d= c("1","2") # this looks like numbers! 
a + d 
# notice that d is looks like number but its not! 
# you can check this
str ( d )
str ( a )

# you can try to convert this. 
as.numeric ( d ) + a



# vectors needs to be of the same type! 

## lists 
# however unlike vector a list can have different types 

x = list ( a = c(1,2,3), b = c(TRUE,FALSE), c= c ( "cats","dogs"))
x[["a"]]


# you can create a 2 dimentional atomic vector represented as a matrix 
set.seed(123) # make this repeatable 
m1 <- matrix( rnorm(15), ncol = 3, nrow = 5) # specify how many rows and columns you need. 
# you can name your columns and rows. 
colnames ( m1 ) = c("cond_1","cond_2","cond_3" )
row.names ( m1 ) = c( "g1","g2","g3","g4","g5")

# you can change the value(s) of a specific cell 
m1["g1","cond_2"] = 10
m1 


# The power of R:  you can do math on any dimension you want.  
rowSums(m1)
colSums(m1)

# NOTE ONE SAYING IN R IS THAT THERE is yet another way to do it! 

# you can add your function to each row 
apply(m1, 1, mean)
# multiple 2 to each number in the matrix 
apply(m1, 2, function(x) {2*x} ) 

# dataframe 
# similar to matrix but it can take mix datatype. In other words 1 column can be character another numbers
df= data.frame ( x = c(1:3), y=c("cond_1","cond_2","cond_3"), stringsAsFactors = F )
str ( df )

# you can convert a matrix to a data.frame 
df = data.frame ( m1 )
# add to column 
df$gender = c("M","F","F","M","M")
df$pid = paste0 ( "sample_" , seq ( 1, 5, 1) )

# you can copy data.frame just like any other objects 
df2 = df 

# you can add rows
df=rbind ( df , data.frame ( cond_1=.2 
                             ,cond_2=1
                             ,cond_3=-.3
                             ,gender="M"
                             ,pid="sample_1") )


df = cbind ( df, data.frame ( random = letters[1:nrow (df )] ) )

# your first package 
## increasing ability to manipulate dataframe with 
library(dplyr)
 
# suppose you want to round all the numbers 

df

df %>% 
  mutate_if(is.numeric, round, 3 )

#  apply a function 
df %>%
  mutate(cond_3mod = log2 ( cond_3 + 1 )  )

# apply on a subset 
df %>% 
  mutate_if(is.numeric, round, 3 ) %>%
  filter ( gender == "M")

# aggregate 
df

df %>% group_by(pid ) %>%
  summarise(
    random2 = paste(unique ( random ), collapse = "," ) ,   
    mean = mean(cond_3)
  ) %>%
  data.frame()

# lots lots more.  If you need to manipulate google dplry ! 




### the concept of wide vs long 
df = df2 
df 
# difference 
# suppose a - c represents different conditions and you want to reformat this such that 
# gender, pid, cond_value 


library ( reshape2)



# most importantly you need to specify the variables you want to keep 
# Specify id.vars: the variables to keep but not split apart on
# in this case its gender, pid 
dmelt = reshape2::melt ( df, id.vars = c("pid","gender"))
# your variable are the columns that were melted. 
# you can reorder 
dmelt = dmelt [ order ( dmelt$pid, dmelt$variable), ]

# take a look at the dataframe and you will see that the rows for each condition have collapse
head ( dmelt )

# now we can convert long back to wide 

# "pid" and "gender" are columns we want to keep the same
# "variable" is the column that contains the names of the new column we want to expand 
# "value" holds the measurements
dcast(dmelt, pid + gender ~ variable, value.var="value")
# compare that with df 
df



str ( dmelt )
# notice there is a row that reads Factor. This could be your best friend or worse enemy. 


library ( ggplot2 )
# works only with dataframes. 
# utilizes factors for order 
# build by layers 
# aes function tells ggplot what information you want

g= ggplot(dmelt, aes(x=pid, y=value)) 

g 
# notice that you get a blank plot.  This is because ggplot now expects you add another layer telling it what type of plot you want. 

# now lets add a dot layer 
g  + geom_point( )

# however you may want to specify by color what the conditions are
# again add to aes 
g  + geom_point(aes ( col=variable) )
# lets change the size a bit 
g  + geom_point(aes ( col=variable) , size=15) 
# lets add shape to define male or female 
g2 = g  + geom_point(aes ( col=variable, shape=gender),  size=12) 

g2

# lets use our own colors 
g2 + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))
# you can literally just use normal names 
g2 + scale_color_manual(values=c("green", "steelblue", "red"))

# you can use a build-in pallete 
library(RColorBrewer)

# color blind friendly 
display.brewer.all(colorblindFriendly = TRUE)
display.brewer.pal(n = 8, name = 'Dark2')
g3 = g2 + scale_color_brewer(palette="Dark2")
g3

# themes can be apply to apply overall look to plot 
g3 + theme_minimal()  
g3 + theme_classic() 
g3 + theme_bw()

# you can customize the theme elements 
g4=g3 + theme_classic() + theme( 
  axis.text.x=element_text(size=12, angle = 90, vjust=.5) # x axis and tilt the labels
  ,axis.text.y = element_text(size=12) 
  , legend.position="bottom" # move the legend 
  , axis.title = element_text( size= 20 )
  ) + 
  xlab ( "" ) + # remove x title 
  ylab ( "condition value ") + # change y
  ggtitle ( "first plot ") # add or change title 

g4

# recall that you can layer ggplot however you want. 

g4 +  geom_boxplot(alpha = 0.8)


# but how does ggplot know what the order is? 
# use factors 
dmelt$pid = factor ( as.character(dmelt$pid), levels= rev ( unique ( dmelt$pid) )   )
# you can see the order 
dmelt$pid
ggplot(dmelt, aes(x=pid, y=value))  + geom_point(aes ( col=variable) , size=10)

#### there is a lot more here but we will stop here. 
# google ggplot plot x? 

# last basic lesson here 
# the concept of loops in R 

for ( i in df$cond_1){
  i = i + 3
  print ( i )
}



## ok so now lets finish this lesson off by talking about the scalibility issue
# lets use Marcela's ewing sarcoma cell line as an example 
# import key 

ekey = readRDS(  "~/workshop.RNA-seq/seq/testout/ews.key.rds" )
ekey 

# get a list of files 
base="~/workshop.RNA-seq/seq/testout/testloop/"
out = "~/workshop.RNA-seq/scripts/output/"
files = list.files( base, pattern = "fastq")

files 

cmd = "sbatch ~/workshop.RNA-seq/scripts/workflow1.sh"
# Recall this is the the function to invoke the job 
# sbatch workflow1.sh 
# where is our raw directory, 
# where do you want your outputs
# name of your file up to R1
# what you want to call it. 

# also we only need just R1 as mentioned previously 
files = files [ grepl ( "R1", files )]

files 


log = ""
for ( f in files ){
  # lets first grab the tube number and leave if it does not match our key 
  RNAseq.id  = stringr::str_match(f, "(.*?)_")[, 2] # recall our naming schema the RNAseq.id  is always the beginning! 
  # check if this is one the sample we need to analyze 
  if ( isTRUE ( any ( ekey$RNAseq.id %in% RNAseq.id ) )){
    # now lets get the base name - recall the script will create input by attaching .fastq.gz at the end so all we need is everything before that 
    # in1="${base}${inputfile}.fastq.gz
    basename = stringr::str_match(f, "(.*).fastq.gz")[, 2]
    
    # now can output and submit the job immediately or save it as a single file 
    cmthis = paste ( cmd, base, out, basename, RNAseq.id )
    # you can do a system call here 
    # or save it 
    log = c( log, c( cmthis, "\n") )
  }
}

# you can now save this to file and then just invoke it with sh xxx 
cat ( log )
#write(log , paste0( "new_loop.sh"), append = F)

# done you now have a fully functional scalable workflow for RNAseq with all the QC you need. 
# now lets go back to the lecture and see what we got. 


### ok RNAseq analysis workflow 
library ( limma )
library ( edgeR )
library ( data.table )
# lets load some cheat codes! 
source ( "~/workshop.RNA-seq/scripts/cheat.R" ) 

# import raw data - I saved this previously in an rds object 
# this is 2D vs 3D ANKR35
data = readRDS( "~/workshop.RNA-seq/rnaseq/counts/ank.rds")

# lets study what we just imported 
names ( data )

key = data$key

# let check out the key 
key 
# so as you can see there are two factors to consider here. 
# plate ( 2D vs 3D )
# genetics ( ANKRD35 vs GFP )


match.data <- data$data # lets import the raw counts 
head ( match.data )
# importantly we need to make sure that the names of the key and the dataframe order is the same! 
match.data =  match.data [ , key$pid]

# sanity check is always necessary and should be done often 
all.equal( colnames ( match.data), key$pid ) # look at what happens when I reverse the pid 
all.equal( colnames ( match.data), rev ( key$pid ) )

dim ( match.data )
dim ( key )

# now lets import this into edgeR 
y.match <- DGEList( match.data)



## remove low counts since this could interfere with some of the statistical approximations
#  here half the total sample size must have > 1 raw to be included 
threshold <- ncol(match.data)/2
keep <- rowSums( cpm(y.match)  > 1) >= threshold


# sanity, however because this was a  Ankrd35 perturbation you may expect low counts so here we force the keep Ankrd35
# we are going to force keep Ankrd35 even if its low.
keep["Ankrd35"] # see this would had been removed! 
keep["Ankrd35"] = TRUE
keep["Ankrd35"] 

# before 
dim(y.match)
# after 
y.match  <- y.match [keep, ]
dim(y.match)

## normalisation using trimmed mean of M-values (TMM) (Robinson and Oshlack 2010)
## in additional to library scaling this will ensure we don't get some strange rna composition issues 
y.match <- calcNormFactors( y.match, method = "TMM" )


## we include another normalization step so that we can input and use limma.  This step will also help stabalize the variance 
# however before we do its better to setup our model matrix. This will be use for linear model
group = as.factor ( key$combo )
design <- model.matrix(~0  + group )
# clean up the name a bit 
colnames ( design ) = gsub ( "group","", colnames ( design ))
voom <- voom(y.match, design = design, plot=T) 

# now lets get the cpm 
cpm <-  data.frame (  cpm(y.match, log=T, prior.count=1)  ) 

# ok so now lets take a look at our alignment and normalization 
# we've already seen that the mean-variance trend looks good however lets take a look at the density plots 
# get good colors 
getPalette = colorRampPalette(brewer.pal(9, "Set1")) # expand color pallete

temp = melt (cpm ) # remember that ggplot only does long form 
colnames(temp) = c("sample", "log2.cpm")
g = ggplot(temp, aes(log2.cpm, fill=sample, colours=sample)) + geom_density(alpha=0.1) + 
  theme_bw()  + scale_color_manual( values= getPalette(12))
plot ( g )
# density plots are good way to see if there are any outlier samples. 

boxplot(cpm, las=2, col="lightblue", main="log2 (cpm + 1) after normalization")

# now lets take a look at some of its qc values 
head ( data$qc)
# here we map important values like library size.., total reads, unique mapped reads etc. 
# here is an example 
# but first we need to melt this. 

# Specify id.vars: the variables to keep but not split apart on
# in this case its gender, pid 
qc = reshape2::melt ( data$qc, id.vars = c("sample","dataset"))
# lets select out some to plot 
qc = qc[qc$variable %in% c( "Library.size..Mil.", "Mapped.reads", "Total.reads","Uniquely.mapped.reads"), ]


ggplot(qc[qc$variable=="Total.reads", ], aes(x=sample , y=value)) +
  geom_bar( stat = "identity", position = position_dodge(), colour='grey' ) +
  theme(legend.position="none", legend.title=element_blank(), legend.key = element_blank(),
        axis.text.x = element_text(angle = 90, size=8),
        axis.text.y = element_text(size=8)
  ) 

# however this is not very super useful; suppose we want to split the groups up between 2d and 3d 
# first merge the info 
qc = merge ( qc, key, by.x= "sample", by.y="pid")

ggplot(qc[qc$variable=="Total.reads", ], aes(x=sample , y=value)) +
  geom_bar( aes ( fill=plate), stat = "identity", position = position_dodge() ) +
  theme(legend.position="none", legend.title=element_blank(), legend.key = element_blank(),
        axis.text.x = element_text(angle = 90, size=8),
        axis.text.y = element_text(size=8)
  ) 

# that is a little better but lets rank it now.  This is when factor is very important 
qc = qc[ order ( qc$plate, qc$variable, qc$value), ]

ggplot(qc[qc$variable=="Total.reads", ], aes(x=sample , y=value)) +
  geom_bar( aes ( fill=plate), stat = "identity", position = position_dodge() ) +
  theme(legend.position="none", legend.title=element_blank(), legend.key = element_blank(),
        axis.text.x = element_text(angle = 90, size=8),
        axis.text.y = element_text(size=8)
  ) 

# so notice the ranking did nothing! 
# this is why factors are so important. 
qc$sample = factor ( as.character ( qc$sample), levels= unique ( qc$sample))


ggplot(qc[qc$variable=="Total.reads", ], aes(x=sample , y=value)) +
  geom_bar( aes ( fill=plate), stat = "identity", position = position_dodge() ) +
  theme(legend.position="none", legend.title=element_blank(), legend.key = element_blank(),
        axis.text.x = element_text(angle = 90, size=8),
        axis.text.y = element_text(size=8)
  ) + scale_fill_manual ( values=getPalette(12))


### however we can do better by plotting everything all at once 

ggplot(qc, aes(x=sample , y=value)) +
  geom_bar( aes ( fill=plate), stat = "identity", position = position_dodge() ) +
  theme(legend.position="none", legend.title=element_blank(), legend.key = element_blank(),
        axis.text.x = element_text(angle = 90, size=8),
        axis.text.y = element_text(size=8)
  ) + scale_fill_manual ( values=getPalette(12))+ facet_wrap(~variable, scales = "free_y" )


# ok so these samles have some pretty good qc's so now lets take a look at the exeriment QC! 
# lets look for batch effect and/or study our experimental distribution 

library ( FactoMineR )
library ( factoextra )
library ( ggrepel )

res.pca <- PCA( t( cpm ), graph = FALSE) ## Principal Component Analysis 
# plot the PCA
p <- fviz_pca_ind(res.pca, label="none", addEllipses=TRUE, ellipse.level=0.95,  axes = c(1, 2)) # convert this to ggplot object
p2 <- p +geom_point(aes(colour= key$plate, shape= key$genetic ), size=8.0, alpha=.5  ) # add geom pts

# add pid annotation so we can see outliers 
pca.unc = p2 + geom_text_repel(aes(label = key$pid ), arrow = arrow(length = unit(0.01, 'npc')), box.padding = unit(1.5, 'lines')   ) +
  scale_colour_discrete(name  ="Plate") + scale_shape_discrete ( name="Genetics") + ggtitle("PCA uncorrected")

# this is as perfect as you can hope for in an experiment! 
# see that the 2D and 3D defines the majority of the variance ( x axis Dimension 1 )
# however the D2 also shows that the seperation is also occurs at the genetic level.

# this plots shows that microenviroment drives most of the variance 
# with a secondary effect from genetics 

pca.unc 

# from here we can explore what the main genes are driving the variance. 
# lesson here is that DGE is not the only method! 
fviz_contrib(res.pca, choice = "var", axes = 1,  top = 20)
fviz_contrib(res.pca, choice = "var", axes = 1:2,  top = 20)

# here is an example of what happens when there is a batch effect 

# here I created a scenerio where we want to study the difference between aspirin and sugar pills 
# this is from the same persons 


bad = readRDS( "~/workshop.RNA-seq/rnaseq/counts/batc.rds")
bad.data = bad$new.data
bad.key = bad$new.key


res.pca2 <- PCA( t( bad.data ), graph = FALSE)
# plot the PCA
p <- fviz_pca_ind(res.pca2, label="none", addEllipses=TRUE, ellipse.level=0.95,  axes = c(1, 2))
p2 <- p +geom_point(aes(colour= bad.key$drug, shape= factor(bad.key$lane) ), size=5.0  )
pca.bad = p2 + geom_text_repel(aes(label = bad.key$pid ), arrow = arrow(length = unit(0.01, 'npc')), box.padding = unit(1.5, 'lines')   ) +
  scale_colour_discrete(name  ="Drug") + scale_shape_discrete ( name="Lane") + ggtitle("PCA uncorrected")

pca.bad 
# so as you can see the lane difference is driving most of the differences! not the sugar or aspirin that is secondary! 
# here there is clearly and observable difference from something that is NON-biological and thus the batch most be removed! 

library (sva)
combat <- ComBat(dat=as.matrix( bad.data ), batch=factor ( bad.key$lane ), mod=NULL, par.prior=TRUE, prior.plots=FALSE)

combat.log = data.frame ( log2 ( 2^combat + 1 ), stringsAsFactors = FALSE)

# plot the PCA, corrcted to see 
# wow! batch has been removed 
p <- fviz_pca_ind( PCA( t( combat), graph = FALSE) , label="none", addEllipses=TRUE, ellipse.level=0.95,  axes = c(1, 2))
p2 <- p +geom_point(aes(colour= bad.key$drug, shape= factor(bad.key$lane) ), size=5.0  )
p2 + geom_text_repel(aes(label = bad.key$pid ), arrow = arrow(length = unit(0.01, 'npc')), box.padding = unit(1.5, 'lines')   ) +
  scale_colour_discrete(name  ="Drug") + scale_shape_discrete ( name="Lane") + ggtitle("Corrected for cell line")

# warning!!!! when running the linear model do not use the batch corrected values, instead run that as covariate. 

bad.key
# so in your design you would add the batch as a covariate 
 model.matrix(~0  + bad.key$drug + bad.key$lane)
# however off note since pid_1 was persumbably the same person we need to use a pair-wise comparison - this will increse your power
 model.matrix(~0  + bad.key$drug + bad.key$lane + bad.key$pid)


# ok finally now we run our linear models and extract DGEs
design 
# first lets ask what we want to get from it. 
# 1. compare all D3 to D2 
# 2. cmpare within D3: ANKRD35 vs GFP 
# 3. cmpare within D2: ANKRD35 vs GFP 

## the way these models works regardless if its edgeR or limma is the concept of a contrast.
# you use this to guide the model to extract what comparisons you want. 

contrast.matrix.all <- makeContrasts( 
  D3 = D3.ANKRD35 - D3.GFP # first variable is always relative to that. 
  ,D2 = D2.ANKRD35 - D2.GFP # first variable is always relative to that. 
  , D3D2 =   ( D3.ANKRD35 + D3.GFP )/2 -  ( D2.ANKRD35 + D2.GFP )/2 # this is a more complicated and we "add" the variables we want but we need to divide by 2 
  , levels=design 
)


fit  <- lmFit(voom, design ) # Fit linear model for each gene given a series of arrays
fit.r <- contrasts.fit(fit, contrast.matrix.all) # Given a fitted linear modelcompute estimated coefficients and standard errors for a given set of contrasts.
fit.r <- eBayes(fit.r) #  compute moderated t-statistics, moderated F-statistic, and log-odds of differential expression by empirical Bayes

# now to get the DGE genes is as simple as getting the coef you named above 

D3D2 <- topTable(fit.r, coef="D3D2", n=Inf )
D2 <- topTable(fit.r, coef="D2", n=Inf )
D3 <- topTable(fit.r, coef="D3", n=Inf )

hist.D3D2<- do.histo ( D3D2 )
hist.D3<- do.histo ( D3 )
hist.D2<- do.histo(D2)

multiplot( hist.D3D2$logfc, hist.D3D2$fdr, ncol=1)



D3D2.final  <- merge( D3D2 , cpm  , by="row.names" )
colnames( D3D2.final  )[1] <- "gene"
rsub <- colnames(D3D2.final )[1:7]


# now we do some post QC and explorations 

pv = .05
fdr = .05
logfc = 1.5 


color.env2 =  brewer.pal(n = length ( unique(key$plate )) + 1, name = "Set3")
names ( color.env2 )= unique ( key$plate)


this.color = color.env2[key$plate]

post.D3D2 <- plot.post ( D3D2.final, "plate", "plate" , key , rsub, exp.group="plate", exp.this="D3", normal.this= "D2", sample.id="pid", GENE_SYMBOL = "gene", fdr=fdr, fold_thres = logfc, title1 = "2D vs 3D ", top10 = NA)

go.D3D2 <- get.go (
  post.D3D2$data[post.D3D2$data$class != "no-change", ] , data$mouse ,g.name="gene",
  bg = row.names( post.D3D2$data ), species = "Mm"
) 









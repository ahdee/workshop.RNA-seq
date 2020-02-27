# R data structure
## think of it a series of "objects"
# Most basic is the atomic vector 


## there are three types here 

1 # numeric 
"cat" # character 
TRUE # logical 

# you can store them individually or together 
a=10 
b=c(2,3,3, 10, 10000)

# you can name your vector 
names ( b ) = c("a","b","c","d","e")
b["e"]



# you can print and do stuff to them 
b[1:3]
b + 12 

# logical vectors can be build from comparison 
a > b 
a >= b 
a == b 
a != b
a %in% b 

# you can use the results of your logical vector to subset things 

b[ a > b  ]
# you can apply logical operations 

b[ a>b & b<3]
b[ b > 100 | b<3]
b[ ! ( b > 100 | b<3 )]

# you cannot coerce a numeric type into a character and vice versa 
c = c ( 'cat', 'dog')
c + a
as.numeric ( c )

# sometimes things can be deceptive 
d= c("1","2")
a + d 
# notice that d is looks like number but its not! 
# you can check this
str ( d )
str ( a )

# you can try to convert this. 
as.numeric ( d ) + a



# last all vectors needs to be of the same type! 

# lists 
# however unlike vector a list can have different types 

x = list ( a = c(1,2,3), b = c(TRUE,FALSE), c= c ( "cats","dogs"))
x[["a"]]


# you can create a 2 dimentional atomic vector represented as a matrix 
set.seed(123) # make this repeatable 
m1 <- matrix( rnorm(15), ncol = 3, nrow = 5)
# you can name your columns and rows. 
colnames ( m1 ) = c("cond_1","cond_2","cond_3" )
row.names ( m1 ) = c( "g1","g2","g3","g4","g5")

# you can change the values 
m1["g1","cond_2"] = 10
m1 
# you can do math 
rowSums(m1)
colSums(m1)

# you can add your function to each row 
apply(m1, 1, mean)
# multiple 2 to each number in the matrix 
apply(m1, 2, function(x) {2*x} ) 

# dataframe 
# similar to matrix but it can take mix datatype. In other words 1 row can be character another numbers
df= data.frame ( x = c(1:3), y=c("cond_1","cond_2","cond_3"), stringsAsFactors = F )
str ( df )
df = data.frame ( m1 )
# add to column 
df$gender = c("M","F","F","M","M")
df$pid = paste0 ( "sample_" , seq ( 1, 5, 1) )

# difference between wide and long 
# suppose a - c represents different conditions and you want to reformat this such that 
# gender, pid, cond_value 

# your first package 
library ( reshape2)
# most importantly you need to specify the variables you want to keep 
# Specify id.vars: the variables to keep but not split apart on
# in this case its gender, pid 
dmelt = reshape2::melt ( df, id.vars = c("pid","gender"))
# your variable are the columns that were melted. 
# you can reorder 
dmelt = dmelt [ order ( dmelt$pid, dmelt$variable), ]




# "pid" and "gender" are columns we want to keep the same
# "variable" is the column that contains the names of the new column to put things in
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

g2 + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))
# you can literally just use normal names 
g2 + scale_color_manual(values=c("green", "steelblue", "red"))

# you can use a build-in pallete 
library(RColorBrewer)

# color blind friendly 
display.brewer.all(colorblindFriendly = TRUE)
display.brewer.pal(n = 8, name = 'Dark2')
g3 = g2 + scale_color_brewer(palette="Dark2")

# themes can be apply to apply overall look to plot 
g3 + theme_minimal()  
g3 + theme_classic() 
g3 + theme_bw()

# you can customize the theme elements 
g4=g3 + theme_classic() + theme( 
  axis.text.x=element_text(size=12, angle = 90, vjust=.5)
  ,axis.text.y = element_text(size=12) 
  , legend.position="bottom"
  , axis.title = element_text( size= 20 )
  ) + 
  xlab ( "" ) + # remove x title 
  ylab ( "condition value ") + # change y
  ggtitle ( "first plot ") # add or change title 


# recall that you can layer ggplot however you want. 

g4 +  geom_boxplot(alpha = 0.8)

#### there is a lot more here but we will stop here. 





# lets load some packages to read in file and manipulate data 
library ( readr) 
library ( data.table )
library ( dplyr )
library ( reshape2)





output = "~/workshop.RNA-seq/scripts/work/"
# check if this exists 
dir.exists( output )
if ( ! dir.exists( output ) ){
  dir.create( output )
}

dir.exists( output )


script <- read_file("~/workshop.RNA-seq/scripts/loop.sh") # read in script 
cat ( script ) # print out file 

scriptnew = gsub ( "LOOP", "10 11 12 CAT", script)


write(scriptnew , paste0(output,  "new_loop.sh"), append = F)


# where is the raw data? 
# the idea here is that we are going to loop through each one 







data = readRDS( "batc.rds" )
new.data = data$new.data
new.key = data$new.key


res.pca <- PCA( t( new.data ), graph = FALSE)
# plot the PCA
p <- fviz_pca_ind(res.pca, label="none", addEllipses=TRUE, ellipse.level=0.95,  axes = c(1, 2))
p2 <- p +geom_point(aes(colour= new.key$drug, shape= factor(new.key$lane) ), size=5.0  )
pca.unc = p2 + geom_text_repel(aes(label = new.key$pid ), arrow = arrow(length = unit(0.01, 'npc')), box.padding = unit(1.5, 'lines')   ) +
  scale_colour_discrete(name  ="Drug") + scale_shape_discrete ( name="Lane") + ggtitle("PCA uncorrected")

# base on above and discussion we find that its crucial to correct for cell type. here we use a null to be as unbias as possible. 




combat <- ComBat(dat=as.matrix( new.data ), batch=factor ( new.key$lane ), mod=NULL, par.prior=TRUE, prior.plots=FALSE)
#combat = setcolorder(combat,key$tube)
#colnames ( v.data )= key$pid
#combat <- ComBat(dat=as.matrix( voom ), batch=factor ( key$cell ), mod=NULL, par.prior=TRUE, prior.plots=FALSE)
combat.log = data.frame ( log2 ( 2^combat + 1 ), stringsAsFactors = FALSE)

# plot the PCA, corrcted to see 
res.pca <- PCA( t( combat), graph = FALSE)

p <- fviz_pca_ind(res.pca, label="none", addEllipses=TRUE, ellipse.level=0.95,  axes = c(1, 2))
p2 <- p +geom_point(aes(colour= new.key$drug, shape= factor(new.key$lane) ), size=5.0  )
pca.cor = p2 + geom_text_repel(aes(label = new.key$pid ), arrow = arrow(length = unit(0.01, 'npc')), box.padding = unit(1.5, 'lines')   ) +
  scale_colour_discrete(name  ="Drug") + scale_shape_discrete ( name="Lane") + ggtitle("Corrected for cell line")

library(ggplot2)
fte_theme <- function() {
  
  # Generate the colors for the chart procedurally with RColorBrewer
  palette <- brewer.pal("Greys", n=9)
  color.background = palette[2]
  color.grid.major = palette[3]
  color.axis.text = palette[6]
  color.axis.title = palette[7]
  color.title = palette[9]
  
  # Begin construction of chart
  theme_bw(base_size=9) +
    
    # Set the entire chart region to a light gray color
    theme(panel.background=element_rect(fill=color.background, color=color.background)) +
    theme(plot.background=element_rect(fill=color.background, color=color.background)) +
    theme(panel.border=element_rect(color=color.background)) +
    
    # Format the grid
    theme(panel.grid.major=element_line(color=color.grid.major,size=.25)) +
    theme(panel.grid.minor=element_blank()) +
    theme(axis.ticks=element_blank()) +
    
    # Format the legend, but hide by default
    theme(legend.position="none") +
    theme(legend.background = element_rect(fill=color.background)) +
    theme(legend.text = element_text(size=7,color=color.axis.title)) +
    
    # Set title and axis labels, and format these and tick marks
    theme(plot.title=element_text(color=color.title, size=10, vjust=1.25)) +
    theme(axis.text.x=element_text(size=7,color=color.axis.text)) +
    theme(axis.text.y=element_text(size=7,color=color.axis.text)) +
    theme(axis.title.x=element_text(size=8,color=color.axis.title, vjust=0)) +
    theme(axis.title.y=element_text(size=8,color=color.axis.title, vjust=1.25)) +
    
    # Plot margins
    theme(plot.margin = unit(c(0.35, 0.2, 0.3, 0.35), "cm"))
}


theme_538 <- theme(panel.background=element_rect(fill="#F0F0F0"), 
                   plot.background=element_rect(fill="#F0F0F0"), 
                   panel.border=element_rect(colour="#F0F0F0"),
                   panel.grid.major=element_line(colour="#D0D0D0",size=.75), axis.ticks=element_blank(),
                   legend.position="none",  plot.title=element_text(face="bold",hjust=-.08,vjust=2,colour="#3C3C3C",size=20),
                   axis.text.x=element_text(size=11,colour="#535353",face="bold"),
                   axis.text.y=element_text(size=11,colour="#535353",face="bold"),
                   axis.title.y=element_text(size=11,colour="#535353",face="bold",vjust=1.5),
                   axis.title.x=element_text(size=11,colour="#535353",face="bold",vjust=-.5),
                   plot.margin = unit(c(1, 1, .5, .7), "cm")) 

theme_538_legend <- theme(panel.background=element_rect(fill="#F0F0F0"), 
                   plot.background=element_rect(fill="#F0F0F0"), 
                   panel.border=element_rect(colour="#F0F0F0"),
                   panel.grid.major=element_line(colour="#D0D0D0",size=.75), axis.ticks=element_blank(),
                   legend.position="right",  plot.title=element_text(face="bold",hjust=-.08,vjust=2,colour="#3C3C3C",size=20),
                   axis.text.x=element_text(size=11,colour="#535353",face="bold"),
                   axis.text.y=element_text(size=11,colour="#535353",face="bold"),
                   axis.title.y=element_text(size=11,colour="#535353",face="bold",vjust=1.5),
                   axis.title.x=element_text(size=11,colour="#535353",face="bold",vjust=-.5),
                   plot.margin = unit(c(1, 1, .5, .7), "cm")) 

theme_538_bar <- theme(panel.background=element_rect(fill="#F0F0F0"), 
                   plot.background=element_rect(fill="#F0F0F0"), 
                  # panel.border=element_rect(colour="#F0F0F0"),
                   panel.grid.major=element_line(colour="#D0D0D0",size=.75), axis.ticks=element_blank(),
                   legend.position="none",  plot.title=element_text(face="bold",hjust=-.08,vjust=2,colour="#3C3C3C",size=20),
                   axis.text.x=element_text(size=11,colour="#535353",face="bold", angle = 45, hjust = 1),
                   axis.text.y=element_text(size=11,colour="#535353",face="bold"),
                   axis.title.y=element_text(size=11,colour="#535353",face="bold",vjust=1.5),
                   axis.title.x=element_text(size=11,colour="#535353",face="bold",vjust=-.5),
                   plot.margin = unit(c(1, 1, .5, .7), "cm")) 
# Set the default theme
# Modify it with theme()
tt <- theme(axis.text = element_text(size=16, colour=NULL),
            axis.title = element_text(size = 20, colour= NULL),
            plot.title = element_text(size = 30),
            # axis.line = element_line(size = 1, colour = "black"),
            axis.text = element_text(colour = "blue"),
            axis.ticks = element_line(size = 2),
            legend.background = element_rect(fill = "grey"),
            legend.position = c(0.13, 0.88),
            panel.grid.major = element_line(colour = "grey80"),
            panel.grid.minor = element_blank(),
            panel.background = element_rect(fill = "white"),
            panel.border = element_rect(colour = "black", fill=NA, size=1),
            text = element_text(size = 12, family = 'Times New Roman')
)

bcn_data_theme <- tt

# Calculates the geodesic distance between two points specified by radian latitude/longitude using the
# Haversine formula (hf)
# http://www.r-bloggers.com/great-circle-distance-calculations-in-r/
gcd.hf <- function(long1, lat1, long2, lat2) {
  R <- 6371 # Earth mean radius [km]
  delta.long <- (long2 - long1)
  delta.lat <- (lat2 - lat1)
  a <- sin(delta.lat/2)^2 + cos(lat1) * cos(lat2) * sin(delta.long/2)^2
  c <- 2 * asin(min(1,sqrt(a)))
  d = R * c
  return(d) # Distance in km
}

#####
# Define simple mathematical function
# for getting euclidean distance in kilometers
#####
get_distance <- function(lon1, 
                         lat1, 
                         lon2, 
                         lat2){
  rad <- pi/180
  a1 <- lat1 * rad
  a2 <- lon1 * rad
  b1 <- lat2 * rad
  b2 <- lon2 * rad
  dlon <- b2 - a2
  dlat <- b1 - a1
  a <- (sin(dlat/2))^2 + cos(a1) * cos(b1) * (sin(dlon/2))^2
  c <- 2 * atan2(sqrt(a), sqrt(1 - a))
  R <- 6378.145
  d <- R * c
  return(d)
}

#####
# CONFIDENCE INTERVALS ON PROPORTIONS
##### 

##### CONFIDENCE INTERVALS ON PROPORTIONS
# (https://aghaynes.wordpress.com/2014/04/09/calculating-confidence-intervals-for-proportions/)
simpasym <- function(n, p, z=1.96, cc=TRUE){
  out <- list()
  if(cc){
    out$lb <- p - z*sqrt((p*(1-p))/n) - 0.5/n
    out$ub <- p + z*sqrt((p*(1-p))/n) + 0.5/n
  } else {
    out$lb <- p - z*sqrt((p*(1-p))/n)
    out$ub <- p + z*sqrt((p*(1-p))/n)
  }
  out
}

##### CONFIDENCE INTERVALS ON PROPORTIONS
# (https://aghaynes.wordpress.com/2014/04/09/calculating-confidence-intervals-for-proportions/)
simpasym <- function(n, p, z=1.96, cc=TRUE){
  out <- list()
  if(cc){
    out$lb <- p - z*sqrt((p*(1-p))/n) - 0.5/n
    out$ub <- p + z*sqrt((p*(1-p))/n) + 0.5/n
  } else {
    out$lb <- p - z*sqrt((p*(1-p))/n)
    out$ub <- p + z*sqrt((p*(1-p))/n)
  }
  out
}

##### FARENHEIT TO CELSIUS
# http://swcarpentry.github.io/r-novice-inflammation/08-making-packages-R.html
f2c <- function(temp) {
  #Converts Fahrenheit to Celsius using fahr_to_kelvin() and kelvin_to_celsius()
  fahr_to_kelvin <- function(temp) {
    #Converts Fahrenheit to Kelvin
    kelvin <- ((temp - 32) * (5/9)) + 273.15
    kelvin
  }
  kelvin_to_celsius <- function(temp) {
    #Converts Kelvin to Celsius
    Celsius <- temp - 273.15
    Celsius
  }
  temp_k <- fahr_to_kelvin(temp)
  result <- kelvin_to_celsius(temp_k)
  result
}

#####
# MULTIPLOT
#####
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
# Front-end needs ---------------------------------------------------------
# Package install and load
  #install.packages('akima')
  library(akima)

# Data read and write -----------------------------------------------------
# Read in all of the data files for the PNR Run to make simulation results
  # Make an empty list to hold the data
    l = length(list.files(pattern = "BayesRes*"))
    d = vector(mode='list',
      length = l)
  # Set up a progress meter
    pb <- txtProgressBar(min = 0, max = l, style = 3, char = "+")
  # Read in each text file and add to corresponding element of list 'd'
    for(i in 1:l){
      d[[i]] = read.csv(list.files(pattern = "BayesRes")[i])
      Sys.sleep(0.001)
      setTxtProgressBar(pb, i)
    }
  # Close the progress meter
    close(pb)

  # Put all of the data into a single df
    pdata = do.call(rbind, lapply(d, data.frame))

  # Save it as an R object
    save(pdata, file = 'pdata.rda')

  # Read in the data file.
    load("pdata.rda")

  # Look at the names to make sure it's all there
    names(pdata)
    
  # Round off detection
    pdata$Pd = round(pdata$Pd, 2)
    
# Contour plots of abundance ----------------------------------------------
# Make a dataframe that can be ordered for the contour plot
	dat = data.frame(x=pdata$Nm,	y=pdata$Pd,	z=pdata$SE)

  # Order the data frame and remove NA values
  	dat = na.omit(dat[with(dat, order(x, y)), ])

  # Interpolate values of 'z' at evenly spaced values of x and y.
  	im = with(dat,
  	          interp(x, y, z, 
  	                 duplicate='mean',
  	                 nx=10,#length(unique(dat$x)),
  	                 ny=10))#length(unique(dat$y))))

# Make the contour plot
    par(mar=c(5, 5.2, 1, 10))
  # filled.contour is the function that actually makes the contour plot
	  filled.contour(
	    im$x,                                 # The variable to be displayed on the x-axis
	    im$y,                                 # The variable to be displayed on the y-axis
	    im$z,                                 # The response variable you wish to plot
	    levels = c(seq(0,max(im$z, na.rm = T),.005)),
	    col=(gray.colors(11)),             # Could also choose 'grey.colors' or 'topo.colors'. If you want the ramp to go the other way, just delete the 'rev'. Note that you will need to change the 20 in parentheses to match the number of levels that you actually have or want to display.
	    main = '',                            # I don't like in-figure titles. You can add one, though. You will, however, need to change the 'mar' argument in the call to par above.
	    ylim=c(min(dat$y),max(dat$y)),                         # Set max and min of y-axis to your data range
	    xlim=c(min(dat$x),max(dat$x)),                         # Set max and min of x-axis to your data range
	    xlab="Number released",               # Change the words in the quotes to change the x-axis label
	    cex.lab=1.5,                          # This makes the labels 1.5x larger than default
	    plot.axes = {                         # This argument tells R to print the axes, but increas the size
	      contour(                            # This is the line that adds the contour lines
	        im$x,                             # The variable to be displayed on the x-axis
	        im$y,                             # The variable to be displayed on the y-axis
	        im$z,                             # The response variable you wish to plot
	        levels = c(seq(0,max(im$z, na.rm = T),0.01)),                     # This number needs to match the one in 'col' on line 102
	        drawlabels = FALSE,               # The labels are realy ugly
	        col = c(rep(rgb(0,0,0, alpha=0.05), 12),
	                'black',
	                rep(rep(rgb(0,0,0, alpha=0.05), 13))
	                ),
	        lwd = c(rep(1,12),2,rep(1,12)),
	        lty = c(rep(1,12),2,rep(1,12)),
	        add = TRUE                        # Add the lines to the current plot
	      );                                  # Close the call to the contour line function
	      axis(1, cex.axis=1.25);             # X axis tick marks & tick labels
	      axis(2, cex.axis=1.25)              # Y axis tick marks & tick labels
	    }                                     # Close the argument plot.axes
	  )                                       # Close the call to filled.contour
	# Finally, add a label for the y-axis
    mtext(side = 2, "Detection probability (p)", line=4, cex.lab=1.5, cex=1.5)

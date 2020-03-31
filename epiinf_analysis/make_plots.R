library(bdskytools)
library(lubridate)

dates <- ymd(read.table("dates.txt", header=F)[[1]])

plotDateSkyline <- function(times, hpd_matrix, endDate,
                            add=FALSE,
                            addWithSecondAxis=FALSE,
                            xlim=NA, ylim=NA,
                            ylab=NA,
                            col=rgb(1,0.5,0),
                            ...) {

    dates <- endDate - times*365

    if (length(ylim)==1 && is.na(ylim))
        ylim <- c(0, max(hpd_matrix[3,]))

    if (length(xlim)==1 && is.na(xlim))
        xlim <- c(min(dates), endDate)

    if (!add) {
        if (addWithSecondAxis) {
            par(new=TRUE)
            plot(1, type='n', xlim=xlim, ylim=ylim, xlab=NA, ylab=NA, xaxt="n", yaxt="n", ...)
            axis(side=4, col=col, col.axis=col)
            mtext(ylab, side=4, line=par("mgp")[1], col=col)
        } else {
            par(new=FALSE)
            plot(1, type='n', xlim=xlim, ylim=ylim, xlab=NA, ylab=NA, xaxt="n", yaxt="n", ...)
            axis(side=2, col=col, col.axis=col)
            axis.Date(1, x=dates)
            mtext(ylab, side=2, line=par("mgp")[1], col=col)
        }
    }

    c <- col2rgb(col)/255
    col_faded <- rgb(c[1], c[2], c[3], 0.5)

    polygon(c(dates, rev(dates)), c(hpd_matrix[1,], rev(hpd_matrix[3,])),
            col=col_faded, border=NA)

    lines(dates, hpd_matrix[2,], lwd=2, col=col)
}

gridSkylineNotEquidistant <- function(sky, origin, endtimes, timegrid) {
    grid <- matrix(NA, nrow=dim(sky)[1], ncol=length(timegrid))

    for (tidx in 1:length(timegrid)) {
        t <- timegrid[tidx]
        interval <- which.min(endtimes>t)
        grid[,tidx] <- sky[,interval]
    }

    return(grid)
}

### BDSKY plots

for (bu in c("36.5", "30.4")) {

    lf <- readLogfile(paste0("Results/BD_epi.clock_8e-4.bu_",bu,".1.log"))

    Re_sky <- getSkylineSubset(lf, "reproductiveNumber")
    Re_hpd <- getMatrixHPD(Re_sky)

    timegrid <- seq(0,max(lf$origin), length.out=101)
    Re_gridded <- gridSkyline(Re_sky, lf$origin, timegrid)
    Re_gridded_hpd <- getMatrixHPD(Re_gridded)

    s_period <- decimal_date(max(dates))-decimal_date(min(dates))
    s_sky <- getSkylineSubset(lf, "samplingProportion")
    s_gridded <- gridSkylineNotEquidistant(s_sky, lf$origin, c(s_period, 0), timegrid)
    s_gridded_hpd <- getMatrixHPD(s_gridded)

    ## plotSkyline(max(dates)-timegrid*365, Re_gridded_hpd, type="smooth")
    png(paste0("figures/skyline_bu_",bu,".png"), width=800, height=600, pointsize=20) 
    par(mar=c(2,3,3,3), mgp=c(2, 0.5, 0))
    plotDateSkyline(timegrid, Re_gridded_hpd, max(dates),
                    ylab=expression(R[e]),
                    main=paste0("Clock rate 8e-4/s/y, Become uninfectious rate ", bu, "/y"))

    plotDateSkyline(timegrid, s_gridded_hpd, max(dates),
                    col=rgb(0.5,0.5,1.0),
                    ylab=expression(s),
                    ylim=c(0,1),
                    addWithSecondAxis=TRUE)
    dev.off()
}

### Trajectory plots


source("~/code/beast_and_friends/EpiInf/scripts/plotTraj.R")

png("figures/trajectories_bu_36.5.png", width=800, height=600, pointsize=20)
par(mar=c(3,3,2,1), mgp=c(2,0.5,0))
plotTraj(dir("Results", "BD_epi.clock_8e-4.bu_36.5.*.traj", full.names=T),
         log='y', ylim=c(1,1e4),
         presentTime = decimal_date(max(dates)),
         timesAreCalendarYears = TRUE,
         main="Clock rate 8e-4/s/y, Become uninfectious rate 36.5/y")
legend("topleft", inset=0.05, c("Median", "95% CI"), lwd=2, lty=c(1,2))
dev.off()

png("figures/trajectories_bu_30.4.png", width=800, height=600, pointsize=20)
plotTraj(dir("Results", "BD_epi.clock_8e-4.bu_30.4.*.traj", full.names=T),
         log='y', ylim=c(1,1e4),
         presentTime = decimal_date(max(dates)),
         timesAreCalendarYears = TRUE,
         main="Clock rate 8e-4/s/y, Become uninfectious rate 30.4/y")
legend("topleft", inset=0.05, c("Median", "95% CI"), lwd=2, lty=c(1,2))
dev.off()
